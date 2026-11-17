/*
 * east_tm: exact results for the discrete-time (MCS) softened East dynamics.
 *
 *   mode finite : exact moments of K and E for the TPS trajectory ensemble
 *                 with tobs sweeps (use these to check east_tps)
 *   mode asym   : long-time limit from the largest eigenvalue of the tilted
 *                 sweep transfer matrix (use these to check east_clone)
 *   mode pk     : exact distribution P(K) of the TPS ensemble at each s, and
 *                 optionally the equal-weight point (keywords kcut, hist_bins,
 *                 smooth_res as in east_mbar); small L * tobs only
 *
 * usage: ./bin/east_tm [input_file]      (default input file: system.inp)
 */
#include "input.h"
#include "peaks.h"
#include "tm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const double *logp; /* ln P(K) at s = 0 */
    long Kmax;
    double Klo, Khi;
    int nbins;
} PkCtx;

/* normalized probabilities P_s(K) into p */
static void pk_at(const PkCtx *c, double s, double *p)
{
    double mx = -INFINITY, sum = 0.0;
    long K;
    for (K = 0; K <= c->Kmax; K++)
        if (c->logp[K] - s * K > mx) mx = c->logp[K] - s * K;
    for (K = 0; K <= c->Kmax; K++) {
        p[K] = exp(c->logp[K] - s * K - mx);
        sum += p[K];
    }
    for (K = 0; K <= c->Kmax; K++) p[K] /= sum;
}

static double pk_prob_below(double s, double Kcut, void *ctx)
{
    const PkCtx *c = (const PkCtx *)ctx;
    double *p = (double *)malloc(sizeof(double) * (size_t)(c->Kmax + 1)), acc = 0.0;
    long K;
    pk_at(c, s, p);
    for (K = 0; K <= c->Kmax && K < Kcut; K++) acc += p[K];
    free(p);
    return acc;
}

static void pk_hist(double s, double *h, void *ctx)
{
    const PkCtx *c = (const PkCtx *)ctx;
    double *p = (double *)malloc(sizeof(double) * (size_t)(c->Kmax + 1));
    const double w = (c->Khi - c->Klo) / c->nbins;
    long K;
    int b;
    pk_at(c, s, p);
    for (b = 0; b < c->nbins; b++) h[b] = 0.0;
    for (K = 0; K <= c->Kmax; K++) {
        b = (int)floor((K - c->Klo) / w);
        if (b >= 0 && b < c->nbins) h[b] += p[K] / w;
    }
    free(p);
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "system.inp";
    const char *mode = "finite", *out = NULL, *basis = "k0", *v;
    Input in;
    int L = 0, tobs = 0, ns = 0, ng = 0, is, ig, hist_bins = 60;
    double smooth_res = 30.0;
    int use_kcut = 0, auto_cut = 0;
    double T = 1.0, eps = 0.0, fd_h = 1e-4, kcut = NAN;
    double *ss = NULL, *gs = NULL;
    TMSolver tm;
    FILE *f = stdout;

    if (argc > 1 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        printf("usage: ./bin/east_tm [input_file]   (default: system.inp)\n"
               "keywords: mode (finite|asym|pk) L T eps s g tobs basis (k0|full) fd_h output\n"
               "          pk mode: kcut (value|auto) hist_bins smooth_res\n");
        return 0;
    }
    if (input_read(&in, path) != 0) return 1;
    if (input_str(&in, "smooth")) {
        fprintf(stderr, "error: keyword 'smooth' was replaced by 'smooth_res'\n");
        return 1;
    }
    if ((v = input_str(&in, "mode"))) mode = v;
    if (input_int(&in, "L", 1, &L) < 0 || input_double(&in, "T", 0, &T) < 0 ||
        input_double(&in, "eps", 1, &eps) < 0 || input_list(&in, "s", 1, &ss, &ns) < 0 ||
        input_list(&in, "g", 1, &gs, &ng) < 0 || input_double(&in, "fd_h", 0, &fd_h) < 0) {
        input_free(&in);
        return 1;
    }
    if ((!strcmp(mode, "finite") || !strcmp(mode, "pk")) && input_int(&in, "tobs", 1, &tobs) < 0) {
        input_free(&in);
        return 1;
    }
    if (!strcmp(mode, "pk")) {
        if (input_int(&in, "hist_bins", 0, &hist_bins) < 0 ||
            input_double(&in, "smooth_res", 0, &smooth_res) < 0) {
            input_free(&in);
            return 1;
        }
        if ((v = input_str(&in, "kcut"))) {
            use_kcut = 1;
            if (!strcmp(v, "auto")) auto_cut = 1;
            else if (input_double(&in, "kcut", 1, &kcut) < 0) return 1;
        }
    }
    if ((v = input_str(&in, "basis"))) basis = v;
    if ((v = input_str(&in, "output"))) out = v;
    input_warn_unused(&in);
    if (strcmp(mode, "finite") && strcmp(mode, "asym") && strcmp(mode, "pk")) {
        fprintf(stderr, "error: mode must be finite, asym or pk\n");
        input_free(&in);
        return 1;
    }
    if (strcmp(basis, "k0") && strcmp(basis, "full")) {
        fprintf(stderr, "error: basis must be k0 or full\n");
        input_free(&in);
        return 1;
    }
    if (L > 20) {
        fprintf(stderr, "error: L = %d is too large for exact enumeration (use L <= 20)\n", L);
        input_free(&in);
        return 1;
    }
    /* the asymptotic strict-East calculation lives in the non-empty sector */
    if (tm_init(&tm, L, T, eps, !strcmp(basis, "k0"), !strcmp(mode, "asym") && eps == 0.0) !=
        0) {
        fprintf(stderr, "error: cannot initialize L=%d T=%g eps=%g\n", L, T, eps);
        input_free(&in);
        return 1;
    }
    if (out) {
        f = fopen(out, "w");
        if (!f) {
            fprintf(stderr, "error: cannot write %s\n", out);
            tm_free(&tm);
            input_free(&in);
            return 1;
        }
    }
    fprintf(f, "# soft_east discrete-time exact (%s), basis %s, dim %ld\n", mode, basis,
            tm.basis.dim);
    fprintf(f, "# L %d  T %.12g  c %.12g  eps %.12g\n", L, T, tm.model.c, eps);
    fprintf(f, "# MCS = L random-site attempts; weight exp(-s K - g E)\n");

    if (!strcmp(mode, "finite")) {
        fprintf(f, "# tobs %d sweeps; E sums the tobs+1 slices; x_0 from equilibrium\n", tobs);
        fprintf(f, "# k = <K>/(L tobs), e = <E>/(L (tobs+1)), chi_XY = Cov(X,Y)/(L tobs)\n");
        fprintf(f, "# %10s %12s %20s %18s %18s %18s %18s %18s %16s %16s %16s %16s %16s\n", "g",
                "s", "lnZ", "<K>", "<E>", "VarK", "VarE", "CovKE", "k", "e", "chi_KK",
                "chi_EE", "chi_KE");
        for (ig = 0; ig < ng; ig++) {
            for (is = 0; is < ns; is++) {
                TMFinite r;
                const double nk = (double)L * tobs, ne = (double)L * (tobs + 1);
                tm_finite(&tm, ss[is], gs[ig], tobs, &r);
                fprintf(f, "%12.6g %12.8g %20.12e %18.10e %18.10e %18.10e %18.10e %18.10e "
                           "%16.10e %16.10e %16.10e %16.10e %16.10e\n",
                        gs[ig], ss[is], r.lnZ, r.meanK, r.meanE, r.varK, r.varE, r.covKE,
                        r.meanK / nk, r.meanE / ne, r.varK / nk, r.varE / nk, r.covKE / nk);
                fflush(f);
            }
            if (ng > 1) fprintf(f, "\n");
        }
    } else if (!strcmp(mode, "pk")) {
        const long Kmax = (long)L * tobs;
        double *lp = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
        double *p = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
        fprintf(f, "# tobs %d sweeps; P(K) of the TPS ensemble (x_0 from equilibrium)\n", tobs);
        fprintf(f, "# columns: g s K k=K/(L tobs) P(K); blocks separated by blank lines\n");
        for (ig = 0; ig < ng; ig++) {
            PkCtx c;
            tm_pk(&tm, gs[ig], tobs, lp);
            c.logp = lp;
            c.Kmax = Kmax;
            c.Klo = -0.5;
            c.Khi = Kmax + 0.5;
            c.nbins = hist_bins;
            for (is = 0; is < ns; is++) {
                long K;
                pk_at(&c, ss[is], p);
                for (K = 0; K <= Kmax; K++)
                    fprintf(f, "%12.6g %12.8g %8ld %12.8f %16.8e\n", gs[ig], ss[is], K,
                            (double)K / Kmax, p[K]);
                fprintf(f, "\n\n");
            }
            if (use_kcut) {
                SplitProblem sp;
                double smin = ss[0], smax = ss[0], sstar, Kc = NAN;
                for (is = 1; is < ns; is++) {
                    if (ss[is] < smin) smin = ss[is];
                    if (ss[is] > smax) smax = ss[is];
                }
                sp.prob_below = pk_prob_below;
                sp.hist = pk_hist;
                sp.ctx = &c;
                sp.Klo = c.Klo;
                sp.Khi = c.Khi;
                sp.nbins = hist_bins;
                sp.smooth_res = smooth_res;
                if (auto_cut) {
                    double best = -1.0, s0 = ss[0];
                    for (is = 0; is < ns; is++) {
                        double m1 = 0.0, m2 = 0.0;
                        long K;
                        pk_at(&c, ss[is], p);
                        for (K = 0; K <= Kmax; K++) {
                            m1 += p[K] * K;
                            m2 += p[K] * (double)K * K;
                        }
                        if (m2 - m1 * m1 > best) {
                            best = m2 - m1 * m1;
                            s0 = ss[is];
                        }
                    }
                    sstar = equal_weight_auto(&sp, s0, smin, smax, &Kc);
                } else {
                    Kc = kcut * Kmax;
                    sstar = equal_weight_fixed(&sp, Kc, smin, smax);
                }
                fprintf(f, "# s_equal g %.6g s* %.10f kcut %.8f%s\n", gs[ig], sstar, Kc / Kmax,
                        isnan(sstar) ? "  (no valley or not bracketed by s)" : "");
            }
            fflush(f);
        }
        free(lp);
        free(p);
    } else {
        if (eps == 0.0)
            fprintf(f, "# eps = 0: non-empty sector; the empty state has theta = 0\n");
        fprintf(f, "# theta = ln lambda_1 per MCS (whole system); k, e per site per MCS\n");
        fprintf(f, "# chi_KK = -dk/ds by central difference with h = %g (0: skipped)\n", fd_h);
        fprintf(f, "# %10s %12s %20s %14s %18s %18s %18s %9s %3s\n", "g", "s", "theta",
                "lam2/lam1", "k", "e", "chi_KK", "resid", "ok");
        for (ig = 0; ig < ng; ig++) {
            for (is = 0; is < ns; is++) {
                TMAsym r, rp, rm;
                double chi = NAN;
                tm_asym(&tm, ss[is], gs[ig], &r);
                if (fd_h > 0.0) {
                    tm_asym(&tm, ss[is] + fd_h, gs[ig], &rp);
                    tm_asym(&tm, ss[is] - fd_h, gs[ig], &rm);
                    chi = -(rp.k - rm.k) / (2.0 * fd_h);
                }
                fprintf(f, "%12.6g %12.8g %20.12e %14.8f %18.10e %18.10e %18.10e %9.1e %3d\n",
                        gs[ig], ss[is], r.theta, r.lam_ratio, r.k, r.e, chi, r.resid, r.ok);
                fflush(f);
            }
            if (ng > 1) fprintf(f, "\n");
        }
    }
    if (f != stdout) fclose(f);
    tm_free(&tm);
    free(ss);
    free(gs);
    input_free(&in);
    return 0;
}
