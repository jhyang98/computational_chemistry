/*
 * Statistical validation of the discrete-time simulations against exact
 * transfer-matrix results (fixed seeds, so the output is reproducible).
 *
 *  1. one-sweep kernel: histogram of x_1 from fixed x_0 vs exact row of M^L
 *     (chi-square), and the mean number of flips
 *  2. TPS: <K>, <E>, Var K, Var E vs exact finite-tobs values
 *  3. cloning: theta_d vs exact largest eigenvalue; k, e vs exact
 *  4. MBAR over TPS windows: k, e, chi at intermediate s vs exact
 *  5. umbrella sampling in K + MBAR: P(K < K_c), moments at several s and
 *     the equal-weight point vs exact P(K) from the transfer matrix; <K> and
 *     <E> within each peak at the exact equal-weight point
 *  6. the same with replica-exchange umbrella sampling (east_reus); also
 *     swap acceptance, round trips and overlap consistency
 *
 * A check passes when |z| = |estimate - exact| / error <= 4.
 */
#include "clone.h"
#include "dyn.h"
#include "mbar.h"
#include "peaks.h"
#include "reus.h"
#include "rng.h"
#include "tm.h"
#include "tps.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ZMAX 4.0

static int n_pass = 0, n_fail = 0;

static void check_true(const char *name, int cond)
{
    if (cond) n_pass++;
    else n_fail++;
    printf("[%s] %s\n", cond ? "PASS" : "FAIL", name);
}

static void zcheck(const char *name, double got, double err, double exact)
{
    double z = (got - exact) / err;
    int ok = fabs(z) <= ZMAX && !isnan(z);
    if (ok) n_pass++;
    else n_fail++;
    printf("[%s] %-52s est % .8e +- %.2e  exact % .8e  z %+6.2f\n", ok ? "PASS" : "FAIL", name, got,
           err, exact, z);
}

/* exact mean number of flips in one sweep from configuration x (full basis) */
static double exact_mean_flips(TMSolver *t, uint32_t x)
{
    const long n = t->basis.dim;
    const int L = t->model.L;
    double *a = (double *)calloc((size_t)n, sizeof(double));
    double *as = (double *)calloc((size_t)n, sizeof(double));
    double *Aa = (double *)malloc(sizeof(double) * (size_t)n);
    double *Aas = (double *)malloc(sizeof(double) * (size_t)n);
    double tot = 0.0, tots = 0.0;
    long i;
    int k;
    a[t->basis.index[x]] = 1.0 / t->psi0[t->basis.index[x]];
    for (k = 0; k < L; k++) {
        tgen_apply_A(&t->tg, a, Aa);
        tgen_apply_A(&t->tg, as, Aas);
        for (i = 0; i < n; i++) {
            double r = t->tg.esc[i];
            double na = a[i] + (Aa[i] - r * a[i]) / L;
            double nas = as[i] + (Aas[i] - r * as[i]) / L - Aa[i] / L; /* M' = -A/L at s=0 */
            a[i] = na;
            as[i] = nas;
        }
    }
    for (i = 0; i < n; i++) {
        tot += t->psi0[i] * a[i];
        tots += t->psi0[i] * as[i];
    }
    free(a);
    free(as);
    free(Aa);
    free(Aas);
    return -tots / tot;
}

static void validate_kernel(void)
{
    const int L = 5;
    const double eps = 0.2;
    const long nsamp = 2000000;
    uint32_t xs[] = {0u, 5u, 26u, 31u};
    TMSolver tm;
    Dyn d;
    Rng r;
    int ix;
    printf("\n== 1. one-sweep kernel (L = %d, eps = %g, %ld sweeps per start) ==\n", L, eps, nsamp);
    tm_init(&tm, L, 1.0, eps, 0, 0);
    dyn_init(&d, L, 1.0, eps);
    rng_seed(&r, 20261106ULL);
    for (ix = 0; ix < 4; ix++) {
        double prob[32], cnt[32], chi2 = 0.0, other_e = 0.0, other_o = 0.0, sumk = 0.0, sumk2 = 0.0;
        unsigned char x[5];
        long m;
        int q, dof = -1, i;
        char name[128];
        tm_mcs_row(&tm, xs[ix], prob);
        memset(cnt, 0, sizeof(cnt));
        for (m = 0; m < nsamp; m++) {
            int n = 0, k;
            uint32_t y = 0;
            for (i = 0; i < L; i++) {
                x[i] = (unsigned char)((xs[ix] >> i) & 1u);
                n += x[i];
            }
            k = dyn_mcs(&d, x, &n, &r);
            for (i = 0; i < L; i++) y |= (uint32_t)x[i] << i;
            cnt[y] += 1.0;
            sumk += k;
            sumk2 += (double)k * k;
        }
        for (q = 0; q < 32; q++) {
            double ex = prob[q] * nsamp;
            if (ex >= 20.0) {
                chi2 += (cnt[q] - ex) * (cnt[q] - ex) / ex;
                dof++;
            } else {
                other_e += ex;
                other_o += cnt[q];
            }
        }
        if (other_e >= 20.0) {
            chi2 += (other_o - other_e) * (other_o - other_e) / other_e;
            dof++;
        }
        {
            /* Wilson-Hilferty: chi2/dof -> approximately normal */
            double kk = dof, z = (pow(chi2 / kk, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * kk))) /
                                 sqrt(2.0 / (9.0 * kk));
            int ok = z <= ZMAX;
            if (ok) n_pass++;
            else n_fail++;
            printf("[%s] x0=%-2u chi2 = %8.2f for %2d dof (z %+5.2f)\n", ok ? "PASS" : "FAIL",
                   (unsigned)xs[ix], chi2, dof, z);
        }
        {
            double mk = sumk / nsamp, vk = sumk2 / nsamp - mk * mk;
            sprintf(name, "mean flips per sweep from x0=%u", (unsigned)xs[ix]);
            zcheck(name, mk, sqrt(vk / nsamp), exact_mean_flips(&tm, xs[ix]));
        }
    }
    tm_free(&tm);
}

static void validate_tps(void)
{
    struct {
        double eps, s, g;
    } cs[] = {{0.1, 0.0, 0.0}, {0.1, 0.03, 0.05}, {0.01, 0.02, 0.0}, {0.01, -0.02, 0.1}};
    const int L = 8, tobs = 40;
    int i;
    printf("\n== 2. TPS vs exact finite-tobs moments (L = %d, tobs = %d) ==\n", L, tobs);
    for (i = 0; i < 4; i++) {
        TPSParams p;
        TPSResult r;
        TMSolver tm;
        TMFinite ex;
        char name[128];
        tps_default_params(&p);
        p.L = L;
        p.tobs = tobs;
        p.eps = cs[i].eps;
        p.s = cs[i].s;
        p.g = cs[i].g;
        p.seed = 1000 + i;
        p.n_equil = 5000;
        p.n_prod = 600000;
        p.sample_every = 10;
        p.nblocks = 50;
        tps_run(&p, &r);
        tm_init(&tm, L, 1.0, p.eps, 1, 0);
        tm_finite(&tm, p.s, p.g, tobs, &ex);
        sprintf(name, "<K>   eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.st.mean_x, r.st.err_mean_x, ex.meanK);
        sprintf(name, "<E>   eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.st.mean_y, r.st.err_mean_y, ex.meanE);
        sprintf(name, "Var K eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.st.var_x, r.st.err_var_x, ex.varK);
        sprintf(name, "Var E eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.st.var_y, r.st.err_var_y, ex.varE);
        sprintf(name, "Cov   eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.st.cov, r.st.err_cov, ex.covKE);
        tps_result_free(&r);
        tm_free(&tm);
    }
}

static void validate_clone(void)
{
    struct {
        double eps, s, g;
    } cs[] = {{0.1, 0.05, 0.0}, {0.1, 0.0, 0.1}, {0.01, 0.02, 0.05}, {0.01, -0.05, 0.0}};
    const int L = 8;
    int i;
    printf("\n== 3. cloning vs exact long-time limit (L = %d) ==\n", L);
    for (i = 0; i < 4; i++) {
        CloneParams p;
        CloneResult r;
        TMSolver tm;
        TMAsym ex;
        char name[128];
        clone_default_params(&p);
        p.L = L;
        p.eps = cs[i].eps;
        p.s = cs[i].s;
        p.g = cs[i].g;
        p.n_clones = 2000;
        p.tobs = 1000;
        p.t_skip = 100;
        p.t_tail = 300;
        p.n_runs = 10;
        p.seed = 3000 + i;
        clone_run(&p, &r, NULL, NULL, NULL);
        tm_init(&tm, L, 1.0, p.eps, 1, 0);
        tm_asym(&tm, p.s, p.g, &ex);
        sprintf(name, "theta_d eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.theta, r.theta_err, ex.theta);
        sprintf(name, "k_d     eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.k, r.k_err, ex.k);
        sprintf(name, "e_d     eps=%g s=%g g=%g", p.eps, p.s, p.g);
        zcheck(name, r.e, r.e_err, ex.e);
        tm_free(&tm);
    }
}

static void validate_mbar(void)
{
    const int L = 8, tobs = 40, nwin = 5, nboot = 100, block = 100;
    const double eps = 0.01, g = 0.0;
    double sw[5] = {0.0, 0.01, 0.02, 0.03, 0.04}, st[4] = {0.005, 0.015, 0.025, 0.035};
    TPSResult res[5];
    MBAR m;
    TMSolver tm;
    Rng r;
    long ntot = 0, off = 0;
    int w, i, b;
    double acc1[4][3], acc2[4][3];
    printf("\n== 4. MBAR over TPS windows vs exact (L = %d, tobs = %d, eps = %g) ==\n", L, tobs,
           eps);
    for (w = 0; w < nwin; w++) {
        TPSParams p;
        tps_default_params(&p);
        p.L = L;
        p.tobs = tobs;
        p.eps = eps;
        p.s = sw[w];
        p.g = g;
        p.seed = 5000 + w;
        p.n_equil = 5000;
        p.n_prod = 400000;
        p.sample_every = 10;
        tps_run(&p, &res[w]);
        ntot += res[w].nsamp;
    }
    mbar_alloc(&m, nwin, ntot);
    for (w = 0; w < nwin; w++) {
        m.ws[w] = sw[w];
        m.wg[w] = g;
        m.N[w] = res[w].nsamp;
        memcpy(m.K + off, res[w].K, sizeof(double) * (size_t)res[w].nsamp);
        memcpy(m.E + off, res[w].E, sizeof(double) * (size_t)res[w].nsamp);
        off += res[w].nsamp;
        tps_result_free(&res[w]);
    }
    mbar_guess(&m);
    mbar_solve(&m, 1e-10, 200);
    memset(acc1, 0, sizeof(acc1));
    memset(acc2, 0, sizeof(acc2));
    rng_seed(&r, 777ULL);
    for (b = 0; b < nboot; b++) {
        MBAR mb;
        mbar_resample(&m, &mb, block, &r);
        mbar_solve(&mb, 1e-10, 200);
        for (i = 0; i < 4; i++) {
            MBARObs o;
            double v[3];
            int q;
            mbar_observe(&mb, st[i], g, &o);
            v[0] = o.meanK;
            v[1] = o.meanE;
            v[2] = o.varK;
            for (q = 0; q < 3; q++) {
                acc1[i][q] += v[q];
                acc2[i][q] += v[q] * v[q];
            }
        }
        mbar_free(&mb);
    }
    tm_init(&tm, L, 1.0, eps, 1, 0);
    for (i = 0; i < 4; i++) {
        MBARObs o;
        TMFinite ex;
        double err[3];
        int q;
        char name[128];
        for (q = 0; q < 3; q++) {
            double mu = acc1[i][q] / nboot;
            err[q] = sqrt((acc2[i][q] / nboot - mu * mu) * nboot / (nboot - 1.0));
        }
        mbar_observe(&m, st[i], g, &o);
        tm_finite(&tm, st[i], g, tobs, &ex);
        sprintf(name, "MBAR <K>  at s=%g (between windows)", st[i]);
        zcheck(name, o.meanK, err[0], ex.meanK);
        sprintf(name, "MBAR <E>  at s=%g (between windows)", st[i]);
        zcheck(name, o.meanE, err[1], ex.meanE);
        sprintf(name, "MBAR VarK at s=%g (between windows)", st[i]);
        zcheck(name, o.varK, err[2], ex.varK);
    }
    tm_free(&tm);
    mbar_free(&m);
}

/* ---- 5. umbrella sampling ---- */

typedef struct {
    const double *logp;
    long Kmax;
    double Klo, Khi;
    int nbins;
} ExactPk;

static void exact_ps(const ExactPk *c, double s, double *p)
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

static double exact_below(double s, double Kcut, void *ctx)
{
    const ExactPk *c = (const ExactPk *)ctx;
    double *p = (double *)malloc(sizeof(double) * (size_t)(c->Kmax + 1)), acc = 0.0;
    long K;
    exact_ps(c, s, p);
    for (K = 0; K <= c->Kmax && K < Kcut; K++) acc += p[K];
    free(p);
    return acc;
}

static void exact_hist(double s, double *h, void *ctx)
{
    const ExactPk *c = (const ExactPk *)ctx;
    double *p = (double *)malloc(sizeof(double) * (size_t)(c->Kmax + 1));
    const double w = (c->Khi - c->Klo) / c->nbins;
    long K;
    int b;
    exact_ps(c, s, p);
    for (b = 0; b < c->nbins; b++) h[b] = 0.0;
    for (K = 0; K <= c->Kmax; K++) {
        b = (int)floor((K - c->Klo) / w);
        if (b >= 0 && b < c->nbins) h[b] += p[K] / w;
    }
    free(p);
}

/* Exact averages of K and E on the two sides of Kcut at (s, g = 0):
   out = {<K>_lo, <K>_hi, <E>_lo, <E>_hi}, lo = K < Kcut. The E averages use
   <E>_side = <E>(s = 0, g = 0) - d/dg ln sum_side P_g(K) exp(-s K) at g = 0
   (central difference with step h). */
static double side_logsum(const double *lp, long Kmax, double s, double Kcut, int hi)
{
    double mx = -INFINITY, sum = 0.0;
    long K;
    for (K = 0; K <= Kmax; K++)
        if ((K >= Kcut) == hi && lp[K] - s * K > mx) mx = lp[K] - s * K;
    for (K = 0; K <= Kmax; K++)
        if ((K >= Kcut) == hi) sum += exp(lp[K] - s * K - mx);
    return mx + log(sum);
}

static void exact_split(TMSolver *tm, const double *lp0, long Kmax, int tobs, double s,
                        double Kcut, double *out)
{
    const double h = 1e-5; /* truncation error ~ h^2: 1e-7 relative here */
    double *lpp = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
    double *lpm = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
    double *p = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
    double w[2] = {0.0, 0.0}, kk[2] = {0.0, 0.0};
    TMFinite f0;
    ExactPk c;
    long K;
    int side;
    c.logp = lp0;
    c.Kmax = Kmax;
    exact_ps(&c, s, p);
    for (K = 0; K <= Kmax; K++) {
        side = K >= Kcut;
        w[side] += p[K];
        kk[side] += p[K] * K;
    }
    out[0] = kk[0] / w[0];
    out[1] = kk[1] / w[1];
    tm_finite(tm, 0.0, 0.0, tobs, &f0);
    tm_pk(tm, h, tobs, lpp);
    tm_pk(tm, -h, tobs, lpm);
    for (side = 0; side < 2; side++)
        out[2 + side] = f0.meanE - (side_logsum(lpp, Kmax, s, Kcut, side) -
                                    side_logsum(lpm, Kmax, s, Kcut, side)) / (2.0 * h);
    free(lpp);
    free(lpm);
    free(p);
}

typedef struct {
    const MBAR *m;
} MbarPk;

static double mbar_below_cb(double s, double Kcut, void *ctx)
{
    return mbar_prob_below(((const MbarPk *)ctx)->m, s, 0.0, Kcut);
}

static void mbar_hist_cb(double s, double *h, void *ctx);
static double g_Klo, g_Khi;
static int g_nb;
static void mbar_hist_cb(double s, double *h, void *ctx)
{
    mbar_histogram(((const MbarPk *)ctx)->m, s, 0.0, g_Klo, g_Khi, g_nb, h);
}

static void umbrella_quantities(const MBAR *m, const double *sv, int ns, const double *Kc, int nc,
                                const SplitProblem *tmpl, double s_split, double K_split,
                                double *out)
{
    /* out: <K>(s), VarK(s) for each s; P(K<Kc) at sv[1] for each Kc; s*_fixed; s*_auto;
       <K>, <E> below and above K_split at s_split */
    MbarPk ctx;
    SplitProblem sp = *tmpl;
    int i, q = 0;
    double Kauto;
    ctx.m = m;
    sp.prob_below = mbar_below_cb;
    sp.hist = mbar_hist_cb;
    sp.ctx = &ctx;
    for (i = 0; i < ns; i++) {
        MBARObs o;
        mbar_observe(m, sv[i], 0.0, &o);
        out[q++] = o.meanK;
        out[q++] = o.varK;
    }
    for (i = 0; i < nc; i++) out[q++] = mbar_prob_below(m, sv[1], 0.0, Kc[i]);
    out[q++] = equal_weight_fixed(&sp, Kc[0], 0.0, 0.06);
    out[q++] = equal_weight_auto(&sp, sv[1], 0.0, 0.06, &Kauto);
    {
        MBARSplit o;
        mbar_observe_split(m, s_split, 0.0, K_split, &o);
        out[q++] = o.meanK_lo;
        out[q++] = o.meanK_hi;
        out[q++] = o.meanE_lo;
        out[q++] = o.meanE_hi;
    }
}

static void validate_umbrella(int use_reus)
{
    const int L = 8, tobs = 100, nwin = 8, nboot = 60, block = 200;
    const double eps = 0.003, s_run = 0.03, kspan = 0.25;
    const double delta = kspan / (nwin + 1);
    const long Kmax = (long)L * tobs;
    double sv[4] = {0.0, 0.03, 0.045, 0.06};
    double Kc[3];
    TPSResult res[8];
    MBAR m;
    TMSolver tm;
    ExactPk ec;
    SplitProblem sp;
    Rng r;
    double *lp = (double *)malloc(sizeof(double) * (size_t)(Kmax + 1));
    double est[4 * 2 + 3 + 2 + 4], acc1[4 * 2 + 3 + 2 + 4], acc2[4 * 2 + 3 + 2 + 4];
    double ex[4 * 2 + 3 + 2 + 4], s_split = NAN;
    int cnt[4 * 2 + 3 + 2 + 4];
    const int nq = 4 * 2 + 3 + 2 + 4;
    long ntot = 0, off = 0;
    int w, i, b, q;
    char name[128];

    const char *tag = use_reus ? "REUS" : "umbrella";
    printf("\n== %d. %s in K (L = %d, tobs = %d, eps = %g, s = %g, %d windows) ==\n",
           use_reus ? 6 : 5, use_reus ? "replica-exchange umbrella sampling" : "umbrella sampling",
           L, tobs, eps, s_run, nwin);

    /* exact reference */
    tm_init(&tm, L, 1.0, eps, 1, 0);
    tm_pk(&tm, 0.0, tobs, lp);
    ec.logp = lp;
    ec.Kmax = Kmax;
    ec.Klo = -0.5;
    ec.Khi = 200.5; /* same histogram range for exact and MBAR */
    ec.nbins = 40;
    g_Klo = ec.Klo;
    g_Khi = ec.Khi;
    g_nb = ec.nbins;
    sp.prob_below = exact_below;
    sp.hist = exact_hist;
    sp.ctx = &ec;
    sp.Klo = ec.Klo;
    sp.Khi = ec.Khi;
    sp.nbins = ec.nbins;
    sp.smooth_res = 20.0; /* 40 bins -> moving average over 3 bins */
    {
        double Kauto = NAN;
        double s_auto = equal_weight_auto(&sp, sv[1], 0.0, 0.06, &Kauto);
        Kc[0] = Kauto; /* the exact valley, also used for the fixed-cut test */
        Kc[1] = 60.5;
        Kc[2] = 120.5;
        q = 0;
        for (i = 0; i < 4; i++) {
            TMFinite f;
            tm_finite(&tm, sv[i], 0.0, tobs, &f);
            ex[q++] = f.meanK;
            ex[q++] = f.varK;
        }
        for (i = 0; i < 3; i++) ex[q++] = exact_below(sv[1], Kc[i], &ec);
        ex[q++] = equal_weight_fixed(&sp, Kc[0], 0.0, 0.06);
        ex[q++] = s_auto;
        s_split = s_auto;
        exact_split(&tm, lp, Kmax, tobs, s_split, Kc[0], ex + q);
        printf("     exact: valley cut K = %.1f, s*(auto) = %.6f\n", Kauto, s_auto);
        printf("     exact at s*: <K> %.3f | %.3f, <E> %.3f | %.3f (K < cut | K >= cut)\n",
               ex[q], ex[q + 1], ex[q + 2], ex[q + 3]);
    }

    if (!use_reus) {
        /* umbrella windows, each started from the empty trajectory */
        for (w = 0; w < nwin; w++) {
            TPSParams p;
            tps_default_params(&p);
            p.L = L;
            p.tobs = tobs;
            p.eps = eps;
            p.s = s_run;
            p.seed = 7000 + w;
            p.n_equil = 2000;
            p.n_prod = 300000;
            p.sample_every = 5;
            p.init = "empty";
            p.umb_kmin = w * delta;
            p.umb_kmax = w * delta + 2.0 * delta;
            if (tps_run(&p, &res[w]) != 0) {
                n_fail++;
                printf("[FAIL] umbrella window %d did not run\n", w);
                return;
            }
            ntot += res[w].nsamp;
        }
        mbar_alloc(&m, nwin, ntot);
        for (w = 0; w < nwin; w++) {
            m.ws[w] = s_run;
            m.wg[w] = 0.0;
            m.Kmin[w] = (double)res[w].Kmin;
            m.Kmax[w] = (double)res[w].Kmax;
            m.N[w] = res[w].nsamp;
            memcpy(m.K + off, res[w].K, sizeof(double) * (size_t)res[w].nsamp);
            memcpy(m.E + off, res[w].E, sizeof(double) * (size_t)res[w].nsamp);
            off += res[w].nsamp;
            tps_result_free(&res[w]);
        }
    } else {
        REUSParams p;
        REUSResult rr;
        long trips = 0;
        double ksmax = 0.0, accmin = 1.0;
        reus_default_params(&p);
        p.L = L;
        p.tobs = tobs;
        p.eps = eps;
        p.s = s_run;
        p.g = 0.0;
        p.seed = 9100;
        p.nwin = nwin;
        p.k_lo = 0.0;
        p.k_hi = kspan;
        p.n_equil = 2000;
        p.n_prod = 300000;
        p.swap_every = 5;
        p.sample_every = 5;
        p.init = "empty";
        if (reus_run(&p, &rr) != 0) {
            n_fail++;
            printf("[FAIL] REUS did not run\n");
            return;
        }
        mbar_alloc(&m, nwin, (long)nwin * rr.nsamp);
        for (w = 0; w < nwin; w++) {
            long n1, n2;
            double ks;
            m.ws[w] = s_run;
            m.wg[w] = 0.0;
            m.Kmin[w] = (double)rr.Kmin[w];
            m.Kmax[w] = (double)rr.Kmax[w];
            m.N[w] = rr.nsamp;
            memcpy(m.K + off, rr.K + (size_t)w * rr.nsamp, sizeof(double) * (size_t)rr.nsamp);
            memcpy(m.E + off, rr.E + (size_t)w * rr.nsamp, sizeof(double) * (size_t)rr.nsamp);
            off += rr.nsamp;
            trips += rr.round_trips[w];
            if (w + 1 < nwin) {
                double a = (double)rr.swap_acc[w] / rr.swap_try[w];
                ks = reus_overlap_ks(&rr, w, &n1, &n2);
                if (ks > ksmax) ksmax = ks;
                if (a < accmin) accmin = a;
            }
        }
        sprintf(name, "REUS round trips > 20 (got %ld)", trips);
        check_true(name, trips > 20);
        sprintf(name, "REUS min swap acceptance > 0.01 (got %.3f)", accmin);
        check_true(name, accmin > 0.01);
        sprintf(name, "REUS max overlap KS distance < 0.1 (got %.3f)", ksmax);
        check_true(name, ksmax < 0.1);
        reus_result_free(&rr);
    }
    mbar_guess(&m);
    if (mbar_solve(&m, 1e-10, 500) != 0) {
        n_fail++;
        printf("[FAIL] MBAR over %s windows did not converge\n", tag);
    }
    umbrella_quantities(&m, sv, 4, Kc, 3, &sp, s_split, Kc[0], est);

    memset(acc1, 0, sizeof(acc1));
    memset(acc2, 0, sizeof(acc2));
    memset(cnt, 0, sizeof(cnt));
    rng_seed(&r, 4242ULL);
    for (b = 0; b < nboot; b++) {
        MBAR mb;
        double v[4 * 2 + 3 + 2 + 4];
        mbar_resample(&m, &mb, block, &r);
        mbar_solve(&mb, 1e-10, 500);
        umbrella_quantities(&mb, sv, 4, Kc, 3, &sp, s_split, Kc[0], v);
        for (q = 0; q < nq; q++)
            if (!isnan(v[q])) {
                acc1[q] += v[q];
                acc2[q] += v[q] * v[q];
                cnt[q]++;
            }
        mbar_free(&mb);
    }
    q = 0;
    for (i = 0; i < 4; i++) {
        double e1, e2;
        int c1 = cnt[q], c2 = cnt[q + 1];
        e1 = sqrt((acc2[q] / c1 - (acc1[q] / c1) * (acc1[q] / c1)) * c1 / (c1 - 1.0));
        e2 = sqrt((acc2[q + 1] / c2 - (acc1[q + 1] / c2) * (acc1[q + 1] / c2)) * c2 / (c2 - 1.0));
        sprintf(name, "%s <K>  at s=%g", tag, sv[i]);
        zcheck(name, est[q], e1, ex[q]);
        sprintf(name, "%s VarK at s=%g", tag, sv[i]);
        zcheck(name, est[q + 1], e2, ex[q + 1]);
        q += 2;
    }
    for (i = 0; i < 3; i++, q++) {
        double e = sqrt((acc2[q] / cnt[q] - (acc1[q] / cnt[q]) * (acc1[q] / cnt[q])) * cnt[q] /
                        (cnt[q] - 1.0));
        sprintf(name, "%s P(K < %.1f) at s=%g", tag, Kc[i], sv[1]);
        zcheck(name, est[q], e, ex[q]);
    }
    for (i = 0; i < 2; i++, q++) {
        double e = (cnt[q] > 1) ? sqrt((acc2[q] / cnt[q] - (acc1[q] / cnt[q]) * (acc1[q] / cnt[q])) *
                                       cnt[q] / (cnt[q] - 1.0))
                                : NAN;
        sprintf(name, i == 0 ? "%s equal-weight s* (fixed cut at exact valley)"
                             : "%s equal-weight s* (automatic valley cut)",
                tag);
        zcheck(name, est[q], e, ex[q]);
    }
    for (i = 0; i < 4; i++, q++) {
        static const char *lab[4] = {"<K> inactive peak", "<K> active peak", "<E> inactive peak",
                                     "<E> active peak"};
        double e = sqrt((acc2[q] / cnt[q] - (acc1[q] / cnt[q]) * (acc1[q] / cnt[q])) * cnt[q] /
                        (cnt[q] - 1.0));
        sprintf(name, "%s %s at exact s*, cut", tag, lab[i]);
        zcheck(name, est[q], e, ex[q]);
    }
    tm_free(&tm);
    mbar_free(&m);
    free(lp);
}

int main(void)
{
    printf("statistical validation (fixed seeds); pass if |z| <= %.0f\n", ZMAX);
    validate_kernel();
    validate_tps();
    validate_clone();
    validate_mbar();
    validate_umbrella(0);
    validate_umbrella(1);
    printf("\n== summary: %d passed, %d failed ==\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
