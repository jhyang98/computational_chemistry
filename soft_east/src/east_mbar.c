/*
 * east_mbar: combine TPS sample files (east_tps) with MBAR. Files may come
 * from different s, g and from umbrella windows in K.
 *
 * usage: ./bin/east_mbar [input_file]    (default input file: system.inp)
 *
 * keywords
 *   window      <file>     one line per TPS sample file (repeat)
 *   target_s    <list>     s values to evaluate (no umbrella window)
 *   target_g    <list>     g values to evaluate
 *   skip        <int>      samples dropped from the start of every file (0)
 *   output      <file>     results (default: stdout)
 *   kcut        <value|auto>
 *                          optional equal-weight point at each target g:
 *                          value: s where P(k < kcut) = 1/2, k = K/(L tobs)
 *                          auto : the cut is the valley between the two peaks of
 *                                 the reweighted P(k), found self-consistently
 *                                 at the equal-weight point (see peaks.h)
 *                          also reports the averages of k and e within each
 *                          peak (k < kcut and k >= kcut) at that point
 *   hist_bins   <int>      bins for histograms and the automatic cut (60)
 *   smooth_res  <value>    smoothing resolution m for the automatic cut (30):
 *                          moving average over the neighbouring
 *                          round(hist_bins / (2 m)) bins on each side, i.e. a
 *                          window of about 1/m of the histogram range
 *   hist_s      <list>     optional: s values for reweighted P(k)
 *   hist_output <file>     histogram file (default: mbar_hist.dat)
 *   n_boot      <int>      moving-block bootstrap replicas for error bars (0)
 *   block       <int>      bootstrap block length in samples (10)
 *   seed        <int>      bootstrap random seed (0: clock)
 *   tol         <value>    MBAR convergence, max |grad_k|/N_k (1e-10)
 */
#include "input.h"
#include "mbar.h"
#include "peaks.h"
#include "rng.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char path[512];
    int L, tobs;
    double T, eps, s, g;
    double Kmin, Kmax; /* umbrella window, +-inf if none */
    long n;
    double *K, *E;
} WinFile;

static int read_window(const char *path, long skip, WinFile *w)
{
    FILE *f = fopen(path, "r");
    char line[1024];
    long cap = 1024, cnt = 0, seen = 0;
    int have[6] = {0, 0, 0, 0, 0, 0};
    memset(w, 0, sizeof(*w));
    strncpy(w->path, path, sizeof(w->path) - 1);
    w->T = 1.0;
    w->Kmin = -INFINITY;
    w->Kmax = INFINITY;
    if (!f) {
        fprintf(stderr, "error: cannot open window file '%s'\n", path);
        return -1;
    }
    w->K = (double *)malloc(sizeof(double) * (size_t)cap);
    w->E = (double *)malloc(sizeof(double) * (size_t)cap);
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            char key[64];
            double val;
            if (sscanf(line, "# param %63s %lf", key, &val) == 2) {
                if (!strcmp(key, "L")) { w->L = (int)val; have[0] = 1; }
                else if (!strcmp(key, "tobs")) { w->tobs = (int)val; have[1] = 1; }
                else if (!strcmp(key, "T")) { w->T = val; have[2] = 1; }
                else if (!strcmp(key, "eps")) { w->eps = val; have[3] = 1; }
                else if (!strcmp(key, "s")) { w->s = val; have[4] = 1; }
                else if (!strcmp(key, "g")) { w->g = val; have[5] = 1; }
                else if (!strcmp(key, "umb_Kmin")) w->Kmin = val;
                else if (!strcmp(key, "umb_Kmax")) w->Kmax = val;
            }
            continue;
        }
        {
            long mv;
            double K, E;
            if (sscanf(line, "%ld %lf %lf", &mv, &K, &E) != 3) continue;
            seen++;
            if (seen <= skip) continue;
            if (cnt == cap) {
                cap *= 2;
                w->K = (double *)realloc(w->K, sizeof(double) * (size_t)cap);
                w->E = (double *)realloc(w->E, sizeof(double) * (size_t)cap);
            }
            w->K[cnt] = K;
            w->E[cnt] = E;
            cnt++;
        }
    }
    fclose(f);
    w->n = cnt;
    if (!have[0] || !have[1] || !have[3] || !have[4] || !have[5]) {
        fprintf(stderr, "error: '%s' lacks '# param' lines for L, tobs, eps, s, g\n", path);
        return -1;
    }
    if (cnt < 2) {
        fprintf(stderr, "error: '%s' has fewer than 2 samples after skip\n", path);
        return -1;
    }
    return 0;
}

/* Group the files into distinct states (s, g, window) and fill an MBAR object. */
static int build_mbar(WinFile *wf, int nf, MBAR *m)
{
    int nstate = 0, i, j;
    int *state_of = (int *)malloc(sizeof(int) * (size_t)nf);
    int *first = (int *)malloc(sizeof(int) * (size_t)nf);
    long *Ns = (long *)calloc((size_t)nf, sizeof(long)), ntot = 0, *off;
    for (i = 0; i < nf; i++) {
        for (j = 0; j < nstate; j++) {
            const WinFile *a = &wf[first[j]];
            if (a->s == wf[i].s && a->g == wf[i].g && a->Kmin == wf[i].Kmin &&
                a->Kmax == wf[i].Kmax)
                break;
        }
        if (j == nstate) first[nstate++] = i;
        state_of[i] = j;
        Ns[j] += wf[i].n;
        ntot += wf[i].n;
    }
    if (mbar_alloc(m, nstate, ntot) != 0) {
        free(state_of);
        free(first);
        free(Ns);
        return -1;
    }
    off = (long *)calloc((size_t)nstate, sizeof(long));
    for (j = 1; j < nstate; j++) off[j] = off[j - 1] + Ns[j - 1];
    for (j = 0; j < nstate; j++) {
        m->ws[j] = wf[first[j]].s;
        m->wg[j] = wf[first[j]].g;
        m->Kmin[j] = wf[first[j]].Kmin;
        m->Kmax[j] = wf[first[j]].Kmax;
        m->N[j] = Ns[j];
    }
    for (i = 0; i < nf; i++) {
        int st = state_of[i];
        memcpy(m->K + off[st], wf[i].K, sizeof(double) * (size_t)wf[i].n);
        memcpy(m->E + off[st], wf[i].E, sizeof(double) * (size_t)wf[i].n);
        off[st] += wf[i].n;
    }
    free(off);
    free(state_of);
    free(first);
    free(Ns);
    return 0;
}

/* ---- equal-weight point through the shared peaks.h routines ---- */

typedef struct {
    const MBAR *m;
    double g;
    double Klo, Khi;
    int nbins;
} MbarCtx;

static double cb_prob_below(double s, double Kcut, void *ctx)
{
    const MbarCtx *c = (const MbarCtx *)ctx;
    return mbar_prob_below(c->m, s, c->g, Kcut);
}

static void cb_hist(double s, double *h, void *ctx)
{
    const MbarCtx *c = (const MbarCtx *)ctx;
    mbar_histogram(c->m, s, c->g, c->Klo, c->Khi, c->nbins, h);
}

/* s* (and Kcut) at target g; auto_cut = 1: valley, else fixed Kcut */
static double equal_weight(const MBAR *m, double g, int auto_cut, double Kcut_fixed,
                           const double *ts, int nts, double Klo, double Khi, int nbins,
                           double smooth_res, double *Kcut_out)
{
    MbarCtx c;
    SplitProblem sp;
    double smin = ts[0], smax = ts[0];
    int i;
    for (i = 1; i < nts; i++) {
        if (ts[i] < smin) smin = ts[i];
        if (ts[i] > smax) smax = ts[i];
    }
    c.m = m;
    c.g = g;
    c.Klo = Klo;
    c.Khi = Khi;
    c.nbins = nbins;
    sp.prob_below = cb_prob_below;
    sp.hist = cb_hist;
    sp.ctx = &c;
    sp.Klo = Klo;
    sp.Khi = Khi;
    sp.nbins = nbins;
    sp.smooth_res = smooth_res;
    if (!auto_cut) {
        *Kcut_out = Kcut_fixed;
        return equal_weight_fixed(&sp, Kcut_fixed, smin, smax);
    }
    {
        /* start from the chi_KK maximum on the target grid */
        double best = -1.0, s0 = ts[0];
        for (i = 0; i < nts; i++) {
            MBARObs o;
            mbar_observe(m, ts[i], g, &o);
            if (o.varK > best) {
                best = o.varK;
                s0 = ts[i];
            }
        }
        return equal_weight_auto(&sp, s0, smin, smax, Kcut_out);
    }
}

typedef struct {
    double k, e, ckk, cee, cke;
} Derived;

static void derive(const MBARObs *o, int L, int tobs, Derived *d)
{
    const double nk = (double)L * tobs, ne = (double)L * (tobs + 1);
    d->k = o->meanK / nk;
    d->e = o->meanE / ne;
    d->ckk = o->varK / nk;
    d->cee = o->varE / nk;
    d->cke = o->covKE / nk;
}

/* Averages within the two peaks at (s*, g), split at Kcut:
   v = {k_in, k_act, dk, e_in, e_act, de, sdk_in, sdk_act}. */
#define NPEAK 8
#define NPEAK_ERR 6
static void peak_values(const MBAR *m, double sstar, double g, double Kcut, int L, int tobs,
                        double *v)
{
    MBARSplit sp;
    const double nk = (double)L * tobs, ne = (double)L * (tobs + 1);
    int q;
    if (isnan(sstar) || isnan(Kcut)) {
        for (q = 0; q < NPEAK; q++) v[q] = NAN;
        return;
    }
    mbar_observe_split(m, sstar, g, Kcut, &sp);
    v[0] = sp.meanK_lo / nk;
    v[1] = sp.meanK_hi / nk;
    v[2] = v[1] - v[0];
    v[3] = sp.meanE_lo / ne;
    v[4] = sp.meanE_hi / ne;
    v[5] = v[4] - v[3];
    v[6] = sqrt(sp.varK_lo) / nk;
    v[7] = sqrt(sp.varK_hi) / nk;
}

static double std_from(double s1, double s2, int n)
{
    double mu = s1 / n, var = s2 / n - mu * mu;
    return (var > 0.0 && n > 1) ? sqrt(var * n / (n - 1.0)) : 0.0;
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "system.inp";
    const char *out = NULL, *hist_out = "mbar_hist.dat", *v;
    Input in;
    int nf, i, j, nts = 0, ntg = 0, nhs = 0, hist_bins = 60, n_boot = 0, block = 10;
    double smooth_res = 30.0;
    int use_kcut = 0, auto_cut = 0;
    long skip = 0;
    double *ts = NULL, *tg = NULL, *hs = NULL, kcut = NAN, tol = 1e-10;
    double Klo, Khi;
    uint64_t seed = 0;
    WinFile *wf;
    MBAR m;
    FILE *f = stdout;
    int L, tobs, npt;
    double *acc1 = NULL, *acc2 = NULL, *sacc1 = NULL, *sacc2 = NULL;
    double *pacc1 = NULL, *pacc2 = NULL;
    int *scount = NULL, *pcount = NULL;

    if (argc > 1 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        printf("usage: ./bin/east_mbar [input_file]   (default: system.inp)\n"
               "keywords: window target_s target_g skip output kcut hist_bins smooth_res\n"
               "          hist_s hist_output n_boot block seed tol\n");
        return 0;
    }
    if (input_read(&in, path) != 0) return 1;
    if (input_str(&in, "smooth")) {
        fprintf(stderr, "error: keyword 'smooth' was replaced by 'smooth_res' (resolution m; the "
                        "moving average spans about 1/m of the histogram range)\n");
        return 1;
    }
    nf = input_count(&in, "window");
    if (nf < 1) {
        fprintf(stderr, "error: no 'window' entries in %s\n", path);
        return 1;
    }
    if (input_list(&in, "target_s", 1, &ts, &nts) < 0 ||
        input_list(&in, "target_g", 1, &tg, &ntg) < 0 || input_long(&in, "skip", 0, &skip) < 0 ||
        input_int(&in, "hist_bins", 0, &hist_bins) < 0 ||
        input_double(&in, "smooth_res", 0, &smooth_res) < 0 || input_int(&in, "n_boot", 0, &n_boot) < 0 ||
        input_int(&in, "block", 0, &block) < 0 || input_u64(&in, "seed", 0, &seed) < 0 ||
        input_double(&in, "tol", 0, &tol) < 0) {
        input_free(&in);
        return 1;
    }
    if ((v = input_str(&in, "kcut"))) {
        use_kcut = 1;
        if (!strcmp(v, "auto")) auto_cut = 1;
        else if (input_double(&in, "kcut", 1, &kcut) < 0) return 1;
    }
    if (input_str(&in, "hist_s") && input_list(&in, "hist_s", 0, &hs, &nhs) < 0) {
        input_free(&in);
        return 1;
    }
    if ((v = input_str(&in, "output"))) out = v;
    if ((v = input_str(&in, "hist_output"))) hist_out = v;

    wf = (WinFile *)calloc((size_t)nf, sizeof(WinFile));
    for (i = 0; i < nf; i++)
        if (read_window(input_nth(&in, "window", i), skip, &wf[i]) != 0) return 1;
    input_warn_unused(&in);
    L = wf[0].L;
    tobs = wf[0].tobs;
    for (i = 1; i < nf; i++) {
        if (wf[i].L != L || wf[i].tobs != tobs || wf[i].eps != wf[0].eps || wf[i].T != wf[0].T) {
            fprintf(stderr, "error: '%s' has different L, tobs, eps or T from '%s'\n", wf[i].path,
                    wf[0].path);
            return 1;
        }
    }
    if (build_mbar(wf, nf, &m) != 0) {
        fprintf(stderr, "error: cannot set up MBAR\n");
        return 1;
    }
    mbar_guess(&m);
    if (mbar_solve(&m, tol, 500) != 0)
        fprintf(stderr, "warning: MBAR did not converge (max |grad|/N = %.2e)\n", m.grad_max);

    /* histogram range covering all samples, bins centred on integers at the ends */
    Klo = Khi = m.K[0];
    for (i = 1; i < m.ntot; i++) {
        if (m.K[i] < Klo) Klo = m.K[i];
        if (m.K[i] > Khi) Khi = m.K[i];
    }
    Klo -= 0.5;
    Khi += 0.5;

    if (out) {
        f = fopen(out, "w");
        if (!f) {
            fprintf(stderr, "error: cannot write %s\n", out);
            return 1;
        }
    }
    fprintf(f, "# soft_east MBAR over %d files, %d states, %ld samples\n", nf, m.nwin, m.ntot);
    fprintf(f, "# L %d  tobs %d  T %.12g  eps %.12g  skip %ld\n", L, tobs, wf[0].T, wf[0].eps,
            skip);
    fprintf(f, "# Newton iterations %d, max |grad_k|/N_k = %.2e\n", m.iters, m.grad_max);
    fprintf(f, "# state: s g Kmin Kmax N f <K>_raw <K>_mbar (both within the state's window)\n");
    {
        long off = 0;
        for (j = 0; j < m.nwin; j++) {
            MBARObs o;
            double raw = 0.0;
            long n;
            for (n = 0; n < m.N[j]; n++) raw += m.K[off + n];
            raw /= (double)m.N[j];
            off += m.N[j];
            mbar_observe_window(&m, m.ws[j], m.wg[j], m.Kmin[j], m.Kmax[j], &o);
            fprintf(f, "#   %12.8g %10.6g %8.0f %8.0f %8ld %16.8e %14.6e %14.6e\n", m.ws[j], m.wg[j],
                    isfinite(m.Kmin[j]) ? m.Kmin[j] : -1.0, isfinite(m.Kmax[j]) ? m.Kmax[j] : -1.0,
                    m.N[j], m.f[j], raw, o.meanK);
        }
    }
    fprintf(f, "# (Kmin, Kmax = -1: no umbrella window)\n");
    fprintf(f, "# k = <K>/(L tobs), e = <E>/(L (tobs+1)), chi_XY = Cov(X,Y)/(L tobs)\n");
    fprintf(f, "# ESS = Kish effective sample size of the reweighted samples (small = unreliable)\n");

    npt = nts * ntg;
    if (n_boot > 0) {
        Rng r;
        int b;
        acc1 = (double *)calloc((size_t)npt * 5, sizeof(double));
        acc2 = (double *)calloc((size_t)npt * 5, sizeof(double));
        sacc1 = (double *)calloc((size_t)ntg, sizeof(double));
        sacc2 = (double *)calloc((size_t)ntg, sizeof(double));
        scount = (int *)calloc((size_t)ntg, sizeof(int));
        pacc1 = (double *)calloc((size_t)ntg * NPEAK_ERR, sizeof(double));
        pacc2 = (double *)calloc((size_t)ntg * NPEAK_ERR, sizeof(double));
        pcount = (int *)calloc((size_t)ntg, sizeof(int));
        rng_seed_auto(&r, seed);
        for (b = 0; b < n_boot; b++) {
            MBAR mb;
            if (mbar_resample(&m, &mb, block, &r) != 0) {
                fprintf(stderr, "error: out of memory in bootstrap\n");
                return 1;
            }
            mbar_solve(&mb, tol, 500);
            for (j = 0; j < ntg; j++) {
                int k;
                for (k = 0; k < nts; k++) {
                    MBARObs o;
                    Derived d;
                    int p = j * nts + k, q;
                    double vals[5];
                    mbar_observe(&mb, ts[k], tg[j], &o);
                    derive(&o, L, tobs, &d);
                    vals[0] = d.k;
                    vals[1] = d.e;
                    vals[2] = d.ckk;
                    vals[3] = d.cee;
                    vals[4] = d.cke;
                    for (q = 0; q < 5; q++) {
                        acc1[p * 5 + q] += vals[q];
                        acc2[p * 5 + q] += vals[q] * vals[q];
                    }
                }
                if (use_kcut) {
                    double Kc;
                    double sstar = equal_weight(&mb, tg[j], auto_cut, kcut * L * tobs, ts, nts,
                                                Klo, Khi, hist_bins, smooth_res, &Kc);
                    if (!isnan(sstar)) {
                        double pv[NPEAK];
                        int q, finite = 1;
                        sacc1[j] += sstar;
                        sacc2[j] += sstar * sstar;
                        scount[j]++;
                        peak_values(&mb, sstar, tg[j], Kc, L, tobs, pv);
                        for (q = 0; q < NPEAK_ERR; q++)
                            if (!isfinite(pv[q])) finite = 0;
                        if (finite) {
                            for (q = 0; q < NPEAK_ERR; q++) {
                                pacc1[j * NPEAK_ERR + q] += pv[q];
                                pacc2[j * NPEAK_ERR + q] += pv[q] * pv[q];
                            }
                            pcount[j]++;
                        }
                    }
                }
            }
            mbar_free(&mb);
        }
        fprintf(f, "# errors: std over %d moving-block bootstrap replicas (block %d samples)\n",
                n_boot, block);
    }

    fprintf(f, "# %10s %12s %16s %16s %16s %16s %16s %12s %12s", "g", "s", "k", "e", "chi_KK",
            "chi_EE", "chi_KE", "rho", "ESS");
    if (n_boot > 0)
        fprintf(f, " %12s %12s %12s %12s %12s", "err_k", "err_e", "err_chiKK", "err_chiEE",
                "err_chiKE");
    fprintf(f, "\n");
    for (j = 0; j < ntg; j++) {
        int k;
        for (k = 0; k < nts; k++) {
            MBARObs o;
            Derived d;
            int p = j * nts + k, q;
            mbar_observe(&m, ts[k], tg[j], &o);
            derive(&o, L, tobs, &d);
            fprintf(f, "%12.6g %12.8g %16.10e %16.10e %16.10e %16.10e %16.10e %12.8f %12.1f",
                    tg[j], ts[k], d.k, d.e, d.ckk, d.cee, d.cke, d.cke / sqrt(d.ckk * d.cee),
                    o.ess);
            if (n_boot > 0)
                for (q = 0; q < 5; q++)
                    fprintf(f, " %12.4e", std_from(acc1[p * 5 + q], acc2[p * 5 + q], n_boot));
            fprintf(f, "\n");
        }
        if (ntg > 1) fprintf(f, "\n");
    }
    if (use_kcut) {
        if (auto_cut)
            fprintf(f, "# equal-weight point: cut at the valley of P(k), P(k < cut) = 1/2\n"
                       "#   %d bins of width dk = %.6g; smoothing resolution m = %g: moving average "
                       "over %d bins (dk window = %.6g)\n",
                    hist_bins, (Khi - Klo) / hist_bins / ((double)L * tobs), smooth_res,
                    2 * smooth_halfwidth(hist_bins, smooth_res) + 1,
                    (2 * smooth_halfwidth(hist_bins, smooth_res) + 1) * (Khi - Klo) / hist_bins /
                        ((double)L * tobs));
        else
            fprintf(f, "# equal-weight point: P(k < %.6g) = 1/2\n", kcut);
        {
            double *sst = (double *)malloc(sizeof(double) * (size_t)ntg);
            double *kcs = (double *)malloc(sizeof(double) * (size_t)ntg);
            for (j = 0; j < ntg; j++) {
                double Kc = NAN;
                double sstar = equal_weight(&m, tg[j], auto_cut, kcut * L * tobs, ts, nts, Klo,
                                            Khi, hist_bins, smooth_res, &Kc);
                sst[j] = sstar;
                kcs[j] = Kc;
                fprintf(f, "# s_equal g %.6g s* %.10f kcut %.8f", tg[j], sstar,
                        Kc / ((double)L * tobs));
                if (n_boot > 0)
                    fprintf(f, " err %.3e (%d/%d replicas)",
                            scount[j] > 1 ? std_from(sacc1[j], sacc2[j], scount[j]) : NAN,
                            scount[j], n_boot);
                if (isnan(sstar))
                    fprintf(f, "  (%s)", auto_cut ? "no two-peak valley or not bracketed by target_s"
                                                  : "not bracketed by target_s");
                fprintf(f, "\n");
            }
            fprintf(f, "# peaks: averages within each peak at the equal-weight point s*\n"
                       "#   inactive: k < kcut, active: k >= kcut; dk = k_act - k_in, "
                       "de = e_act - e_in\n"
                       "#   sdk_in, sdk_act: standard deviation of k within the peak\n");
            if (n_boot > 0)
                fprintf(f, "#   err_*: std over the bootstrap replicas, each at its own s* and "
                           "kcut; nrep = replicas used\n");
            fprintf(f, "#   columns: g s* kcut k_in k_act dk e_in e_act de sdk_in sdk_act");
            if (n_boot > 0)
                fprintf(f, " err_k_in err_k_act err_dk err_e_in err_e_act err_de nrep");
            fprintf(f, "\n");
            for (j = 0; j < ntg; j++) {
                double pv[NPEAK];
                int q;
                peak_values(&m, sst[j], tg[j], kcs[j], L, tobs, pv);
                fprintf(f, "# peaks %.6g %.10f %.8f", tg[j], sst[j], kcs[j] / ((double)L * tobs));
                for (q = 0; q < NPEAK; q++) fprintf(f, " %.8e", pv[q]);
                if (n_boot > 0) {
                    for (q = 0; q < NPEAK_ERR; q++)
                        fprintf(f, " %.3e",
                                pcount[j] > 1 ? std_from(pacc1[j * NPEAK_ERR + q],
                                                         pacc2[j * NPEAK_ERR + q], pcount[j])
                                              : NAN);
                    fprintf(f, " %d", pcount[j]);
                }
                fprintf(f, "\n");
            }
            free(sst);
            free(kcs);
        }
    }
    if (f != stdout) fclose(f);

    if (nhs > 0) {
        FILE *h = fopen(hist_out, "w");
        double *hh = (double *)malloc(sizeof(double) * (size_t)hist_bins);
        const double nk = (double)L * tobs, width = (Khi - Klo) / hist_bins;
        int b;
        if (!h) {
            fprintf(stderr, "error: cannot write %s\n", hist_out);
            return 1;
        }
        fprintf(h, "# reweighted activity distribution P(k), k = K/(L tobs)\n");
        fprintf(h, "# columns: g s k P(k) (density in k); blocks separated by blank lines\n");
        for (j = 0; j < ntg; j++) {
            int k;
            for (k = 0; k < nhs; k++) {
                mbar_histogram(&m, hs[k], tg[j], Klo, Khi, hist_bins, hh);
                for (b = 0; b < hist_bins; b++)
                    fprintf(h, "%12.6g %12.8g %14.8e %14.8e\n", tg[j], hs[k],
                            (Klo + (b + 0.5) * width) / nk, hh[b] * nk);
                fprintf(h, "\n\n");
            }
        }
        fclose(h);
        free(hh);
    }

    for (i = 0; i < nf; i++) {
        free(wf[i].K);
        free(wf[i].E);
    }
    free(wf);
    free(acc1);
    free(acc2);
    free(sacc1);
    free(sacc2);
    free(scount);
    free(pacc1);
    free(pacc2);
    free(pcount);
    mbar_free(&m);
    free(ts);
    free(tg);
    free(hs);
    input_free(&in);
    return 0;
}
