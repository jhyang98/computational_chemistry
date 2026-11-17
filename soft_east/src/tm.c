#include "tm.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NWORK 16

static double vdot(long n, const double *a, const double *b)
{
    long i;
    double s = 0.0;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

int tm_init(TMSolver *t, int L, double T, double eps, int use_translation, int exclude_empty)
{
    long n, j;
    memset(t, 0, sizeof(*t));
    if (model_init(&t->model, L, T, eps) != 0) return -1;
    if (basis_build(&t->basis, L, use_translation, exclude_empty) != 0) return -1;
    if (tgen_build(&t->tg, &t->model, &t->basis) != 0) {
        basis_free(&t->basis);
        return -1;
    }
    n = t->basis.dim;
    t->psi0 = (double *)malloc(sizeof(double) * (size_t)n);
    t->work = (double *)malloc(sizeof(double) * (size_t)n * NWORK);
    if (!t->psi0 || !t->work) {
        tm_free(t);
        return -1;
    }
    for (j = 0; j < n; j++) {
        double nn = t->tg.occ[j];
        double c = t->model.c;
        double logpi = nn * log(c) + (L - nn) * log(1.0 - c);
        t->psi0[j] = sqrt((double)t->basis.period[j]) * exp(0.5 * logpi);
    }
    lanczos_default_opts(&t->lopt);
    t->lopt.tol = 1e-12;
    return 0;
}

void tm_free(TMSolver *t)
{
    tgen_free(&t->tg);
    basis_free(&t->basis);
    free(t->psi0);
    free(t->work);
    t->psi0 = t->work = NULL;
}

/* y = M_s x (single attempt) */
static void apply_M(const TMSolver *t, double s, const double *x, double *y)
{
    const long n = t->basis.dim;
    const double invL = 1.0 / t->model.L;
    long i;
    tgen_apply(&t->tg, s, 0.0, x, y);
    for (i = 0; i < n; i++) y[i] = x[i] + invL * y[i];
}

int tm_finite(TMSolver *t, double s, double g, int tobs, TMFinite *out)
{
    const long n = t->basis.dim;
    const int L = t->model.L;
    const double es = exp(-s), invL = 1.0 / L;
    double *u = t->work, *us = u + n, *ug = us + n, *uss = ug + n, *ugg = uss + n,
           *usg = ugg + n;
    double *Au = usg + n, *Aus = Au + n, *Aug = Aus + n, *Auss = Aug + n, *Augg = Auss + n,
           *Ausg = Augg + n;
    double *d = Ausg + n;
    double lnscale = 0.0, Z, Zs, Zg, Zss, Zgg, Zsg;
    long i;
    int step, a;

    for (i = 0; i < n; i++) d[i] = exp(-g * t->tg.occ[i]);

    /* initial slice: u = D psi0 and its g-derivatives */
    for (i = 0; i < n; i++) {
        double ni = t->tg.occ[i];
        u[i] = d[i] * t->psi0[i];
        us[i] = 0.0;
        ug[i] = -ni * u[i];
        uss[i] = 0.0;
        ugg[i] = ni * ni * u[i];
        usg[i] = 0.0;
    }

    for (step = 0; step < tobs; step++) {
        double scale = 0.0;
        for (a = 0; a < L; a++) {
            tgen_apply_A(&t->tg, u, Au);
            tgen_apply_A(&t->tg, us, Aus);
            tgen_apply_A(&t->tg, ug, Aug);
            tgen_apply_A(&t->tg, uss, Auss);
            tgen_apply_A(&t->tg, ugg, Augg);
            tgen_apply_A(&t->tg, usg, Ausg);
            for (i = 0; i < n; i++) {
                const double r = t->tg.esc[i];
                /* M v = v + (es A v - r v)/L ; M' v = -es A v/L ; M'' v = es A v/L */
                double nu = u[i] + invL * (es * Au[i] - r * u[i]);
                double nus = us[i] + invL * (es * Aus[i] - r * us[i]) - invL * es * Au[i];
                double nug = ug[i] + invL * (es * Aug[i] - r * ug[i]);
                double nuss = uss[i] + invL * (es * Auss[i] - r * uss[i]) -
                              2.0 * invL * es * Aus[i] + invL * es * Au[i];
                double nugg = ugg[i] + invL * (es * Augg[i] - r * ugg[i]);
                double nusg = usg[i] + invL * (es * Ausg[i] - r * usg[i]) - invL * es * Aug[i];
                u[i] = nu;
                us[i] = nus;
                ug[i] = nug;
                uss[i] = nuss;
                ugg[i] = nugg;
                usg[i] = nusg;
            }
        }
        /* energy weight of the new slice: D, D' = -n D, D'' = n^2 D */
        for (i = 0; i < n; i++) {
            const double ni = t->tg.occ[i], di = d[i];
            double nugg = di * (ugg[i] - 2.0 * ni * ug[i] + ni * ni * u[i]);
            double nug = di * (ug[i] - ni * u[i]);
            double nusg = di * (usg[i] - ni * us[i]);
            u[i] *= di;
            us[i] *= di;
            uss[i] *= di;
            ug[i] = nug;
            ugg[i] = nugg;
            usg[i] = nusg;
            if (fabs(u[i]) > scale) scale = fabs(u[i]);
        }
        if (scale > 0.0) {
            for (i = 0; i < n; i++) {
                u[i] /= scale;
                us[i] /= scale;
                ug[i] /= scale;
                uss[i] /= scale;
                ugg[i] /= scale;
                usg[i] /= scale;
            }
            lnscale += log(scale);
        }
    }
    Z = vdot(n, t->psi0, u);
    Zs = vdot(n, t->psi0, us) / Z;
    Zg = vdot(n, t->psi0, ug) / Z;
    Zss = vdot(n, t->psi0, uss) / Z;
    Zgg = vdot(n, t->psi0, ugg) / Z;
    Zsg = vdot(n, t->psi0, usg) / Z;
    out->lnZ = log(Z) + lnscale;
    out->meanK = -Zs;
    out->meanE = -Zg;
    out->varK = Zss - Zs * Zs;
    out->varE = Zgg - Zg * Zg;
    out->covKE = Zsg - Zs * Zg;
    return 0;
}

/* S x = D^{1/2} M^L D^{1/2} x */
static void tm_matvec(const double *x, double *y, void *ctx)
{
    TMSolver *t = (TMSolver *)ctx;
    const long n = t->basis.dim;
    double *a = t->work + 12 * n, *b = t->work + 13 * n;
    long i;
    int k;
    for (i = 0; i < n; i++) a[i] = exp(-0.5 * t->cur_g * t->tg.occ[i]) * x[i];
    for (k = 0; k < t->model.L; k++) {
        apply_M(t, t->cur_s, a, b);
        memcpy(a, b, sizeof(double) * (size_t)n);
    }
    for (i = 0; i < n; i++) y[i] = exp(-0.5 * t->cur_g * t->tg.occ[i]) * a[i];
}

int tm_asym(TMSolver *t, double s, double g, TMAsym *out)
{
    const long n = t->basis.dim;
    const int L = t->model.L;
    double evals[2], resid[2];
    double *evecs = (double *)malloc(sizeof(double) * (size_t)n * 2);
    double *ws = (double *)malloc(sizeof(double) * (size_t)n * (size_t)L);
    double *Aw = t->work;
    LanczosInfo info;
    double dlam = 0.0, en = 0.0;
    long i;
    int j;

    if (!evecs || !ws) {
        free(evecs);
        free(ws);
        return -1;
    }
    t->cur_s = s;
    t->cur_g = g;
    if (lanczos_top(n, tm_matvec, t, t->psi0, (n > 1) ? 2 : 1, evals, evecs, resid, &t->lopt,
                    &info) != 0) {
        free(evecs);
        free(ws);
        return -1;
    }
    out->theta = log(evals[0]);
    out->lam_ratio = (n > 1) ? evals[1] / evals[0] : 0.0;
    out->resid = resid[0];
    out->ok = info.converged;

    /* e_d = <v|n|v>/L */
    for (i = 0; i < n; i++) en += t->tg.occ[i] * evecs[i] * evecs[i];
    out->e = en / L;

    /* d lambda/ds = sum_j <M^j w| M'_s |M^{L-1-j} w>, w = D^{1/2} v, M' = -exp(-s) A/L */
    for (i = 0; i < n; i++) ws[i] = exp(-0.5 * g * t->tg.occ[i]) * evecs[i];
    for (j = 1; j < L; j++) apply_M(t, s, ws + (size_t)(j - 1) * n, ws + (size_t)j * n);
    for (j = 0; j < L; j++) {
        tgen_apply_A(&t->tg, ws + (size_t)(L - 1 - j) * n, Aw);
        dlam += vdot(n, ws + (size_t)j * n, Aw);
    }
    dlam *= -exp(-s) / L;
    out->k = -(dlam / evals[0]) / L;
    free(evecs);
    free(ws);
    return 0;
}

int tm_mcs_row(TMSolver *t, uint32_t x, double *prob)
{
    const long n = t->basis.dim;
    const int L = t->model.L;
    double *a = t->work, *b = t->work + n;
    long i;
    int k;
    if (t->basis.use_translation || t->basis.exclude_empty) return -1;
    memset(a, 0, sizeof(double) * (size_t)n);
    a[t->basis.index[x]] = 1.0 / t->psi0[t->basis.index[x]];
    for (k = 0; k < L; k++) {
        apply_M(t, 0.0, a, b);
        memcpy(a, b, sizeof(double) * (size_t)n);
    }
    for (i = 0; i < n; i++) prob[t->basis.rep[i]] = t->psi0[i] * a[i];
    return 0;
}

int tm_pk(TMSolver *t, double g, int tobs, double *logp)
{
    const long n = t->basis.dim;
    const int L = t->model.L;
    const long Kmax = (long)L * tobs;
    const double invL = 1.0 / L;
    double *u = (double *)calloc((size_t)(Kmax + 1) * (size_t)n, sizeof(double));
    double *Au = (double *)malloc(sizeof(double) * (size_t)n);
    double *d = (double *)malloc(sizeof(double) * (size_t)n);
    double lnscale = 0.0, norm;
    long K, i, top = 0; /* highest K that can be nonzero */
    int step, a;

    if (!u || !Au || !d) {
        free(u);
        free(Au);
        free(d);
        return -1;
    }
    for (i = 0; i < n; i++) {
        d[i] = exp(-g * t->tg.occ[i]);
        u[i] = d[i] * t->psi0[i]; /* K = 0, slice 0 weighted */
    }
    for (step = 0; step < tobs; step++) {
        double scale = 0.0;
        for (a = 0; a < L; a++) {
            /* u_K <- (I - R/L) u_K + (A/L) u_{K-1}; descending K updates in place */
            for (K = top + 1; K >= 0; K--) {
                double *uk = u + (size_t)K * n;
                if (K > Kmax) continue;
                for (i = 0; i < n; i++) uk[i] *= 1.0 - invL * t->tg.esc[i];
                if (K > 0) {
                    tgen_apply_A(&t->tg, u + (size_t)(K - 1) * n, Au);
                    for (i = 0; i < n; i++) uk[i] += invL * Au[i];
                }
            }
            if (top < Kmax) top++;
        }
        for (K = 0; K <= top; K++) {
            double *uk = u + (size_t)K * n;
            for (i = 0; i < n; i++) {
                uk[i] *= d[i];
                if (uk[i] > scale) scale = uk[i];
            }
        }
        if (scale > 0.0) {
            for (K = 0; K <= top; K++) {
                double *uk = u + (size_t)K * n;
                for (i = 0; i < n; i++) uk[i] /= scale;
            }
            lnscale += log(scale);
        }
    }
    norm = 0.0;
    for (K = 0; K <= Kmax; K++) {
        double p = 0.0;
        const double *uk = u + (size_t)K * n;
        for (i = 0; i < n; i++) p += t->psi0[i] * uk[i];
        logp[K] = (p > 0.0) ? log(p) : -INFINITY;
        norm += p;
    }
    for (K = 0; K <= Kmax; K++) logp[K] -= log(norm);
    (void)lnscale;
    free(u);
    free(Au);
    free(d);
    return 0;
}
