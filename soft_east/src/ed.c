#include "ed.h"
#include "cg.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void ed_matvec(const double *x, double *y, void *ctx)
{
    EDSolver *ed = (EDSolver *)ctx;
    tgen_apply(&ed->tg, ed->cur_s, ed->cur_g, x, y);
}

static double vdot(long n, const double *a, const double *b)
{
    long i;
    double s = 0.0;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

int ed_init(EDSolver *ed, int L, double T, double eps, int use_translation)
{
    long n;
    memset(ed, 0, sizeof(*ed));
    if (model_init(&ed->model, L, T, eps) != 0) return -1;
    /* strict East: the empty configuration is disconnected, drop it */
    if (basis_build(&ed->basis, L, use_translation, eps == 0.0) != 0) return -1;
    if (tgen_build(&ed->tg, &ed->model, &ed->basis) != 0) {
        basis_free(&ed->basis);
        return -1;
    }
    n = ed->basis.dim;
    lanczos_default_opts(&ed->lopt);
    ed->cg_tol = 1e-10;
    ed->cg_maxit = 20000;
    ed->fd_h = 1e-4;
    ed->v = (double *)malloc(sizeof(double) * (size_t)n);
    ed->evecs = (double *)malloc(sizeof(double) * (size_t)n * 2);
    ed->w1 = (double *)malloc(sizeof(double) * (size_t)n);
    ed->w2 = (double *)malloc(sizeof(double) * (size_t)n);
    ed->w3 = (double *)malloc(sizeof(double) * (size_t)n);
    ed->w4 = (double *)malloc(sizeof(double) * (size_t)n);
    if (!ed->v || !ed->evecs || !ed->w1 || !ed->w2 || !ed->w3 || !ed->w4) {
        ed_free(ed);
        return -1;
    }
    ed_reset_guess(ed);
    return 0;
}

void ed_free(EDSolver *ed)
{
    tgen_free(&ed->tg);
    basis_free(&ed->basis);
    free(ed->v);
    free(ed->evecs);
    free(ed->w1);
    free(ed->w2);
    free(ed->w3);
    free(ed->w4);
    ed->v = ed->evecs = ed->w1 = ed->w2 = ed->w3 = ed->w4 = NULL;
}

void ed_reset_guess(EDSolver *ed)
{
    const double c = ed->model.c;
    const int L = ed->model.L;
    long j;
    for (j = 0; j < ed->basis.dim; j++) {
        double n = ed->tg.occ[j];
        double logpi = n * log(c) + (L - n) * log(1.0 - c);
        ed->v[j] = sqrt((double)ed->basis.period[j]) * exp(0.5 * logpi);
    }
}

/* Top two eigenpairs plus Hellmann-Feynman first derivatives. */
static int solve_top(EDSolver *ed, double s, double g, EDPoint *out)
{
    const long n = ed->basis.dim;
    const int L = ed->model.L;
    double evals[2], resid[2];
    LanczosInfo info;
    double *v;
    double vmin, vmax;
    long i;

    ed->cur_s = s;
    ed->cur_g = g;
    if (lanczos_top(n, ed_matvec, ed, ed->v, 2, evals, ed->evecs, resid, &ed->lopt, &info) != 0)
        return -1;
    v = ed->evecs;
    memcpy(ed->v, v, sizeof(double) * (size_t)n);

    out->s = s;
    out->g = g;
    out->theta = evals[0];
    out->theta2 = (n > 1) ? evals[1] : NAN;
    out->resid = resid[0];
    out->resid2 = (n > 1) ? resid[1] : 0.0;
    out->nmatvec = info.nmatvec;
    out->ok = info.converged;

    /* Hellmann-Feynman: dH/ds = -exp(-s) A, dH/dg = -diag(n) */
    tgen_apply_A(&ed->tg, v, ed->w1);
    out->k = exp(-s) * vdot(n, v, ed->w1) / L;
    {
        double acc = 0.0;
        for (i = 0; i < n; i++) acc += ed->tg.occ[i] * v[i] * v[i];
        out->e = acc / L;
    }

    vmin = vmax = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] < vmin) vmin = v[i];
        if (v[i] > vmax) vmax = v[i];
    }
    out->min_ratio = (vmax != 0.0) ? vmin / vmax : 0.0;

    out->chi_kk = out->chi_ke = out->chi_ee = out->rho = NAN;
    out->cg_iters = 0;
    out->cg_relres = 0.0;
    return 0;
}

static void finish_rho(EDPoint *p)
{
    if (p->chi_kk > 0.0 && p->chi_ee > 0.0)
        p->rho = p->chi_ke / sqrt(p->chi_kk * p->chi_ee);
    else
        p->rho = NAN;
}

int ed_point(EDSolver *ed, double s, double g, ChiMode chi, EDPoint *out)
{
    const long n = ed->basis.dim;
    const int L = ed->model.L;
    long i;

    if (solve_top(ed, s, g, out) != 0) return -1;
    if (chi == CHI_NONE) return 0;

    if (chi == CHI_LR) {
        /* second-order perturbation theory:
           d2theta/dadb = <v|H_ab|v> + 2 <H_a v| R |H_b v>,
           R = reduced resolvent (theta0 - H)^{-1} on the complement of v. */
        const double *v = ed->v;
        double es = exp(-s);
        double *bs = ed->w1, *bg = ed->w2, *xs = ed->w3, *xg = ed->w4;
        double vAv, c1, c2;
        CGInfo ci1, ci2;

        tgen_apply_A(&ed->tg, v, bs);
        vAv = vdot(n, v, bs);
        for (i = 0; i < n; i++) {
            bs[i] = -es * bs[i];
            bg[i] = -ed->tg.occ[i] * v[i];
        }
        c1 = vdot(n, v, bs);
        c2 = vdot(n, v, bg);
        for (i = 0; i < n; i++) {
            bs[i] -= c1 * v[i];
            bg[i] -= c2 * v[i];
        }
        ed->cur_s = s;
        ed->cur_g = g;
        if (cg_reduced_resolvent(n, ed_matvec, ed, out->theta, v, bs, xs, ed->cg_tol,
                                 ed->cg_maxit, &ci1) != 0)
            return -1;
        if (cg_reduced_resolvent(n, ed_matvec, ed, out->theta, v, bg, xg, ed->cg_tol,
                                 ed->cg_maxit, &ci2) != 0)
            return -1;
        out->chi_kk = (es * vAv + 2.0 * vdot(n, bs, xs)) / L;
        out->chi_ee = (2.0 * vdot(n, bg, xg)) / L;
        out->chi_ke = (vdot(n, bg, xs) + vdot(n, bs, xg)) / L; /* symmetrized */
        out->cg_iters = ci1.iters + ci2.iters;
        out->cg_relres = (ci1.relres > ci2.relres) ? ci1.relres : ci2.relres;
        if (!ci1.converged || !ci2.converged) out->ok = 0;
        finish_rho(out);
        return 0;
    }

    if (chi == CHI_FD) {
        const double h = ed->fd_h;
        double *v_center = (double *)malloc(sizeof(double) * (size_t)n);
        EDPoint sp, sm, gp, gm;
        int ok = out->ok;
        if (!v_center) return -1;
        memcpy(v_center, ed->v, sizeof(double) * (size_t)n);
        if (solve_top(ed, s + h, g, &sp) != 0) { free(v_center); return -1; }
        memcpy(ed->v, v_center, sizeof(double) * (size_t)n);
        if (solve_top(ed, s - h, g, &sm) != 0) { free(v_center); return -1; }
        memcpy(ed->v, v_center, sizeof(double) * (size_t)n);
        if (solve_top(ed, s, g + h, &gp) != 0) { free(v_center); return -1; }
        memcpy(ed->v, v_center, sizeof(double) * (size_t)n);
        if (solve_top(ed, s, g - h, &gm) != 0) { free(v_center); return -1; }
        memcpy(ed->v, v_center, sizeof(double) * (size_t)n);
        free(v_center);
        out->chi_kk = -(sp.k - sm.k) / (2.0 * h);
        out->chi_ee = -(gp.e - gm.e) / (2.0 * h);
        /* dk/dg and de/ds both equal -chi_KE (Maxwell relation); average them */
        out->chi_ke = -0.5 * ((gp.k - gm.k) + (sp.e - sm.e)) / (2.0 * h);
        out->nmatvec += sp.nmatvec + sm.nmatvec + gp.nmatvec + gm.nmatvec;
        out->ok = ok && sp.ok && sm.ok && gp.ok && gm.ok;
        finish_rho(out);
        return 0;
    }
    return -1;
}

int ed_find_peak(EDSolver *ed, double g, const double *sgrid, int ns, double stol,
                 EDPoint *best, int *at_edge)
{
    const double gr = 0.5 * (sqrt(5.0) - 1.0);
    double *chi = NULL;
    double a, b, x1, x2, f1, f2;
    int i, imax = 0;
    EDPoint p;

    *at_edge = 0;
    if (ns < 1) return -1;
    chi = (double *)malloc(sizeof(double) * (size_t)ns);
    if (!chi) return -1;

    for (i = 0; i < ns; i++) {
        if (ed_point(ed, sgrid[i], g, CHI_LR, &p) != 0) {
            free(chi);
            return -1;
        }
        chi[i] = p.chi_kk;
        if (i == 0 || chi[i] > chi[imax]) imax = i;
    }
    if (ns == 1) {
        free(chi);
        return ed_point(ed, sgrid[0], g, CHI_LR, best);
    }
    if (imax == 0 || imax == ns - 1) *at_edge = 1;
    a = sgrid[imax > 0 ? imax - 1 : 0];
    b = sgrid[imax < ns - 1 ? imax + 1 : ns - 1];
    free(chi);

    x1 = b - gr * (b - a);
    x2 = a + gr * (b - a);
    if (ed_point(ed, x1, g, CHI_LR, &p) != 0) return -1;
    f1 = p.chi_kk;
    if (ed_point(ed, x2, g, CHI_LR, &p) != 0) return -1;
    f2 = p.chi_kk;
    while (b - a > stol) {
        if (f1 > f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - gr * (b - a);
            if (ed_point(ed, x1, g, CHI_LR, &p) != 0) return -1;
            f1 = p.chi_kk;
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + gr * (b - a);
            if (ed_point(ed, x2, g, CHI_LR, &p) != 0) return -1;
            f2 = p.chi_kk;
        }
    }
    return ed_point(ed, 0.5 * (a + b), g, CHI_LR, best);
}

int ed_find_strict_crossing(EDSolver *ed, double g, double s0, double stol, EDPoint *out)
{
    const int L = ed->model.L;
    double lo, hi, s, step;
    int it, have_lo = 0, have_hi = 0;
    EDPoint p;

    if (!ed->basis.exclude_empty) return -1;

    /* bracket: theta decreases with s; lo has theta > 0, hi has theta < 0 */
    s = s0;
    if (ed_point(ed, s, g, CHI_NONE, &p) != 0) return -1;
    if (fabs(p.theta) <= 1e-13) {
        *out = p;
        return 0;
    }
    step = 0.05;
    if (p.theta > 0.0) {
        lo = s;
        have_lo = 1;
        for (it = 0; it < 60 && !have_hi; it++) {
            s = lo + step;
            if (ed_point(ed, s, g, CHI_NONE, &p) != 0) return -1;
            if (p.theta < 0.0) { hi = s; have_hi = 1; }
            else { lo = s; step *= 2.0; }
        }
    } else {
        hi = s;
        have_hi = 1;
        for (it = 0; it < 60 && !have_lo; it++) {
            s = hi - step;
            if (ed_point(ed, s, g, CHI_NONE, &p) != 0) return -1;
            if (p.theta > 0.0) { lo = s; have_lo = 1; }
            else { hi = s; step *= 2.0; }
        }
    }
    if (!have_lo || !have_hi) return -1;

    /* safeguarded Newton, d theta/ds = -L k */
    s = 0.5 * (lo + hi);
    for (it = 0; it < 200; it++) {
        double snew, df;
        if (ed_point(ed, s, g, CHI_NONE, &p) != 0) return -1;
        if (fabs(p.theta) <= 1e-13) break;
        if (p.theta > 0.0) lo = s;
        else hi = s;
        df = -L * p.k;
        snew = (df != 0.0) ? s - p.theta / df : 0.5 * (lo + hi);
        if (!(snew > lo && snew < hi)) snew = 0.5 * (lo + hi);
        if (fabs(snew - s) < stol || hi - lo < stol) {
            s = snew;
            if (ed_point(ed, s, g, CHI_NONE, &p) != 0) return -1;
            break;
        }
        s = snew;
    }
    *out = p;
    return 0;
}
