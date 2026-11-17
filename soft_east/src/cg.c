#include "cg.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static double vdot(long n, const double *a, const double *b)
{
    long i;
    double s = 0.0;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static void project_out(long n, const double *v, double *x)
{
    double c = vdot(n, v, x);
    long i;
    for (i = 0; i < n; i++) x[i] -= c * v[i];
}

int cg_reduced_resolvent(long n, matvec_fn H, void *ctx, double theta0, const double *v,
                         const double *b, double *x, double tol, int maxit, CGInfo *info)
{
    double *r = NULL, *p = NULL, *q = NULL;
    double bnorm, rr, rr_new, alpha, beta, pq;
    long i;
    int it;

    info->iters = 0;
    info->relres = 0.0;
    info->converged = 0;
    memset(x, 0, sizeof(double) * (size_t)n);

    r = (double *)malloc(sizeof(double) * (size_t)n);
    p = (double *)malloc(sizeof(double) * (size_t)n);
    q = (double *)malloc(sizeof(double) * (size_t)n);
    if (!r || !p || !q) {
        free(r); free(p); free(q);
        return -1;
    }

    memcpy(r, b, sizeof(double) * (size_t)n);
    project_out(n, v, r);
    bnorm = sqrt(vdot(n, r, r));
    if (bnorm == 0.0) {
        info->converged = 1;
        free(r); free(p); free(q);
        return 0;
    }
    memcpy(p, r, sizeof(double) * (size_t)n);
    rr = bnorm * bnorm;

    for (it = 1; it <= maxit; it++) {
        H(p, q, ctx);
        for (i = 0; i < n; i++) q[i] = theta0 * p[i] - q[i];
        project_out(n, v, q);
        pq = vdot(n, p, q);
        if (!(pq > 0.0)) break; /* loss of positive definiteness */
        alpha = rr / pq;
        for (i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        project_out(n, v, r);
        rr_new = vdot(n, r, r);
        info->iters = it;
        info->relres = sqrt(rr_new) / bnorm;
        if (info->relres <= tol) {
            info->converged = 1;
            break;
        }
        beta = rr_new / rr;
        rr = rr_new;
        for (i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
    }
    project_out(n, v, x);

    /* report the true residual */
    H(x, q, ctx);
    for (i = 0; i < n; i++) q[i] = theta0 * x[i] - q[i];
    project_out(n, v, q);
    memcpy(r, b, sizeof(double) * (size_t)n);
    project_out(n, v, r);
    for (i = 0; i < n; i++) r[i] -= q[i];
    info->relres = sqrt(vdot(n, r, r)) / bnorm;
    if (info->relres > 10.0 * tol) info->converged = 0;

    free(r); free(p); free(q);
    return 0;
}
