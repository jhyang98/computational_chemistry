#include "lanczos.h"
#include "eig_dense.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void lanczos_default_opts(LanczosOpts *o)
{
    o->mmax = 100;
    o->max_restarts = 5000;
    o->tol = 1e-11;
}

static double vdot(long n, const double *a, const double *b)
{
    long i;
    double s = 0.0;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static double vnrm(long n, const double *a)
{
    return sqrt(vdot(n, a, a));
}

/* Deterministic pseudo-random fill in [-0.5, 0.5) (used only on breakdown). */
static void fill_random(long n, double *x, unsigned long long *state)
{
    long i;
    for (i = 0; i < n; i++) {
        *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
        x[i] = (double)(*state >> 11) * (1.0 / 9007199254740992.0) - 0.5;
    }
}

/* Classical Gram-Schmidt, applied twice, against columns 0..k-1 of V.
   Accumulated projection coefficients are written to h (may be NULL). */
static void orthogonalize(long n, int k, const double *V, double *w, double *h, double *tmp)
{
    int pass, i;
    if (h)
        for (i = 0; i < k; i++) h[i] = 0.0;
    for (pass = 0; pass < 2; pass++) {
        for (i = 0; i < k; i++) tmp[i] = vdot(n, V + (size_t)i * (size_t)n, w);
        for (i = 0; i < k; i++) {
            const double *vi = V + (size_t)i * (size_t)n;
            double ci = tmp[i];
            long r;
            for (r = 0; r < n; r++) w[r] -= ci * vi[r];
            if (h) h[i] += ci;
        }
    }
}

int lanczos_top(long n, matvec_fn A, void *ctx, const double *v0, int nev,
                double *evals, double *evecs, double *resid,
                const LanczosOpts *opt, LanczosInfo *info)
{
    int m = opt->mmax;
    int ld, kk, j, i, restart, full = 0, nw;
    double *V = NULL, *H = NULL, *Hc = NULL, *Y = NULL, *th = NULL, *h = NULL, *tmp = NULL;
    double *w = NULL, *U = NULL;
    double anorm = 0.0, beta_last = 0.0, nrm;
    unsigned long long rng = 0x9E3779B97F4A7C15ULL;
    long r;

    info->converged = 0;
    info->nmatvec = 0;
    info->nrestart = 0;

    if (n < 1 || nev < 1) return -1;
    if (m > n) m = (int)n;
    if (m < nev + 2 && n >= nev + 2) m = nev + 2;
    if (m > n) m = (int)n;
    ld = m + 1;

    V = (double *)calloc((size_t)n * (size_t)(m + 1), sizeof(double));
    U = (double *)calloc((size_t)n * (size_t)m, sizeof(double));
    H = (double *)calloc((size_t)ld * (size_t)ld, sizeof(double));
    Hc = (double *)malloc(sizeof(double) * (size_t)m * (size_t)m);
    Y = (double *)malloc(sizeof(double) * (size_t)m * (size_t)m);
    th = (double *)malloc(sizeof(double) * (size_t)m);
    h = (double *)malloc(sizeof(double) * (size_t)(m + 1));
    tmp = (double *)malloc(sizeof(double) * (size_t)(m + 1));
    w = (double *)malloc(sizeof(double) * (size_t)n);
    if (!V || !U || !H || !Hc || !Y || !th || !h || !tmp || !w) {
        free(V); free(U); free(H); free(Hc); free(Y); free(th); free(h); free(tmp); free(w);
        return -1;
    }

    /* starting vector */
    if (v0) memcpy(V, v0, sizeof(double) * (size_t)n);
    else
        for (r = 0; r < n; r++) V[r] = 1.0;
    nrm = vnrm(n, V);
    if (!(nrm > 0.0)) {
        for (r = 0; r < n; r++) V[r] = 1.0;
        nrm = vnrm(n, V);
    }
    for (r = 0; r < n; r++) V[r] /= nrm;

    kk = 0;
    for (restart = 0;; restart++) {
        int mcur = m;
        full = 0;
        for (j = kk; j < m; j++) {
            double *vj = V + (size_t)j * (size_t)n;
            double beta, wn;
            A(vj, w, ctx);
            info->nmatvec++;
            wn = vnrm(n, w);
            if (wn > anorm) anorm = wn;
            orthogonalize(n, j + 1, V, w, h, tmp);
            for (i = 0; i <= j; i++) {
                H[i * ld + j] = h[i];
                H[j * ld + i] = h[i];
            }
            if (j + 1 == n) { /* basis spans the whole space */
                full = 1;
                mcur = j + 1;
                beta_last = 0.0;
                break;
            }
            beta = vnrm(n, w);
            if (beta <= 1e-13 * anorm) {
                /* invariant subspace: continue with a fresh orthogonal direction */
                int tries;
                double *vn = V + (size_t)(j + 1) * (size_t)n;
                for (tries = 0; tries < 10; tries++) {
                    fill_random(n, vn, &rng);
                    orthogonalize(n, j + 1, V, vn, NULL, tmp);
                    nrm = vnrm(n, vn);
                    if (nrm > 1e-8) break;
                }
                for (r = 0; r < n; r++) vn[r] /= nrm;
                H[(j + 1) * ld + j] = 0.0;
                H[j * ld + (j + 1)] = 0.0;
            } else {
                double *vn = V + (size_t)(j + 1) * (size_t)n;
                for (r = 0; r < n; r++) vn[r] = w[r] / beta;
                H[(j + 1) * ld + j] = beta;
                H[j * ld + (j + 1)] = beta;
            }
        }
        if (!full) beta_last = H[m * ld + (m - 1)];

        /* Rayleigh-Ritz on the projected matrix */
        for (i = 0; i < mcur; i++)
            for (j = 0; j < mcur; j++) Hc[i * mcur + j] = H[i * ld + j];
        eig_sym_jacobi(mcur, Hc, th, Y);

        nw = (nev < mcur) ? nev : mcur;
        {
            int est_ok = 1;
            for (i = 0; i < nw; i++) {
                double est = fabs(beta_last * Y[(mcur - 1) * mcur + i]);
                if (est > opt->tol) est_ok = 0;
            }
            if (est_ok || full || restart >= opt->max_restarts) {
                int all_ok = 1;
                for (i = 0; i < nw; i++) {
                    double *x = evecs + (size_t)i * (size_t)n;
                    double sum = 0.0;
                    memset(x, 0, sizeof(double) * (size_t)n);
                    for (j = 0; j < mcur; j++) {
                        double yj = Y[j * mcur + i];
                        const double *vj = V + (size_t)j * (size_t)n;
                        for (r = 0; r < n; r++) x[r] += yj * vj[r];
                    }
                    nrm = vnrm(n, x);
                    for (r = 0; r < n; r++) {
                        x[r] /= nrm;
                        sum += x[r];
                    }
                    if (sum < 0.0)
                        for (r = 0; r < n; r++) x[r] = -x[r];
                    A(x, w, ctx);
                    info->nmatvec++;
                    evals[i] = vdot(n, x, w); /* Rayleigh quotient */
                    for (r = 0; r < n; r++) w[r] -= evals[i] * x[r];
                    resid[i] = vnrm(n, w);
                    if (resid[i] > opt->tol) all_ok = 0;
                }
                for (i = nw; i < nev; i++) {
                    evals[i] = NAN;
                    resid[i] = NAN;
                }
                if (all_ok || full || restart >= opt->max_restarts) {
                    info->converged = all_ok || full;
                    info->nrestart = restart;
                    break;
                }
            }
        }

        /* thick restart: keep the leading Ritz vectors plus the residual direction */
        {
            int keep = nev + (mcur - nev) / 2;
            if (keep > mcur - 1) keep = mcur - 1;
            if (keep < nw) keep = nw;
            for (i = 0; i < keep; i++) {
                double *u = U + (size_t)i * (size_t)n;
                memset(u, 0, sizeof(double) * (size_t)n);
                for (j = 0; j < mcur; j++) {
                    double yj = Y[j * mcur + i];
                    const double *vj = V + (size_t)j * (size_t)n;
                    for (r = 0; r < n; r++) u[r] += yj * vj[r];
                }
            }
            /* continuation vector V[mcur] -> V[keep] */
            memmove(V + (size_t)keep * (size_t)n, V + (size_t)mcur * (size_t)n,
                    sizeof(double) * (size_t)n);
            memcpy(V, U, sizeof(double) * (size_t)n * (size_t)keep);
            memset(H, 0, sizeof(double) * (size_t)ld * (size_t)ld);
            for (i = 0; i < keep; i++) {
                H[i * ld + i] = th[i];
                H[keep * ld + i] = beta_last * Y[(mcur - 1) * mcur + i];
                H[i * ld + keep] = H[keep * ld + i];
            }
            kk = keep;
        }
    }

    free(V); free(U); free(H); free(Hc); free(Y); free(th); free(h); free(tmp); free(w);
    return 0;
}
