/*
 * Thick-restart Lanczos with full reorthogonalization for the largest
 * eigenvalues of a real symmetric operator.
 */
#ifndef LANCZOS_H
#define LANCZOS_H

typedef void (*matvec_fn)(const double *x, double *y, void *ctx);

typedef struct {
    int mmax;          /* maximum Krylov basis size (capped at n) */
    int max_restarts;
    double tol;        /* absolute residual tolerance ||A x - theta x|| */
} LanczosOpts;

typedef struct {
    int converged;     /* 1 if all wanted pairs met the tolerance */
    int nmatvec;
    int nrestart;
} LanczosInfo;

void lanczos_default_opts(LanczosOpts *o);

/*
 * n      : dimension
 * A, ctx : operator y = A x
 * v0     : starting vector (length n), NULL for a uniform vector
 * nev    : number of wanted eigenpairs (largest algebraic)
 * evals  : out, nev eigenvalues in descending order
 * evecs  : out, nev normalized eigenvectors, evecs[k * n + i]
 * resid  : out, explicit residual norms of the returned pairs
 * Returns 0 on success, -1 on allocation failure.
 */
int lanczos_top(long n, matvec_fn A, void *ctx, const double *v0, int nev,
                double *evals, double *evecs, double *resid,
                const LanczosOpts *opt, LanczosInfo *info);

#endif
