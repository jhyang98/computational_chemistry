/*
 * Conjugate gradient for the reduced resolvent of the top eigenpair:
 *
 *   (theta0 - H) x = b,   with b and x orthogonal to v,
 *
 * where (theta0, v) is the largest eigenpair of the symmetric operator H.
 * On the complement of v the matrix theta0 - H is positive definite, so CG
 * applies. Used for second derivatives of theta (second-order perturbation).
 */
#ifndef CG_H
#define CG_H

#include "lanczos.h"

typedef struct {
    int iters;
    double relres;     /* final ||r|| / ||b|| */
    int converged;
} CGInfo;

int cg_reduced_resolvent(long n, matvec_fn H, void *ctx, double theta0, const double *v,
                         const double *b, double *x, double tol, int maxit, CGInfo *info);

#endif
