/*
 * Symmetrized tilted generator of the softened East model.
 *
 * The trajectory weight is exp(-s K - g E), with K the number of spin flips
 * and E = int_0^t sum_i n_i(t') dt'. The tilted generator is
 *
 *   W(s,g)[x', x] = exp(-s) w(x -> x')            (x' != x)
 *   W(s,g)[x,  x] = -r(x) - g n(x)
 *
 * Detailed balance makes it similar to the symmetric matrix
 *
 *   H(s,g) = exp(-s) A - diag(r) - g diag(n),
 *   A[x', x] = sqrt(w(x -> x') w(x' -> x)) = sqrt(c (1 - c)) C_i,
 *
 * which has the same spectrum. Only H is stored, in the chosen basis.
 * theta(s,g) = largest eigenvalue of H (whole system, per unit time).
 */
#ifndef TGEN_H
#define TGEN_H

#include "basis.h"
#include "model.h"

typedef struct {
    long dim;
    long nnz;
    long *rowptr;  /* CSR of the off-diagonal coupling A (symmetric) */
    long *col;
    double *val;
    double *esc;   /* r(x): escape rate (diagonal) */
    double *occ;   /* n(x): number of excitations (diagonal) */
} TGen;

int tgen_build(TGen *t, const Model *m, const Basis *b);
void tgen_free(TGen *t);

/* y = A x */
void tgen_apply_A(const TGen *t, const double *x, double *y);

/* y = H(s,g) x */
void tgen_apply(const TGen *t, double s, double g, const double *x, double *y);

/* Largest |A[i][j] - A[j][i]| (should be at round-off level). */
double tgen_asymmetry(const TGen *t);

#endif
