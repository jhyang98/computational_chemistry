/*
 * Dense symmetric eigensolver (cyclic Jacobi rotations).
 * Used for the projected matrix inside Lanczos and for small-size checks.
 */
#ifndef EIG_DENSE_H
#define EIG_DENSE_H

/*
 * a : n x n symmetric matrix, row-major; destroyed on output.
 * w : eigenvalues, sorted in descending order.
 * z : eigenvectors, row-major, column k is the eigenvector of w[k]
 *     (z[i * n + k] = component i). May be NULL if not needed.
 * Returns the number of sweeps used, or -1 if not converged.
 */
int eig_sym_jacobi(int n, double *a, double *w, double *z);

#endif
