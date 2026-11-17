/*
 * Exact diagonalization driver for the (s, g)-biased softened East model.
 *
 * Conventions (trajectory weight exp(-s K - g E), theta = largest eigenvalue
 * of the tilted generator, whole system, per unit time):
 *
 *   k      =  K / (t L)       = -(1/L) d theta / ds
 *   e      =  E / (t L)       = -(1/L) d theta / dg
 *   chi_KK =  Var K / (t L)   =  (1/L) d^2 theta / ds^2
 *   chi_KE =  Cov(K,E)/(t L)  =  (1/L) d^2 theta / ds dg
 *   chi_EE =  Var E / (t L)   =  (1/L) d^2 theta / dg^2
 *
 * First derivatives use the Hellmann-Feynman theorem; second derivatives use
 * second-order perturbation theory with the reduced resolvent (CG), or
 * central finite differences of the first derivatives (for cross-checks).
 */
#ifndef ED_H
#define ED_H

#include "basis.h"
#include "lanczos.h"
#include "model.h"
#include "tgen.h"

typedef enum { CHI_NONE = 0, CHI_LR = 1, CHI_FD = 2 } ChiMode;

typedef struct {
    double s, g;
    double theta;    /* largest eigenvalue in the basis sector */
    double theta2;   /* second largest eigenvalue in the same sector */
    double k, e;
    double chi_kk, chi_ke, chi_ee, rho;
    double resid;    /* Lanczos residual of the top pair */
    double resid2;   /* Lanczos residual of the second pair */
    double cg_relres;
    double min_ratio; /* min/max component of the top eigenvector (>= 0 for Perron) */
    int nmatvec;
    int cg_iters;
    int ok;           /* 1 if every solver met its tolerance */
} EDPoint;

typedef struct {
    Model model;
    Basis basis;
    TGen tg;
    LanczosOpts lopt;
    double cg_tol;
    int cg_maxit;
    double fd_h;
    double *v;        /* warm-start vector */
    double *evecs;    /* workspace for two eigenvectors */
    double *w1, *w2, *w3, *w4;
    double cur_s, cur_g;
} EDSolver;

int ed_init(EDSolver *ed, int L, double T, double eps, int use_translation);
void ed_free(EDSolver *ed);

/* Reset the warm-start vector to the equilibrium Perron vector sqrt(pi). */
void ed_reset_guess(EDSolver *ed);

/* Solve at (s, g). Updates the warm-start vector. Returns 0 on success. */
int ed_point(EDSolver *ed, double s, double g, ChiMode chi, EDPoint *out);

/*
 * Pseudo-critical point: maximize chi_KK(s) at fixed g.
 * The coarse grid sgrid[0..ns-1] (increasing) is scanned first, then the
 * maximum is refined by golden-section search to |delta s| < stol.
 * *at_edge is set to 1 if the coarse maximum lies on the grid boundary.
 */
int ed_find_peak(EDSolver *ed, double g, const double *sgrid, int ns, double stol,
                 EDPoint *best, int *at_edge);

/*
 * Strict East (eps = 0, empty configuration excluded): coexistence point
 * where the largest eigenvalue of the non-empty sector crosses the empty-state
 * eigenvalue 0, i.e. theta(s*, g) = 0. s0 is the initial guess.
 */
int ed_find_strict_crossing(EDSolver *ed, double g, double s0, double stol, EDPoint *out);

#endif
