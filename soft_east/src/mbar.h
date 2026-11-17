/*
 * MBAR (Shirts & Chodera 2008) for path ensembles biased by exp(-s K - g E),
 * optionally restricted to umbrella windows Kmin <= K <= Kmax.
 *
 * State k samples trajectories with reduced bias
 *   u_k(x) = s_k K(x) + g_k E(x)   if Kmin_k <= K(x) <= Kmax_k,  +inf otherwise,
 * relative to the unbiased path measure. The free energies f_k are found by
 * Newton's method on the convex MBAR objective
 *   F(f) = sum_n ln sum_k N_k exp(f_k - u_k(x_n)) - sum_k N_k f_k,  f_0 = 0.
 * Averages at any unrestricted (s, g) then follow from the weights
 *   w_n(s, g) ~ exp(-s K_n - g E_n) / sum_k N_k exp(f_k - u_k(x_n)).
 */
#ifndef MBAR_H
#define MBAR_H

#include "rng.h"

typedef struct {
    int nwin;
    double *ws, *wg;     /* bias of each state */
    double *Kmin, *Kmax; /* umbrella window of each state (-inf, +inf: none) */
    long *N;             /* samples per state */
    long ntot;
    double *K, *E;       /* samples, grouped by state */
    double *f;           /* free energies, f[0] = 0 */
    double *logden;      /* ln sum_k N_k exp(f_k - u_k(x_n)) */
    int iters;
    double grad_max;     /* final max |grad_k| / N_k */
} MBAR;

typedef struct {
    double meanK, meanE, varK, varE, covKE;
    double ess;          /* Kish effective sample size */
} MBARObs;

/* Averages at (s, g) on the two sides of a cut: lo = K < Kcut, hi = K >= Kcut. */
typedef struct {
    double p_lo;                      /* P(K < Kcut) */
    double meanK_lo, meanE_lo, varK_lo;
    double meanK_hi, meanE_hi, varK_hi; /* NAN on a side without weight */
} MBARSplit;

/* Allocates storage (windows default to none); fill ws, wg, Kmin, Kmax, N, K, E. */
int mbar_alloc(MBAR *m, int nwin, long ntot);
void mbar_free(MBAR *m);

/*
 * Initial guess for f along the state order: neighbouring states with the
 * same (s, g) and overlapping windows use the overlap of their samples;
 * otherwise trapezoidal integration of <K> ds + <E> dg.
 */
void mbar_guess(MBAR *m);

/* Newton iterations from the current f; returns 0 when max |grad|/N < tol. */
int mbar_solve(MBAR *m, double tol, int maxit);

/* Averages at (s, g) without window, or restricted to [Kmin, Kmax]. */
void mbar_observe(const MBAR *m, double s, double g, MBARObs *o);
void mbar_observe_window(const MBAR *m, double s, double g, double Kmin, double Kmax,
                         MBARObs *o);

/* Conditional averages at (s, g) on either side of Kcut (same cut rule as mbar_prob_below). */
void mbar_observe_split(const MBAR *m, double s, double g, double Kcut, MBARSplit *o);

/* Reweighted probability that K < Kcut at (s, g). */
double mbar_prob_below(const MBAR *m, double s, double g, double Kcut);

/* Reweighted density of K at (s, g) on [klo, khi) with nbins bins. */
void mbar_histogram(const MBAR *m, double s, double g, double klo, double khi, int nbins,
                    double *h);

/*
 * Moving-block bootstrap replica: dst gets the same states and sample counts
 * as src, each state's samples redrawn as consecutive blocks of `block`
 * samples with random starting points. dst->f is initialized from src->f;
 * call mbar_solve(dst, ...) afterwards. dst must be freed with mbar_free.
 */
int mbar_resample(const MBAR *src, MBAR *dst, int block, Rng *rng);

#endif
