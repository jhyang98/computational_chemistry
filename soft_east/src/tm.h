/*
 * Exact results for the discrete-time (MCS) dynamics of dyn.h.
 *
 * Single attempt, symmetrized with sqrt(pi) and tilted by exp(-s) per flip:
 *   M_s = I + (exp(-s) A - diag(r)) / L
 * (A and r as in tgen.h). One sweep is M_s^L; the energy weight of a slice is
 * D_g = diag(exp(-g n)). For the TPS ensemble of dyn.h,
 *
 *   Z(s, g) = psi0^T (D_g M_s^L)^{tobs} D_g psi0,   psi0 = sqrt(pi),
 *
 * so <K> = -d ln Z/ds, <E> = -d ln Z/dg, Var K = d2 ln Z/ds2, etc. These are
 * computed by forward-mode differentiation (exact up to round-off).
 *
 * Asymptotically, theta_d(s, g) = ln lambda_max(D^{1/2} M_s^L D^{1/2}) per MCS
 * (whole system); k_d = -(1/L) d theta_d/ds, e_d = -(1/L) d theta_d/dg.
 */
#ifndef TM_H
#define TM_H

#include <stdint.h>

#include "basis.h"
#include "lanczos.h"
#include "model.h"
#include "tgen.h"

typedef struct {
    Model model;
    Basis basis;
    TGen tg;
    double *psi0;      /* sqrt(pi) in the basis */
    double *work;      /* scratch, 16 vectors */
    LanczosOpts lopt;
    double cur_s, cur_g;
} TMSolver;

/* exclude_empty should be 1 only for eps = 0 asymptotic (non-empty sector). */
int tm_init(TMSolver *t, int L, double T, double eps, int use_translation, int exclude_empty);
void tm_free(TMSolver *t);

typedef struct {
    double lnZ;
    double meanK, meanE, varK, varE, covKE;
} TMFinite;

int tm_finite(TMSolver *t, double s, double g, int tobs, TMFinite *out);

typedef struct {
    double theta;       /* ln lambda_1, whole system per MCS */
    double lam_ratio;   /* lambda_2 / lambda_1 */
    double k, e;        /* per site per MCS */
    double resid;
    int ok;
} TMAsym;

int tm_asym(TMSolver *t, double s, double g, TMAsym *out);

/*
 * Exact distribution of K in the TPS ensemble at s = 0 and energy bias g:
 * logp[K] = ln P(K) for K = 0 .. L*tobs (-INFINITY where P(K) = 0 or
 * underflows). Other s follow from ln P_s(K) = ln P(K) - s K + const.
 * Cost grows like (L tobs)^2 * L * dim, so this is for small systems.
 */
int tm_pk(TMSolver *t, double g, int tobs, double *logp);

/*
 * Probability distribution after one unbiased MCS starting from configuration
 * x (requires a full basis: use_translation = 0, exclude_empty = 0).
 * prob has 2^L entries indexed by configuration.
 */
int tm_mcs_row(TMSolver *t, uint32_t x, double *prob);

#endif
