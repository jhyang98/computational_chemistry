/*
 * Transition path sampling of the (s, g)-biased softened East model.
 *
 * Trajectory ensemble: P[x] ~ pi(x_0) prod_t T(x_t | x_{t-1}) exp(-s K - g E),
 * with T the MCS kernel of dyn.h and pi the equilibrium product measure.
 * Moves: forward shooting from a uniform slice tau in [0, tobs-1] or backward
 * shooting from tau in [1, tobs], each with probability 1/2, regenerating the
 * segment with fresh noise from the unbiased kernel. Acceptance:
 * min(1, exp(-s dK - g dE)).
 *
 * Umbrella window (optional): K is restricted to [Kmin, Kmax] with
 * Kmin = ceil(umb_kmin L tobs), Kmax = floor(umb_kmax L tobs). Moves that
 * would leave the window are rejected. If the starting trajectory lies
 * outside, an approach phase first accepts only moves that bring K closer to
 * the window (not part of the sampled ensemble), then n_equil and n_prod
 * moves follow inside the window.
 */
#ifndef TPS_H
#define TPS_H

#include <stdint.h>

#include "dyn.h"
#include "stats.h"

typedef struct {
    int L;
    double T, eps, s, g;
    int tobs;
    uint64_t seed;      /* 0: from the clock */
    long n_equil;       /* shooting moves before sampling */
    long n_prod;        /* shooting moves while sampling */
    int sample_every;   /* moves between recorded samples */
    const char *init;   /* "equilibrium", "empty" or a trajectory file */
    const char *save;   /* final trajectory file, or NULL */
    int nblocks;        /* blocks for jackknife errors */
    long progress;      /* print progress every this many moves (0: never) */
    double umb_kmin;    /* umbrella window in k = K/(L tobs); NAN: no bound */
    double umb_kmax;
    long n_approach;    /* maximum moves to enter the window */
} TPSParams;

typedef struct {
    long nsamp;
    double *K, *E;      /* recorded samples */
    long *move;         /* move index of each sample */
    long n_fwd, n_bwd, acc_fwd, acc_bwd;
    PairStats st;
    uint64_t seed_used;
    int umbrella;       /* 1 if a K window was applied */
    long Kmin, Kmax;    /* window in K (0 .. L tobs without umbrella) */
    long n_approach_used;
} TPSResult;

void tps_default_params(TPSParams *p);

/*
 * One shooting move. With on = 1, moves that would leave [Kmin, Kmax] are
 * rejected, and while the trajectory is outside the window only moves that
 * bring K closer are accepted (approach phase). Returns 1 if accepted.
 */
int tps_shoot(Traj *tr, Shot *sh, const Dyn *d, double s, double g, int on, long Kmin, long Kmax,
              Rng *r, int *dir_out);
int tps_run(const TPSParams *p, TPSResult *res);
void tps_result_free(TPSResult *res);

#endif
