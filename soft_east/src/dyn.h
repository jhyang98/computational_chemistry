/*
 * Discrete-time Monte Carlo dynamics of the softened East model and
 * trajectory storage for path sampling.
 *
 * One Monte Carlo sweep (MCS) = L update attempts. Each attempt picks a site
 * i uniformly at random (with replacement) and flips it with probability
 *
 *   C_i * c        if n_i = 0
 *   C_i * (1 - c)  if n_i = 1,     C_i = n_{i-1} + eps (1 - n_{i-1}),
 *
 * with c = 1 / (1 + exp(1/T)) and periodic boundaries. This is the update of
 * the earlier East-model TPS code (facilitating site = left neighbour)
 * extended by the softness eps. The single-attempt kernel is reversible with
 * respect to the product measure of density c, so a sweep is reversible too,
 * which makes backward shooting with the forward kernel valid.
 *
 * Trajectory conventions (t_obs = number of MCS):
 *   slices x_0, ..., x_{t_obs}
 *   E = sum_{t=0}^{t_obs} n(x_t)       (all t_obs + 1 slices)
 *   K = sum_{t=1}^{t_obs} k_t          (k_t = flips during sweep t-1 -> t)
 */
#ifndef DYN_H
#define DYN_H

#include "rng.h"

typedef struct {
    int L;
    double T, c, eps;
    double pflip[2][2]; /* [n_{i-1}][n_i] -> flip probability per attempt */
} Dyn;

int dyn_init(Dyn *d, int L, double T, double eps);

/* Draw a configuration from the product measure with density c. Returns n. */
int dyn_sample_equilibrium(const Dyn *d, unsigned char *x, Rng *r);

/* One MCS in place. Updates *n (number of excitations); returns flips. */
int dyn_mcs(const Dyn *d, unsigned char *x, int *n, Rng *r);

/*
 * Run nsteps MCS starting from x_start (not modified). Slice j of out
 * (out + j*L) is the configuration after j+1 sweeps; n_out[j] its number of
 * excitations and k_out[j] the flips during that sweep.
 */
void dyn_propagate(const Dyn *d, const unsigned char *x_start, int nsteps, unsigned char *out,
                   int *n_out, int *k_out, Rng *r);

typedef struct {
    int L, tobs;
    unsigned char *x; /* (tobs + 1) * L */
    int *n;           /* n[t], t = 0..tobs */
    int *k;           /* k[t], t = 1..tobs (k[0] = 0) */
    long K, E;
} Traj;

int traj_alloc(Traj *tr, int L, int tobs);
void traj_free(Traj *tr);
void traj_recount(Traj *tr); /* recompute n[], K, E from the slices and k[] */

/* Fresh unbiased trajectory: x_0 from the product measure, then forward. */
void traj_init_equilibrium(Traj *tr, const Dyn *d, Rng *r);

/* Text file: header lines, then one line "k_t slice" per slice. */
int traj_save(const Traj *tr, const char *path);
int traj_load(Traj *tr, const char *path);

/*
 * Shooting move workspace and proposal. A proposal regenerates the segment
 * after (forward) or before (backward) slice tau with fresh noise and stores
 * it in the workspace; dK and dE are the changes of K and E.
 */
typedef struct {
    unsigned char *buf; /* tobs * L */
    int *nbuf, *kbuf;   /* tobs */
    int dir, tau, len;
    long dK, dE;
} Shot;

int shot_alloc(Shot *sh, int L, int tobs);
void shot_free(Shot *sh);
void shot_propose(Shot *sh, const Traj *tr, const Dyn *d, int dir, int tau, Rng *r);
void shot_accept(const Shot *sh, Traj *tr);

#endif
