/*
 * Discrete-time cloning (population dynamics) for the (s, g)-biased softened
 * East model, as in the earlier East-model cloning code, with constant
 * population and systematic resampling.
 *
 * Each step, every clone runs one unbiased MCS and gets the weight
 * exp(-s k - g n), where k is its number of flips in the sweep and n its
 * number of excitations after the sweep. The population is then resampled in
 * proportion to the weights. After discarding t_skip steps,
 *
 *   theta_d = (1/(tobs - t_skip)) sum_t ln <w>_t     (per MCS, whole system)
 *
 * estimates ln lambda_max of the tilted sweep transfer matrix (east_tm asym).
 * k and e come from the histories carried along by the clones: the flips and
 * excitations accumulated in the steps t_skip < t <= tobs - t_tail by the
 * ancestors of the final population. Leaving out the last t_tail steps
 * removes the final-time boundary layer, where the statistics differ from the
 * bulk of the biased trajectories; t_tail should exceed the relaxation time.
 */
#ifndef CLONE_H
#define CLONE_H

#include <stdint.h>

typedef struct {
    int L;
    double T, eps, s, g;
    int n_clones;
    int tobs;      /* cloning steps (MCS) */
    int t_skip;    /* initial steps excluded from the estimates */
    int t_tail;    /* final steps excluded from the k, e histories */
    int n_runs;    /* independent runs for error bars */
    uint64_t seed; /* 0: from the clock */
} CloneParams;

typedef struct {
    double theta, theta_err;
    double k, k_err, e, e_err;
    uint64_t seed_used;
} CloneResult;

void clone_default_params(CloneParams *p);

/* Optional per-run output: theta_runs, k_runs, e_runs may be NULL. */
int clone_run(const CloneParams *p, CloneResult *res, double *theta_runs, double *k_runs,
              double *e_runs);

#endif
