/*
 * Replica-exchange umbrella sampling (REUS) in the activity K.
 *
 * nwin umbrella windows at the same (s, g): window w restricts K to
 *   [ceil(k_lo + w D) L tobs, floor(k_lo + (w+2) D) L tobs],
 *   D = (k_hi - k_lo) / (nwin + 1),
 * so neighbouring windows overlap by half a window. One replica (trajectory)
 * sits in each window. A cycle is swap_every shooting moves per replica
 * (tps_shoot, run in parallel with OpenMP if enabled), followed by swap
 * attempts between neighbouring windows (even pairs on even cycles, odd pairs
 * on odd cycles). A swap of replicas a (in w) and b (in w+1) is accepted with
 *   min(1, exp(u_w(a) + u_{w+1}(b) - u_w(b) - u_{w+1}(a))),
 * u_w(x) = s K + g E inside window w and +inf outside; with equal s and g
 * this means: accepted whenever each trajectory lies in the other's window.
 *
 * Each replica first enters its starting window (approach phase of
 * tps_shoot, at most n_approach moves). A replica that does not get there
 * restarts from a copy of the trajectory of a replica that did, the one whose
 * K is closest to the nearest missing frontier, and approaches again; this is
 * repeated until every window is reached or three rounds in a row bring no
 * new window in. Then n_equil moves
 * with swaps, then n_prod moves with swaps while the (K, E) of the replica in
 * each window is recorded every sample_every moves.
 *
 * Checkpoints (optional): with a checkpoint file set, the full state
 * (trajectories, random-number states, window assignment, counters and the
 * samples so far) is written during and after the approach phase, every
 * checkpoint_every seconds of wall time, and when the run is asked to stop
 * (stop_after seconds elapsed, or SIGTERM / SIGUSR1). reus_run then returns
 * REUS_STOPPED; called again with the same parameters it resumes from the
 * file and gives exactly the same result as an uninterrupted run. The file
 * is written to <file>.tmp first and then renamed.
 */
#ifndef REUS_H
#define REUS_H

#include <stdint.h>

#include "dyn.h"

typedef struct {
    int L;
    double T, eps, s, g;
    int tobs;
    uint64_t seed;     /* 0: from the clock */
    int nwin;
    double k_lo, k_hi; /* k range covered by the windows */
    long n_approach;   /* maximum moves per replica to enter its window */
    long n_equil;      /* moves per replica before sampling */
    long n_prod;       /* moves per replica while sampling */
    int swap_every;    /* moves between swap attempts */
    int sample_every;  /* moves between samples (multiple of swap_every) */
    const char *init;  /* "empty", "equilibrium", or a file pattern with one %d
                          (window index); a missing file falls back to empty */
    long progress;     /* report every this many cycles (0: never) */
    const char *checkpoint;  /* checkpoint file (NULL or "": none) */
    double checkpoint_every; /* seconds of wall time between checkpoints (0: only on stop) */
    double stop_after;       /* write a checkpoint and stop after this many seconds (0: never) */
    long stop_approach_round; /* testing: stop after this many approach retry rounds (0: never) */
    long stop_cycle;         /* testing: checkpoint and stop once this many cycles are done (0: never) */
} REUSParams;

#define REUS_STOPPED 3 /* reus_run: stopped after writing a checkpoint, not finished */

typedef struct {
    int nwin;
    long *Kmin, *Kmax;   /* window bounds in K */
    long nsamp;          /* samples per window */
    double *K, *E;       /* sample i of window w at [w * nsamp + i] */
    int *rep;            /* replica that produced each sample */
    long *moves, *acc;   /* shooting moves per window (production) */
    long *swap_try, *swap_acc; /* per neighbour pair (w, w+1) */
    long *round_trips;   /* per replica: window 0 -> nwin-1 -> 0 */
    long *approach;      /* per replica: moves needed to enter the window (all attempts) */
    int *init_from_file; /* per replica: 1 if started from a trajectory file */
    int *seeded_from;    /* per replica: window whose trajectory it restarted from, -1 none */
    Traj *final;         /* final trajectory in each window */
    uint64_t seed_used;
} REUSResult;

void reus_default_params(REUSParams *p);

/* 0: finished (free res with reus_result_free); REUS_STOPPED: stopped, resume
   later from the checkpoint (res holds nothing); -1: error. */
int reus_run(const REUSParams *p, REUSResult *res);
void reus_result_free(REUSResult *res);

/* Window bounds in K for window w. */
void reus_window(const REUSParams *p, int w, long *Kmin, long *Kmax);

/*
 * Kolmogorov-Smirnov distance between the K distributions of windows w and
 * w+1 restricted to their overlap (should be small when converged); the
 * numbers of samples of each window in the overlap are returned in n1, n2.
 */
double reus_overlap_ks(const REUSResult *res, int w, long *n1, long *n2);

/* Effective number of independent samples of K in window w (block estimate). */
double reus_neff(const REUSResult *res, int w);

#endif
