#include "tps.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void tps_default_params(TPSParams *p)
{
    memset(p, 0, sizeof(*p));
    p->L = 100;
    p->T = 1.0;
    p->eps = 0.0;
    p->s = 0.0;
    p->g = 0.0;
    p->tobs = 400;
    p->seed = 0;
    p->n_equil = 10000;
    p->n_prod = 100000;
    p->sample_every = 10;
    p->init = "equilibrium";
    p->save = NULL;
    p->nblocks = 20;
    p->progress = 0;
    p->umb_kmin = NAN;
    p->umb_kmax = NAN;
    p->n_approach = 1000000;
}

typedef struct {
    int on;
    long Kmin, Kmax;
} Window;

static long window_dist(const Window *w, long K)
{
    if (K < w->Kmin) return w->Kmin - K;
    if (K > w->Kmax) return K - w->Kmax;
    return 0;
}

int tps_shoot(Traj *tr, Shot *sh, const Dyn *d, double s, double g, int on, long Kmin, long Kmax,
              Rng *r, int *dir_out)
{
    int dir = (rng_uniform(r) < 0.5) ? 1 : -1;
    int tau = (dir > 0) ? rng_int(r, tr->tobs) : 1 + rng_int(r, tr->tobs);
    double logacc;
    *dir_out = dir;
    shot_propose(sh, tr, d, dir, tau, r);
    if (on) {
        Window w;
        long dcur, dnew;
        w.on = 1;
        w.Kmin = Kmin;
        w.Kmax = Kmax;
        dcur = window_dist(&w, tr->K);
        dnew = window_dist(&w, tr->K + sh->dK);
        if (dcur > 0) { /* approach phase: only get closer */
            if (dnew < dcur) {
                shot_accept(sh, tr);
                return 1;
            }
            return 0;
        }
        if (dnew > 0) return 0; /* would leave the window */
    }
    logacc = -s * (double)sh->dK - g * (double)sh->dE;
    if (logacc >= 0.0 || rng_uniform(r) < exp(logacc)) {
        shot_accept(sh, tr);
        return 1;
    }
    return 0;
}

static int tps_move(Traj *tr, Shot *sh, const Dyn *d, double s, double g, const Window *w,
                    Rng *r, int *dir_out)
{
    return tps_shoot(tr, sh, d, s, g, w->on, w->Kmin, w->Kmax, r, dir_out);
}

int tps_run(const TPSParams *p, TPSResult *res)
{
    Dyn d;
    Rng r;
    Traj tr;
    Shot sh;
    Window w;
    long m, cap;
    int dir;

    memset(res, 0, sizeof(*res));
    if (p->tobs < 1 || p->sample_every < 1 || p->n_prod < 1) {
        fprintf(stderr, "error: tobs, sample_every and n_prod must be positive\n");
        return -1;
    }
    if (dyn_init(&d, p->L, p->T, p->eps) != 0) {
        fprintf(stderr, "error: invalid L, T or eps\n");
        return -1;
    }
    res->seed_used = rng_seed_auto(&r, p->seed);
    if (traj_alloc(&tr, p->L, p->tobs) != 0 || shot_alloc(&sh, p->L, p->tobs) != 0) {
        fprintf(stderr, "error: out of memory\n");
        return -1;
    }

    if (!p->init || !strcmp(p->init, "equilibrium")) {
        traj_init_equilibrium(&tr, &d, &r);
    } else if (!strcmp(p->init, "empty")) {
        /* every slice empty: a valid (inactive) trajectory when eps > 0 */
        traj_recount(&tr);
    } else if (traj_load(&tr, p->init) != 0) {
        traj_free(&tr);
        shot_free(&sh);
        return -1;
    }
    if (p->eps == 0.0 && tr.E == 0)
        fprintf(stderr, "warning: eps = 0 and the trajectory is empty; shooting cannot leave it\n");

    /* umbrella window */
    {
        const double LT = (double)p->L * p->tobs;
        w.on = !isnan(p->umb_kmin) || !isnan(p->umb_kmax);
        w.Kmin = isnan(p->umb_kmin) ? 0 : (long)ceil(p->umb_kmin * LT - 1e-9);
        w.Kmax = isnan(p->umb_kmax) ? (long)LT : (long)floor(p->umb_kmax * LT + 1e-9);
        if (w.Kmin < 0) w.Kmin = 0;
        if (w.Kmax > (long)LT) w.Kmax = (long)LT;
        res->umbrella = w.on;
        res->Kmin = w.Kmin;
        res->Kmax = w.Kmax;
        if (w.on && w.Kmin > w.Kmax) {
            fprintf(stderr, "error: empty umbrella window (K from %ld to %ld)\n", w.Kmin, w.Kmax);
            traj_free(&tr);
            shot_free(&sh);
            return -1;
        }
        if (w.on) {
            for (m = 0; m < p->n_approach && window_dist(&w, tr.K) > 0; m++) {
                tps_move(&tr, &sh, &d, p->s, p->g, &w, &r, &dir);
                if (p->progress > 0 && (m + 1) % p->progress == 0)
                    fprintf(stderr, "approach %ld: K=%ld (window %ld..%ld)\n", m + 1, tr.K,
                            w.Kmin, w.Kmax);
            }
            res->n_approach_used = m;
            if (window_dist(&w, tr.K) > 0) {
                fprintf(stderr, "error: K=%ld did not reach the window %ld..%ld in %ld moves\n",
                        tr.K, w.Kmin, w.Kmax, p->n_approach);
                traj_free(&tr);
                shot_free(&sh);
                return -1;
            }
        }
    }

    for (m = 0; m < p->n_equil; m++) {
        tps_move(&tr, &sh, &d, p->s, p->g, &w, &r, &dir);
        if (p->progress > 0 && (m + 1) % p->progress == 0)
            fprintf(stderr, "equil %ld: K=%ld E=%ld\n", m + 1, tr.K, tr.E);
    }

    cap = p->n_prod / p->sample_every + 1;
    res->K = (double *)malloc(sizeof(double) * (size_t)cap);
    res->E = (double *)malloc(sizeof(double) * (size_t)cap);
    res->move = (long *)malloc(sizeof(long) * (size_t)cap);
    if (!res->K || !res->E || !res->move) {
        traj_free(&tr);
        shot_free(&sh);
        tps_result_free(res);
        return -1;
    }
    for (m = 0; m < p->n_prod; m++) {
        int acc = tps_move(&tr, &sh, &d, p->s, p->g, &w, &r, &dir);
        if (dir > 0) {
            res->n_fwd++;
            res->acc_fwd += acc;
        } else {
            res->n_bwd++;
            res->acc_bwd += acc;
        }
        if ((m + 1) % p->sample_every == 0) {
            res->K[res->nsamp] = (double)tr.K;
            res->E[res->nsamp] = (double)tr.E;
            res->move[res->nsamp] = m + 1;
            res->nsamp++;
        }
        if (p->progress > 0 && (m + 1) % p->progress == 0)
            fprintf(stderr, "prod %ld: K=%ld E=%ld\n", m + 1, tr.K, tr.E);
    }
    if (res->nsamp >= 2) pair_stats(res->K, res->E, res->nsamp, p->nblocks, &res->st);

    if (p->save && traj_save(&tr, p->save) != 0) {
        traj_free(&tr);
        shot_free(&sh);
        return -1;
    }
    traj_free(&tr);
    shot_free(&sh);
    return 0;
}

void tps_result_free(TPSResult *res)
{
    free(res->K);
    free(res->E);
    free(res->move);
    res->K = res->E = NULL;
    res->move = NULL;
    res->nsamp = 0;
}
