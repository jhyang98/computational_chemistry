#include "reus.h"
#include "rng.h"
#include "tps.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    Traj tr;
    Shot sh;
    Rng rng;
    int win;  /* current window */
    int trip; /* 0: not yet at window 0, 1: left window 0, 2: reached the top since then */
} Replica;

void reus_default_params(REUSParams *p)
{
    memset(p, 0, sizeof(*p));
    p->L = 100;
    p->T = 1.0;
    p->tobs = 600;
    p->nwin = 22;
    p->k_lo = 0.0;
    p->k_hi = 0.08;
    p->n_approach = 1000000;
    p->n_equil = 20000;
    p->n_prod = 1000000;
    p->swap_every = 10;
    p->sample_every = 10;
    p->init = "empty";
    p->checkpoint = NULL;
    p->checkpoint_every = 1800.0;
    p->stop_after = 0.0;
    p->stop_cycle = 0;
}

void reus_window(const REUSParams *p, int w, long *Kmin, long *Kmax)
{
    const double LT = (double)p->L * p->tobs;
    const double D = (p->k_hi - p->k_lo) / (p->nwin + 1);
    long a = (long)ceil((p->k_lo + w * D) * LT - 1e-9);
    long b = (long)floor((p->k_lo + (w + 2) * D) * LT + 1e-9);
    if (a < 0) a = 0;
    if (b > (long)LT) b = (long)LT;
    *Kmin = a;
    *Kmax = b;
}

static int file_exists(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static int init_replica(const REUSParams *p, const Dyn *d, Replica *r, int i, int *from_file)
{
    *from_file = 0;
    if (traj_alloc(&r->tr, p->L, p->tobs) != 0 || shot_alloc(&r->sh, p->L, p->tobs) != 0)
        return -1;
    if (!strcmp(p->init, "empty")) {
        traj_recount(&r->tr);
    } else if (!strcmp(p->init, "equilibrium")) {
        traj_init_equilibrium(&r->tr, d, &r->rng);
    } else {
        char path[1024];
        const char *pc = strchr(p->init, '%');
        if (!pc) {
            snprintf(path, sizeof(path), "%s", p->init); /* one file for every window */
        } else {
            /* exactly one conversion of the form %d or %0Nd */
            const char *q = pc + 1;
            while (*q >= '0' && *q <= '9') q++;
            if (*q != 'd' || strchr(q + 1, '%')) {
                fprintf(stderr, "error: init pattern '%s' must contain one %%d (e.g. w%%02d.traj)\n",
                        p->init);
                return -1;
            }
            snprintf(path, sizeof(path), p->init, i);
        }
        if (file_exists(path) && traj_load(&r->tr, path) == 0) {
            *from_file = 1;
        } else {
            fprintf(stderr, "note: window %d: no usable trajectory '%s', starting empty\n", i,
                    path);
            memset(r->tr.x, 0, (size_t)(p->tobs + 1) * (size_t)p->L);
            memset(r->tr.k, 0, sizeof(int) * (size_t)(p->tobs + 1));
            traj_recount(&r->tr);
        }
    }
    return 0;
}

/* ---- stop requests and checkpoints ---- */

static volatile sig_atomic_t stop_requested = 0;

static void on_stop_signal(int sig)
{
    (void)sig;
    stop_requested = 1;
}

typedef struct {
    char magic[8];
    int L, tobs, nwin, swap_every, sample_every;
    double T, eps, s, g, k_lo, k_hi;
    long n_equil, n_prod, nsamp;
} CkHeader;

static void ck_header(const REUSParams *p, long nsamp, CkHeader *h)
{
    memset(h, 0, sizeof(*h)); /* also clears the padding written to the file */
    memcpy(h->magic, "SEREUS01", 8);
    h->L = p->L;
    h->tobs = p->tobs;
    h->nwin = p->nwin;
    h->swap_every = p->swap_every;
    h->sample_every = p->sample_every;
    h->T = p->T;
    h->eps = p->eps;
    h->s = p->s;
    h->g = p->g;
    h->k_lo = p->k_lo;
    h->k_hi = p->k_hi;
    h->n_equil = p->n_equil;
    h->n_prod = p->n_prod;
    h->nsamp = nsamp;
}

static int ck_io(FILE *f, void *buf, size_t size, size_t n, int writing)
{
    if (n == 0) return 0;
    return (writing ? fwrite(buf, size, n, f) : fread(buf, size, n, f)) == n ? 0 : -1;
}

/* State after c_next cycles: header, seed, master rng, replicas, counters, samples. */
static int ck_body(FILE *f, int writing, const CkHeader *h, REUSResult *res, Replica *rep,
                   int *rep_of_win, Rng *master, long *c_next, long *nfilled)
{
    const int W = h->nwin;
    const size_t xs = (size_t)(h->tobs + 1) * (size_t)h->L, ts = (size_t)(h->tobs + 1);
    int i, err = 0;
    err |= ck_io(f, &res->seed_used, sizeof(uint64_t), 1, writing);
    err |= ck_io(f, master, sizeof(Rng), 1, writing);
    err |= ck_io(f, c_next, sizeof(long), 1, writing);
    err |= ck_io(f, nfilled, sizeof(long), 1, writing);
    if (err || *nfilled < 0 || *nfilled > res->nsamp) return -1;
    for (i = 0; i < W; i++) {
        Replica *r = &rep[i];
        err |= ck_io(f, &r->rng, sizeof(Rng), 1, writing);
        err |= ck_io(f, &r->win, sizeof(int), 1, writing);
        err |= ck_io(f, &r->trip, sizeof(int), 1, writing);
        err |= ck_io(f, &r->tr.K, sizeof(long), 1, writing);
        err |= ck_io(f, &r->tr.E, sizeof(long), 1, writing);
        err |= ck_io(f, r->tr.x, 1, xs, writing);
        err |= ck_io(f, r->tr.n, sizeof(int), ts, writing);
        err |= ck_io(f, r->tr.k, sizeof(int), ts, writing);
    }
    err |= ck_io(f, rep_of_win, sizeof(int), (size_t)W, writing);
    err |= ck_io(f, res->moves, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->acc, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->swap_try, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->swap_acc, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->round_trips, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->approach, sizeof(long), (size_t)W, writing);
    err |= ck_io(f, res->init_from_file, sizeof(int), (size_t)W, writing);
    err |= ck_io(f, res->seeded_from, sizeof(int), (size_t)W, writing);
    for (i = 0; i < W; i++) {
        const size_t o = (size_t)i * (size_t)res->nsamp;
        err |= ck_io(f, res->K + o, sizeof(double), (size_t)*nfilled, writing);
        err |= ck_io(f, res->E + o, sizeof(double), (size_t)*nfilled, writing);
        err |= ck_io(f, res->rep + o, sizeof(int), (size_t)*nfilled, writing);
    }
    return err ? -1 : 0;
}

static int ck_write(const char *path, const REUSParams *p, REUSResult *res, Replica *rep,
                    int *rep_of_win, Rng *master, long c_next, long nfilled)
{
    char tmp[1100];
    CkHeader h;
    FILE *f;
    int err;
    ck_header(p, res->nsamp, &h);
    snprintf(tmp, sizeof(tmp), "%s.tmp", path);
    f = fopen(tmp, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot write checkpoint %s\n", tmp);
        return -1;
    }
    err = ck_io(f, &h, sizeof(h), 1, 1);
    if (!err) err = ck_body(f, 1, &h, res, rep, rep_of_win, master, &c_next, &nfilled);
    if (fclose(f) != 0) err = -1;
    if (err) {
        remove(tmp);
        fprintf(stderr, "error: writing checkpoint %s failed\n", tmp);
        return -1;
    }
#ifdef _WIN32
    remove(path); /* rename does not replace an existing file on Windows */
#endif
    if (rename(tmp, path) != 0) {
        fprintf(stderr, "error: cannot rename %s to %s\n", tmp, path);
        return -1;
    }
    return 0;
}

static int ck_read(const char *path, const REUSParams *p, REUSResult *res, Replica *rep,
                   int *rep_of_win, Rng *master, long *c_next)
{
    CkHeader want, got;
    FILE *f = fopen(path, "rb");
    long nfilled = 0;
    int err;
    if (!f) {
        fprintf(stderr, "error: cannot open checkpoint %s\n", path);
        return -1;
    }
    ck_header(p, res->nsamp, &want);
    memset(&got, 0, sizeof(got));
    err = ck_io(f, &got, sizeof(got), 1, 0);
    if (err || memcmp(got.magic, want.magic, 8) != 0) {
        fprintf(stderr, "error: %s is not a REUS checkpoint of this version\n", path);
        fclose(f);
        return -1;
    }
    if (got.L != want.L || got.tobs != want.tobs || got.nwin != want.nwin ||
        got.swap_every != want.swap_every || got.sample_every != want.sample_every ||
        got.T != want.T || got.eps != want.eps || got.s != want.s || got.g != want.g ||
        got.k_lo != want.k_lo || got.k_hi != want.k_hi || got.n_equil != want.n_equil ||
        got.n_prod != want.n_prod || got.nsamp != want.nsamp) {
        fprintf(stderr,
                "error: checkpoint %s was written with different parameters\n"
                "  checkpoint: L %d T %g eps %g s %g g %g tobs %d nwin %d k [%g, %g] "
                "n_equil %ld n_prod %ld swap_every %d sample_every %d\n"
                "  this run  : L %d T %g eps %g s %g g %g tobs %d nwin %d k [%g, %g] "
                "n_equil %ld n_prod %ld swap_every %d sample_every %d\n"
                "  (remove the checkpoint to start over)\n",
                path, got.L, got.T, got.eps, got.s, got.g, got.tobs, got.nwin, got.k_lo,
                got.k_hi, got.n_equil, got.n_prod, got.swap_every, got.sample_every, want.L,
                want.T, want.eps, want.s, want.g, want.tobs, want.nwin, want.k_lo, want.k_hi,
                want.n_equil, want.n_prod, want.swap_every, want.sample_every);
        fclose(f);
        return -1;
    }
    err = ck_body(f, 0, &got, res, rep, rep_of_win, master, c_next, &nfilled);
    fclose(f);
    if (err) {
        fprintf(stderr, "error: checkpoint %s is truncated or damaged\n", path);
        return -1;
    }
    return 0;
}

static void traj_copy(Traj *dst, const Traj *src)
{
    memcpy(dst->x, src->x, (size_t)(src->tobs + 1) * (size_t)src->L);
    memcpy(dst->n, src->n, sizeof(int) * (size_t)(src->tobs + 1));
    memcpy(dst->k, src->k, sizeof(int) * (size_t)(src->tobs + 1));
    dst->K = src->K;
    dst->E = src->E;
}

static long window_gap(long K, long Kmin, long Kmax)
{
    return K < Kmin ? Kmin - K : (K > Kmax ? K - Kmax : 0);
}

/* Approach phase of tps_shoot until K is in [Kmin, Kmax]; returns the moves used. */
static long approach(Replica *r, const Dyn *d, double s, double g, long Kmin, long Kmax,
                     long nmax)
{
    long m;
    int dir;
    for (m = 0; m < nmax && !stop_requested; m++) {
        if (r->tr.K >= Kmin && r->tr.K <= Kmax) break;
        tps_shoot(&r->tr, &r->sh, d, s, g, 1, Kmin, Kmax, &r->rng, &dir);
    }
    return m;
}

int reus_run(const REUSParams *p, REUSResult *res)
{
    const int W = p->nwin;
    Dyn d;
    Rng master;
    Replica *rep = NULL;
    int *rep_of_win = NULL;
    double *ws = NULL, *wg = NULL;
    long ncyc_eq, ncyc_pr, c, c0 = 0;
    int samp_cyc, i, fail = 0, stopped = 0;
    const int ck_on = p->checkpoint && p->checkpoint[0];
    int resume = 0;
    time_t t_start = time(NULL), t_ck;

    memset(res, 0, sizeof(*res));
    if (W < 2 || p->swap_every < 1 || p->sample_every < p->swap_every ||
        p->sample_every % p->swap_every != 0 || p->n_prod < p->sample_every ||
        !(p->k_hi > p->k_lo)) {
        fprintf(stderr, "error: need nwin >= 2, k_hi > k_lo, swap_every >= 1, sample_every a "
                        "multiple of swap_every, n_prod >= sample_every\n");
        return -1;
    }
    if (dyn_init(&d, p->L, p->T, p->eps) != 0) {
        fprintf(stderr, "error: invalid L, T or eps\n");
        return -1;
    }
    res->seed_used = rng_seed_auto(&master, p->seed);
    res->nwin = W;
    ncyc_eq = p->n_equil / p->swap_every;
    ncyc_pr = p->n_prod / p->swap_every;
    samp_cyc = p->sample_every / p->swap_every;
    res->nsamp = ncyc_pr / samp_cyc;

    res->Kmin = (long *)calloc((size_t)W, sizeof(long));
    res->Kmax = (long *)calloc((size_t)W, sizeof(long));
    res->K = (double *)malloc(sizeof(double) * (size_t)W * (size_t)res->nsamp);
    res->E = (double *)malloc(sizeof(double) * (size_t)W * (size_t)res->nsamp);
    res->rep = (int *)malloc(sizeof(int) * (size_t)W * (size_t)res->nsamp);
    res->moves = (long *)calloc((size_t)W, sizeof(long));
    res->acc = (long *)calloc((size_t)W, sizeof(long));
    res->swap_try = (long *)calloc((size_t)W, sizeof(long));
    res->swap_acc = (long *)calloc((size_t)W, sizeof(long));
    res->round_trips = (long *)calloc((size_t)W, sizeof(long));
    res->approach = (long *)calloc((size_t)W, sizeof(long));
    res->init_from_file = (int *)calloc((size_t)W, sizeof(int));
    res->seeded_from = (int *)malloc(sizeof(int) * (size_t)W);
    res->final = (Traj *)calloc((size_t)W, sizeof(Traj));
    rep = (Replica *)calloc((size_t)W, sizeof(Replica));
    rep_of_win = (int *)malloc(sizeof(int) * (size_t)W);
    ws = (double *)malloc(sizeof(double) * (size_t)W);
    wg = (double *)malloc(sizeof(double) * (size_t)W);
    if (!res->Kmin || !res->Kmax || !res->K || !res->E || !res->rep || !res->moves ||
        !res->acc || !res->swap_try || !res->swap_acc || !res->round_trips || !res->approach ||
        !res->init_from_file || !res->seeded_from || !res->final || !rep || !rep_of_win || !ws ||
        !wg) {
        fprintf(stderr, "error: out of memory\n");
        fail = 1;
        goto done;
    }
    for (i = 0; i < W; i++) {
        reus_window(p, i, &res->Kmin[i], &res->Kmax[i]);
        if (res->Kmin[i] > res->Kmax[i]) {
            fprintf(stderr, "error: window %d is empty\n", i);
            fail = 1;
            goto done;
        }
        ws[i] = p->s;
        wg[i] = p->g;
        res->seeded_from[i] = -1;
        rng_seed(&rep[i].rng, res->seed_used + 0x9E3779B97F4A7C15ULL * (uint64_t)(i + 1));
        rep[i].win = i;
        rep[i].trip = 0;
        rep_of_win[i] = i;
    }
    if (ck_on) {
        FILE *t = fopen(p->checkpoint, "rb");
        resume = (t != NULL);
        if (t) fclose(t);
        stop_requested = 0;
        signal(SIGTERM, on_stop_signal);
#ifdef SIGUSR1
        signal(SIGUSR1, on_stop_signal);
#endif
    }
    for (i = 0; i < W; i++) {
        int bad = resume ? (traj_alloc(&rep[i].tr, p->L, p->tobs) != 0 ||
                            shot_alloc(&rep[i].sh, p->L, p->tobs) != 0)
                         : (init_replica(p, &d, &rep[i], i, &res->init_from_file[i]) != 0);
        if (bad) {
            fprintf(stderr, "error: cannot initialize replica %d\n", i);
            fail = 1;
            goto done;
        }
    }
    if (resume) {
        if (ck_read(p->checkpoint, p, res, rep, rep_of_win, &master, &c0) != 0) {
            fail = 1;
            goto done;
        }
        fprintf(stderr, "east_reus: resumed from %s at cycle %ld of %ld\n", p->checkpoint, c0,
                ncyc_eq + ncyc_pr);
        goto cycles;
    }

    /* approach: every replica enters its own window */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (i = 0; i < W; i++)
        res->approach[i] = approach(&rep[i], &d, ws[i], wg[i], res->Kmin[i], res->Kmax[i],
                                    p->n_approach);

    /*
     * A replica that did not reach its window restarts from a copy of the
     * trajectory of a replica that did: the one whose K is closest to the
     * missing window (neighbouring windows overlap, so this is usually inside
     * it or close by), ties going to the nearest window index. Repeated until
     * every window is reached or three rounds in a row bring no new window in.
     */
    {
        int stall = 0;
        int *in = (int *)malloc(sizeof(int) * (size_t)W);
        int *src = (int *)malloc(sizeof(int) * (size_t)W);
        if (!in || !src) {
            free(in);
            free(src);
            fprintf(stderr, "error: out of memory\n");
            fail = 1;
            goto done;
        }
        for (;;) {
            int nmiss = 0, nin = 0, nnew = 0;
            if (stop_requested) break;
            for (i = 0; i < W; i++) {
                in[i] = rep[i].tr.K >= res->Kmin[i] && rep[i].tr.K <= res->Kmax[i];
                nin += in[i];
                nmiss += !in[i];
            }
            if (nmiss == 0 || nin == 0) break;
            for (i = 0; i < W; i++) {
                int j, best = -1;
                long bgap = 0;
                src[i] = -1;
                if (in[i]) continue;
                for (j = 0; j < W; j++) {
                    long gap;
                    if (!in[j]) continue;
                    gap = window_gap(rep[j].tr.K, res->Kmin[i], res->Kmax[i]);
                    if (best < 0 || gap < bgap || (gap == bgap && abs(j - i) < abs(best - i))) {
                        best = j;
                        bgap = gap;
                    }
                }
                src[i] = best;
                fprintf(stderr,
                        "note: window %d (K %ld..%ld) not reached (K = %ld); restarting from the "
                        "trajectory of window %d (K = %ld)\n",
                        i, res->Kmin[i], res->Kmax[i], rep[i].tr.K, best, rep[best].tr.K);
                traj_copy(&rep[i].tr, &rep[best].tr);
                res->seeded_from[i] = best;
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
            for (i = 0; i < W; i++)
                if (src[i] >= 0)
                    res->approach[i] += approach(&rep[i], &d, ws[i], wg[i], res->Kmin[i],
                                                 res->Kmax[i], p->n_approach);
            for (i = 0; i < W; i++)
                if (src[i] >= 0 && rep[i].tr.K >= res->Kmin[i] && rep[i].tr.K <= res->Kmax[i])
                    nnew++;
            stall = (nnew == 0) ? stall + 1 : 0;
            if (stall == 3) break; /* no progress */
        }
        free(in);
        free(src);
    }
    if (stop_requested) { /* asked to stop before the first checkpoint: start over next time */
        fprintf(stderr, "east_reus: stop requested during the approach phase; no checkpoint\n");
        stopped = 1;
        goto done;
    }
    for (i = 0; i < W; i++) {
        if (rep[i].tr.K < res->Kmin[i] || rep[i].tr.K > res->Kmax[i]) {
            fprintf(stderr, "error: window %d (K %ld..%ld) not reached (K = %ld)", i,
                    res->Kmin[i], res->Kmax[i], rep[i].tr.K);
            if (res->seeded_from[i] >= 0)
                fprintf(stderr, ", also not from the trajectory of window %d",
                        res->seeded_from[i]);
            fprintf(stderr, "; %ld moves per attempt\n", p->n_approach);
            fail = 1;
        }
    }
    if (fail) goto done;
    if (ck_on && ck_write(p->checkpoint, p, res, rep, rep_of_win, &master, 0, 0) != 0) {
        fail = 1;
        goto done;
    }

cycles:
    t_ck = time(NULL);
    for (c = c0; c < ncyc_eq + ncyc_pr; c++) {
        const int prod = (c >= ncyc_eq);
        int w;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (i = 0; i < W; i++) {
            Replica *r = &rep[i];
            const int wi = r->win;
            int m, dir;
            long acc = 0;
            for (m = 0; m < p->swap_every; m++)
                acc += tps_shoot(&r->tr, &r->sh, &d, ws[wi], wg[wi], 1, res->Kmin[wi],
                                 res->Kmax[wi], &r->rng, &dir);
            if (prod) { /* each replica owns a different window index */
                res->moves[wi] += p->swap_every;
                res->acc[wi] += acc;
            }
        }
        /* neighbour swaps */
        for (w = (int)(c % 2); w + 1 < W; w += 2) {
            const int a = rep_of_win[w], b = rep_of_win[w + 1];
            const double Ka = (double)rep[a].tr.K, Ea = (double)rep[a].tr.E;
            const double Kb = (double)rep[b].tr.K, Eb = (double)rep[b].tr.E;
            int ok = Ka >= res->Kmin[w + 1] && Ka <= res->Kmax[w + 1] && Kb >= res->Kmin[w] &&
                     Kb <= res->Kmax[w];
            if (ok) {
                double logacc = (ws[w] * Ka + wg[w] * Ea + ws[w + 1] * Kb + wg[w + 1] * Eb) -
                                (ws[w] * Kb + wg[w] * Eb + ws[w + 1] * Ka + wg[w + 1] * Ea);
                ok = logacc >= 0.0 || rng_uniform(&master) < exp(logacc);
            }
            if (prod) {
                res->swap_try[w]++;
                res->swap_acc[w] += ok;
            }
            if (ok) {
                rep_of_win[w] = b;
                rep_of_win[w + 1] = a;
                rep[a].win = w + 1;
                rep[b].win = w;
            }
        }
        /* round trips: window 0 -> top -> window 0 */
        for (i = 0; i < W; i++) {
            if (rep[i].win == 0) {
                if (rep[i].trip == 2 && prod) res->round_trips[i]++;
                rep[i].trip = 1;
            } else if (rep[i].win == W - 1 && rep[i].trip == 1) {
                rep[i].trip = 2;
            }
        }
        if (prod && (c - ncyc_eq + 1) % samp_cyc == 0) {
            long idx = (c - ncyc_eq + 1) / samp_cyc - 1;
            for (w = 0; w < W; w++) {
                const Replica *r = &rep[rep_of_win[w]];
                res->K[(size_t)w * res->nsamp + idx] = (double)r->tr.K;
                res->E[(size_t)w * res->nsamp + idx] = (double)r->tr.E;
                res->rep[(size_t)w * res->nsamp + idx] = rep_of_win[w];
            }
        }
        if (p->progress > 0 && (c + 1) % p->progress == 0) {
            long trips = 0;
            for (i = 0; i < W; i++) trips += res->round_trips[i];
            fprintf(stderr, "cycle %ld/%ld (%s): round trips so far %ld\n", c + 1,
                    ncyc_eq + ncyc_pr, prod ? "production" : "equilibration", trips);
        }
        if (ck_on && c + 1 < ncyc_eq + ncyc_pr) {
            const time_t now = time(NULL);
            const int stop = stop_requested || (p->stop_cycle > 0 && c + 1 >= p->stop_cycle) ||
                             (p->stop_after > 0.0 && difftime(now, t_start) >= p->stop_after);
            if (stop || (p->checkpoint_every > 0.0 &&
                         difftime(now, t_ck) >= p->checkpoint_every)) {
                const long nfilled = (c + 1 > ncyc_eq) ? (c + 1 - ncyc_eq) / samp_cyc : 0;
                if (ck_write(p->checkpoint, p, res, rep, rep_of_win, &master, c + 1, nfilled) !=
                    0) {
                    fail = 1;
                    goto done;
                }
                t_ck = now;
                if (stop) {
                    fprintf(stderr, "east_reus: checkpoint written at cycle %ld of %ld; stopping\n",
                            c + 1, ncyc_eq + ncyc_pr);
                    stopped = 1;
                    goto done;
                }
            }
        }
    }

    for (i = 0; i < W; i++) {
        const Traj *src = &rep[rep_of_win[i]].tr;
        if (traj_alloc(&res->final[i], p->L, p->tobs) != 0) {
            fail = 1;
            break;
        }
        memcpy(res->final[i].x, src->x, (size_t)(p->tobs + 1) * (size_t)p->L);
        memcpy(res->final[i].n, src->n, sizeof(int) * (size_t)(p->tobs + 1));
        memcpy(res->final[i].k, src->k, sizeof(int) * (size_t)(p->tobs + 1));
        res->final[i].K = src->K;
        res->final[i].E = src->E;
    }

done:
    if (rep) {
        for (i = 0; i < W; i++) {
            traj_free(&rep[i].tr);
            shot_free(&rep[i].sh);
        }
    }
    free(rep);
    free(rep_of_win);
    free(ws);
    free(wg);
    if (ck_on) {
        signal(SIGTERM, SIG_DFL);
#ifdef SIGUSR1
        signal(SIGUSR1, SIG_DFL);
#endif
    }
    if (fail || stopped) {
        reus_result_free(res);
        return fail ? -1 : REUS_STOPPED;
    }
    return 0;
}

void reus_result_free(REUSResult *res)
{
    int i;
    if (res->final)
        for (i = 0; i < res->nwin; i++) traj_free(&res->final[i]);
    free(res->final);
    free(res->Kmin);
    free(res->Kmax);
    free(res->K);
    free(res->E);
    free(res->rep);
    free(res->moves);
    free(res->acc);
    free(res->swap_try);
    free(res->swap_acc);
    free(res->round_trips);
    free(res->approach);
    free(res->init_from_file);
    free(res->seeded_from);
    memset(res, 0, sizeof(*res));
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

double reus_overlap_ks(const REUSResult *res, int w, long *n1, long *n2)
{
    const long lo = res->Kmin[w] > res->Kmin[w + 1] ? res->Kmin[w] : res->Kmin[w + 1];
    const long hi = res->Kmax[w] < res->Kmax[w + 1] ? res->Kmax[w] : res->Kmax[w + 1];
    const long n = res->nsamp;
    double *a = (double *)malloc(sizeof(double) * (size_t)n);
    double *b = (double *)malloc(sizeof(double) * (size_t)n);
    long na = 0, nb = 0, i = 0, j = 0;
    double dmax = 0.0;
    for (i = 0; i < n; i++) {
        double x = res->K[(size_t)w * n + i], y = res->K[(size_t)(w + 1) * n + i];
        if (x >= lo && x <= hi) a[na++] = x;
        if (y >= lo && y <= hi) b[nb++] = y;
    }
    *n1 = na;
    *n2 = nb;
    if (na == 0 || nb == 0) {
        free(a);
        free(b);
        return NAN;
    }
    qsort(a, (size_t)na, sizeof(double), cmp_double);
    qsort(b, (size_t)nb, sizeof(double), cmp_double);
    i = j = 0;
    while (i < na && j < nb) {
        double v = a[i] < b[j] ? a[i] : b[j], dd;
        while (i < na && a[i] <= v) i++;
        while (j < nb && b[j] <= v) j++;
        dd = fabs((double)i / na - (double)j / nb);
        if (dd > dmax) dmax = dd;
    }
    free(a);
    free(b);
    return dmax;
}

double reus_neff(const REUSResult *res, int w)
{
    const long n = res->nsamp;
    const int nb = 20;
    const long B = n / nb;
    const double *x = res->K + (size_t)w * n;
    double m = 0.0, v = 0.0, vb = 0.0, neff;
    long i;
    int b;
    if (B < 1) return NAN;
    for (i = 0; i < nb * B; i++) m += x[i];
    m /= (double)(nb * B);
    for (i = 0; i < nb * B; i++) v += (x[i] - m) * (x[i] - m);
    v /= (double)(nb * B - 1);
    for (b = 0; b < nb; b++) {
        double mb = 0.0;
        for (i = 0; i < B; i++) mb += x[b * B + i];
        mb /= (double)B;
        vb += (mb - m) * (mb - m);
    }
    vb /= (nb - 1);
    if (!(vb > 0.0) || !(v > 0.0)) return (double)n;
    neff = (double)n * v / ((double)B * vb);
    return neff < (double)n ? neff : (double)n;
}
