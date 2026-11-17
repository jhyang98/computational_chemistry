/*
 * east_reus: replica-exchange umbrella sampling in the activity K.
 *
 * usage: ./bin/east_reus [input_file]      (default input file: system.inp)
 *
 * Writes one sample file per window (same format as east_tps, readable by
 * east_mbar), the final trajectory of each window, and a summary with swap
 * acceptance, round trips, overlap consistency and effective sample sizes.
 *
 * Checkpoints: with 'checkpoint <file>', the run can be stopped (stop_after
 * seconds, SIGTERM or SIGUSR1) and resumed by running the same input again.
 * Exit status: 0 finished, 3 stopped with a checkpoint (run again to resume),
 * 1 error. The checkpoint is removed once all output files are written.
 */
#include "input.h"
#include "reus.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void write_summary(FILE *f, const REUSParams *p, const REUSResult *r)
{
    const double LT = (double)p->L * p->tobs;
    long trips = 0;
    int w;
    fprintf(f, "# REUS L=%d T=%g eps=%g s=%g g=%g tobs=%d nwin=%d k=[%g, %g] seed=%llu\n", p->L,
            p->T, p->eps, p->s, p->g, p->tobs, p->nwin, p->k_lo, p->k_hi,
            (unsigned long long)r->seed_used);
    fprintf(f, "# moves per window: n_equil %ld, n_prod %ld; swap every %d, sample every %d\n",
            p->n_equil, p->n_prod, p->swap_every, p->sample_every);
    fprintf(f, "#\n# per window (production)\n");
    fprintf(f, "# file: 1 if started from a trajectory file; approach: moves to enter the window\n"
               "# (all attempts); from: window whose trajectory it restarted from (-1: none)\n");
    fprintf(f, "# %3s %7s %7s %9s %9s %10s %8s %9s %6s %9s %5s\n", "w", "Kmin", "Kmax", "kmin",
            "kmax", "<K>", "acc", "N_eff", "file", "approach", "from");
    for (w = 0; w < r->nwin; w++) {
        double m = 0.0;
        long i;
        for (i = 0; i < r->nsamp; i++) m += r->K[(size_t)w * r->nsamp + i];
        m /= (double)r->nsamp;
        fprintf(f, "  %3d %7ld %7ld %9.6f %9.6f %10.2f %8.4f %9.0f %6d %9ld %5d\n", w, r->Kmin[w],
                r->Kmax[w], r->Kmin[w] / LT, r->Kmax[w] / LT, m,
                r->moves[w] ? (double)r->acc[w] / r->moves[w] : 0.0, reus_neff(r, w),
                r->init_from_file[w], r->approach[w], r->seeded_from[w]);
    }
    fprintf(f, "#\n# neighbour pairs: swap acceptance and KS distance of K in the overlap\n");
    fprintf(f, "# (KS should be small, roughly below 1.4*sqrt(1/n1+1/n2) times sqrt(2 tau))\n");
    fprintf(f, "# %5s %16s %10s %8s %8s %8s\n", "pair", "overlap K", "swap_acc", "KS", "n1", "n2");
    for (w = 0; w + 1 < r->nwin; w++) {
        long n1, n2;
        double ks = reus_overlap_ks(r, w, &n1, &n2);
        long lo = r->Kmin[w + 1], hi = r->Kmax[w];
        fprintf(f, "  %02d-%02d %7ld..%-7ld %10.4f %8.3f %8ld %8ld\n", w, w + 1, lo, hi,
                r->swap_try[w] ? (double)r->swap_acc[w] / r->swap_try[w] : 0.0, ks, n1, n2);
    }
    for (w = 0; w < r->nwin; w++) trips += r->round_trips[w];
    fprintf(f, "#\n# round trips (window 0 -> %d -> 0) during production: %ld in total\n",
            r->nwin - 1, trips);
    fprintf(f, "# per replica:");
    for (w = 0; w < r->nwin; w++) fprintf(f, " %ld", r->round_trips[w]);
    fprintf(f, "\n");
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "system.inp";
    const char *prefix = "reus_", *v;
    Input in;
    REUSParams p;
    REUSResult r;
    int w, save_traj = 1, rc;
    char fname[1024];
    FILE *f;

    if (argc > 1 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        printf("usage: ./bin/east_reus [input_file]   (default: system.inp)\n"
               "keywords: L T eps s g tobs seed nwin k_lo k_hi n_approach n_equil n_prod\n"
               "          swap_every sample_every init progress output_prefix save_traj\n"
               "          checkpoint checkpoint_every stop_after\n"
               "exit status: 0 finished, 3 stopped with a checkpoint (run again to resume), "
               "1 error\n");
        return 0;
    }
    if (input_read(&in, path) != 0) return 1;
    reus_default_params(&p);
    if (input_int(&in, "L", 1, &p.L) < 0 || input_double(&in, "T", 0, &p.T) < 0 ||
        input_double(&in, "eps", 1, &p.eps) < 0 || input_double(&in, "s", 1, &p.s) < 0 ||
        input_double(&in, "g", 1, &p.g) < 0 || input_int(&in, "tobs", 1, &p.tobs) < 0 ||
        input_u64(&in, "seed", 0, &p.seed) < 0 || input_int(&in, "nwin", 1, &p.nwin) < 0 ||
        input_double(&in, "k_lo", 1, &p.k_lo) < 0 || input_double(&in, "k_hi", 1, &p.k_hi) < 0 ||
        input_long(&in, "n_approach", 0, &p.n_approach) < 0 ||
        input_long(&in, "n_equil", 0, &p.n_equil) < 0 ||
        input_long(&in, "n_prod", 1, &p.n_prod) < 0 ||
        input_int(&in, "swap_every", 0, &p.swap_every) < 0 ||
        input_int(&in, "sample_every", 0, &p.sample_every) < 0 ||
        input_long(&in, "progress", 0, &p.progress) < 0 ||
        input_int(&in, "save_traj", 0, &save_traj) < 0 ||
        input_double(&in, "checkpoint_every", 0, &p.checkpoint_every) < 0 ||
        input_double(&in, "stop_after", 0, &p.stop_after) < 0) {
        input_free(&in);
        return 1;
    }
    if ((v = input_str(&in, "init"))) p.init = v;
    if ((v = input_str(&in, "output_prefix"))) prefix = v;
    if ((v = input_str(&in, "checkpoint"))) p.checkpoint = v;
    input_warn_unused(&in);

#ifdef _OPENMP
    fprintf(stderr, "east_reus: OpenMP with up to %d threads\n", omp_get_max_threads());
#endif
    rc = reus_run(&p, &r);
    if (rc == REUS_STOPPED) {
        fprintf(stderr, "east_reus: stopped before the end; run the same input again to resume\n");
        input_free(&in);
        return 3;
    }
    if (rc != 0) {
        input_free(&in);
        return 1;
    }

    for (w = 0; w < r.nwin; w++) {
        const double LT = (double)p.L * p.tobs;
        long i;
        snprintf(fname, sizeof(fname), "%sw%02d.dat", prefix, w);
        f = fopen(fname, "w");
        if (!f) {
            fprintf(stderr, "error: cannot write %s\n", fname);
            reus_result_free(&r);
            input_free(&in);
            return 1;
        }
        fprintf(f, "# soft_east REUS samples, window %d of %d\n", w, r.nwin);
        fprintf(f, "# param L %d\n# param T %.12g\n# param eps %.12g\n# param s %.12g\n"
                   "# param g %.12g\n# param tobs %d\n# param seed %llu\n"
                   "# param n_equil %ld\n# param n_prod %ld\n# param sample_every %d\n"
                   "# param init %s\n",
                p.L, p.T, p.eps, p.s, p.g, p.tobs, (unsigned long long)r.seed_used, p.n_equil,
                p.n_prod, p.sample_every, p.init);
        fprintf(f, "# param umb_kmin %.12g\n# param umb_kmax %.12g\n", r.Kmin[w] / LT,
                r.Kmax[w] / LT);
        fprintf(f, "# param umb_Kmin %ld\n# param umb_Kmax %ld\n", r.Kmin[w], r.Kmax[w]);
        fprintf(f, "# param swap_every %d\n", p.swap_every);
        fprintf(f, "# K = flips in tobs MCS, E = sum of excitations over the tobs+1 slices\n");
        fprintf(f, "# columns: move K E replica\n");
        for (i = 0; i < r.nsamp; i++)
            fprintf(f, "%ld %.0f %.0f %d\n", (i + 1) * (long)p.sample_every,
                    r.K[(size_t)w * r.nsamp + i], r.E[(size_t)w * r.nsamp + i],
                    r.rep[(size_t)w * r.nsamp + i]);
        fclose(f);
        if (save_traj) {
            snprintf(fname, sizeof(fname), "%sw%02d.traj", prefix, w);
            traj_save(&r.final[w], fname);
        }
    }
    snprintf(fname, sizeof(fname), "%ssummary.dat", prefix);
    f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "error: cannot write %s\n", fname);
        reus_result_free(&r);
        input_free(&in);
        return 1;
    }
    write_summary(f, &p, &r);
    fclose(f);
    if (p.checkpoint && p.checkpoint[0]) remove(p.checkpoint); /* all outputs are written */
    write_summary(stdout, &p, &r);
    printf("# files: %sw00.dat ... %sw%02d.dat, %ssummary.dat\n", prefix, prefix, r.nwin - 1,
           prefix);
    reus_result_free(&r);
    input_free(&in);
    return 0;
}
