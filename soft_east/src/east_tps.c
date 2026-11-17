/*
 * east_tps: transition path sampling of the (s, g)-biased softened East model.
 *
 * usage: ./bin/east_tps [input_file]      (default input file: system.inp)
 */
#include "input.h"
#include "tps.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void write_summary(FILE *f, const char *pre, const TPSParams *p, const TPSResult *r)
{
    const PairStats *st = &r->st;
    const double nk = (double)p->L * p->tobs;        /* K normalization */
    const double ne = (double)p->L * (p->tobs + 1);  /* E normalization */
    fprintf(f, "%ssamples %ld, blocks %d\n", pre, r->nsamp, st->nblocks);
    fprintf(f, "%sacceptance forward %.4f backward %.4f total %.4f\n", pre,
            r->n_fwd ? (double)r->acc_fwd / r->n_fwd : 0.0,
            r->n_bwd ? (double)r->acc_bwd / r->n_bwd : 0.0,
            (double)(r->acc_fwd + r->acc_bwd) / (double)(r->n_fwd + r->n_bwd));
    fprintf(f, "%s<K>    = %.8e +- %.2e\n", pre, st->mean_x, st->err_mean_x);
    fprintf(f, "%s<E>    = %.8e +- %.2e\n", pre, st->mean_y, st->err_mean_y);
    fprintf(f, "%sVar K  = %.8e +- %.2e\n", pre, st->var_x, st->err_var_x);
    fprintf(f, "%sVar E  = %.8e +- %.2e\n", pre, st->var_y, st->err_var_y);
    fprintf(f, "%sCov KE = %.8e +- %.2e\n", pre, st->cov, st->err_cov);
    fprintf(f, "%sk = <K>/(L tobs)        = %.8e +- %.2e\n", pre, st->mean_x / nk,
            st->err_mean_x / nk);
    fprintf(f, "%se = <E>/(L (tobs+1))    = %.8e +- %.2e\n", pre, st->mean_y / ne,
            st->err_mean_y / ne);
    fprintf(f, "%schi_KK = Var K/(L tobs) = %.8e +- %.2e\n", pre, st->var_x / nk,
            st->err_var_x / nk);
    fprintf(f, "%schi_EE = Var E/(L tobs) = %.8e +- %.2e\n", pre, st->var_y / nk,
            st->err_var_y / nk);
    fprintf(f, "%schi_KE = Cov/(L tobs)   = %.8e +- %.2e\n", pre, st->cov / nk, st->err_cov / nk);
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "system.inp";
    Input in;
    TPSParams p;
    TPSResult r;
    const char *out = "tps_samples.dat";
    const char *v;
    FILE *f;
    long i;

    if (argc > 1 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        printf("usage: ./bin/east_tps [input_file]   (default: system.inp)\n"
               "keywords: L T eps s g tobs seed n_equil n_prod sample_every init save\n"
               "          output nblocks progress umb_kmin umb_kmax n_approach\n");
        return 0;
    }
    if (input_read(&in, path) != 0) return 1;
    tps_default_params(&p);
    if (input_int(&in, "L", 1, &p.L) < 0 || input_double(&in, "T", 0, &p.T) < 0 ||
        input_double(&in, "eps", 1, &p.eps) < 0 || input_double(&in, "s", 1, &p.s) < 0 ||
        input_double(&in, "g", 1, &p.g) < 0 || input_int(&in, "tobs", 1, &p.tobs) < 0 ||
        input_u64(&in, "seed", 0, &p.seed) < 0 || input_long(&in, "n_equil", 0, &p.n_equil) < 0 ||
        input_long(&in, "n_prod", 1, &p.n_prod) < 0 ||
        input_int(&in, "sample_every", 0, &p.sample_every) < 0 ||
        input_int(&in, "nblocks", 0, &p.nblocks) < 0 ||
        input_long(&in, "progress", 0, &p.progress) < 0 ||
        input_double(&in, "umb_kmin", 0, &p.umb_kmin) < 0 ||
        input_double(&in, "umb_kmax", 0, &p.umb_kmax) < 0 ||
        input_long(&in, "n_approach", 0, &p.n_approach) < 0) {
        input_free(&in);
        return 1;
    }
    if ((v = input_str(&in, "init"))) p.init = v;
    if ((v = input_str(&in, "save"))) p.save = v;
    if ((v = input_str(&in, "output"))) out = v;
    input_warn_unused(&in);

    if (tps_run(&p, &r) != 0) {
        input_free(&in);
        return 1;
    }

    f = fopen(out, "w");
    if (!f) {
        fprintf(stderr, "error: cannot write %s\n", out);
        tps_result_free(&r);
        input_free(&in);
        return 1;
    }
    fprintf(f, "# soft_east TPS samples\n");
    fprintf(f, "# param L %d\n# param T %.12g\n# param eps %.12g\n# param s %.12g\n"
               "# param g %.12g\n# param tobs %d\n# param seed %llu\n"
               "# param n_equil %ld\n# param n_prod %ld\n# param sample_every %d\n"
               "# param init %s\n",
            p.L, p.T, p.eps, p.s, p.g, p.tobs, (unsigned long long)r.seed_used, p.n_equil,
            p.n_prod, p.sample_every, p.init);
    if (r.umbrella) {
        fprintf(f, "# param umb_kmin %.12g\n# param umb_kmax %.12g\n", p.umb_kmin, p.umb_kmax);
        fprintf(f, "# param umb_Kmin %ld\n# param umb_Kmax %ld\n", r.Kmin, r.Kmax);
        fprintf(f, "# umbrella: K restricted to [%ld, %ld]; entered the window after %ld moves\n",
                r.Kmin, r.Kmax, r.n_approach_used);
        fprintf(f, "# averages in the summary are restricted to the window (combine with east_mbar)\n");
    }
    fprintf(f, "# K = flips in tobs MCS, E = sum of excitations over the tobs+1 slices\n");
    fprintf(f, "# columns: move K E\n");
    for (i = 0; i < r.nsamp; i++) fprintf(f, "%ld %.0f %.0f\n", r.move[i], r.K[i], r.E[i]);
    write_summary(f, "# summary ", &p, &r);
    fclose(f);

    printf("TPS L=%d T=%g eps=%g s=%g g=%g tobs=%d seed=%llu -> %s\n", p.L, p.T, p.eps, p.s, p.g,
           p.tobs, (unsigned long long)r.seed_used, out);
    if (r.umbrella)
        printf("  umbrella window K in [%ld, %ld] (k in [%.6g, %.6g]), entered after %ld moves;\n"
               "  averages below are restricted to the window\n",
               r.Kmin, r.Kmax, (double)r.Kmin / ((double)p.L * p.tobs),
               (double)r.Kmax / ((double)p.L * p.tobs), r.n_approach_used);
    write_summary(stdout, "  ", &p, &r);
    tps_result_free(&r);
    input_free(&in);
    return 0;
}
