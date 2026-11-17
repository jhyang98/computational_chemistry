/*
 * east_clone: discrete-time cloning for the (s, g)-biased softened East model.
 *
 * usage: ./bin/east_clone [input_file]   (default input file: system.inp)
 */
#include "clone.h"
#include "input.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "system.inp";
    const char *out = NULL, *v;
    Input in;
    CloneParams p;
    double *ss = NULL, *gs = NULL, *th = NULL, *kr = NULL, *er = NULL;
    int ns = 0, ng = 0, is, ig, i, print_runs = 0;
    FILE *f = stdout;

    if (argc > 1 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        printf("usage: ./bin/east_clone [input_file]   (default: system.inp)\n"
               "keywords: L T eps s g n_clones tobs t_skip t_tail n_runs seed print_runs output\n"
               "          (s and g accept lists)\n");
        return 0;
    }
    if (input_read(&in, path) != 0) return 1;
    clone_default_params(&p);
    if (input_int(&in, "L", 1, &p.L) < 0 || input_double(&in, "T", 0, &p.T) < 0 ||
        input_double(&in, "eps", 1, &p.eps) < 0 || input_list(&in, "s", 1, &ss, &ns) < 0 ||
        input_list(&in, "g", 1, &gs, &ng) < 0 || input_int(&in, "n_clones", 1, &p.n_clones) < 0 ||
        input_int(&in, "tobs", 1, &p.tobs) < 0 || input_int(&in, "t_skip", 0, &p.t_skip) < 0 ||
        input_int(&in, "t_tail", 0, &p.t_tail) < 0 ||
        input_int(&in, "n_runs", 0, &p.n_runs) < 0 || input_u64(&in, "seed", 0, &p.seed) < 0 ||
        input_int(&in, "print_runs", 0, &print_runs) < 0) {
        input_free(&in);
        return 1;
    }
    if ((v = input_str(&in, "output"))) out = v;
    input_warn_unused(&in);
    if (out) {
        f = fopen(out, "w");
        if (!f) {
            fprintf(stderr, "error: cannot write %s\n", out);
            input_free(&in);
            return 1;
        }
    }
    th = (double *)malloc(sizeof(double) * (size_t)p.n_runs);
    kr = (double *)malloc(sizeof(double) * (size_t)p.n_runs);
    er = (double *)malloc(sizeof(double) * (size_t)p.n_runs);

    fprintf(f, "# soft_east discrete-time cloning\n");
    fprintf(f, "# L %d  T %.12g  eps %.12g  n_clones %d  tobs %d  t_skip %d  t_tail %d  n_runs %d\n",
            p.L, p.T, p.eps, p.n_clones, p.tobs, p.t_skip, p.t_tail, p.n_runs);
    fprintf(f, "# theta per MCS (whole system) from steps t_skip < t <= tobs\n");
    fprintf(f, "# k, e per site per MCS from clone histories over t_skip < t <= tobs - t_tail\n");
    fprintf(f, "# errors: standard error over independent runs\n");
    fprintf(f, "# %10s %12s %18s %10s %16s %10s %16s %10s %20s\n", "g", "s", "theta", "err", "k",
            "err", "e", "err", "seed");
    for (ig = 0; ig < ng; ig++) {
        for (is = 0; is < ns; is++) {
            CloneResult r;
            p.s = ss[is];
            p.g = gs[ig];
            if (clone_run(&p, &r, th, kr, er) != 0) {
                input_free(&in);
                return 1;
            }
            /* next point continues the random stream deterministically */
            if (p.seed != 0) p.seed = p.seed * 6364136223846793005ULL + 1442695040888963407ULL;
            fprintf(f, "%12.6g %12.8g %18.10e %10.2e %16.10e %10.2e %16.10e %10.2e %20llu\n", p.g,
                    p.s, r.theta, r.theta_err, r.k, r.k_err, r.e, r.e_err,
                    (unsigned long long)r.seed_used);
            if (print_runs)
                for (i = 0; i < p.n_runs; i++)
                    fprintf(f, "#   run %3d theta %.10e k %.10e e %.10e\n", i, th[i], kr[i], er[i]);
            fflush(f);
        }
    }
    if (f != stdout) fclose(f);
    free(th);
    free(kr);
    free(er);
    free(ss);
    free(gs);
    input_free(&in);
    return 0;
}
