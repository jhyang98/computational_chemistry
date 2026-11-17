/*
 * east_ed: exact diagonalization of the continuous-time tilted generator of
 * the softened East model with activity bias s and energy bias g.
 *
 * Modes:
 *   scan    observables on an (s, g) grid
 *   peak    pseudo-critical point s*(g) = argmax_s chi_KK (eps > 0)
 *   strict  strict East coexistence point theta_nonempty(s*, g) = 0 (eps = 0)
 */
#include "cli_util.h"
#include "ed.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EAST_ED_VERSION "0.1.0"

typedef struct {
    const char *mode;
    double *Ls;   int nL;
    double *eps;  int neps;
    double *gs;   int ng;
    double *ss;   int ns;
    double T;
    int use_translation;
    ChiMode chi;
    double fd_h;
    double stol;
    double side;   /* > 0: report k, e at s* -/+ side in peak mode */
    double tol;
    int mmax;
    double s0;
    int have_s0;
    const char *outpath;
    int verbose;
} Args;

static void usage(void)
{
    fprintf(stderr,
        "east_ed %s - exact diagonalization of the (s,g)-biased softened East model\n"
        "\n"
        "usage: ./bin/east_ed <mode> [options]\n"
        "\n"
        "modes:\n"
        "  scan     observables on a grid of s and g\n"
        "  peak     s*(g) = argmax_s chi_KK at each g (requires eps > 0)\n"
        "  strict   strict East (eps = 0): s*(g) where theta_nonempty = 0\n"
        "\n"
        "lists accept 'a', 'a,b,c' or 'a:step:b'; numbers accept 'p/q'\n"
        "\n"
        "options:\n"
        "  -L <list>        system sizes (default 10)\n"
        "  --T <val>        temperature, J = 1 (default 1)\n"
        "  --eps <list>     softness values (default 0.01; strict mode uses 0)\n"
        "  --g <list>       energy bias values (default 0)\n"
        "  --s <list>       activity bias values; coarse grid in peak mode\n"
        "                   (default -0.05:0.005:0.10)\n"
        "  --basis k0|full  zero-momentum sector (default) or full space\n"
        "  --chi lr|fd|none second derivatives: linear response (default),\n"
        "                   finite differences, or none\n"
        "  --fd-h <val>     finite-difference step (default 1e-4)\n"
        "  --stol <val>     tolerance on s* (default 1e-7)\n"
        "  --side <val>     peak mode: also report k, e at s* -/+ val\n"
        "  --s0 <val>       strict mode: initial guess for s* (default 0)\n"
        "  --tol <val>      Lanczos residual tolerance (default 1e-11)\n"
        "  --mmax <int>     Lanczos basis size (default 100)\n"
        "  --out <file>     write results to file (default stdout)\n"
        "  -v               progress messages on stderr\n",
        EAST_ED_VERSION);
}

static int need_value(int i, int argc, const char *opt)
{
    if (i + 1 >= argc) {
        fprintf(stderr, "error: option %s needs a value\n", opt);
        return -1;
    }
    return 0;
}

static int parse_args(int argc, char **argv, Args *a)
{
    int i;
    memset(a, 0, sizeof(*a));
    if (argc < 2) return -1;
    a->mode = argv[1];
    a->T = 1.0;
    a->use_translation = 1;
    a->chi = CHI_LR;
    a->fd_h = 1e-4;
    a->stol = 1e-7;
    a->side = 0.0;
    a->tol = 1e-11;
    a->mmax = 100;
    a->s0 = 0.0;

    for (i = 2; i < argc; i++) {
        const char *o = argv[i];
        if (!strcmp(o, "-v")) {
            a->verbose = 1;
            continue;
        }
        if (need_value(i, argc, o)) return -1;
        if (!strcmp(o, "-L")) {
            if (parse_list(argv[++i], &a->Ls, &a->nL)) goto bad;
        } else if (!strcmp(o, "--T")) {
            if (parse_number(argv[++i], &a->T)) goto bad;
        } else if (!strcmp(o, "--eps")) {
            if (parse_list(argv[++i], &a->eps, &a->neps)) goto bad;
        } else if (!strcmp(o, "--g")) {
            if (parse_list(argv[++i], &a->gs, &a->ng)) goto bad;
        } else if (!strcmp(o, "--s")) {
            if (parse_list(argv[++i], &a->ss, &a->ns)) goto bad;
        } else if (!strcmp(o, "--basis")) {
            i++;
            if (!strcmp(argv[i], "k0")) a->use_translation = 1;
            else if (!strcmp(argv[i], "full")) a->use_translation = 0;
            else goto bad;
        } else if (!strcmp(o, "--chi")) {
            i++;
            if (!strcmp(argv[i], "lr")) a->chi = CHI_LR;
            else if (!strcmp(argv[i], "fd")) a->chi = CHI_FD;
            else if (!strcmp(argv[i], "none")) a->chi = CHI_NONE;
            else goto bad;
        } else if (!strcmp(o, "--fd-h")) {
            if (parse_number(argv[++i], &a->fd_h)) goto bad;
        } else if (!strcmp(o, "--stol")) {
            if (parse_number(argv[++i], &a->stol)) goto bad;
        } else if (!strcmp(o, "--side")) {
            if (parse_number(argv[++i], &a->side)) goto bad;
        } else if (!strcmp(o, "--s0")) {
            if (parse_number(argv[++i], &a->s0)) goto bad;
            a->have_s0 = 1;
        } else if (!strcmp(o, "--tol")) {
            if (parse_number(argv[++i], &a->tol)) goto bad;
        } else if (!strcmp(o, "--mmax")) {
            double m;
            if (parse_number(argv[++i], &m)) goto bad;
            a->mmax = (int)m;
        } else if (!strcmp(o, "--out")) {
            a->outpath = argv[++i];
        } else {
            fprintf(stderr, "error: unknown option %s\n", o);
            return -1;
        }
        continue;
    bad:
        fprintf(stderr, "error: bad value for %s: %s\n", o, argv[i]);
        return -1;
    }

    if (!a->Ls) {
        a->Ls = (double *)malloc(sizeof(double));
        a->Ls[0] = 10;
        a->nL = 1;
    }
    if (!a->eps) {
        a->eps = (double *)malloc(sizeof(double));
        a->eps[0] = strcmp(a->mode, "strict") ? 0.01 : 0.0;
        a->neps = 1;
    }
    if (!a->gs) {
        a->gs = (double *)malloc(sizeof(double));
        a->gs[0] = 0.0;
        a->ng = 1;
    }
    if (!a->ss) {
        if (parse_list("-0.05:0.005:0.10", &a->ss, &a->ns)) return -1;
    }
    return 0;
}

static void print_header(FILE *f, int argc, char **argv, const Args *a)
{
    int i;
    fprintf(f, "# east_ed %s\n# command:", EAST_ED_VERSION);
    for (i = 0; i < argc; i++) fprintf(f, " %s", argv[i]);
    fprintf(f, "\n# T = %.12g, c = %.12g, basis = %s\n", a->T, 1.0 / (1.0 + exp(1.0 / a->T)),
            a->use_translation ? "k0 (zero-momentum sector)" : "full");
    fprintf(f, "# weight exp(-s K - g E); theta = largest eigenvalue (whole system, per unit time)\n");
    fprintf(f, "# k = -(1/L) dtheta/ds, e = -(1/L) dtheta/dg, chi_XY = (1/L) d2theta/dXdY\n");
}

static int init_solver(EDSolver *ed, int L, double T, double eps, const Args *a)
{
    if (ed_init(ed, L, T, eps, a->use_translation) != 0) {
        fprintf(stderr, "error: cannot initialize L=%d T=%g eps=%g\n", L, T, eps);
        return -1;
    }
    ed->lopt.tol = a->tol;
    ed->lopt.mmax = a->mmax;
    ed->fd_h = a->fd_h;
    return 0;
}

static void warn_point(const EDPoint *p, int L, double eps)
{
    if (!p->ok)
        fprintf(stderr, "warning: L=%d eps=%g g=%g s=%g: solver tolerance not met "
                        "(resid=%.2e, cg_relres=%.2e)\n",
                L, eps, p->g, p->s, p->resid, p->cg_relres);
    if (p->min_ratio < -1e-8)
        fprintf(stderr, "warning: L=%d eps=%g g=%g s=%g: top eigenvector not positive "
                        "(min/max=%.2e)\n",
                L, eps, p->g, p->s, p->min_ratio);
}

static int run_scan(FILE *f, const Args *a)
{
    int il, ie, ig, is;
    fprintf(f, "# %3s %10s %12s %12s %20s %20s %14s %18s %18s %18s %18s %18s %12s %9s %7s %7s %3s\n",
            "L", "eps", "g", "s", "theta", "theta2", "gap", "k", "e",
            "chi_KK", "chi_KE", "chi_EE", "rho", "resid", "nmv", "cg_it", "ok");
    for (il = 0; il < a->nL; il++) {
        for (ie = 0; ie < a->neps; ie++) {
            EDSolver ed;
            int L = (int)a->Ls[il];
            if (init_solver(&ed, L, a->T, a->eps[ie], a)) return -1;
            if (a->verbose)
                fprintf(stderr, "L=%d eps=%g: basis dim %ld\n", L, a->eps[ie], ed.basis.dim);
            for (ig = 0; ig < a->ng; ig++) {
                for (is = 0; is < a->ns; is++) {
                    EDPoint p;
                    if (ed_point(&ed, a->ss[is], a->gs[ig], a->chi, &p) != 0) {
                        fprintf(stderr, "error: solver failed\n");
                        ed_free(&ed);
                        return -1;
                    }
                    warn_point(&p, L, a->eps[ie]);
                    fprintf(f, "%5d %10.6g %12.6g %12.8g %20.12e %20.12e %14.6e %18.10e %18.10e "
                               "%18.10e %18.10e %18.10e %12.8f %9.1e %7d %7d %3d\n",
                            L, a->eps[ie], p.g, p.s, p.theta, p.theta2, p.theta - p.theta2, p.k,
                            p.e, p.chi_kk, p.chi_ke, p.chi_ee, p.rho, p.resid, p.nmatvec,
                            p.cg_iters, p.ok);
                }
                fprintf(f, "\n");
                fflush(f);
            }
            ed_free(&ed);
        }
    }
    return 0;
}

static int run_peak(FILE *f, const Args *a)
{
    int il, ie, ig;
    fprintf(f, "# edge = 1: the largest chi_KK on the s grid is at its boundary; s* is then\n"
               "#           not defined and the row is nan (the result depends on the s range)\n");
    if (a->side > 0.0)
        fprintf(f, "# active = s* - %g, inactive = s* + %g; dk = k_active - k_inactive, "
                   "de = e_active - e_inactive\n",
                a->side, a->side);
    fprintf(f, "# %3s %10s %12s %16s %18s %18s %14s %16s %16s %16s %16s %12s %4s %3s",
            "L", "eps", "g", "s_star", "chi_KK_max", "theta", "gap", "k", "e",
            "chi_KE", "chi_EE", "rho", "edge", "ok");
    if (a->side > 0.0)
        fprintf(f, " %14s %14s %14s %14s %14s %14s", "k_active", "e_active", "k_inactive",
                "e_inactive", "dk", "de");
    fprintf(f, "\n");

    for (il = 0; il < a->nL; il++) {
        for (ie = 0; ie < a->neps; ie++) {
            EDSolver ed;
            int L = (int)a->Ls[il];
            if (a->eps[ie] <= 0.0) {
                fprintf(stderr, "error: peak mode needs eps > 0 (use strict mode for eps = 0)\n");
                return -1;
            }
            if (init_solver(&ed, L, a->T, a->eps[ie], a)) return -1;
            for (ig = 0; ig < a->ng; ig++) {
                EDPoint p;
                int edge = 0;
                ed_reset_guess(&ed);
                if (ed_find_peak(&ed, a->gs[ig], a->ss, a->ns, a->stol, &p, &edge) != 0) {
                    fprintf(stderr, "error: peak search failed\n");
                    ed_free(&ed);
                    return -1;
                }
                warn_point(&p, L, a->eps[ie]);
                if (edge) {
                    /* no interior maximum of chi_KK on the s grid: s* is not defined */
                    fprintf(stderr, "note: L=%d eps=%g g=%g: chi_KK maximum at the edge of the "
                                    "s grid; s* reported as nan\n",
                            L, a->eps[ie], a->gs[ig]);
                    fprintf(f, "%5d %10.6g %12.6g %16s %18s %18s %14s %16s %16s %16s %16s %12s "
                               "%4d %3d",
                            L, a->eps[ie], a->gs[ig], "nan", "nan", "nan", "nan", "nan", "nan",
                            "nan", "nan", "nan", edge, p.ok);
                    if (a->side > 0.0)
                        fprintf(f, " %14s %14s %14s %14s %14s %14s", "nan", "nan", "nan", "nan",
                                "nan", "nan");
                    fprintf(f, "\n");
                    fflush(f);
                    continue;
                }
                if (a->verbose)
                    fprintf(stderr, "L=%d eps=%g g=%g: s*=%.8f chi=%.6e\n", L, a->eps[ie],
                            a->gs[ig], p.s, p.chi_kk);
                fprintf(f, "%5d %10.6g %12.6g %16.10f %18.10e %18.10e %14.6e %16.10e %16.10e "
                           "%16.10e %16.10e %12.8f %4d %3d",
                        L, a->eps[ie], p.g, p.s, p.chi_kk, p.theta, p.theta - p.theta2, p.k, p.e,
                        p.chi_ke, p.chi_ee, p.rho, edge, p.ok);
                if (a->side > 0.0) {
                    EDPoint pa, pi;
                    double sstar = p.s;
                    if (ed_point(&ed, sstar - a->side, a->gs[ig], CHI_NONE, &pa) != 0 ||
                        ed_point(&ed, sstar + a->side, a->gs[ig], CHI_NONE, &pi) != 0) {
                        ed_free(&ed);
                        return -1;
                    }
                    fprintf(f, " %14.8e %14.8e %14.8e %14.8e %14.8e %14.8e", pa.k, pa.e, pi.k,
                            pi.e, pa.k - pi.k, pa.e - pi.e);
                }
                fprintf(f, "\n");
                fflush(f);
            }
            ed_free(&ed);
        }
    }
    return 0;
}

static int run_strict(FILE *f, const Args *a)
{
    int il, ie, ig;
    for (ie = 0; ie < a->neps; ie++) {
        if (a->eps[ie] != 0.0) {
            fprintf(stderr, "error: strict mode is defined for eps = 0 only\n");
            return -1;
        }
    }
    fprintf(f, "# coexistence with the empty state (eigenvalue 0): theta_nonempty(s*, g) = 0\n");
    fprintf(f, "# active branch at s*: k, e; inactive (empty) branch has k = e = 0\n");
    fprintf(f, "# slope = -e/k = ds*/dg (Clausius-Clapeyron-like relation with k = e = 0 inactive)\n");
    fprintf(f, "# s_exact: product-state solution, exp(-2s*) = [(1+2g)^2 - (1-2c)^2] / [4c(1-c)]\n");
    fprintf(f, "# %3s %12s %16s %16s %16s %16s %14s %14s %9s %3s\n", "L", "g", "s_star", "s_exact",
            "k_active", "e_active", "slope", "gap_nonempty", "resid", "ok");
    for (il = 0; il < a->nL; il++) {
        EDSolver ed;
        int L = (int)a->Ls[il];
        double guess = a->s0;
        if (init_solver(&ed, L, a->T, 0.0, a)) return -1;
        for (ig = 0; ig < a->ng; ig++) {
            EDPoint p;
            const double c = ed.model.c, g = a->gs[ig];
            double num = (1.0 + 2.0 * g) * (1.0 + 2.0 * g) - (1.0 - 2.0 * c) * (1.0 - 2.0 * c);
            double s_exact = (num > 0.0) ? -0.5 * log(num / (4.0 * c * (1.0 - c))) : NAN;
            if (ed_find_strict_crossing(&ed, g, guess, a->stol, &p) != 0) {
                fprintf(stderr, "error: strict crossing not found for L=%d g=%g\n", L, g);
                ed_free(&ed);
                return -1;
            }
            warn_point(&p, L, 0.0);
            fprintf(f, "%5d %12.6g %16.10f %16.10f %16.10e %16.10e %14.8f %14.6e %9.1e %3d\n", L,
                    p.g, p.s, s_exact, p.k, p.e, -p.e / p.k, p.theta - p.theta2, p.resid, p.ok);
            fflush(f);
            guess = p.s;
        }
        ed_free(&ed);
    }
    return 0;
}

int main(int argc, char **argv)
{
    Args a;
    FILE *f = stdout;
    int rc;

    if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        usage();
        return argc < 2 ? 1 : 0;
    }
    if (parse_args(argc, argv, &a) != 0) {
        usage();
        return 1;
    }
    if (a.outpath) {
        f = fopen(a.outpath, "w");
        if (!f) {
            fprintf(stderr, "error: cannot open %s\n", a.outpath);
            return 1;
        }
    }
    print_header(f, argc, argv, &a);

    if (!strcmp(a.mode, "scan") || !strcmp(a.mode, "point")) rc = run_scan(f, &a);
    else if (!strcmp(a.mode, "peak")) rc = run_peak(f, &a);
    else if (!strcmp(a.mode, "strict")) rc = run_strict(f, &a);
    else {
        fprintf(stderr, "error: unknown mode %s\n", a.mode);
        usage();
        rc = -1;
    }

    if (f != stdout) fclose(f);
    free(a.Ls);
    free(a.eps);
    free(a.gs);
    free(a.ss);
    return rc == 0 ? 0 : 1;
}
