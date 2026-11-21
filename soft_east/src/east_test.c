/*
 * Validation suite for east_ed.
 *
 *  1. sparse coupling matrix is symmetric; zero-momentum sector dimensions
 *  2. s = g = 0: theta = 0, k and e equal their equilibrium values
 *  3. eps = 1 (non-interacting Glauber spins): closed-form theta and derivatives
 *  4. zero-momentum sector agrees with the full configuration space
 *  5. Lanczos + CG agree with dense diagonalization + sum over states
 *  6. Hellmann-Feynman and linear-response derivatives agree with finite
 *     differences; Maxwell relation dk/dg = de/ds = -chi_KE
 *  7. strict East: s*(g) equals the closed-form, L-independent line
 *  8. top eigenvector is positive (Perron-Frobenius)
 *  9. softened East: exact product-state eigenvector on the line
 *     mu(s) = g / (1 - eps)
 * 10. discrete-time (MCS) transfer matrix: brute-force enumeration of all
 *     attempt sequences, equilibrium values, derivatives vs finite
 *     differences, k0 vs full basis, long-time limit, eps = 1 closed form,
 *     one-sweep kernel normalization
 * 11. REUS start-up: windows that the approach phase cannot reach are seeded
 *     from the trajectory of a reached window (fixed seed, small system)
 * 12. REUS checkpoints: a run stopped twice and resumed from its checkpoint
 *     gives exactly the same samples and counters as an uninterrupted run;
 *     a checkpoint with other parameters is refused
 */
#include "ed.h"
#include "eig_dense.h"
#include "reus.h"
#include "tm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n_pass = 0, n_fail = 0;

static void check(const char *name, double got, double expect, double tol)
{
    double err = fabs(got - expect);
    double scale = fabs(expect) > 1.0 ? fabs(expect) : 1.0;
    int ok = (err <= tol * scale) && !isnan(got);
    if (ok) n_pass++;
    else n_fail++;
    printf("[%s] %-58s got % .12e  expect % .12e  err %.1e\n", ok ? "PASS" : "FAIL", name, got,
           expect, err);
}

static void check_true(const char *name, int cond)
{
    if (cond) n_pass++;
    else n_fail++;
    printf("[%s] %s\n", cond ? "PASS" : "FAIL", name);
}

/* ---------- closed form for eps = 1 (per-site two-state process) ---------- */

static double glauber_lambda(double c, double s, double g)
{
    double tau = -(1.0 + g);
    double det = c * (1.0 - c + g) - exp(-2.0 * s) * c * (1.0 - c);
    return 0.5 * (tau + sqrt(tau * tau - 4.0 * det));
}

/* derivatives of lambda by Richardson-extrapolated central differences */
static double d1(double c, double s, double g, int ds, double h)
{
    double fp = ds ? glauber_lambda(c, s + h, g) : glauber_lambda(c, s, g + h);
    double fm = ds ? glauber_lambda(c, s - h, g) : glauber_lambda(c, s, g - h);
    return (fp - fm) / (2.0 * h);
}

static double d1r(double c, double s, double g, int ds)
{
    double h = 1e-3;
    return (4.0 * d1(c, s, g, ds, h / 2) - d1(c, s, g, ds, h)) / 3.0;
}

static double d2(double c, double s, double g, int which, double h)
{
    /* which: 0 = ss, 1 = gg, 2 = sg */
    double f0 = glauber_lambda(c, s, g);
    if (which == 0)
        return (glauber_lambda(c, s + h, g) - 2.0 * f0 + glauber_lambda(c, s - h, g)) / (h * h);
    if (which == 1)
        return (glauber_lambda(c, s, g + h) - 2.0 * f0 + glauber_lambda(c, s, g - h)) / (h * h);
    return (glauber_lambda(c, s + h, g + h) - glauber_lambda(c, s + h, g - h) -
            glauber_lambda(c, s - h, g + h) + glauber_lambda(c, s - h, g - h)) /
           (4.0 * h * h);
}

static double d2r(double c, double s, double g, int which)
{
    double h = 2e-3;
    return (4.0 * d2(c, s, g, which, h / 2) - d2(c, s, g, which, h)) / 3.0;
}

/* ---------- dense reference ---------- */

typedef struct {
    double theta, theta2, k, e, chi_kk, chi_ke, chi_ee;
} DenseRef;

static int dense_reference(const EDSolver *ed, double s, double g, DenseRef *out)
{
    const TGen *t = &ed->tg;
    const int n = (int)t->dim;
    const int L = ed->model.L;
    const double es = exp(-s);
    double *H = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    double *Aden = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    double *w = (double *)malloc(sizeof(double) * (size_t)n);
    double *z = (double *)malloc(sizeof(double) * (size_t)n * (size_t)n);
    double *Av = (double *)malloc(sizeof(double) * (size_t)n);
    double vAv = 0.0, sum_ss = 0.0, sum_gg = 0.0, sum_sg = 0.0, e = 0.0;
    int i, j, m;
    long p;

    if (!H || !Aden || !w || !z || !Av) return -1;
    for (i = 0; i < n; i++) {
        for (p = t->rowptr[i]; p < t->rowptr[i + 1]; p++) {
            Aden[i * n + t->col[p]] += t->val[p];
            H[i * n + t->col[p]] += es * t->val[p];
        }
        H[i * n + i] -= t->esc[i] + g * t->occ[i];
    }
    if (eig_sym_jacobi(n, H, w, z) < 0) return -1;

    /* top eigenvector v = z[:,0] */
    for (i = 0; i < n; i++) {
        double acc = 0.0;
        for (j = 0; j < n; j++) acc += Aden[i * n + j] * z[j * n + 0];
        Av[i] = acc;
        vAv += z[i * n + 0] * acc;
        e += t->occ[i] * z[i * n + 0] * z[i * n + 0];
    }
    for (m = 1; m < n; m++) {
        double ms = 0.0, mg = 0.0; /* <m|H_s|0>, <m|H_g|0> */
        double dth = w[0] - w[m];
        for (i = 0; i < n; i++) {
            ms += z[i * n + m] * (-es * Av[i]);
            mg += z[i * n + m] * (-t->occ[i] * z[i * n + 0]);
        }
        sum_ss += ms * ms / dth;
        sum_gg += mg * mg / dth;
        sum_sg += ms * mg / dth;
    }
    out->theta = w[0];
    out->theta2 = (n > 1) ? w[1] : NAN;
    out->k = es * vAv / L;
    out->e = e / L;
    out->chi_kk = (es * vAv + 2.0 * sum_ss) / L;
    out->chi_ee = 2.0 * sum_gg / L;
    out->chi_ke = 2.0 * sum_sg / L;
    free(H);
    free(Aden);
    free(w);
    free(z);
    free(Av);
    return 0;
}

/* ---------- tests ---------- */

static long necklaces(int L)
{
    /* (1/L) sum_{d | L} phi(d) 2^{L/d} */
    long total = 0;
    int d;
    for (d = 1; d <= L; d++) {
        if (L % d == 0) {
            int phi = 0, k;
            for (k = 1; k <= d; k++) {
                int a = k, b = d;
                while (b) {
                    int r = a % b;
                    a = b;
                    b = r;
                }
                if (a == 1) phi++;
            }
            total += (long)phi * (1L << (L / d));
        }
    }
    return total / L;
}

static void test_structure(void)
{
    int Ls[] = {6, 8, 10, 12, 14};
    int i;
    char name[128];
    printf("\n== 1. structure ==\n");
    for (i = 0; i < 5; i++) {
        EDSolver ed;
        ed_init(&ed, Ls[i], 1.0, 0.01, 1);
        sprintf(name, "k0 dimension L=%d", Ls[i]);
        check(name, (double)ed.basis.dim, (double)necklaces(Ls[i]), 0.0);
        sprintf(name, "A symmetric (k0, eps=0.01) L=%d", Ls[i]);
        check(name, tgen_asymmetry(&ed.tg), 0.0, 1e-14);
        ed_free(&ed);
        ed_init(&ed, Ls[i], 1.0, 0.0, 1);
        sprintf(name, "k0 dimension, empty excluded L=%d", Ls[i]);
        check(name, (double)ed.basis.dim, (double)necklaces(Ls[i]) - 1.0, 0.0);
        sprintf(name, "A symmetric (k0, eps=0) L=%d", Ls[i]);
        check(name, tgen_asymmetry(&ed.tg), 0.0, 1e-14);
        ed_free(&ed);
    }
}

static void test_equilibrium(void)
{
    double epss[] = {0.0, 0.0003, 0.003, 0.01, 0.3, 1.0};
    int Ls[] = {6, 9, 12};
    int ie, il;
    char name[128];
    printf("\n== 2. equilibrium (s = g = 0) ==\n");
    for (il = 0; il < 3; il++) {
        for (ie = 0; ie < 6; ie++) {
            EDSolver ed;
            EDPoint p;
            int L = Ls[il];
            double eps = epss[ie];
            double c, kx, ex;
            ed_init(&ed, L, 1.0, eps, 1);
            c = ed.model.c;
            ed_point(&ed, 0.0, 0.0, CHI_NONE, &p);
            if (eps > 0.0) {
                kx = 2.0 * c * (1.0 - c) * (c + eps * (1.0 - c));
                ex = c;
            } else { /* product measure conditioned on a non-empty configuration */
                double pne = 1.0 - pow(1.0 - c, L);
                kx = 2.0 * c * c * (1.0 - c) / pne;
                ex = c / pne;
            }
            sprintf(name, "theta(0,0)=0       L=%d eps=%g", L, eps);
            check(name, p.theta, 0.0, 1e-12);
            sprintf(name, "k(0,0)             L=%d eps=%g", L, eps);
            check(name, p.k, kx, 1e-10);
            sprintf(name, "e(0,0)             L=%d eps=%g", L, eps);
            check(name, p.e, ex, 1e-10);
            ed_free(&ed);
        }
    }
    /* other temperature */
    {
        EDSolver ed;
        EDPoint p;
        double c, kx;
        ed_init(&ed, 10, 4.0 / 3.0, 0.01, 1);
        c = ed.model.c;
        ed_point(&ed, 0.0, 0.0, CHI_NONE, &p);
        kx = 2.0 * c * (1.0 - c) * (c + 0.01 * (1.0 - c));
        check("c(T=4/3) = 1/(1+exp(0.75))", c, 1.0 / (1.0 + exp(0.75)), 1e-15);
        check("k(0,0) T=4/3 L=10 eps=0.01", p.k, kx, 1e-10);
        check("e(0,0) T=4/3 L=10 eps=0.01", p.e, c, 1e-10);
        ed_free(&ed);
    }
}

static void test_glauber(void)
{
    double pts[][2] = {{0.0, 0.0}, {0.2, 0.0}, {-0.3, 0.1}, {0.1, -0.15}, {0.5, 0.15}};
    int Ls[] = {6, 9};
    int il, ip, ib;
    char name[128];
    printf("\n== 3. eps = 1 (independent spins, closed form) ==\n");
    for (il = 0; il < 2; il++) {
        for (ib = 0; ib < 2; ib++) {
            EDSolver ed;
            int L = Ls[il];
            ed_init(&ed, L, 1.0, 1.0, ib == 0);
            for (ip = 0; ip < 5; ip++) {
                double s = pts[ip][0], g = pts[ip][1];
                double c = ed.model.c;
                EDPoint p;
                const char *bn = ib == 0 ? "k0" : "full";
                ed_point(&ed, s, g, CHI_LR, &p);
                sprintf(name, "theta  L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.theta, L * glauber_lambda(c, s, g), 1e-11);
                sprintf(name, "k      L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.k, -d1r(c, s, g, 1), 1e-9);
                sprintf(name, "e      L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.e, -d1r(c, s, g, 0), 1e-9);
                sprintf(name, "chi_KK L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.chi_kk, d2r(c, s, g, 0), 1e-7);
                sprintf(name, "chi_EE L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.chi_ee, d2r(c, s, g, 1), 1e-7);
                sprintf(name, "chi_KE L=%d %s s=%g g=%g", L, bn, s, g);
                check(name, p.chi_ke, d2r(c, s, g, 2), 1e-7);
            }
            ed_free(&ed);
        }
    }
}

static void test_k0_vs_full(void)
{
    double pts[][2] = {{0.0, 0.0}, {0.03, 0.0}, {0.017, 0.0}, {-0.1, 0.1}, {0.2, -0.1}};
    double epss[] = {0.0, 0.003, 0.01};
    int ie, ip;
    char name[128];
    printf("\n== 4. zero-momentum sector vs full space (L = 8) ==\n");
    for (ie = 0; ie < 3; ie++) {
        EDSolver a, b;
        ed_init(&a, 8, 1.0, epss[ie], 1);
        ed_init(&b, 8, 1.0, epss[ie], 0);
        for (ip = 0; ip < 5; ip++) {
            EDPoint pa, pb;
            double s = pts[ip][0], g = pts[ip][1];
            ed_point(&a, s, g, CHI_LR, &pa);
            ed_point(&b, s, g, CHI_LR, &pb);
            sprintf(name, "theta  eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.theta, pb.theta, 1e-11);
            sprintf(name, "k      eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.k, pb.k, 1e-9);
            sprintf(name, "e      eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.e, pb.e, 1e-9);
            sprintf(name, "chi_KK eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.chi_kk, pb.chi_kk, 1e-7);
            sprintf(name, "chi_KE eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.chi_ke, pb.chi_ke, 1e-7);
            sprintf(name, "chi_EE eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, pa.chi_ee, pb.chi_ee, 1e-7);
        }
        ed_free(&a);
        ed_free(&b);
    }
}

static void test_dense(void)
{
    double pts[][2] = {{0.0, 0.0}, {0.0166, 0.0}, {0.0175, 0.0}, {-0.15, 0.1}, {0.3, -0.1},
                       {0.1, -0.05}};
    double epss[] = {0.0, 0.003, 0.01};
    int Ls[] = {8, 10};
    int il, ie, ip;
    char name[128];
    printf("\n== 5. Lanczos + CG vs dense diagonalization (k0 sector) ==\n");
    for (il = 0; il < 2; il++) {
        for (ie = 0; ie < 3; ie++) {
            EDSolver ed;
            ed_init(&ed, Ls[il], 1.0, epss[ie], 1);
            for (ip = 0; ip < 6; ip++) {
                EDPoint p;
                DenseRef r;
                double s = pts[ip][0], g = pts[ip][1];
                memset(&r, 0, sizeof(r));
                ed_point(&ed, s, g, CHI_LR, &p);
                if (dense_reference(&ed, s, g, &r) != 0) {
                    check_true("dense reference failed", 0);
                    continue;
                }
                sprintf(name, "theta  L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.theta, r.theta, 1e-11);
                sprintf(name, "theta2 L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.theta2, r.theta2, 1e-9);
                sprintf(name, "k      L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.k, r.k, 1e-9);
                sprintf(name, "e      L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.e, r.e, 1e-9);
                sprintf(name, "chi_KK L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.chi_kk, r.chi_kk, 1e-7);
                sprintf(name, "chi_KE L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.chi_ke, r.chi_ke, 1e-7);
                sprintf(name, "chi_EE L=%d eps=%g s=%g g=%g", Ls[il], epss[ie], s, g);
                check(name, p.chi_ee, r.chi_ee, 1e-7);
            }
            ed_free(&ed);
        }
    }
}

static void test_finite_differences(void)
{
    double pts[][2] = {{0.0, 0.0}, {0.02, 0.05}, {-0.08, 0.1}, {0.12, -0.1}};
    double epss[] = {0.003, 0.01, 0.1};
    int ie, ip;
    const int L = 10;
    char name[128];
    printf("\n== 6. derivatives vs finite differences (L = 10) ==\n");
    for (ie = 0; ie < 3; ie++) {
        EDSolver ed;
        ed_init(&ed, L, 1.0, epss[ie], 1);
        for (ip = 0; ip < 4; ip++) {
            /* central differences with steps h and h/2, Richardson-extrapolated:
               near the crossover theta varies on a scale ~1e-3, so a plain
               O(h^2) difference is not accurate enough for a strict test */
            double s = pts[ip][0], g = pts[ip][1];
            double hs[2] = {2e-5, 1e-5};
            double kfd[2], efd[2], dkdg[2], deds[2], ckk[2], cee[2];
            EDPoint p0;
            int ih;
            ed_point(&ed, s, g, CHI_LR, &p0);
            for (ih = 0; ih < 2; ih++) {
                const double h = hs[ih];
                EDPoint a, b, cc, d, pfd;
                ed_point(&ed, s + h, g, CHI_NONE, &a);
                ed_point(&ed, s - h, g, CHI_NONE, &b);
                ed_point(&ed, s, g + h, CHI_NONE, &cc);
                ed_point(&ed, s, g - h, CHI_NONE, &d);
                kfd[ih] = -(a.theta - b.theta) / (2.0 * h * L);
                efd[ih] = -(cc.theta - d.theta) / (2.0 * h * L);
                dkdg[ih] = (cc.k - d.k) / (2.0 * h);
                deds[ih] = (a.e - b.e) / (2.0 * h);
                ed.fd_h = h;
                ed_point(&ed, s, g, CHI_FD, &pfd);
                ckk[ih] = pfd.chi_kk;
                cee[ih] = pfd.chi_ee;
            }
#define RICH(x) ((4.0 * (x)[1] - (x)[0]) / 3.0)
            sprintf(name, "k = -(1/L) dtheta/ds  eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, p0.k, RICH(kfd), 1e-8);
            sprintf(name, "e = -(1/L) dtheta/dg  eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, p0.e, RICH(efd), 1e-8);
            sprintf(name, "Maxwell dk/dg = de/ds  eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, RICH(dkdg), RICH(deds), 1e-6);
            sprintf(name, "dk/dg = -chi_KE        eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, RICH(dkdg), -p0.chi_ke, 1e-6);
            sprintf(name, "chi_KK LR vs FD        eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, p0.chi_kk, RICH(ckk), 1e-6);
            sprintf(name, "chi_EE LR vs FD        eps=%g s=%g g=%g", epss[ie], s, g);
            check(name, p0.chi_ee, RICH(cee), 1e-6);
#undef RICH
        }
        ed_free(&ed);
    }
}

/*
 * Product-state solution. With h(s) = exp(-s) sqrt(c(1-c)) sigma^x - diag(c, 1-c)
 * acting on one site and (a, b) its top eigenvector with eigenvalue mu(s),
 * the product state psi = (a, b)^{(x) L} satisfies
 *   H psi = [mu (1 - eps) - g] N psi + mu eps L psi.
 * On the line mu(s) = g / (1 - eps) it is the Perron vector with
 *   theta = L eps g / (1 - eps),
 *   e = rho = b^2 / (a^2 + b^2),
 *   k = exp(-s) sqrt(c(1-c)) [rho + eps (1 - rho)] 2 a b / (a^2 + b^2).
 * For eps = 0 this line is the strict East coexistence line (theta = 0 =
 * eigenvalue of the empty state), independent of L:
 *   exp(-2 s*) = [(1 + 2 g)^2 - (1 - 2 c)^2] / [4 c (1 - c)].
 */
static double product_line_s(double c, double eps, double g)
{
    double gg = g / (1.0 - eps);
    double num = (1.0 + 2.0 * gg) * (1.0 + 2.0 * gg) - (1.0 - 2.0 * c) * (1.0 - 2.0 * c);
    return -0.5 * log(num / (4.0 * c * (1.0 - c)));
}

static void test_strict(void)
{
    int Ls[] = {8, 10, 12, 14};
    double gs[] = {-0.15, -0.1, 0.0, 0.05, 0.15};
    int il, ig;
    char name[128];
    printf("\n== 7. strict East crossing vs closed form ==\n");
    for (il = 0; il < 4; il++) {
        EDSolver ed;
        EDPoint p;
        ed_init(&ed, Ls[il], 1.0, 0.0, 1);
        for (ig = 0; ig < 5; ig++) {
            double g = gs[ig];
            ed_find_strict_crossing(&ed, g, 0.05, 1e-11, &p);
            sprintf(name, "strict s*(g=%g) L=%d", g, Ls[il]);
            check(name, p.s, product_line_s(ed.model.c, 0.0, g), 1e-9);
            sprintf(name, "slope -e/k = ds*/dg (FD of closed form) g=%g L=%d", g, Ls[il]);
            check(name, -p.e / p.k,
                  (product_line_s(ed.model.c, 0.0, g + 1e-5) -
                   product_line_s(ed.model.c, 0.0, g - 1e-5)) / 2e-5,
                  1e-7);
        }
        ed_free(&ed);
    }
}

static void test_product_line(void)
{
    double epss[] = {0.003, 0.01, 0.1};
    double gs[] = {-0.1, 0.0, 0.1};
    int Ls[] = {8, 12};
    int il, ie, ig;
    char name[128];
    printf("\n== 9. softened East: exact product-state line ==\n");
    for (il = 0; il < 2; il++) {
        for (ie = 0; ie < 3; ie++) {
            EDSolver ed;
            const double eps = epss[ie];
            const int L = Ls[il];
            ed_init(&ed, L, 1.0, eps, 1);
            for (ig = 0; ig < 3; ig++) {
                const double c = ed.model.c, q = sqrt(c * (1.0 - c)), g = gs[ig];
                double s = product_line_s(c, eps, g);
                double mu = g / (1.0 - eps);
                double a = 1.0, b = (mu + c) / (exp(-s) * q);
                double nn = a * a + b * b, rho = b * b / nn;
                double kx = exp(-s) * q * (rho + eps * (1.0 - rho)) * 2.0 * a * b / nn;
                EDPoint p;
                ed_point(&ed, s, g, CHI_NONE, &p);
                sprintf(name, "theta = L eps g/(1-eps)  L=%d eps=%g g=%g", L, eps, g);
                check(name, p.theta, L * eps * g / (1.0 - eps), 1e-11);
                sprintf(name, "e = rho                  L=%d eps=%g g=%g", L, eps, g);
                check(name, p.e, rho, 1e-10);
                sprintf(name, "k (product state)        L=%d eps=%g g=%g", L, eps, g);
                check(name, p.k, kx, 1e-10);
            }
            ed_free(&ed);
        }
    }
}

static void test_perron(void)
{
    double pts[][2] = {{0.0, 0.0}, {0.5, 0.0}, {-0.5, 0.0}, {0.02, 0.15}, {0.0, -0.15}};
    int ip;
    EDSolver ed;
    char name[128];
    printf("\n== 8. Perron vector positivity (L = 12, eps = 0.003) ==\n");
    ed_init(&ed, 12, 1.0, 0.003, 1);
    for (ip = 0; ip < 5; ip++) {
        EDPoint p;
        ed_point(&ed, pts[ip][0], pts[ip][1], CHI_NONE, &p);
        sprintf(name, "min/max component >= 0 at s=%g g=%g (%.3e)", pts[ip][0], pts[ip][1],
                p.min_ratio);
        check_true(name, p.min_ratio >= -1e-10);
        sprintf(name, "gap > 0 at s=%g g=%g", pts[ip][0], pts[ip][1]);
        check_true(name, p.theta - p.theta2 > 0.0);
    }
    ed_free(&ed);
}

/* ---------- discrete-time transfer matrix ---------- */

typedef struct {
    int L, nattempt;
    double c, eps, s, g;
    double Z, ZK, ZE;
} Enum;

/* brute-force sum over every attempt sequence (site choice and flip outcome) */
static void enum_rec(Enum *en, uint32_t x, int a, double w, int K, int E)
{
    const int L = en->L;
    int i;
    if (a > 0 && a % L == 0) { /* end of a sweep: energy weight of the new slice */
        int n = 0, q;
        for (q = 0; q < L; q++) n += (x >> q) & 1u;
        w *= exp(-en->g * n);
        E += n;
    }
    if (a == en->nattempt) {
        en->Z += w;
        en->ZK += w * K;
        en->ZE += w * E;
        return;
    }
    for (i = 0; i < L; i++) {
        int left = (x >> (i == 0 ? L - 1 : i - 1)) & 1u;
        int ni = (x >> i) & 1u;
        double p = (left ? 1.0 : en->eps) * (ni ? 1.0 - en->c : en->c);
        if (p > 0.0) enum_rec(en, x ^ (1u << i), a + 1, w * p / L * exp(-en->s), K + 1, E);
        if (p < 1.0) enum_rec(en, x, a + 1, w * (1.0 - p) / L, K, E);
    }
}

static void test_transfer_matrix(void)
{
    char name[160];
    int ie;
    double epss[] = {0.0, 0.01, 0.1};
    printf("\n== 10. discrete-time transfer matrix ==\n");

    /* brute-force enumeration of all attempt sequences: L = 3, tobs = 2 */
    {
        double pts[][2] = {{0.0, 0.0}, {0.3, 0.0}, {-0.2, 0.15}, {0.1, -0.3}};
        int ip;
        for (ie = 0; ie < 3; ie++) {
            TMSolver tm;
            tm_init(&tm, 3, 1.0, epss[ie], 0, 0);
            for (ip = 0; ip < 4; ip++) {
                Enum en;
                TMFinite r;
                uint32_t x0;
                memset(&en, 0, sizeof(en));
                en.L = 3;
                en.nattempt = 3 * 2;
                en.c = tm.model.c;
                en.eps = epss[ie];
                en.s = pts[ip][0];
                en.g = pts[ip][1];
                for (x0 = 0; x0 < 8u; x0++) {
                    int n = conf_popcount(x0);
                    double pi = pow(en.c, n) * pow(1.0 - en.c, 3 - n);
                    enum_rec(&en, x0, 0, pi * exp(-en.g * n), 0, n);
                }
                tm_finite(&tm, en.s, en.g, 2, &r);
                sprintf(name, "lnZ vs enumeration  L=3 tobs=2 eps=%g s=%g g=%g", en.eps, en.s, en.g);
                check(name, r.lnZ, log(en.Z), 1e-12);
                sprintf(name, "<K> vs enumeration  L=3 tobs=2 eps=%g s=%g g=%g", en.eps, en.s, en.g);
                check(name, r.meanK, en.ZK / en.Z, 1e-12);
                sprintf(name, "<E> vs enumeration  L=3 tobs=2 eps=%g s=%g g=%g", en.eps, en.s, en.g);
                check(name, r.meanE, en.ZE / en.Z, 1e-12);
            }
            tm_free(&tm);
        }
    }

    /* finite time: equilibrium values and derivatives */
    for (ie = 0; ie < 3; ie++) {
        TMSolver tm, tf;
        TMFinite r0, r, rf;
        const int L = 8, t = 30;
        double c, kx;
        tm_init(&tm, L, 1.0, epss[ie], 1, 0);
        tm_init(&tf, L, 1.0, epss[ie], 0, 0);
        c = tm.model.c;
        kx = 2.0 * c * (1.0 - c) * (c + epss[ie] * (1.0 - c));
        tm_finite(&tm, 0.0, 0.0, t, &r0);
        sprintf(name, "finite lnZ(0,0) = 0              eps=%g", epss[ie]);
        check(name, r0.lnZ, 0.0, 1e-12);
        sprintf(name, "finite <K>(0,0) = t L k_eq       eps=%g", epss[ie]);
        check(name, r0.meanK, t * L * kx, 1e-11);
        sprintf(name, "finite <E>(0,0) = (t+1) L c      eps=%g", epss[ie]);
        check(name, r0.meanE, (t + 1) * L * c, 1e-11);
        {
            const double s = 0.03, g = 0.05;
            double hs[2] = {2e-4, 1e-4}, fK[2], fE[2], fV[2], fC[2], fVE[2];
            int ih;
            tm_finite(&tm, s, g, t, &r);
            tm_finite(&tf, s, g, t, &rf);
            for (ih = 0; ih < 2; ih++) {
                TMFinite a, b, cc, d;
                double h = hs[ih];
                tm_finite(&tm, s + h, g, t, &a);
                tm_finite(&tm, s - h, g, t, &b);
                tm_finite(&tm, s, g + h, t, &cc);
                tm_finite(&tm, s, g - h, t, &d);
                fK[ih] = -(a.lnZ - b.lnZ) / (2 * h);
                fE[ih] = -(cc.lnZ - d.lnZ) / (2 * h);
                fV[ih] = -(a.meanK - b.meanK) / (2 * h);
                fVE[ih] = -(cc.meanE - d.meanE) / (2 * h);
                fC[ih] = -(cc.meanK - d.meanK) / (2 * h);
            }
#define RICH(x) ((4.0 * (x)[1] - (x)[0]) / 3.0)
            sprintf(name, "finite <K> = -dlnZ/ds            eps=%g", epss[ie]);
            check(name, r.meanK, RICH(fK), 1e-8);
            sprintf(name, "finite <E> = -dlnZ/dg            eps=%g", epss[ie]);
            check(name, r.meanE, RICH(fE), 1e-8);
            sprintf(name, "finite VarK = -d<K>/ds           eps=%g", epss[ie]);
            check(name, r.varK, RICH(fV), 1e-7);
            sprintf(name, "finite VarE = -d<E>/dg           eps=%g", epss[ie]);
            check(name, r.varE, RICH(fVE), 1e-7);
            sprintf(name, "finite CovKE = -d<K>/dg          eps=%g", epss[ie]);
            check(name, r.covKE, RICH(fC), 1e-7);
#undef RICH
            sprintf(name, "finite k0 vs full basis lnZ      eps=%g", epss[ie]);
            check(name, r.lnZ, rf.lnZ, 1e-11);
            sprintf(name, "finite k0 vs full basis VarK     eps=%g", epss[ie]);
            check(name, r.varK, rf.varK, 1e-10);
        }
        tm_free(&tm);
        tm_free(&tf);
    }

    /* asymptotic */
    for (ie = 0; ie < 3; ie++) {
        const double eps = epss[ie];
        TMSolver tm, tf;
        TMAsym a0, a, af;
        const int L = 8;
        double c, kx, ex, pne;
        tm_init(&tm, L, 1.0, eps, 1, eps == 0.0);
        tm_init(&tf, L, 1.0, eps, 0, eps == 0.0);
        c = tm.model.c;
        pne = (eps == 0.0) ? 1.0 - pow(1.0 - c, L) : 1.0;
        kx = 2.0 * c * (1.0 - c) * (c + eps * (1.0 - c)) / pne;
        ex = c / pne;
        tm_asym(&tm, 0.0, 0.0, &a0);
        sprintf(name, "asym theta_d(0,0) = 0            eps=%g", eps);
        check(name, a0.theta, 0.0, 1e-12);
        sprintf(name, "asym k_d(0,0) equilibrium        eps=%g", eps);
        check(name, a0.k, kx, 1e-10);
        sprintf(name, "asym e_d(0,0) equilibrium        eps=%g", eps);
        check(name, a0.e, ex, 1e-10);
        {
            const double s = 0.03, g = 0.05;
            double hs[2] = {2e-4, 1e-4}, fk[2], fe[2];
            int ih;
            tm_asym(&tm, s, g, &a);
            tm_asym(&tf, s, g, &af);
            for (ih = 0; ih < 2; ih++) {
                TMAsym p, q, u, w;
                double h = hs[ih];
                tm_asym(&tm, s + h, g, &p);
                tm_asym(&tm, s - h, g, &q);
                tm_asym(&tm, s, g + h, &u);
                tm_asym(&tm, s, g - h, &w);
                fk[ih] = -(p.theta - q.theta) / (2 * h * L);
                fe[ih] = -(u.theta - w.theta) / (2 * h * L);
            }
            sprintf(name, "asym k_d = -(1/L) dtheta/ds      eps=%g", eps);
            check(name, a.k, (4.0 * fk[1] - fk[0]) / 3.0, 1e-8);
            sprintf(name, "asym e_d = -(1/L) dtheta/dg      eps=%g", eps);
            check(name, a.e, (4.0 * fe[1] - fe[0]) / 3.0, 1e-8);
            sprintf(name, "asym k0 vs full basis theta      eps=%g", eps);
            check(name, a.theta, af.theta, 1e-11);
        }
        if (eps > 0.0) {
            /* long-time limit of the finite-time free energy */
            TMSolver t2;
            TMFinite f1, f2;
            tm_init(&t2, L, 1.0, eps, 1, 0);
            tm_finite(&t2, 0.05, 0.02, 300, &f1);
            tm_finite(&t2, 0.05, 0.02, 600, &f2);
            tm_asym(&tm, 0.05, 0.02, &a);
            sprintf(name, "[lnZ(600)-lnZ(300)]/300 = theta_d eps=%g", eps);
            check(name, (f2.lnZ - f1.lnZ) / 300.0, a.theta, 1e-9);
            tm_free(&t2);
        }
        tm_free(&tm);
        tm_free(&tf);
    }

    /* eps = 1, g = 0: M = I + H/L with H a sum of independent sites, so
       theta_d = L ln(1 + mu(s)), mu = top eigenvalue of the one-site h(s) */
    {
        double sv[] = {0.1, -0.2, 0.5};
        int i;
        TMSolver tm;
        tm_init(&tm, 8, 1.0, 1.0, 1, 0);
        for (i = 0; i < 3; i++) {
            TMAsym a;
            double c = tm.model.c, s = sv[i];
            double mu = 0.5 * (-1.0 + sqrt((1 - 2 * c) * (1 - 2 * c) +
                                           4.0 * exp(-2.0 * s) * c * (1 - c)));
            tm_asym(&tm, s, 0.0, &a);
            sprintf(name, "asym eps=1 g=0: theta_d = L ln(1+mu(s)), s=%g", s);
            check(name, a.theta, 8.0 * log(1.0 + mu), 1e-11);
        }
        tm_free(&tm);
    }

    /* one-sweep kernel rows are probability distributions */
    {
        TMSolver tm;
        double prob[1 << 6];
        uint32_t xs[] = {0u, 5u, 42u, 63u};
        int i, q;
        tm_init(&tm, 6, 1.0, 0.05, 0, 0);
        for (i = 0; i < 4; i++) {
            double sum = 0.0, mn = 1.0;
            tm_mcs_row(&tm, xs[i], prob);
            for (q = 0; q < 64; q++) {
                sum += prob[q];
                if (prob[q] < mn) mn = prob[q];
            }
            sprintf(name, "MCS row sums to 1 from x=%u", (unsigned)xs[i]);
            check(name, sum, 1.0, 1e-13);
            sprintf(name, "MCS row nonnegative from x=%u", (unsigned)xs[i]);
            check_true(name, mn > -1e-15);
        }
        tm_free(&tm);
    }
}

static void reus_small_params(REUSParams *p)
{
    reus_default_params(p);
    p->L = 8;
    p->tobs = 20;
    p->eps = 0.05;
    p->s = 0.0;
    p->g = 0.0;
    p->seed = 1;
    p->nwin = 6;
    p->k_lo = 0.0;
    p->k_hi = 0.3;
    p->n_approach = 5;
    p->n_equil = 20;
    p->n_prod = 200;
    p->swap_every = 5;
    p->sample_every = 5;
    p->init = "empty";
}

static void test_reus_seeding(void)
{
    REUSParams p;
    REUSResult r;
    int w, rc, nseeded = 0, in_all = 1, final_in = 1, valid_src = 1;
    long i;
    printf("\n== 11. REUS start-up: seeding unreached windows from reached ones ==\n");
    reus_small_params(&p); /* n_approach = 5: too few for the upper windows from empty */
    rc = reus_run(&p, &r);
    check_true("REUS runs although some windows need seeding", rc == 0);
    if (rc != 0) return;
    for (w = 0; w < r.nwin; w++) {
        if (r.seeded_from[w] >= 0) nseeded++;
        if (r.seeded_from[w] == w || r.seeded_from[w] < -1 || r.seeded_from[w] >= r.nwin)
            valid_src = 0;
        if (r.final[w].K < r.Kmin[w] || r.final[w].K > r.Kmax[w]) final_in = 0;
        for (i = 0; i < r.nsamp; i++) {
            double K = r.K[(size_t)w * r.nsamp + i];
            if (K < r.Kmin[w] || K > r.Kmax[w]) in_all = 0;
        }
    }
    check_true("at least one window was seeded from another window", nseeded > 0);
    check_true("seeding source is another valid window", valid_src);
    check_true("all samples lie inside their windows", in_all);
    check_true("final trajectories lie inside their windows", final_in);
    reus_result_free(&r);
}

static int same_long(const long *a, const long *b, int n)
{
    return memcmp(a, b, sizeof(long) * (size_t)n) == 0;
}

static void test_reus_checkpoint(void)
{
    const char *ck = "east_test_reus.ckpt";
    REUSParams p;
    REUSResult a, b;
    int rc1, rc2, rc3, rc4, w, same;
    printf("\n== 12. REUS checkpoint and resume ==\n");
    remove(ck);
    reus_small_params(&p);
    if (reus_run(&p, &a) != 0) {
        check_true("uninterrupted reference run", 0);
        return;
    }
    p.checkpoint = ck;
    p.checkpoint_every = 0.0;
    p.stop_approach_round = 1;
    rc1 = reus_run(&p, &b);
    p.stop_approach_round = 0;
    rc2 = reus_run(&p, &b);
    check_true("approach phase stops with a checkpoint", rc1 == REUS_STOPPED);
    same = rc2 == 0 && a.nsamp == b.nsamp &&
           memcmp(a.K, b.K, sizeof(double) * (size_t)a.nwin * (size_t)a.nsamp) == 0 &&
           memcmp(a.E, b.E, sizeof(double) * (size_t)a.nwin * (size_t)a.nsamp) == 0 &&
           memcmp(a.rep, b.rep, sizeof(int) * (size_t)a.nwin * (size_t)a.nsamp) == 0 &&
           same_long(a.approach, b.approach, a.nwin) &&
           memcmp(a.seeded_from, b.seeded_from, sizeof(int) * (size_t)a.nwin) == 0;
    check_true("approach checkpoint resumes bit-identically", same);
    if (rc2 == 0) reus_result_free(&b);
    remove(ck);

    p.checkpoint = ck;
    p.checkpoint_every = 0.0;
    p.stop_cycle = 3; /* during equilibration (4 cycles) */
    rc1 = reus_run(&p, &b);
    p.stop_cycle = 25; /* during production */
    rc2 = reus_run(&p, &b);
    p.stop_cycle = 0;
    rc3 = reus_run(&p, &b);
    check_true("first part stops with a checkpoint (REUS_STOPPED)", rc1 == REUS_STOPPED);
    check_true("second part resumes and stops again", rc2 == REUS_STOPPED);
    check_true("third part resumes and finishes", rc3 == 0);
    if (rc3 != 0) {
        reus_result_free(&a);
        remove(ck);
        return;
    }
    same = a.nsamp == b.nsamp &&
           memcmp(a.K, b.K, sizeof(double) * (size_t)a.nwin * (size_t)a.nsamp) == 0 &&
           memcmp(a.E, b.E, sizeof(double) * (size_t)a.nwin * (size_t)a.nsamp) == 0 &&
           memcmp(a.rep, b.rep, sizeof(int) * (size_t)a.nwin * (size_t)a.nsamp) == 0;
    check_true("resumed run: samples (K, E, replica) identical", same);
    same = same_long(a.moves, b.moves, a.nwin) && same_long(a.acc, b.acc, a.nwin) &&
           same_long(a.swap_try, b.swap_try, a.nwin) && same_long(a.swap_acc, b.swap_acc, a.nwin) &&
           same_long(a.round_trips, b.round_trips, a.nwin) &&
           same_long(a.approach, b.approach, a.nwin) &&
           memcmp(a.seeded_from, b.seeded_from, sizeof(int) * (size_t)a.nwin) == 0 &&
           a.seed_used == b.seed_used;
    check_true("resumed run: counters, approach and seeding identical", same);
    same = 1;
    for (w = 0; w < a.nwin; w++)
        if (a.final[w].K != b.final[w].K || a.final[w].E != b.final[w].E ||
            memcmp(a.final[w].x, b.final[w].x, (size_t)(p.tobs + 1) * (size_t)p.L) != 0)
            same = 0;
    check_true("resumed run: final trajectories identical", same);
    reus_result_free(&b);

    /* a checkpoint with other parameters must be refused */
    p.stop_cycle = 3;
    remove(ck);
    rc1 = reus_run(&p, &b);
    p.stop_cycle = 0;
    p.s = 0.01;
    rc4 = reus_run(&p, &b);
    check_true("checkpoint with different s is refused", rc1 == REUS_STOPPED && rc4 == -1);
    if (rc4 == 0) reus_result_free(&b);
    reus_result_free(&a);
    remove(ck);
}

int main(void)
{
    test_structure();
    test_equilibrium();
    test_glauber();
    test_k0_vs_full();
    test_dense();
    test_finite_differences();
    test_strict();
    test_perron();
    test_product_line();
    test_transfer_matrix();
    test_reus_seeding();
    test_reus_checkpoint();
    printf("\n== summary: %d passed, %d failed ==\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
