/*
 * east_edan: analysis of east_ed peak tables for the hypotheses.
 *
 *   h2     shift of the coexistence line: r = (s*_soft - s*_strict) / eps for
 *          every L, extrapolated to L -> infinity from the three largest L
 *   h1     finite-size growth at s*: alpha = dln chi_KK / dln L and
 *          beta = dln (L gap) / dln L for successive L
 *   cross  eps at which alpha crosses given thresholds and beta crosses 0,
 *          for every g and pair of successive L
 *   h3     two-level estimates of the jumps at s*:
 *          dk2 = sqrt(2 chi_KK gap / L), de2 = sqrt(2 chi_EE gap / L)
 *
 * The strict line is s*_strict(g) = -(1/2) ln{[(1+2g)^2 - (1-2c)^2] / [4c(1-c)]},
 * with c read from the table header.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXROW 20000
#define MAXV 256

typedef struct {
    int L;
    double eps, g, s, chi, theta, gap, k, e, chi_ke, chi_ee, rho;
    int edge, ok;
} Row;

typedef struct {
    int L;
    double g, k, e;
} StrictRow;

static Row rows[MAXROW];
static int nrow = 0;
static StrictRow srows[MAXROW];
static int nsrow = 0;
static double c_val = NAN;

static int Ls[MAXV], nL = 0;
static double epss[MAXV], gs[MAXV];
static int neps = 0, ng = 0;

static int same(double a, double b) { return fabs(a - b) <= 1e-12 * (1.0 + fabs(a) + fabs(b)); }

static void pnum(FILE *f, const char *fmt, int width, double x)
{
    if (isnan(x))
        fprintf(f, " %*s", width, "nan");
    else
        fprintf(f, fmt, x);
}

static int read_c(const char *line, double *c)
{
    const char *p = strstr(line, "c = ");
    if (strncmp(line, "# T =", 5) != 0 || !p) return 0;
    *c = strtod(p + 4, NULL);
    return 1;
}

static int read_peak(const char *path)
{
    FILE *f = fopen(path, "r");
    char line[4096];
    int nread = 0;
    if (!f) {
        fprintf(stderr, "error: cannot open %s\n", path);
        return -1;
    }
    while (fgets(line, sizeof line, f)) {
        Row r;
        char *p = line, *q;
        double v[14];
        int i, j, dup = 0;
        if (line[0] == '#') {
            double c;
            if (read_c(line, &c)) {
                if (!isnan(c_val) && fabs(c - c_val) > 1e-9) {
                    fprintf(stderr, "error: %s has c = %.10g, other tables c = %.10g\n", path, c,
                            c_val);
                    fclose(f);
                    return -1;
                }
                c_val = c;
            }
            continue;
        }
        for (i = 0; i < 14; i++) {
            v[i] = strtod(p, &q);
            if (q == p) break;
            p = q;
        }
        if (i == 0) continue; /* blank line */
        if (i < 14) {
            fprintf(stderr, "error: %s: expected 14 columns (peak table): %s", path, line);
            fclose(f);
            return -1;
        }
        r.L = (int)v[0];
        r.eps = v[1];
        r.g = v[2];
        r.s = v[3];
        r.chi = v[4];
        r.theta = v[5];
        r.gap = v[6];
        r.k = v[7];
        r.e = v[8];
        r.chi_ke = v[9];
        r.chi_ee = v[10];
        r.rho = v[11];
        r.edge = (int)v[12];
        r.ok = (int)v[13];
        for (j = 0; j < nrow; j++)
            if (rows[j].L == r.L && same(rows[j].eps, r.eps) && same(rows[j].g, r.g)) {
                if (!(fabs(rows[j].s - r.s) < 1e-6) && !(isnan(rows[j].s) && isnan(r.s)))
                    fprintf(stderr, "warning: L=%d eps=%g g=%g appears twice with s* %.8g, %.8g; "
                                    "the first is used\n",
                            r.L, r.eps, r.g, rows[j].s, r.s);
                dup = 1;
                break;
            }
        if (dup) continue;
        if (nrow >= MAXROW) {
            fprintf(stderr, "error: too many rows\n");
            fclose(f);
            return -1;
        }
        rows[nrow++] = r;
        nread++;
    }
    fclose(f);
    if (nread == 0) fprintf(stderr, "warning: no rows read from %s\n", path);
    return 0;
}

static int read_strict(const char *path)
{
    FILE *f = fopen(path, "r");
    char line[4096];
    if (!f) {
        fprintf(stderr, "error: cannot open %s\n", path);
        return -1;
    }
    while (fgets(line, sizeof line, f)) {
        char *p = line, *q;
        double v[6];
        int i, j, dup = 0;
        if (line[0] == '#') continue;
        for (i = 0; i < 6; i++) {
            v[i] = strtod(p, &q);
            if (q == p) break;
            p = q;
        }
        if (i == 0) continue;
        if (i < 6) {
            fprintf(stderr, "error: %s: not a strict table: %s", path, line);
            fclose(f);
            return -1;
        }
        for (j = 0; j < nsrow; j++)
            if (srows[j].L == (int)v[0] && same(srows[j].g, v[1])) dup = 1;
        if (dup || nsrow >= MAXROW) continue;
        srows[nsrow].L = (int)v[0];
        srows[nsrow].g = v[1];
        srows[nsrow].k = v[4];
        srows[nsrow].e = v[5];
        nsrow++;
    }
    fclose(f);
    return 0;
}

static int cmp_int(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }
static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void add_unique_d(double *arr, int *n, double x)
{
    int i;
    for (i = 0; i < *n; i++)
        if (same(arr[i], x)) return;
    if (*n < MAXV) arr[(*n)++] = x;
}

static void build_lists(void)
{
    int i, j;
    for (i = 0; i < nrow; i++) {
        int found = 0;
        for (j = 0; j < nL; j++)
            if (Ls[j] == rows[i].L) found = 1;
        if (!found && nL < MAXV) Ls[nL++] = rows[i].L;
        add_unique_d(epss, &neps, rows[i].eps);
        add_unique_d(gs, &ng, rows[i].g);
    }
    qsort(Ls, (size_t)nL, sizeof(int), cmp_int);
    qsort(epss, (size_t)neps, sizeof(double), cmp_dbl);
    qsort(gs, (size_t)ng, sizeof(double), cmp_dbl);
}

static const Row *find(int L, double eps, double g)
{
    int i;
    for (i = 0; i < nrow; i++)
        if (rows[i].L == L && same(rows[i].eps, eps) && same(rows[i].g, g)) return &rows[i];
    return NULL;
}

static const StrictRow *find_strict(int L, double g)
{
    int i;
    for (i = 0; i < nsrow; i++)
        if (srows[i].L == L && same(srows[i].g, g)) return &srows[i];
    return NULL;
}

static double s_strict(double g)
{
    double c = c_val;
    double a = ((1.0 + 2.0 * g) * (1.0 + 2.0 * g) - (1.0 - 2.0 * c) * (1.0 - 2.0 * c)) /
               (4.0 * c * (1.0 - c));
    return a > 0.0 ? -0.5 * log(a) + 0.0 : NAN; /* + 0.0: no negative zero */
}

/* x(L) -> infinity from three equally spaced sizes, geometric increments (Aitken). */
static double extrap_geo(const int *L, const double *x)
{
    double d1 = x[1] - x[0], d2 = x[2] - x[1], q;
    if (L[1] - L[0] != L[2] - L[1]) return NAN;
    if (!(d1 * d2 > 0.0) || fabs(d2) >= fabs(d1)) return NAN;
    q = d2 / d1;
    return x[2] + d2 * q / (1.0 - q);
}

/* x(L) = x_inf - A L^-p through three points; NAN if no p > 0 fits. */
static double extrap_pow(const int *L, const double *x, double *p_out)
{
    double d1 = x[1] - x[0], d2 = x[2] - x[1], rho, lo = 1e-3, hi = 50.0, p, A;
    int it;
    *p_out = NAN;
    if (!(d1 * d2 > 0.0)) return NAN;
    rho = d2 / d1;
#define RATIO(pp)                                                                                  \
    ((pow(L[1], -(pp)) - pow(L[2], -(pp))) / (pow(L[0], -(pp)) - pow(L[1], -(pp))))
    /* RATIO decreases from its p -> 0 limit to 0 as p grows */
    if (!(rho < RATIO(lo)) || !(rho > RATIO(hi))) return NAN;
    for (it = 0; it < 200; it++) {
        double mid = 0.5 * (lo + hi);
        if (RATIO(mid) > rho)
            lo = mid;
        else
            hi = mid;
    }
    p = 0.5 * (lo + hi);
#undef RATIO
    A = d2 / (pow(L[1], -p) - pow(L[2], -p));
    *p_out = p;
    return x[2] + A * pow(L[2], -p);
}

static int usable(const Row *r) { return r && r->ok && !r->edge && !isnan(r->s); }

static void run_h2(FILE *f)
{
    int ig, ie, i;
    double pred = 1.0 / (2.0 * c_val * (1.0 - c_val));
    fprintf(f, "# h2: r = (s*_soft - s*_strict) / eps; hypothesis 2 predicts r = 1/[2c(1-c)] = "
               "%.6f\n",
            pred);
    fprintf(f, "# r_geo: geometric (Aitken) extrapolation from the three largest L (equal steps)\n"
               "# r_pow: power-law fit r(L) = r_inf - A L^-p through the same three points\n"
               "# (nan where the three values are not monotone with shrinking steps)\n");
    fprintf(f, "#      eps          g   s_strict");
    for (i = 0; i < nL; i++) fprintf(f, "     r(L=%2d)", Ls[i]);
    fprintf(f, " Lmax        r_geo        r_pow        p     r_pred\n");
    for (ig = 0; ig < ng; ig++)
        for (ie = 0; ie < neps; ie++) {
            int Lx[MAXV], n = 0, any = 0;
            double rx[MAXV], geo = NAN, pw = NAN, p = NAN, ss = s_strict(gs[ig]);
            for (i = 0; i < nL; i++)
                if (find(Ls[i], epss[ie], gs[ig])) any = 1;
            if (!any) continue;
            fprintf(f, "%10.6g %10.6g %10.6f", epss[ie], gs[ig], ss);
            for (i = 0; i < nL; i++) {
                const Row *r = find(Ls[i], epss[ie], gs[ig]);
                double v = usable(r) ? (r->s - ss) / epss[ie] : NAN;
                pnum(f, " %11.6f", 11, v);
                if (!isnan(v)) {
                    Lx[n] = Ls[i];
                    rx[n++] = v;
                }
            }
            if (n >= 3) {
                geo = extrap_geo(Lx + n - 3, rx + n - 3);
                pw = extrap_pow(Lx + n - 3, rx + n - 3, &p);
            }
            fprintf(f, " %4d", n > 0 ? Lx[n - 1] : 0);
            pnum(f, " %12.6f", 12, geo);
            pnum(f, " %12.6f", 12, pw);
            pnum(f, " %8.3f", 8, p);
            fprintf(f, " %10.6f\n", pred);
        }
}

static double alpha_of(const Row *a, const Row *b)
{
    if (!usable(a) || !usable(b)) return NAN;
    return log(b->chi / a->chi) / log((double)b->L / a->L);
}

static double beta_of(const Row *a, const Row *b)
{
    if (!usable(a) || !usable(b)) return NAN;
    return log((b->L * b->gap) / (a->L * a->gap)) / log((double)b->L / a->L);
}

static void run_h1(FILE *f)
{
    int ig, ie, i;
    fprintf(f, "# h1: alpha = dln chi_KK_max / dln L, beta = dln (L gap) / dln L between "
               "successive L\n"
               "# first-order-like: chi grows (alpha >= 1), L gap falls (beta < 0);\n"
               "# crossover: alpha -> 0, beta -> 1; beta = 0 if gap ~ 1/L\n");
    fprintf(f, "#      eps          g");
    for (i = 0; i + 1 < nL; i++) fprintf(f, "  alpha%02d-%02d", Ls[i], Ls[i + 1]);
    for (i = 0; i + 1 < nL; i++) fprintf(f, "   beta%02d-%02d", Ls[i], Ls[i + 1]);
    fprintf(f, "\n");
    for (ig = 0; ig < ng; ig++)
        for (ie = 0; ie < neps; ie++) {
            int any = 0;
            for (i = 0; i < nL; i++)
                if (find(Ls[i], epss[ie], gs[ig])) any = 1;
            if (!any) continue;
            fprintf(f, "%10.6g %10.6g", epss[ie], gs[ig]);
            for (i = 0; i + 1 < nL; i++)
                pnum(f, " %12.5f", 12,
                     alpha_of(find(Ls[i], epss[ie], gs[ig]), find(Ls[i + 1], epss[ie], gs[ig])));
            for (i = 0; i + 1 < nL; i++)
                pnum(f, " %12.5f", 12,
                     beta_of(find(Ls[i], epss[ie], gs[ig]), find(Ls[i + 1], epss[ie], gs[ig])));
            fprintf(f, "\n");
        }
}

/* Crossing of y(eps) - y0 on the sorted eps grid, interpolated in ln eps.
   where: 0 found (first crossing), 1 none and y > y0 everywhere, -1 none and
   y < y0 everywhere, 2 more than one crossing (the first is returned),
   3 fewer than two points. */
static double crossing(const double *e, const double *y, int n, double y0, int *where, int *npts)
{
    int i, m = 0, nx = 0, above = 0, below = 0;
    double ee[MAXV], yy[MAXV], first = NAN;
    for (i = 0; i < n; i++)
        if (!isnan(y[i])) {
            ee[m] = e[i];
            yy[m++] = y[i];
        }
    *npts = m;
    if (m < 2) {
        *where = 3;
        return NAN;
    }
    for (i = 0; i < m; i++) {
        if (yy[i] > y0) above++;
        if (yy[i] < y0) below++;
    }
    for (i = 0; i + 1 < m; i++) {
        double a = yy[i] - y0, b = yy[i + 1] - y0;
        if (a == 0.0 || a * b < 0.0) {
            double t = (a == 0.0) ? 0.0 : a / (a - b);
            double x = exp(log(ee[i]) + t * (log(ee[i + 1]) - log(ee[i])));
            if (nx == 0) first = x;
            nx++;
        }
    }
    if (nx == 0) {
        *where = below == 0 ? 1 : -1;
        return NAN;
    }
    *where = nx == 1 ? 0 : 2;
    return first;
}

static void run_cross(FILE *f, const double *thr, int nthr)
{
    int ig, i, it;
    fprintf(f, "# cross: eps_x where alpha(eps) = threshold or beta(eps) = 0, interpolated in "
               "ln eps,\n"
               "# for every g and pair of successive L (eps grid of the tables)\n"
               "# where: 0 single crossing; 2 several crossings (smallest eps given);\n"
               "#        1 none, indicator above the level at every eps; -1 none, below at "
               "every eps;\n"
               "#        3 fewer than two eps values\n");
    fprintf(f, "#          g  L1  L2 indicator     level        eps_x where npts\n");
    for (ig = 0; ig < ng; ig++)
        for (i = 0; i + 1 < nL; i++) {
            double al[MAXV], be[MAXV], x;
            int ie, where, npts;
            for (ie = 0; ie < neps; ie++) {
                const Row *a = find(Ls[i], epss[ie], gs[ig]), *b = find(Ls[i + 1], epss[ie], gs[ig]);
                al[ie] = alpha_of(a, b);
                be[ie] = beta_of(a, b);
            }
            for (it = 0; it < nthr; it++) {
                x = crossing(epss, al, neps, thr[it], &where, &npts);
                if (where == 3) continue;
                fprintf(f, "%12.6g %3d %3d %9s %9.4g", gs[ig], Ls[i], Ls[i + 1], "alpha", thr[it]);
                pnum(f, " %12.5e", 12, x);
                fprintf(f, " %5d %4d\n", where, npts);
            }
            x = crossing(epss, be, neps, 0.0, &where, &npts);
            if (where == 3) continue;
            fprintf(f, "%12.6g %3d %3d %9s %9.4g", gs[ig], Ls[i], Ls[i + 1], "beta", 0.0);
            pnum(f, " %12.5e", 12, x);
            fprintf(f, " %5d %4d\n", where, npts);
        }
}

static void run_h3(FILE *f)
{
    int ig, ie, i;
    fprintf(f, "# h3: two-level estimates at s*: dk2 = sqrt(2 chi_KK gap / L), "
               "de2 = sqrt(2 chi_EE gap / L)\n"
               "# (exact for two mixing levels; rho = chi_KE / sqrt(chi_KK chi_EE) = 1 there)\n"
               "# k_strict, e_strict: active branch on the strict line (strict table, if "
               "given); dk_ratio = dk2 / k_strict\n");
    fprintf(f, "#   L        eps          g          dk2          de2          rho     "
               "k_strict     e_strict     dk_ratio     de_ratio\n");
    for (ig = 0; ig < ng; ig++)
        for (ie = 0; ie < neps; ie++)
            for (i = 0; i < nL; i++) {
                const Row *r = find(Ls[i], epss[ie], gs[ig]);
                const StrictRow *sr = find_strict(Ls[i], gs[ig]);
                double dk = NAN, de = NAN, rho = NAN, ks = NAN, es = NAN;
                if (!r) continue;
                if (usable(r)) {
                    dk = sqrt(2.0 * r->chi * r->gap / r->L);
                    de = sqrt(2.0 * r->chi_ee * r->gap / r->L);
                    rho = r->rho;
                }
                if (sr) {
                    ks = sr->k;
                    es = sr->e;
                }
                fprintf(f, "%5d %10.6g %10.6g", r->L, epss[ie], gs[ig]);
                pnum(f, " %12.6f", 12, dk);
                pnum(f, " %12.6f", 12, de);
                pnum(f, " %12.8f", 12, rho);
                pnum(f, " %12.6f", 12, ks);
                pnum(f, " %12.6f", 12, es);
                pnum(f, " %12.6f", 12, dk / ks);
                pnum(f, " %12.6f", 12, de / es);
                fprintf(f, "\n");
            }
}

static void usage(void)
{
    fprintf(stderr,
            "usage: ./bin/east_edan <mode> [options] <peak tables ...>\n\n"
            "modes:\n"
            "  h2     r = (s*_soft - s*_strict)/eps for every L, extrapolated in L\n"
            "  h1     alpha = dln chi_KK/dln L and beta = dln(L gap)/dln L, successive L\n"
            "  cross  eps where alpha crosses the thresholds and beta crosses 0\n"
            "  h3     two-level jump estimates dk2, de2 at s*\n\n"
            "options:\n"
            "  --alpha <a,b,...>  thresholds for cross (default 0.25,0.5,1)\n"
            "  --strict <file>    strict table (east_ed strict) for h3 ratios; may be repeated\n"
            "  --out <file>       output file (default stdout)\n\n"
            "Rows with edge = 1, ok = 0 or s* = nan are treated as missing. A point\n"
            "present in several tables is taken from the first one given.\n");
}

int main(int argc, char **argv)
{
    const char *mode, *out = NULL, *strict[16];
    double thr[MAXV] = {0.25, 0.5, 1.0};
    int nthr = 3, i, nfile = 0, nstrict = 0;
    FILE *f = stdout;
    if (argc < 3) {
        usage();
        return 1;
    }
    mode = argv[1];
    if (strcmp(mode, "h1") && strcmp(mode, "h2") && strcmp(mode, "h3") && strcmp(mode, "cross")) {
        fprintf(stderr, "error: unknown mode '%s'\n", mode);
        usage();
        return 1;
    }
    for (i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--out") && i + 1 < argc)
            out = argv[++i];
        else if (!strcmp(argv[i], "--strict") && i + 1 < argc) {
            if (nstrict == 16) {
                fprintf(stderr, "error: too many --strict tables\n");
                return 1;
            }
            strict[nstrict++] = argv[++i];
        }
        else if (!strcmp(argv[i], "--alpha") && i + 1 < argc) {
            char *p = argv[++i], *q;
            nthr = 0;
            while (*p && nthr < MAXV) {
                thr[nthr] = strtod(p, &q);
                if (q == p) {
                    fprintf(stderr, "error: bad --alpha list '%s'\n", argv[i]);
                    return 1;
                }
                nthr++;
                p = (*q == ',') ? q + 1 : q;
            }
        } else if (argv[i][0] == '-' && argv[i][1] == '-') {
            fprintf(stderr, "error: unknown option %s\n", argv[i]);
            usage();
            return 1;
        } else {
            if (read_peak(argv[i]) != 0) return 1;
            nfile++;
        }
    }
    if (nfile == 0 || nrow == 0) {
        fprintf(stderr, "error: no peak table rows\n");
        return 1;
    }
    if (isnan(c_val)) {
        fprintf(stderr, "error: no '# T = .., c = ..' header line found\n");
        return 1;
    }
    for (i = 0; i < nstrict; i++)
        if (read_strict(strict[i]) != 0) return 1;
    build_lists();
    if (out) {
        f = fopen(out, "w");
        if (!f) {
            fprintf(stderr, "error: cannot write %s\n", out);
            return 1;
        }
    }
    fprintf(f, "# east_edan %s, c = %.11f, %d rows,", mode, c_val, nrow);
    fprintf(f, " L =");
    for (i = 0; i < nL; i++) fprintf(f, " %d", Ls[i]);
    fprintf(f, "\n# input:");
    for (i = 2; i < argc; i++) fprintf(f, " %s", argv[i]);
    fprintf(f, "\n");
    if (!strcmp(mode, "h2"))
        run_h2(f);
    else if (!strcmp(mode, "h1"))
        run_h1(f);
    else if (!strcmp(mode, "cross"))
        run_cross(f, thr, nthr);
    else
        run_h3(f);
    if (f != stdout) fclose(f);
    return 0;
}
