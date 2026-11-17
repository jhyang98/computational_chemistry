#include "stats.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

typedef struct {
    double n, sx, sy, sxx, syy, sxy;
} Sums;

static void moments(const Sums *s, double shift_x, double shift_y, double *mx, double *my,
                    double *vx, double *vy, double *cxy)
{
    /* sums are of shifted data; shifts restore the means */
    double ax = s->sx / s->n, ay = s->sy / s->n;
    double f = s->n / (s->n - 1.0);
    *mx = ax + shift_x;
    *my = ay + shift_y;
    *vx = f * (s->sxx / s->n - ax * ax);
    *vy = f * (s->syy / s->n - ay * ay);
    *cxy = f * (s->sxy / s->n - ax * ay);
}

void pair_stats(const double *x, const double *y, long n, int nblocks, PairStats *o)
{
    Sums tot = {0, 0, 0, 0, 0, 0};
    Sums *blk;
    double sx0 = x[0], sy0 = y[0]; /* shift for numerical stability */
    long i, bsize;
    int b, B;
    double jm[5] = {0, 0, 0, 0, 0}, jv[5] = {0, 0, 0, 0, 0};
    double full[5];

    o->n = n;
    if (nblocks > n) nblocks = (int)n;
    if (nblocks < 2) nblocks = 2;
    bsize = n / nblocks;
    B = nblocks;

    blk = (Sums *)calloc((size_t)B, sizeof(Sums));
    for (i = 0; i < n; i++) {
        double dx = x[i] - sx0, dy = y[i] - sy0;
        long bi = bsize > 0 ? i / bsize : 0;
        tot.n += 1;
        tot.sx += dx;
        tot.sy += dy;
        tot.sxx += dx * dx;
        tot.syy += dy * dy;
        tot.sxy += dx * dy;
        if (bi < B) { /* trailing remainder is used in the totals only */
            blk[bi].n += 1;
            blk[bi].sx += dx;
            blk[bi].sy += dy;
            blk[bi].sxx += dx * dx;
            blk[bi].syy += dy * dy;
            blk[bi].sxy += dx * dy;
        }
    }
    moments(&tot, sx0, sy0, &full[0], &full[1], &full[2], &full[3], &full[4]);
    o->mean_x = full[0];
    o->mean_y = full[1];
    o->var_x = full[2];
    o->var_y = full[3];
    o->cov = full[4];

    /* leave-one-block-out jackknife over the blocked part */
    {
        Sums all = {0, 0, 0, 0, 0, 0};
        double (*est)[5] = (double (*)[5])malloc(sizeof(double) * 5 * (size_t)B);
        int q;
        for (b = 0; b < B; b++) {
            all.n += blk[b].n;
            all.sx += blk[b].sx;
            all.sy += blk[b].sy;
            all.sxx += blk[b].sxx;
            all.syy += blk[b].syy;
            all.sxy += blk[b].sxy;
        }
        for (b = 0; b < B; b++) {
            Sums s;
            s.n = all.n - blk[b].n;
            s.sx = all.sx - blk[b].sx;
            s.sy = all.sy - blk[b].sy;
            s.sxx = all.sxx - blk[b].sxx;
            s.syy = all.syy - blk[b].syy;
            s.sxy = all.sxy - blk[b].sxy;
            moments(&s, sx0, sy0, &est[b][0], &est[b][1], &est[b][2], &est[b][3], &est[b][4]);
            for (q = 0; q < 5; q++) jm[q] += est[b][q] / B;
        }
        for (b = 0; b < B; b++)
            for (q = 0; q < 5; q++) jv[q] += (est[b][q] - jm[q]) * (est[b][q] - jm[q]);
        for (q = 0; q < 5; q++) jv[q] = sqrt((B - 1.0) / B * jv[q]);
        free(est);
    }
    o->nblocks = B;
    o->err_mean_x = jv[0];
    o->err_mean_y = jv[1];
    o->err_var_x = jv[2];
    o->err_var_y = jv[3];
    o->err_cov = jv[4];
    free(blk);
}

double log_sum_exp(const double *a, long n)
{
    double m = -DBL_MAX, s = 0.0;
    long i;
    for (i = 0; i < n; i++)
        if (a[i] > m) m = a[i];
    if (m == -DBL_MAX) return -INFINITY;
    for (i = 0; i < n; i++) s += exp(a[i] - m);
    return m + log(s);
}
