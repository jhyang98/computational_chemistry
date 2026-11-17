/*
 * Means, variances and covariance of a correlated pair series (x_n, y_n)
 * with block-jackknife error bars.
 */
#ifndef STATS_H
#define STATS_H

typedef struct {
    long n;
    int nblocks;
    double mean_x, mean_y, var_x, var_y, cov;
    double err_mean_x, err_mean_y, err_var_x, err_var_y, err_cov;
} PairStats;

/* nblocks is reduced to n if n < nblocks; needs n >= 2. */
void pair_stats(const double *x, const double *y, long n, int nblocks, PairStats *o);

/* log(sum exp(a_i)) */
double log_sum_exp(const double *a, long n);

#endif
