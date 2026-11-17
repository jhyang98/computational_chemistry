/*
 * Two-peak analysis of an activity distribution P_s(K) and the equal-weight
 * coexistence point.
 *
 * The distribution is supplied through callbacks so that the same code works
 * for exact P(K) (east_tm) and for MBAR-reweighted samples (east_mbar).
 *
 * Valley: the histogram (nbins bins on [Klo, Khi)) is smoothed with a moving
 * average over 2h+1 bins, h = round(nbins / (2 smooth_res)), so the averaging
 * window spans about (Khi - Klo) / smooth_res whatever the number of bins
 * (smooth_res <= 0: no smoothing). The highest local maximum is taken as
 * the first peak; the second peak is the local maximum j that maximizes the
 * depth of the valley between them, depth = min(h_j, h_max) - min(h between).
 * A split is accepted only if depth >= 0.1 * min(h_j, h_max). The valley is
 * the lowest smoothed bin between the two peaks, and Kcut its bin centre.
 *
 * Equal weight with a fixed Kcut: the s in [slo, shi] where
 * P_s(K < Kcut) = 1/2 (bisection; P_s(K < Kcut) increases with s).
 *
 * Equal weight with an automatic Kcut: start at s0, find the valley of P_s0,
 * solve the fixed-Kcut problem, recompute the valley at the new s, and repeat
 * until s and Kcut no longer change (if the valley alternates between two
 * neighbouring bins, the midpoint of the two solutions is returned).
 */
#ifndef PEAKS_H
#define PEAKS_H

typedef struct {
    double (*prob_below)(double s, double Kcut, void *ctx); /* P_s(K < Kcut) */
    void (*hist)(double s, double *h, void *ctx);           /* binned P_s on [Klo, Khi) */
    void *ctx;
    double Klo, Khi;
    int nbins;
    double smooth_res; /* smoothing resolution m: window ~ (Khi - Klo) / m */
} SplitProblem;

/* Half-width h = round(nbins / (2 m)) of the moving average (0 if m <= 0). */
int smooth_halfwidth(int nbins, double smooth_res);

/* Valley bin of h smoothed with half-width `half` bins, or -1 if no acceptable second peak. */
int split_two_peaks(const double *h, int nb, int half);

/* NAN if P_s(K < Kcut) - 1/2 does not change sign on [slo, shi]. */
double equal_weight_fixed(const SplitProblem *sp, double Kcut, double slo, double shi);

/* NAN if no valley is found or the iteration does not settle; *Kcut_out set on success. */
double equal_weight_auto(const SplitProblem *sp, double s0, double slo, double shi,
                         double *Kcut_out);

#endif
