#include "peaks.h"

#include <math.h>
#include <stdlib.h>

int smooth_halfwidth(int nbins, double smooth_res)
{
    if (!(smooth_res > 0.0)) return 0;
    return (int)floor(nbins / (2.0 * smooth_res) + 0.5);
}

int split_two_peaks(const double *h, int nb, int half)
{
    double *sm;
    int i, j, gmax = 0, best = -1, valley = -1;
    double best_depth = 0.0;
    if (nb < 3) return -1;
    sm = (double *)malloc(sizeof(double) * (size_t)nb);
    if (!sm) return -1;
    for (i = 0; i < nb; i++) {
        double acc = 0.0;
        int c = 0;
        for (j = i - half; j <= i + half; j++)
            if (j >= 0 && j < nb) {
                acc += h[j];
                c++;
            }
        sm[i] = acc / c;
    }
    for (i = 1; i < nb; i++)
        if (sm[i] > sm[gmax]) gmax = i;
    for (i = 0; i < nb; i++) {
        int is_max;
        double lo, vmin, depth;
        int a, b, k, vk;
        if (i == gmax) continue;
        is_max = (i == 0 || sm[i] > sm[i - 1]) && (i == nb - 1 || sm[i] >= sm[i + 1]);
        if (!is_max) continue;
        a = i < gmax ? i : gmax;
        b = i < gmax ? gmax : i;
        vmin = sm[a];
        vk = a;
        for (k = a; k <= b; k++)
            if (sm[k] < vmin) {
                vmin = sm[k];
                vk = k;
            }
        lo = sm[i] < sm[gmax] ? sm[i] : sm[gmax];
        depth = lo - vmin;
        if (depth > best_depth && depth >= 0.1 * lo) {
            best_depth = depth;
            best = i;
            valley = vk;
        }
    }
    free(sm);
    return (best >= 0) ? valley : -1;
}

double equal_weight_fixed(const SplitProblem *sp, double Kcut, double slo, double shi)
{
    double plo = sp->prob_below(slo, Kcut, sp->ctx) - 0.5;
    double phi = sp->prob_below(shi, Kcut, sp->ctx) - 0.5;
    int it;
    if (plo * phi > 0.0) return NAN;
    for (it = 0; it < 60 && shi - slo > 1e-8; it++) {
        double mid = 0.5 * (slo + shi);
        double pm = sp->prob_below(mid, Kcut, sp->ctx) - 0.5;
        if ((pm < 0.0) == (plo < 0.0)) {
            slo = mid;
            plo = pm;
        } else {
            shi = mid;
        }
    }
    return 0.5 * (slo + shi);
}

double equal_weight_auto(const SplitProblem *sp, double s0, double slo, double shi,
                         double *Kcut_out)
{
    double *h = (double *)malloc(sizeof(double) * (size_t)sp->nbins);
    const double w = (sp->Khi - sp->Klo) / sp->nbins;
    double s = s0, s_prev = NAN, Kcut = NAN, result = NAN;
    int it, b_prev = -2, b_prev2 = -3;
    if (!h) return NAN;
    for (it = 0; it < 20; it++) {
        int b;
        double snew;
        sp->hist(s, h, sp->ctx);
        b = split_two_peaks(h, sp->nbins, smooth_halfwidth(sp->nbins, sp->smooth_res));
        if (b < 0) break;
        Kcut = sp->Klo + (b + 0.5) * w;
        snew = equal_weight_fixed(sp, Kcut, slo, shi);
        if (isnan(snew)) break;
        if (b == b_prev || fabs(snew - s) < 1e-8) { /* settled */
            result = snew;
            break;
        }
        if (b == b_prev2) { /* the valley alternates between two bins: take the middle */
            result = 0.5 * (snew + s_prev);
            Kcut = sp->Klo + (0.5 * (b + b_prev) + 0.5) * w;
            break;
        }
        b_prev2 = b_prev;
        b_prev = b;
        s_prev = s;
        s = snew;
    }
    free(h);
    if (!isnan(result)) *Kcut_out = Kcut;
    return result;
}
