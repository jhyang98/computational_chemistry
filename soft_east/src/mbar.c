#include "mbar.h"
#include "rng.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int mbar_alloc(MBAR *m, int nwin, long ntot)
{
    int k;
    memset(m, 0, sizeof(*m));
    m->nwin = nwin;
    m->ntot = ntot;
    m->ws = (double *)calloc((size_t)nwin, sizeof(double));
    m->wg = (double *)calloc((size_t)nwin, sizeof(double));
    m->Kmin = (double *)malloc(sizeof(double) * (size_t)nwin);
    m->Kmax = (double *)malloc(sizeof(double) * (size_t)nwin);
    m->N = (long *)calloc((size_t)nwin, sizeof(long));
    m->f = (double *)calloc((size_t)nwin, sizeof(double));
    m->K = (double *)calloc((size_t)(ntot > 0 ? ntot : 1), sizeof(double));
    m->E = (double *)calloc((size_t)(ntot > 0 ? ntot : 1), sizeof(double));
    m->logden = (double *)calloc((size_t)(ntot > 0 ? ntot : 1), sizeof(double));
    if (!m->ws || !m->wg || !m->Kmin || !m->Kmax || !m->N || !m->f || !m->K || !m->E ||
        !m->logden) {
        mbar_free(m);
        return -1;
    }
    for (k = 0; k < nwin; k++) {
        m->Kmin[k] = -INFINITY;
        m->Kmax[k] = INFINITY;
    }
    return 0;
}

void mbar_free(MBAR *m)
{
    free(m->ws);
    free(m->wg);
    free(m->Kmin);
    free(m->Kmax);
    free(m->N);
    free(m->f);
    free(m->K);
    free(m->E);
    free(m->logden);
    memset(m, 0, sizeof(*m));
}

static int inside(const MBAR *m, int k, long n)
{
    return m->K[n] >= m->Kmin[k] && m->K[n] <= m->Kmax[k];
}

/* logden_n and the objective F(f) */
static double objective(const MBAR *m, const double *f, double *logden)
{
    const int W = m->nwin;
    double F = 0.0;
    long n;
    int k;
    double *a = (double *)malloc(sizeof(double) * (size_t)W);
    for (n = 0; n < m->ntot; n++) {
        double mx = -DBL_MAX, sum = 0.0;
        for (k = 0; k < W; k++) {
            if (!inside(m, k, n)) {
                a[k] = -INFINITY;
                continue;
            }
            a[k] = log((double)m->N[k]) + f[k] - m->ws[k] * m->K[n] - m->wg[k] * m->E[n];
            if (a[k] > mx) mx = a[k];
        }
        for (k = 0; k < W; k++)
            if (a[k] > -INFINITY) sum += exp(a[k] - mx);
        logden[n] = mx + log(sum);
        F += logden[n];
    }
    for (k = 0; k < W; k++) F -= (double)m->N[k] * f[k];
    free(a);
    return F;
}

/* Gaussian elimination with partial pivoting; a is n x n row-major, b overwritten by x. */
static int solve_linear(int n, double *a, double *b)
{
    int i, j, k;
    for (k = 0; k < n; k++) {
        int p = k;
        double mx = fabs(a[k * n + k]);
        for (i = k + 1; i < n; i++)
            if (fabs(a[i * n + k]) > mx) {
                mx = fabs(a[i * n + k]);
                p = i;
            }
        if (mx == 0.0) return -1;
        if (p != k) {
            for (j = 0; j < n; j++) {
                double t = a[k * n + j];
                a[k * n + j] = a[p * n + j];
                a[p * n + j] = t;
            }
            {
                double t = b[k];
                b[k] = b[p];
                b[p] = t;
            }
        }
        for (i = k + 1; i < n; i++) {
            double fct = a[i * n + k] / a[k * n + k];
            for (j = k; j < n; j++) a[i * n + j] -= fct * a[k * n + j];
            b[i] -= fct * b[k];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (j = i + 1; j < n; j++) s -= a[i * n + j] * b[j];
        b[i] = s / a[i * n + i];
    }
    return 0;
}

void mbar_guess(MBAR *m)
{
    const int W = m->nwin;
    long *off = (long *)calloc((size_t)W, sizeof(long));
    double *mK = (double *)calloc((size_t)W, sizeof(double));
    double *mE = (double *)calloc((size_t)W, sizeof(double));
    int k;
    for (k = 1; k < W; k++) off[k] = off[k - 1] + m->N[k - 1];
    for (k = 0; k < W; k++) {
        long n;
        for (n = 0; n < m->N[k]; n++) {
            mK[k] += m->K[off[k] + n];
            mE[k] += m->E[off[k] + n];
        }
        mK[k] /= (double)m->N[k];
        mE[k] /= (double)m->N[k];
    }
    m->f[0] = 0.0;
    for (k = 1; k < W; k++) {
        const int j = k - 1;
        double lo = m->Kmin[j] > m->Kmin[k] ? m->Kmin[j] : m->Kmin[k];
        double hi = m->Kmax[j] < m->Kmax[k] ? m->Kmax[j] : m->Kmax[k];
        int windowed = isfinite(m->Kmin[k]) || isfinite(m->Kmax[k]) || isfinite(m->Kmin[j]) ||
                       isfinite(m->Kmax[j]);
        if (windowed && m->ws[j] == m->ws[k] && m->wg[j] == m->wg[k] && lo <= hi) {
            /* Z_k / Z_j = (fraction of j's samples in the overlap) / (fraction of k's) */
            long n, cj = 0, ck = 0;
            for (n = 0; n < m->N[j]; n++)
                if (m->K[off[j] + n] >= lo && m->K[off[j] + n] <= hi) cj++;
            for (n = 0; n < m->N[k]; n++)
                if (m->K[off[k] + n] >= lo && m->K[off[k] + n] <= hi) ck++;
            if (cj > 0 && ck > 0) {
                double a = (double)cj / m->N[j], b = (double)ck / m->N[k];
                m->f[k] = m->f[j] - log(a / b);
                continue;
            }
            m->f[k] = m->f[j];
            continue;
        }
        /* f = -ln Z + const, df = <K> ds + <E> dg */
        m->f[k] = m->f[j] + 0.5 * (mK[k] + mK[j]) * (m->ws[k] - m->ws[j]) +
                  0.5 * (mE[k] + mE[j]) * (m->wg[k] - m->wg[j]);
    }
    free(off);
    free(mK);
    free(mE);
}

int mbar_solve(MBAR *m, double tol, int maxit)
{
    const int W = m->nwin, D = W - 1;
    double *grad = (double *)malloc(sizeof(double) * (size_t)W);
    double *H = (double *)malloc(sizeof(double) * (size_t)(D > 0 ? D * D : 1));
    double *step = (double *)malloc(sizeof(double) * (size_t)(D > 0 ? D : 1));
    double *ftry = (double *)malloc(sizeof(double) * (size_t)W);
    double *ldtry = (double *)malloc(sizeof(double) * (size_t)m->ntot);
    double *sw = (double *)malloc(sizeof(double) * (size_t)W);
    double *wn = (double *)malloc(sizeof(double) * (size_t)W);
    int *idx = (int *)malloc(sizeof(int) * (size_t)W);
    double F;
    int it, k, l, rc = -1;

    if (!grad || !H || !step || !ftry || !ldtry || !sw || !wn || !idx) goto done;
    m->f[0] = 0.0;
    F = objective(m, m->f, m->logden);
    if (W == 1) {
        rc = 0;
        m->iters = 0;
        m->grad_max = 0.0;
        goto done;
    }
    for (it = 0; it < maxit; it++) {
        long n;
        double gmax = 0.0, slope = 0.0, alpha = 1.0;
        memset(sw, 0, sizeof(double) * (size_t)W);
        memset(H, 0, sizeof(double) * (size_t)(D * D));
        for (n = 0; n < m->ntot; n++) {
            int nz = 0, a, b;
            for (k = 0; k < W; k++) {
                if (!inside(m, k, n)) continue;
                wn[k] = exp(m->f[k] - m->ws[k] * m->K[n] - m->wg[k] * m->E[n] - m->logden[n]);
                sw[k] += wn[k];
                idx[nz++] = k;
            }
            for (a = 0; a < nz; a++) {
                k = idx[a];
                if (k == 0) continue;
                for (b = 0; b < nz; b++) {
                    l = idx[b];
                    if (l == 0) continue;
                    H[(k - 1) * D + (l - 1)] -= (double)m->N[k] * (double)m->N[l] * wn[k] * wn[l];
                }
            }
        }
        for (k = 1; k < W; k++) {
            grad[k] = (double)m->N[k] * sw[k] - (double)m->N[k];
            H[(k - 1) * D + (k - 1)] += (double)m->N[k] * sw[k];
            if (fabs(grad[k]) / (double)m->N[k] > gmax) gmax = fabs(grad[k]) / (double)m->N[k];
        }
        m->grad_max = gmax;
        m->iters = it;
        if (gmax < tol) {
            rc = 0;
            break;
        }
        for (k = 1; k < W; k++) step[k - 1] = -grad[k];
        if (solve_linear(D, H, step) != 0) break;
        for (k = 1; k < W; k++) slope += grad[k] * step[k - 1];
        if (slope >= 0.0) { /* not a descent direction: fall back to gradient */
            for (k = 1; k < W; k++) step[k - 1] = -grad[k] / (double)m->N[k];
            slope = 0.0;
            for (k = 1; k < W; k++) slope += grad[k] * step[k - 1];
        }
        for (;;) {
            double Ft;
            ftry[0] = 0.0;
            for (k = 1; k < W; k++) ftry[k] = m->f[k] + alpha * step[k - 1];
            Ft = objective(m, ftry, ldtry);
            /* Near the minimum the decrease of F is below its round-off
               (F is a sum over all samples), so take full Newton steps there. */
            if (gmax < 1e-4 || Ft <= F + 1e-4 * alpha * slope || alpha < 1e-12) {
                memcpy(m->f, ftry, sizeof(double) * (size_t)W);
                memcpy(m->logden, ldtry, sizeof(double) * (size_t)m->ntot);
                F = Ft;
                break;
            }
            alpha *= 0.5;
        }
    }
done:
    free(grad);
    free(H);
    free(step);
    free(ftry);
    free(ldtry);
    free(sw);
    free(wn);
    free(idx);
    return rc;
}

/* normalized log-weights at (s, g) restricted to [Kmin, Kmax] (-inf outside) */
static void target_logw(const MBAR *m, double s, double g, double Kmin, double Kmax, double *lw)
{
    double mx = -DBL_MAX, sum = 0.0;
    long n;
    for (n = 0; n < m->ntot; n++) {
        if (m->K[n] < Kmin || m->K[n] > Kmax) {
            lw[n] = -INFINITY;
            continue;
        }
        lw[n] = -s * m->K[n] - g * m->E[n] - m->logden[n];
        if (lw[n] > mx) mx = lw[n];
    }
    for (n = 0; n < m->ntot; n++)
        if (lw[n] > -INFINITY) sum += exp(lw[n] - mx);
    for (n = 0; n < m->ntot; n++)
        if (lw[n] > -INFINITY) lw[n] -= mx + log(sum);
}

void mbar_observe_window(const MBAR *m, double s, double g, double Kmin, double Kmax,
                         MBARObs *o)
{
    double *lw = (double *)malloc(sizeof(double) * (size_t)m->ntot);
    double mK = 0, mE = 0, vK = 0, vE = 0, cKE = 0, sw2 = 0;
    long n;
    target_logw(m, s, g, Kmin, Kmax, lw);
    for (n = 0; n < m->ntot; n++) {
        double w = exp(lw[n]);
        mK += w * m->K[n];
        mE += w * m->E[n];
        sw2 += w * w;
    }
    for (n = 0; n < m->ntot; n++) {
        double w = exp(lw[n]), dK = m->K[n] - mK, dE = m->E[n] - mE;
        vK += w * dK * dK;
        vE += w * dE * dE;
        cKE += w * dK * dE;
    }
    o->meanK = mK;
    o->meanE = mE;
    o->varK = vK;
    o->varE = vE;
    o->covKE = cKE;
    o->ess = 1.0 / sw2;
    free(lw);
}

void mbar_observe(const MBAR *m, double s, double g, MBARObs *o)
{
    mbar_observe_window(m, s, g, -INFINITY, INFINITY, o);
}

void mbar_observe_split(const MBAR *m, double s, double g, double Kcut, MBARSplit *o)
{
    double *lw = (double *)malloc(sizeof(double) * (size_t)m->ntot);
    double w0 = 0, k0 = 0, e0 = 0, v0 = 0, w1 = 0, k1 = 0, e1 = 0, v1 = 0;
    long n;
    target_logw(m, s, g, -INFINITY, INFINITY, lw);
    for (n = 0; n < m->ntot; n++) {
        double w = exp(lw[n]);
        if (m->K[n] < Kcut) {
            w0 += w;
            k0 += w * m->K[n];
            e0 += w * m->E[n];
        } else {
            w1 += w;
            k1 += w * m->K[n];
            e1 += w * m->E[n];
        }
    }
    k0 = w0 > 0.0 ? k0 / w0 : NAN;
    e0 = w0 > 0.0 ? e0 / w0 : NAN;
    k1 = w1 > 0.0 ? k1 / w1 : NAN;
    e1 = w1 > 0.0 ? e1 / w1 : NAN;
    for (n = 0; n < m->ntot; n++) {
        double w = exp(lw[n]);
        if (m->K[n] < Kcut)
            v0 += w * (m->K[n] - k0) * (m->K[n] - k0);
        else
            v1 += w * (m->K[n] - k1) * (m->K[n] - k1);
    }
    o->p_lo = w0 / (w0 + w1);
    o->meanK_lo = k0;
    o->meanE_lo = e0;
    o->varK_lo = w0 > 0.0 ? v0 / w0 : NAN;
    o->meanK_hi = k1;
    o->meanE_hi = e1;
    o->varK_hi = w1 > 0.0 ? v1 / w1 : NAN;
    free(lw);
}

double mbar_prob_below(const MBAR *m, double s, double g, double Kcut)
{
    double *lw = (double *)malloc(sizeof(double) * (size_t)m->ntot);
    double p = 0.0;
    long n;
    target_logw(m, s, g, -INFINITY, INFINITY, lw);
    for (n = 0; n < m->ntot; n++)
        if (m->K[n] < Kcut) p += exp(lw[n]);
    free(lw);
    return p;
}

void mbar_histogram(const MBAR *m, double s, double g, double klo, double khi, int nbins,
                    double *h)
{
    double *lw = (double *)malloc(sizeof(double) * (size_t)m->ntot);
    double width = (khi - klo) / nbins;
    long n;
    int b;
    for (b = 0; b < nbins; b++) h[b] = 0.0;
    target_logw(m, s, g, -INFINITY, INFINITY, lw);
    for (n = 0; n < m->ntot; n++) {
        int bi = (int)floor((m->K[n] - klo) / width);
        if (bi >= 0 && bi < nbins) h[bi] += exp(lw[n]);
    }
    for (b = 0; b < nbins; b++) h[b] /= width; /* probability density in K */
    free(lw);
}

int mbar_resample(const MBAR *src, MBAR *dst, int block, Rng *rng)
{
    long off = 0;
    int k;
    if (mbar_alloc(dst, src->nwin, src->ntot) != 0) return -1;
    for (k = 0; k < src->nwin; k++) {
        long n = src->N[k], filled = 0;
        long bl = (block < n) ? block : n;
        if (bl < 1) bl = 1;
        dst->ws[k] = src->ws[k];
        dst->wg[k] = src->wg[k];
        dst->Kmin[k] = src->Kmin[k];
        dst->Kmax[k] = src->Kmax[k];
        dst->N[k] = n;
        dst->f[k] = src->f[k];
        while (filled < n) {
            long start = (long)(rng_uniform(rng) * (double)(n - bl + 1)), q;
            for (q = 0; q < bl && filled < n; q++, filled++) {
                dst->K[off + filled] = src->K[off + start + q];
                dst->E[off + filled] = src->E[off + start + q];
            }
        }
        off += n;
    }
    return 0;
}
