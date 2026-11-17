#include "tgen.h"

#include <math.h>
#include <stdlib.h>

int tgen_build(TGen *t, const Model *m, const Basis *b)
{
    const int L = m->L;
    const double amp0 = sqrt(m->c * (1.0 - m->c));
    long n = b->dim;
    long j, cap, nnz = 0;
    long *tmp_col = NULL;
    double *tmp_val = NULL;
    int k;

    t->dim = n;
    t->nnz = 0;
    t->rowptr = (long *)malloc(sizeof(long) * (size_t)(n + 1));
    t->esc = (double *)malloc(sizeof(double) * (size_t)n);
    t->occ = (double *)malloc(sizeof(double) * (size_t)n);
    cap = n * (long)L;
    t->col = (long *)malloc(sizeof(long) * (size_t)(cap > 0 ? cap : 1));
    t->val = (double *)malloc(sizeof(double) * (size_t)(cap > 0 ? cap : 1));
    tmp_col = (long *)malloc(sizeof(long) * (size_t)L);
    tmp_val = (double *)malloc(sizeof(double) * (size_t)L);
    if (!t->rowptr || !t->esc || !t->occ || !t->col || !t->val || !tmp_col || !tmp_val) {
        free(tmp_col);
        free(tmp_val);
        tgen_free(t);
        return -1;
    }

    /*
     * Row j of A in the symmetry-adapted basis:
     *   A[j', j] = sqrt(p_j / p_j') * sum_{x' in orbit(j')} A[x', rep_j]
     * Because A is symmetric, filling "column j" gives row j as well.
     */
    for (j = 0; j < n; j++) {
        uint32_t x = b->rep[j];
        int nt = 0;
        t->rowptr[j] = nnz;
        t->esc[j] = model_escape_rate(m, x);
        t->occ[j] = (double)conf_popcount(x);
        for (k = 0; k < L; k++) {
            double a = amp0 * model_constraint(m, x, k);
            uint32_t xp;
            long jp;
            int q, found = 0;
            if (a == 0.0) continue;
            xp = x ^ (1u << k);
            jp = b->index[xp];
            if (jp < 0) continue; /* only the excluded empty state */
            a *= sqrt((double)b->period[j] / (double)b->period[jp]);
            for (q = 0; q < nt; q++) {
                if (tmp_col[q] == jp) {
                    tmp_val[q] += a;
                    found = 1;
                    break;
                }
            }
            if (!found) {
                tmp_col[nt] = jp;
                tmp_val[nt] = a;
                nt++;
            }
        }
        for (k = 0; k < nt; k++) {
            t->col[nnz] = tmp_col[k];
            t->val[nnz] = tmp_val[k];
            nnz++;
        }
    }
    t->rowptr[n] = nnz;
    t->nnz = nnz;
    free(tmp_col);
    free(tmp_val);
    return 0;
}

void tgen_free(TGen *t)
{
    free(t->rowptr);
    free(t->col);
    free(t->val);
    free(t->esc);
    free(t->occ);
    t->rowptr = NULL;
    t->col = NULL;
    t->val = NULL;
    t->esc = NULL;
    t->occ = NULL;
    t->dim = 0;
    t->nnz = 0;
}

void tgen_apply_A(const TGen *t, const double *x, double *y)
{
    long i, p;
    for (i = 0; i < t->dim; i++) {
        double acc = 0.0;
        for (p = t->rowptr[i]; p < t->rowptr[i + 1]; p++) acc += t->val[p] * x[t->col[p]];
        y[i] = acc;
    }
}

void tgen_apply(const TGen *t, double s, double g, const double *x, double *y)
{
    const double es = exp(-s);
    long i, p;
    for (i = 0; i < t->dim; i++) {
        double acc = 0.0;
        for (p = t->rowptr[i]; p < t->rowptr[i + 1]; p++) acc += t->val[p] * x[t->col[p]];
        y[i] = es * acc - (t->esc[i] + g * t->occ[i]) * x[i];
    }
}

double tgen_asymmetry(const TGen *t)
{
    long i, p, q;
    double worst = 0.0;
    for (i = 0; i < t->dim; i++) {
        for (p = t->rowptr[i]; p < t->rowptr[i + 1]; p++) {
            long j = t->col[p];
            double aji = 0.0;
            double d;
            for (q = t->rowptr[j]; q < t->rowptr[j + 1]; q++) {
                if (t->col[q] == i) aji += t->val[q];
            }
            d = fabs(t->val[p] - aji);
            if (d > worst) worst = d;
        }
    }
    return worst;
}
