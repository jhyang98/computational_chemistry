#include "clone.h"
#include "dyn.h"
#include "rng.h"
#include "stats.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void clone_default_params(CloneParams *p)
{
    memset(p, 0, sizeof(*p));
    p->L = 100;
    p->T = 1.0;
    p->n_clones = 1000;
    p->tobs = 600;
    p->t_skip = 100;
    p->t_tail = 100;
    p->n_runs = 10;
}

typedef struct {
    unsigned char *x; /* N * L */
    int *n;
    long *Kh, *Eh;
} Population;

static int pop_alloc(Population *P, int N, int L)
{
    P->x = (unsigned char *)malloc((size_t)N * (size_t)L);
    P->n = (int *)malloc(sizeof(int) * (size_t)N);
    P->Kh = (long *)malloc(sizeof(long) * (size_t)N);
    P->Eh = (long *)malloc(sizeof(long) * (size_t)N);
    return (P->x && P->n && P->Kh && P->Eh) ? 0 : -1;
}

static void pop_free(Population *P)
{
    free(P->x);
    free(P->n);
    free(P->Kh);
    free(P->Eh);
}

static void mean_sem(const double *a, int n, double *mean, double *sem)
{
    double m = 0.0, v = 0.0;
    int i;
    for (i = 0; i < n; i++) m += a[i];
    m /= n;
    for (i = 0; i < n; i++) v += (a[i] - m) * (a[i] - m);
    *mean = m;
    *sem = (n > 1) ? sqrt(v / (n - 1) / n) : NAN;
}

static int one_run(const CloneParams *p, const Dyn *d, Rng *r, double *theta, double *k,
                   double *e)
{
    const int N = p->n_clones, L = p->L;
    Population A, B;
    double *lw = (double *)malloc(sizeof(double) * (size_t)N);
    double *cw = (double *)malloc(sizeof(double) * (size_t)N);
    int *kk = (int *)malloc(sizeof(int) * (size_t)N);
    double lnZ = 0.0;
    int j, t;

    if (!lw || !cw || !kk || pop_alloc(&A, N, L) || pop_alloc(&B, N, L)) {
        free(lw);
        free(cw);
        free(kk);
        return -1;
    }
    for (j = 0; j < N; j++) {
        A.n[j] = dyn_sample_equilibrium(d, A.x + (size_t)j * L, r);
        A.Kh[j] = 0;
        A.Eh[j] = 0;
    }

    for (t = 1; t <= p->tobs; t++) {
        double lmax, total, u;
        int m, src;
        for (j = 0; j < N; j++) {
            kk[j] = dyn_mcs(d, A.x + (size_t)j * L, &A.n[j], r);
            lw[j] = -p->s * kk[j] - p->g * A.n[j];
            if (t > p->t_skip && t <= p->tobs - p->t_tail) {
                A.Kh[j] += kk[j];
                A.Eh[j] += A.n[j];
            }
        }
        lmax = lw[0];
        for (j = 1; j < N; j++)
            if (lw[j] > lmax) lmax = lw[j];
        total = 0.0;
        for (j = 0; j < N; j++) {
            total += exp(lw[j] - lmax);
            cw[j] = total;
        }
        if (t > p->t_skip) lnZ += lmax + log(total / N);

        /* systematic resampling */
        u = rng_uniform(r);
        src = 0;
        for (m = 0; m < N; m++) {
            double target = (u + m) / N * total;
            while (src < N - 1 && cw[src] < target) src++;
            memcpy(B.x + (size_t)m * L, A.x + (size_t)src * L, (size_t)L);
            B.n[m] = A.n[src];
            B.Kh[m] = A.Kh[src];
            B.Eh[m] = A.Eh[src];
        }
        {
            Population tmp = A;
            A = B;
            B = tmp;
        }
    }
    {
        double sk = 0.0, se = 0.0, tt = (double)(p->tobs - p->t_skip);
        double th = (double)(p->tobs - p->t_tail - p->t_skip);
        for (j = 0; j < N; j++) {
            sk += (double)A.Kh[j];
            se += (double)A.Eh[j];
        }
        *theta = lnZ / tt;
        *k = sk / N / (L * th);
        *e = se / N / (L * th);
    }
    pop_free(&A);
    pop_free(&B);
    free(lw);
    free(cw);
    free(kk);
    return 0;
}

int clone_run(const CloneParams *p, CloneResult *res, double *theta_runs, double *k_runs,
              double *e_runs)
{
    Dyn d;
    Rng r;
    double *th, *kv, *ev;
    int i;

    if (p->n_clones < 2 || p->t_skip < 0 || p->t_tail < 0 || p->tobs <= p->t_skip + p->t_tail ||
        p->n_runs < 1) {
        fprintf(stderr, "error: need n_clones >= 2, t_skip >= 0, t_tail >= 0, "
                        "t_skip + t_tail < tobs, n_runs >= 1\n");
        return -1;
    }
    if (dyn_init(&d, p->L, p->T, p->eps) != 0) {
        fprintf(stderr, "error: invalid L, T or eps\n");
        return -1;
    }
    res->seed_used = rng_seed_auto(&r, p->seed);
    th = (double *)malloc(sizeof(double) * (size_t)p->n_runs);
    kv = (double *)malloc(sizeof(double) * (size_t)p->n_runs);
    ev = (double *)malloc(sizeof(double) * (size_t)p->n_runs);
    if (!th || !kv || !ev) {
        free(th);
        free(kv);
        free(ev);
        return -1;
    }
    for (i = 0; i < p->n_runs; i++) {
        if (one_run(p, &d, &r, &th[i], &kv[i], &ev[i]) != 0) {
            free(th);
            free(kv);
            free(ev);
            return -1;
        }
        if (theta_runs) theta_runs[i] = th[i];
        if (k_runs) k_runs[i] = kv[i];
        if (e_runs) e_runs[i] = ev[i];
    }
    mean_sem(th, p->n_runs, &res->theta, &res->theta_err);
    mean_sem(kv, p->n_runs, &res->k, &res->k_err);
    mean_sem(ev, p->n_runs, &res->e, &res->e_err);
    free(th);
    free(kv);
    free(ev);
    return 0;
}
