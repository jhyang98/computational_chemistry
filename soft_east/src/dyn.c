#include "dyn.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int dyn_init(Dyn *d, int L, double T, double eps)
{
    if (L < 2 || !(T > 0.0) || eps < 0.0 || eps > 1.0) return -1;
    d->L = L;
    d->T = T;
    d->c = 1.0 / (1.0 + exp(1.0 / T));
    d->eps = eps;
    d->pflip[1][0] = d->c;
    d->pflip[1][1] = 1.0 - d->c;
    d->pflip[0][0] = eps * d->c;
    d->pflip[0][1] = eps * (1.0 - d->c);
    return 0;
}

int dyn_sample_equilibrium(const Dyn *d, unsigned char *x, Rng *r)
{
    int i, n = 0;
    for (i = 0; i < d->L; i++) {
        x[i] = (unsigned char)(rng_uniform(r) < d->c);
        n += x[i];
    }
    return n;
}

int dyn_mcs(const Dyn *d, unsigned char *x, int *n, Rng *r)
{
    const int L = d->L;
    int a, flips = 0;
    for (a = 0; a < L; a++) {
        int i = rng_int(r, L);
        int left = x[i == 0 ? L - 1 : i - 1];
        double p = d->pflip[left][x[i]];
        if (p > 0.0 && rng_uniform(r) < p) {
            *n += x[i] ? -1 : 1;
            x[i] ^= 1u;
            flips++;
        }
    }
    return flips;
}

void dyn_propagate(const Dyn *d, const unsigned char *x_start, int nsteps, unsigned char *out,
                   int *n_out, int *k_out, Rng *r)
{
    const int L = d->L;
    const unsigned char *prev = x_start;
    int j, i, n = 0;
    for (i = 0; i < L; i++) n += x_start[i];
    for (j = 0; j < nsteps; j++) {
        unsigned char *cur = out + (size_t)j * (size_t)L;
        memcpy(cur, prev, (size_t)L);
        k_out[j] = dyn_mcs(d, cur, &n, r);
        n_out[j] = n;
        prev = cur;
    }
}

int traj_alloc(Traj *tr, int L, int tobs)
{
    tr->L = L;
    tr->tobs = tobs;
    tr->x = (unsigned char *)calloc((size_t)(tobs + 1) * (size_t)L, 1);
    tr->n = (int *)calloc((size_t)(tobs + 1), sizeof(int));
    tr->k = (int *)calloc((size_t)(tobs + 1), sizeof(int));
    tr->K = tr->E = 0;
    if (!tr->x || !tr->n || !tr->k) {
        traj_free(tr);
        return -1;
    }
    return 0;
}

void traj_free(Traj *tr)
{
    free(tr->x);
    free(tr->n);
    free(tr->k);
    tr->x = NULL;
    tr->n = tr->k = NULL;
}

void traj_recount(Traj *tr)
{
    int t, i;
    tr->K = 0;
    tr->E = 0;
    tr->k[0] = 0;
    for (t = 0; t <= tr->tobs; t++) {
        const unsigned char *s = tr->x + (size_t)t * (size_t)tr->L;
        int n = 0;
        for (i = 0; i < tr->L; i++) n += s[i];
        tr->n[t] = n;
        tr->E += n;
        if (t > 0) tr->K += tr->k[t];
    }
}

void traj_init_equilibrium(Traj *tr, const Dyn *d, Rng *r)
{
    dyn_sample_equilibrium(d, tr->x, r);
    if (tr->tobs > 0)
        dyn_propagate(d, tr->x, tr->tobs, tr->x + tr->L, tr->n + 1, tr->k + 1, r);
    traj_recount(tr);
}

int traj_save(const Traj *tr, const char *path)
{
    FILE *f = fopen(path, "w");
    int t, i;
    char *line;
    if (!f) {
        fprintf(stderr, "error: cannot write trajectory file '%s'\n", path);
        return -1;
    }
    line = (char *)malloc((size_t)tr->L + 1);
    if (!line) {
        fclose(f);
        return -1;
    }
    fprintf(f, "# soft_east trajectory: one line per slice, 'k_t config'\n");
    fprintf(f, "# k_t = flips during sweep t-1 -> t (0 for t = 0)\n");
    fprintf(f, "L %d\ntobs %d\nK %ld\nE %ld\n", tr->L, tr->tobs, tr->K, tr->E);
    for (t = 0; t <= tr->tobs; t++) {
        const unsigned char *s = tr->x + (size_t)t * (size_t)tr->L;
        for (i = 0; i < tr->L; i++) line[i] = (char)('0' + s[i]);
        line[tr->L] = '\0';
        fprintf(f, "%d %s\n", tr->k[t], line);
    }
    free(line);
    fclose(f);
    return 0;
}

int traj_load(Traj *tr, const char *path)
{
    FILE *f = fopen(path, "r");
    char buf[4096];
    int L = -1, tobs = -1, t = 0, i;
    long K = -1, E = -1;
    if (!f) {
        fprintf(stderr, "error: cannot open trajectory file '%s'\n", path);
        return -1;
    }
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '#') continue;
        if (!strncmp(buf, "L ", 2)) {
            L = atoi(buf + 2);
            continue;
        }
        if (!strncmp(buf, "tobs ", 5)) {
            tobs = atoi(buf + 5);
            continue;
        }
        if (!strncmp(buf, "K ", 2)) {
            K = atol(buf + 2);
            continue;
        }
        if (!strncmp(buf, "E ", 2)) {
            E = atol(buf + 2);
            continue;
        }
        if (L != tr->L || tobs != tr->tobs) {
            fprintf(stderr, "error: trajectory file '%s' has L=%d tobs=%d, expected L=%d tobs=%d\n",
                    path, L, tobs, tr->L, tr->tobs);
            fclose(f);
            return -1;
        }
        {
            int kt;
            char *p = buf;
            char *sp = strchr(p, ' ');
            if (!sp || t > tobs) {
                fprintf(stderr, "error: malformed trajectory file '%s'\n", path);
                fclose(f);
                return -1;
            }
            kt = atoi(p);
            sp++;
            for (i = 0; i < L; i++) {
                if (sp[i] != '0' && sp[i] != '1') {
                    fprintf(stderr, "error: malformed slice %d in '%s'\n", t, path);
                    fclose(f);
                    return -1;
                }
                tr->x[(size_t)t * (size_t)L + (size_t)i] = (unsigned char)(sp[i] - '0');
            }
            tr->k[t] = (t == 0) ? 0 : kt;
            t++;
        }
    }
    fclose(f);
    if (t != tobs + 1) {
        fprintf(stderr, "error: trajectory file '%s' has %d slices, expected %d\n", path, t,
                tobs + 1);
        return -1;
    }
    traj_recount(tr);
    if ((K >= 0 && K != tr->K) || (E >= 0 && E != tr->E)) {
        fprintf(stderr, "error: trajectory file '%s' totals do not match its slices\n", path);
        return -1;
    }
    return 0;
}

int shot_alloc(Shot *sh, int L, int tobs)
{
    int m = tobs > 0 ? tobs : 1;
    sh->buf = (unsigned char *)malloc((size_t)m * (size_t)L);
    sh->nbuf = (int *)malloc(sizeof(int) * (size_t)m);
    sh->kbuf = (int *)malloc(sizeof(int) * (size_t)m);
    if (!sh->buf || !sh->nbuf || !sh->kbuf) {
        shot_free(sh);
        return -1;
    }
    return 0;
}

void shot_free(Shot *sh)
{
    free(sh->buf);
    free(sh->nbuf);
    free(sh->kbuf);
    sh->buf = NULL;
    sh->nbuf = sh->kbuf = NULL;
}

void shot_propose(Shot *sh, const Traj *tr, const Dyn *d, int dir, int tau, Rng *r)
{
    const int L = tr->L, tobs = tr->tobs;
    long oldK = 0, oldE = 0, newK = 0, newE = 0;
    int t, j;
    sh->dir = dir;
    sh->tau = tau;
    if (dir > 0) {
        /* forward: regenerate slices tau+1 .. tobs */
        sh->len = tobs - tau;
        dyn_propagate(d, tr->x + (size_t)tau * (size_t)L, sh->len, sh->buf, sh->nbuf, sh->kbuf,
                      r);
        for (t = tau + 1; t <= tobs; t++) {
            oldE += tr->n[t];
            oldK += tr->k[t];
        }
    } else {
        /* backward: regenerate slices tau-1 .. 0 (buffer slot j holds slice tau-1-j) */
        sh->len = tau;
        dyn_propagate(d, tr->x + (size_t)tau * (size_t)L, sh->len, sh->buf, sh->nbuf, sh->kbuf,
                      r);
        for (t = 0; t < tau; t++) oldE += tr->n[t];
        for (t = 1; t <= tau; t++) oldK += tr->k[t];
    }
    for (j = 0; j < sh->len; j++) {
        newE += sh->nbuf[j];
        newK += sh->kbuf[j];
    }
    sh->dK = newK - oldK;
    sh->dE = newE - oldE;
}

void shot_accept(const Shot *sh, Traj *tr)
{
    const int L = tr->L;
    int j;
    if (sh->dir > 0) {
        for (j = 0; j < sh->len; j++) {
            int t = sh->tau + 1 + j;
            memcpy(tr->x + (size_t)t * (size_t)L, sh->buf + (size_t)j * (size_t)L, (size_t)L);
            tr->n[t] = sh->nbuf[j];
            tr->k[t] = sh->kbuf[j];
        }
    } else {
        for (j = 0; j < sh->len; j++) {
            int t = sh->tau - 1 - j; /* new slice t; transition (t, t+1) */
            memcpy(tr->x + (size_t)t * (size_t)L, sh->buf + (size_t)j * (size_t)L, (size_t)L);
            tr->n[t] = sh->nbuf[j];
            tr->k[t + 1] = sh->kbuf[j];
        }
    }
    tr->K += sh->dK;
    tr->E += sh->dE;
}
