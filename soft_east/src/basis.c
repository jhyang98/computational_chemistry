#include "basis.h"
#include "model.h"

#include <stdlib.h>

int basis_build(Basis *b, int L, int use_translation, int exclude_empty)
{
    uint32_t nconf, x;
    long count = 0;

    b->L = L;
    b->use_translation = use_translation;
    b->exclude_empty = exclude_empty;
    b->rep = NULL;
    b->period = NULL;
    b->index = NULL;
    b->dim = 0;

    if (L < 2 || L > MODEL_LMAX) return -1;
    nconf = 1u << L;

    b->index = (long *)malloc(sizeof(long) * nconf);
    b->rep = (uint32_t *)malloc(sizeof(uint32_t) * nconf);
    b->period = (int *)malloc(sizeof(int) * nconf);
    if (!b->index || !b->rep || !b->period) {
        basis_free(b);
        return -1;
    }
    for (x = 0; x < nconf; x++) b->index[x] = -1;

    for (x = 0; x < nconf; x++) {
        if (b->index[x] != -1) continue; /* already assigned to an orbit */
        if (exclude_empty && x == 0u) continue;
        if (use_translation) {
            /* x is the smallest member of its orbit: all smaller ones were
               visited before and their orbits fully assigned. */
            uint32_t y = x;
            int p = 0;
            do {
                b->index[y] = count;
                y = conf_rotate(y, L);
                p++;
            } while (y != x);
            b->rep[count] = x;
            b->period[count] = p;
        } else {
            b->index[x] = count;
            b->rep[count] = x;
            b->period[count] = 1;
        }
        count++;
    }
    b->dim = count;

    /* shrink the per-state arrays */
    {
        uint32_t *r2 = (uint32_t *)realloc(b->rep, sizeof(uint32_t) * (size_t)count);
        int *p2 = (int *)realloc(b->period, sizeof(int) * (size_t)count);
        if (r2) b->rep = r2;
        if (p2) b->period = p2;
    }
    return 0;
}

void basis_free(Basis *b)
{
    free(b->rep);
    free(b->period);
    free(b->index);
    b->rep = NULL;
    b->period = NULL;
    b->index = NULL;
    b->dim = 0;
}
