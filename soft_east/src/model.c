#include "model.h"

#include <math.h>

int model_init(Model *m, int L, double T, double eps)
{
    if (L < 2 || L > MODEL_LMAX) return -1;
    if (!(T > 0.0)) return -1;
    if (eps < 0.0 || eps > 1.0) return -1;
    m->L = L;
    m->T = T;
    m->c = 1.0 / (1.0 + exp(1.0 / T));
    m->eps = eps;
    return 0;
}

int conf_popcount(uint32_t x)
{
    int n = 0;
    while (x) {
        x &= x - 1u;
        n++;
    }
    return n;
}

double model_constraint(const Model *m, uint32_t x, int i)
{
    int left = (i == 0) ? m->L - 1 : i - 1;
    int nl = (int)((x >> left) & 1u);
    return nl ? 1.0 : m->eps;
}

double model_escape_rate(const Model *m, uint32_t x)
{
    double r = 0.0;
    int i;
    for (i = 0; i < m->L; i++) {
        int ni = (int)((x >> i) & 1u);
        double ci = model_constraint(m, x, i);
        r += ci * (ni ? (1.0 - m->c) : m->c);
    }
    return r;
}

uint32_t conf_rotate(uint32_t x, int L)
{
    uint32_t mask = (L >= 32) ? 0xffffffffu : ((1u << L) - 1u);
    return ((x << 1) | (x >> (L - 1))) & mask;
}
