#include "rng.h"

#include <time.h>

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t rotl(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

void rng_seed(Rng *r, uint64_t seed)
{
    uint64_t x = seed;
    int i;
    for (i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
}

uint64_t rng_seed_auto(Rng *r, uint64_t seed)
{
    if (seed == 0) {
        uint64_t t = (uint64_t)time(NULL);
        uint64_t c = (uint64_t)clock();
        uint64_t x = t * 0x9E3779B97F4A7C15ULL ^ (c << 17) ^ (uint64_t)(uintptr_t)&seed;
        seed = splitmix64(&x);
        if (seed == 0) seed = 1;
    }
    rng_seed(r, seed);
    return seed;
}

uint64_t rng_next(Rng *r)
{
    uint64_t *s = r->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

double rng_uniform(Rng *r)
{
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

int rng_int(Rng *r, int n)
{
    /* multiply-shift: bias below 2^-32 relative, negligible for lattice sizes */
    uint64_t x = rng_next(r) >> 32;
    return (int)((x * (uint64_t)n) >> 32);
}
