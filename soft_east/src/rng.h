/*
 * xoshiro256** pseudo-random generator (Blackman & Vigna), seeded with
 * splitmix64. Portable replacement for rand(), whose RAND_MAX is only 32767
 * on some platforms.
 */
#ifndef RNG_H
#define RNG_H

#include <stdint.h>

typedef struct {
    uint64_t s[4];
} Rng;

void rng_seed(Rng *r, uint64_t seed);

/* Seed from the clock when seed == 0; returns the seed actually used. */
uint64_t rng_seed_auto(Rng *r, uint64_t seed);

uint64_t rng_next(Rng *r);

/* Uniform double in [0, 1) with 53 random bits. */
double rng_uniform(Rng *r);

/* Uniform integer in [0, n), n > 0. */
int rng_int(Rng *r, int n);

#endif
