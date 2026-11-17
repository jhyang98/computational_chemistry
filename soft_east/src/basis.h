/*
 * Basis of the configuration space.
 *
 * With use_translation = 1 the basis spans the zero-momentum (translation
 * invariant) sector: one state per orbit of the translation group,
 *   |r~> = p_r^{-1/2} sum_{j < p_r} T^j |r>,
 * where r is the smallest configuration of the orbit and p_r its size.
 * The Perron eigenvector of the tilted generator is translation invariant,
 * so the largest eigenvalue always lies in this sector.
 *
 * With use_translation = 0 every configuration is its own basis state.
 *
 * With exclude_empty = 1 the all-empty configuration is removed. This is used
 * for the strict East model (eps = 0), where the empty configuration is
 * disconnected from the rest and has eigenvalue 0 for all (s, g).
 */
#ifndef BASIS_H
#define BASIS_H

#include <stdint.h>

typedef struct {
    int L;
    int use_translation;
    int exclude_empty;
    long dim;          /* number of basis states */
    uint32_t *rep;     /* representative configuration of each state */
    int *period;       /* orbit size (1 without translation) */
    long *index;       /* configuration -> basis index, -1 if not in basis */
} Basis;

int basis_build(Basis *b, int L, int use_translation, int exclude_empty);
void basis_free(Basis *b);

#endif
