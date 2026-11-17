/*
 * Softened East model on a periodic ring of L sites.
 *
 * Site i is excited (n_i = 1) or not (n_i = 0); bit i of a configuration word
 * stores n_i. The facilitating neighbour of site i is its left neighbour i-1
 * (periodic). Rates (J = 1):
 *
 *   w_i(0 -> 1) = c       * C_i
 *   w_i(1 -> 0) = (1 - c) * C_i
 *   C_i         = n_{i-1} + eps * (1 - n_{i-1})
 *   c           = 1 / (1 + exp(1/T))
 *
 * Since C_i does not depend on n_i, detailed balance holds with respect to the
 * product measure with density c.
 */
#ifndef MODEL_H
#define MODEL_H

#include <stdint.h>

#define MODEL_LMAX 26

typedef struct {
    int L;       /* number of sites */
    double T;    /* temperature in units of J */
    double c;    /* equilibrium excitation density */
    double eps;  /* softness of the kinetic constraint */
} Model;

/* Returns 0 on success, -1 on invalid parameters. */
int model_init(Model *m, int L, double T, double eps);

int conf_popcount(uint32_t x);

/* Kinetic constraint C_i of configuration x. */
double model_constraint(const Model *m, uint32_t x, int i);

/* Total escape rate r(x) = sum_i C_i [c (1 - n_i) + (1 - c) n_i]. */
double model_escape_rate(const Model *m, uint32_t x);

/* Rotate a configuration by one site (periodic translation). */
uint32_t conf_rotate(uint32_t x, int L);

#endif
