// ---------------------- include/poly.h ----------------------
#ifndef POLY_H
#define POLY_H

#include <pbc/pbc.h>

// Given m points (tau_i, eval_i), compute coeffs of polynomial deg<m
// via naive Lagrange interpolation
// allocates out[0..m-1]

/**
 * Given:
 *   - tau[0..m-1]: the distinct points τₖ
 *   - eval[0..m-1]: the values A_j(τₖ) (or B_j, or C_j)
 *   - m            : number of constraints
 *
 * Compute the unique polynomial
 *   P(x) = ∑_{k=0..m-1} eval[k] · ℓ_k(x)
 * of degree < m, where ℓ_k are the Lagrange basis polynomials.
 *
 * The output is placed in out[0..m-1], i.e. the coefficients
 * P(x) = out[0] + out[1]·x + … + out[m-1]·x^{m-1}.
 *
 * Preconditions:
 *   - out, tau, eval have been allocated to length m
 *   - pairing is already initialized
 */

void lagrange_interpolation(element_t *out, element_t *tau, element_t *eval,
                            int m, pairing_t pairing);

#endif // POLY_H
