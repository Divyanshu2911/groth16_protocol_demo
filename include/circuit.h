// ---------------------- include/circuit.h ----------------------
#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <pbc/pbc.h>

// Structure to hold R1CS in sparse form
typedef struct {
    int n_vars, n_cons;
    // Each constraint: A, B, C are arrays of length n_vars
    // stored as lists of (index,value)
    // For simplicity: dense matrix
    element_t **A, **B, **C; // dimensions [n_cons][n_vars]
} r1cs_t;

// Build R1CS for polynomial evaluation: y = \sum_{i=0}^d a_i x^i
// Inputs: degree d; coeffs a[0..d]; evaluation point x; claimed y
// Outputs: r1cs struct, wire vector w (length n_vars)
void build_r1cs(int d, element_t *coeffs, element_t x, element_t y,
                r1cs_t *r1cs, element_t **wires, pairing_t pairing);

#endif // CIRCUIT_H