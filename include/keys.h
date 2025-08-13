// ---------------------- include/keys.h ----------------------
#ifndef KEYS_H
#define KEYS_H

#include <pbc/pbc.h>
#include "circuit.h"

// Proving and verifying key structures
typedef struct {
    element_t *A_query; // G1 elements length = n_vars
    element_t *B_query; // G2 elements length = n_vars
    element_t *C_query; // G1 elements length = n_vars
} pk_t;

typedef struct {
    element_t *g2_tau;    // G2^{τ^i}
    element_t g2_gamma;   // G2^{γ}
    element_t g2;         // base
} vk_t;

// Generate keys from R1CS and tau-powers
void keygen(const r1cs_t *r1cs, int deg, pairing_t pairing,
            pk_t *pk, vk_t *vk);

#endif // KEYS_H