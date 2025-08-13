// src/circuit.c
#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>
#include "../include/circuit.h"

void build_r1cs(int d,
                element_t *coeffs,
                element_t x,
                element_t y,
                r1cs_t *r1cs,
                element_t **wires,
                pairing_t pairing) {
    int n_pow  = d + 1;        // x^0 … x^d
    int n_sum  = d + 1;        // partial sums s₀ … s_d
    int n_vars = n_pow + n_sum;  
    int n_cons = d        // xᶦ = x·xᶦ⁻¹
                 + d      // sᵢ = sᵢ₋₁ + aᵢ·xᶦ
                 + 1;     // final check

    // 1) allocate and init wires
    *wires = malloc(sizeof(element_t) * n_vars);
    for(int i = 0; i < n_vars; i++) {
        element_init_Zr((*wires)[i], pairing);
    }

    // 2) compute powers: w₀=1, w₁=x, w₂=x², …
    element_set1((*wires)[0]);
    element_set((*wires)[1], x);
    for(int i = 2; i < n_pow; i++) {
        element_mul((*wires)[i], (*wires)[i-1], x);
    }

    // 3) compute partial sums: s₀ = a₀·w₀, sᵢ = sᵢ₋₁ + aᵢ·wᵢ
    int off = n_pow;
    element_mul((*wires)[off + 0], coeffs[0], (*wires)[0]);
    for(int i = 1; i < n_sum; i++) {
        element_mul((*wires)[off + i], coeffs[i], (*wires)[i]);
        element_add((*wires)[off + i],
                    (*wires)[off + i],
                    (*wires)[off + i - 1]);
    }

    // 4) allocate & zero‐init R1CS matrices
    r1cs->n_vars = n_vars;
    r1cs->n_cons = n_cons;
    r1cs->A = malloc(sizeof(element_t*) * n_cons);
    r1cs->B = malloc(sizeof(element_t*) * n_cons);
    r1cs->C = malloc(sizeof(element_t*) * n_cons);
    for(int c = 0; c < n_cons; c++) {
        r1cs->A[c] = malloc(sizeof(element_t) * n_vars);
        r1cs->B[c] = malloc(sizeof(element_t) * n_vars);
        r1cs->C[c] = malloc(sizeof(element_t) * n_vars);
        for(int v = 0; v < n_vars; v++) {
            element_init_Zr(r1cs->A[c][v], pairing);
            element_init_Zr(r1cs->B[c][v], pairing);
            element_init_Zr(r1cs->C[c][v], pairing);
            element_set0(r1cs->A[c][v]);
            element_set0(r1cs->B[c][v]);
            element_set0(r1cs->C[c][v]);
        }
    }

    // 5) fill constraints
    int ci = 0;
    //   a) x^i constraints: w_i * x = w_{i+1}
    for(int i = 1; i < n_pow; i++, ci++) {
        element_set1(r1cs->A[ci][i]);
        element_set1(r1cs->B[ci][1]);
        element_set1(r1cs->C[ci][i+1]);
    }
    //   b) sum constraints: s_{i-1} + (a_i·x^i) = s_i
    for(int i = 1; i < n_sum; i++, ci++) {
        element_set1(r1cs->A[ci][off + i - 1]); // s_{i-1}
        element_set1(r1cs->A[ci][off + i]);     // a_i·x^i
        element_set1(r1cs->B[ci][0]);           // ×1
        element_set1(r1cs->C[ci][off + i]);     // s_i
    }
    //   c) final check: s_d * 1 = y
    element_set1(r1cs->A[ci][off + d]);
    element_set1(r1cs->B[ci][0]);
    element_set(r1cs->C[ci][0], y);
}
