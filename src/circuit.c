// src/circuit.c
#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>
#include "../include/circuit.h"

// src/circuit.c  (replace the whole build_r1cs with this)
void build_r1cs(int d,
                element_t *coeffs,
                element_t x,
                element_t y,
                r1cs_t *r1cs,
                element_t **wires,
                pairing_t pairing)
{
    // Layout:
    // wires[0..d]     : w_i = x^i   (w0=1, w1=x, …, wd=x^d)
    // wires[d+1..2d+1]: s_i partial sums (s0..sd)
    int n_pow = d + 1;
    int n_sum = d + 1;
    int n_vars = n_pow + n_sum;

    // Constraints:
    // 1) s0 = a0 * w0                      -> 1
    // 2) w_i * w1 = w_{i+1} for i=1..d-1   -> (d-1)
    // 3) s_i = s_{i-1} + a_i * w_i for i=1..d -> d
    // 4) s_d * 1 = y                       -> 1
    int n_cons = (d > 0 ? (2 * d + 1) : 2);

    // ---- wires (w^i and partial sums) ----
    *wires = malloc(sizeof(element_t) * n_vars);
    for (int i = 0; i < n_vars; i++)
        element_init_Zr((*wires)[i], pairing);

    // powers
    element_set1((*wires)[0]); // w0 = 1
    if (n_pow >= 2)
        element_set((*wires)[1], x); // w1 = x
    for (int i = 2; i < n_pow; i++)
    {
        element_mul((*wires)[i], (*wires)[i - 1], x); // wi = w_{i-1} * x
    }

    // partial sums (computed just for printing/debug)
    int off = n_pow; // start of s_i region
    // s0 = a0 * w0
    element_mul((*wires)[off + 0], coeffs[0], (*wires)[0]);
    for (int i = 1; i < n_sum; i++)
    { // i = 1..d
        element_t term;
        element_init_Zr(term, pairing);
        element_mul(term, coeffs[i], (*wires)[i]);                   // a_i * w_i
        element_add((*wires)[off + i], (*wires)[off + i - 1], term); // s_i = s_{i-1} + term
        element_clear(term);
    }

    // ---- allocate R1CS A,B,C ----
    r1cs->n_vars = n_vars;
    r1cs->n_cons = n_cons;
    r1cs->A = malloc(sizeof(element_t *) * n_cons);
    r1cs->B = malloc(sizeof(element_t *) * n_cons);
    r1cs->C = malloc(sizeof(element_t *) * n_cons);
    for (int c = 0; c < n_cons; c++)
    {
        r1cs->A[c] = malloc(sizeof(element_t) * n_vars);
        r1cs->B[c] = malloc(sizeof(element_t) * n_vars);
        r1cs->C[c] = malloc(sizeof(element_t) * n_vars);
        for (int v = 0; v < n_vars; v++)
        {
            element_init_Zr(r1cs->A[c][v], pairing);
            element_set0(r1cs->A[c][v]);
            element_init_Zr(r1cs->B[c][v], pairing);
            element_set0(r1cs->B[c][v]);
            element_init_Zr(r1cs->C[c][v], pairing);
            element_set0(r1cs->C[c][v]);
        }
    }

    // ---- fill constraints ----
    int ci = 0;

    // (1) s0 = a0 * w0  =>  (a0·w0) * 1 = s0
    element_set(r1cs->A[ci][0], coeffs[0]); // a0 · w0
    element_set1(r1cs->B[ci][0]);           // × 1 (wire w0 == 1)
    element_set1(r1cs->C[ci][off + 0]);     // = s0
    ci++;

    // (2) power chain: w_i * w1 = w_{i+1} for i=1..d-1
    for (int i = 1; i <= d - 1; i++, ci++)
    {
        element_set1(r1cs->A[ci][i]);     // w_i
        element_set1(r1cs->B[ci][1]);     // × w1 (which is x)
        element_set1(r1cs->C[ci][i + 1]); // = w_{i+1}
    }

    // (3) sums: s_i = s_{i-1} + a_i·w_i   encoded as  (s_{i-1} + a_i·w_i) * 1 = s_i
    for (int i = 1; i <= d; i++, ci++)
    {
        element_set1(r1cs->A[ci][off + i - 1]); // s_{i-1}
        element_set(r1cs->A[ci][i], coeffs[i]); // + a_i · w_i
        element_set1(r1cs->B[ci][0]);           // × 1
        element_set1(r1cs->C[ci][off + i]);     // = s_i
    }

    // (4) final check: s_d * 1 = y (y lives as constant via w0 == 1)
    element_set1(r1cs->A[ci][off + d]); // s_d
    element_set1(r1cs->B[ci][0]);       // × 1
    element_set(r1cs->C[ci][0], y);     // = y (since w0 == 1)
}
