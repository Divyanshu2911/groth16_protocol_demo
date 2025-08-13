// src/interpolate.c
#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>
#include "../include/circuit.h"
#include "../include/poly.h"

// Helper: multiply two polynomials a,b of degree< m,
// truncating to degree < m: result in res (all length m).
static void poly_mul(
    element_t *res,
    element_t *a,
    element_t *b,
    int m,
    pairing_t pairing
) {
    // temp polynomial
    element_t *tmp = malloc(sizeof(element_t)*m);
    for(int i=0;i<m;i++){
        element_init_Zr(tmp[i], pairing);
        element_set0(tmp[i]);
    }
    for(int i=0;i<m;i++){
        if (!element_is0(a[i])) {
            for(int j=0; j+i < m; j++){
                if (!element_is0(b[j])) {
                    element_t t; element_init_Zr(t, pairing);
                    element_mul(t, a[i], b[j]);
                    element_add(tmp[i+j], tmp[i+j], t);
                    element_clear(t);
                }
            }
        }
    }
    // copy back
    for(int i=0;i<m;i++){
        element_set(res[i], tmp[i]);
        element_clear(tmp[i]);
    }
    free(tmp);
}

void lagrange_interpolation(
    element_t *out,
    element_t *tau,
    element_t *eval,
    int m,
    pairing_t pairing
) {
    // Initialize out[] = 0
    for(int i=0;i<m;i++){
        element_init_Zr(out[i], pairing);
        element_set0(out[i]);
    }

    // Workspace for basis polynomial
    element_t *basis = malloc(sizeof(element_t)*m);
    element_t *next_basis = malloc(sizeof(element_t)*m);

    for(int k=0;k<m;k++){
        // 1) build ℓ_k numerator = ∏_{j≠k} (x - τ[j])
        // start with basis = [1,0,0,...]
        for(int i=0;i<m;i++){
            element_init_Zr(basis[i], pairing);
            element_set0(basis[i]);
        }
        element_set1(basis[0]);

        for(int j=0;j<m;j++){
            if (j==k) continue;
            // multiply basis by (x - τ[j])
            // i.e. polynomial [ -τ[j], 1 ]
            for(int i=0;i<m;i++){
                element_init_Zr(next_basis[i], pairing);
                element_set0(next_basis[i]);
            }
            // next_basis[i+1] += basis[i]
            for(int i=0;i<m-1;i++){
                if (!element_is0(basis[i])) {
                    element_add(next_basis[i+1], next_basis[i+1], basis[i]);
                }
            }
            // next_basis[i] += basis[i] * (-τ[j])
            for(int i=0;i<m;i++){
                if (!element_is0(basis[i])) {
                    element_t neg; element_init_Zr(neg, pairing);
                    element_neg(neg, tau[j]);
                    element_mul(neg, neg, basis[i]);
                    element_add(next_basis[i], next_basis[i], neg);
                    element_clear(neg);
                }
            }
            // swap next_basis→basis
            for(int i=0;i<m;i++){
                element_set(basis[i], next_basis[i]);
                element_clear(next_basis[i]);
            }
        }

        // 2) compute denominator denom = ∏_{j≠k} (τ[k] - τ[j])
        element_t denom; element_init_Zr(denom, pairing);
        element_set1(denom);
        for(int j=0;j<m;j++){
            if (j==k) continue;
            element_t diff; element_init_Zr(diff, pairing);
            element_sub(diff, tau[k], tau[j]);
            element_mul(denom, denom, diff);
            element_clear(diff);
        }
        // invert denom
        element_t inv_denom; element_init_Zr(inv_denom, pairing);
        element_invert(inv_denom, denom);
        element_clear(denom);

        // 3) scale basis by eval[k] * inv_denom
        element_t scale; element_init_Zr(scale, pairing);
        element_mul(scale, eval[k], inv_denom);
        element_clear(inv_denom);

        for(int i=0;i<m;i++){
            element_t tmp; element_init_Zr(tmp, pairing);
            element_mul(tmp, basis[i], scale);
            element_add(out[i], out[i], tmp);
            element_clear(tmp);
        }
        element_clear(scale);

        // clear basis
        for(int i=0;i<m;i++){
            element_clear(basis[i]);
        }
    }

    free(basis);
    free(next_basis);
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s path/to/a.param [degree of x in y = f(x)] x y [comma separated coefficients starting from x^0]\n"
            "  (same args as build_circuit)\n", argv[0]);
        return 1;
    }

    // --- init pairing (as in build_circuit) ---
    FILE *fp = fopen(argv[1],"r");
    fseek(fp,0,SEEK_END); long sz=ftell(fp); fseek(fp,0,SEEK_SET);
    char *buf = malloc(sz+1); 
    if (fread(buf, 1, sz, fp) != sz) {
        perror("fread a.param");
        return 1;
    }
    buf[sz] = '\0';
    fclose(fp);
    pbc_param_t params; pbc_param_init_set_buf(params,buf,sz+1); free(buf);
    pairing_t pairing; pairing_init_pbc_param(pairing, params);

    // --- parse inputs & build R1CS ---
    int i, argi = 2;
    int d = atoi(argv[argi++]);
    element_t x, y;
    element_init_Zr(x,pairing); element_init_Zr(y,pairing);
    element_set_str(x,argv[argi++],10);
    element_set_str(y,argv[argi++],10);
    element_t *coeffs = malloc((d+1)*sizeof(element_t));
    for(i=0;i<=d;i++){
        element_init_Zr(coeffs[i],pairing);
        element_set_str(coeffs[i],argv[argi++],10);
    }
    r1cs_t r1cs; element_t *wires;
    build_r1cs(d, coeffs, x, y, &r1cs, &wires, pairing);

    int m = r1cs.n_cons;
    int n = r1cs.n_vars;

    // --- pick tau points = [1,2,…,m] in Zr ---
    element_t *tau = malloc(sizeof(element_t)*m);
    for(i=0;i<m;i++){
        element_init_Zr(tau[i],pairing);
        element_set_si(tau[i], i+1);
    }

    // --- for each variable j, gather A_eval[k] = r1cs.A[k][j] ---
    element_t *A_eval = malloc(sizeof(element_t)*m);
    element_t *B_eval = malloc(sizeof(element_t)*m);
    element_t *C_eval = malloc(sizeof(element_t)*m);

    element_t *polyA = malloc(sizeof(element_t)*m);
    element_t *polyB = malloc(sizeof(element_t)*m);
    element_t *polyC = malloc(sizeof(element_t)*m);

    for(int j=0;j<n;j++){
        // extract column j
        for(int k=0;k<m;k++){
            element_init_Zr(A_eval[k],pairing);
            element_set(A_eval[k], r1cs.A[k][j]);
            element_init_Zr(B_eval[k],pairing);
            element_set(B_eval[k], r1cs.B[k][j]);
            element_init_Zr(C_eval[k],pairing);
            element_set(C_eval[k], r1cs.C[k][j]);
        }
        // interp: get polynomials of degree < m
        lagrange_interpolation(polyA, tau, A_eval, m, pairing);
        lagrange_interpolation(polyB, tau, B_eval, m, pairing);
        lagrange_interpolation(polyC, tau, C_eval, m, pairing);

        // Print out the coefficient vectors:
        printf("Variable %d → A_j(x) coeffs:", j);
        for(i=0;i<m;i++) element_printf(" %B", polyA[i]);
        printf("\n");
        printf("           B_j(x) coeffs:");
        for(i=0;i<m;i++) element_printf(" %B", polyB[i]);
        printf("\n");
        printf("           C_j(x) coeffs:");
        for(i=0;i<m;i++) element_printf(" %B", polyC[i]);
        printf("\n\n");

        // clear eval & poly arrays for next j
        for(i=0;i<m;i++){
            element_clear(A_eval[i]);
            element_clear(B_eval[i]);
            element_clear(C_eval[i]);
            element_clear(polyA[i]);
            element_clear(polyB[i]);
            element_clear(polyC[i]);
        }
    }

    // Cleanup r1cs, wires, tau, etc. (omitted for brevity)
    return 0;
}
