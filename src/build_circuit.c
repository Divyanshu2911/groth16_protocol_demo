// ---------------------- src/build_circuit.c ----------------------
#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>
#include "../include/circuit.h"

/*
void build_r1cs(int d, element_t *coeffs, element_t x, element_t y,
                r1cs_t *r1cs, element_t **wires, pairing_t pairing) {
    int n_pow = d + 1;
    int n_sum = d + 1;
    int n_vars = n_pow + n_sum; // wire0(1) + x^i + s_i
    int n_cons = d  // constraints for x^i = x * x^{i-1}
                + d  // constraints for s_i = s_{i-1} + a_i*x^i
                + 1; // final check s_d - y = 0

    // Allocate wires
    *wires = malloc(sizeof(element_t) * n_vars);
    for(int i=0;i<n_vars;i++) element_init_Zr((*wires)[i], pairing);

    // Compute wires: w0=1, w1=x, w2=x*x,...
    element_set1((*wires)[0]);
    element_set((*wires)[1], x);
    for(int i=2;i<n_pow;i++){
        element_mul((*wires)[i], (*wires)[i-1], x);
    }
    // Partial sums: s0 = a0*w0, s_i = s_{i-1} + a_i*w_i
    int off = n_pow;
    element_mul((*wires)[off+0], coeffs[0], (*wires)[0]);
    for(int i=1;i<n_sum;i++){
        element_mul((*wires)[off+i], coeffs[i], (*wires)[i]);
        element_add((*wires)[off+i], (*wires)[off+i], (*wires)[off+i-1]);
    }

    // Allocate R1CS matrices (dense for demo)
    r1cs->n_vars = n_vars;
    r1cs->n_cons = n_cons;
    r1cs->A = malloc(sizeof(element_t*) * n_cons);
    r1cs->B = malloc(sizeof(element_t*) * n_cons);
    r1cs->C = malloc(sizeof(element_t*) * n_cons);
    for(int i=0;i<n_cons;i++){
        r1cs->A[i] = malloc(sizeof(element_t) * n_vars);
        r1cs->B[i] = malloc(sizeof(element_t) * n_vars);
        r1cs->C[i] = malloc(sizeof(element_t) * n_vars);
        for(int j=0;j<n_vars;j++){
            element_init_Zr(r1cs->A[i][j], pairing);
            element_init_Zr(r1cs->B[i][j], pairing);
            element_init_Zr(r1cs->C[i][j], pairing);
            element_set0(r1cs->A[i][j]);
            element_set0(r1cs->B[i][j]);
            element_set0(r1cs->C[i][j]);
        }
    }
    int ci = 0;
    // 1) x^i constraints: for i=1..d: w_i * x = w_{i+1}
    for(int i=1;i<d+1;i++){
        element_set1(r1cs->A[ci][i]); // w_i
        element_set1(r1cs->B[ci][1]); // x = wire index 1
        element_set1(r1cs->C[ci][i+1]);
        ci++;
    }
    // 2) sum constraints: for i=1..d: s_{i-1} + a_i*w_i = s_i
    for(int i=1;i<n_sum;i++){
        element_set1(r1cs->A[ci][off+i-1]);
        element_set1(r1cs->A[ci][n_pow + i]); // a_i*w_i is precomputed as wire
        element_set1(r1cs->B[ci][0]);       // multiply by 1
        element_set1(r1cs->C[ci][off+i]);   // s_i
        ci++;
    }
    // 3) final: s_d - y = 0  => s_d * 1 = y
    element_set1(r1cs->A[ci][off + d]); // s_d
    element_set1(r1cs->B[ci][0]);       // 1
    // move y into C[ci][0]
    element_set(r1cs->C[ci][0], y);
}
*/

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s a.param d x y a0 a1 … ad\n"
            "  a.param : PBC parameter file (any path)\n"
            "  d              : degree of the polynomial\n"
            "  x              : evaluation point (in Zr)\n"
            "  y              : claimed result y = ∑ a_i·x^i\n"
            "  a0…ad          : coefficients in Zr\n",
            argv[0]
        );
        return 1;
    }

    // 1) Load the params file into a buffer
    const char *params_path = argv[1];
    FILE *fp = fopen(params_path, "r");
    if (!fp) { perror("fopen a.param"); return 1; }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *param_buf = malloc(sz + 1);
    if (fread(param_buf, 1, sz, fp) != sz) {
        perror("fread a.param"); return 1;
    }
    param_buf[sz] = '\0';
    fclose(fp);

    // 2) Init PBC param and pairing
    pbc_param_t params;
    pbc_param_init_set_buf(params, param_buf, sz + 1);
    free(param_buf);

    pairing_t pairing;
    pairing_init_pbc_param(pairing, params);

    // 3) Parse degree, x, y
    int argi = 2;
    int d = atoi(argv[argi++]);
    element_t x, y;
    element_init_Zr(x, pairing);
    element_init_Zr(y, pairing);
    element_set_str(x, argv[argi++], 10);
    element_set_str(y, argv[argi++], 10);

    // 4) Parse coefficients a0…ad
    element_t *coeffs = malloc((d + 1) * sizeof(element_t));
    for (int i = 0; i <= d; i++) {
        element_init_Zr(coeffs[i], pairing);
        element_set_str(coeffs[i], argv[argi++], 10);
    }

    // 5) Build R1CS and wires
    r1cs_t r1cs;
    element_t *wires;
    build_r1cs(d, coeffs, x, y, &r1cs, &wires, pairing);

    printf("Built R1CS: %d variables, %d constraints\n",
           r1cs.n_vars, r1cs.n_cons);

    // 6) Print all wire values
    for (int i = 0; i < r1cs.n_vars; i++) {
        element_printf("  wire[%d] = %B\n", i, wires[i]);
    }

    // Cleanup omitted for brevity
    return 0;
}