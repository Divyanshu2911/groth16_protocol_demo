// src/keygen.c
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <pbc/pbc.h>
#include "circuit.h"
#include "poly.h"
#include "keys.h"
#include "fmt.h"

// Horner evaluation: out = coeffs[0] + coeffs[1]*t + ... (degree < m)
static void poly_eval(element_t out, element_t *coeffs, int m, element_t t, pairing_t pairing)
{
    (void)pairing; // unused here but kept for consistent signature
    element_set0(out);
    for (int i = m - 1; i >= 0; i--)
    {
        element_mul(out, out, t);         // out = out * t
        element_add(out, out, coeffs[i]); // out += coeffs[i]
    }
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        fprintf(stderr, "Usage: %s pairing.params d x y a0…ad\n", argv[0]);
        return 1;
    }

    fmt_init(1, stdout);
    fmt_banner("Key Generation (demo)");

    // --- init pairing from params file ---
    FILE *fp = fopen(argv[1], "r");
    if (!fp)
    {
        fprintf(stderr, "Error opening '%s': %s\n", argv[1], strerror(errno));
        return 1;
    }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *buf = (char *)malloc(sz + 1);
    size_t rd = fread(buf, 1, sz, fp);
    fclose(fp);
    if (rd != (size_t)sz)
    {
        fprintf(stderr, "Short read: %zu of %ld\n", rd, sz);
        free(buf);
        return 1;
    }
    buf[sz] = '\0';
    pbc_param_t params;
    pbc_param_init_set_buf(params, buf, sz + 1);
    free(buf);
    pairing_t pairing;
    pairing_init_pbc_param(pairing, params);

    // --- parse inputs ---
    int argi = 2;
    int d = atoi(argv[argi++]);

    element_t x, y;
    element_init_Zr(x, pairing);
    element_init_Zr(y, pairing);
    element_set_str(x, argv[argi++], 10);
    element_set_str(y, argv[argi++], 10);

    element_t *coeffs = (element_t *)malloc((d + 1) * sizeof(element_t));
    for (int i = 0; i <= d; i++)
    {
        element_init_Zr(coeffs[i], pairing);
        element_set_str(coeffs[i], argv[argi++], 10);
    }

    // --- R1CS and wires (for sizes & witness) ---
    r1cs_t r;
    element_t *wires;
    build_r1cs(d, coeffs, x, y, &r, &wires, pairing);
    int m = r.n_cons;
    int n = r.n_vars;
    fmt_kv_i("constraints (m)", m);
    fmt_kv_i("variables (n)", n);

    // --- interpolation points: τ_k = 1..m ---
    element_t *tau_pts = (element_t *)malloc(sizeof(element_t) * m);
    for (int i = 0; i < m; i++)
    {
        element_init_Zr(tau_pts[i], pairing);
        element_set_si(tau_pts[i], i + 1);
    }

    // --- sample toxic waste τ and choose bases g1,g2 ---
    element_t tau_secret;
    element_init_Zr(tau_secret, pairing);
    element_random(tau_secret);
    fmt_kv_e("tau", tau_secret);

    element_t g1, g2;
    element_init_G1(g1, pairing);
    element_random(g1);
    element_init_G2(g2, pairing);
    element_random(g2);
    fmt_kv_e("g1", g1);
    fmt_kv_e("g2", g2);

    // --- print g2^{tau^i} for i = 0..m-1 ---
    element_t tau_pow;
    element_init_Zr(tau_pow, pairing);
    element_set1(tau_pow);
    element_t g2_pow;
    element_init_G2(g2_pow, pairing);
    fmt_sub("g2^{tau^i}");
    for (int i = 0; i < m; i++)
    {
        if (i > 0)
            element_mul(tau_pow, tau_pow, tau_secret);
        element_pow_zn(g2_pow, g2, tau_pow);
        printf("  i=%d : ", i);
        element_printf("%B\n", g2_pow);
    }

    // --- temp arrays for per-column interpolation ---
    element_t *Ae = (element_t *)malloc(sizeof(element_t) * m);
    element_t *Be = (element_t *)malloc(sizeof(element_t) * m);
    element_t *Ce = (element_t *)malloc(sizeof(element_t) * m);

    element_t *pA = (element_t *)malloc(sizeof(element_t) * m);
    element_t *pB = (element_t *)malloc(sizeof(element_t) * m);
    element_t *pC = (element_t *)malloc(sizeof(element_t) * m);

    // --- accumulators for aggregates A(τ), B(τ), C(τ) ---
    element_t Aagg, Bagg, Cagg;
    element_init_Zr(Aagg, pairing);
    element_set0(Aagg);
    element_init_Zr(Bagg, pairing);
    element_set0(Bagg);
    element_init_Zr(Cagg, pairing);
    element_set0(Cagg);

    // --- per-variable query scalars and elements ---
    element_t valA, valB, valC;
    element_init_Zr(valA, pairing);
    element_init_Zr(valB, pairing);
    element_init_Zr(valC, pairing);

    element_t AqueryG1, BqueryG2, CqueryG1;
    element_init_G1(AqueryG1, pairing);
    element_init_G2(BqueryG2, pairing);
    element_init_G1(CqueryG1, pairing);

    fmt_sub("Per-variable queries");
    for (int j = 0; j < n; j++)
    {
        // column j evaluations at the m points
        for (int k = 0; k < m; k++)
        {
            element_init_Zr(Ae[k], pairing);
            element_set(Ae[k], r.A[k][j]);
            element_init_Zr(Be[k], pairing);
            element_set(Be[k], r.B[k][j]);
            element_init_Zr(Ce[k], pairing);
            element_set(Ce[k], r.C[k][j]);
        }
        // interpolate coefficient vectors (degree < m)
        lagrange_interpolation(pA, tau_pts, Ae, m, pairing);
        lagrange_interpolation(pB, tau_pts, Be, m, pairing);
        lagrange_interpolation(pC, tau_pts, Ce, m, pairing);

        // evaluate at secret tau
        poly_eval(valA, pA, m, tau_secret, pairing);
        poly_eval(valB, pB, m, tau_secret, pairing);
        poly_eval(valC, pC, m, tau_secret, pairing);

        // queries in the exponent
        element_pow_zn(AqueryG1, g1, valA);
        element_pow_zn(BqueryG2, g2, valB);
        element_pow_zn(CqueryG1, g1, valC);

        // print
        printf("  var %d\n", j);
        fmt_kv_e("    A_j(τ)", valA);
        fmt_kv_e("    B_j(τ)", valB);
        fmt_kv_e("    C_j(τ)", valC);
        fmt_kv_e("    A_query[G1]", AqueryG1);
        fmt_kv_e("    B_query[G2]", BqueryG2);
        fmt_kv_e("    C_query[G1]", CqueryG1);

        // aggregate: A(τ) += w_j * A_j(τ), etc.
        element_t t;
        element_init_Zr(t, pairing);
        element_mul(t, wires[j], valA);
        element_add(Aagg, Aagg, t);
        element_mul(t, wires[j], valB);
        element_add(Bagg, Bagg, t);
        element_mul(t, wires[j], valC);
        element_add(Cagg, Cagg, t);
        element_clear(t);

        // clear temporaries for this j
        for (int k = 0; k < m; k++)
        {
            element_clear(Ae[k]);
            element_clear(Be[k]);
            element_clear(Ce[k]);
            element_clear(pA[k]);
            element_clear(pB[k]);
            element_clear(pC[k]);
        }
    }

    fmt_hr();
    fmt_sub("Aggregate check (QAP divisibility)");

    // A(τ), B(τ), C(τ)
    fmt_kv_e("A(τ)", Aagg);
    fmt_kv_e("B(τ)", Bagg);
    fmt_kv_e("C(τ)", Cagg);

    // D = A(τ)·B(τ) − C(τ)
    element_t D;
    element_init_Zr(D, pairing);
    element_mul(D, Aagg, Bagg);
    element_t tmp;
    element_init_Zr(tmp, pairing);
    element_set(tmp, Cagg);
    element_sub(D, D, tmp);
    fmt_kv_e("A(τ)·B(τ) − C(τ)", D);

    // Z(τ) = ∏_{k=1..m} (τ − τ_k), with τ_k = {1..m}
    element_t Ztau;
    element_init_Zr(Ztau, pairing);
    element_set1(Ztau);
    for (int k = 0; k < m; k++)
    {
        element_t diff;
        element_init_Zr(diff, pairing);
        element_sub(diff, tau_secret, tau_pts[k]);
        element_mul(Ztau, Ztau, diff);
        element_clear(diff);
    }
    fmt_kv_e("Z(τ)", Ztau);

    // H(τ) = (A·B − C) / Z(τ)
    element_t invZ;
    element_init_Zr(invZ, pairing);
    element_invert(invZ, Ztau);
    element_t Htau;
    element_init_Zr(Htau, pairing);
    element_mul(Htau, D, invZ);
    fmt_kv_e("H(τ)", Htau);

    // sanity: check H(τ)*Z(τ) == D
    element_t check;
    element_init_Zr(check, pairing);
    element_mul(check, Htau, Ztau);
    int ok = (element_cmp(check, D) == 0);
    fmt_kv_s("division exact", ok ? "yes" : "NO (unexpected)");

    // --- cleanup (minimal; expand if you like) ---
    element_clear(check);
    element_clear(Htau);
    element_clear(invZ);
    element_clear(Ztau);
    element_clear(tmp);
    element_clear(D);
    element_clear(Aagg);
    element_clear(Bagg);
    element_clear(Cagg);
    element_clear(AqueryG1);
    element_clear(BqueryG2);
    element_clear(CqueryG1);
    element_clear(valA);
    element_clear(valB);
    element_clear(valC);
    for (int i = 0; i < m; i++)
        element_clear(tau_pts[i]);
    free(tau_pts);
    free(Ae);
    free(Be);
    free(Ce);
    free(pA);
    free(pB);
    free(pC);
    // (r1cs matrices and wires clearing omitted in this demo)
    return 0;
}
