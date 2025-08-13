// src/verifier.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <arpa/inet.h>   // ntohl
#include <pbc/pbc.h>
#include "fmt.h"

// robust readers -------------------------------------------------
static int read_exact(FILE *f, void *p, size_t n) {
    unsigned char *b = (unsigned char*)p;
    size_t rtot = 0;
    while (rtot < n) {
        size_t r = fread(b + rtot, 1, n - rtot, f);
        if (r == 0) return -1;
        rtot += r;
    }
    return 0;
}
static int read_u32(FILE *f, uint32_t *out_host) {
    uint32_t v_be;
    if (read_exact(f, &v_be, 4) != 0) return -1;
    *out_host = ntohl(v_be);
    return 0;
}
static int read_elem_G1(FILE *f, pairing_t pairing, element_t out) {
    uint32_t len; if (read_u32(f, &len) != 0) return -1;
    unsigned char *buf = (unsigned char*)malloc(len);
    if (!buf) return -1;
    if (read_exact(f, buf, len) != 0) { free(buf); return -1; }
    element_init_G1(out, pairing);
    element_from_bytes(out, buf);
    free(buf);
    return 0;
}
static int read_elem_G2(FILE *f, pairing_t pairing, element_t out) {
    uint32_t len; if (read_u32(f, &len) != 0) return -1;
    unsigned char *buf = (unsigned char*)malloc(len);
    if (!buf) return -1;
    if (read_exact(f, buf, len) != 0) { free(buf); return -1; }
    element_init_G2(out, pairing);
    element_from_bytes(out, buf);
    free(buf);
    return 0;
}
// ----------------------------------------------------------------

// Compute Z(x) = ∏_{k=1}^m (x - k) coefficients in F_r (degree m)
static void compute_Z_coeffs(element_t *coef, int m, pairing_t pairing) {
    // coef has length m+1, initialize to 0; start with P(x) = 1
    for (int i = 0; i <= m; i++) { element_init_Zr(coef[i], pairing); element_set0(coef[i]); }
    element_set1(coef[0]); // degree 0: 1

    // temp buffer
    element_t *next = (element_t*)malloc(sizeof(element_t)*(m+1));
    for (int step = 1; step <= m; step++) {
        for (int i = 0; i <= m; i++) { element_init_Zr(next[i], pairing); element_set0(next[i]); }
        // multiply current poly by (x - step)
        // next[i+1] += coef[i]      (× x)
        for (int i = 0; i <= m-1; i++) element_add(next[i+1], next[i+1], coef[i]);
        // next[i] += coef[i] * (-step)
        element_t negk; element_init_Zr(negk, pairing); element_set_si(negk, -step);
        for (int i = 0; i <= m; i++) {
            if (!element_is0(coef[i])) {
                element_t t; element_init_Zr(t, pairing);
                element_mul(t, coef[i], negk);
                element_add(next[i], next[i], t);
                element_clear(t);
            }
        }
        element_clear(negk);
        // swap into coef
        for (int i = 0; i <= m; i++) { element_set(coef[i], next[i]); element_clear(next[i]); }
    }
    free(next);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s pairing.params [proof.bin]\n", argv[0]);
        return 1;
    }
    const char *proof_path = (argc > 2 ? argv[2] : "proof_demo.bin");

    fmt_init(1, stdout);
    fmt_banner("Verifier (QAP-aware demo)");

    // --- load pairing params ---
    FILE *fp = fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "Error opening '%s': %s\n", argv[1], strerror(errno)); return 1; }
    fseek(fp, 0, SEEK_END); long sz = ftell(fp); fseek(fp, 0, SEEK_SET);
    char *pbuf = (char*)malloc(sz + 1);
    size_t rd = fread(pbuf, 1, sz, fp); fclose(fp);
    if (rd != (size_t)sz) { fprintf(stderr, "Short read: %zu of %ld\n", rd, sz); free(pbuf); return 1; }
    pbuf[sz] = '\0';
    pbc_param_t params; pbc_param_init_set_buf(params, pbuf, sz + 1); free(pbuf);
    pairing_t pairing; pairing_init_pbc_param(pairing, params);

    // --- read proof file contents ---
    FILE *pf = fopen(proof_path, "rb");
    if (!pf) { fprintf(stderr, "Error opening '%s': %s\n", proof_path, strerror(errno)); return 1; }

    element_t piA, piB, piC, piH, g2;
    if (read_elem_G1(pf, pairing, piA) != 0) { fprintf(stderr, "read piA failed\n"); fclose(pf); return 1; }
    if (read_elem_G2(pf, pairing, piB) != 0) { fprintf(stderr, "read piB failed\n"); fclose(pf); return 1; }
    if (read_elem_G1(pf, pairing, piC) != 0) { fprintf(stderr, "read piC failed\n"); fclose(pf); return 1; }
    if (read_elem_G1(pf, pairing, piH) != 0) { fprintf(stderr, "read piH failed\n"); fclose(pf); return 1; }
    if (read_elem_G2(pf, pairing, g2)  != 0) { fprintf(stderr, "read g2 failed\n");  fclose(pf); return 1; }

    uint32_t m;
    if (read_u32(pf, &m) != 0) { fprintf(stderr, "read m failed\n"); fclose(pf); return 1; }

    element_t *g2_tau = (element_t*)malloc(sizeof(element_t)*(m+1));
    for (uint32_t i = 0; i <= m; i++) {
        if (read_elem_G2(pf, pairing, g2_tau[i]) != 0) { fprintf(stderr, "read g2^tau^i failed\n"); fclose(pf); return 1; }
    }
    fclose(pf);

    fmt_kv_s("proof file", proof_path);
    fmt_sub("Proof elements (decoded)");
    fmt_kv_e("piA (G1)", piA);
    fmt_kv_e("piB (G2)", piB);
    fmt_kv_e("piC (G1)", piC);
    fmt_kv_e("piH (G1)", piH);
    fmt_kv_e("g2  (G2)", g2);
    fmt_kv_i("m", m);

    // --- build g2^{Z(τ)} from g2^{τ^i} and Z(x) coefficients ---
    element_t *coef = (element_t*)malloc(sizeof(element_t)*(m+1));
    compute_Z_coeffs(coef, (int)m, pairing);

    element_t g2Z; element_init_G2(g2Z, pairing); element_set0(g2Z); // identity
    element_t tmp; element_init_G2(tmp, pairing);
    for (uint32_t i = 0; i <= m; i++) {
        // tmp = (g2^{τ^i})^{coef[i]}
        element_pow_zn(tmp, g2_tau[i], coef[i]);
        element_add(g2Z, g2Z, tmp); // additive notation → sums exponents
    }
    fmt_sub("Reconstructed");
    fmt_kv_e("g2^{Z(τ)}", g2Z);

    // --- pairing check: e(piA,piB) ?= e(piC,g2) * e(piH, g2^{Z(τ)}) ---
    fmt_sub("Pairing check");
    element_t L, R, T;
    element_init_GT(L, pairing); element_init_GT(R, pairing); element_init_GT(T, pairing);
    pairing_apply(L, piA, piB, pairing);
    pairing_apply(R, piC, g2, pairing);
    pairing_apply(T, piH, g2Z, pairing);
    element_mul(R, R, T);

    fmt_kv_e("e(piA, piB)", L);
    fmt_kv_e("RHS",          R);

    int ok = (element_cmp(L, R) == 0);
    fmt_kv_s("result", ok ? "ACCEPT" : "REJECT");

    // cleanup (demo)
    element_clear(L); element_clear(R); element_clear(T);
    element_clear(g2Z); element_clear(tmp);
    for (uint32_t i = 0; i <= m; i++) element_clear(g2_tau[i]);
    free(g2_tau);
    for (uint32_t i = 0; i <= m; i++) element_clear(coef[i]);
    free(coef);
    element_clear(piA); element_clear(piB); element_clear(piC); element_clear(piH); element_clear(g2);
    return ok ? 0 : 1;
}
