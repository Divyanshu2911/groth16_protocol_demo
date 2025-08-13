// src/prover.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <arpa/inet.h>   // htonl
#include <pbc/pbc.h>
#include "../include/circuit.h"
#include "../include/poly.h"
#include "../include/fmt.h"

// Horner eval: out = sum_{i=0..m-1} coeffs[i] * t^i
static void poly_eval(element_t out, element_t *coeffs, int m, element_t t) {
    element_set0(out);
    for (int i = m - 1; i >= 0; i--) {
        element_mul(out, out, t);
        element_add(out, out, coeffs[i]);
    }
}

static int write_exact(FILE *f, const void *p, size_t n) {
    const unsigned char *b = (const unsigned char*)p;
    size_t w = 0;
    while (w < n) {
        size_t r = fwrite(b + w, 1, n - w, f);
        if (r == 0) return -1;
        w += r;
    }
    return 0;
}

static int write_u32(FILE *f, uint32_t v_host) {
    uint32_t v = htonl(v_host);
    return write_exact(f, &v, 4);
}

static int write_elem(FILE *f, element_t e) {
    int len = element_length_in_bytes(e);
    unsigned char *buf = (unsigned char*)malloc(len);
    element_to_bytes(buf, e);
    int ok = (write_u32(f, (uint32_t)len) == 0 && write_exact(f, buf, (size_t)len) == 0);
    free(buf);
    return ok ? 0 : -1;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s pairing.params d x y a0…ad\n", argv[0]);
        return 1;
    }

    fmt_init(1, stdout);
    fmt_banner("Prover (QAP-aware demo)");

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

    // --- parse polynomial inputs ---
    int argi = 2;
    int d = atoi(argv[argi++]);
    element_t x, y; element_init_Zr(x, pairing); element_init_Zr(y, pairing);
    element_set_str(x, argv[argi++], 10);
    element_set_str(y, argv[argi++], 10);
    element_t *coeffs = (element_t*)malloc((d+1) * sizeof(element_t));
    for (int i = 0; i <= d; i++) { element_init_Zr(coeffs[i], pairing); element_set_str(coeffs[i], argv[argi++], 10); }

    // --- build R1CS & wires (witness) ---
    r1cs_t r; element_t *wires;
    build_r1cs(d, coeffs, x, y, &r, &wires, pairing);
    int m = r.n_cons, n = r.n_vars;
    fmt_kv_i("constraints (m)", m);
    fmt_kv_i("variables (n)", n);

    // --- interpolation points τ_k = 1..m ---
    element_t *tau_pts = (element_t*)malloc(sizeof(element_t)*m);
    for (int i = 0; i < m; i++) { element_init_Zr(tau_pts[i], pairing); element_set_si(tau_pts[i], i+1); }

    // --- toxic waste tau (demo) + bases ---
    element_t tau_secret; element_init_Zr(tau_secret, pairing); element_random(tau_secret);
    element_t g1, g2; element_init_G1(g1, pairing); element_random(g1);
    element_init_G2(g2, pairing); element_random(g2);
    fmt_kv_e("tau", tau_secret);
    fmt_kv_e("g1", g1);
    fmt_kv_e("g2", g2);

    // --- column polys and aggregates ---
    element_t *Ae = (element_t*)malloc(sizeof(element_t)*m);
    element_t *Be = (element_t*)malloc(sizeof(element_t)*m);
    element_t *Ce = (element_t*)malloc(sizeof(element_t)*m);
    element_t *pA = (element_t*)malloc(sizeof(element_t)*m);
    element_t *pB = (element_t*)malloc(sizeof(element_t)*m);
    element_t *pC = (element_t*)malloc(sizeof(element_t)*m);

    element_t valA, valB, valC;
    element_init_Zr(valA, pairing); element_init_Zr(valB, pairing); element_init_Zr(valC, pairing);

    element_t Aagg, Bagg, Cagg;
    element_init_Zr(Aagg, pairing); element_set0(Aagg);
    element_init_Zr(Bagg, pairing); element_set0(Bagg);
    element_init_Zr(Cagg, pairing); element_set0(Cagg);

    fmt_sub("Per-variable scalars and aggregation");
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < m; k++) {
            element_init_Zr(Ae[k], pairing); element_set(Ae[k], r.A[k][j]);
            element_init_Zr(Be[k], pairing); element_set(Be[k], r.B[k][j]);
            element_init_Zr(Ce[k], pairing); element_set(Ce[k], r.C[k][j]);
        }
        lagrange_interpolation(pA, tau_pts, Ae, m, pairing);
        lagrange_interpolation(pB, tau_pts, Be, m, pairing);
        lagrange_interpolation(pC, tau_pts, Ce, m, pairing);

        poly_eval(valA, pA, m, tau_secret);
        poly_eval(valB, pB, m, tau_secret);
        poly_eval(valC, pC, m, tau_secret);

        // accumulate A(τ),B(τ),C(τ)
        element_t t; element_init_Zr(t, pairing);
        element_mul(t, wires[j], valA); element_add(Aagg, Aagg, t);
        element_mul(t, wires[j], valB); element_add(Bagg, Bagg, t);
        element_mul(t, wires[j], valC); element_add(Cagg, Cagg, t);
        element_clear(t);

        for (int k = 0; k < m; k++) { element_clear(Ae[k]); element_clear(Be[k]); element_clear(Ce[k]); element_clear(pA[k]); element_clear(pB[k]); element_clear(pC[k]); }
    }

    fmt_kv_e("A(τ)", Aagg);
    fmt_kv_e("B(τ)", Bagg);
    fmt_kv_e("C(τ)", Cagg);

    // --- compute Z(τ) and H(τ) = (A·B − C)/Z(τ) ---
    element_t Ztau; element_init_Zr(Ztau, pairing); element_set1(Ztau);
    for (int k = 0; k < m; k++) {
        element_t diff; element_init_Zr(diff, pairing);
        element_sub(diff, tau_secret, tau_pts[k]);
        element_mul(Ztau, Ztau, diff);
        element_clear(diff);
    }
    element_t D; element_init_Zr(D, pairing);
    element_mul(D, Aagg, Bagg); element_sub(D, D, Cagg);

    element_t invZ; element_init_Zr(invZ, pairing); element_invert(invZ, Ztau);
    element_t Htau; element_init_Zr(Htau, pairing); element_mul(Htau, D, invZ);

    fmt_sub("QAP divisibility (prover)");
    fmt_kv_e("Z(τ)", Ztau);
    fmt_kv_e("H(τ)", Htau);

    // --- proof elements: piA, piB, piC, piH ---
    element_t piA, piB, piC, piH;
    element_init_G1(piA, pairing);
    element_init_G2(piB, pairing);
    element_init_G1(piC, pairing);
    element_init_G1(piH, pairing);

    element_pow_zn(piA, g1, Aagg);
    element_pow_zn(piB, g2, Bagg);
    element_pow_zn(piC, g1, Cagg);
    element_pow_zn(piH, g1, Htau);

    fmt_sub("Proof elements");
    fmt_kv_e("piA (G1)", piA);
    fmt_kv_e("piB (G2)", piB);
    fmt_kv_e("piC (G1)", piC);
    fmt_kv_e("piH (G1)", piH);

    // --- publish g2^{τ^i} for i=0..m (lets verifier build g2^{Z(τ)} ) ---
    element_t *g2_tau = (element_t*)malloc(sizeof(element_t)*(m+1));
    element_t tp; element_init_Zr(tp, pairing); element_set1(tp);
    for (int i = 0; i <= m; i++) {
        element_init_G2(g2_tau[i], pairing);
        element_pow_zn(g2_tau[i], g2, tp);    // g2^{τ^i}
        element_mul(tp, tp, tau_secret);      // τ^{i+1}
    }

    // --- serialize proof: piA,piB,piC,piH,g2, m, g2^{τ^0..τ^m} ---
    const char *proof_path = "proof_demo.bin";
    FILE *pf = fopen(proof_path, "wb");
    if (!pf) { fprintf(stderr, "Error opening '%s' for write: %s\n", proof_path, strerror(errno)); return 1; }

    int ok = 0
      || write_elem(pf, piA)
      || write_elem(pf, piB)
      || write_elem(pf, piC)
      || write_elem(pf, piH)
      || write_elem(pf, g2)
      || write_u32(pf, (uint32_t)m);
    if (ok) { fprintf(stderr, "Error writing proof header\n"); fclose(pf); return 1; }

    for (int i = 0; i <= m; i++) if (write_elem(pf, g2_tau[i])) { fprintf(stderr, "Error writing g2^tau^i\n"); fclose(pf); return 1; }
    fclose(pf);
    fmt_kv_s("proof file", proof_path);

    // preview pairing that SHOULD pass with g2^{Z(τ)} and piH (verifier reconstructs g2^{Z(τ)})
    fmt_sub("Preview done. (Verifier will do final check)");

    // Cleanup minimal (demo)
    element_clear(invZ); element_clear(Htau); element_clear(D); element_clear(Ztau);
    element_clear(Aagg); element_clear(Bagg); element_clear(Cagg);
    element_clear(valA); element_clear(valB); element_clear(valC);
    element_clear(piA); element_clear(piB); element_clear(piC); element_clear(piH);
    element_clear(g1); element_clear(g2); element_clear(tau_secret);
    element_clear(tp);
    for (int i = 0; i < m; i++) element_clear(tau_pts[i]);
    for (int i = 0; i <= m; i++) element_clear(g2_tau[i]);
    free(tau_pts); free(g2_tau);
    free(Ae); free(Be); free(Ce); free(pA); free(pB); free(pC);
    // (r1cs matrices and wires clearing omitted)
    return 0;
}
