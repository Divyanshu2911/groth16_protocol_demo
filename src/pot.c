// ---------------------- src/pot.c ----------------------
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <pbc/pbc.h>
#include "pot.h"
#include "fmt.h"

void generate_pot(int deg, pairing_t pairing)
{
    fmt_banner("Powers of Tau");
    // 1) sample secret tau
    element_t tau;
    element_init_Zr(tau, pairing);
    element_random(tau);
    fmt_kv_e("tau (Zr)", tau);

    // 2) compute tau^i in Zr: tp[0]=1; tp[i]=tp[i-1]*tau
    element_t *tp = malloc((deg + 1) * sizeof(element_t));
    for (int i = 0; i <= deg; i++)
        element_init_Zr(tp[i], pairing);
    element_set1(tp[0]);
    for (int i = 1; i <= deg; i++)
        element_mul(tp[i], tp[i - 1], tau);

    fmt_sub("tau^i (Zr)");
    for (int i = 0; i <= deg; i++)
    {
        printf("  i=%d : ", i);
        element_printf("%B\n", tp[i]);
    }

    // 3) pick generators g1 in G1, g2 in G2
    element_t g1, g2;
    element_init_G1(g1, pairing);
    element_random(g1);
    element_init_G2(g2, pairing);
    element_random(g2);
    fmt_kv_e("g1 (G1)", g1);
    fmt_kv_e("g2 (G2)", g2);

    // 4) exponentiate: g1^{tau^i}, g2^{tau^i}
    element_t tmpG1, tmpG2;
    element_init_G1(tmpG1, pairing);
    element_init_G2(tmpG2, pairing);
    fmt_sub("G1 powers");
    for (int i = 0; i <= deg; i++)
    {
        element_pow_zn(tmpG1, g1, tp[i]);
        printf("  i=%d : ", i);
        element_printf("%B\n", tmpG1);
    }
    fmt_sub("G2 powers");
    for (int i = 0; i <= deg; i++)
    {
        element_pow_zn(tmpG2, g2, tp[i]);
        printf("  i=%d : ", i);
        element_printf("%B\n", tmpG2);
    }

    // cleanup
    element_clear(tmpG1);
    element_clear(tmpG2);
    element_clear(g1);
    element_clear(g2);
    for (int i = 0; i <= deg; i++)
        element_clear(tp[i]);
    free(tp);
    element_clear(tau);
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s pairing.params [deg]", argv[0]);
        return 1;
    }
    int deg = (argc > 2 ? atoi(argv[2]) : 8);
    FILE *fp = fopen(argv[1], "r");
    if (!fp)
    {
        fprintf(stderr, "Error opening '%s': %s", argv[1], strerror(errno));
        return 1;
    }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *buf = malloc(sz + 1);
    size_t rd = fread(buf, 1, sz, fp);
    fclose(fp);
    if (rd != (size_t)sz)
    {
        fprintf(stderr, "Short read: %zu of %ld", rd, sz);
        free(buf);
        return 1;
    }
    buf[sz] = '\0';
    pbc_param_t params;
    pbc_param_init_set_buf(params, buf, sz + 1);
    free(buf);
    pairing_t pairing;
    pairing_init_pbc_param(pairing, params);
    generate_pot(deg, pairing);
    return 0;
}