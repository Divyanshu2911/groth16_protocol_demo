// ---------------------- src/pot.c ----------------------
#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>
#include "../include/pot.h"

void generate_pot(int deg, pairing_t pairing) {
    // 1) sample secret tau
    element_t tau; element_init_Zr(tau, pairing);
    element_random(tau);
    printf("secret tau (in Zr): %B\n", tau);

    // 2) compute tau^i in Zr
    element_t *tp = malloc((deg+1)*sizeof(element_t));
    for(int i=0;i<=deg;i++){
        element_init_Zr(tp[i], pairing);
        if(i==0) element_set1(tp[i]);
        else     element_pow_zn(tp[i], tau, tp[i-1]);
        element_printf("tau^%d = %B\n", i, tp[i]);
    }

    // 3) pick generators g1 in G1, g2 in G2
    element_t g1, g2, tmp;
    element_init_G1(g1, pairing); element_random(g1);
    element_init_G2(g2, pairing); element_random(g2);
    element_init_G1(tmp, pairing);
    printf("generator g1: %B\n", g1);
    printf("generator g2: %B\n", g2);

    // 4) exponentiate: g1^{tau^i}, g2^{tau^i}
    for(int i=0;i<=deg;i++){
        element_pow_zn(tmp, g1, tp[i]);
        element_printf("g1^{tau^%d} = %B\n", i, tmp);
    }
    element_init_G2(tmp, pairing);
    for(int i=0;i<=deg;i++){
        element_pow_zn(tmp, g2, tp[i]);
        element_printf("g2^{tau^%d} = %B\n", i, tmp);
    }

    // cleanup
    for(int i=0;i<=deg;i++) element_clear(tp[i]); free(tp);
    element_clear(tau);
    element_clear(g1); element_clear(g2); element_clear(tmp);
}

int main(int argc, char **argv) {
    if(argc<2){ fprintf(stderr,"Usage: %s path/to/a.param [deg]\n",argv[0]); return 1; }
    int deg = (argc>2? atoi(argv[2]) : 8);
    FILE *fp = fopen(argv[1],"r"); fseek(fp,0,SEEK_END); long sz=ftell(fp); fseek(fp,0,SEEK_SET);
    char *buf=malloc(sz+1); fread(buf,1,sz,fp); buf[sz]='\0'; fclose(fp);
    pbc_param_t params; pbc_param_init_set_buf(params,buf,sz+1); free(buf);
    pairing_t pairing; pairing_init_pbc_param(pairing,params);
    generate_pot(deg, pairing);
    return 0;
}