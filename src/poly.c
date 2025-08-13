// src/poly.c
#include <pbc/pbc.h>
#include <stdlib.h>
#include "../include/poly.h"

void lagrange_interpolation(element_t *out, element_t *tau, element_t *eval,
                            int m, pairing_t pairing)
{
    // init out to 0
    for (int i = 0; i < m; i++)
    {
        element_init_Zr(out[i], pairing);
        element_set0(out[i]);
    }

    element_t *basis = malloc(sizeof(element_t) * m);
    element_t *nb = malloc(sizeof(element_t) * m);

    for (int k = 0; k < m; k++)
    {
        // basis = 1 (degree 0 polynomial)
        for (int i = 0; i < m; i++)
        {
            element_init_Zr(basis[i], pairing);
            element_set0(basis[i]);
        }
        element_set1(basis[0]);

        // multiply ∏_{j≠k} (x − τ[j])
        for (int j = 0; j < m; j++)
        {
            if (j == k)
                continue;

            for (int i = 0; i < m; i++)
            {
                element_init_Zr(nb[i], pairing);
                element_set0(nb[i]);
            }

            // nb[i+1] += basis[i]  (× x)
            for (int i = 0; i < m - 1; i++)
            {
                element_add(nb[i + 1], nb[i + 1], basis[i]);
            }
            // nb[i] += basis[i] * (−τ[j])
            for (int i = 0; i < m; i++)
            {
                if (!element_is0(basis[i]))
                {
                    element_t tmp;
                    element_init_Zr(tmp, pairing);
                    element_neg(tmp, tau[j]);
                    element_mul(tmp, tmp, basis[i]);
                    element_add(nb[i], nb[i], tmp);
                    element_clear(tmp);
                }
            }

            for (int i = 0; i < m; i++)
            {
                element_set(basis[i], nb[i]);
                element_clear(nb[i]);
            }
        }

        // denom = ∏_{j≠k} (τ[k] − τ[j])
        element_t denom;
        element_init_Zr(denom, pairing);
        element_set1(denom);
        for (int j = 0; j < m; j++)
            if (j != k)
            {
                element_t diff;
                element_init_Zr(diff, pairing);
                element_sub(diff, tau[k], tau[j]);
                element_mul(denom, denom, diff);
                element_clear(diff);
            }
        element_t inv;
        element_init_Zr(inv, pairing);
        element_invert(inv, denom);
        element_clear(denom);

        // scale = eval[k] * inv
        element_t scale;
        element_init_Zr(scale, pairing);
        element_mul(scale, eval[k], inv);
        element_clear(inv);

        // out += scale * basis
        for (int i = 0; i < m; i++)
        {
            element_t tmp;
            element_init_Zr(tmp, pairing);
            element_mul(tmp, basis[i], scale);
            element_add(out[i], out[i], tmp);
            element_clear(tmp);
        }
        element_clear(scale);

        for (int i = 0; i < m; i++)
            element_clear(basis[i]);
    }

    free(basis);
    free(nb);
}
