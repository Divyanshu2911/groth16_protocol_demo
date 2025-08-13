// ---------------------- include/pot.h ----------------------
#ifndef POT_H
#define POT_H

#include <pbc/pbc.h>

// Compute and print powers of tau up to degree 'deg'
// tau is random in Zr
void generate_pot(int deg, pairing_t pairing);

#endif // POT_H