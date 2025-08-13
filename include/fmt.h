#ifndef FMT_H
#define FMT_H

#include <stdio.h>
#include <pbc/pbc.h>

/* Initialize formatting.
 * use_color: 1=ANSI colors, 0=no color (or set NO_COLOR env var).
 * out: currently unused (PBC prints to stdout); keep stdout for now.
 */
void fmt_init(int use_color, FILE *out);

void fmt_hr(void);                  // horizontal rule
void fmt_banner(const char *title); // big section header
void fmt_sub(const char *title);    // subsection header

void fmt_kv_s(const char *k, const char *v);
void fmt_kv_i(const char *k, long long v);
void fmt_kv_e(const char *k, element_t e); // pretty-print a field/group element

/* Print a vector of elements with indices:
 * title: printed once above the list.
 */
void fmt_vec_e(const char *title, element_t *arr, int n);

#endif
