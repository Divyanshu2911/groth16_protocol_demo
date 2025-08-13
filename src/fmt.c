#include "../include/fmt.h"
#include <stdlib.h>

static int g_color = 1;

/* ANSI codes (only used when g_color==1) */
#define C_RESET "\033[0m"
#define C_BANNER "\033[1;36m"
#define C_SUB "\033[1;33m"
#define C_KEY "\033[1;37m"
#define C_DIM "\033[2;37m"

static const char *maybe(int on, const char *s) { return on ? s : ""; }

void fmt_init(int use_color, FILE *out)
{
    (void)out;
    /* Honor NO_COLOR env if set. */
    const char *nc = getenv("NO_COLOR");
    g_color = (nc && *nc) ? 0 : use_color;
    setvbuf(stdout, NULL, _IOLBF, 0); // line-buffered for nicer streaming
}

void fmt_hr(void)
{
    printf("%s────────────────────────────────────────────────────────%s\n",
           maybe(g_color, C_DIM), maybe(g_color, C_RESET));
}

void fmt_banner(const char *title)
{
    printf("\n%s== %s ==%s\n", maybe(g_color, C_BANNER), title, maybe(g_color, C_RESET));
    fmt_hr();
}

void fmt_sub(const char *title)
{
    printf("%s-- %s --%s\n", maybe(g_color, C_SUB), title, maybe(g_color, C_RESET));
}

void fmt_kv_s(const char *k, const char *v)
{
    printf("%s%-18s%s : %s\n", maybe(g_color, C_KEY), k, maybe(g_color, C_RESET), v);
}

void fmt_kv_i(const char *k, long long v)
{
    printf("%s%-18s%s : %lld\n", maybe(g_color, C_KEY), k, maybe(g_color, C_RESET), v);
}

void fmt_kv_e(const char *k, element_t e)
{
    printf("%s%-18s%s : ", maybe(g_color, C_KEY), k, maybe(g_color, C_RESET));
    element_printf("%B\n", e);
}

void fmt_vec_e(const char *title, element_t *arr, int n)
{
    if (title && *title)
    {
        printf("%s%s%s\n", maybe(g_color, C_SUB), title, maybe(g_color, C_RESET));
    }
    for (int i = 0; i < n; i++)
    {
        printf("  [%d] = ", i);
        element_printf("%B\n", arr[i]);
    }
}
