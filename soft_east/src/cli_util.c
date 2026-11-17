#include "cli_util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int parse_number(const char *str, double *val)
{
    char *end = NULL;
    double a = strtod(str, &end);
    if (end == str) return -1;
    if (*end == '/') {
        const char *p = end + 1;
        double b = strtod(p, &end);
        if (end == p || *end != '\0' || b == 0.0) return -1;
        *val = a / b;
        return 0;
    }
    if (*end != '\0') return -1;
    *val = a;
    return 0;
}

static int count_char(const char *s, char ch)
{
    int n = 0;
    for (; *s; s++)
        if (*s == ch) n++;
    return n;
}

int parse_list(const char *str, double **out, int *n)
{
    char *buf = NULL;
    double *vals = NULL;
    int cnt;

    *out = NULL;
    *n = 0;
    if (!str || !*str) return -1;
    buf = (char *)malloc(strlen(str) + 1);
    if (!buf) return -1;
    strcpy(buf, str);

    if (count_char(buf, ':') == 2) {
        char *p1 = strchr(buf, ':');
        char *p2 = strchr(p1 + 1, ':');
        double a, d, b, scale;
        int i;
        *p1 = '\0';
        *p2 = '\0';
        if (parse_number(buf, &a) || parse_number(p1 + 1, &d) || parse_number(p2 + 1, &b) ||
            d == 0.0 || (b - a) / d < -1e-9) {
            free(buf);
            return -1;
        }
        cnt = (int)floor((b - a) / d + 1e-9) + 1;
        vals = (double *)malloc(sizeof(double) * (size_t)cnt);
        if (!vals) {
            free(buf);
            return -1;
        }
        scale = fabs(a) + fabs(b) + fabs(d);
        for (i = 0; i < cnt; i++) {
            double v = a + i * d;
            if (fabs(v) < 1e-12 * scale) v = 0.0;
            vals[i] = v;
        }
    } else if (count_char(buf, ':') == 0) {
        char *tok;
        int i = 0;
        cnt = count_char(buf, ',') + 1;
        vals = (double *)malloc(sizeof(double) * (size_t)cnt);
        if (!vals) {
            free(buf);
            return -1;
        }
        tok = strtok(buf, ",");
        while (tok) {
            if (parse_number(tok, &vals[i]) != 0) {
                free(vals);
                free(buf);
                return -1;
            }
            i++;
            tok = strtok(NULL, ",");
        }
        if (i != cnt) {
            free(vals);
            free(buf);
            return -1;
        }
    } else {
        free(buf);
        return -1;
    }
    free(buf);
    *out = vals;
    *n = cnt;
    return 0;
}
