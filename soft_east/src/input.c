#include "input.h"
#include "cli_util.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *trim(char *s)
{
    char *end;
    while (*s && isspace((unsigned char)*s)) s++;
    end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
    return s;
}

int input_read(Input *in, const char *path)
{
    FILE *f;
    char line[2048];
    int lineno = 0, cap = 32;

    in->e = NULL;
    in->n = 0;
    strncpy(in->path, path, sizeof(in->path) - 1);
    in->path[sizeof(in->path) - 1] = '\0';

    f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "error: cannot open input file '%s'\n", path);
        return -1;
    }
    in->e = (InputEntry *)malloc(sizeof(InputEntry) * (size_t)cap);
    if (!in->e) {
        fclose(f);
        return -1;
    }
    while (fgets(line, sizeof(line), f)) {
        char *p, *hash, *key, *val;
        lineno++;
        hash = strchr(line, '#');
        if (hash) *hash = '\0';
        p = trim(line);
        if (!*p) continue;
        key = p;
        while (*p && !isspace((unsigned char)*p) && *p != '=') p++;
        if (*p) {
            *p = '\0';
            p++;
        }
        val = trim(p);
        if (*val == '=') val = trim(val + 1); /* allow "key = value" */
        if (!*val) {
            fprintf(stderr, "error: %s:%d: keyword '%s' has no value\n", path, lineno, key);
            fclose(f);
            input_free(in);
            return -1;
        }
        if (strlen(key) >= sizeof(in->e[0].key) || strlen(val) >= sizeof(in->e[0].val)) {
            fprintf(stderr, "error: %s:%d: line too long\n", path, lineno);
            fclose(f);
            input_free(in);
            return -1;
        }
        if (in->n == cap) {
            InputEntry *ne;
            cap *= 2;
            ne = (InputEntry *)realloc(in->e, sizeof(InputEntry) * (size_t)cap);
            if (!ne) {
                fclose(f);
                input_free(in);
                return -1;
            }
            in->e = ne;
        }
        strcpy(in->e[in->n].key, key);
        strcpy(in->e[in->n].val, val);
        in->e[in->n].line = lineno;
        in->e[in->n].used = 0;
        in->n++;
    }
    fclose(f);
    return 0;
}

void input_free(Input *in)
{
    free(in->e);
    in->e = NULL;
    in->n = 0;
}

const char *input_str(Input *in, const char *key)
{
    int i;
    for (i = in->n - 1; i >= 0; i--) {
        if (!strcmp(in->e[i].key, key)) {
            in->e[i].used = 1;
            return in->e[i].val;
        }
    }
    return NULL;
}

int input_count(Input *in, const char *key)
{
    int i, c = 0;
    for (i = 0; i < in->n; i++)
        if (!strcmp(in->e[i].key, key)) c++;
    return c;
}

const char *input_nth(Input *in, const char *key, int idx)
{
    int i, c = 0;
    for (i = 0; i < in->n; i++) {
        if (!strcmp(in->e[i].key, key)) {
            if (c == idx) {
                in->e[i].used = 1;
                return in->e[i].val;
            }
            c++;
        }
    }
    return NULL;
}

static int missing(Input *in, const char *key, int required)
{
    if (required) {
        fprintf(stderr, "error: %s: required keyword '%s' is missing\n", in->path, key);
        return -1;
    }
    return 1;
}

static void bad(Input *in, const char *key, const char *val)
{
    fprintf(stderr, "error: %s: bad value for '%s': %s\n", in->path, key, val);
}

int input_double(Input *in, const char *key, int required, double *out)
{
    const char *v = input_str(in, key);
    if (!v) return missing(in, key, required);
    if (parse_number(v, out) != 0) {
        bad(in, key, v);
        return -1;
    }
    return 0;
}

int input_long(Input *in, const char *key, int required, long *out)
{
    const char *v = input_str(in, key);
    char *end = NULL;
    double d;
    long x;
    if (!v) return missing(in, key, required);
    x = strtol(v, &end, 10);
    if (end != v && *end == '\0') {
        *out = x;
        return 0;
    }
    /* accept forms like 1e6 */
    if (parse_number(v, &d) == 0 && d == (double)(long)d) {
        *out = (long)d;
        return 0;
    }
    bad(in, key, v);
    return -1;
}

int input_int(Input *in, const char *key, int required, int *out)
{
    long x = *out;
    int rc = input_long(in, key, required, &x);
    if (rc == 0) *out = (int)x;
    return rc;
}

int input_u64(Input *in, const char *key, int required, uint64_t *out)
{
    const char *v = input_str(in, key);
    char *end = NULL;
    unsigned long long x;
    if (!v) return missing(in, key, required);
    x = strtoull(v, &end, 10);
    if (end == v || *end != '\0') {
        bad(in, key, v);
        return -1;
    }
    *out = (uint64_t)x;
    return 0;
}

int input_list(Input *in, const char *key, int required, double **vals, int *n)
{
    const char *v = input_str(in, key);
    if (!v) return missing(in, key, required);
    if (parse_list(v, vals, n) != 0) {
        bad(in, key, v);
        return -1;
    }
    return 0;
}

void input_warn_unused(const Input *in)
{
    int i;
    for (i = 0; i < in->n; i++)
        if (!in->e[i].used)
            fprintf(stderr, "warning: %s:%d: keyword '%s' is not used by this program\n",
                    in->path, in->e[i].line, in->e[i].key);
}
