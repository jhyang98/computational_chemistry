/*
 * Keyword input file (system.inp).
 *
 *   # comment
 *   key   value          # trailing comment
 *
 * One keyword per line; the value is the rest of the line. A keyword may be
 * repeated where a program accepts several entries (e.g. "window" in MBAR).
 * Numeric values accept fractions "p/q"; list values accept "a,b,c" or
 * "a:step:b".
 */
#ifndef INPUT_H
#define INPUT_H

#include <stdint.h>

typedef struct {
    char key[64];
    char val[1024];
    int line;
    int used;
} InputEntry;

typedef struct {
    InputEntry *e;
    int n;
    char path[512];
} Input;

/* Returns 0 on success, -1 if the file cannot be read or is malformed. */
int input_read(Input *in, const char *path);
void input_free(Input *in);

/* Last occurrence of key, or NULL. Marks it as used. */
const char *input_str(Input *in, const char *key);

/* Number of occurrences of key, and the idx-th one. */
int input_count(Input *in, const char *key);
const char *input_nth(Input *in, const char *key, int idx);

/*
 * Typed getters. If the key is absent: return 1 and leave *out at the
 * default when required == 0, or print an error and return -1 when
 * required == 1. Return 0 when found and parsed, -1 on a parse error.
 */
int input_int(Input *in, const char *key, int required, int *out);
int input_long(Input *in, const char *key, int required, long *out);
int input_u64(Input *in, const char *key, int required, uint64_t *out);
int input_double(Input *in, const char *key, int required, double *out);
int input_list(Input *in, const char *key, int required, double **vals, int *n);

/* Print a warning for every keyword that was never read. */
void input_warn_unused(const Input *in);

#endif
