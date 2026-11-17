#ifndef CLI_UTIL_H
#define CLI_UTIL_H

/*
 * Parse a list of numbers:
 *   "a"          single value
 *   "a,b,c"      explicit list
 *   "a:d:b"      range a, a+d, ..., up to b (inclusive within round-off)
 * Single values may be written as fractions "p/q".
 * On success *out is malloc'ed with *n entries; returns 0. Returns -1 on error.
 */
int parse_list(const char *str, double **out, int *n);

/* Parse one number, allowing a fraction "p/q". Returns 0 on success. */
int parse_number(const char *str, double *val);

#endif
