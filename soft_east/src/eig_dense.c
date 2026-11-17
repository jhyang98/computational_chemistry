#include "eig_dense.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define JACOBI_MAX_SWEEPS 100

static void sort_desc(int n, double *w, double *z)
{
    int i, j, k;
    for (i = 0; i < n - 1; i++) {
        int best = i;
        for (j = i + 1; j < n; j++)
            if (w[j] > w[best]) best = j;
        if (best != i) {
            double t = w[i];
            w[i] = w[best];
            w[best] = t;
            if (z) {
                for (k = 0; k < n; k++) {
                    double u = z[k * n + i];
                    z[k * n + i] = z[k * n + best];
                    z[k * n + best] = u;
                }
            }
        }
    }
}

int eig_sym_jacobi(int n, double *a, double *w, double *z)
{
    int sweep, p, q, r;

    if (z) {
        memset(z, 0, sizeof(double) * (size_t)n * (size_t)n);
        for (p = 0; p < n; p++) z[p * n + p] = 1.0;
    }
    if (n == 1) {
        w[0] = a[0];
        return 0;
    }

    for (sweep = 0; sweep < JACOBI_MAX_SWEEPS; sweep++) {
        double off = 0.0, diag = 0.0;
        for (p = 0; p < n; p++) {
            diag += a[p * n + p] * a[p * n + p];
            for (q = p + 1; q < n; q++) off += a[p * n + q] * a[p * n + q];
        }
        if (off == 0.0 || off <= 1e-32 * diag) {
            for (p = 0; p < n; p++) w[p] = a[p * n + p];
            sort_desc(n, w, z);
            return sweep;
        }
        for (p = 0; p < n - 1; p++) {
            for (q = p + 1; q < n; q++) {
                double apq = a[p * n + q];
                double app, aqq, theta, t, c, s;
                if (apq == 0.0) continue;
                app = a[p * n + p];
                aqq = a[q * n + q];
                /* skip rotations that cannot change anything at double precision */
                if (fabs(apq) < 1e-300 ||
                    (fabs(apq) < 1e-18 * (fabs(app) + fabs(aqq)) && sweep > 3)) {
                    a[p * n + q] = 0.0;
                    a[q * n + p] = 0.0;
                    continue;
                }
                theta = (aqq - app) / (2.0 * apq);
                if (fabs(theta) > 1e150) {
                    t = 0.5 / theta;
                } else {
                    t = 1.0 / (fabs(theta) + sqrt(theta * theta + 1.0));
                    if (theta < 0.0) t = -t;
                }
                c = 1.0 / sqrt(t * t + 1.0);
                s = t * c;
                for (r = 0; r < n; r++) {
                    double arp, arq;
                    if (r == p || r == q) continue;
                    arp = a[r * n + p];
                    arq = a[r * n + q];
                    a[r * n + p] = a[p * n + r] = c * arp - s * arq;
                    a[r * n + q] = a[q * n + r] = s * arp + c * arq;
                }
                a[p * n + p] = app - t * apq;
                a[q * n + q] = aqq + t * apq;
                a[p * n + q] = a[q * n + p] = 0.0;
                if (z) {
                    for (r = 0; r < n; r++) {
                        double zrp = z[r * n + p];
                        double zrq = z[r * n + q];
                        z[r * n + p] = c * zrp - s * zrq;
                        z[r * n + q] = s * zrp + c * zrq;
                    }
                }
            }
        }
    }
    for (p = 0; p < n; p++) w[p] = a[p * n + p];
    sort_desc(n, w, z);
    return -1;
}
