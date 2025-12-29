void mmm(double *a, double *b, double *c, int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 2)
        for (j = 0; j < n; j += 2)
            for (k = 0; k < n; k += 2) {
                c[i * n + j] += a[i * n + k] * b[k * n + j] + a[i * n + (k + 1)] * b[(k + 1) * n + j];
                c[(i + 1) * n + j] += a[(i + 1) * n + k] * b[k * n + j] + a[(i + 1) * n + (k + 1)] * b[(k + 1) * n + j];
                c[i * n + (j + 1)] += a[i * n + k] * b[k * n + (j + 1)] + a[i * n + (k + 1)] * b[(k + 1) * n + (j + 1)];
                c[(i + 1) * n + (j + 1)] += a[(i + 1) * n + k] * b[k * n + (j + 1)] + a[(i + 1) * n + (k + 1)] * b[(k + 1) * n + (j + 1)];
            }
}
