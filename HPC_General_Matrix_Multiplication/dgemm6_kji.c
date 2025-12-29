void dgemm6_kji(double *C, double *A, double *B, int n)
{
    int k, j, i;
    for (k = 0; k < n; k++)
        for (j = 0; j < n; j++)
            for (i = 0; i < n; i++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm6_kji2(double *C, double *A, double *B, int n)
{
    int k, kk, j, jj, i, ii;
    int b = 1;  // Adjust block size
    for (k = 0; k < n; k += b)
        for (j = 0; j < n; j += b)
            for (i = 0; i < n; i += b)
                for (kk = k; kk < k + b; kk++)
                    for (jj = j; jj < j + b; jj++)
                        for (ii = i; ii < i + b; ii++)
                            C[ii * n + jj] += A[ii * n + kk] * B[kk * n + jj];
}
