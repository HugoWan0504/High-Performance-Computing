void dgemm6_ikj(double *C, double *A, double *B, int n)
{
    int i, k, j;
    for (i = 0; i < n; i++)
        for (k = 0; k < n; k++)
            for (j = 0; j < n; j++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm6_ikj2(double *C, double *A, double *B, int n)
{
    int i, ii, k, kk, j, jj;
    int b = 1;  // Adjust block size
    for (i = 0; i < n; i += b)
        for (k = 0; k < n; k += b)
            for (j = 0; j < n; j += b)
                for (ii = i; ii < i + b; ii++)
                    for (kk = k; kk < k + b; kk++)
                        for (jj = j; jj < j + b; jj++)
                            C[ii * n + jj] += A[ii * n + kk] * B[kk * n + jj];
}
