void dgemm6_kij(double *C, double *A, double *B, int n)
{
    int k, i, j;
    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm6_kij2(double *C, double *A, double *B, int n)
{
    int k, kk, i, ii, j, jj;
    int b = 1;  // Adjust block size
    for (k = 0; k < n; k += b)
        for (i = 0; i < n; i += b)
            for (j = 0; j < n; j += b)
                for (kk = k; kk < k + b; kk++)
                    for (ii = i; ii < i + b; ii++)
                        for (jj = j; jj < j + b; jj++)
                            C[ii * n + jj] += A[ii * n + kk] * B[kk * n + jj];
}
