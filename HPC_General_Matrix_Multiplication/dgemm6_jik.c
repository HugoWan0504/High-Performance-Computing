void dgemm6_jik(double *C, double *A, double *B, int n)
{
    int j, i, k;
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            for (k = 0; k < n; k++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm6_jik2(double *C, double *A, double *B, int n)
{
    int j, jj, i, ii, k, kk;
    int b = 1;  // Adjust block size
    for (j = 0; j < n; j += b)
        for (i = 0; i < n; i += b)
            for (k = 0; k < n; k += b)
                for (jj = j; jj < j + b; jj++)
                    for (ii = i; ii < i + b; ii++)
                        for (kk = k; kk < k + b; kk++)
                            C[ii * n + jj] += A[ii * n + kk] * B[kk * n + jj];
}
