void dgemm6_jki(double *C, double *A, double *B, int n)
{
    int j, k, i;
    for (j = 0; j < n; j++)
        for (k = 0; k < n; k++)
            for (i = 0; i < n; i++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm6_jki2(double *C, double *A, double *B, int n)
{
    int j, jj, k, kk, i, ii;
    int b = 1;  // Adjust block size
    for (j = 0; j < n; j += b)
        for (k = 0; k < n; k += b)
            for (i = 0; i < n; i += b)
                for (jj = j; jj < j + b; jj++)
                    for (kk = k; kk < k + b; kk++)
                        for (ii = i; ii < i + b; ii++)
                            C[ii * n + jj] += A[ii * n + kk] * B[kk * n + jj];
}
