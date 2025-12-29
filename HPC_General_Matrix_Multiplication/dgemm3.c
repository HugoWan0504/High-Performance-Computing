void dgemm3(double *C, double *A, double *B, int n)
{
    int i, j, k;
    for (i = 0; i < n; i += 2)
        for (j = 0; j < n; j += 2)
            for (k = 0; k < n; k += 2) {

                // Load 16 registers for a 2x2 block of C
                register double c00 = C[i * n + j];
                register double c01 = C[i * n + j + 1];
                register double c10 = C[(i + 1) * n + j];
                register double c11 = C[(i + 1) * n + j + 1];

                register double c20 = C[i * n + (j + 2)];
                register double c21 = C[i * n + (j + 3)];
                register double c30 = C[(i + 1) * n + (j + 2)];
                register double c31 = C[(i + 1) * n + (j + 3)];

                // Perform the multiplication and update C registers
                for (int kk = k; kk < k + 2; kk++) {
                    register double a0 = A[i * n + kk];
                    register double a1 = A[(i + 1) * n + kk];

                    register double b0 = B[kk * n + j];
                    register double b1 = B[kk * n + j + 1];

                    c00 += a0 * b0;
                    c01 += a0 * b1;
                    c10 += a1 * b0;
                    c11 += a1 * b1;

                    register double b2 = B[kk * n + (j + 2)];
                    register double b3 = B[kk * n + (j + 3)];

                    c20 += a0 * b2;
                    c21 += a0 * b3;
                    c30 += a1 * b2;
                    c31 += a1 * b3;
                }

                // Store results back to C
                C[i * n + j] = c00;
                C[i * n + j + 1] = c01;
                C[(i + 1) * n + j] = c10;
                C[(i + 1) * n + j + 1] = c11;

                C[i * n + (j + 2)] = c20;
                C[i * n + (j + 3)] = c21;
                C[(i + 1) * n + (j + 2)] = c30;
                C[(i + 1) * n + (j + 3)] = c31;
            }
}
