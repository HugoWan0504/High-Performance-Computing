void dgemm7(double *C, double *A, double *B, int n) 
{
    int i, j, k, ii, jj, kk;
    
    int b = 64;  // This can be optimized based on the cache size
    
    for (i = 0; i < n; i += b)
        for (j = 0; j < n; j += b)
            for (k = 0; k < n; k += b)

                // Process sub-blocks of A, B, and C
                for (ii = i; ii < i + b && ii < n; ii += 2)
                    for (jj = j; jj < j + b && jj < n; jj += 2) {
                        // Load a 2x2 block of C into registers
                        register double c00 = C[ii * n + jj];
                        register double c01 = C[ii * n + jj + 1];
                        register double c10 = C[(ii + 1) * n + jj];
                        register double c11 = C[(ii + 1) * n + jj + 1];

                        for (kk = k; kk < k + b && kk < n; kk += 2) {
                            // Load elements of A and B into registers for the 2x2 multiplication
                            register double a00 = A[ii * n + kk];
                            register double a01 = A[ii * n + kk + 1];
                            register double a10 = A[(ii + 1) * n + kk];
                            register double a11 = A[(ii + 1) * n + kk + 1];

                            register double b00 = B[kk * n + jj];
                            register double b01 = B[kk * n + jj + 1];
                            register double b10 = B[(kk + 1) * n + jj];
                            register double b11 = B[(kk + 1) * n + jj + 1];

                            // Perform 2x2 block multiplication and accumulate in registers
                            c00 += a00 * b00 + a01 * b10;
                            c01 += a00 * b01 + a01 * b11;
                            c10 += a10 * b00 + a11 * b10;
                            c11 += a10 * b01 + a11 * b11;
                        }

                        // Store the results back into C
                        C[ii * n + jj] = c00;
                        C[ii * n + jj + 1] = c01;
                        C[(ii + 1) * n + jj] = c10;
                        C[(ii + 1) * n + jj + 1] = c11;
                    }
}
