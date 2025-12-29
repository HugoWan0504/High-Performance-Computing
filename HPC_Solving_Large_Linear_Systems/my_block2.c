#ifndef __MY_BLOCK_C__
#define __MY_BLOCK_C__

#include "include.h"
#include <math.h>

// Second Try
void mydgemm(double *A, double *B, int n, int bid, int b) {
    int i, j, k;
    int start_i = (bid + 1) * b;
    int end_i = n;
    int start_j = (bid + 1) * b;
    int end_j = n;
    int start_k = bid * b;
    int end_k = (bid + 1) * b;

    for (i = start_i; i < end_i; i += 4) {
        for (j = start_j; j < end_j; j += 4) {
            for (k = start_k; k < end_k; ++k) {
                for (int ii = i; ii < i + 4 && ii < end_i; ++ii) {
                    for (int jj = j; jj < j + 4 && jj < end_j; ++jj) {
                        A[ii * n + jj] -= A[ii * n + k] * B[k * n + jj];
                    }
                }
            }
        }
    }
}


int mydgetrf_block(double *A, int *ipiv, int n) {
    int b = 4096;
    int k, i, j;

    for (k = 0; k < n; k += b) {
        int end_k = (k + b < n) ? k + b : n;

        if (mydgetrf(A + k * n + k, ipiv + k, end_k - k) == 0) {
            return 0;
        }

        for (i = k + b; i < n; i += b) {
            int end_i = (i + b < n) ? i + b : n;

            for (j = k; j < end_k; ++j) {
                for (int ii = i; ii < end_i; ++ii) {
                    A[ii * n + j] /= A[j * n + j];
                }
            }

            for (int ii = i; ii < end_i; ++ii) {
                for (int jj = end_k; jj < n; ++jj) {
                    for (j = k; j < end_k; ++j) {
                        A[ii * n + jj] -= A[ii * n + j] * A[j * n + jj];
                    }
                }
            }
        }
    }

    return 1;
}

void my_block_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf_block(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif