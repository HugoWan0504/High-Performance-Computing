#ifndef __MY_BLOCK_C__
#define __MY_BLOCK_C__

#include "include.h"

// First Try
void mydgemm(double *A,double *B,int n,int bid,int b)
{
    //TODO
    //Implement a matrix multiplication here following dgemm7 in HW1
    //The first matrix is A[(bid+1)*b:n,bid*b:(bid+1)*b]
    //The second matrix is B[bid*b:(bid+1)*b,(bid+1)*b:n]
    //b is the block size for dgetrf
    int i, j, k;
    for (i = (bid + 1) * b; i < n; i++) {
        for (j = (bid + 1) * b; j < n; j++) {
            for (k = bid * b; k < (bid + 1) * b; k++) {
                A[i * n + j] -= A[i * n + k] * B[k * n + j];
            }
        }
    }
}

int mydgetrf_block(double *A,int *ipiv,int n)
{
    int b = 64;  // Can be changed
    //TODO
    int k, i;
    
    // Factorize the current block (k:k+b-1) x (k:k+b-1)
    for (k = 0; k < n; k += b) {
        
        if (mydgetrf(A + k * n + k, ipiv + k, b) == 0) {
            return 0;
        }

        // Updatek
        for (i = k + b; i < n; i += b) {
            mydgemm(A, A, n, k / b, b);
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