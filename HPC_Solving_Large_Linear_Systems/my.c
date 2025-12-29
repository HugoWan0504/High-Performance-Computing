#ifndef __MY_C__
#define __MY_C__

#include "include.h"

// Q2
int mydgetrf(double *A,int *ipiv,int n)
{
    //TODO
    //The return value (an integer) can be 0 or 1
    //If 0, the matrix is irreducible and the result will be ignored
    //If 1, the result is valid
    int i, j, k, maxind;
    double max, temp;
    for (k = 0; k < n; k++) {
        maxind = k;
        max = fabs(A[k * n + k]);
        for (i = k + 1; i < n; i++) {
            if (fabs(A[i * n + k]) > max) {
                max = fabs(A[i * n + k]);
                maxind = i;
            }
        }

        if (max == 0) {
            return 0;
        }

        if (maxind != k) {
            for (j = 0; j < n; j++) {
                temp = A[k * n + j];
                A[k * n + j] = A[maxind * n + j];
                A[maxind * n + j] = temp;
            }
            int temp_pivot = ipiv[k];
            ipiv[k] = ipiv[maxind];
            ipiv[maxind] = temp_pivot;
        }

        for (i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];
            for (j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
    return 1;
}

// Q2
void mydtrsv(char UPLO,double *A,double *B,int n,int *ipiv)
{
    //TODO
    int i, j;
    double temp;

    double *B_pivoted = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        B_pivoted[i] = B[ipiv[i]];
    }

    if (UPLO == 'L') {
        for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
                B_pivoted[i] -= A[i * n + j] * B_pivoted[j];
            }
            B_pivoted[i] /= A[i * n + i];
        }
    } else if (UPLO == 'U') {
        for (i = n - 1; i >= 0; i--) {
            for (j = i + 1; j < n; j++) {
                B_pivoted[i] -= A[i * n + j] * B_pivoted[j];
            }
            B_pivoted[i] /= A[i * n + i];
        }
    }

    for (i = 0; i < n; i++) {
        B[i] = B_pivoted[i];
    }
    free(B_pivoted);
}


/*
// Q3   Not Needed But It's Here
void blocked_gepp(double *A, int n, int b) {
    for (int ib = 0; ib < n; ib += b) {
        int end = ib + b > n ? n : ib + b;

        // Perform LU on block A(ib:end, ib:end)
        mydgetrf(A, ib, end, n); 

        for (int i = end; i < n; i++) {
            for (int j = ib; j < end; j++) {
                A[i * n + j] /= A[j * n + j];
            }
        }

        for (int i = ib; i < end; i++) {
            for (int j = end; j < n; j++) {
                A[i * n + j] -= A[i * n + ib] * A[ib * n + j];
            }
        }

        for (int i = end; i < n; i++) {
            for (int j = end; j < n; j++) {
                for (int k = ib; k < end; k++) {
                    A[i * n + j] -= A[i * n + k] * A[k * n + j];
                }
            }
        }
    }
}
*/

// Given
void my_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif