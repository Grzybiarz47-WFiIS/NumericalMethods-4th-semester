#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"
#include "ludcmp.c"
#include "lubksb.c"

#define N 3

void prepare(float ** M)
{
	M[1][1] = 1.0;
	M[1][2] = 2.0;
	M[1][3] = 3.0;
	M[2][1] = 4.0;
	M[2][2] = 5.0;
	M[2][3] = 6.0;
	M[3][1] = 7.0;
	M[3][2] = 8.0;
	M[3][3] = 9.0;
}
void rotate(float ** M, float ** M1, int * indx)
{
	float * b = vector(1, N);
	for(int i = 1; i <= N; i++)	
	{
		for(int i = 1; i <= N; i++)
			b[i] = 0.0;
		b[i] = 1.0;
		lubksb(M, N, indx, b);
		for(int j = 1; j <= N; j++)
			M1[j][i] = b[j];
	}
	free_vector(b, 1, N);
}
void print(float ** A)
{
	for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
		{
			printf("%g ", A[i][j]);
			if(j == N)
				printf("\n");
		}
	printf("////////////////////////////\n");
}
int main()
{
	float ** A, ** B, ** CA, ** CB;
	int * indx = ivector(1, N);
	float * ptr = malloc(sizeof(float));
	A = matrix(1, N, 1, N);
	B = matrix(1, N, 1, N);
	CA = matrix(1, N, 1, N);
	CB = matrix(1, N, 1, N);
	////////////////////
	prepare(A);
	prepare(B);
	prepare(CA);
	prepare(CB);
	B[1][1] = 1.1;
	CB[1][1] = 1.1;
	////////////////////
	ludcmp(A, N, indx, ptr);
	ludcmp(B, N, indx, ptr);
	////////////////////
	float ** A1, ** B1;
	A1 = matrix(1, N, 1, N);
	B1 = matrix(1, N, 1, N);
	rotate(A, A1, indx);
	rotate(B, B1, indx);
	////////////////////
	float ** A2 = matrix(1, N, 1, N);
	float ** B2 = matrix(1, N, 1, N);
	print(CA);
	print(A1);
	print(CB);
	print(B1);
	for(int i = 1; i <= N; i++)
	{
		for(int j = 1; j <= N; j++)
		{
			A2[i][j] = 0.0;
			B2[i][j] = 0.0;
			for(int k = 1; k <= N; k++)
			{
				A2[i][j] += A1[i][k] * CA[k][j];
				B2[i][j] += B1[i][k] * CB[k][j]; 
			}
		}
	}
	print(A2);
	print(B2);
	////////////////////
	free_matrix(A, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);
	free_matrix(A1, 1, N, 1, N);
	free_matrix(B1, 1, N, 1, N);
	free_matrix(A2, 1, N, 1, N);
	free_matrix(B2, 1, N, 1, N);
	free_matrix(CA, 1, N, 1, N);
	free_matrix(CB, 1, N, 1, N);
	free_ivector(indx, 1, N);
	free(ptr);
	
	return 0;
}

