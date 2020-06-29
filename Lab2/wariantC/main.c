#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"
#include "ludcmp.c"
#include "lubksb.c"

#define N 50
#define Xa (1.0/2.0)
#define Xb 2.0
#define h ((2*Xb)/((float)(N-1)))

float p(float x)
{
	if(x < -Xa || x > Xa || x == 0)
		return 0.0;
	else if(x < 0)
		return -1.0;
	return 1.0;
}
int main()
{
	float **M, *b, *x;
	//int *v;
	

	M = matrix(1, N, 1, N);
	//v = ivector(1, N);
	b = vector(1, N);
	x = vector(1, N);
	//float *ptr = malloc(sizeof(float));
	
	for(int i = 1; i <= N; ++i)
	{
		b[i] = 0.0;
		for (int j = 1; j <= N; ++j)
			M[i][j] = 0.0;
	}
	///////////////////////////////
	
	x[1] = -Xb;
	x[N] = Xb;
	b[1] = 0.0;
	b[N] = 0.0;
	M[1][1] = 1.0;
	M[1][2] = 0.0;
	M[N][N-1] = 0.0;
	M[N][N] = 1.0;
	for(int i = 2; i < N; i++) //serce
	{
		x[i] = -Xb + h*((float)(i-1));
		b[i] = h*h*p(x[i]);
		M[i][i-1] = M[i][i+1] = 1.0;
		M[i][i] = -2.0;
	}
	///////////////////////////////
	
	//ludcmp(M, N, v, ptr); pierwszy sposob
	//lubksb(M, N, v, b);  
	
	////////////////////////////////
	
	for(int i = 2; i <= N; i++) //drugi sposob
	{
		M[i][i-1] = M[i][i-1]/M[i-1][i-1];
		M[i][i] = M[i][i] - M[i][i-1]*M[i-1][i];
		b[i] = b[i] - M[i][i-1]*b[i-1];
	}
	b[N] /= M[N][N];
	for(int i = N - 1; i >= 1; i--)
		b[i] = (b[i] - M[i][i+1]*b[i+1])/M[i][i];
	
	///////////////////////////////
	for(int i = 1; i <= N; i++)
		printf("%g %g\n", x[i], b[i]);
	///////////////////////////////	
	free_matrix(M, 1, N, 1, N);
	free_vector(b, 1, N);
	free_vector(x, 1, N);
	//free_ivector(v, 1, N);
	//free(ptr);

	return 0;
}

