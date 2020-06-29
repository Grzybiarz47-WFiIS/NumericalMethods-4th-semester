#include <stdio.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"
#include "tred2.c"
#include "pythag.c"
#include "tqli.c"

#define iter 8
#define N 7

void printVector(float* v, const char* str)
{
	FILE *fp = fopen("out.dat", str);
	for(int i = 1; i <= N; i++)
		fprintf(fp, "%g\n", v[i]);
	fprintf(fp, "\n");
	fclose(fp);
}
void prepareVector(float *v, float x)
{
	for(int i = 1; i <= N; i++)
		v[i] = x;
}
float** createMatrix()
{
	float **M = matrix(1, N, 1, N);
	for(int i = 1; i <= 7; i++)
		for(int j = 1; j <= N; j++)
			M[i][j] = sqrt(i + j);
	return M;
}
void print(float **M)
{
	for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
		{
			printf("%g ", M[i][j]);
			if(j == N)
				printf("\n");
		}
		printf("\n");
}
void product(float** M, float* v1, float* v2)
{
	prepareVector(v2, 0.);
	for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
			v2[i] += M[i][j]*v1[j];
}
float scalar(float* v1, float* v2)
{
	float res = 0;
	for(int i = 1; i <= N; i++)
		res += v1[i]*v2[i];
	return res;
}
void redHot(float** M, float* v, float lambda)
{
	for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
			M[i][j] -= lambda*v[i]*v[j];
}

int main()
{
	float *d, *e, *f;
	d = vector(1, N);
	e = vector(1, N);
	f = vector(1, N);
	float **A, **B;
	A = createMatrix();
	B = createMatrix();
	///////////////////// - 1 sposob
	tred2(A, N, d, e);
	tqli(d, e, N, A);
	printVector(d, "w");
	//print(A);
	///////////////////// - 2 sposob
	for(int i = 1; i <= N; i++)
	{
		float lambda;
		prepareVector(d, 1.);
		for(int j = 1; j <= iter; j++)
		{
			product(B, d, e);
			lambda = scalar(e, d) / scalar(d, d);
			float e_norm = sqrt(scalar(e, e));
			for(int k = 1; k <= N; k++)
			{
				e[k] /= e_norm;
				d[k] = e[k];
			}
		}
		redHot(B, d, lambda);
		//print(B);
		//printVector(d, "a");
		f[i] = lambda; 
	}
	printVector(f, "a");
	///////////////////// 
	free_vector(d, 1, N);
	free_vector(e, 1, N);
	free_vector(f, 1, N);
	free_matrix(A, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);
	return 0;
}
