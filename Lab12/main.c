#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define a 0.
#define b 1.
#define n 8

double f(double x)
{
	return log(pow(x, 3.) + 3*pow(x, 2.) + x + 0.1)*sin(18*x);
}

void prepareValues(double* fVal, int N, double h)
{
	for(int i = 0; i <= N; i++)
		fVal[i] = f(a + i*h);
}

void Simpson(double** D)
{
	for(int i = 0; i <= n; i++)
	{
		int N = pow(2, i + 1); 
		double h = (b - a)/N;
		double sum = 0.;
		
		double* fVal = malloc((N + 1)*sizeof(double));
		prepareValues(fVal, N, h);
		
		for(int i = 0; i < N/2; i++)
			sum += (h/3)*(fVal[2*i] + 4*fVal[2*i + 1] + fVal[2*i + 2]);
		
		D[i][0] = sum; 
		
		free(fVal);
	}
}

void Milne(double** D)
{
	for(int i = 0; i <= n; i++)
	{
		int N = pow(2, i + 2); 
		double h = (b - a)/N;
		double sum = 0.;
		
		double* fVal = malloc((N + 1)*sizeof(double));
		prepareValues(fVal, N, h);
		
		for(int i = 0; i < N/4; i++)
			sum += (4*h/90)*(7*fVal[4*i] + 32*fVal[4*i + 1] + 12*fVal[4*i + 2] + 32*fVal[4*i + 3] + 7*fVal[4*i + 4]);
		
		D[i][0] = sum; 
		
		free(fVal);
	}
}

void findIntegral(void(*method)(double**), FILE* fp)
{
	double ** D = malloc((n + 1)*sizeof(double*));
	for(int i = 0; i <= n; i++)
		D[i] = malloc((n + 1)*sizeof(double));
	
	method(D);
	
	for(int i = 1; i <= n; i++)
		for(int j = 1; j <= i; j++)
			D[i][j] = (pow(4, j)*D[i][j - 1] - D[i - 1][j - 1])/(pow(4, j) - 1);
	
	for(int i = 0; i <= n; i++)
		fprintf(fp, "%.10g %.10g\n", D[i][0], D[i][i]);
	fprintf(fp, "\n\n");
	
	for(int i = 0; i <= n; i++)
		free(D[i]);
	free(D);
}


int main()
{
	FILE* fp = fopen("out.dat", "w");
	findIntegral(Simpson, fp);
	findIntegral(Milne, fp);
	fclose(fp);
	return 0;
}



