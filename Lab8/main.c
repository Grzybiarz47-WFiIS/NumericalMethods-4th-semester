#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

#define xmax 5.0
#define xmin -5.0
#define alpha 0.0
#define beta 0.0
#define step 0.01

float f1(float x)
{
	return 1.0 / (1.0 + x*x);
}

float f2(float x)
{
	return cos(2*x);	
}

float** createNodes(int n)
{
	float** xw = matrix(1, n, 1, 1);
	float h = (xmax - xmin) / (n - 1);
	for(int i = 1; i <= n; i++)
		xw[i][1] = xmin + (i - 1)*h;
	return xw;
}

float** createValues(int n, float (*f)(float), float** xw)
{
	float** yw = matrix(1, n, 1, 1);
	for(int i = 1; i <= n; i++)
		yw[i][1] = f(xw[i][1]);
	return yw;
}

float countM(float x3, float x2, float x1, float y3, float y2, float y1)
{
	float h1 = x2 - x1;
	float h2 = x3 - x2;
	return (6.0 / (h1 + h2))*(((y3 - y2) / h2) - ((y2 - y1) / h1)); 
}

float countLambda(float x3, float x2, float x1)
{
	float h1 = x2 - x1;
	float h2 = x3 - x2;
	return h2 / (h1 + h2);
}

float countNum(float x, float (*f)(float))
{
	return (f(x + step) + f(x - step) - 2*f(x)) / (step*step);
}

void wyzM(float** xw, float** yw, float** m, int n)
{
	float** M = matrix(1, n, 1, n);

	for(int i = 1; i <= n; i++)
		for(int j = 1; j <= n; j++)
			M[i][j] = 0.0;

	M[1][1] = M[n][n] = 1.0;
	m[1][1] = alpha;
	m[n][1] = beta;

	for(int i = 2; i <= n - 1; i++)
	{
		m[i][1] = countM(xw[i + 1][1], xw[i][1], xw[i - 1][1], yw[i + 1][1], yw[i][1], yw[i - 1][1]);
		
		float lambda = countLambda(xw[i + 1][1], xw[i][1], xw[i - 1][1]);
		
		M[i][i] = 2.0;
		M[i][i - 1] = 1.0 - lambda;
		M[i][i + 1] = lambda;
	}
	
	gaussj(M, n, m, 1);

	free_matrix(M, 1, n, 1, n);
}

float wyzSx(float** xw, float** yw, float** m, int n, float x)
{
	int a = n;
	for(int i = 2; i <= n; i++)
		if(x >= xw[i - 1][1] && x <= xw[i][1])
		{
			a = i;
			break;
		}
	
	float h = xw[a][1] - xw[a - 1][1];
	float A = ((yw[a][1] - yw[a - 1][1]) / h) - (h / 6.0)*(m[a][1] - m[a - 1][1]);
	float B = yw[a - 1][1] - m[a - 1][1]*((h*h) / 6.0);
	
	return m[a - 1][1]*(pow(xw[a][1] - x, 3) / (6.0*h)) + m[a][1]*(pow(x - xw[a - 1][1], 3) / (6.0*h)) + A*(x - xw[a - 1][1]) + B;
}

void interpolation(int n, float (*f)(float), FILE* fp, int flag)
{
	float** xw = createNodes(n);
	float** yw = createValues(n, f, xw);
	float** m = matrix(1, n, 1, 1);

	wyzM(xw, yw, m, n);
	
	if(flag)
	{
		for(float x = xmin; x <= xmax; x+= step)
			fprintf(fp, "%g %g\n", x, wyzSx(xw, yw, m, n, x));
		fprintf(fp, "\n\n");
	}
	else
	{
		for(int i = 1; i <= n; i++)
			fprintf(fp, "%g %g %g\n", xw[i][1], m[i][1], countNum(xw[i][1], f));
		fprintf(fp, "\n\n");
	}
	
	free_matrix(xw, 1, n, 1, 1);
	free_matrix(yw, 1, n, 1, 1);
	free_matrix(m, 1, n, 1, 1);
}

int main()
{	
	FILE* fp;
	fp = fopen("f1.dat", "w");
	interpolation(5, f1, fp, 1);
	interpolation(8, f1, fp, 1);
	interpolation(21, f1, fp, 1);
	fclose(fp);
	fp = fopen("f2.dat", "w");
	interpolation(5, f2, fp, 1);
	interpolation(8, f2, fp, 1);
	interpolation(21, f2, fp, 1);
	fclose(fp);
	fp = fopen("pochodne.dat", "w");
	interpolation(10, f1, fp, 0);
	fclose(fp);
	return 0;
}






