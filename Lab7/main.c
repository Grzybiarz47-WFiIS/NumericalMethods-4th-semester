#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define xmin -5.
#define xmax 5.
#define maxNodes 20
#define step 5
#define stepW 0.01

typedef unsigned int uint;

double f(double x)
{
	return 1 / (1 + x*x);
}

void equal(double* xm, const uint size)
{
	double h = (xmax - xmin) / size;
	for(int i = 0; i <= size; i++)
		xm[i] = xmin + i*h;
}

void chebyshev(double* xm, const uint size)
{
	for(int i = 0; i <= size; i++)
		xm[i] = 0.5*((xmin - xmax)*cos((M_PI*(2*i + 1))/(2*size + 2)) + (xmin + xmax));
}

double Newton(double** fm, double* xm, double x, const uint size)
{
	double res = 0.0;
	for(int j = 0; j <= size; j++)
	{
		double factor = 1.0;
		for(int i = 0; i < j; i++)
			factor *= (x - xm[i]);
		res += fm[j][j]*factor;
	}
	return res;
}
void interpolation(void (*prepare)(double*, const uint), FILE *fp)
{
	for(int n = step; n <= maxNodes; n += step)
	{
		double* xm = malloc((n + 1)*sizeof(double));
		prepare(xm, n);
		double ** fm = malloc((n + 1)*sizeof(double*));
		for(int i = 0; i <= n; i++)
		{	
			fm[i] = malloc((n + 1)*sizeof(double));
			fm[i][0] = f(xm[i]);
		}
		//////////////////////////
		for(int j = 1; j <= n; j++)
			for(int i = j; i <= n; i++)
				fm[i][j] = (fm[i][j - 1] - fm[i - 1][j - 1])/(xm[i] - xm[i - j]);
		//////////////////////////
		for(double x = xmin; x <= xmax; x += stepW)
				fprintf(fp, "%g %g\n", x, Newton(fm, xm, x, n));
		fprintf(fp, "\n\n");
		//////////////////////////
		for(int i = 0; i <= n; i++)
			free(fm[i]);
		free(fm);
		free(xm);
	}
}

int main()
{
	FILE* fp;
	fp = fopen("zad_1.dat", "w");
	interpolation(equal, fp);
	fclose(fp);
	fp = fopen("zad_2.dat", "w");
	interpolation(chebyshev, fp);
	fclose(fp);
	return 0;
}
