#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 5
#define ITmax 30

void print(FILE* fp, double L, double iter, double x, double R1, double R2)
{
	fprintf(fp, "%3g %3g %12g %12g %12g\n", L, iter, x, R1, R2);
}
double licz_r(double* a, double* b, int n, double x)
{
	b[n] = 0.0;
	for(int k = n - 1; k >= 0; k--)
		b[k] = a[k + 1] + x*b[k + 1];
	return (a[0] + x*b[0]);
}
int main()
{
	double * a = malloc(sizeof(double) * (N + 1));
	double * b = malloc(sizeof(double) * (N + 1));
	double * c = malloc(sizeof(double) * N);
	a[0] = 240.0;
	a[1] = -196.0;
	a[2] = -92.0;
	a[3] = 33.0;
	a[4] = 14.0;
	a[5] = 1.0;
	FILE * fp;
	fp = fopen("out.dat","w");
	for(int i = 1; i <= N; i++)
	{
		int n = N - i + 1;
		double x0 = 0.0;
		for(int iter = 1; iter <= ITmax; iter++)
		{
			double R1 = licz_r(a, b, n, x0);
			double R2 = licz_r(b, c, n - 1, x0);
			double x1 = x0 - (R1 / R2);
			print(fp, i, iter, x1, R1, R2);
			if(fabs(x1 - x0) < 1E-7)
				break;
			x0 = x1;
		}
		for(int j = 0; j <= (n - 1); j++)
			a[j] = b[j];
		fprintf(fp, "\n");
	}
	fclose(fp);
	free(a);
	free(b);
	free(c); 
	return 0;
}

