#include "math.h"
#include "stdio.h"

#include "nrutil.h"
#include "nrutil.c"
#include "gauleg.c"
#include "gaulag.c"
#include "gammln.c"
#include "gauher.c"

#define alpha 0.0

double f1(double x)
{
	return 1.0/(x*sqrt(x*x - 1));
}

double f2(double x)
{
	return log(fabs(x))*exp(-x*x);
}

double f3(double x)
{
	return sin(2*x)*exp(-3*x);
}

void GauLeg(FILE* fp, double(*f)(double), const int N, const float a, const float b, const double c)
{
	for(int n = 2; n <= N; ++n)
	{
		float* x = vector(1, n);
		float* w = vector(1, n);
		gauleg(a, b, x, w, n);
		
		double sum = 0.0;
		for(int i = 1; i <= n; ++i)
			sum += w[i]*f(x[i]);
		
		fprintf(fp, "%d %g\n", n, fabs(c - sum));
			
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	fprintf(fp, "\n\n");
}

void GauLag(FILE* fp, double(*f)(double), const int N, const double c)
{
	for(int n = 2; n <= N; ++n)
	{
		float* x = vector(1, n);
		float* w = vector(1, n);
		gaulag(x, w, n, alpha);
		
		double sum = 0.0;
		for(int i = 1; i <= n; ++i)
			sum += (w[i]*f(x[i]))/(exp(-x[i]));
		
		fprintf(fp, "%d %g\n", n, fabs(c - sum));
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	fprintf(fp, "\n\n");
}

void GauHer(FILE* fp, double(*f)(double), const int N, const double c)
{
	for(int n = 2; n <= N; n += 2)
	{
		float* x = vector(1, n);
		float* w = vector(1, n);
		gauher(x, w, n);
		
		double sum = 0.0;
		for(int i = 1; i <= n; ++i)
			sum += (w[i]*f(x[i]))/(2*exp(-x[i]*x[i]));
		
		fprintf(fp, "%d %g\n", n, fabs(c - sum));
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	fprintf(fp, "\n\n");
}


int main()
{
	FILE* fp = fopen("out.dat", "w");
	
	//podpunkt a
	GauLeg(fp, f1, 100, 1, 2, M_PI/3);
	
	//podpunkt b
	GauHer(fp, f2, 100, -0.8700577);
	GauLeg(fp, f2, 100, 0, 5, -0.8700577);
	
	//podpunkt c
	GauLag(fp, f3, 20, 2.0/13.0);
	
	fclose(fp);
	return 0;
}

