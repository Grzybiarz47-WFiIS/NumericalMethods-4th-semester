#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_fft_complex.h>

#define T 1.0
#define Tmax (3.0*T)
#define sigma (T/20.0)
#define w (2*M_PI/T)
#define stride 1

double delta()
{
	return (rand() / (RAND_MAX + 1.0)) - 0.5;
}

double funf(double t)
{
	return sin(w*t) + sin(2*w*t) + sin(3*w*t) + delta();
}

double fung(double t)
{
	return (1.0 / (sigma*sqrt(2.0*M_PI)))*exp((-t*t) / (2.0*sigma*sigma));
}

double* createValues(double(*fun)(double), const int N, const double dt)
{
	double* arr = calloc(2*N, sizeof(double));
	for(int i = 0; i < N; i++)
		arr[2*i] = fun(dt*i);
	return arr;
}

double findMax(double* f, const int N)
{
	double max = 0.0;
	for(int i = 0; i < N; i++)
		if(max < fabs(f[2*i]))
			max = fabs(f[2*i]);
	return max / 2.5;
}

void print(FILE* fp, double* f, const int N, const double dt)
{
	for(int i = 0; i < N; i++)
		fprintf(fp, "%g %g\n", dt*i, f[2*i]);
	fprintf(fp, "\n\n");
}

void FFT(const int k, const char* file)
{
	FILE* fp = fopen(file, "w");
	
	const int N = (int)pow(2, k);
	const double dt = Tmax / (double)N;
	
	double * f = createValues(funf, N, dt);
	double * g1 = createValues(fung, N, dt);
	double * g2 = createValues(fung, N, dt);

	print(fp, f, N, dt);

	gsl_fft_complex_radix2_forward(f, stride, N);
	gsl_fft_complex_radix2_forward(g1, stride, N);
	gsl_fft_complex_radix2_backward(g2, stride, N);
	
	for(int i = 0; i < N; i++)
	{
		double a1 = f[2*i];
		double b1 = f[2*i + 1];
		double a2 = g1[2*i] + g2[2*i];
		double b2 = g1[2*i + 1] + g2[2*i + 1];
		f[2*i] = a1*a2 - b1*b2;
		f[2*i + 1] = a1*b2 + a2*b1;
	}	
	
	gsl_fft_complex_radix2_backward(f, stride, N);
	
	double fmax = findMax(f, N);
	for(int i = 0; i < N; i++)
		f[2*i] /= fmax;
	
	print(fp, f, N, dt);
	
	free(f);
    free(g1);
    free(g2);
    
    fclose(fp);
}

int main()
{
	srand((unsigned int)time(NULL));
	FFT(8, "k8.dat");
	FFT(10, "k10.dat");
	FFT(12, "k12.dat");
	return 0;
}

