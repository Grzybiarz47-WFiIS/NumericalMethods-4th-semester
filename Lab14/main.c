#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define X0 10
#define eps 14.06
typedef unsigned long ul;

double countMI(double* x, const ul N)
{
	double res = 0.0;
	for(int i = 0; i < N; ++i)
		res += x[i];
	return res / N;
}

double countSIGMA(double* x, const ul N, const double mi)
{
	double res = 0.0;
	for(int i = 0; i < N; ++i)
		res += pow(x[i] - mi, 2.0);
	return sqrt(res / N);
}

double F(const double x, const double mi, const double delta)
{
	if(x > mi)
		return (-1/(delta*delta))*(x*x/2 - mi*x + mi*mi) + x/delta;
	else
		return (-1/(delta*delta))*(-x*x/2 + mi*x) + x/delta;
}

double countDX(const double xmin, const double xmax, const ul K)
{
	return (xmax - xmin)/K; 
}

void printX(FILE* fp, double* x, const ul N)
{
	for(int i = 0; i < N - 1; ++i)
		fprintf(fp, "%g %g\n", x[i], x[i + 1]);
	fprintf(fp, "\n\n");
}

void printHist1(FILE* fp, int* hist, const double xmin, const ul K, const ul N, const double dx)
{
	for(int j = 0; j < K; ++j)	
		fprintf(fp, "%g %g\n", (2*xmin + (2*j + 1)*dx)/2.0, (double)(hist[j])/N);
	fprintf(fp, "\n\n");
}

void printHist2(FILE* fp, int* hist, double* p, const double xmin, const ul K, const ul N, const double dx)
{
	for(int j = 0; j < K; ++j)
		fprintf(fp, "%g %g %g\n", (2*xmin + (2*j + 1)*dx)/2.0, (double)(hist[j])/N, p[j]);
	fprintf(fp, "\n\n");
}

double genU1()
{
	static ul X = X0;
	static ul m = pow(2, 15.0);
	static ul a = 123;
	static ul c = 1;
	
	X = (a*X + c) % m;
	return X/(m + 1.0);
}

double genU2()
{
	static ul X = X0;
	static ul m = pow(2, 32.0);
	static ul a = 69069;
	static ul c = 1;
	
	X = (a*X + c) % m;
	return X/(m + 1.0);
}

void genRand(FILE* fpU, FILE* fpUHist, const ul N, const ul K, double(*gen)())
{
	double* x = malloc(sizeof(double)*N);
	int* hist = calloc(sizeof(int), K);

	const double dx = countDX(0, 1, K);

	for(int i = 0; i < N; ++i)
	{
		x[i] = gen();
		int index = x[i]/dx;
		hist[index]++;
	}
	
	printX(fpU, x, N);
	
	const double mi = countMI(x, N);
	const double sigma = countSIGMA(x, N, mi);
	printf("%g %g\n", mi, sigma);
	
	printHist1(fpUHist, hist, 0, K, N, dx);
	
	free(x);
	free(hist);
}

void genT(FILE* fpTHist, const ul N, const double mi, const double delta, const ul K)
{
	double* p = malloc(sizeof(double)*K);
	int* hist = calloc(sizeof(int), K);
	
	const double dx = countDX(mi - delta, mi + delta, K);
	
	for(int i = 0; i < N; ++i)
	{
		double x = mi + (genU2() + genU2() - 1)*delta;
		int index = (x - (mi - delta))/dx;
		hist[index]++;
	}
		
	double stat = 0.0;
	for(int j = 0; j < K; ++j)
	{
		p[j] =  F(mi - delta + (j + 1)*dx, mi, delta) - F(mi - delta + j*dx, mi, delta);
		stat += pow(hist[j] - N*p[j], 2.0)/(N*p[j]);
	}
	printf("%g\n", stat);
	
	if(stat < eps)
		printf("hipoteza nie zostala odrzucona na poziomie istotnosci alpha = 0.05\n");
	else
		printf("hipoteza zostala odrzucona na poziomie istotnosci alpha = 0.05\n");
	
	printHist2(fpTHist, hist, p, mi - delta, K, N, dx);

	free(p);
	free(hist);
}

int main()
{
	//czesc 1
	FILE* fpU = fopen("U.dat", "w"), *fpUHist = fopen("U_hist.dat", "w"); 
	genRand(fpU, fpUHist, 1E4, 12, genU1);
	genRand(fpU, fpUHist, 1E4, 12, genU2);
	fclose(fpU);
	fclose(fpUHist);
	
	//czesc 2
	FILE* fpTHist = fopen("T_hist.dat", "w");
	genT(fpTHist, 1E3, 4, 3, 10);
	fclose(fpTHist);
	
	return 0;
}
