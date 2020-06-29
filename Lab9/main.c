#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define xmin -4.0
#define xmax 4.0
#define x0 2.0
#define w 1.0
#define N 201
#define m 50
#define k 7
#define frand() ((double)rand())/(RAND_MAX+1.0)

const double sigma = (xmax - xmin) / 16.0;
const double step = (xmax - xmin) / (N - 1);

double f(double x)
{
	double a = sin((14.0*M_PI*x)/(xmax - xmin));
	double b = exp(-pow(x - x0, 2.0)/(2.0*sigma*sigma));
	double c = exp(-pow(x + x0, 2.0)/(2.0*sigma*sigma));
	return a*(b + c);
}

double Crand()
{
	return ((frand() - 0.5) / 5.0);
}

double* createNodes()
{
	double* x = malloc(N*sizeof(double));
	for(int i = 0; i < N; i++)
		x[i] = xmin + i*step;
	return x;
}

double* createYszum(double* x)
{
	double* y = malloc(N*sizeof(double));
	for(int i = 0; i < N; i++)
		y[i] = f(x[i]) + Crand();
	return y;
}

double findAlpha(double* x, double* phi)
{
	double num = 0.0;
	double den = 0.0;
	for(int i = 0; i < N; i++)
	{
		num += x[i]*phi[i]*phi[i];
		den += phi[i]*phi[i];
	}
	return num / den;
}

double findBeta(double* x, double* phi1, double* phi2)
{
	double num = 0.0;
	double den = 0.0;
	for(int i = 0; i < N; i++)
	{
		num += x[i]*phi1[i]*phi2[i];
		den += phi1[i]*phi1[i];
	}
	return num / den;
}

double** createGram(double* x)
{
	double** phi = malloc((m + 1)*sizeof(double*));
	for(int i = 0; i <= m; i++)
		phi[i] = malloc(N*sizeof(double));
	
	for(int i = 0; i < N; i++)
		phi[0][i] = 1.0;
	
	double alpha = findAlpha(x, phi[0]);
	double beta;
	for(int i = 0; i < N; i++)
		phi[1][i] = (x[i] - alpha)*phi[0][i];
	
	for(int j = 1; j < m; j++)
	{
		alpha = findAlpha(x, phi[j]);
		beta = findBeta(x, phi[j - 1], phi[j]);
		for(int i = 0; i < N; i++)
			phi[j + 1][i] = (x[i] - alpha)*phi[j][i] - beta*phi[j - 1][i];
	}
	return phi;
}

double* findC(double* y, double** Gram)
{
	double* c = malloc((m + 1)*sizeof(double));
	for(int j = 0; j <= m; j++)
	{
		c[j] = 0.0;
		for(int i = 0; i < N; i++)
			c[j] += y[i]*Gram[j][i];
	}		
	return c;
}

double* findS(double** Gram)
{
	double* s = malloc((m + 1)*sizeof(double));
	for(int j = 0; j <= m; j++)
	{
		s[j] = 0.0;
		for(int i = 0; i < N; i++)
			s[j] += Gram[j][i]*Gram[j][i];
	}		
	return s;
}

void printPolynomials(double* x, double** Gram, FILE* fp)
{
	for(int i = 0; i < N; i++)
	{
		fprintf(fp, "%g ", x[i]);
		for(int j = 0; j < k; j++)
			fprintf(fp, "%g ", Gram[j][i] / Gram[j][0]);
		fprintf(fp, "\n");
	}
}

void printValues(double* x, double* y, FILE* fp)
{
	for(int i = 0; i < N; i++)
		fprintf(fp, "%g %g\n", x[i], y[i]);
}

void printApprox(const int A, double* c, double* s, double* x, double** Gram, FILE* fp)
{
	for(int i = 0; i < N; i++)
	{
		double approx = 0.0;
		for(int j = 0; j <= A; j++)
			approx += Gram[j][i]*(c[j] / s[j]);
		fprintf(fp, "%g %g\n", x[i], approx);
	}
	fprintf(fp, "\n\n");
}

int main()
{
	srand((unsigned int)time(NULL));
	
	double* x = createNodes();
	double* yszum = createYszum(x); 
	double** Gram = createGram(x);
	double* c = findC(yszum, Gram);
	double* s = findS(Gram);
	
	FILE * fp;
	fp = fopen("Gram.dat", "w");
	printPolynomials(x, Gram, fp);
	fclose(fp);
	
	fp = fopen("pkt.dat", "w");
	printValues(x, yszum, fp);
	fclose(fp);
	
	fp = fopen("approx.dat", "w");
	printApprox(10, c, s, x, Gram, fp);
	printApprox(30, c, s, x, Gram, fp);
	printApprox(50, c, s, x, Gram, fp);
	fclose(fp);
	
	free(x);
	free(yszum);
	for(int i = 0; i <= m; i++)
		free(Gram[i]);
	free(Gram);
	free(c);
	free(s);
	return 0;
}
