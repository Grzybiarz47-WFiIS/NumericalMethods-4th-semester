#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

#define w 1.
#define N 2000
#define h 0.02

void complete(double * v0, double * v1, double * v2, double * b, double beta, double omega, double F)
{
	double tmp1 = beta*h + 1;
	double tmp2 = w*w*h*h - 2 - beta*h;
	for(int i = 2; i <= N; ++i)
	{
		v0[i] = tmp1;
		v1[i] = tmp2;
		v2[i] = 1;
		b[i] = F*sin(omega*h*(i - 1))*h*h;
	}
}
void jakobi(double * v0, double * v1, double * v2, double * b, double * xs, double * xn)
{
	for(int it = 1; it < 1E+6; ++it)
	{
		xn[0] = b[0] / v0[0];
		xn[1] = (b[1] - v1[1]*xs[0]) / v0[1];
		double sn = xn[0]*xn[0] + xn[1]*xn[1];
		double ss = xs[0]*xs[0] + xs[1]*xs[1];
		for(int i = 2; i <= N; ++i)
		{
			xn[i] = (b[i] - v1[i]*xs[i - 1] - v2[i]*xs[i - 2]) / v0[i];
			sn += xn[i]*xn[i];
			ss += xs[i]*xs[i];
		}
		if(fabs(sn - ss) < 1E-6)
		{
			//printf("Iteracje: %d\n", it);
			break;
		}
		for(int i = 0; i <= N; ++i)
			xs[i] = xn[i];
	}
	for(int i = 0; i <= N; ++i)
		printf("%g %g\n", h*i, xn[i]);
	printf("\n\n");
}
int main()
{
	double * vec_0 = malloc(sizeof(double) * (N + 1)); //inicjalizacja
	double * vec_1 = malloc(sizeof(double) * (N + 1));
	double * vec_2 = malloc(sizeof(double) * (N + 1));
	double * b = malloc(sizeof(double) * (N + 1));
	double * xs = calloc(N + 1, sizeof(double));
	double * xn = malloc(sizeof(double) * (N + 1));
	//////////
	vec_0[0] = vec_0[1] = 1; //poczatkowe wartosci
	vec_1[0] = 0;
	vec_1[1] = -1;
	vec_2[0] = vec_2[1] = 0;
	b[0] = 1;
	b[1] = 0;
	////////// 
	complete(vec_0, vec_1, vec_2, b, 0.0, 0.8, 0.0); //pierwsza wersja
	jakobi(vec_0, vec_1, vec_2, b, xs, xn);
	////////// 
	complete(vec_0, vec_1, vec_2, b, 0.4, 0.8, 0.0); //druga wersja
	jakobi(vec_0, vec_1, vec_2, b, xs, xn);
	////////// 
	complete(vec_0, vec_1, vec_2, b, 0.4, 0.8, 0.1); //trzecia wersja
	jakobi(vec_0, vec_1, vec_2, b, xs, xn);
	//////////
	free(vec_0); //zwalnianie pamieci
	free(vec_1);
	free(vec_2);
	free(b);
	free(xs);
	free(xn);
	return 0;
}
