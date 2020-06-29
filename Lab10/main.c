#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define iter 20
#define step 100
#define N 200
#define x0 5.0
#define y0 5.0
#define arrSize 3
#define xmax 10.0
#define xmin -10.0
#define ymax 10.0
#define ymin -10.0

const int toPrint[arrSize] = {0, 7, 20};

double d_rand(const double min, const double max)
{
	double r = (double)rand() / RAND_MAX;
	r = r * (max - min) + min;
	return r;
}

double f(const double x, const double y)
{
	double a = sin(x)*sin(y);
	double b = exp(- pow(x + M_PI / 2.0, 2) - pow(y - M_PI / 2.0, 2));
	return a - b;
}

double find_deltax(double x)
{
	double deltax;
	do{
		deltax = d_rand(-1, 1);
	}while(x + deltax > xmax || x + deltax < xmin);
	return deltax;
}

double find_deltay(double y)
{
	double deltay;
	do{
		deltay = d_rand(-1, 1);
	}while(y + deltay > ymax || y + deltay < ymin);
	return deltay;
}

void print(double* x, double* y, FILE* fp)
{
	for(int i = 0; i < N; i++)
		fprintf(fp, "%g %g\n", x[i], y[i]);
	fprintf(fp, "\n\n");
}

int main()
{
	srand((unsigned int)time(NULL));
	
	double* x = malloc(sizeof(double)*N);
	double* y = malloc(sizeof(double)*N);
	
	FILE* fpw0, *fpT;
	fpw0 = fopen("w0.dat", "w");
	fpT = fopen("T.dat", "w");
	
	for(int i = 0; i < N; ++i)
	{
		x[i] = x0;
		y[i] = y0;
	}
	
	for(int it = 0; it <= iter; ++it)
	{
		double T = 10 / pow(2, it);
		for(int k = 0; k < step; ++k)
		{
			for(int i = 0; i < N; ++i)
			{
				double deltax = find_deltax(x[i]);
				double deltay = find_deltay(y[i]);
				if(f(x[i] + deltax, y[i] + deltay) < f(x[i], y[i]))
				{
					x[i] += deltax;
					y[i] += deltay;
				}
				else if(d_rand(0, 1) < exp(-(f(x[i] + deltax, y[i] + deltay) - f(x[i], y[i])) / T))
				{
					x[i] += deltax;
					y[i] += deltay;
				}
			}
			fprintf(fpw0, "%g\n", f(x[0], y[0]));
		}
		for(int i = 0; i < arrSize; ++i)
			if(toPrint[i] == it)
				print(x, y, fpT);
	}
	double min = 2E+9;
	double resx = 0.0;
	double resy = 0.0;
	for(int i = 0; i < N; i++)
		if(min > f(x[i], y[i]))
		{
			resx = x[i];
			resy = y[i];
			min = f(resx, resy);
		}
	printf("Minimum %g w punkcie %g %g\n", min, resx, resy);
	
	fclose(fpw0);
	fclose(fpT);
	
	free(x);
	free(y);
	return 0;
}
