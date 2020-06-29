#include <stdio.h>

#include "nrutil.h"
#include "nrutil.c"
#include "tred2.c"
#include "pythag.c"
#include "tqli.c"

#define t -0.021
#define m 10
#define nx 20
#define ny 20
#define N (nx*ny)

int main()
{
	float ** H; //alokacja pamieci - macierz
	H = matrix(1, N, 1, N);
	float * d, * e; //alokacja pamieci - wektory
	int * indx;
	d = vector(1, N);
	e = vector(1, N);
	indx = ivector(1, N);
	for(int i = 1; i <= nx; i++) //wypelnianie macierzy H wg wzoru
	{					
		for(int j = 1; j <= ny; j++)
		{
			int l = j + (i - 1)*ny;
			for(int k = 1; k <= N; k++)
				H[l][k] = 0.;
			if(i > 1) 
				H[l][l - ny] = t; 
			if(i < nx)
				H[l][l + ny] = t;
			if(j > 1)
				H[l][l - 1] = t; 
			if(j < ny) 
				H[l][l + 1] = t; 
			H[l][l]= -4*t;
		}
	}
	for(int i = 1; i <= N; i++)
		indx[i] = i;
 	tred2(H, N, d, e); //redukcja do macierzy trojdiagonalnej
	tqli(d, e, N, H); //diagonalizacja macierzy T (zapisanej w d i e)
	for(int l = 1; l <= N - 1; l++) //sortowanie wg wzoru
	{
		for(int k = N; k >= l + 1; k--)
		{
			float e1 = d[k - 1];
			float e2 = d[k];
			int l1 = indx[k - 1];
			int l2 = indx[k];
			if(e2 < e1)
			{ 
				d[k] = e1;
				d[k - 1] = e2;
				indx[k] = l1;
				indx[k - 1] = l2;
			}
		}
	}
	FILE * fp;
	fp=fopen("out.dat","w");
	for(int i = 1; i <= nx; i++) //zapis do pliku wektorów własnych wg wzoru
	{
		for(int j = 1; j <= ny; j++)
		{
			int l = j + (i - 1)*ny;
			fprintf(fp, "%6d %6d ", i, j);
			for(int k = 1; k <= m; k++)
				fprintf(fp," %12.6f ", H[l][indx[k]]);
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	for(int i = 1; i <= m; i++) //wypisywanie m wartosci wlasnych
		fprintf(fp, "%g\n", d[i]);
	fclose(fp);
	free_matrix(H, 1, N, 1, N); //zwalnianie pamieci
	free_vector(d, 1, N);
	free_vector(e, 1, N);
	free_ivector(indx, 1, N);
	return 0;
}
