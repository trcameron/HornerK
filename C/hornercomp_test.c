#include "horner.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpfr.h>
#include <float.h>
/* rand double complex */
double complex rand_dble_cmplx(const double a,const double b){
	double real = (b-a)*((double)rand()/RAND_MAX) + a;
	double imag = (b-a)*((double)rand()/RAND_MAX) + a;
	return real + imag*I;
}
/* rand double complex poly */
double complex* rand_poly_cmplx(const unsigned int deg){
	double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
	for(unsigned int i=0; i<=deg; ++i){
		poly[i] = rand_dble_cmplx(-1,1);
	}
	return poly;
}
/* Main Function */
int main(int argc,char **argv){
	srand(time(NULL));
	clock_t start, stop;
	double hornerKtime, hornerComptime, hornerMPFRtime;
	double complex *poly, resK, resComp, resMPFR, x;
	FILE *myfile = fopen("../csv_files/hornercomp_test.csv","w+");
	fprintf(myfile,"k, deg, HornerK Time, HornerComp Time, HornerMPFR Time\n");
	const int itnum = 100;
	const int k = 2;
	const int prec = (64*k) - round(4*log2(64*k)) + 13;
	for(unsigned int deg=20; deg<=81920; deg*=2){
		hornerKtime = 0; hornerComptime = 0; hornerMPFRtime = 0;
		for(unsigned int it=0; it<itnum; ++it){
			poly = rand_poly_cmplx(deg);
			x = rand_dble_cmplx(-1,1);
			x = x/cabs(x);

			start = clock();
			resK = hornerK_cmplx(poly,x,deg,k);
			stop = clock();
			hornerKtime += (double)(stop-start)/CLOCKS_PER_SEC;

			start = clock();
			resComp = hornerComp_cmplx(poly,x,deg);
			stop = clock();
			hornerComptime += (double)(stop-start)/CLOCKS_PER_SEC;

			start = clock();
			resMPFR = hornerMPFR_cmplx(poly,x,deg,prec);
			stop = clock();
			hornerMPFRtime += (double)(stop-start)/CLOCKS_PER_SEC;

			free(poly);
		}
		fprintf(myfile,"%d, %d, %.4e, %.4e, %.4e\n",k,deg,hornerKtime/itnum,hornerComptime/itnum,hornerMPFRtime/itnum);
	}
	fclose(myfile);
	return 0;
}