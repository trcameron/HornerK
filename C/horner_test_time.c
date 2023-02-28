#include "horner.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpfr.h>
#include <float.h>
/* random double */
double rand_dble(const double a,const double b){
	return (b-a)*((double)rand()/RAND_MAX) + a;
}
/* random polynomial */
double* rand_poly(const unsigned int deg){
	double* poly = (double*)malloc((deg+1)*sizeof(double));
	for(unsigned int i=0; i<=deg; ++i){
		poly[i] = rand_dble(-1,1);
	}
	return poly;
}
/* Main Function */
int main(int argc,char **argv){
	srand(time(NULL));
	clock_t start, stop;
	double hornerKtime, hornerMPFRtime, *poly, resK, resMPFR, x;
	int prec;
	const int itnum = 100;
	FILE *myfile = fopen("../csv_files/horner_test_time.csv","w+");
	fprintf(myfile,"k, prec, HornerK Time, HornerMPFR, Ratio\n");
	for(unsigned int k=2; k<=8; ++k){
		prec = (64*k) - round(4*log2(64*k)) + 13;
		hornerKtime = 0; hornerMPFRtime = 0;
		for(unsigned int deg=20; deg<=81920; deg*=2){
			for(unsigned int it=0; it<itnum; ++it){
				poly = rand_poly(deg);
				x = rand_dble(-1,1);

				start = clock();
				resK = hornerK(poly,x,deg,k);
				stop = clock();
				hornerKtime += (double)(stop-start)/CLOCKS_PER_SEC;

				start = clock();
				resMPFR = hornerMPFR(poly,x,deg,prec);
				stop = clock();
				hornerMPFRtime += (double)(stop-start)/CLOCKS_PER_SEC;

				free(poly);
			}
		}
		fprintf(myfile,"%d, %d, %.4e, %.4e, %.4e\n",k,prec,hornerKtime/(13*itnum),hornerMPFRtime/(13*itnum),hornerMPFRtime/hornerKtime);
	}
	fclose(myfile);
	
	return 0;
}