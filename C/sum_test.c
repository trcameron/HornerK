#include "sum.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
double rand_dble(const double a,const double b){
	return (b-a)*((double)rand()/RAND_MAX) + a;
}
double* rand_vec(const unsigned int n){
	double* p = (double*)malloc(n*sizeof(double));
	for(unsigned int i=0; i<n; ++i){
		p[i] = rand_dble(-1,1);
	}
	return p;
}
/* Main Function */
int main(int argc,char **argv){
	srand(time(NULL));
	unsigned int prec;
	double errk, sumk_time, sumMPFR_time, resk, resMPFR;
	double* p;
	clock_t start, stop;
	for(unsigned int k=2; k<=16; ++k){
		printf("k = %d\n",k);
		prec = (64*k) - round(4*log2(64*k)) + 13;
		printf("prec = %d\n",prec);
		errk = 0; sumk_time = 0; sumMPFR_time = 0; 
		for(unsigned int n=100; n<=2500; n+=10){
			p = rand_vec(n);
			
			start = clock();
			resk = sumk(p,n,k);
			stop = clock();
			sumk_time += (double)(stop-start)/CLOCKS_PER_SEC;
			
			start = clock();
			resMPFR = sumMPFR(p,n,prec);
			stop = clock();
			sumMPFR_time += (double)(stop-start)/CLOCKS_PER_SEC;
			
			errk += fabs(resk - resMPFR);
			
			free(p);
		}
		printf("\t sumk_time = %.4e\n",sumk_time);
		printf("\t sumMPFR_time = %.4e\n",sumMPFR_time);
		printf("\t errk = %.4e\n",errk);
	}
	return 0;
}