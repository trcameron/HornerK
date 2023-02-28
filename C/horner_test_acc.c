#include "horner.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpfr.h>
#include <float.h>
/* binomial coefficient */
long long int binom_coeff(const unsigned int n,const unsigned int k){
	long long int c[n+1][k+1];
	for(unsigned int i=0; i<=n; ++i){
		for(unsigned int j=0; j<=fmin(i,k); ++j){
			if(j==0 || j==i){
				c[i][j] = 1;
			}
			else{
				c[i][j] = c[i-1][j-1] + c[i-1][j];
			}
		}
	}
	return c[n][k];
}
/* binomial polynomial */
double* binom_poly(const unsigned int deg){
	double* poly = (double*)malloc((deg+1)*sizeof(double));
	for(unsigned int k=0; k<=deg; ++k){
		poly[k] = pow(-1,deg-k)*(double)binom_coeff(deg,k);
	}
	return poly;
}
/* Main Function */
int main(int argc,char **argv){
	const double x = 220.0/219.0;
	// printf("x: %.17f\n",x);
	double cond, err, *poly, resK, resMPFR;
	const unsigned int prec = 1024;
	FILE *myfile = fopen("../csv_files/horner_test_acc.csv","w+");
	fprintf(myfile,"deg, cond, prec, err1, err2, err3, err4, err5, err6, err7, err8\n");
	for(unsigned int deg=2; deg<=50; ++deg){
		cond = pow(fabs((1+x)/(1-x)),deg);
		poly = binom_poly(deg);

		fprintf(myfile,"%d, %.4e, %d",deg,cond,(int)ceil(log2(cond)));

		resMPFR = hornerMPFR(poly,x,deg,prec);
		resK = horner(poly,x,deg);
		err = fabs(resK - resMPFR)/fabs(resMPFR);
		if(err > 1.0){
			err = 1.0;
		}
		else if(err < DBL_EPSILON){
			err = DBL_EPSILON;
		}
		fprintf(myfile,", %.4e",err);
		for(unsigned int k=2; k<=8; ++k){
			resK = hornerK(poly,x,deg,k);
			err = fabs(resK - resMPFR)/fabs(resMPFR);
			if(err > 1.0){
				err = 1.0;
			}
			else if(err < DBL_EPSILON){
				err = DBL_EPSILON;
			}
			fprintf(myfile,", %.4e",err);
		}
		fprintf(myfile,"\n");
		free(poly);
	}
	fclose(myfile);

	return 0;
}