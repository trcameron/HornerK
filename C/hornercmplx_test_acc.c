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
/* power of I */
double complex ipow(const unsigned int n){
	if(n%4 == 0){
		return 1;
	}
	if(n%4 == 1){
		return I;
	}
	if(n%4 == 2){
		return -1;
	}
	return -I;
}
/* binomial polynomial */
double complex* binom_poly_cmplx(const unsigned int deg){
	double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
	for(unsigned int k=0; k<=deg; ++k){
		poly[k] = pow(-1,deg-k)*ipow(deg-k)*(double)binom_coeff(deg,k);
	}
	return poly;
}
/* Main Function */
int main(int argc,char **argv){
	const double complex x = I*220.0/219.0, y = I;
	double cond, err;
	double complex *poly, resK, resMPFR;
	const unsigned int prec = 1024;
	FILE *myfile = fopen("../csv_files/hornercmplx_test_acc.csv","w+");
	fprintf(myfile,"deg, cond, prec, err1, err2, err3, err4, err5, err6, err7, err8\n");
	for(unsigned int deg=2; deg<=50; ++deg){
		cond = pow((1+cabs(x))/cabs(x-y),deg);
		poly = binom_poly_cmplx(deg);

		fprintf(myfile,"%d, %.4e, %d",deg,cond,(int)ceil(log2(cond)));

		resMPFR = hornerMPFR_cmplx(poly,x,deg,prec);
		resK = horner_cmplx(poly,x,deg);
		err = cabs(resK - resMPFR)/cabs(resMPFR);
		if(err > 1.0){
			err = 1.0;
		}
		else if(err < DBL_EPSILON){
			err = DBL_EPSILON;
		}
		fprintf(myfile,", %.4e",err);
		for(unsigned int k=2; k<=8; ++k){
			resK = hornerK_cmplx(poly,x,deg,k);
			err = cabs(resK - resMPFR)/cabs(resMPFR);
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