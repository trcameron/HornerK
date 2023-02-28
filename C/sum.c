#include "sum.h"
#include "math.h"
#include <mpc.h>
#include <mpfr.h>
#include <stdlib.h>
/* Sumk (summation in k-fold precision rounded to the working precision) */
double sumk(double* p,const unsigned int n,const unsigned int k){
	for(unsigned int i=1; i<k; ++i){
		vec_sum(p,n);
	}
	double res = p[0];
	for(unsigned int i=1; i<n; ++i){
		res += p[i];
	}
	return res;
}
/* Sumk Cmplx (complex summation in k-fold precision) */
double complex sumk_cmplx(double complex* p,const unsigned int n,const unsigned int k){
	struct eft* eft_arr = (struct eft*)malloc(2*sizeof(struct eft));
	for(unsigned int i=1; i<k; ++i){
		vec_sum_cmplx(p,eft_arr,n);
	}
	free(eft_arr);
	double complex res = p[0];
	for(unsigned int i=1; i<n; ++i){
		res += p[i];
	}
	return res;
}
/* Sumkk (summation in k-fold precision stored in k-parts) */
double* sumkk(double* p,const unsigned int n,const unsigned int k){
	double* res = (double*)malloc(k*sizeof(double));
	for(unsigned int i=0; i<k-1; ++i){
		vec_sum(p,n-i);
		res[i] = p[n-i-1];
	}
	res[k-1] = p[0];
	for(unsigned int i=1; i<n-k+1; ++i){
		res[k-1] += p[i];
	}
	return res;
}
/* Sumkk Cmplx (complex summation in k-fold precision stored in k-parts) */
double complex* sumkk_cmplx(double complex* p,const unsigned int n,const unsigned int k){
	double complex* res = (double complex*)malloc(k*sizeof(double complex));
	struct eft* eft_arr = (struct eft*)malloc(2*sizeof(struct eft));
	for(unsigned int i=0; i<k-1; ++i){
		vec_sum_cmplx(p,eft_arr,n-i);
		res[i] = p[n-i-1];
	}
	free(eft_arr);
	res[k-1] = p[0];
	for(unsigned int i=1; i<n-k+1; ++i){
		res[k-1] += p[i];
	}
	return res;
}
/* Sum MPFR (sum in multi-precision using MPFR) */
double sumMPFR(double* p,const unsigned int n,const unsigned int prec){
	mpfr_t res_mp;
	mpfr_init2(res_mp,prec);
	mpfr_set_d(res_mp,0,MPFR_RNDN);
	for(unsigned int i=0; i<n; ++i){
		mpfr_add_d(res_mp,res_mp,p[i],MPFR_RNDN);
	}
	return mpfr_get_d(res_mp,MPFR_RNDN);
}
/* Sort */
void sort(double *p)
{
	double max, temp;
	int ind, i, j;
	for(i=0; i<3; ++i){
		max = fabs(p[i]);
		ind = i;
		for(j=i+1; j<4; ++j){
			temp = fabs(p[j]);
			if(temp > max){
				max = temp;
				ind = j;
			}
		}
		if(ind != i){
			temp = p[i];
			p[i] = p[ind];
			p[ind] = temp;
		}
	}
}
/* Priest Summation */
double priest_sum(double *p){
	sort(p);
	double s = p[0], c = 0, y, u, t, v, z;
	for(int i=1; i<4; ++i){
		y = c + p[i];
		u = p[i] - (y - c);
		t = y + s;
		v = y - (t - s);
		z = u + v;
		s = t + z;
		c = z - (s - t);
	}
	return s;
}
/* Priest Summation Cmplx */
double complex priest_sum_cmplx(const double complex *p){
	double realp[4] = {creal(p[0]),creal(p[1]),creal(p[2]),creal(p[3])};
	double imagp[4] = {cimag(p[0]),cimag(p[1]),cimag(p[2]),cimag(p[3])};
	return priest_sum(realp)+I*priest_sum(imagp);
}