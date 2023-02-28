#include "eft.h"
#include "math.h"
#include <stdlib.h>
/* Two Sum (error free transformation of sum operation )*/
void two_sum(const double a,const double b,struct eft* res){
	res->fl_res = a + b;
	double t = res->fl_res - a;
	res->fl_err = (a - (res->fl_res - t)) + (b - t);
}
/* Two Sum Cmplx (error free transformation of complex sum operation) */
void two_sum_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_sum* res){
	two_sum(creal(a),creal(b),eft_arr);
	two_sum(cimag(a),cimag(b),eft_arr+1);
	res->fl_res = eft_arr[0].fl_res + I*eft_arr[1].fl_res;
	res->fl_err = eft_arr[0].fl_err + I*eft_arr[1].fl_err;
}
/* Two Product (error free transformation of product operation)*/
void two_prod(const double a,const double b,struct eft* res){
	res->fl_res = a*b;
	res->fl_err = fma(a,b,-res->fl_res);
}
/* Two Product Cmplx (error free transformation of complex product operation)*/
void two_prod_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_prod* res){
	two_prod(creal(a),creal(b),eft_arr);
	two_prod(cimag(a),cimag(b),eft_arr+1);
	two_prod(creal(a),cimag(b),eft_arr+2);
	two_prod(cimag(a),creal(b),eft_arr+3);
	two_sum(eft_arr[0].fl_res,-eft_arr[1].fl_res,eft_arr+4);
	two_sum(eft_arr[2].fl_res,eft_arr[3].fl_res,eft_arr+5);
	res->fl_res = eft_arr[4].fl_res + I*eft_arr[5].fl_res;
	res->fl_err1 = eft_arr[0].fl_err + I*eft_arr[2].fl_err;
	res->fl_err2 = -eft_arr[1].fl_err + I*eft_arr[3].fl_err;
	res->fl_err3 = eft_arr[4].fl_err + I*eft_arr[5].fl_err;
}
/* Vector Sum (distillation algorithm) */
void vec_sum(double* p,const unsigned int n){
	struct eft res;
	for(unsigned int i=1; i<n; ++i){
		two_sum(p[i],p[i-1],&res);
		p[i] = res.fl_res;
		p[i-1] = res.fl_err;
	}
}
/* Vector Sum Cmplx (distillation algorithm in complex arithmetic) */
void vec_sum_cmplx(double complex* p,struct eft* eft_arr,const unsigned int n){
	struct eft_cmplx_sum res;
	for(unsigned int i=1; i<n; ++i){
		two_sum_cmplx(p[i],p[i-1],eft_arr,&res);
		p[i] = res.fl_res;
		p[i-1] = res.fl_err;
	}
}