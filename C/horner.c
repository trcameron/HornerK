#include "horner.h"
#include <mpc.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include "sum.h"
/* Horner (Horner's method in double precision) */
double horner(const double* poly,const double x,const unsigned int deg){
	double res = poly[deg];
	for(int i=deg-1; i>=0; --i){
		res = poly[i] + x*res;
	}
	return res;
}
/* Horner Cmplx (Horner's method in complex double precision) */
double complex horner_cmplx(const double complex* poly,const double complex x,const unsigned int deg){
	double complex res = poly[deg];
	for(int i=deg-1; i>=0; --i){
		res = poly[i] + x*res;
	}
	return res;
}
/* HornerK Step (Horner's method step in k-fold precision and stored in k-parts) */
void hornerK_step(const double a,const double b,double* x,double* e,const unsigned int k){
	const unsigned int n = 2*k;
	double p,h;
	struct eft res;
	two_prod(b,x[0],&res);
	p = res.fl_res; e[0] = res.fl_err;
	for(unsigned int i=1; i<k; ++i){
		two_prod(b,x[i],&res);
		h = res.fl_res; e[i] = res.fl_err;
		two_sum(p,h,&res);
		p = res.fl_res; e[k+i-1] = res.fl_err;
	}
	two_sum(p,a,&res);
	x[0] = res.fl_res; e[n-1] = res.fl_err;
	for(unsigned int i=0; i<k-2; ++i){
		vec_sum(e,n-i);
		x[i+1] = e[n-i-1];
	}
	x[k-1] = e[0];
	for(unsigned int i=1; i<k+2; ++i){
		x[k-1] += e[i];
	}
}
/* HornerK Step Cmplx (Complex Horner's method step in k-fold precision and stored in k-parts) */
void hornerK_step_cmplx(const double complex a,const double complex b,double complex* x,double complex* e,struct eft* eft_arr,const unsigned int k){
	double complex p,h;
	struct eft_cmplx_sum sum;
	struct eft_cmplx_prod prod;
	two_prod_cmplx(b,x[0],eft_arr,&prod);
	p = prod.fl_res; e[0] = prod.fl_err1; e[1] = prod.fl_err2; e[2] = prod.fl_err3;
	for(unsigned int i=1; i<k; ++i){
		two_prod_cmplx(b,x[i],eft_arr,&prod);
		h = prod.fl_res; e[3*i] = prod.fl_err1; e[3*i+1] = prod.fl_err2; e[3*i+2] = prod.fl_err3;
		two_sum_cmplx(p,h,eft_arr,&sum);
		p = sum.fl_res; e[3*k+i-1] = sum.fl_err;
	}
	two_sum_cmplx(p,a,eft_arr,&sum);
	x[0] = sum.fl_res; e[4*k-1] = sum.fl_err;
	for(unsigned int i=0; i<k-2; ++i){
		vec_sum_cmplx(e,eft_arr,4*k-i);
		x[i+1] = e[4*k-i-1];
	}
	x[k-1] = e[0];
	for(unsigned int i=1; i<3*k+2; ++i){
		x[k-1] += e[i];
	}
}
/* Hornerk (Horner's method in k-fold precision and rounded into the working precision) */
double hornerK(const double* poly,const double x,const unsigned int deg,const unsigned int k){
	double* coll = (double*)malloc(k*sizeof(double));
	double* e = (double*)malloc(2*k*sizeof(double));
	coll[0] = poly[deg];
	memset(coll+1,0,(k-1)*sizeof(double));
	for(int i=deg-1; i>=0; --i){
		hornerK_step(poly[i],x,coll,e,k);
	}
	double res = sumk(coll,k,k);
	free(coll);
	free(e);
	return res;
}
/* Hornerk Cmplx (Complex Horner's method in k-fold precision and rounded into the working precision) */
double complex hornerK_cmplx(const double complex* poly,const double complex x,const unsigned int deg,const unsigned int k){
	double complex* coll = (double complex*)malloc(k*sizeof(double complex));
	double complex* e = (double complex*)malloc(4*k*sizeof(double complex));
	struct eft* eft_arr = (struct eft*)malloc(6*sizeof(struct eft));
	coll[0] = poly[deg];
	memset(coll+1,0,(k-1)*sizeof(double complex));
	for(int i=deg-1; i>=0; --i){
		hornerK_step_cmplx(poly[i],x,coll,e,eft_arr,k);
	}
	double complex res = sumk_cmplx(coll,k,k);
	free(coll);
	free(e);
	free(eft_arr);
	return res;
}
/* HornerMPFR (Horner's method in multi-precision using MPFR) */
double hornerMPFR(const double* poly,const double x,const unsigned int deg,const unsigned int prec){
	mpfr_t res_mp;
	mpfr_init2(res_mp,prec);
	mpfr_set_d(res_mp,poly[deg],MPFR_RNDN);

	for(int i=deg-1; i>=0; --i){
		mpfr_mul_d(res_mp,res_mp,x,MPFR_RNDN);
		mpfr_add_d(res_mp,res_mp,poly[i],MPFR_RNDN);
	}
	return mpfr_get_d(res_mp,MPFR_RNDN);
}
/* HornerMPFR Cmplx (Complex Horner's method in multi-precision using MPFR) */
double complex hornerMPFR_cmplx(const double complex* poly,const double complex x,const unsigned int deg,const unsigned int prec){
	mpc_t* poly_mp = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
	mpc_t x_mp, res_mp;
	mpc_init2(x_mp,prec); mpc_init2(res_mp,prec);
	mpc_set_dc(x_mp,x,MPC_RNDNN); mpc_set_dc(res_mp,poly[deg],MPC_RNDNN);
	for(unsigned int i=0; i<=deg; ++i){
		mpc_init2(poly_mp[i],prec);
		mpc_set_dc(poly_mp[i],poly[i],MPC_RNDNN);
	}
	for(int i=deg-1; i>=0; --i){
		mpc_mul(res_mp,res_mp,x_mp,MPC_RNDNN);
		mpc_add(res_mp,poly_mp[i],res_mp,MPC_RNDNN);
	}
	double complex res = mpc_get_dc(res_mp,MPC_RNDNN);
	mpc_clear(x_mp); mpc_clear(res_mp);
	for(unsigned int i=0; i<=deg; ++i){
		mpc_clear(poly_mp[i]);
	}
	free(poly_mp);
	return res;
}
/* HornerComp Cmplx (compensated Horner's method in double precision and rounded into working preccision) */
double complex hornerComp_cmplx(const double complex* poly,const double complex x,const unsigned int deg){
	struct eft eft_arr[6];
	struct eft_cmplx_sum tsc;
	struct eft_cmplx_prod tpc;
	double complex e = 0;
	double complex p[4];
	double complex res = poly[deg];
	for(int i=deg-1; i>=0; --i){
		two_prod_cmplx(res,x,eft_arr,&tpc);
		two_sum_cmplx(tpc.fl_res,poly[i],eft_arr,&tsc);
		res = tsc.fl_res;
		p[0] = tpc.fl_err1; p[1] = tpc.fl_err2; p[2] = tpc.fl_err3; p[3] = tsc.fl_err;
		e = e*x + priest_sum_cmplx(p);
	}
	res += e;
	return res;
}