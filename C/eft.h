#ifndef EFT
#define EFT
#include <complex.h>
/* EFT Data Structure */
struct eft{
	double fl_res, fl_err;
};
/* EFT Cmplx Data Structure for Sum */
struct eft_cmplx_sum{
	double complex fl_res, fl_err;
};
/* EFT Cmplx Data Structure for Product */
struct eft_cmplx_prod{
	double complex fl_res, fl_err1, fl_err2, fl_err3;
};
/* Two Sum (error free transformation of sum operation) */
void two_sum(const double a,const double b,struct eft* res);
/* Two Sum Cmplx (error free transformation of complex sum operation) */
void two_sum_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_sum* res);
/* Two Product (error free transformation of product operation)*/
void two_prod(const double a,const double b,struct eft* res);
/* Two Product Cmplx (error free transformation of complex product operation)*/
void two_prod_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_prod* res);
/* Vector Sum (distillation algorithm) */
void vec_sum(double* p,const unsigned int n);
/* Vector Sum Cmplx (distillation algorithm in complex arithmetic) */
void vec_sum_cmplx(double complex* p,struct eft* eft_arr,const unsigned int n);
#endif