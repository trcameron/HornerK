#ifndef SUM
#define SUM
#include "eft.h"
/* Sumk (summation in k-fold precision) */
double sumk(double* p,const unsigned int n,const unsigned int k);
/* Sumk Cmplx (complex summation in k-fold precision) */
double complex sumk_cmplx(double complex* p,const unsigned int n,const unsigned int k);
/* Sumkk (summation in k-fold precision stored in k-parts) */
double* sumkk(double* p,const unsigned int n,const unsigned int k);
/* Sumkk Cmplx (complex summation in k-fold precision stored in k-parts) */
double complex* sumkk_cmplx(double complex* p,const unsigned int n,const unsigned int k);
/* Sum MPFR (sum in multi-precision using MPFR) */
double sumMPFR(double* p,const unsigned int n,const unsigned int prec);
/* Sort */
void sort(double *p);
/* Priest Summation */
double priest_sum(double *p);
/* Priest Summation Cmplx */
double complex priest_sum_cmplx(const double complex *p);
#endif