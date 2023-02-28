#ifndef HORNER
#define HORNER
#include <complex.h>
#include "eft.h"
/* Horner (Horner's method in double precision) */
double horner(const double* poly,const double x,const unsigned int deg);
/* Horner Cmplx (Horner's method in complex double precision) */
double complex horner_cmplx(const double complex* poly,const double complex x,const unsigned int deg);
/* HornerK Step (Horner's method step in k-fold precision and stored in k-parts) */
void hornerK_step(const double a,const double b,double* x,double* e,const unsigned int k);
/* HornerK Step Cmplx (Complex Horner's method step in k-fold precision and stored in k-parts) */
void hornerK_step_cmplx(const double complex a,const double complex b,double complex* x,double complex* e,struct eft* eft_arr,const unsigned int k);
/* Hornerk (Horner's method in k-fold precision and rounded into the working precision) */
double hornerK(const double* poly,const double x,const unsigned int deg,const unsigned int k);
/* Hornerk Cmplx (Complex Horner's method in k-fold precision and rounded into the working precision) */
double complex hornerK_cmplx(const double complex* poly,const double complex x,const unsigned int deg,const unsigned int k);
/* HornerMPFR (Horner's method in multi-precision using MPFR) */
double hornerMPFR(const double* poly,const double x,const unsigned int deg,const unsigned int prec);
/* HornerMPFR Cmplx (Complex Horner's method in multi-precision using MPFR) */
double complex hornerMPFR_cmplx(const double complex* poly,const double complex x,const unsigned int deg,const unsigned int prec);
/* HornerComp Cmplx (compensated Horner's method in double precision and rounded into working preccision) */
double complex hornerComp_cmplx(const double complex* poly,const double complex x,const unsigned int deg);
#endif