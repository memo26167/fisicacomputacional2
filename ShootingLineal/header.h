#ifndef BIB_H
#define BIB_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
typedef double (*FC)(double) ;
typedef void (*RKF)(double*,double,double*,FC ,FC ,FC );
typedef struct problemaRKF {
  RKF f;
  FC p;
  FC q;
  FC r;
} pRKF;

void rk(double *,int,int,double, pRKF, gsl_matrix *,gsl_vector *);
void condIniciales(gsl_matrix *, double *, int);
void metodoRKF(gsl_matrix *, gsl_vector *, int, int,double ,pRKF);
void fEqDiferencial(gsl_vector *,double , gsl_vector *,int , pRKF);
void ecuacion1(double *,double,double *,FC ,FC,FC);
void ecuacion2(double *,double,double *,FC ,FC,FC);
void shooting(double , double , double *, double, pRKF, FILE* );
#endif
