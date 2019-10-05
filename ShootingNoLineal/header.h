#ifndef BIB_H
#define BIB_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef double (*FC)(double,double,double) ;
typedef void (*RKF)(double*,double,double*,FC ,FC ,FC );
typedef struct problemaRKF {
  FC fn;//funcion de la edo
  FC fy;//derivada parcial con respecto a y de la funcion de la edo
  FC fyp;//derivada parcial con respecto a y' de la funcion de la edo
  RKF f;//funcion del sistema de edo's
} pRKF;
void fEqDiferencial(gsl_vector *,double, gsl_vector*,int, pRKF);
void metodoRKF(gsl_matrix *, gsl_vector *, int, int,double, pRKF);
void condIniciales(gsl_matrix *, double *, int);
void rk(double*,int, int , double, pRKF, gsl_matrix *,gsl_vector *);
void shooting(double, double, double *,double, pRKF , FILE*);
#endif
