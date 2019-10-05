#ifndef BIB_H
#define BIB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void condIniciales(gsl_matrix *, double *, int);
void fEqDiferencial(gsl_vector *,double, gsl_vector *,int , void (*)(double*,double,double*));
void metodoRKF(gsl_matrix *, gsl_vector *, int, int,double,double,void (*)(double *, double, double *));
void rk(double *,int , int, double, double,void (*)(double *, double, double *), FILE *);
#endif
