#ifndef BIB_H
#define BIB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


void diferencial(double *, double*,double ,int , double);
double potencial(double);
void numerov(int,double,double,double,int,double,double,double,FILE*);

#endif
