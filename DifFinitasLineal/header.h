#ifndef BIB_H
#define BIB_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
typedef double (*FC)(double) ;
typedef struct problemaEDOFinitas{
  FC p;
  FC q;
  FC r;
  double* cond;
  int n;
  double xa;
  double xb;
  gsl_vector *w;
  double *x;
  double h;
}pEDOF;

void metodoDifFin(pEDOF);
void diferenciasFinitas(pEDOF,double,double, FILE*);
#endif
