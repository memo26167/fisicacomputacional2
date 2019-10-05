#ifndef BIB_H
#define BIB_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef double (*FC)(double,double,double) ;
typedef void (*RKF)(double*,double,double*,FC ,FC ,FC );
typedef struct problemaDiferenciasFinitasNoLineal {
  FC fn;//funcion de la edo
  FC fy;//derivada parcial con respecto a y de la funcion de la edo
  FC fyp;//derivada parcial con respecto a y' de la funcion de la edo
  double *cond;// condiciones de frontera, primero w0 y luego wn+1
  int n;//cantidad de particiones
  double h;//paso
  gsl_vector *x;// vector de variable independiente
  gsl_vector *w;// vector solucion a iterar
} pDFNL;

double diferenciasStep(pDFNL);
void difFinitasNoLineales(pDFNL,double, double,double, FILE*);

#endif
