#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct EcuacionDiferencialParcialParabolica{
  int nx;
  int ny;
  int nt;
  double tol;
  double D;
  double** m_inicial;//condicion inicial
  double cond_alpha[4];
  double cond_beta[4];

  double (*f1)(double x ,double* parametros);
  double (*f2)(double y ,double* parametros);
  double (*f3)(double x ,double* parametros);
  double (*f4)(double y ,double* parametros);

  double (*g1)(double x ,double* parametros);
  double (*g2)(double y ,double* parametros);
  double (*g3)(double x ,double* parametros);
  double (*g4)(double y ,double* parametros);
}ecp;
