#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct EcuacionDiferencialParcialParabolica{
  /*Parametros del problema (no relacionados a la EDP, sino del contexto)*/
  double *parametros;
  double D;
  double** m_inicial;//condicion inicial

  double xi;
  double xf;
  double yi;
  double yf;

  double (*fuente)(double x,double y ,double* parametros);

  /* condiciones de borde */
  //coeficientes
  double alpha[4];
  double beta[4];

  //funciones
  double (*f1)(double y ,double* parametros);
  double (*f2)(double y ,double* parametros);
  double (*f3)(double x ,double* parametros);
  double (*f4)(double x ,double* parametros);

  double (*g1)(double y ,double* parametros);
  double (*g2)(double y ,double* parametros);
  double (*g3)(double x ,double* parametros);
  double (*g4)(double x ,double* parametros);

  /*Parametros del método*/
  int nx;
  int ny;
  int nt;
  double tol;

}ecp;

void resolver_metodo(ecp);
