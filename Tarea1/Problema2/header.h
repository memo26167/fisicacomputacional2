#ifndef BIB_H
#define BIB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef double (*FC)(double);//funcion común
typedef double (*F2V)(double,double);//funcion dos variables
typedef struct problemaEliptico
{
  int n;//cantidad de puntos en la direccion x
  int m;//cantidad de puntos en la direccion y
  //vector con sol, que simula una matriz nXm cuya primera componente es x segunda y
  gsl_vector *phi;
  gsl_vector *x;
  gsl_vector *y;
  double xa;
  double xb;
  double ya;
  double yb;
  
  F2V S;  //funcion fuente de la ecuacion diferencial

  FC f0; //funcion del borde inferior de y, y=0,  que depende de x
  FC f1; //funcion del borde inferior de x, x=0,  que depende de y
  FC f2; //funcion del borde superior de y, y=Ly, que depende de x
  FC f3; //funcion del borde superior de x, x=Lx, que depende de y

  //vectores de condiciones de borde
  /* 0 Primer componente,  borde inferior de y, y=0
   * 1 Segundo componente, borde inferior de x, x=0
   * 2 Tercer componente,  borde superior de y, y=Ly
   * 3 Segundo componente, borde superior de x, x=Lx
   */
  double * alpha;
  double * beta;
  double * gamma;
}pEl;

void diferenciasFinitas(pEl,FILE*);
#endif
