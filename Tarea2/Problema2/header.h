#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct EcHiperbolica
{
 /*Parametros del problema (no relacionados a la EDP, sino del contexto)*/
  double *parametros;
  double D;

  double xi;
  double xf;

  // DEFINIR TIEMPOS INICIALES
  double ti;
  double tf;

  /* condiciones de borde */
  //coeficientes
  double alpha[2];
  double beta[2];
  double gamma[2];

  //funciones de condiciones de borde
  double (*alpha1)(double t,double* parametros);
  double (*alpha2)(double t,double* parametros);

  double (*beta1)(double t,double* parametros);
  double (*beta2)(double t,double* parametros);

  double (*gamma1)(double t,double* parametros);
  double (*gamma2)(double t,double* parametros);
  
  // funciones de condiciones iniciales
  double (*f)(double x ,double* parametros);
  double (*g)(double x ,double* parametros);
  
  /*Parametros del método*/
  int nx;
  int nt;
  double lambda;

}hiper;

void resolver_metodo(hiper, FILE*);
void rellenar_A(double **, hiper);
