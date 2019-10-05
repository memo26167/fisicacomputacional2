#include <stdio.h>
#include <stdlib.h>
#include "header.h"


void ecuacion1(double *k,double x, double *y, FC p, FC q, FC r)
{
  k[0]=y[1];
  k[1]=p(x)*y[1]+q(x)*y[0]+r(x);
}

void ecuacion2(double *k,double x, double *y,FC p, FC q,FC r)
{
  k[0]=y[1];
  k[1]=p(x)*y[1]+q(x)*y[0];
  x=r(x); // para sacar warnings
}

/* Realiza el método de Shooting Lineal */
void shooting(double * cond_fronteras, double paso_h, pRKF problema,gsl_vector *sol_x,gsl_vector *solucion_frontera, int num_it, FILE* fp)
{
  /* Se definen los problemas de valor inicial
   * Se asignan las funciones de
   * -EDO de problema de valor inicial
   * -Funcion P del problema de frontera
   * -Funcion Q del problema de frontera
   * -Funcion R del problema de frontera
   * -Condiciones iniciales del problema de valor inicial
   * -Matriz solucion del problema de valor inicial
   */

  // Problema de valor inicial 1
  pRKF edo1;
  edo1.f=&ecuacion1;
  edo1.p=problema.p;
  edo1.q=problema.q;
  edo1.r=problema.r;
  double edo1y0[2]={cond_fronteras[0],0};
  gsl_matrix *sol_y1=gsl_matrix_alloc(num_it,2);
  // Se resuelve con el metodo RKF
  rk(edo1y0,2,num_it,paso_h,edo1,sol_y1,sol_x);

  // Problema de valor inicial 2
  pRKF edo2;
  edo2.f=&ecuacion2;
  edo2.p=problema.p;
  edo2.q=problema.q;
  edo2.r=problema.r;
  double edo2y0[2]={0,1};
  gsl_matrix *sol_y2=gsl_matrix_alloc(num_it,2);
  // Se resuelve con el metodo RKF
  rk(edo2y0,2,num_it,paso_h,edo2,sol_y2,sol_x);

  // Se calcula la solucion del problema de frontera
  double escalar=(cond_fronteras[1]-gsl_matrix_get(sol_y1,num_it-1,0))/gsl_matrix_get(sol_y2,num_it-1,0);
  double aux=0;
  for (int i = 0; i < num_it; ++i) {
    aux=gsl_matrix_get(sol_y1,i,0)+escalar*gsl_matrix_get(sol_y2,i,0);
    gsl_vector_set(solucion_frontera,i,aux);
  }
  
  // Se libera 
  gsl_matrix_free(sol_y1);
  gsl_matrix_free(sol_y2);
}
