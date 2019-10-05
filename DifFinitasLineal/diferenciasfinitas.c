#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void diferenciasFinitas(pEDOF prob,double xa,double xb, FILE* archivo)
{
  // Se completan variables del problema
  double h=(xb-xa)/(prob.n+1);
  prob.h=h;
  gsl_vector *sol_y=gsl_vector_alloc(prob.n);
  prob.w=sol_y;
  
  double x[prob.n];
  x[0]=h+xa;
  for (int i = 0; i < prob.n-1; ++i) {
    x[i+1]=x[i]+h;
  }
  prob.x=x;

  // Se resuelve el metodo
  metodoDifFin(prob);
  
  // Se imprime en el archivo
  fprintf(archivo,"%E %E\n", xa, prob.cond[0]);
  for (int i = 0; i < prob.n; ++i) {
    fprintf(archivo,"%E %E\n",prob.x[i],gsl_vector_get(prob.w,i));
  }
  fprintf(archivo,"%E %E\n", xb, prob.cond[1]);
  gsl_vector_free(sol_y);
}

void metodoDifFin(pEDOF prob)
{
  /* Este metodo genera el sistema tridiagonal
   * que se resuelve con el solver tridiagonal gsl
   * creando tres vectores de la matriz A que son
   * d, diagonal
   * e, diagonal superior
   * f, diagonal inferior
   * Ademas de crear un vector b
   * Y usar el vector solucion prob.w
   */
 
  int n=prob.n;
  double h=prob.h;
  gsl_vector *d=gsl_vector_alloc(n);
  gsl_vector *e=gsl_vector_alloc(n-1);
  gsl_vector *f=gsl_vector_alloc(n-1);
  gsl_vector *b=gsl_vector_alloc(n);
  double aux=0;

  // Se definen los valores de los vectores d,e,f y b
  for (int i = 0; i < n; ++i) {
    gsl_vector_set(d,i,2+h*h*prob.q(prob.x[i]));
  }
  for (int i = 0; i < n-1; ++i) {
    gsl_vector_set(e,i,-1+h/2*prob.p(prob.x[i]));
    gsl_vector_set(f,i,-1-h/2*prob.p(prob.x[i+1]));
  }
  aux=-h*h*prob.r(prob.x[0])+(1+h/2*prob.p(prob.x[0]))*prob.cond[0];
  gsl_vector_set(b,0,aux);
  aux=-h*h*prob.r(prob.x[n-1])+(1-h/2*prob.p(prob.x[n-1]))*prob.cond[1];
  gsl_vector_set(b,n-1,aux);
  for (int i = 1; i < n-1; ++i) {
    gsl_vector_set(b,i,-h*h*prob.r(prob.x[i]));
  }

  // Se resuelve el sistema
  gsl_linalg_solve_tridiag(d,e,f,b,prob.w);

  gsl_vector_free(d);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
}
