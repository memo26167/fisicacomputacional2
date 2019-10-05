#include <stdio.h>
#include <stdlib.h>
#include "header.h"
/*Se modifica el problema,
  Antes cond[1]=wn+1
  Ahora cond[1]=wn'
 */


double diferenciasStep(pDFNL prob)
{
  /* Se escribe el sistema tridiagona, Jv=b, con J=[f,d,e] donde
   * d es el vector diagonal
   * e es el vector superior
   * f es el vector inferior
   * b es el vector a la derecha de el sistema  resolver
   * v es el vector solucion
   */
  
  //salida de la funcion
  double error=0;
  //variables del problema
  double h=prob.h;
  int n=prob.n;
  gsl_vector *jacob_d=gsl_vector_alloc(n);
  gsl_vector *jacob_e=gsl_vector_alloc(n-1);
  gsl_vector *jacob_f=gsl_vector_alloc(n-1);
  gsl_vector *vector_b=gsl_vector_alloc(prob.n);
  gsl_vector *vector_v=gsl_vector_alloc(prob.n);
  gsl_vector_set_zero(jacob_d);
  gsl_vector_set_zero(jacob_e);
  gsl_vector_set_zero(jacob_f);
  gsl_vector_set_zero(vector_b);
  gsl_vector_set_zero(vector_v);
  
  double aux=0;
  double xi=0; //variable que almacena el valor de x en un componente i
  double wi=0; //variable que almacena el valor de w en un componente i
  double wpi=0;//variable que almacena el valor de w' en un componente i

  // Primero definimos extremos dependientes de w0 y wn+1
  //i=1 (enrealidad i=0, pero se usa el i del libro burden) (requiere w0)
  xi=gsl_vector_get(prob.x,0);
  wi=gsl_vector_get(prob.w,0);
  wpi=( gsl_vector_get(prob.w,1)-prob.cond[0])/(2*h);//(w2-w0)/(2h)
  // * vector d
  aux=2+h*h*prob.fy(xi,wi,wpi);
  gsl_vector_set(jacob_d, 0,aux);
  // * vector e
  aux=-1+h/2*prob.fyp(xi,wi,wpi);
  gsl_vector_set(jacob_e, 0,aux);
  // * vector b
  aux= -(-gsl_vector_get(prob.w,1)+2*wi-prob.cond[0]+h*h*prob.fn(xi,wi,wpi));
  gsl_vector_set(vector_b,0,aux);
    
  //i=n (enrealidad i=n-1)(requiere wn+1)
  xi=gsl_vector_get(prob.x,n-1);
  wi=gsl_vector_get(prob.w,n-1);
  wpi=prob.cond[1];//wn'=(wn+1-wn-1)/(2h), entonces wn+1=2*h*wn'+wn-1
  // * vector d
  aux=2+h*h*prob.fy(xi,wi,wpi);
  gsl_vector_set(jacob_d, n-1,aux);
  // * vector f
  aux=-1-h/2*prob.fyp(xi,wi,wpi);
  gsl_vector_set(jacob_f, n-2,aux);
  // * vector b
  aux= -(-2*h*prob.cond[1]+2*wi-2*gsl_vector_get(prob.w,n-2)+h*h*prob.fn(xi,wi,wpi));
  gsl_vector_set(vector_b,n-1,aux);
  
  for (int i = 1; i <= n-2; ++i) {
    //parte desde 1 porq ya se calculo el comp. 0 de vector d
    //y termina en n-2 porq ya se calculo el comp n-1 del vector d
    xi=gsl_vector_get(prob.x,i);
    wi=gsl_vector_get(prob.w,i);
    wpi=(gsl_vector_get(prob.w,i+1)-gsl_vector_get(prob.w,i-1))/(2*h);
    //vector d
    aux=2+h*h*prob.fy(xi,wi,wpi);
    gsl_vector_set(jacob_d,i,aux);
    //vector e
    aux=-1+h/2*prob.fyp(xi,wi,wpi);
    gsl_vector_set(jacob_e,i,aux);
    //vector f
    aux=-1-h/2*prob.fyp(xi,wi,wpi);
    gsl_vector_set(jacob_f,i-1,aux);
    //vector b
    aux= -(-gsl_vector_get(prob.w,i+1)+2*wi-gsl_vector_get(prob.w,i-1)+h*h*prob.fn(xi,wi,wpi));
    gsl_vector_set(vector_b,i,aux);
  }

  gsl_linalg_solve_tridiag(jacob_d,jacob_e,jacob_f,vector_b,vector_v);
  gsl_vector_add(prob.w,vector_v);// w(k+1)=w(k)+v

  
  error=gsl_blas_dnrm2(vector_v);
  
  gsl_vector_free(jacob_d);
  gsl_vector_free(jacob_e);
  gsl_vector_free(jacob_f);
  gsl_vector_free(vector_b);
  gsl_vector_free(vector_v);
  
  return error;
}

void difFinitasNoLineales(pDFNL prob, double xa, double xb,double tol, FILE* archivo)
{
  gsl_rng * r= gsl_rng_alloc (gsl_rng_taus);
  prob.h=(xb-xa)/(prob.n);
  double aux=0;
  // definir vector x e w
  gsl_vector *x=gsl_vector_alloc(prob.n);
  gsl_vector *w=gsl_vector_alloc(prob.n);
  prob.x=x;
  prob.w=w;
  for (int i = 0; i < prob.n; ++i) {
    //inicializar x
    gsl_vector_set(prob.x,i,xa+(i+1)*prob.h);
    aux=gsl_ran_gaussian(r,1);
    //inicializar w
    gsl_vector_set(prob.w,i,aux);
  }
  
  double error=1000;
  while(error>tol){
    error=diferenciasStep(prob);
  }
  
  gsl_rng_free (r);
  
  // Imprimir archivo
  fprintf(archivo, "%E %E \n", xa,prob.cond[0]);
  for (int i = 0; i < prob.n; ++i) {
    fprintf(archivo, "%E %E \n",gsl_vector_get(prob.x,i),gsl_vector_get(prob.w,i));
  }
}

