#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void diferenciasFinitas(pEl prob ,FILE *archivo)
{
  double n=prob.n;
  double m=prob.m;
  double deltax=(prob.xb-prob.xa)/(n-1);
  double deltay=(prob.yb-prob.ya)/(m-1);

  /* Inicializar variables */
  gsl_vector *sol_phi=gsl_vector_alloc(n*m);
  gsl_vector_set_zero(sol_phi);
  prob.phi=sol_phi;

  gsl_vector *x=gsl_vector_alloc(n);
  gsl_vector *y=gsl_vector_alloc(m);
  gsl_vector_set_zero(x);
  gsl_vector_set_zero(y);
  
  for (int i = 0; i < n; ++i) {
    gsl_vector_set(x,i,deltax*i+prob.xa);
  }
  for (int j = 0; j < m; ++j) {
    gsl_vector_set(y,j,deltay*j+prob.ya);
  }
  prob.x=x;
  prob.y=y;
  
  gsl_matrix *sys_eq=gsl_matrix_alloc(n*m,n*m);//sistema de ecuaciones
  gsl_matrix_set_zero(sys_eq);//sistema de ecuaciones en cero
  gsl_vector *sys_b=gsl_vector_alloc(n*m);
  gsl_vector_set_zero(sys_b);
  double gamma[4];
  for (int i = 0; i < 4; ++i) {
    gamma[i]=0;
    if(prob.alpha[i]!=0){
      gamma[i]=prob.beta[i]/prob.alpha[i];;
    }
  }
  int k=0;

  /* Rellenar Sistema de ecuaciones de diferencias finitas */
  double aux;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      k=i+j*n;
      if(i==0){

	gsl_matrix_set(sys_eq,k,k,1);
	aux=gamma[1]*prob.f1(gsl_vector_get(y,j));
	gsl_vector_set(sys_b,k,aux);
      }
      else if(i==n-1){

	gsl_matrix_set(sys_eq,k,k,1);
	aux=gamma[3]*prob.f3(gsl_vector_get(y,j)); 
	gsl_vector_set(sys_b,k,aux);
      }
      else if(j==0){

	gsl_matrix_set(sys_eq,k,k,1);
	aux=gamma[0]*prob.f0(gsl_vector_get(x,i));
	gsl_vector_set(sys_b,k,aux);
      }
      else if(j==m-1){
	gsl_matrix_set(sys_eq,k,k,1);
	aux=gamma[2]*prob.f2(gsl_vector_get(x,i));
	gsl_vector_set(sys_b,k,aux);
      }
      else{//puntos fuera de los bordes
	gsl_matrix_set(sys_eq,k,k,2*(pow(deltax,-2)+pow(deltay,-2)));//i,j	
	gsl_matrix_set(sys_eq,k,k-1,-pow(deltax,-2));//i-1,j
	gsl_matrix_set(sys_eq,k,k+1,-pow(deltax,-2));//i+1,j
	gsl_matrix_set(sys_eq,k,k-n,-pow(deltay,-2));//i,j-1
	gsl_matrix_set(sys_eq,k,k+n,-pow(deltay,-2));//i,j+1
	gsl_vector_set(sys_b,k,prob.S(gsl_vector_get(x,i),gsl_vector_get(y,j)));
      }
    }
  }//fori

  /* for (int i = 0; i < n*m; ++i) { */
  /*   for (int j = 0; j < n*m; ++j) { */
  /*     printf("%f ", gsl_matrix_get(sys_eq,i,j)); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* for (int i = 0; i < n; ++i) { */
  /*   printf("%f \n",gsl_vector_get(sys_b,i)); */
  /* } */
  
  gsl_permutation * p = gsl_permutation_alloc (n*m);
  int signum;
  gsl_linalg_LU_decomp (sys_eq, p, &signum);
  gsl_linalg_LU_solve (sys_eq, p, sys_b, prob.phi);
  gsl_permutation_free (p);
  gsl_matrix_free(sys_eq);
  
  // Imprimir archivo
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      fprintf(archivo, "%E ", gsl_vector_get(x,i));
      fprintf(archivo, "%E ", gsl_vector_get(y,j));
      fprintf(archivo, "%E\n",gsl_vector_get(sol_phi,j*n+i));
    }
    fprintf(archivo, "\n");
  }
  
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(sys_b);
}
