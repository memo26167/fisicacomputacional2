#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void resuelveSistema(double ** A, double * b, double * x, int n)
{
  gsl_matrix *m=gsl_matrix_alloc(n,n);
  gsl_permutation *p=gsl_permutation_alloc(n);
  int signum=0;
  int i=0;

  gsl_vector *b_vec=gsl_vector_alloc(n);
  gsl_vector *x_sol=gsl_vector_alloc(n);


  for (i = 0; i < n; ++i) {
    gsl_vector_set(b_vec,i,b[i]);
    
    for (int j = 0; j < n; ++j) {
      gsl_matrix_set(m, i, j, A[i][j]);
      
    }
  }

  gsl_linalg_LU_decomp(m,p,&signum);//Realiza descomposición LU
  gsl_linalg_LU_solve(m,p,b_vec,x_sol);//Resuelve el sistema

  for (int i = 0; i < n; ++i) {
    x[i] = gsl_vector_get(x_sol, i);
  }
  
  gsl_matrix_free(m);
  gsl_vector_free(b_vec);
  gsl_permutation_free(p);  
}


void resuelveTridiagonal(double* sol_x, double* diagonal,double* superior, double* inferior, double* vector_b, double n){

  gsl_vector *b    = gsl_vector_alloc(n);
  gsl_vector *x    = gsl_vector_alloc(n);
  
  gsl_vector *diag = gsl_vector_alloc(n);
  gsl_vector *sup  = gsl_vector_alloc(n-1);
  gsl_vector *inf  = gsl_vector_alloc(n-1);

  for (int i = 0; i < n; ++i) {
    gsl_vector_set(diag, i, diagonal[i]);
    gsl_vector_set(b   , i, vector_b[i]);
  }
  
  for (int i = 0; i < n-1; ++i) {
    gsl_vector_set(sup, i, superior[i]);
    gsl_vector_set(inf, i, inferior[i]);
  }
  
  gsl_linalg_solve_tridiag(diag,sup,inf,b,x);

  for (int i = 0; i < n; ++i) {
    sol_x[i]=gsl_vector_get(x,i);
  }

  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(diag);
  gsl_vector_free(sup);
  gsl_vector_free(inf);
}
