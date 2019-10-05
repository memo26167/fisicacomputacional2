#include <stdio.h>
#include <stdlib.h>
#include "header.h"


void normalizar(gsl_vector* y, int n){
  double A=0;
  double aux;
  // Evitar numeros muy pequeños o muy grandes
  // como 1e-300, dividiremos la solucion por su promedio antes de normalizar
  double sum=0;
  for (int i = 0; i < n; ++i) {
    sum=gsl_vector_get(y, i)+sum;
  }
  sum=sum/n;
  for (int i = 0; i < n; ++i) {
    aux=gsl_vector_get(y,i);
    gsl_vector_set(y,i,aux/sum); 
  }
  // Calcular la norma de cada punto y acumular en A
  for (int i = 0; i < n; ++i) {
    A=pow(gsl_vector_get(y,i),2)+A;
  }
  // Dividir en la norma
  for (int i = 0; i < n; ++i) {
    aux=gsl_vector_get(y,i);
    gsl_vector_set(y, i, aux/sqrt(A));
  }
  
}
