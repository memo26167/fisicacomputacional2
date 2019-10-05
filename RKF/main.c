#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/*
 Este es el sistema de ecuaciones diferenciales, que es la funcion F para el metodo.
 Toma como entrada la variable independiente (o tiempo)
 y la solucion de la ecuacion diferencial en un determinado instante, 
 y de salida entrega un vector del mismo tamaño del sistema de eq diferenciales
 por medio de un puntero dado en la entrada.
 */

void funcion(double* k, double x, double *y)
{
  k[0]=y[1];//siempre quedará asi para una edo de 2do orden
  k[1]=-y[0]+y[1];

  //para evitar warnings
  x=x+1;
}

int main(void)
{
  double y0=0;
  double y0p=10;
  double y[2]={y0,y0p};// cond. iniciales

  // Puntero que guarda la funcion F del sismtema de ecuaciones diferenciales
  void (*puntero_funcion)(double*, double, double *)=&funcion;

  /* No se controlará el tiempo,sino que el número de iteraciones.
     Esto debido a que como el paso de tiempo h cambiará, no se puede
     determinar el tiempo del último componente.
     Entonces se preferirá determinar un total de iteraciones, que dará
     el tamaño del vector solucion.
   */
  int num_it=10000;
  int num_eq=2;

  // Paso de tiempo
  double paso_h=0.5;

  // Tolerancia de error;
  double tol=pow(10,-3);

  // Se imprime la solucion en un archivo para su análisis posterior
  FILE *archivo = fopen("salida.dat", "w");

  rk(y,num_eq,num_it,paso_h,tol,puntero_funcion,archivo);
  
  fclose(archivo);
  return 0;
}

/* void funcionEqDiferencial(gsl_vector *k,double x, gsl_vector *y) */
/* { */
/*   //gsl_vector_set(k,0,gsl_vector_get(y,1)); */
/*   //gsl_vector_set(k,1,-9.8+pow(65.351,-3)*pow(gsl_vector_get(y,0),2)*exp(-pow(10.53,-5)*gsl_vector_get(y,0))); */
/*   gsl_vector_set(k,0,-gsl_vector_get(y,1)); */
/*   gsl_vector_set(k,1,gsl_vector_get(y,0)); */

/*   //para evitar warnings */
/*   x=x+1; */
/* } */
