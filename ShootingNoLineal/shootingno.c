#include <stdio.h>
#include <stdlib.h>
#include "header.h"


void funcionSistema(double *k,double x, double *y, FC fn, FC fy, FC fyp)
{ //y[0] es y, y[1] es y'
  //y[2] es z, y[3] es z'
  //Primer sistema, que resuelve problema de valor inicial con parametro tk
  k[0]=y[1];
  k[1]=fn(x,y[0],y[1]);
  //Segundo sistema, que obtiene valores de la derivada de y respecto a tk 
  k[2]=y[3];
  k[3]=fy(x,y[0],y[1])*y[2]+fyp(x,y[0],y[1])*y[3];
}

void shooting(double xi, double xf, double * cond_fronteras, double paso_h, pRKF problema, FILE* fp)
{
  // Se redondea hacia arriba para obtener la mayor cantidad de componentes
  int num_it=ceil((xf-xi)/paso_h);
  // Se recalcula el paso para asegurar que se obtenga el tiempo final
  paso_h=(xf-xi)/num_it;
  ++num_it;
  
  //Se define un vector para almacenar la solucion para el problema de frontera
  gsl_vector *solucion_frontera=gsl_vector_alloc(num_it);
  gsl_vector_set_zero(solucion_frontera);
  // Vector con contenido del valor de la variable independiente de la iteracion
  // que da solucion al problema de el paso adaptativo, puede representar el tiempo
  gsl_vector *sol_x=gsl_vector_alloc(num_it);
  gsl_vector_set_zero(sol_x);
  gsl_vector_set(sol_x,0,xi);

  /* Se definen el problema de valor inicial
   * Se asigna
   * -Parametro de condicion inicial
   * -Funcion f vectorial del sistema de EDO de problema de valor inicial
   * -Funcion fn escalar (funcion del problema de frontera)
   * -Funcion fy escalar (derivada parcial de F con respecto a y)
   * -Funcion fyp escalar (derivada parcial de F con respecto a y')
   * -Condiciones iniciales del problema de valor inicial
   * -Matriz solucion del problema de valor inicial
   */
  
  problema.f=funcionSistema;
  gsl_matrix *sol_y=gsl_matrix_alloc(num_it,4);
  double tka=10;
  double tks=1;
  double cond_iniciales[4]={cond_fronteras[0],tka,0.0,1.0};
  double aux=0;
  while (fabs(tka-tks)>pow(10,-3)){
    tka=tks;
    cond_iniciales[1]=tka;  
    // Se resuelve con el metodo RKF
    rk(cond_iniciales,4,num_it,paso_h,problema,sol_y,sol_x);
    // Se calcula la solucion del problema de frontera

    aux=(gsl_matrix_get(sol_y,num_it-1,0)-cond_fronteras[1])/gsl_matrix_get(sol_y,num_it-1,2);
    tks=tka-aux;
    printf("%f\n", tka);
  }
  aux=0;
  for (int i = 0; i < num_it; ++i) {
    aux=gsl_matrix_get(sol_y,i,0);
    gsl_vector_set(solucion_frontera,i,aux);
  }

  for (int i = 0; i < num_it; ++i) {
    fprintf(fp, "%E ", gsl_vector_get(sol_x,i));
    fprintf(fp, "%E", gsl_vector_get(solucion_frontera,i));
    fprintf(fp,"\n");
  }
  
  // Se libera 
  gsl_matrix_free(sol_y);
  gsl_vector_free(solucion_frontera);
  gsl_vector_free(sol_x);
}
