#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* 
 * Este código realiza el método de Runge-Kutta-Fehlberg
 * de orden 5, con correcciones de Cash y Karp.
 * Toma el problem RKF y mediante el metodo
 * fEqDiferencial, lo trata como un problema a resolver del tipo Runge-Kutta estandar
 */
void rk
(double * y0,int num_eq, int num_it, double paso_h, pRKF funcion, gsl_matrix *sol_y,gsl_vector *sol_x)
{
  // Matriz con contenido de la solucion y
  // su primra componetente es la iteracion y su segunda es el num de el component de la solucion
  gsl_matrix_set_zero(sol_y); 

  // Condiciones iniciales
  condIniciales(sol_y,y0,num_eq);
  
  // Se obtiene la solucion del metodo llamando a la funcion RKF
  // por medio de la matriz sol_y y el vector sol_x
  metodoRKF(sol_y, sol_x, num_it,num_eq,paso_h,funcion);
}

void condIniciales(gsl_matrix * y, double *y0, int num_eq)
{
  for (int i = 0; i < num_eq; ++i) {
    gsl_matrix_set(y,0,i,y0[i]);
    gsl_matrix_set(y,0,i,y0[i]);
  }  
}

/* 
 * Esta función realiza el método de RKF
 * Realiza la iteracion del método RK5 junto al paso adaptativo
 */
void metodoRKF(gsl_matrix *y, gsl_vector *x, int num_it, int num_eq,double paso_h_anterior, pRKF funcion)
{
  // Se define la tabla de valores de Cash y Karp
  double tabla_a[6]={0 ,1.0/5.0 ,3.0/10.0 ,3.0/5.0 ,1.0 ,7.0/8.0};
  double tabla_b[6][5]={
		  {0              ,0           ,0             ,0                ,0},
		  {1.0/5.0        ,0           ,0             ,0                ,0},
		  {3.0/40.0       ,9.0/40.0    ,0             ,0                ,0},
		  {3.0/10.0       ,-9.0/10.0   ,6.0/5.0       ,0                ,0},
		  {-11.0/54.0     ,5.0/2.0     ,-70.0/27.0    ,35.0/27.0        ,0},
		  {1631.0/55296.0 ,175.0/512.0 ,575.0/13824.0 ,44275.0/110592.0 ,253.0/4096.0}};
  double tabla_c[6]={37.0/378.0 ,0 ,250.0/621.0 ,125.0/594.0 ,0 ,512.0/1771.0};
  
  // Almacena los vectores K  del metodo RKF
  gsl_matrix *vectores_k=gsl_matrix_alloc(6,num_eq);

  // Almacena la diferencia componente a componente
  // entre la solucion de grado 5 y de grado 4, mediante las tablas
  gsl_vector *vector_error=gsl_vector_alloc(num_eq);
  
  // El paso h que se adaptará en cada iteracion
  double paso_h;

  // variables acumulacion y auxilio, se reciclaran
  gsl_vector *acum=gsl_vector_alloc(num_eq);
  gsl_vector *yAnt=gsl_vector_alloc(num_eq);
  gsl_vector *kv=gsl_vector_alloc(num_eq);
  gsl_vector *ku=gsl_vector_alloc(num_eq);
  
  for (int t=0; t < num_it-1; ++t) {
    //asignar valores a la variable independiente (paso h antiguo)
    gsl_vector_set(x,t+1, paso_h+gsl_vector_get(x,t));
    gsl_matrix_get_row(yAnt,y,t);
    paso_h=paso_h_anterior;
    gsl_matrix_set_zero(vectores_k);
    gsl_vector_set_zero(acum);
    gsl_vector_set_zero(ku);
    gsl_vector_set_zero(kv);
    //Calcular K0
    gsl_vector_set_zero(kv);
    fEqDiferencial(kv,gsl_vector_get(x,t),yAnt,num_eq,funcion);
    gsl_vector_scale(kv,paso_h);
    gsl_matrix_set_row(vectores_k,0,kv);

    //Calcular el resto de Ki's
    for (int i = 1; i < 6; ++i) {
      gsl_vector_set_zero(ku);
      gsl_vector_set_zero(acum);
      for (int j = 0; j < i; ++j) {
	gsl_matrix_get_row(kv,vectores_k,j);
	gsl_vector_scale(kv,tabla_b[i][j]);
	gsl_vector_add(acum,kv);
      }
      gsl_vector_add(acum,yAnt);
      fEqDiferencial(ku,gsl_vector_get(x,t)+tabla_a[i]*paso_h,acum,num_eq,funcion);
      gsl_vector_scale(ku,paso_h);
      gsl_matrix_set_row(vectores_k,i,ku);
    }

    //Calcular Soluciones (paso h antiguo)
    gsl_vector_set_zero(acum);
    for (int i = 0; i<6; ++i) {
      gsl_matrix_get_row(kv,vectores_k,i);
      gsl_vector_scale(kv,tabla_c[i]);
      gsl_vector_add(acum,kv);
    }
    gsl_vector_add(acum,yAnt);
    gsl_matrix_set_row(y,t+1,acum);
  }

  //Terminar algoritmo
  gsl_matrix_free(vectores_k);
  gsl_vector_free(vector_error);
  gsl_vector_free(acum);
  gsl_vector_free(ku);
  gsl_vector_free(kv);
  gsl_vector_free(yAnt);
}

/* Esta metodo toma la funcion f de la ecuaciones diferencial en terminos
 * de vectores tipo puntero double
 * y lo transforma a vectores tipo gsl_vector
 */
void fEqDiferencial(gsl_vector *sol_k,double x, gsl_vector *sol_y,int num_eq, pRKF funcion)
{
  double k[num_eq];
  double y[num_eq];
  for (int i = 0; i < num_eq; ++i) {
    y[i]=gsl_vector_get(sol_y,i);
    k[i]=0;
  }
  funcion.f(k,x,y,funcion.fn,funcion.fy,funcion.fyp);
  for (int i = 0; i < num_eq; ++i) {
    gsl_vector_set(sol_k,i,k[i]);
  }
}
