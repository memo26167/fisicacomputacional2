/* Este programa resuelve Ecuacion Diferencial Parcial Hiperbolica de segundo orden, en una dimension espacial,
 * es decir:
 * 
 * ∂²U(x,t)/∂t² - k² ∂²U(x,t)/∂x² = 0
 * 
 * Sujeta a las condiciones de borde
 * 
 * Para el tiempo
 *    U(x,0)  = f(x)
 * ∂U(x,0)/∂t = g(x)
 * 
 * Para el espacio
 * U(0,t) = α(t)
 * U(L,t) = β(t)
 *
 * En donde se utiliza diferencias finitas para discretizar el problema
 *
 * ( U[x,t+1] -2U[x,t] + U[x,t-1] )/ Δt² - k²( U[x+1,t] -2U[x,t] + U[x-1,t] )/ Δx² = 0
 *
 * Que se puede resumir en la ecuacion vectorial:
 * 
 * U⃗[t+1]= A U⃗[t] - U⃗[t-1]
 * Sin embargo, si t=1, entonces se debe ocupar la condicion de borde.
 *
 * ∂U(x,0)/∂t = g(x)
 * (U[x,1]-U[x,0])/Δt = g(x)
 * Entonces
 * U[x,1] = g(x)Δt + U[x,0]
 * 
 * I
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"


double alpha1(double t, double* parametros);
double alpha2(double t, double* parametros);
double beta1(double t, double* parametros);
double beta2(double t, double* parametros);
double gamma1(double t, double* parametros);
double gamma2(double t, double* parametros);
double f(double x, double* parametros);
double g(double x, double* parametros);


int main()
{
  
  hiper problema;

  problema.parametros=malloc(sizeof(double)*1);
  problema.parametros[0]=0.9;//p0
  problema.D=1;//el omega cuadrado

  problema.xi=0;
  problema.xf=1;

  // DEFINIR TIEMPOS INICIALES
  problema.ti=0;
  problema.tf=5;//5;
  
  /* condiciones de borde */
  //coeficientes

  // TUBO CERRADO
  // * x=0
  problema.alpha[0]=1;
  problema.beta[0]=0;
  problema.gamma[0]=1;//problema.parametros[0];
  
  // * x=L
  problema.alpha[1]=0;
  problema.beta[1]=1;
  problema.gamma[0]=1;//problema.parametros[0];

  //funciones de condiciones de borde

  //FALTA DETERMINAR LAS FUNCIONES,TAMAÑOS DE GRILLA, ETC
  problema.alpha1=&alpha1;
  problema.alpha2=&alpha2;

  problema.beta1=&beta1;
  problema.beta2=&beta2;

  problema.gamma1=&gamma1;
  problema.gamma2=&gamma2;
  
  // funciones de condiciones iniciales
  problema.f=&f;
  problema.g=&g;
  
  /*Parametros del método*/
  problema.nx=100;

  FILE *archivo = fopen("datos.dat", "w");
  resolver_metodo( problema, archivo);
  fclose(archivo);
  free(archivo);
  
  return 0;
}

double alpha1(double t, double* parametros)
{
  double y=parametros[0]*t*0+1;
  return y;
}

double alpha2(double t, double* parametros)
{
  double y=parametros[0]*t*0+1;
  return y;
}

double beta1(double t, double* parametros)
{
  double y=1+t*0*parametros[0];
  return y;
}

double beta2(double t, double* parametros)
{
  double y=1+t*0*parametros[0];
  return y;
}

double gamma1(double t, double* parametros)
{
  double y=parametros[0]+t*0;
  return y;
}

double gamma2(double t, double* parametros)
{
  double y=parametros[0]+t*0;
  return y;
}

double f(double x, double* parametros)
{
  double pi=3.1415;
  double y=parametros[0]*cos(2*pi*x);
  return y;
}

double g(double x, double* parametros)
{
  double y=0*x*parametros[0];
  return y;
}
