#include <stdio.h>
#include <stdlib.h>
#include "header.h"


/* Este programa resuelve el problema de frontera de segundo orden
 * y''=p(x)y'+q(x)y+r(x)
 * s.t.  y(a)=alpha
 *       y(b)=beta
 *       a<=x<=b
 * Mediante la resolucion de un sistema de ecuaciones lineales que
 * discretizan este problema mediante formulas de derivada numerica centrada
 * La discretizacion divide el intevalo [a,b] en n+1 sub intervalos
 * donde el paso h=(a-b)/(n+1)
 * Si anotamos wi=y(xi), entonces
 * -(1+h/2 p(xi))wi-1 + (2+h² q(xi))wi - (1-h/2 p(xi))wi+1 = -h² r(xi)
 * w0 = alpha
 * wn+1 = beta
 *
 * w será un vector con el contenido de y(x) que resuelve el problema
 */

double funcionP(double x)
{
  double p;
  p=-2/x*0;
  return p;
}

double funcionQ(double x)
{
  double q;
  q=-2*5+x*0;
  return q;
}

double funcionR(double x)
{
  double r;
  r=x*0+1e-200;
  return r;
}

int main(void)
{
  // Valores del intervalo de a<x<b
  double xa=0;
  double xb=1;
  int particiones=1000;

  double alpha=0;
  double beta=0;
  double cond[2]={alpha,beta};// condiciones de frontera
  
  pEDOF problema_frontera;
  
  problema_frontera.p=&funcionP;
  problema_frontera.q=&funcionQ;
  problema_frontera.r=&funcionR;
  problema_frontera.cond=cond;
  problema_frontera.n=particiones;
  // a<x<b
  
  // Se abre el archivo para imprimir los datos y se realiza el metodo
  FILE *archivo = fopen("salida.dat", "w");
  // Se termina el metodo y se cierra el archivo
  diferenciasFinitas(problema_frontera,xa,xb,archivo);
  
  fclose(archivo);
  return 0;
}
