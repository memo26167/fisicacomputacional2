#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* Este programa resuelve el problema de frontera de segundo orden
 * y''=p(x)y'+q(x)y+r(x)
 * s.t.  y(a)=alpha
 *       y(b)=beta
 *       a<=x<=b
 * Mediante la resolucion de dos problemas de valor inicial junto a RK4
 * 
 * y1''=p(x)y1'+q(x)y1+r(x)     y2''=p(x)y2'+q(x)y2
 * s.t.  y1(a)=alpha            s.t.  y2(a)=0
 *       y1'(a)=0                     y2'(a)=1
 *       a<=x<=b                      a<=x<=b
 *
 * Para luego obtener la solucion general del problema de frontera
 * con la ecuacion
 * 
 * y(x)=y1(x)+(beta-y1(b))/y2(b)*y2(x)
 */

double funcionP(double x)
{
  double p;
  p=0*x;
  return p;
}

double funcionQ(double x)
{
  double q;
  q=0*x-2*(60-50*exp(-pow(x-0.5,2)/0.08));
  return q;
}

double funcionR(double x)
{
  double r;
  r=0*x+1e-300;
  return r;
}

int main(void)
{
  // Valores del intervalo de a<x<b
  double intval_a=0;
  double intval_b=1;
  
  double ya=0;
  double yb=0;
  double cond[2]={ya,yb};// condiciones de frontera
  pRKF problema_frontera;
  problema_frontera.p=&funcionP;
  problema_frontera.q=&funcionQ;
  problema_frontera.r=&funcionR;
  double paso_h=0.0001;

  // Se abre el archivo para imprimir los datos y se realiza el metodo
  FILE *archivo = fopen("salida.dat", "w");
  shooting(intval_a, intval_b, cond,paso_h, problema_frontera, archivo);
  // Se termina el metodo y se cierra el archivo
  fclose(archivo);
  return 0;
}
