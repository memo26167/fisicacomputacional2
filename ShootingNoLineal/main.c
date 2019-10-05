#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* Este programa resuelve el problema de frontera de segundo orden
 * y''=f(x,y,y')
 * s.t.  y(a)=alpha
 *       y(b)=beta
 *       a<=x<=b
 * Mediante la resolucion de dos problemas de valor inicial junto a RK5
 * 
 * y''=f(x,y,y')     z''=df/dy*z+df/dy' *z'
 * s.t.  y(a)=alpha  s.t.  z(a)=0
 *       y'(a)=tk          z'(a)=1
 *       a<=x<=b           a<=x<=b
 *
 * Para luego obtener un tk actualizado con el ecuacion
 * tk+1=tk - (y(b)-beta)/z(b) 
 * El y obtenido es la solucion
 */

double funcion(double x,double y, double yp)
{
  double salida;
  salida=-2/x*yp;
  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

//derivada parcial de la funcion con respecto a y
double funcion_y(double x,double y, double yp)
{
  double salida;
  salida=0;
  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

//derivada parcial de la funcion con respecto a y'
double funcion_yp(double x,double y, double yp)
{
  double salida;
  salida=-2/x;
  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

int main(void)
{
  pRKF problema;
  problema.fn=&funcion;
  problema.fy=&funcion_y;
  problema.fyp=&funcion_yp;
  double alpha=110;
  double beta=0;
  double cond_frontera[2]={alpha,beta};
  double paso_h=0.01;
  double a=2;//xi
  double b=4;//xf
  // Se abre el archivo para imprimir los datos y se realiza el metodo
  FILE *archivo = fopen("salida.dat", "w");
  
  shooting(a,b,cond_frontera,paso_h,problema,archivo);
  fclose(archivo);
  return 0;
}
