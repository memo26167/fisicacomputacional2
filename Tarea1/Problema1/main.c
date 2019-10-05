#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* Este programa resuelve el problema no lineal, de frontera de segundo orden
 * y''=f(xi,y(xi),y'(xi))
 * s.t.  y(a)=alpha
 *       y(b)=beta
 *       a<=x<=b
 * Mediante la resolucion de un sistema de ecuaciones lineales que
 * discretizan este problema mediante formulas de derivada numerica centrada,
 * y la aplicacion de el metodo de Newton Multivariable con aplicacion
 * de la matriz Jacobiana.
 * La discretizacion divide el intevalo [a,b] en n+1 sub intervalos
 * donde el paso h=(a-b)/(n+1)
 * Si anotamos wi=y(xi), entonces
 * -wi+1 + 2wi - wi-1 + h²f(xi,wi,(wi+1-wi-1)/(2h))=0
 * w0 = alpha
 * wn+1 = beta
 * La matriz jacobiana J
 * J(w1,..,wn)ij={ -1+h/2*fyp* (xi,wi,(wi+1-wi-1)/(2h)), para i=j-1, j=2..N
                    2+h²/2*fy* (xi,wi,(wi+1-wi-1)/(2h)), para i=j  , j=1..N
		   -1-h/2*fyp* (xi,wi,(wi+1-wi-1)/(2h)), para i=j+1, j=1..N-1
 * El vector b
 * b=-[ -w2   + 2w1 - w0   + h² f(x1,w1, (w2-w0)/(2h) ), ... ,
        -wi+1 + 2wi - wi-1 + h² f(xi,wi, (wi+1-wi-1)/(2h) ), ... ,
	-wn+1 + 2wn - wn-1 + h² f(xn,wn, (wn+1-wn-1)/(2h) )]
      
 * v será un vector con el contenido de y(x) que resuelve el problema
 * Jv=b
 * Se resuelve este problema de forma iterativa,
 * en la que despues de calcular v, se obtiene el nuevo w tal que
 * w=w+v
 * Despues de n iteraciones, w será el vector solucion al problema
*/

double funcion(double x,double y, double yp)
{
  double salida;
  /*Funcion*/
  double k_termic=120;
  double area=1.5e-4;
  double perimetro=0.106;
  double h_convec=10;
  double temp_inf=293;
  double stefan_boltzmann=5.670373e-8;
  double aux1=perimetro*h_convec/(area*k_termic)*(y-temp_inf);
  double aux2=perimetro*stefan_boltzmann/(area*k_termic)*(pow(y,4)-pow(temp_inf,4));
  salida=aux1+aux2;
  /* */

  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

//derivada parcial de la funcion con respecto a y
double funcion_y(double x,double y, double yp)
{
  double salida;
  /*Funcion*/
  double k_termic=120;
  double area=1.5e-4;
  double perimetro=0.106;
  double h_convec=10;
  double stefan_boltzmann=5.670373e-8;
  double aux1=perimetro*h_convec/(area*k_termic);
  double aux2=4*perimetro*stefan_boltzmann/(area*k_termic)*pow(y,3);
  salida=aux1+aux2;
  /* */
  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

//derivada parcial de la funcion con respecto a y'
double funcion_yp(double x,double y, double yp)
{
  double salida;
  salida=0;
  //confundir warnings
  x=yp*0+x;
  x=y*0+x;
  return salida;
}

int main(void)
{
  // Valores del intervalo de a<x<b
  double xa=0;
  double xb=0.1;
  int particiones=1000;

  // condiciones de frontera
  double alpha=773.15;
  double beta=0;
  double cond[2]={alpha,beta};

  // definicion del problema de frontera
  pDFNL problema;
  problema.fn=&funcion;
  problema.fy=&funcion_y;
  problema.fyp=&funcion_yp;
  problema.cond=cond;
  problema.n=particiones;

  // Se abre el archivo para imprimir los datos y se realiza el metodo
  FILE *archivo = fopen("salida.dat", "w");
  difFinitasNoLineales(problema,xa,xb,1e-3,archivo);
  // Se termina el metodo y se cierra el archivo
  fclose(archivo);
  
  return 0;
}


