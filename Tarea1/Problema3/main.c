#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

/* Este programa realiza metodo de Numeov
 * Para resolver el problema de autovalores
 * De Schrödinger
 * 
 */

double potencial(double x){
  double pot;
  pot=50*exp(-pow(x-0.5,2)/0.08);
  return pot;
}

int main(void)
{
  int N=100;
  int numero_autofunciones=5;
  
  double xa=0;
  double xb=1;
  double paso_h=(xb-xa)/N;
  double ya0=0;
  double ya1=2.0/N;//
  double nmax=1000;
  double tol=1e-10;

  FILE *archivo = fopen("salida.dat", "w");

  numerov(numero_autofunciones,xa,ya0,ya1,N,paso_h,nmax,tol,archivo);
  
  fclose(archivo);	
  return 0;  
}




