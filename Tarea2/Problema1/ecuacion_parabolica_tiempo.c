#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void resolver_metodo(ecp p)
{
  //Se definen las variables del problema
  int n=p.nx;
  int m=p.ny;

  double delta_x=(p.xf-p.xi)/(p.nx-1);
  double delta_y=(p.yf-p.yi)/(p.ny-1);
  double delta_t=0.3*pow(delta_x,2)/p.D;
  double gammax=p.D*delta_t/(pow(delta_x,2));
  double gammay=p.D*delta_t/(pow(delta_y,2));  

  double mA[n][m];//Matriz Actual
  double mS[n][m];//Matriz Siguiente

  //Grilla;
  double x[n];
  double y[m];
  for (int i = 0; i < n; ++i) {
    x[i]=p.xi+delta_x*i;
  }
  for (int j = 0; j < m; ++j) {
    y[j]=p.yi+delta_y*j;
  }
  
  //condicion inicial
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      mA[i][j]=p.m_inicial[i][j];
    }    
  }

  //puntos siguientes fuera de los bordes
  for (int i = 1; i < n-1; ++i) {
    for (int j = 1; j < m-1; ++j) {
      mS[i][j]=(1 - 2*gammax - 2*gammay)*mA[i][j]
	+ gammax*(mA[i+1][j] + mA[i-1][j])
	+ gammay*(mA[i][j+1] + mA[i][j-1])
	+ p.fuente(x[i],y[j],p.parametros)*delta_t;
    }
  }
  
  /* Bordes sin esquinas */

  double phi0j=0;
  double phiNxj=0;
  double phii0=0;
  double phiiNy=0;
  //Borde izquierdo i=-1
  for (int j = 1; j < m-1; ++j) {
    phi0j=mA[1][j]
      -2*delta_x /( p.beta[0] * p.g1(y[j],p.parametros) )
      * ( p.f1(y[j],p.parametros) - p.alpha[0]*mA[0][j] );
    
    mS[0][j]=(1 - 2*gammax - 2*gammay)*mA[0][j]
	+ gammax*(mA[0+1][j] + phi0j)
	+ gammay*(mA[0][j+1] + mA[0][j-1])
	+ p.fuente(x[0],y[j],p.parametros)*delta_t;
  }

  //Borde derecho i=n
  for (int j = 1; j < m-1; ++j) {
    phiNxj=mA[n-2][j]
      -2*delta_x /( p.beta[1] * p.g2(y[j],p.parametros) )
      * ( p.f2(y[j],p.parametros) - p.alpha[1]*mA[n-1][j] );
    
    mS[n-1][j]=(1 - 2*gammax - 2*gammay)*mA[n-1][j]
	+ gammax*(phiNxj + mA[n-2][j])
	+ gammay*(mA[n-1][j+1] + mA[n-1][j-1])
	+ p.fuente(x[n-1],y[j],p.parametros)*delta_t;
  }

  //Borde derecho j=-1
  //Borde superior j=n
}
