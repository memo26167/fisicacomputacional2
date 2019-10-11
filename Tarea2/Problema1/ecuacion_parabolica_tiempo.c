#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void resolver_metodo(ecp p, FILE* archivo)
{
  //Se definen las variables del problema
  int n=p.nx;
  int m=p.ny;
  
  double delta_x=(p.xf-p.xi)/(p.nx-1);
  double delta_y=(p.yf-p.yi)/(p.ny-1);
  double delta_t=0.02*pow(delta_x,2)/p.D;
  p.nt=(p.tf-p.ti)/(delta_t)+1;
  int count_t=0;
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

  
  printf("%f %f %f %d\n", p.ti,p.tf,delta_t,p.nt);

  for (double t = p.ti; t < p.tf; t+=delta_t , count_t+=1 ) {//iteracion de puntos
    // Puntos a un tiempo siguiente, fuera de los bordes
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
    //Borde izquierdo i=0
    for (int j = 1; j < m-1; ++j) {
      phi0j=mA[1][j]
	-2*delta_x /( p.beta[0] * p.g1(t,y[j],p.parametros) )
	* ( p.f1(t,y[j],p.parametros) - p.alpha[0]*mA[0][j] );

      //Dirichlet
      if(p.beta[0]==0){
	phi0j= p.f1(t,y[j],p.parametros)/p.alpha[0];
      }
      
      mS[0][j]=(1 - 2*gammax - 2*gammay)*mA[0][j]
	+ gammax*(mA[0+1][j] + phi0j)
	+ gammay*(mA[0][j+1] + mA[0][j-1])
	+ p.fuente(x[0],y[j],p.parametros)*delta_t;
    }

    //Borde derecho i=n
    for (int j = 1; j < m-1; ++j) {
      phiNxj=mA[n-2][j]
	+2*delta_x /( p.beta[1] * p.g2(t,y[j],p.parametros) )
	* ( p.f2(t,y[j],p.parametros) - p.alpha[1]*mA[n-1][j] );

      //Dirichlet
      if(p.beta[1]==0){
	phiNxj= p.f2(t,y[j],p.parametros)/p.alpha[1];
      }
      
      mS[n-1][j]=(1 - 2*gammax - 2*gammay)*mA[n-1][j]
	+ gammax*(phiNxj + mA[n-2][j])
	+ gammay*(mA[n-1][j+1] + mA[n-1][j-1])
	+ p.fuente(x[n-1],y[j],p.parametros)*delta_t;
    }

    //Borde inferior j=0
    for (int i = 1; i < n-1; ++i) {
      phii0=mA[i][1]
	-2*delta_y /( p.beta[2] * p.g3(t,x[i],p.parametros) )
	* ( p.f3(t,x[i],p.parametros) - p.alpha[2]*mA[i][0] );

      //Dirichlet
      if(p.beta[2]==0){
	phii0=p.f3(t,x[i],p.parametros) / p.alpha[2];
      }
      
      mS[i][0]=(1 - 2*gammax - 2*gammay)*mA[i][0]
	+ gammax*(mA[i+1][0] + mA[i-1][0])
	+ gammay*(mA[i][0+1] + phii0)
	+ p.fuente(x[i],y[0],p.parametros)*delta_t;
    }
  
    //Borde superior j=n
    for (int i = 1; i < n-1; ++i) {
      phiiNy=mA[i][m-2]
	+2*delta_y /( p.beta[3] * p.g4(t,x[i],p.parametros) )
	* ( p.f4(t,x[i],p.parametros) - p.alpha[3]*mA[i][m-1] );

      //Dirichlet
      if (p.beta[3]==0){
	phiiNy=p.f4(t,x[i],p.parametros) / p.alpha[3];
      }
      
      mS[i][m-1]=(1 - 2*gammax - 2*gammay)*mA[i][m-1]
	+ gammax*(mA[i+1][m-1] + mA[i-1][m-1])
	+ gammay*(phiiNy + mA[i][m-2])
	+ p.fuente(x[i],y[m-1],p.parametros)*delta_t;
    }

    /* Esquinas */

    //Inferior Izquierda, i=0 j=0
    double phi01=mA[1][0]
      -2*delta_x /( p.beta[0] * p.g1(t,y[0],p.parametros) )
      * ( p.f1(t,y[0],p.parametros) - p.alpha[0]*mA[0][0] );

    //Dirichlet
    if(p.beta[0]==0){
      phi01= p.f1(t,y[0],p.parametros) / p.alpha[0];
    }
    
    double phi10=mA[0][1]
      -2*delta_y /( p.beta[2] * p.g3(t,x[0],p.parametros) )
      * ( p.f3(t,x[0],p.parametros) - p.alpha[2]*mA[0][0] );

    //Dirichlet
    if(p.beta[2]==0){
      phi10=p.f3(t,x[0],p.parametros) / p.alpha[2];
    }
    
    mS[0][0]=(1 - 2*gammax - 2*gammay)*mA[0][0]
      + gammax*(mA[0+1][0] + phi01)
      + gammay*(mA[0][0+1] + phi10)
      + p.fuente(x[0],y[0],p.parametros)*delta_t;

    //Superior Izquierda, i=0 j=m-1
    double phi0Ny=mA[1][m-1]
      -2*delta_x /( p.beta[0] * p.g1(t,y[0],p.parametros) )
      * ( p.f1(t,y[0],p.parametros) - p.alpha[0]*mA[0][0] );

    //Dirichlet
    if(p.beta[0]==0){
      phi0Ny=p.f1(t,y[0],p.parametros) / p.alpha[0];
    }
    
    double phi1Ny=mA[0][m-2]
      +2*delta_y /( p.beta[3] * p.g4(t,x[0],p.parametros) )
      * ( p.f4(t,x[0],p.parametros) - p.alpha[3]*mA[0][m-1] );

    //Dirichlet
    if(p.beta[3]==0){
      phi1Ny=p.f4(t,x[0],p.parametros) / p.alpha[3] ;
    }
    
    mS[0][m-1]=(1 - 2*gammax - 2*gammay)*mA[0][m-1]
      + gammax*(mA[0+1][m-1] + phi0Ny)
      + gammay*(phi1Ny + mA[0][m-2])
      + p.fuente(x[0],y[m-1],p.parametros)*delta_t;

    //Inferior Derecha, i=n-1 j=0
    double phiNx0=mA[n-1][1]
      -2*delta_y /( p.beta[2] * p.g3(t,x[n-1],p.parametros) )
      * ( p.f3(t,x[n-1],p.parametros) - p.alpha[2]*mA[n-1][0] );

    //Dirichlet
    if(p.beta[2]==0){
      phiNx0=p.f3(t,x[n-1],p.parametros) / p.alpha[2] ;
    }
    
    double phiNx1=mA[n-2][0]
      +2*delta_x /( p.beta[1] * p.g2(t,y[0],p.parametros) )
      * ( p.f2(t,y[0],p.parametros) - p.alpha[1]*mA[n-1][0] );

    //Dirichlet
    if(p.beta[1]==0){
      phiNx1=p.f2(t,y[0],p.parametros) / p.alpha[1] ;
    }
    
    mS[n-1][0]=(1 - 2*gammax - 2*gammay)*mA[n-1][0]
      + gammax*(phiNx1 + mA[n-2][0])
      + gammay*(mA[n-1][0+1] + phiNx0)
      + p.fuente(x[n-1],y[0],p.parametros)*delta_t;

    //Superior Derecha, i=n-1 j=m-1
    double phiNx1Ny=mA[n-2][m-1]
      +2*delta_x /( p.beta[1] * p.g2(t,y[m-1],p.parametros) )
      * ( p.f2(t,y[m-1],p.parametros) - p.alpha[1]*mA[n-1][m-1] );

    //Dirichlet
    if(p.beta[1]==0){
      phiNx1Ny=p.f2(t,y[m-1],p.parametros) / p.alpha[1] ;
    }
    
    double phiNxNy1=mA[n-1][m-2]
      +2*delta_y /( p.beta[3] * p.g4(t,x[n-1],p.parametros) )
      * ( p.f4(t,x[n-1],p.parametros) - p.alpha[3]*mA[n-1][m-1] );

    //Dirichlet
    if(p.beta[3]==0){
      phiNxNy1=p.f4(t,x[n-1],p.parametros) / p.alpha[3];
    }
    
    mS[n-1][m-1]=(1 - 2*gammax - 2*gammay)*mA[n-1][m-1]
      + gammax*(phiNx1Ny + mA[n-2][m-1])
      + gammay*(phiNxNy1 + mA[n-1][m-2])
      + p.fuente(x[n-1],y[m-1],p.parametros)*delta_t;

    // Imprimir archivo
    if(count_t%30==0){
      for (int i = 0; i < n; ++i) {
	for (int j = 0; j < m; ++j) {
	  fprintf(archivo, "%E ", x[i]);
	  fprintf(archivo, "%E ", y[j]);
	  fprintf(archivo, "%E\n",mS[i][j]);
	}
      }
      fprintf(archivo, "\n");
    }


    // Actualizar
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
	mA[i][j]=mS[i][j];
      }
    }

  }
  
    
}
