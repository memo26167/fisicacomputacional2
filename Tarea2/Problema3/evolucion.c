#include <stdlib.h>
#include <stdio.h>
#include "header.h"


void evolucionFTBS(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  
  for (int i = 1; i < nx; ++i) {
    u_siguiente[i] = (1.0-gamma)*u_anterior[i]+gamma*u_anterior[i-1];
  }
  u_siguiente[0] = u0;
  for (int i = 0; i < nx; ++i) {
    u_anterior[i] = u_siguiente[i];
  }
}

void evolucionFTCS(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  for (int i = 1; i < nx-1; ++i) {
    u_siguiente[i] = u_anterior[i] - gamma/2.0*(u_anterior[i+1] - u_anterior[i-1]);
  }
  u_siguiente[0] = u0;
  u_siguiente[nx-1] = 2.0*u_siguiente[nx-2] - u_siguiente[nx-3];
  
  for (int i = 0; i < nx; ++i) {
    u_anterior[i] = u_siguiente[i];
  }
}

void evolucionLAX(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  for (int i = 1; i < nx-1; ++i) {
    u_siguiente[i] = (u_anterior[i+1]+u_anterior[i-1])/2.0
      - gamma/2.0*(u_anterior[i+1]-u_anterior[i-1]);
    /* u_siguiente[i] = u_anterior[i] */
    /*   - gamma/2.0*(u_anterior[i+1] - u_anterior[i-1]) */
    /*   + pow(gamma,2.0)/2.0*(u_anterior[i+1] - 2*u_anterior[i] + u_anterior[i-1]); */
  }
  u_siguiente[0] = u0;
  u_siguiente[nx-1] = 2.0*u_siguiente[nx-2] - u_siguiente[nx-3];
  
  for (int i = 0; i < nx; ++i) {
    u_anterior[i] = u_siguiente[i];
  }
}

void evolucionBTCS(double*u_anterior,int nx,double gamma,double u0)
{
  double diag[nx];
  double u_diag[nx-1]; //diagonal superior
  double l_diag[nx-1]; //diagon inferior
  double b[nx];
  double* u_siguiente=malloc(sizeof(double) * nx);

  for (int i = 0; i < nx-1; ++i) {
    diag[i]=1;
  }
  diag[nx-1]=gamma+1;

  for (int i = 0; i < nx-2; ++i) {
    l_diag[i]=-gamma/2;
  }
  l_diag[nx-2]=-gamma;
  
  for (int i = 1; i < nx-1; ++i) {
    u_diag[i]=gamma/2;
  }
  u_diag[0]=0;

  for (int i = 1; i < nx-1; ++i) {
    b[i]= u_anterior[i];
  }
  b[0]=u0;

  //resolver sistema tridiagonal
  resuelveTridiagonal(u_siguiente,diag,u_diag,l_diag,b,nx);

  for (int i = 0; i < nx; ++i) {
    u_anterior[i] = u_siguiente[i];
  }
  free(u_siguiente);
}
