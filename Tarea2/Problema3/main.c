#include <stdio.h>
#include <stdlib.h>
#include "header.h"

int main()
{
  //Variables del método
  int    nx = 501;
  double xi = 0.0;
  double xf = 100.0;
  //minimo en el cual se dice que la onda llegó al destino
  double tol_onda = 5e-3; 
  double a  = 1.0;

  double dx = (xf-xi)/(nx-1);
  double dt = 1.0/10000.0;
  double gamma = dt/dx*a;
  double u0 = 0;

  double x[nx];
  double *u_FTBS = malloc(sizeof(double) * nx);
  double *u_FTCS = malloc(sizeof(double) * nx);
  double *u_BTCS = malloc(sizeof(double) * nx);
  double *u_LAX  = malloc(sizeof(double) * nx);
  
  
  for (int i = 0; i < nx; ++i) {
    x[i] = dx*i;
   }

  rellenarUinicial(u_FTBS, x, nx, u0);
  rellenarUinicial(u_FTCS, x, nx, u0);
  rellenarUinicial(u_BTCS, x, nx, u0);
  rellenarUinicial(u_LAX,  x, nx, u0);

  FILE*archivo=fopen("datos.dat","w");

  // Posicion en x, en la que se detendrá la evolucion
  double porcentaje=25;
  int pos_corte =(nx-2)*porcentaje/100;
  while(u_FTBS[pos_corte]>=0 && u_FTBS[pos_corte]<tol_onda){
    evolucionFTBS(u_FTBS, nx, gamma, u0);
    evolucionFTCS(u_FTCS, nx, gamma, u0);
    evolucionBTCS(u_BTCS, nx, gamma, u0);
  }
  dt=0.1;
  gamma = dt/dx*a;
  while(u_LAX[pos_corte]>=0 && u_LAX[pos_corte]<tol_onda){
    evolucionLAX (u_LAX,  nx, gamma, u0);
  }
  
  for (int i = 0; i < nx; ++i) {
    fprintf(archivo, "%E ",   x[i]);
    fprintf(archivo, "%E ",   u_FTBS[i]);
    fprintf(archivo, "%E ",   u_FTCS[i]);
    fprintf(archivo, "%E ",   u_BTCS[i]);
    fprintf(archivo, "%E \n", u_LAX[i]);
  }
  fprintf(archivo, "\n");

  free(u_FTBS);
  free(u_BTCS);
  free(u_FTCS);
  free(u_LAX);
  fclose(archivo);
  return 0;
}

void rellenarUinicial(double*u,double*x,int nx,double u0)
{
  for (int i = 0; i < nx; ++i) {
    if(x[i]>=0.0 && x[i]<=2.0){
      u[i] = 1;
    }
    else{
      u[i] = 0;
    }
  }
  u[0] = u0;
}




