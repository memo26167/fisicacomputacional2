#include <stdio.h>
#include <stdlib.h>

void rellenarUinicial(double*U,double*x,int nx,double u0);
void evolucionFTBS(double*u_anterior,int nx,double gamma,double u0);
void evolucionFTCS(double*u_anterior,int nx,double gamma,double u0);
void evolucionLAX(double*u_anterior,int nx,double gamma,double u0);


int main()
{

  //Variables del método
  int nx=1001;
  double xi=0.0;
  double xf=100.0;
  //minimo en el cual se dice que la onda llegó al destino
  double tol_onda=5e-3; 
  double a=1.0;

  double dx=(xf-xi)/(nx-1);
  double dt=1.0/200.0;
  double gamma=dt/dx*a;
  double u0=0;
  
  double x[nx];
  double *u_FTBS= malloc(sizeof(double) * nx);
  double *u_FTCS= malloc(sizeof(double) * nx);
  double *u_BTCS= malloc(sizeof(double) * nx);
  double *u_LAX= malloc(sizeof(double) * nx);
  
  for (int i = 0; i < nx; ++i) {
    x[i]=dx*i;
   }

  rellenarUinicial(u_FTBS,x,nx,u0);
  rellenarUinicial(u_FTCS,x,nx,u0);
  rellenarUinicial(u_BTCS,x,nx,u0);
  rellenarUinicial(u_LAX,x,nx,u0);

  FILE*archivo=fopen("datos.dat","w");

  int pos_corte=(nx-1)*0.4;
  while(u_FTBS[pos_corte]>=0 && u_FTBS[pos_corte]<tol_onda){
    evolucionFTBS(u_FTBS,nx,gamma,u0);
    evolucionFTCS(u_FTCS,nx,gamma,u0);
    evolucionLAX(u_LAX,nx,gamma,u0);
  }
  //FALTA BTCS

  for (int i = 0; i < nx; ++i) {
    fprintf(archivo,"%E ",x[i]);
    fprintf(archivo,"%E ", u_FTBS[i]);
    fprintf(archivo,"%E ", u_FTCS[i]);
    fprintf(archivo,"%E \n", u_LAX[i]);
  }
  fprintf(archivo, "\n");

  fclose(archivo);
  return 0;
}

void rellenarUinicial(double*u,double*x,int nx,double u0)
{
  for (int i = 0; i < nx; ++i) {
    if(x[i]>=0.0 && x[i]<=2.0){
      u[i]=1;
    }
    else{
      u[i]=0;
    }
  }
  u[0]=u0;
}

void evolucionFTBS(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  
  for (int i = 1; i < nx; ++i) {
    u_siguiente[i]=(1.0-gamma)*u_anterior[i]+gamma*u_anterior[i-1];
  }
  u_siguiente[0]=u0;
  for (int i = 0; i < nx; ++i) {
    u_anterior[i]=u_siguiente[i];
  }
}

void evolucionFTCS(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  for (int i = 1; i < nx-1; ++i) {
    u_siguiente[i]=-gamma/2.0*(u_anterior[i+1]-u_anterior[i-1])+u_anterior[i];
  }
  u_siguiente[0]=u0;
  u_siguiente[nx-1]=2.0*u_siguiente[nx]-u_siguiente[nx-1];
  
  for (int i = 0; i < nx; ++i) {
    u_anterior[i]=u_siguiente[i];
  }
}

void evolucionLAX(double*u_anterior,int nx,double gamma,double u0)
{
  double u_siguiente[nx];
  for (int i = 1; i < nx-1; ++i) {
    u_siguiente[i]=(u_anterior[i+1]+u_anterior[i-1])/2.0  -gamma/2.0*(u_anterior[i+1]-u_anterior[i-1]);
  }
  u_siguiente[0]=u0;
  u_siguiente[nx-1]=2.0*u_siguiente[nx]-u_siguiente[nx-1];
  
  for (int i = 0; i < nx; ++i) {
    u_anterior[i]=u_siguiente[i];
  }
}
