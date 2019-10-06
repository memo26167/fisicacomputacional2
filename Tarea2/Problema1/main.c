
/*  Se resolverá el problema de ecuacion diferencial paracial Parabolica
 * dependiente del tiempo.
 * Se tiene la ecuacion
 * 
 * dΦ(x,y,t)/dt = D∇²Φ(x,y,t) +S(x,y)
 *
 * dΦ(x,y,t)/dt = D(∂²Φ/∂²x + ∂²Φ/∂²y) +S(x,y)
 *
 * Donde S(x,y) es una fuente
 *
 * Se genera una cuadricula de puntos, con separacion Δx ligado a i, Δy ligado a j y Δt ligado a n (sin cuadricula) 
 * en donde, en los puntos alejados de los bordes
 * se utiliza diferencia central para aproximar la derivada parcial
 * 
 * y la ecuacion queda:
 * 
 * Φ[n+1,i,j]=Φ[n,i,j]+D*Δt*(Φ[n,i-1,j] -2Φ[n,i,j] + Φ[n,i+1,j])/(Δx)² + (Φ[n,i,j-1] -2Φ[n,i,j] + Φ[n,i,j+1])/(Δy)² + S[i,j]*Δt
 *
 * si se define Cx=d*Δt/(Δx)² y Cy=d*Δt/(Δy)²,la ecuacion se puede reescribir
 *
 * Φ[n+1,i,j] = (1-2Cx-2Cy)*Φ[n,i,j] + Cx*(Φ[n,i-1,j] + Φ[n,i+1,j]) + Cy*(Φ[n,i,j-1] + Φ[n,i,j+1]) + S[i,j]*Δt
 * 
 * Para los puntos en bordes, se utiliza la condicion de borde mixta:
 * 
 * Para x=0,  α1Φ + β1*g1(y)* ∂Φ/∂x = γ1*f1(y) BORDE IZQUIERDO
 * Para x=Lx, α2Φ + β2*g2(y)* ∂Φ/∂x = γ2*f2(y) BORDE DERECHO
 * Para y=0,  α3Φ + β3*g3(x)* ∂Φ/∂y = γ3*f3(x) BORDE INFERIOR
 * Para y=Ly, α4Φ + β4*g4(x)* ∂Φ/∂y = γ4*f4(x) BORDE SUPERIOR
 *
 * Si β==0{
 *  β=1e-300
 * }
 * que permite escribir la condicion de tipo Dirichlet aproximada sin realizar grandes desajustes
 * 
 * Esto permite escribir la ecuacion de arriba para cada borde o esquina considerado, la insercion de puntos fantasmas (no existentes en la cuadricula):
 *
 * Φ[0,j]   = Φ[2,j]    -2Δx /β1 * (γ1 - α1*Φ[1,j] ) 
 * Φ[Nx+1,j]= Φ[Nx-1,j] -2Δx /β2 * (γ2 - α2*Φ[Nx,j])
 * Φ[i,0]   = Φ[i,2]    -2Δy /β3 * (γ3 - α3*Φ[i,1] )
 * Φ[i,Ny+1]= Φ[i,Ny-1] -2Δy /β4 * (γ4 - α4*Φ[i,Ny])
 *
 * Donde está βi, enrealidad debe estar βi*gi
 * Donde está γi, enrealidad debe estar γi*fi, 
 *  sin embargo, fi integrará γi
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"

// Funcion fuente de la EDP
double fuente(double x,double y ,double* parametros);
// Funciones de las condiciones de borde
double f1(double y, double* parametros);
double f2(double y, double* parametros);
double f3(double x, double* parametros);
double f4(double x, double* parametros);
double g1(double y, double* parametros);
double g2(double y, double* parametros);
double g3(double x, double* parametros);
double g4(double x, double* parametros);

int main(void)
{
  ecp problema;//se define el problem
  /*Parametros del problema*/

  //constantes del problema a resolver, depende del contexto
  problema.parametros= malloc(sizeof(double)*4);
  problema.parametros[0]=28; //T infinito
  problema.parametros[1]=5; //h 
  problema.parametros[2]=401; //k
  problema.parametros[3]=113e-6; 
  //constante de la EDP
  problema.D=1;
  
  //parametros del metodo
  problema.nx=100;
  problema.ny=100;
  problema.xi=0;
  problema.xf=1;
  problema.yi=0;
  problema.yf=1;
  //condicion inicial
  problema.m_inicial= malloc(sizeof(double*)*problema.nx);
  for (int i = 0; i < problema.nx; ++i) {
    problema.m_inicial[i]= malloc(sizeof(double)*problema.ny);
  }
  for (int i = 0; i < problema.nx; ++i) {
    for (int j = 0; j < problema.ny; ++j) {
      problema.m_inicial[i][j]=28;
    }
  }
  
  //funcion fuente de la EDP
  problema.fuente=&fuente;

  /* Condiciones de borde */
  problema.alpha[0]=1;
  problema.alpha[1]=1;
  problema.alpha[2]=0;
  problema.alpha[3]=0;

  problema.beta[0]=0;
  problema.beta[1]=-problema.parametros[2];
  problema.beta[2]=1;
  problema.beta[3]=1;
  
  problema.f1=&f1;
  problema.f2=&f2;
  problema.f3=&f3;
  problema.f4=&f4;
  problema.g1=&g1;
  problema.g2=&g2;
  problema.g3=&g3;
  problema.g4=&g4;

  resolver_metodo(problema);
    
  return 0;
}

double fuente(double x,double y ,double* parametros){
  double S;
  S=x*y*parametros[0]*0;
  return S;
}

double f1(double y, double* parametros){
  y=y*parametros[0]*0+1;
  return y;
}

double f2(double y, double* parametros){
  y=y*parametros[0]*0+1;
  return y;
}
double f3(double x, double* parametros){
  x=x*parametros[0]*0+1;
  return x;
}
double f4(double x, double* parametros){
  x=x*parametros[0]*0+1;
  return x;
}
double g1(double y, double* parametros){
  y=y*parametros[0]*0+1;
  return y;
}
double g2(double y, double* parametros){
  y=y*parametros[0]*0+1;
  return y;
}
double g3(double x, double* parametros){
  x=x*parametros[0]*0+1;
  return x;
}
double g4(double x, double* parametros){
  x=x*parametros[0]*0+1;
  return x;
}

