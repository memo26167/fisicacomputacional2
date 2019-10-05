
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
 * Para y=0,  α1Φ + β1*g1(x)* ∂Φ/∂y = γ1*f1(x)
 * Para x=0,  α2Φ + β2*g2(y)* ∂Φ/∂x = γ2*f2(y)
 * Para y=Ly, α3Φ + β3*g3(x)* ∂Φ/∂y = γ3*f3(x)
 * Para x=Lx, α4Φ + β4*g4(y)* ∂Φ/∂x = γ4*f4(y)
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
 * Donde está γi, enrealidad debe estar γi*fi
 */


#include <stdio.h>
#include <stdlib.h>
#include "header.h"

int main(void)
{
  ecp problema;//se define el problem
  problema.nx=100;
  problema.ny=100;
  problema.tol=1e-4;
  problema.D=1;
  problema.m_inicial= malloc(sizeof(double*)*problema.nx);
  for (int i = 0; i < problema.nx; ++i) {
    problema.m_inicial[i]= malloc(sizeof(double)*problema.ny);
  }
  
  for (int i = 0; i < problema.nx; ++i) {
    for (int j = 0; j < problema.ny; ++j) {
      problema.m_inicial[i][j]=28;
    }
  }
  
  return 0;
}



