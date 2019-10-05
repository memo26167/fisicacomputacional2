#include <stdio.h>
#include <stdlib.h>
#include "header.h"
/* Este programa resuelve la ecuacion diferencial parcial eliptica
 * en geometria cartesiana bidimensional,que es de la forma 
 * 
 * ∇²Φ(x,y)=∂²Φ/∂²x+∂²Φ/∂²y=-S(x,y)
 *
 * Donde S(x,y) es una fuente
 *
 * Se genera una cuadricula de puntos, con separacion Δx e Δy, 
 * en donde, en los puntos alejados de los bordes
 * se utiliza diferencia central para aproximar la derivada parcial
 * 
 * ∂²Φ/∂²x= (Φ[i-1,j] -2Φ[i,j] + Φ[i+1,j])/(Δx)²
 * 
 * y la ecuacion queda:
 * 
 * (Φ[i-1,j] -2Φ[i,j] + Φ[i+1,j])/(Δx)² + (Φ[i,j-1] -2Φ[i,j] + Φ[i,j+1])/(Δy)² = S[i,j]
 *
 * Para los puntos en bordes, se utiliza la condicion de borde mixta:
 * 
 * Para y=0,  ∂Φ/∂y + α0Φ= β0* f0(x)
 * Para x=0,  ∂Φ/∂x + α1Φ= β1* f1(y)
 * Para y=Ly, ∂Φ/∂y + α2Φ= β2* f2(x)
 * Para x=Lx, ∂Φ/∂x + α3Φ= β3* f3(y)
 *
 * Si α!=0{
 *  γ=β/α;
 * }
 * Sino{
 * γ=0;
 * }
 *
 * Esto permite escribir la ecuacion de arriba de la forma
 * (Φ[i-1,0] -2Φ[i,0] + Φ[i+1,0])/(Δx)² + 2(Φ[i,1] -Φ[i,0])/(Δy)² + 2αΦ[i,0]/(Δy)= S[i,j] + 2γαf(x)/Δy
 *
 * Si α =0,        entonces
 * ∂Φ/∂y = β * f(x)
 *
 * Si α > 1E+200, entonces
 * Φ = β/α * f(x) = γ*f(x)
 *
 * Por lo que se debe entregar dos vectores, α y β de condiciones de borde
 * Primer componente,  borde inferior de y, y=0
 * Segundo componente, borde inferior de x, x=0
 * Tercer componente,  borde superior de y, y=Ly
 * Segundo componente, borde superior de x, x=Lx
 *
 * Si se da el segundo caso, β debe ser calculado de manera que
 * β=γ*α
 */

double f0(double x)
{
  double salida;
  salida=x*(10-x);
  return salida;
}

double f1(double y)
{
  double salida;
  salida=y*(10-y);
  return salida;
}

double f2(double x)
{
  double salida;
  salida=0*x;
  return salida;
}

double f3(double y)
{
  double salida;
  salida=0*y;
  return salida;
}

double funcion_fuente(double x, double y)
{
  double salida;
  salida=100*exp(-(pow(x-5,2)/2+pow(y-5,2)/2))/1.04;
  return salida;
}

int main(void)
{
  double vector_alpha[4]={1,1,1,1};
  double vector_beta[4]={1,1,1,1};
  pEl problema;

  problema.n=50;
  problema.m=50;
  problema.S=&funcion_fuente;
  problema.f0=&f0;
  problema.f1=&f1;
  problema.f2=&f2;
  problema.f3=&f3;
  problema.alpha=vector_alpha;
  problema.beta=vector_beta;
  // xa <= x <= xb
  problema.xa=0; 
  problema.xb=10;
  // ya <= y <= yb
  problema.ya=0;
  problema.yb=10;
  
  FILE *archivo = fopen("salida.dat", "w");
  diferenciasFinitas(problema,archivo);
  fclose(archivo);
  
  return 0;
}

