#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void rellenarVectorA(hiper p, double *vA,double*x);


// Esta funcion resuelve el problema de la ecuacion diferencial hiperbolica de segundo orden
// con una dimension en el tiempo y una en el espacio
void resolver_metodo(hiper p, FILE* archivo)
{
  //Se definen las variables del problema
  
  double delta_x=(p.xf-p.xi)/(p.nx-1);
  p.lambda=1.0;
  double lambda2=pow(p.lambda,2.0); //lambda cuadrado
  double delta_t=0.5*sqrt(p.lambda)*delta_x/p.D;
  p.nt=(p.tf-p.ti)/(delta_t); //1020 da bien
  
  double vP[p.nx];//Vector Previo
  double vA[p.nx];//Vetor Actual
  double vS[p.nx];//Vector Siguiente
    
  //Grilla;
  double x[p.nx];
  for (int i = 0; i < p.nx; ++i) {
    x[i]=p.xi+delta_x*i;
  }

  //Inicializar matrices
  for (int i = 0; i < p.nx; ++i) {
    vP[i]=0;
    vA[i]=0;
    vS[i]=0;
  }
  
  //condicion inicial t=0
  for (int i = 0; i < p.nx; ++i) {
    vP[i]=p.f(x[i],p.parametros);
  }

  rellenarVectorA(p,vA,x);
  
  printf("%f %f %f %d\n", p.ti,p.tf,delta_t,p.nt);

  double t=0;
  for (int time = 0; time < p.nt; time+=1) {//iteracion de puntos
    t+=delta_t;
    // Puntos a un tiempo siguiente, fuera de los bordes
    for (int i = 1; i < p.nx-1; ++i) {
      vS[i]=2.0*(1.0-lambda2)*vA[i]+lambda2*(vA[i+1]+vA[i-1])-vP[i];
    }
    
    /* Bordes */
  
    //discretizacion de la derivada
    vS[0]=2.0*delta_x/p.beta[0]*( p.alpha[0]*vS[1] - p.gamma[0]) + vS[2];
    
    if(p.beta[0]==0){
      vS[0]=p.gamma[0]/ p.alpha[0];
    }
     
    vS[p.nx-1]=-2.0*delta_x/p.beta[1]*( p.alpha[1]*vS[p.nx-2] - p.gamma[1]) + vS[p.nx-3];

    if(p.beta[1]==0){
      vS[p.nx-1]=p.gamma[1]/p.alpha[1];
    }
    
    // Actualizar
    for (int i = 0; i < p.nx; ++i) {
      vP[i]=vA[i];//el actual ahora es el previo
      vA[i]=vS[i];//el siguiente ahora es el actual
    }

  }
  
  // Imprimir archivo
  for (int i = 0; i < p.nx; ++i) {
    fprintf(archivo, "%E ", x[i]);
    fprintf(archivo, "%E\n",vP[i]);
  }
  fprintf(archivo, "\n");

}

void rellenarVectorA(hiper p, double *vA,double*x)
{
  double delta_x=(p.xf-p.xi)/(p.nx-1);
  p.lambda=1.0;
  double lambda2=pow(p.lambda,2.0); //lambda cuadrado
  double delta_t=0.5*sqrt(p.lambda)*delta_x/p.D;
  /* puntos fantasmas */

  //condicion justo despues de la inicial, t=0+delta_t
  for (int i = 1; i < p.nx-1; ++i) {
    vA[i]=(1.0-lambda2)*p.f(x[i],p.parametros)
      + lambda2/2.0 * ( p.f(x[i+1],p.parametros) + p.f(x[i-1],p.parametros) )
      + delta_t*p.g(x[i],p.parametros);
  }

  //puntos en los bordes de vA

  //Izquierdo
    vA[0]=2.0*delta_x/p.beta[0]*( p.alpha[0]*vA[1] - p.gamma[0]) + vA[2];
  //Dirichlet
  if(p.beta[0]==0.0){
    vA[0]=p.gamma[0] / p.alpha[0];
  }

  //Fantasma Derecho
  vA[p.nx-1]=-2.0*delta_x/p.beta[1] *( p.alpha[1]*vA[p.nx-2] - p.gamma[1]) + vA[p.nx-3];
  //Dirichlet
  if(p.beta[1]==0.0){
    vA[p.nx-1]=p.gamma[1]/ p.alpha[1];
  }
}
				   
