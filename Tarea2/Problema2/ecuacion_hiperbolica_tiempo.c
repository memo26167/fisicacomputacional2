#include <stdio.h>
#include <stdlib.h>
#include "header.h"

// Esta funcion resuelve el problema de la ecuacion diferencial hiperbolica de segundo orden
// con una dimension en el tiempo y una en el espacio
void resolver_metodo(hiper p, FILE* archivo)
{
  //Se definen las variables del problema
  
  double delta_x=(p.xf-p.xi)/(p.nx-1);
  p.lambda=0.1;
  double delta_t=p.lambda*delta_x/p.D;
  p.nt=(p.tf-p.ti)/(delta_t)+1;
  int count_t=0;
  
  double vP[p.nx];//Vector Previo
  double vA[p.nx];//Vetor Actual
  double vS[p.nx];//Vector Siguiente
  
  //matriz para evolucionar sistema
  double ** A;
  A=malloc(sizeof(double*)*p.nx);
  for (int i = 0; i < p.nx; ++i) {
    A[i]= malloc(sizeof(double)*p.nx);
  }
  rellenar_A(A,p);
  //variable para poder realizar la multiplicacion Aij*vAj
  double acum;
  
  //Grilla;
  double x[p.nx];
  for (int i = 0; i < p.nx; ++i) {
    x[i]=p.xi+delta_x*i;
  }
  
  //condicion inicial t=0
  for (int i = 0; i < p.nx; ++i) {
    vP[i]=p.f(x[i],p.parametros);
  }

  //condicion justo despues de la inicial, t=0+delta_t
  for (int i = 0; i < p.nx; ++i) {
    acum=0;
    for (int j = 0; j < p.nx; ++j) {
      acum=acum+A[i][j]*vP[j];
    }
    vA[i]=2*acum+4*delta_t*p.g(x[i],p.parametros);
  }

  
  printf("%f %f %f %d\n", p.ti,p.tf,delta_t,p.nt);

  for (double t = p.ti; t < p.tf; t+=delta_t , count_t+=1 ) {//iteracion de puntos
    // Puntos a un tiempo siguiente, fuera de los bordes
    for (int i = 1; i < p.nx-1; ++i) {
      acum=0;
      for (int j = 0; j < p.nx; ++j) {
	acum=acum+A[i][j]*vA[j];
      }
      vS[i]=acum-vP[i];      
    }
    
    /* Bordes */

    /* puntos fantasmas */
    double v0=0;
    double vN=0;

    //discretizacion de la derivada
    v0=2*delta_x/(p.beta[1]*p.beta1(t,p.parametros))
      *( p.alpha[1]*p.alpha1(t,p.parametros)*vA[0] - p.gamma[1]*p.gamma1(t,p.parametros))
      + vA[1];

    if(p.beta[1]==0){
      v0=p.gamma[1]*p.gamma1(t,p.parametros) / (p.alpha[1]*p.alpha1(t,p.parametros) );
    }
     
    vN=-2*delta_x/(p.beta[2]*p.beta2(t,p.parametros))
      *( p.alpha[2]*p.alpha2(t,p.parametros)*vA[p.nx-1] - p.gamma[2]*p.gamma2(t,p.parametros))
      + vA[p.nx-2];

    if(p.beta[2]==0){
      vN=p.gamma[2]*p.gamma2(t,p.parametros) / (p.alpha[2]*p.alpha2(t,p.parametros) );
    }
    
    //Borde izquierdo i=0
    vS[0]= A[0][0]*vA[0] + A[0][1]*( vA[1] - v0 ) - vP[0];

    //Borde derecho i=n
    vS[p.nx-1]= A[p.nx-1][p.nx-1]*vA[p.nx-1] + A[p.nx-1][p.nx-2]*( vN- vA[p.nx-2] ) - vP[p.nx-1];

    // Imprimir archivo
    for (int i = 0; i < p.nx; ++i) {
      fprintf(archivo, "%E ", x[i]);
      fprintf(archivo, "%E\n",vS[i]);
    }
    fprintf(archivo, "\n");
 
    // Actualizar
    for (int i = 0; i < p.nx; ++i) {
      vP[i]=vA[i];//el actual ahora es el previo
      vA[i]=vS[i];//el siguiente ahora es el actual
      vS[i]=0;//el nuevo siguiente está por calcularse
    }

  }
  free(A);
}


// Esta funcion rellena la matriz A utilizada por la funcion resolver_metodo
void rellenar_A(double **A, hiper p)
{
  //inicializamos
  for (int i = 0; i < p.nx; ++i) {
    for (int j = 0; j < p.nx; ++j) {
      A[i][j]=0;
    } 
  }
  
  for (int i = 1; i < p.nx-1; ++i) {
    A[i][i]=2*(1-pow(p.lambda,2));
    A[i][i-1]=-pow(p.lambda,2);
    A[i][i+1]=pow(p.lambda,2);
  }
  
  A[0][0]=2*(1-pow(p.lambda,2));
  A[0][1]=pow(p.lambda,2);

  A[p.nx-1][p.nx-1]=2*(1-pow(p.lambda,2));
  A[p.nx-1][p.nx-2]=-pow(p.lambda,2);  
}
