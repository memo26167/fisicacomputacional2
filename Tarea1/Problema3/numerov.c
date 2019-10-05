#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void diferencial(double *x, double*y,double energia,int N, double paso_h)
{
  double aux1=0.0;
  double aux2=0.0;
  double aux3=0.0;
  for (int i = 1; i < N+1; ++i) {
    x[i+1]=x[i]+paso_h;
    aux3=( 1.0+1.0/(6.0*N*N)*(energia-potencial(x[i+1])) );
    aux1=2.0*(1.0-5.0/(6.0*N*N)*(energia-potencial(x[i])) )*y[i];
    aux2=-(1.0+1.0/(6.0*N*N)*(energia-potencial(x[i-1])))*y[i-1];
    y[i+1]=(aux1+aux2)/aux3; 
  }
  double A=0;
  for (int i = 0; i < N+1; ++i) {
    A=A+y[i];
  }
  if(A==0){
    A=1;
  }
  for (int i = 0; i < N+1; ++i) {
    y[i]=y[i]/fabs(A);
  }
  //normalizacion
  for (int i = 0; i < N+1; ++i) {
    A=A+pow(y[i],2.0);
  }
  if(A==0){
    A=1;
  }
  for (int i = 0; i < N+1; ++i) {
    y[i]=y[i]/sqrt(fabs(A));
  }
  
}

void numerov(int numero_autofunciones,double xa,double ya0,double ya1,int N,double paso_h,double nmax,double tol,FILE*archivo)
{
  
  double En[numero_autofunciones];
  double y[N+1];
  double x[N+1];
  double ea=1.0;
  double eb=2.0;
  double ec=3.0;
  double phia=0;
  double phib=0;
  double phic=0;
  double iteraciones=0;
  
  for (int i = 0; i <numero_autofunciones; ++i) {
    x[0]=xa;  
    x[1]=x[0]+paso_h;
    y[0]=ya0;
    y[1]=ya1;
    diferencial(x,y,ea,N,paso_h);
    phia=y[N];
    iteraciones=0;

    while(phia*phib>=0 && iteraciones<=nmax){
      iteraciones=iteraciones+1;
      eb=eb+1;
      diferencial(x,y,eb,N,paso_h);
      phib=y[N];
    }
    iteraciones=0;
    


    while(fabs(y[N])>tol && iteraciones <nmax ) {
      iteraciones=iteraciones+1;
      ec=(ea+eb)/2.0;
      diferencial(x,y,ec,N,paso_h);
      phic=y[N];
      if(phia*phic<0){
	eb=ec;
	phib=phic;
      }
      else{
	phia=phic;
	ea=ec;
      }
    }
 
    iteraciones=0;
    En[i]=ec;
    //printf("%f\n", En[i]);
    phib=0;

    fprintf(archivo,"\"En=%f\" \n",En[i]);
         for (int i = 0; i < N+1; ++i) {
      fprintf(archivo, "%E %E", x[i],y[i]);
      fprintf(archivo,"\n");
    }
    fprintf(archivo, "\n\n");
    ea=En[i]+1;
    eb=ea+1;
  }  
}
