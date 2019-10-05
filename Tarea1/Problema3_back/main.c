#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* Este programa resuelve el problema de frontera de segundo orden
 * y''=p(x)y'+q(x)y+r(x)
 * s.t.  y(a)=alpha
 *       y(b)=beta
 *       a<=x<=b
 * Mediante la resolucion de dos problemas de valor inicial junto a RK4
 * 
 * y1''=p(x)y1'+q(x)y1+r(x)     y2''=p(x)y2'+q(x)y2
 * s.t.  y1(a)=alpha            s.t.  y2(a)=0
 *       y1'(a)=0                     y2'(a)=1
 *       a<=x<=b                      a<=x<=b
 *
 * Para luego obtener la solucion general del problema de frontera
 * con la ecuacion
 * 
 * y(x)=y1(x)+(beta-y1(b))/y2(b)*y2(x)
 */

double energia=0;// metodo KISS

double funcionP(double x)
{
  double p;
  p=0*x;
  return p;
}

double funcionQ(double x)
{
  double q;
  double potencialV=50*exp(-pow(x-0.5,2)/0.08);//50*exp(-pow(x-0.5,2.0)/0.08)*0;//LINEA CAMBIADA
  q=0*x-2*(energia-potencialV);
  return q;
}

double funcionR(double x)
{
  double r;
  r=0*x+1e-300;
  return r;
}


int main(void)
{
  // Valores del intervalo de a<x<b
  double intval_a=0;
  double intval_b=1;
  
  double ya=0;
  double yb=0;
  double cond[2]={ya,yb};// condiciones de frontera
  pRKF problema_frontera;
  problema_frontera.p=&funcionP;
  problema_frontera.q=&funcionQ;
  problema_frontera.r=&funcionR;
  double paso_h=0.01;

  //ALGORITMO DE NUMEROV
  //energias y soluciones iniciales
  int numero_autofunciones=5;
  double ea=0;
  double eb=0;
  double ec;
  double phia=0;
  double phib=0;
  double phic=0;
  double En[numero_autofunciones];
  double iteraciones=0;
  double nmax=5000;
  double tol=1e-10;


  // Se redondea hacia arriba para obtener la mayor cantidad de componentes
  int num_it=ceil((intval_b-intval_a)/paso_h);
  // Se recalcula el paso para asegurar que se obtenga el tiempo final
  paso_h=(intval_b-intval_a)/num_it;
  ++num_it;
  
  //Se define un vector para almacenar la solucion para el problema de frontera
  gsl_vector *solucion_frontera=gsl_vector_alloc(num_it);
  gsl_vector_set_zero(solucion_frontera);
  // Vector con contenido del valor de la variable independiente de la iteracion
  // que da solucion al problema de el paso adaptativo, puede representar el tiempo
  gsl_vector *sol_x=gsl_vector_alloc(num_it);
  gsl_vector_set_zero(sol_x);
  gsl_vector_set(sol_x,0,intval_a);

  FILE *archivo = fopen("salida.dat", "w");

  for (int i = 0; i <numero_autofunciones; ++i) {  
    energia=ea;
    shooting(cond,paso_h, problema_frontera,sol_x,solucion_frontera,num_it,archivo);
    phia=gsl_vector_get(solucion_frontera,num_it-4);
    phia=copysign(1.0,phia);
    iteraciones=0;

    while(phia*phib>=0 && iteraciones<nmax){
      iteraciones=iteraciones+1;
      eb=eb+0.1; 
      energia=eb; 
      shooting(cond,paso_h, problema_frontera,sol_x,solucion_frontera,num_it,archivo);
      phib=gsl_vector_get(solucion_frontera,num_it-4);
      phib=copysign(1.0,phib);
    }
    iteraciones=0;
    ++eb;
    
    while(fabs(ea-eb)>tol && iteraciones <nmax ){
      iteraciones=iteraciones+1;
      ec=(ea+eb)/2;
      energia=ec;
      shooting(cond,paso_h, problema_frontera,sol_x,solucion_frontera,num_it,archivo);
      phic=gsl_vector_get(solucion_frontera,num_it-4);
      phic=copysign(1.0,phic);
      if(phia*phic<0){
    	phib=phic;
    	eb=ec;
      }
      else{
    	phia=phic;
    	ea=ec;
      }
    }
    En[i]=ec;
    printf("%f\n", En[i]);
    phib=0;
    fprintf(archivo,"\"En=%f\" \n",En[i]);
    normalizar(solucion_frontera,num_it);
    
    for (int i = 0; i < num_it; ++i) {
      fprintf(archivo, "%E ", gsl_vector_get(sol_x,i));
      fprintf(archivo, "%E", gsl_vector_get(solucion_frontera,i));
      fprintf(archivo,"\n");
    }
    fprintf(archivo, "\n\n");
    ++ea;
  }

  fclose(archivo);
  return 0;
}
