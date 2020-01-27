#include <stdio.h>
#include <stdlib.h>
#include "header.h"


double funcionDistribucion(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
// La distribucion de probabilidad no normalizada, de que ocurra un estado siguiente dado el estado actual
{
  double p = 0;
  double diferencia_energia = funcionEnergia(estado_siguiente, num_moleculas, parametros)
                              - funcionEnergia(estado_actual, num_moleculas, parametros);
  // parametros[0] es beta = 1 / (kB * T)
  p = exp(-parametros[0]*diferencia_energia);
  return p;
}


double probabilidadAceptacion(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
{
  double prob = funcionDistribucion(estado_actual,estado_siguiente,num_moleculas,parametros);
  if (prob >= 1.0){
    prob = 1.0;
  }
  return prob;
}


// Esta funcion contiene el algoritmo de Metropolis
void siguienteEstado(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
{ 
  // parametros[1] es dmax
  // Volumen=1
  
  double dmax = parametros[1];
  double probabilidad;
  double random;
  //Posible estado siguiente

  while(1){
    for (int i = 0; i < num_moleculas; ++i) {
      for (int j = 0; j < 3; ++j) {
	random = 2*gsl_rng_uniform(generador_uniforme)-1;
	estado_siguiente[i][j] = estado_actual[i][j] + dmax*random;
      }
    }
    probabilidad = probabilidadAceptacion(estado_actual,estado_siguiente,num_moleculas,parametros); 

    if( probabilidad == 1.0 ){
      break; //aceptamos el nuevo estado
    }
    else{
      if( probabilidad >= gsl_rng_uniform(generador_uniforme)){
	break; //aceptamos el nuevo estado
      }
    }
    
  }  
}
