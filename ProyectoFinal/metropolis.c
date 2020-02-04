#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/* TODO: Revisar esta seccion de código*/


double probabilidadSiguiente(double **estado, double **estado_actual, int num_moleculas, double* parametros)
{
  double prob = 1;
  
  //diferencia de las posiciones en x ,y ,z
  double dif_posiciones[3];
  double dif_cuadrado;

  double diametro_cuadrado = parametros[2];
  double aux = 0;//valor auxiliar
  
  for (int i = 0; i < num_moleculas-1; ++i) {
    for (int j = i+1; j < num_moleculas; ++j) {
      
      /* TODO: Condiciones periodicas de borde */
      dif_posiciones[0] = estado[i][0] - estado[j][0];// - largo_cajaX*round( (estado[i][0] - estado[j][0])/largo_cajaX )
      dif_posiciones[1] = estado[i][1] - estado[j][1];
      dif_posiciones[2] = estado[i][2] - estado[j][2];
      dif_cuadrado = pow(dif_posiciones[0],2) + pow(dif_posiciones[1],2) + pow(dif_posiciones[2],2);

      // Si almenos dos pelotitas solapan, la probabilidad es inmediatamente 0 y el sistema no se permite
      if(dif_cuadrado < diametro_cuadrado){
	aux = 1;
	break;
      }
    }
    if(aux == 1){
      break;
    }
  }

  if(aux == 1){
    prob = 0;
  }

  return prob;
}


// Esta funcion contiene el algoritmo de Metropolis
void metropolis
(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros, double* prob_actual)
{ 
  // Volumen=1
  
  double dmax = parametros[1];
  double probabilidad_aceptacion;
  double prob_siguiente;
  double perturbacion;
  double random;
  int aux = 0;
  int hola = 0;

  
  while(aux == 0){
    ++hola;
    printf("%d\n", hola);
    for (int i = 0; i < num_moleculas; ++i) {
      for (int j = 0; j < 3; ++j) {

	// El número random estará entre -1 y 1
	perturbacion = 2*gsl_rng_uniform(generador_uniforme)-1;

	//Posible estado siguiente
	estado_siguiente[i][j] = estado_actual[i][j] + dmax*perturbacion;
      }
    }
    
    prob_siguiente = probabilidadSiguiente(estado_siguiente, estado_actual, num_moleculas, parametros);
    
    probabilidad_aceptacion = prob_siguiente/ *prob_actual;

    // Aceptamos estado siguiente si pA =!
    if( probabilidad_aceptacion == 1.0 ){
      aux = 1;
    }
    else{
      // Tomamos numero aleatorio entre 0 y 1
      random = gsl_rng_uniform(generador_uniforme);

      // Aceptamos estado siguiente si pA >= r
      if( probabilidad_aceptacion >= random){
	aux = 1;
      }
    }  
  }
  *prob_actual = prob_siguiente;
}
