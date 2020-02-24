#include <stdio.h>
#include <stdlib.h>
#include "header.h"


void metropolis
(double **estado, int num_moleculas, double* parametros, double* aceptados)
/* Programa que realiza la simulación de metropolis calculando el siguiente estado.
 
   Entradas: El estado actual
             El número de moleculas
	     El vector de parámetros
	     
   Salidas:  El estado siguiente (se modifica el estado actual)
	     El número de movimientos aceptados (del total de moleculas).
*/
{
  // el desplazamiento máximo de una molecula
  double diametro = parametros[0];
  double dmax = parametros[1];
  

  // posiciones y distancias
  double xi, yi, zi;
  double dx, dy, dz;
  double d;

  // distancia minima entre moleculas, 10 para un valor alto
  double min;
  
  for (int i = 0; i < num_moleculas; ++i) {
    
    min = 10.0;

    /* Tomamos una posicion al azar */
    xi = estado[i][0] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1.0);
    yi = estado[i][1] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1.0);
    zi = estado[i][2] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1.0);

    /* Condicion de Borde Periodica */
    xi = xi - floor(xi);
    yi = yi - floor(yi);
    zi = zi - floor(zi);

    for (int j = 0; j < num_moleculas; ++j) {
      if (i != j) {

	/* Distancia mínima entre la particula i y la imagen más cercana de j */

	dx = (xi - estado[j][0]);
	dx = dx - round(dx);

	dy = (yi - estado[j][1]);
	dy = dy - round(dy);
	
	dz = (zi - estado[j][2]);
	dz = dz - round(dz);
	
	d  = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

	if ( d < min ){
	  min = d;
	}
		
      }
    }

    // Aceptamos el estado si se cumple la condicion
    if (min > diametro) {
      estado[i][0] = xi;
      estado[i][1] = yi;
      estado[i][2] = zi;
      *aceptados = *aceptados +1;
    }

  }  
  
}
