#include <stdio.h>
#include <stdlib.h>
#include "header.h"


void metropolis2
(double **estado, double **estado_siguiente, int num_moleculas, double* parametros, double* aceptados)
{
  
  double dmax = parametros[1];
  double diametro_cuadrado = parametros[2];
  
  // copiar estado actual con estado siguiente

  copiarVectores(estado_siguiente, estado, num_moleculas, 3);

  double xi, yi, zi;
  double dx, dy, dz;
  double d2;
  
  // el minimo representa la distancia minima al cuadrado
  double min = 10;
  
  for (int i = 0; i < num_moleculas; ++i) {
    min = 10;
    
    xi = estado[i][0] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1);
    yi = estado[i][1] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1);
    zi = estado[i][2] + dmax*( 2*gsl_rng_uniform(generador_uniforme) - 1);

    /* Condicion de Borde Periodica */

    xi = xi - floor(xi);
    yi = yi - floor(yi);
    zi = zi - floor(zi);

    for (int j = 0; j < num_moleculas-1; ++j) {
      if (i != j) {

	/* Distancia mínima entre la particula i, j */
	dx = (xi - estado[j][0]) - floor((xi - estado[j][0])); 
	dy = (yi - estado[j][1]) - floor((yi - estado[j][1]));
	dz = (zi - estado[j][2]) - floor((zi - estado[j][2]));

	// Trabajamos con el cuadrado de la distancia para evitar
	// calcular raices
	d2  = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

	if ( d2 < min ){
	  min = d2;
	}	
      }
    }

    if (min > diametro_cuadrado) {
      estado_siguiente[i][0] = xi;
      estado_siguiente[i][1] = yi;
      estado_siguiente[i][2] = zi;
      *aceptados = *aceptados +1;
    }

      
  }  
}
