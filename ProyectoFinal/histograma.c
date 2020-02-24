/* Archivo que se utiliza para construir la distribucion g(r)
   en base a histogramas
   
   Se tienen dos funciones.
   
   void histograma, cuya funcion es calcular el histograma de un sistema i_hist 

   void analizarHistograma, cuya funcion es promediar los histogramas y escalarlos 
                    en base al inverso de la distancia al cuadrado mas algunos parametros
*/

#include <stdlib.h>
#include <stdio.h>
#include "header.h"


void histograma
(double ** estado, int num_moleculas, double* hist, int num_hist, double* parametros)
/* Este programa construye el histograma de un sistema de particulas (estado)
   Entradas: El estado de las particulas
             El numero de moleculas
	     El vector (matriz) histograma
	     El numero de particiones del histograma
	     El vector de parámetros
   Salidas: El histograma
*/
{
  double diametro = parametros[0];

  // El rango final del histograma será tres veces el diametro
  double rango_final = diametro*3.0;
  double rango_inicial = diametro*1.0;
  // La distancia diferencial con la cual se construira el histograma
  double diff_hist = (rango_final - rango_inicial)/( (double) num_hist - 1.0);
  
  // el índice que correspondería para una cierta distancia
  double dist_k;
  int k;

  // las distancias entre particulas
  double dx, dy, dz;
  double d;
  
  for (int i = 0; i < num_moleculas-1; ++i) {
    for (int j = i+1; j < num_moleculas; ++j) {

      // Analizamos las imágenes de las moleculas de las celdas contiguas (no necesariamente las mas cercanas)      
      for (double l = -1; l <= 1; ++l) {
	for (double m = -1; m <= 1; ++m) {
	  for (double n = -1; n <= 1; ++n) {

	    dx = (estado[i][0] - estado[j][0]) + l;
            
	    dy = (estado[i][1] - estado[j][1]) + m;
	    
	    dz = (estado[i][2] - estado[j][2]) + n;
	    
	    d = sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));

	    // Calculamos el índice en el cual debemos incrementar el histograma 
	    dist_k = (d-rango_inicial)/diff_hist;
	    k = floor(dist_k);
	    
	    if ( 0 <= k && k < num_hist ){
	      hist[k] = hist[k] + 1.0;
	    }

	    
	  }
	}
      }      
    }
  }
}


void analizarHistograma(double * hist, int num_hist, int num_moleculas, double** prom_hist, double* parametros)
/* Este programa construye la distribucion g(r) de los sistemas de particulas
   Entradas: Los histogramas
             El numero de particiones del histograma
	     El número de moleculas
	     El vector de parámetros

   Salidas:  El histograma promedio
*/
{
  double diametro = parametros[0];

  // El rango final del histograma será tres veces el diametro
  double rango_final = diametro*3.0;
  double rango_inicial = diametro*1.0;
  // La distancia diferencial con la cual se construira el histograma
  double diff_hist = (rango_final - rango_inicial)/( (double) num_hist - 1.0);

  // Iteramos en las particiones del histograma
  for (int i = 0; i < num_hist; ++i) {
    
    prom_hist[i][0] = ((double)i)*diff_hist + rango_inicial;

    prom_hist[i][1] = hist[i] + prom_hist[i][1];
    
    /* Se promedia la acumulacion y se multiplica por dos, 
       ya que la acumulacion sólo considera 1 vez un par de moleculas 
       (se debe considerar 2 veces).
       Recordar que parametros[2] es el número de histogramas aceptadas evitando el estado no termalizado */
    prom_hist[i][1] = 2.0 *prom_hist[i][1] /(parametros[2]);

    // Se arregla la métrica de la distribucion
    prom_hist[i][1] = prom_hist[i][1] /(pow(prom_hist[i][0], 2.0));

    // se normaliza la distribucion
    prom_hist[i][1] =  prom_hist[i][1] /pow(num_moleculas, 2.0)
                       /(4.0*3.14159265359*diff_hist);
    
    // Se normaliza la distancia
    prom_hist[i][0] = prom_hist[i][0] / diametro ;
  }
}
