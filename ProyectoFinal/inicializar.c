
#include <stdlib.h>
#include <stdio.h>
#include "header.h"

void inicializar(double ** estado, int num_moleculas, double* parametros)
/* Programa que inicializa el sistme en una red FCC,
   Entradas: estado del sistema, diametro de pelotas, 
   Salidas:  estado
   Consideraciones: Volumen del sistema = 1*/
{
  // R es el radio de una molecula
  double R = parametros[0]/2.0;

  // Tamaño de la celda unitaria
  double a = 2.0*sqrt(2.0)*R;
  
  // Celda unitaria FCC
  double U[4][3]={
    { 0, 0, 0},
    { a/2, a/2, 0},
    { 0, a/2, a/2},
    { a/2, 0, a/2}};

  
  // Cada celda unitaria tiene 4 átomos
  // Calculo de celdas unitarias por cada eje
  double cpj = pow(num_moleculas/4.0, 1.0/3.0);

  int imax = floor(cpj);
  int jmax = imax;
  int kmax = imax;
  if ( imax*jmax*kmax*4 < num_moleculas ) { imax = ceil(cpj); }
  if ( imax*jmax*kmax*4 < num_moleculas ) { jmax = ceil(cpj); }
  if ( imax*jmax*kmax*4 < num_moleculas ) { kmax = ceil(cpj); }

  // indice del número de moleculas
  int n = 0;

  // la coordenada espacial de una celda a probar
  double x = 0;
  double y = 0;
  double z = 0;

  for (int k = 0; k < kmax; ++k) {
    for (int j = 0; j < jmax; ++j) {
      for (int i = 0; i < imax; ++i) {

	// calculo de la ubicacion de la posible siguiente celda
	x = a*(double)i;
	y = a*(double)j;
	z = a*(double)k;

	// iteramos en la celda unitaria U
	for (int l = 0;  l < 4; ++ l) {
	  /* si la ubicacion de la molecula de la celda unitaria 
	     esta dentro del volumen total, aceptamos*/

	  if( x + U[l][0]< 1 ){ // - 2*R ){
	    if( y + U[l][1]< 1 ){
	      if( z + U[l][2]< 1 ){
		
		if (n < num_moleculas) {
		  
		  
		  estado[n][0] = x + U[l][0];
		  estado[n][1] = y + U[l][1];
		  estado[n][2] = z + U[l][2];
		  ++n;
		  
		  //cierre de if
		}}}}
	}
      }
    }
  }
}
