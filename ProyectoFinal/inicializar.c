/* Programa que inicializa el sistme en una red HEXAGONAL,
   Entradas: estado del sistema, diametro de pelotas, 
   Salidas:  estado
   Consideraciones: Volumen del sistema = 1*/

#include "header.h"

void inicializar(double ** estado, int num_moleculas, double* parametros)
{
  int imax = ceil( pow(num_moleculas,1.0/3.0));
  int jmax = ceil( pow(num_moleculas,1.0/3.0));
  int kmax = floor( pow(num_moleculas,1.0/3.0));
  /* printf("%d", imax); */
  /* printf("%d", jmax); */
  /* printf("%d", kmax); */
  if(imax*jmax*kmax < num_moleculas){
    kmax = ceil( pow(num_moleculas,1.0/3.0));
  }
  
  // n es el índice de la molecula
  int n = 0;

  // R es el radio de una molecula
  double R = parametros[3];

  // auxiliar para romper for
  
  int aux = 0;
  for (int i = 0; i < imax; ++i) {
    for (int j = 0; j < jmax; ++j) {
      for (int k = 0; k < kmax; ++k) {
        estado[n][0] = R*(1.0 + 2.0*sqrt(6.0)/3.0*k);//z;
	estado[n][1] = R*(1.0 + sqrt(3.0)*(j + k%2/3) );//x;
	estado[n][2] = R*(2*(i+1) - (j + k)% 2);//y;

	++n;
	
	if(n == num_moleculas){
	  aux = 1;
	}

	//rompemos primero for
	if(aux == 1){
	  break;
	}	
      }

      //rompemos segundo for
      if(aux == 1){
	break;
      } 
    }

    //rompemos tercer for
    if(aux == 1){
      break;
    }    
  }
  
  
}
