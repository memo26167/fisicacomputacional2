 #include <stdio.h>
#include <stdlib.h>
#include "header.h"

void imprimirEstado(double ** estado,int num_moleculas, FILE* archivo, double* parametros)
/* Este programa imprime el estado final de la simulacion del algoritmo de metropolis 
   En formato XYZ para que lo acepte OVITO
*/
{
  double radio = parametros[0]/2;
  //*** CONFIGURACION DE ARCHIVO XYZ ***

  //numero de particulas
  fprintf(archivo, "%d\n", num_moleculas);

  //tamaño de caja
  fprintf(archivo,"Lattice=\"1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0\" ");

  fprintf(archivo,"Properties=id:I:1:pos:R:3:radius:R:1 \n");
  for(int i = 0; i < num_moleculas; ++i){
    //id de la molecula
    fprintf(archivo,"%d ", i);
    
    //coordenada x y z de la molecula
    fprintf(archivo,"%E %E %E ", estado[i][0], estado[i][1], estado[i][2]);
    
    //Radio
    fprintf(archivo,"%E ", radio);
    
    //final de linea
    fprintf(archivo, "\n");
  }
  
}


void imprimirDatosGrafico(double** matriz, int dimension, FILE* archivo)
/* Imprime matriz de dimensiones dimension*2 */
{
  for (int i = 0; i < dimension; ++i) {
    fprintf(archivo, "%E %E\n", matriz[i][0], matriz[i][1]);
  }
  
}

void copiarVectores(double ** a, double** b, int size_x, int size_y)
// Funcion que copia la matriz "b" en la "a"
{
  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      a[i][j] = b[i][j];
    }
  }
}
