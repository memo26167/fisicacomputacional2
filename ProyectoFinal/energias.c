#include <stdio.h>
#include <stdlib.h>
#include "header.h"


double energiaCinetica(double **estado,int num_moleculas, double* parametros);

double energiaPotencial(double **estado,int num_moleculas, double* parametros)
// Potencial de Lennard Jones con sigma y epsilon = 1
{
  double potencial = 0;

  //diferencia de las posiciones en x ,y ,z
  double dif_posiciones[3];
  double dif_cuadrado;

  double rcorte_cuadrado = 12.25;//3.5 al cuadrado
  double aux = 0;//valor auxiliar
  
  for (int i = 0; i < num_moleculas-1; ++i) {
    for (int j = i+1; j < num_moleculas; ++j) {
      dif_posiciones[0] = estado[i][0] - estado[j][0];// - largo_cajaX*round( (estado[i][0] - estado[j][0])/largo_cajaX )
      dif_posiciones[1] = estado[i][1] - estado[j][1];
      dif_posiciones[2] = estado[i][2] - estado[j][2];
      dif_cuadrado = pow(dif_posiciones[0],2) + pow(dif_posiciones[1],2) + pow(dif_posiciones[2],2);

      if(dif_cuadrado < rcorte_cuadrado){
	aux = 4*( 1.0/pow(dif_cuadrado,6) - 1.0/pow(dif_cuadrado,3) );
      }
      else{
	aux = 0;
      }
      potencial = potencial + aux;
    }
  }
  return potencial;
}


double funcionEnergia(double **estado,int num_moleculas, double* parametros)
{
  return energiaCinetica(estado, num_moleculas, parametros) + energiaPotencial(estado, num_moleculas, parametros);
}
