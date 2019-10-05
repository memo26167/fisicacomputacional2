#include <stdio.h>
#include <stdlib.h>
#include "header.h"

void resolver_metodo(ecp problema){
  //Se definen las variables del problema
  int n=problema.nx;
  int m=problema.ny;
  
  double delta_x=problema.tol;
  double delta_y=problema.tol;
  double delta_t=0.3*pow(delta_x,2)/problema.D;
  double gammax=problema.D*delta_t/(pow(delta_x,2));
  double gammay=problema.D*delta_t/(pow(delta_y,2));  



  

  double mA[n][m];//Matriz Actual
  double mS[n][m];//Matriz Siguiente
  
  //condicion inicial
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      mA[i][j]=problema.m_inicial[i][j];
    }    
  }

  //puntos siguientes fuera de los bordes
  for (int i = 1; i < n-1; ++i) {
    for (int j = 1; j < m-1; ++j) {
      mS[i][j]=(1 - 2*gammax - 2*gammay)*mA[i][j]
	+ gammax*(mA[i+1][j] + mA[i-1][j])
	+ gammay*(mA[i][j+1] + mA[i][j-1]);
    }
  }
  
  /* Bordes sin esquinas */
  //Borde inferior
  for (int i = 1; i < n-1; ++i) {
  }
  
}
