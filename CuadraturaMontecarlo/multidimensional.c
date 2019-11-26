/* Programa que realiza la cuadratura de montecarlo para resolver integrales
 * multidimensionales de mecanica estadística.
 *
 * Sean n átomos sujetos a un potencial U(r) dentro de una caja centrada en (0,0,0)
 *  de tamaño L.
 * Determinar promedios de observables (e.potencial promedio, etc).
 * 
 * O = int ( O(r) * e^(-beta U(r)) dr) / int ( e^(-beta U(r)) dr)
 *
 * O= I1 / I2
 *
 * Se toman m sistemas aleatorios, cada uno compuesto de n átomos que son representados
 *  por sus tres coordenadas, en total son 3n*m numeros aleatorios
 * 
 * 
 * Se aproxima I2 = V^n/m * sum g(Rj)
 * donde Rj es un sistema aleatorio 
 */


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

double funcionPotencial(double rx,double ry,double rz){
  double U = pow(rx,2)+pow(ry,2)+pow(rz,2);
  return U;
}

double funcionAIntegrar(double** R, int n, double parametro) //funcion g, R de tamaño [n,3]
{
  double y = 0;
  double sum = 0;
  for (int i = 0; i < n; ++i) {
    sum = funcionPotencial(R[i][0],R[i][1],R[i][2]) + sum;
  }
  /* TODO: ERROR, La suma es tan grande que exp(sum) es practicamente 0 */
  //printf("%f\n", sum);
  y = exp( -parametro*sum );
  return y;
}

int main(void )
{

  int num_moleculas = 20;
  int num_sistemas = 1e7;
  //volumen del espacio de configuracion total (incluyendo todas las moleculas)
  double volumen = pow(1, num_moleculas);
  double largo = 1; //caja cubica
  double beta = 2.0;
  
  gsl_rng * generador = gsl_rng_alloc (gsl_rng_taus);

  //las variables aleatoreas no pueden ser guardadas en un vector simple tan grande
  //ya que el comando linux ulimit -s dicta cuanta memoria stacked le da al programa
  
  double** sistemaR = malloc(sizeof(double*) * num_moleculas);
  for (int i = 0; i < num_moleculas; ++i) {
    sistemaR[i] = malloc(sizeof(double) * 3);
  }

  double fun_eval = 0;
  double sum = 0;

  for (int i = 0; i < num_sistemas; ++i) {
    // Obtenemos variable aleatoria
    for (int i = 0; i < num_moleculas; ++i) {
      for (int j = 0; j < 3; ++j) {
	sistemaR[i][j] = gsl_rng_uniform(generador)*largo-largo/2.0;
      }
      //printf("%f %f %f \n", sistemaR[i][0], sistemaR[i][1], sistemaR[i][2]);
    }
    //Realizamos evaluacion de funcion y suamos para obtener la integral
    fun_eval = funcionAIntegrar(sistemaR, num_moleculas, beta);
    //printf("%f\n", fun_eval);
    sum = sum + fun_eval;
  }

  double integral = volumen * sum / num_sistemas;

  printf("Integral = %E\n", integral);
  
  return 0;
}
