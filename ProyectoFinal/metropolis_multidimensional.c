/* Programa que realiza la cuadratura de montecarlo para resolver integrales
 * multidimensionales de mecanica estadística mediante algoritmo de metropolis.
 *
 * Sean n átomos con hamiltoniano H
 * Objetivo: Determinar promedios de observables (e.potencial promedio, etc).
 * 
 * O = int ( O(x) * q(x)dx) 
 *
 * Donde q(x) es la distribucion de probabilidad del sistema en el estado x
 * Sea p(x) = Z q(x), siendo p(x) la distribucion no normalizada
 * Problema, no conocemos q(x), sólo p(x)
 * 
 * Algoritmo de metropolis.
 * Objetivo: Crear serie de estados {x1,x2..xn},
    proponiendo un xt=xn + perturbacion que sigan la distribucion q(x)
 * Se define la probabilidad de aceptacion de xt como:
 *
 * Pr(A|xt,xn) = min( p(xt)/p(xn), 1) 
 *
 * Se elije x0 aleatoriamente y se generan M sistemas aleatorios,
    cada uno compuesto de N átomos que son representados
    por sus tres coordenadas espaciales
    , en total son 3N*M numeros aleatorios
    
    
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"

 
int main(int argc, char *argv[])
{
  int num_moleculas = 20;
  int num_sistemas = 1e7;
  double beta = 2.0;
  double dmax = 0.5;
  double parametros[2] = {beta, dmax};
  
  generador_uniforme = gsl_rng_alloc(gsl_rng_taus);

  //las variables aleatoreas no pueden ser guardadas en un vector simple tan grande
  //ya que el comando linux ulimit -s dicta cuanta memoria stacked le da al programa
  //por lo tanto, utilizamos malloc
  
  double** estado_actual = malloc(sizeof(double*) * num_moleculas);
  double** estado_siguiente = malloc(sizeof(double*) * num_moleculas);
  for (int i = 0; i < num_moleculas; ++i) {
    estado_actual[i] = malloc(sizeof(double) * 3);
    estado_siguiente[i] = malloc(sizeof(double) * 3);
  }

  /* TODO: Inicializar estados */
  
  double fun_eval = 0;
  double sum = 0;

  for (int i = 0; i < num_sistemas; ++i) {
    // Obtenemos variable aleatoria

    siguienteEstado(estado_actual, estado_siguiente, num_moleculas, parametros);
    
    //Realizamos evaluacion de funcion y suamos para obtener la integral
    fun_eval = funcionEnergia(estado_siguiente, num_moleculas, parametros);
    
    sum = sum + fun_eval;
  }

  double integral = sum / num_sistemas;
  
  printf("Integral = %E\n", integral);
  
  return 0;
}

