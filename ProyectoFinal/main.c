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

 
int main(/*int argc, char *argv[]*/)
{
  int num_moleculas = 24;
  int num_sistemas = 1e2;
  double aceptados = 1;
  double ratio_aceptacion;
  
  double beta = 2.0;
  double dmax = 0.5;
  double sigma = sqrt(14.0/256.0);

  /* Vector de parametros del sistema */
  double parametros[4] = {beta, dmax, pow(sigma,2),sigma/2};
  // 0 beta = 1/(kB T)
  // 1 dmax distantcia maxima de random
  // 2 pow(sigma,2) diametro al cuadrado
  // 3 sigma/2 radio


  /* Variables de histograma */
  // EL numero de particiones del histograma
  int num_hist = 1000;
  
  double ** hist = malloc(sizeof(double*) * num_hist);
  for (int i = 0; i < num_hist; ++i) {
    hist[i] = malloc(sizeof(double) * num_sistemas);
  }

  double ** prom_hist = malloc(sizeof(double*) * num_hist);
  for (int i = 0; i < num_hist; ++i) {
    prom_hist[i] = malloc(sizeof(double) * 2);
  }
  
  FILE *archivo;
  archivo = fopen("coordenadas.xyz", "w");

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

  int numero_estados = 0;
  inicializar(estado_actual, num_moleculas, parametros);
  
  for (int i = 0; i < num_sistemas; ++i) {
    // Realizamos algoritmo de Metropolis para cada sistema

    metropolis2(estado_actual, estado_siguiente, num_moleculas, parametros, &aceptados);

    /* Analizamos cuantos estados estamos aceptando */
    numero_estados = numero_estados + num_moleculas;
    
    ratio_aceptacion = aceptados / numero_estados;

    if(ratio_aceptacion > 0.55){
      parametros[1] = parametros[1]*1.05;
    }
    if(ratio_aceptacion < 0.45){
      parametros[1] = parametros[1]*0.95;
    }

    
    if(i % 10){
      printf("i %d ratio %E dmax %E\n", i,  ratio_aceptacion, parametros[1]);
      imprimirEstado(estado_siguiente, num_moleculas, archivo, parametros);
    }
    
    copiarVectores(estado_actual, estado_siguiente, num_moleculas, 3);

    /* TODO: Histograma */
    histograma(estado_actual, num_moleculas, hist, i, num_hist);
  }

  /* TODO: Analizar sistema */
  analizarHistograma(hist, num_sistemas, num_hist, prom_hist);

  /* TODO: Salida del sistema */
  imprimirEstado(estado_actual, num_moleculas, archivo, parametros);

  FILE *histograma;
  histograma = fopen("distribucion_radial.dat", "w");
  imprimirDatosGrafico(prom_hist, num_hist, histograma);
  fclose(histograma);
  
  fclose(archivo);
  free(estado_actual);
  free(estado_siguiente);
  free(hist);
  free(prom_hist);
  
  return 0;
}

