/* Programa que realiza la cuadratura de montecarlo para resolver integrales
 * multidimensionales de mecanica estadística mediante algoritmo de metropolis.
 *
 * Sean n átomos con hamiltoniano H
 * Objetivo: Determinar promedios de observables (presion, etc).
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
 * Al cabo de cierta cantidad de sistemas, se espera que el sistema
   siga la distribucion q(x), esto es termalizacion    
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"

 
int main(/*int argc, char *argv[]*/)
{
  int num_moleculas = 24;  
  int num_sistemas = 1e6;//debe ser 1e6 
  int hist_aceptados = 0;
  double aceptados = 1;
  double ratio_aceptacion;
  
  double dmax = 0.5;
  double sigma = sqrt(14.0/256.0);

  
  /* Vector de parametros del sistema */
  double parametros[3] = {sigma, dmax , hist_aceptados};
  // 0 sigma diametro
  // 1 dmax desplazamiento maximo de movimientos aleatorios
  // 2 numero de histogramas aceptados (evitar tomarlos hasta alcanzar el equilibrio)

  
  /* Variables de histograma */

  // EL numero de particiones del histograma
  int num_hist = 1000;

  // Los histogramas
  double * hist = malloc(sizeof(double) * num_hist);

  // El histograma promedio con su coordenada de distancia y de distribucion g(r)
  double ** prom_hist = malloc(sizeof(double*) * num_hist);
  for (int i = 0; i < num_hist; ++i) {
    prom_hist[i] = malloc(sizeof(double) * 2);
  }

  // Archivo de salida del sistema en formato XYZ
  FILE *archivo;
  archivo = fopen("coordenadas.xyz", "w");

  // El generador de números aleatorios
  generador_uniforme = gsl_rng_alloc(gsl_rng_taus);

  /* El vector de estado */
  
  // Las variables aleatoreas no pueden ser guardadas en un vector simple tan grande
  // ya que el comando linux ulimit -s dicta cuanta memoria stacked le da al programa.
  // Por lo tanto, utilizamos malloc.
  double** estado= malloc(sizeof(double*) * num_moleculas);
  for (int i = 0; i < num_moleculas; ++i) {
    estado[i] = malloc(sizeof(double) * 3);
  }

  /******** PROCEDIMIENTOS ********/

  /**** Inicializacion ****/
  inicializar(estado, num_moleculas, parametros);

  
  // Realizamos algoritmo de Metropolis para cada sistema
  for (int i = 0; i < num_sistemas; ++i) {
    printf("%f \n", (double)i/(double)num_sistemas*100.0);
    aceptados = 0;
    
    /**** Metropolis ****/
    metropolis(estado, num_moleculas, parametros, &aceptados);
 
    ratio_aceptacion = aceptados / num_moleculas;    

    // Se modifica el desplazamiento maximo en funcion del ratio de aceptacion
    if(ratio_aceptacion > 0.55){
      parametros[1] = parametros[1]*1.05;
    }
    if(ratio_aceptacion < 0.45){
      parametros[1] = parametros[1]*0.95;
    }

    /**** Histograma ****/
    if(i>num_sistemas*0.25){
      parametros[2] = parametros[2] + 1;
      histograma(estado, num_moleculas, hist, num_hist, parametros);
    }

    /* Posible Salida XYZ para analizar en OVITO */
    /* if(1%10){ */
    /*   imprimirEstado(estado, num_moleculas, archivo, parametros); */
    /* } */
  }

  /**** Analisis de datos ****/
  analizarHistograma(hist, num_hist, num_moleculas, prom_hist, parametros);

  char titulo_archivo[100];
  sprintf(titulo_archivo, "distribucion_radial_%d.dat", num_moleculas);
  printf("%s\n", titulo_archivo);
  
  FILE *histograma;
  histograma = fopen(titulo_archivo, "w");
  imprimirDatosGrafico(prom_hist, num_hist, histograma);

  /**** Liberar y terminar programa ****/
  fclose(histograma);
  fclose(archivo);
  free(estado);
  free(hist);
  free(prom_hist);
  
  return 0;
}

