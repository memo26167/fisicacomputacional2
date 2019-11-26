/* Programa que realiza la cuadratura de montecarlo para resolver integrales
 * multidimensionales de mecanica estadística mediante algoritmo de metropolis.
 *
 * Sean n átomos con hamiltoniano H
 * Determinar promedios de observables (e.potencial promedio, etc).
 * 
 * O = int ( O(x) * q(x)dx) 
 *
 * Donde q(x) es la distribucion de probabilidad del sistema en el estado x
 * Sea p(x) = Z q(x), siendo p(x) la distribucion no normalizada
 * 
 * Algoritmo de metropolis.
 * Objetivo: Crear serie de estados {x1,x2..xn},
    proponiendo un xt=xn + perturbacion 
 * Se define la probabilidad de aceptacion de xt como
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
#include <gsl/gsl_rng.h>
#include <math.h>

gsl_rng * generador_uniforme;

double funcionEnergia(double **estado,int num_moleculas, double* parametros);


double funcionDistribucion(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
// La distribucion de probabilidad no normalizada, de que ocurra un estado siguiente dado el estado actual
{
  double p = 0;
  double diferencia_energia = funcionEnergia(estado_siguiente, num_moleculas, parametros)
                              - funcionEnergia(estado_actual, num_moleculas, parametros);
  // parametros[0] es beta = 1 / (kB * T)
  p = exp(-parametros[0]*diferencia_energia);
  return p;
}


double probabilidadAceptacion(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
{
  double prob = funcionDistribucion(estado_actual,estado_siguiente,num_moleculas,parametros);
  if (prob>= 1.0){
    prob = 1.0;
  }
  return prob;
}


// Esta funcion contiene el algoritmo de Metropolis
void siguienteEstado(double **estado_actual, double **estado_siguiente, int num_moleculas, double* parametros)
{ 
  // parametros[1] es dmax
  // Volumen=1
  
  double dmax = parametros[1];
  double probabilidad;
  //Posible estado siguiente
  while(1){
    for (int i = 0; i < num_moleculas; ++i) {
      for (int j = 0; j < 3; ++j) {
	estado_siguiente[i][j] = estado_actual[i][j] + dmax*gsl_rng_uniform(generador_uniforme);
      }
    }
    probabilidad = probabilidadAceptacion(estado_actual,estado_siguiente,num_moleculas,parametros); 

    if( probabilidad == 1.0 ){
      break; //aceptamos el nuevo estado
    }
    else{
      if( probabilidad >= gsl_rng_uniform(generador_uniforme)){
	break; //aceptamos el nuevo estado
      }
    }
    
  }  
}

int main(void )
{

  int num_moleculas = 20;
  int num_sistemas = 1e7;
  double beta = 2.0;
  double dmax = 0.5;
  double parametros[2] = {beta, dmax};
  
  generador_uniforme = gsl_rng_alloc(gsl_rng_taus);

  //las variables aleatoreas no pueden ser guardadas en un vector simple tan grande
  //ya que el comando linux ulimit -s dicta cuanta memoria stacked le da al programa
  
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
