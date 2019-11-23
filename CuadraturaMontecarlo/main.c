/* Programa que realiza la cuadratura de montecarlo para resolver integrales
 * unidimensionales.
 *
 * Sea I = int_a^b ( f(x)dx )
 * Se toman N numeros aleatorios entre a y b
 * Se aproxima I = (b-a) / N * sum( f(xi) ) 
 */


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>


double funcionAIntegrar(double x)
{
  double y = x;
  return y;
}

int main(void )
{

  double a = 0.0;
  double b = 1.0;
  
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

  int N = 1e7;
  //las variables aleatoreas no pueden ser guardadas en un vector simple tan grande
  //ya que el comando linux ulimit -s dicta cuanta memoria stacked le da al programa
  
  double var_rand = 0.0;
  double fun_eval = 0;
  double sum = 0;

  for (int i = 0; i < N; ++i) {
    // Obtenemos variable aleatoria
    var_rand = gsl_rng_uniform(r)/(b-a); 
    fun_eval = funcionAIntegrar(var_rand);
    sum = sum + fun_eval;
  }

  double integral = (b-a)*sum/N;

  printf("Integral = %f\n", integral);
  
  return 0;
}
