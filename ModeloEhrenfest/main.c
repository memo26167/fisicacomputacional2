/* Programa que simula numericamente el modelo de Ehrenfest
   En donde se tienen N moleculas, y dos cajas A y B.
   Inicialmente la caja B tiene N moleculas y la caja A 0 molecula
   Se toma una molecula de forma aleatoria, se reconoce la caja y se cambia de caja.

   Autor: Guillermo Fonseca*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#define N 100
#define Pasos 200
#define Print_Usuario 0

void modeloEhrenfest(void){
    /* Variable caja contendra cada molecula
       Se le asignara un 1 si la molecula respectiva existe*/
    double caja_A[N];

    /* Variable para exportar un gráfico de distribucion de probabilidad estacionaria */
    double histograma[N];
    
    for (int i = 0; i < N; ++i) {
	caja_A[i] = 0;
	histograma[i] = 0;
    }

    gsl_rng * generador = gsl_rng_alloc (gsl_rng_taus);

    /* TODO: 
       - Inicializar contadores
       - realizar estadisticad4
    */
    int num_pasos = Pasos;
    int num_moleculas = 0;
    int x = 0;

    

    for (int i = 0; i < num_pasos; ++i) {
	gsl_rng_set(generador, 2*x+1+i);
	
	/* Se toma un numero random entero x entre 0 y N-1
	   y se cambia el valor de caja_A[x] a 1 si inicialmente era 0 o a 0 si inicialmente era 1 */
	x = gsl_rng_uniform(generador)*N;
	
	if (caja_A[x] == 1) {caja_A[x] = 0;}
	else {caja_A[x] = 1;}
	
	/* Se refresca el valor de moleculas */
	num_moleculas = 0;
	/* Se cuenta el número de moleculas en la caja */
	for (int j = 0; j < N; ++j) { num_moleculas += caja_A[j];}
	
	/* Se le suma 1 al histograma para la cantidad de moleculas num_moleculas */
	++histograma[num_moleculas-1];

	if (!Print_Usuario){
	    printf("%d %d\n", i, num_moleculas);
	}
	
    }
    
    double suma = 0;

    for (int i = 0; i < N; ++i) {
	histograma[i] = histograma[i]/num_pasos;
    }

    if (!Print_Usuario){
	printf("\n");
	for (int i = 0; i < N; ++i) {
	    printf("%d %f\n", i+1, histograma[i]);
	}
    }


    if (Print_Usuario){
	printf("La probabilidad de que a %d paso\n", num_pasos);
	for (int i = 0; i < N; ++i) {
	    suma += histograma[i];
	    printf("existan %d moleculas en la caja es de : %f\n", i+1, histograma[i]);	
	}
	printf("suma de probabilidades suma = %f\n", suma);
    }
}

int main(void)
{
    modeloEhrenfest();
    
    return 0;
}

