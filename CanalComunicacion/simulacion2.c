/*
 * Calcula distribucion estacionaria y probabilidad de primer retorno (NO MATRICIAL) para el canal de comunicacion binaria,
 * en donde los estados posibles son 0 y 1
 * Autor: Guillermo Fonseca 
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

#define P 0.1 // probabilidad de transicion de estado 0 -> 1
#define Q 0.01 // probabilidad de transicion de estado 1 -> 0


void multiplicacion_matricial(double** a, double** b, double** c, int n){
    // A of size nxn
    double suma = 0;
    
    double A[n][n];

    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    suma = 0;
	    for (int k = 0; k < n; ++k) {
		suma += b[i][k]*c[k][j];
	    }
	    A[i][j] = suma;
	}
    }
    
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    a[i][j] = A[i][j];
	}
    }
    //        printf("a = %f %f \n    %f %f\n", a[0][0], a[0][1], a[1][0], a[1][1]);
    
}

void tiempo_recurrencia(void){
    // cuantas recurrencia se realizarán
    int num_sist = 2000;
    
    //las probabilidades de transicion de estado x hasta el paso n
    double prob_transicion_0[num_sist];
    double prob_transicion_1[num_sist];
    
    // matriz de transicion
    double ** m_t = malloc(sizeof(double*) * 2);

    // matriz de transicion eleavado a n
    double ** m_n = malloc(sizeof(double*) * 2);
    
    for (int i = 0; i < 2; ++i) {
	m_t[i] = malloc(sizeof(double) * 2);
	m_n[i] = malloc(sizeof(double) * 2);
    }
    m_t[0][0] = 1 - P;
    m_t[0][1] = P;
    m_t[1][0] = Q;
    m_t[1][1] = 1 - Q;

    m_n[0][0] = 1 - P;
    m_n[0][1] = P;
    m_n[1][0] = Q;
    m_n[1][1] = 1 - Q;


    // sacar probabilidad de transición de un estado al mismo estado, a n pasos
    prob_transicion_0[0] = m_t[0][0];
    prob_transicion_1[0] = m_t[1][1];
    for (int i = 1; i < num_sist; ++i) {
	multiplicacion_matricial(m_n, m_n, m_t, 2);
	prob_transicion_0[i] = m_n[0][0];
	prob_transicion_1[i] = m_n[1][1];
    }

 
    /* Calculo de probabilidad de retorno 
      Mediante formula de recurrencia 
      Que sencillamente es resolver problema diagonal inferior
      O sustitucion hacia atrás
    */ 
    
    // probabilidades de primer retorno a paso n, del estado k
    double ret_0[num_sist-1];
    double ret_1[num_sist-1];
    
    double sum = 0;//para acumular la formula de recurrencia
       

    ret_0[0] = prob_transicion_0[0];
    ret_1[0] = prob_transicion_1[0];

    for (int i = 1; i < num_sist - 1; ++i) {
	sum = 0;
        for (int k = 0; k < i; ++k) {
	    sum += ret_0[k]*prob_transicion_0[i-k-1];
	}
	ret_0[i] = (prob_transicion_0[i] - sum);

	sum = 0;
        for (int k = 0; k < i; ++k) {
	    sum += ret_1[k]*prob_transicion_1[i-k-1];
	}
	ret_1[i] = (prob_transicion_1[i] - sum);
    }

    // Calculo probabilidad de retorno y tiempo promedio de retorno

    double retorno_0, retorno_1;
    double tiempo_0, tiempo_1;
    for (int i = 0; i < num_sist - 1; ++i) {
	retorno_0 += ret_0[i];
	retorno_1 += ret_1[i];
	tiempo_0 += (i+1)*ret_0[i];
	tiempo_1 += (i+1)*ret_1[i];
    }
 
    printf("estado 0 prob retorno: %f tiempo retorno %f\n", retorno_0, tiempo_0);
    printf("estado 1 prob retorno: %f tiempo retorno %f\n", retorno_1, tiempo_1);

    free(m_t);
    free(m_n);
}


void cadenaMarkov(void){
    // Usando Park-Miller
    gsl_rng *generador = gsl_rng_alloc(gsl_rng_minstd);

    // Parametros de la cadena
    int numero_semillas = 1000000;
    int largo_cadena = 100;
    int estado_inicial = 0;
    int estado = 0;
  
    // Contadores para realizar estadística
    int cont_0 = 0;// contador general de todos los 0
    int cont_1 = 0;// contador general de todos los 1
  
    // un número aleatorio entre 0 y 1, que servira para evaluar la probabilidad
    double num_random = 0;
  
    // Iteraremos, cambiando cada vez las semillas para realizar distintos generadores, de los cuales
    // realizaremos cadenas de markov.
    estado = estado_inicial;
    for (int s = 0; s < numero_semillas; ++s) {
	gsl_rng_set(generador, 3*s+2);

	for (int n = 0; n < largo_cadena; ++n) {
	    /* GENERAR NUEVO ESTADO*/
      
	    num_random = gsl_rng_uniform(generador);
	    // ***Se le aplica la matriz de transicion al estado 0***
	    if (estado == 0) {
		cont_0 += 1; // se cuenta la cantidad de 0
	
		// si num_random < P, se cambia de estado,
		if(num_random < P){
		    estado = 1;
		}
		// caso contrario, se mantiene en el mismo estado.
	    }   
	    // ***Se le aplica la matriz de transicion al estado 1***
	    else{
		cont_1 += 1; // se cuenta la cantidad de 0	

		// si num_random < Q, se cambia de estado,
		if(num_random < Q){
		    estado = 0;
		}
		// caso contrario, se mantiene en el mismo estado.	
	    }

      
	}
    }

    printf("conteo 0s: %d conteo 1s: %d\n", cont_0, cont_1);
    printf("prob 0s: %f prob 1s: %f\n",
	   (double)cont_0/(numero_semillas*largo_cadena),
	   (double)cont_1/(numero_semillas*largo_cadena) );

    printf("teorico prob 0s: %f prob 1s: %f\n",
	   Q/(P+Q),
	   P/(P+Q));
    
    gsl_rng_free(generador);
}

int main(void)
{
    //cadenaMarkov();
    tiempo_recurrencia();
    return 0;
}

