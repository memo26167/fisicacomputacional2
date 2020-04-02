/*
 * Calcula distribucion estacionaria y probabilidad de primer retorno (NO MATRICIAL) para el canal de comunicacion binaria,
 * en donde los estados posibles son 0 y 1
 * Autor: Guillermo Fonseca 
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

#define P 0.5 // probabilidad de transicion de estado 0 -> 1
#define Q 0.5 // probabilidad de transicion de estado 1 -> 0

void cadenaMarkov(void){
    // Usando Park-Miller
    gsl_rng *generador = gsl_rng_alloc(gsl_rng_minstd);

    // Parametros de la cadena
    int numero_semillas = 1000000;
    int largo_cadena = 1000;
    printf("Numero de semillas %d Largo de cadena %d \n", numero_semillas, largo_cadena);
    
    int estado_inicial = 0;
    int estado = 0;
    
    // Contadores para realizar estadística
    int cont_0 = 0;// contador general de todos los 0
    int cont_1 = 0;// contador general de todos los 1

    int cont_i0 = 0;// Numero de cadenas totales con estado inicial 0
    int cont_i1 = 0;// Numero de cadenas totales con estado inicial 1

    
    int bool = 1;// 1 true, 0 false, para ver si se realizó retorno o no
    
    
    int permanencia_0 = 0;// contador del número de pasos que permanece 0 (para calcular retornos de 1)
    int permanencia_1 = 0;// contador del número de pasos que permanece 1 (para calcular retornos de 0)

    double prob_retorno0 = 0; // probabilidad de que se retorne al estado 0
    double prob_retorno1 = 0; // probabilidad de que se retorne al estado 1
    double tiempo_retorno0 = 0; // numero de pasos que transucrren en promedio para pasar al estado 0
    double tiempo_retorno1 = 0; // numero de pasos que transucrren en promedio para pasar al estado 1

    /* El tiempo de retorno es muy sensible a cambios en el número de cadena 
       - Al aumentar el número de semillas, se promedia más el tiempo de retorno
       - Al aumentar el número de cadena, se acerca más al tiempo de retorno teórico
       - Al disminuir el número de cadena, se aleja del valor teorico, y no importa
         si aumentamos el número de semilla, el valor se erá afectado.
       Esto se debe a las probabilidades de transicion, si estas son muy bajas, necesitaremos
       cadenas mas largas para que se produzcan los retornos.
    */
    
    // un número aleatorio entre 0 y 1, que servira para evaluar la probabilidad
    double num_random = 0;
  
    /* Iteraremos, cambiando cada vez las semillas para realizar distintos generadores, de los cuales
     realizaremos cadenas de markov.
    */
    
    for (int s = 0; s < numero_semillas; ++s) {
	gsl_rng_set(generador, 3*s+2);
	estado = estado_inicial;

	
	if(estado_inicial == 0){++cont_i0;} // contar el numero de cadenas que empiezan en 0
	if(estado_inicial == 1){++cont_i1;} // contar el numero de cadenas que empiezan en 1
	
	bool = 1;
	permanencia_0 = 0;
	permanencia_1 = 0;

	 
	for (int n = 0; n < largo_cadena; ++n) {
	    /* GENERAR NUEVO ESTADO*/
	    
	    num_random = gsl_rng_uniform(generador);

	    // ***Se le aplica la matriz de transicion al estado 0***
	    if (estado == 0) {
		// se cuenta la cantidad de 0
		cont_0 += 1; 		
		// si num_random < P, se cambia de estado,
		if(num_random < P){estado = 1;}
		// caso contrario, se mantiene en el mismo estado.
	    }
	    
	    // ***Se le aplica la matriz de transicion al estado 1, similar al caso anterior***
	    else{
		cont_1 += 1; 			
		if(num_random < Q){estado = 0;}
	    }

	    /* Realizar contabilidad para probabilidad de retorno */
	    if (estado_inicial == 0){
		// Se suma la cantidad de veces que se permanece en 1 (fuera de 0)
		// Apenas se vuelva al estado 0, bool = 0 y nos quedaremos con ese valor de permanencia
		++permanencia_1;
		
		if (estado == 0 && bool == 1) {
		    bool = 0; // no volver a este if
		    // Sumamos las cadenas que se devuelven a 0
		    ++prob_retorno0;
		    // Sumamos el tiempo que se demora una cadena a volver a 0
		    tiempo_retorno0 += permanencia_1;
		}
	    }
	    // Similar al caso anterior, pero para 1
	    if (estado_inicial == 1){
		++permanencia_0;
		if (estado == 1 && bool == 1) {
		    bool = 0;
		    ++prob_retorno1;
		    tiempo_retorno1 += permanencia_0;
		}
	    }
	}
	/* Actualizar estado inicial */
	estado_inicial = estado;
    }
    
    
    prob_retorno0 = prob_retorno0/cont_i0;
    prob_retorno1 = prob_retorno1/cont_i1;
    tiempo_retorno0 = tiempo_retorno0/cont_i0;
    tiempo_retorno1 = tiempo_retorno1/cont_i1;

    printf("prob 0s: %f prob 1s: %f\n",
	   (double)cont_0/(numero_semillas*largo_cadena),
	   (double)cont_1/(numero_semillas*largo_cadena) );
    printf("prob_retorno0 %f tiempo_retorno0 %f\n", prob_retorno0, tiempo_retorno0);
    printf("prob_retorno1 %f tiempo_retorno1 %f\n", prob_retorno1, tiempo_retorno1);
    
    gsl_rng_free(generador);
}

int main(void)
{
    cadenaMarkov();
    return 0;
}


