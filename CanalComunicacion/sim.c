/**
   Calcula distribución estacionaria y prob. de primer retorno matricial (n=1,...,5) para el canal de comunicación binaria.
   Autor: Christopher Paredes 21/03/20
 */
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

#define P 0.1
#define Q 0.01
#define N 100000 //cantidad de pasos (n)
#define NS 10000 //cantidad de semillas
#define E_INIT 1
double p00_n(int n);
void distribucion_estacionaria(void);
void prob_primer_retorno(void);

double p00_n(int n){
	return (Q+pow((1-P-Q),n)*P)/(P+Q);	
}

//caso n=5 (por generalizar)
void prob_primer_retorno(void){
	int i=0,s=0;
	gsl_vector *b=gsl_vector_alloc(5);
	gsl_vector *f=gsl_vector_alloc(5);
	gsl_matrix *A=gsl_matrix_alloc(5,5);	
	gsl_permutation *p =gsl_permutation_alloc(5);
	double sum=0;
	gsl_matrix_set_zero(A);	

	for(i=0;i<5;++i){
		gsl_vector_set(b,i,p00_n(i+1));
	}

	for(i=0;i<5;++i){//diagonal
		gsl_matrix_set(A,i,i,p00_n(0));
	}

	gsl_matrix_set(A,1,0,p00_n(1));
	gsl_matrix_set(A,2,1,p00_n(1));
	gsl_matrix_set(A,3,2,p00_n(1));
	gsl_matrix_set(A,4,3,p00_n(1));

	gsl_matrix_set(A,2,0,p00_n(2));
	gsl_matrix_set(A,3,1,p00_n(2));
	gsl_matrix_set(A,4,2,p00_n(2));

	gsl_matrix_set(A,3,0,p00_n(3));
	gsl_matrix_set(A,4,1,p00_n(3));

	gsl_matrix_set(A,4,0,p00_n(4));

	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,f);
	double mu=0;
	double f0_i;
	for(i=0;i<5;++i){
		f0_i=gsl_vector_get(f,i);
		printf("f00^(%i):%f\n",i+1,f0_i);
		sum=sum+f0_i;
		mu=mu+(i+1)*f0_i;		
	}
	printf("f_0:%f\n",sum);
	printf("mu_0:%f\n",mu);	
}

// Simulación de NS cadenas de markov para obtener la distribución estacionaria
void distribucion_estacionaria(void){
	int i=0,s=0;
	gsl_rng *r=gsl_rng_alloc(gsl_rng_minstd);
	int E=0;
	double x=0.0;
	int cont_1s=0;
	int cont_0s=0;   
	
	for(s=0;s<NS;++s){
		E=E_INIT;
		gsl_rng_set(r,3*s+2);
		
		for(i=0;i<N;++i){
			x=gsl_rng_uniform(r);
			if(E==0){
				if(x<P){
					E=1;
				}
			}else{
				if(x<Q){
					E=0;
				}
			}

			if(i==(N-1)){
				if(E==0){
					cont_0s++;
				}else{
					cont_1s++;
				}
			}
		}
	}
	
	printf("conteo 0s:%i conteo 1s:%i\n",cont_0s,cont_1s);
	printf("q1:%.10f q2:%.10f\n",(double)cont_0s/NS,(double)cont_1s/NS);

	gsl_rng_free(r);	
}

int main(int argc, char *argv[]){
	distribucion_estacionaria();
	prob_primer_retorno();
    return 0;
}
