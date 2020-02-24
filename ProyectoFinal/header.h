#ifndef BIB_H
#define BIB_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <float.h>
#include <string.h>

gsl_rng * generador_uniforme;

void inicializar(double ** estado, int num_moleculas, double* parametros);

void metropolis
(double **estado, int num_moleculas, double* parametros, double* aceptados);

void imprimirEstado(double ** estado,int num_moleculas, FILE* archivo, double* parametros);

void copiarVectores(double ** actual, double** siguiente, int size_x, int size_y);

void histograma(double ** estado, int num_moleculas, double* hist, int num_hist, double* parametros);

void analizarHistograma
(double * hist, int num_hist, int num_moleculas, double** prom_hist, double* parametros);

void imprimirDatosGrafico(double** matriz, int dimension, FILE* archivo);
#endif
