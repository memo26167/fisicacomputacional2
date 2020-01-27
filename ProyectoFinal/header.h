#ifndef BIB_H
#define BIB_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

gsl_rng * generador_uniforme;

double energiaCinetica(double **estado,int num_moleculas, double* parametros);
double energiaPotencial(double **estado,int num_moleculas, double* parametros);
double funcionEnergia(double **estado,int num_moleculas, double* parametros);



#endif
