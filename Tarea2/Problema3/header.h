#ifndef BIB_H
#define BIB_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

void rellenarUinicial(double*U,double*x,int nx,double u0);
void evolucionFTBS(double*u_anterior,int nx,double gamma,double u0);
void evolucionFTCS(double*u_anterior,int nx,double gamma,double u0);
void evolucionLAX(double*u_anterior,int nx,double gamma,double u0);
void evolucionBTCS(double*u_anterior,int nx,double gamma,double u0);

void resuelveSistema(double ** A, double * b, double * x, int n);
void resuelveTridiagonal(double* sol_x, double* diagonal,double* superior, double* inferior, double* vector_b, double n);

#endif
