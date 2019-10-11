/* Este programa resuelve Ecuacion Diferencial Parcial Hiperbolica de segundo orden, en una dimension espacial,
 * es decir:
 * 
 * ∂²U(x,t)/∂t² - k² ∂²U(x,t)/∂x² = 0
 * 
 * Sujeta a las condiciones de borde
 * 
 * Para el tiempo
 *    U(x,0)  = f(x)
 * ∂U(x,0)/∂t = g(x)
 * 
 * Para el espacio
 * U(0,t) = α(t)
 * U(L,t) = β(t)
 *
 * En donde se utiliza diferencias finitas para discretizar el problema
 *
 * ( U[x,t+1] -2U[x,t] + U[x,t-1] )/ Δt² - k²( U[x+1,t] -2U[x,t] + U[x-1,t] )/ Δx² = 0
 *
 * Que se puede resumir en la ecuacion vectorial:
 * 
 * U⃗[t+1]= A U⃗[t] - U⃗[t-1]
 * Sin embargo, si t=1, entonces se debe ocupar la condicion de borde.
 *
 * ∂U(x,0)/∂t = g(x)
 * (U[x,1]-U[x,0])/Δt = g(x)
 * Entonces
 * U[x,1] = g(x)Δt + U[x,0]
 * 
 * I
 */

#include <stdio.h>
#include <stdlib.h>




int main()
{

  
  return 0;
}
