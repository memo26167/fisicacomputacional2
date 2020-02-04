#include <stdlib.h>
#include <stdio.h>
#include "header.h"

void copiarVectores(double ** actual, double** siguiente, int size_x, int size_y)
{
  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      actual[i][j] = siguiente[i][j];
    }
  }
}
