#ifndef LU
#define LU
#include <iostream>

void LU_Decomposition(double **a, double **l, double **u, const int n);

double *LUsolver(double **l, double **u, double *b, const int n);
#endif