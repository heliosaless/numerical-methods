#ifndef HOUSEHOLDER
#define HOUSEHOLDER
#include "matrix.h"

double **houseHolder_(double **a, const int i, const int n);
double **houseHolder(double **a, double **&h, const int n);

#endif