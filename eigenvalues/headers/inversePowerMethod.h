#ifndef INVERSE
#define INVERSE
#include <iostream>
#include <math.h>
#include "matrix.h"
#include "LU.h"

void inversePowerMethod(double **A, double *vo,
                        double &lambda, const int n, const double eps);
#endif