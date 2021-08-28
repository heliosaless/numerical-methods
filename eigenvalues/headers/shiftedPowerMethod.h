#ifndef SHIFTED
#define SHIFTED

#include <iostream>
#include <math.h>
#include "matrix.h"
#include "inversePowerMethod.h"

void shiftedPowerMethod(double **A, double *vo,
                        double &lambda, const int n,
                        const double eps, const int U);

#endif