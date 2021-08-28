#ifndef QR
#define QR
#include "matrix.h"
double squareSum(double **a, const int n);

void qrDecomp(double **a, double **&q, double **&r, const int n);

void qrMethod(double **a, double **&a_new, double **&p, const int n, const double epls);

#endif