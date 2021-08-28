#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <math.h>
#include <iomanip>

void matrixAlloc(double **&a, const int p1, const int p2);

void matrixDealloc(double **a, const int p1);

void matrixAssign(double **a, double **b, const int p1, const int p2);

void matrixPrint(double **a, const int p1, const int p2);

double **matrixMultiply(double **a, double **b, const int p11, const int p21, const int p12, const int p22);

double norm(double *a, const int n);

void normalize(double *a, const int n);

void assign(double *s, double *t, const int n);

double dotProduct(double *a, double *b, const int n);
#endif