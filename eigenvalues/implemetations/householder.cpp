#include "../headers/householder.h"
using namespace std;

double **houseHolder_(double **a, const int i, const int n)
{
    double *w = new double[n];
    double *w_ = new double[n];
    double *N = new double[n];

    for (int k = 0; k < n; k++)
    {
        w[k] = 0;
        w_[k] = 0;
        N[k] = 0;
    }

    for (int k = i + 1; k < n; ++k)
        w[k] = a[k][i];

    w_[i + 1] = norm(w, n);
     
    for (int k = 0; k < n; ++k)
        N[k] = w[k] - w_[k];

    normalize(N, n);

    double **h;
    matrixAlloc(h, n, n);
    double **n1;
    matrixAlloc(n1, n, 1);
    double **n1T;
    matrixAlloc(n1T, 1, n);

    for (int k = 0; k < n; ++k)
    {
        n1[k][0] = N[k];
        n1T[0][k] = N[k];
    }

    double **aux = matrixMultiply(n1, n1T, n, 1, 1, n);

    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
        {
            if (k == j)
                h[k][j] = 1 - 2 * aux[k][j];
            else
                h[k][j] = -2 * aux[k][j];
        }

    //matrixPrint(h, n, n);
    matrixDealloc(n1, n);
    matrixDealloc(n1T, 1);
    matrixDealloc(aux, n);
    delete[] w;
    delete[] w_;
    delete[] N;
    return h;
}

double **houseHolder(double **a, double **&h, const int n)
{
    double **h_;
    double **hi;
    double **a_New_;
    double **a_New;
    double **a_Old;

    matrixAlloc(a_Old, n, n);

    matrixAlloc(h, n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i == j)
                h[i][j] = 1;

    matrixAssign(a_Old, a, n, n);

    for (int i = 0; i < n - 2; ++i)
    {
        hi = houseHolder_(a_Old, i, n);
        //matrixPrint(hi, n, n);

        a_New_ = matrixMultiply(a_Old, hi, n, n, n, n);
        a_New = matrixMultiply(hi, a_New_, n, n, n, n);
        matrixDealloc(a_New_, n);

        matrixAssign(a_Old, a_New, n, n);
        h_ = matrixMultiply(h, hi, n, n, n, n);
        matrixAssign(h, h_, n, n);

        matrixDealloc(h_, n);
        matrixDealloc(hi, n);
        if (i < n - 3)
            matrixDealloc(a_New, n);
    }

    matrixDealloc(a_Old, n);

    return a_New;
}