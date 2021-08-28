#include "../headers/qr.h"

using namespace std;

double squareSum(double **a, const int n){
    double sum = 0;
    for (int j = 0; j < n-1; ++j)
        for (int i = j+1; i < n; ++i)
            sum += a[i][j]*a[i][j];
        
    return sum;
}

double ** calculateHj(double **a, const int i, const int n){
    double *w = new double[n];
    double *w_ = new double[n];
    double *N = new double[n];

    for (int k = 0; k < n; k++)
    {
        w[k] = 0;
        w_[k] = 0;
        N[k] = 0;
    }

    for (int k = i; k < n; ++k)
        w[k] = a[k][i];

    w_[i] = norm(w, n);

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

    matrixDealloc(n1, n);
    matrixDealloc(n1T, 1);
    matrixDealloc(aux, n);
    delete[] w;
    delete[] w_;
    delete[] N;
    return h;
}

void qrDecomp(double **a, double **&q, double **&r, const int n){
    double **a_new; matrixAlloc(a_new, n,n);
    double **a_old; matrixAlloc(a_old, n,n);
    
    matrixAlloc(q, n, n);
    matrixAlloc(r, n, n);

    for (int i = 0; i < n; ++i)
        q[i][i] = 1;

    double **hj;

    matrixAssign(a_new, a, n, n);
    for (int j = 0; j < n-1; ++j)
    {
        matrixAssign(a_old, a_new, n,n);
        hj = calculateHj(a_old, j, n);
        
        double **a_new_ = matrixMultiply(hj,a_old, n,n,n,n);
        matrixAssign(a_new, a_new_, n, n);

        double **q_ = matrixMultiply(q, hj, n, n, n, n);
        matrixAssign(q, q_, n, n);

        matrixDealloc(q_, n);
        matrixDealloc(a_new_, n);
        matrixDealloc(hj, n);
    }

    matrixAssign(r, a_new, n, n);
    matrixDealloc(a_new, n);
    matrixDealloc(a_old, n);
}

void qrMethod(double **a, double **&a_new, double **&p, const int n, const double epls){
    
    matrixAlloc(a_new, n , n);
    double **a_old; matrixAlloc(a_old, n , n);

    double **q, **r;

    matrixAlloc(p, n, n);
    for (int i = 0; i < n; ++i)
        p[i][i] = 1;    
    
    double *vector = new double[n];
    double val = 1;

    matrixAssign(a_new, a, n, n);
    int iter = 0;
    while(val > epls){
        //cout << "Iteration: " << ++iter << endl;
        matrixAssign(a_old, a_new, n, n);

        qrDecomp(a_old, q, r, n);
        /*
        
        cout << "\nQ: \n";
        matrixPrint(q, n, n);
        cout << "\nR: \n";
        matrixPrint(r, n, n);
        
        */
        double **a_new_ = matrixMultiply(r,q,n,n,n,n);
        matrixAssign(a_new, a_new_, n, n);
        //matrixPrint(a_new, n, n);

        double **p_ = matrixMultiply(p, q, n, n, n, n);
        matrixAssign(p, p_, n, n);
        val = squareSum(a_new, n);


        matrixDealloc(a_new_, n);
        matrixDealloc(p_, n);
    }

    matrixDealloc(a_old, n);
} 
