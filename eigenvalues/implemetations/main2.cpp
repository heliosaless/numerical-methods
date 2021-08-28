#include "../headers/qr.h"
#include "../headers/householder.h"
#include <iostream>

using namespace std;


int main(int argc, char const *argv[])
{
    const int n = 5;

    double **A = new double *[n];
    for (int i = 0; i < n; ++i)
        A[i] = new double[n];

    double A_[n][n] = {
        {40, 8, 4, 2, 1},
        {8, 30, 12, 6, 2},
        {4, 12, 20, 1, 2},
        {2, 6, 1, 25, 4},
        {1, 2, 2, 4, 5}};

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = A_[i][j];

    cout << "A:\n";
    matrixPrint(A, n, n);

    double **a_new, **p;
    qrMethod(A, a_new, p, n, 10e-7);

    cout << "Eigenvalues:\n";
    matrixPrint(a_new, n, n);
    cout << "Eigenvectors:\n";
    matrixPrint(p, n, n);


    double **h;
    double **tridiag = houseHolder(A, h, n);
    
    cout << "Tridiag: \n";
    matrixPrint(tridiag, n, n);

    double **a_new2, **p2;
    qrMethod(tridiag, a_new2, p2, n, 10e-7);

    cout << "Eigenvalues:\n";
    matrixPrint(a_new2, n, n);
    cout << "Tridiagonal Eigenvectors:\n";
    matrixPrint(p2, n, n);

    double **e = matrixMultiply(h, p2, n, n, n, n);
    cout << "A Eigenvectors: \n";
    matrixPrint(e, n, n);


    return 0;
}
