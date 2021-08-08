#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

void matrixAlloc(double **&a, const int p1, const int p2){
    a = new double*[p1];
    for (int i = 0; i < p1; ++i){
        a[i] = new double[p2];   
        for (int j = 0; j < p2; ++j)
            a[i][j] = 0;        
    }    
}

void matrixDealloc(double **a, const int p1){
    for (int i = 0; i < p1; ++i)
        delete[] a[i];
    delete[] a;       
}

void matrixAssign(double **a, double **b, const int p1, const int p2){
    for (int i = 0; i < p1; ++i)
        for (int j = 0; j < p2; ++j)
            a[i][j] = b[i][j];
}

void matrixPrint(double **a, const int p1, const int p2){
    for (int i = 0; i < p1; ++i){
        for (int j = 0; j < p2; ++j)
            cout << fixed << setprecision(3) <<  a[i][j] << "\t ";
        cout << endl;
    }
    cout << endl;
}

double** matrixMultiply(double **a, double **b, const int p11, const int p21, const int p12, const int p22){
    if(p12 != p21) {cout << "multiplication error" << endl; return nullptr;}
    double **c;
    matrixAlloc(c, p11, p22);

    double sum;
    for (int i = 0; i < p11; ++i)
        for (int j = 0; j < p22; ++j){
            sum = 0.;
            for (int k = 0; k < p12; ++k) sum += a[i][k]*b[k][j];

            c[i][j] = sum;
        }

    return c;
}

double norm(double *a, const int n){
    double sum = 0.;
    for (int i = 0; i < n; ++i){
        sum += a[i]*a[i];
    }
    return sqrt(sum);
}

void normalize(double *a, const int n){
    double s = norm(a, n);
    for (int i = 0; i < n; ++i)
        a[i] /= s;     
}

double** houseHolder_(double **a, const int i, const int n){
    double *w = new double[n];
    double *w_ = new double[n];
    double *N = new double[n];

    for (int k = 0; k < n; k++)
    {
        w[k]  = 0;
        w_[k] = 0;
        N[k]  = 0;
    }

    for (int k = i+1; k < n; ++k)
        w[k] = a[k][i];

    w_[i+1] = norm(w, n);

    for (int k = 0; k < n; ++k)
        N[k] = w[k] - w_[k];

    normalize(N, n);

    double **h; matrixAlloc(h, n, n);
    double **n1; matrixAlloc(n1, n, 1);
    double **n1T; matrixAlloc(n1T, 1, n);

    for (int k = 0; k < n; ++k)
    {
        n1[k][0] = N[k];
        n1T[0][k] = N[k];
    }


    double **aux = matrixMultiply(n1, n1T, n, 1, 1, n); 

    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++){
            if(k==j) h[k][j] = 1  - 2*aux[k][j];
            else h[k][j] = - 2*aux[k][j];
        }


    matrixDealloc(n1,n); matrixDealloc(n1T, 1); matrixDealloc(aux, n);
    delete[] w; delete[] w_; delete[] N;
    return h;
}

double** houseHolder(double **a, double **&h, const int n){
    double **h_;
    double **hi;
    double **a_New_;
    double **a_New;
    double **a_Old;

    matrixAlloc(a_Old, n, n);

    matrixAlloc(h, n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if(i==j) h[i][j] = 1;            

    matrixAssign(a_Old, a, n, n);

    for (int i = 0; i < n-2; ++i)
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
        if(i < n-3) matrixDealloc(a_New, n); 
    }

    matrixDealloc(a_Old, n);

    return a_New;
}


int main(int argc, char const *argv[])
{
    const int n = 5;
    double A[n][n] = {{40, 8, 4, 2, 1},
                    {8, 30, 12, 6, 2},
                    {4, 12, 20, 1, 2},
                    {2, 6, 1, 25, 4},
                    {1, 2, 2, 4, 5}};

    double **a;
    matrixAlloc(a, n, n);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            a[i][j] = A[i][j];

    double **h;

    matrixPrint(a, n, n);
    double **t = houseHolder(a, h, n);

    matrixPrint(t, n, n);
    matrixPrint(h, n, n);

    return 0;
}
