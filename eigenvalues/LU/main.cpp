#include <iostream>
using namespace std;

void LU_Decomposition(double **a, double **l, double **u, const int n){

    double aux;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            if(i==j){ l[i][j] = 1; u[i][j] = 0;}
            else{
                l[i][j] = 0;
                u[i][j] = 0;
            }   
        }

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i <= j; ++i)
        {
            aux = 0.;
            for (int k = 0; k < i; ++k)
                aux += (l[i][k] * u[k][j]);
            u[i][j] = a[i][j] - aux;
        }

        for (int i = j+1; i < n; ++i)
        {
            aux = 0.;
            for (int k = 0; k < j; ++k)
                aux += (l[i][k] * u[k][j]);

            l[i][j] = (a[i][j] - aux) / u[j][j];
        }
    }

}


double* LUsolver(double **l, double **u, double *b, const int n){

    double aux;
    double *y = new double[n];
    y[0] = b[0];

    for (int i = 1; i < n; ++i)
    {
        aux = 0.;
        for (int k = 0; k < i; ++k)
        {
            aux += l[i][k]*y[k];
        }
        y[i] = b[i] - aux;
    }


    double *x = new double[n];
    x[n-1] = y[n-1]/u[n-1][n-1];

    for (int i = n-2; i >= 0; --i)
    {
        aux = 0;
        for (int k = i+1; k < n; ++k)
        {
            aux += u[i][k] * x[k];
        }
        x[i] = (y[i] - aux)/u[i][i];
    }

    delete[] y;
    return x;
}



int main(int argc, char const *argv[])
{
    const int n = 3;
    double A_[n][n] = {{1., 1., 1.}, {4., 3., -1.}, {3., 5., 3.}};

    double **a = new double*[n];
    double **l = new double*[n];
    double **u = new double*[n];


    for (int i = 0; i < n; ++i){
        a[i] = new double[n];
        l[i] = new double[n];
        u[i] = new double[n];
        for (int j = 0; j < n; ++j)
            a[i][j] = A_[i][j];
    }

    LU_Decomposition(a, l, u, n);

    double *b = new double[n];
    b[0] = 1.; b[1] = 6.; b[2] = 4.;

    double *x = LUsolver(l, u, b, n);
    cout << "{" << x[0] << "," << x[1] << "," << x[2] << "}";


    for (int i = 0; i < n; ++i) {delete[] a[i]; delete[] l[i]; delete[] u[i];}    
    delete[] a;delete[] l; delete[] u; delete[] b; delete[] x;

    return 0;
}