#include "../headers/matrix.h"
#include "../headers/householder.h"
using namespace std;

void matrixAlloc(double **&a, const int p1, const int p2)
{
    a = new double *[p1];
    for (int i = 0; i < p1; ++i)
    {
        a[i] = new double[p2];
        for (int j = 0; j < p2; ++j)
            a[i][j] = 0;
    }
}

void matrixDealloc(double **a, const int p1)
{
    for (int i = 0; i < p1; ++i)
        delete[] a[i];
    delete[] a;
}

void matrixAssign(double **a, double **b, const int p1, const int p2)
{
    for (int i = 0; i < p1; ++i)
        for (int j = 0; j < p2; ++j)
            a[i][j] = b[i][j];
}

void matrixPrint(double **a, const int p1, const int p2)
{
    for (int i = 0; i < p1; ++i)
    {
        for (int j = 0; j < p2; ++j)
            cout << fixed << setprecision(3) << a[i][j] << "\t ";
        cout << endl;
    }
    cout << endl;
}

double **matrixMultiply(double **a, double **b, const int p11, const int p21, const int p12, const int p22)
{
    if (p12 != p21)
    {
        cout << "multiplication error" << endl;
        return nullptr;
    }
    double **c;
    matrixAlloc(c, p11, p22);

    double sum;
    for (int i = 0; i < p11; ++i)
        for (int j = 0; j < p22; ++j)
        {
            sum = 0.;
            for (int k = 0; k < p12; ++k)
                sum += a[i][k] * b[k][j];

            c[i][j] = sum;
        }

    return c;
}

double norm(double *a, const int n)
{
    double sum = 0.;
    for (int i = 0; i < n; ++i)
    {
        sum += a[i] * a[i];
    }
    return sqrt(sum);
}

void normalize(double *a, const int n)
{
    double s = norm(a, n);
    for (int i = 0; i < n; ++i)
        a[i] /= s;
}

void assign(double *s, double *t, const int n)
{
    for (int i = 0; i < n; ++i)
        s[i] = t[i];
}

double dotProduct(double *a, double *b, const int n)
{
    double result_ = 0;
    for (int i = 0; i < n; i++)
        result_ += a[i] * b[i];
    return result_;
}
/*
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
*/
