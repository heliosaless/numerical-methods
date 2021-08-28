#include "../headers/powerMethod.h"
#include "../headers/matrix.h"

using namespace std;

void powerMethod(double **A, double *vo, double &lamb_new, const int n, const double eplison)
{
    double *vnew = new double[n];
    assign(vnew, vo, n);
    double *vold = new double[n];
    double lamb_old = 0.;
    lamb_new = 1.;

    while (abs(lamb_new - lamb_old) / lamb_new > eplison)
    {
        lamb_old = lamb_new;

        normalize(vnew, n);
        assign(vold, vnew, n);

        for (int i = 0; i < n; ++i)
        {
            double sum = 0.;
            for (int j = 0; j < n; ++j)
            {
                sum += A[i][j] * vold[j];
            }
            vnew[i] = sum;
        }

        lamb_new = dotProduct(vold, vnew, n);
    }

    assign(vo, vold, n);
    delete[] vnew;
    delete[] vold;
}
/*
int main(int argc, char const *argv[])
{
    const int n = 5;

    double **A = new double *[n];
    for (int i = 0; i < n; ++i)
        A[i] = new double[n];

    double vo[n] = {1, 1, 1, 1, 1};
    double lambd;
    //double A_[n][n] = {{2 / 3., 1 / 3., 1 / 3.}, {1 / 3., 4 / 3., 1 / 3.}, {1 / 3., 1 / 3., 2 / 3.}};
    //double A_[n][n] = {{5, 2, 1}, {2, 3, 1}, {1, 1, 2}};
    double A_[n][n] = {{40, 8, 4, 2, 1}, {8, 30, 12, 6, 2}, {4, 12, 20, 1, 2}, {2, 6, 1, 25, 4}, {1, 2, 2, 4, 5}};

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = A_[i][j];

    powerMethod(A, vo, lambd, n, 10e-7);

    cout << lambd << endl;
    cout << "{";
    for (int i = 0; i < n; ++i)
        if (i != n - 1)
            cout << vo[i] << ",";
        else
            cout << vo[i];
    cout << "}" << endl;

    return 0;
}
*/