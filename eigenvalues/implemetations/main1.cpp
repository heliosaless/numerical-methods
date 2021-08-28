#include "../headers/inversePowerMethod.h"
#include "../headers/shiftedPowerMethod.h"
#include "../headers/powerMethod.h"

using namespace std;

void print(double lambda, double *vo, const int n){
    cout << lambda << endl;
    cout << "{";
    for (int i = 0; i < n; ++i)
        if (i != n - 1)
            cout << vo[i] << ",";
        else
            cout << vo[i];
    cout << "}" << endl;
}

int main(int argc, char const *argv[])
{
    const int n = 5;

    double **A = new double *[n];
    for (int i = 0; i < n; ++i)
        A[i] = new double[n];

    double vo1[n] = {1, 1, 1};
    double vo2[n] = {1, 1, 1};
    double vo3[n] = {1, 1, 1};
    double vo4[n] = {1, 1, 1};
    double vo5[n] = {1, 1, 1};

    double lambd1;
    double lambd2;
    double lambd3;
    double lambd4;
    double lambd5;

    //double A_[n][n] = {{2 / 3., 1 / 3., 1 / 3.}, {1 / 3., 4 / 3., 1 / 3.}, {1 / 3., 1 / 3., 2 / 3.}};
    //double A_[n][n] = {{5, 2, 1}, {2, 3, 1}, {1, 1, 2}};
    double A_[n][n] = {
        {40, 8, 4, 2, 1},
        {8, 30, 12, 6, 2},
        {4, 12, 20, 1, 2},
        {2, 6, 1, 25, 4},
        {1, 2, 2, 4, 5}};
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = A_[i][j];


    powerMethod(A, vo1, lambd1, n, 10e-6);
    print(lambd1, vo1, n);

    shiftedPowerMethod(A, vo2, lambd2, n, 10e-6, 30);
    print(lambd2, vo2, n);

    shiftedPowerMethod(A, vo3, lambd3, n, 10e-6, 20);
    print(lambd3, vo3, n);
    
    shiftedPowerMethod(A, vo4, lambd4, n, 10e-6, 10);
    print(lambd4, vo4, n);

    inversePowerMethod(A, vo5, lambd5, n, 10e-6);
    print(lambd5, vo5, n);

    return 0;
}