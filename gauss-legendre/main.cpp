#include <iostream>
#include <math.h> 
using namespace std;

/*
*	Nome: Helio Matheus Sales Silva
*	Matr : 400800
*	

*/

double f(double x){
	return pow(sin(2*x) + 4*x*x + 3*x, 2);
}

double x(double alpha, double xi, double xf){
	return ((xi+xf)/2) + ((xf-xi)*alpha/2);
}

double gl2point(double xi, double xf){
	double alpha1 = -0.57735026919;
	double alpha2 = 0.57735026919;
	double fx1 = f(x(alpha1, xi, xf));
	double fx2 = f(x(alpha2, xi, xf));

	double w1 = 1;
	double w2 = 1;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	return ((xf-xi)/2) * i;
}

double gl3point(double xi, double xf){
	double alpha1 = -0.77459666924;
	double alpha2 = 0;
	double alpha3 = 0.77459666924;
	double fx1 = f(x(alpha1, xi, xf));
	double fx2 = f(x(alpha2, xi, xf));
	double fx3 = f(x(alpha3, xi, xf));
	double w1 = 0.55555555555;
	double w2 = 0.88888888888;
	double w3 = 0.55555555555;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;

	return ((xf-xi)/2) * i;
}

double gl4point(double xi, double xf){
	double alpha1 = -0.86113631159;
	double alpha2 = -0.33998104358;
	double alpha3 = 0.33998104358;
	double alpha4 = 0.86113631159;
	double fx1 = f(x(alpha1, xi, xf));
	double fx2 = f(x(alpha2, xi, xf));
	double fx3 = f(x(alpha3, xi, xf));
	double fx4 = f(x(alpha3, xi, xf));
	double w1 = 0.3478548451374539;
	double w2 = 0.652145;
	double w3 = 0.652145;
	double w4 = 0.3478548451374539;

	double i = 0.0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;
	i += fx4*w4;

	return ((xf-xi)/2) * i;
}

double integrate(double xi, double xf, int points){
	if(points == 2) return gl2point(xi,xf);
	if(points == 3) return gl3point(xi,xf);
	if(points == 4) return gl4point(xi,xf);
}
 
double run(int a, int b, double eplison, int points){
	double old_ = 0, new_ = 0, error = 10000;
	int n = 0;

	while(error >= eplison){
		n = n + 1;
		double delx = (b-a)/(double)n;
		new_ = 0;

		for (int k = 0; k < n; k++){
			double xi = a + k*delx;
			double xf = xi + delx;
			new_ += integrate(xi, xf, points);
		}

		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}

	cout <<"Usando " << points << " pontos: Temos " << n << " iteracoes e " << "obtemos valor final " << new_ << endl;
}

int main(int argc, char const *argv[])
{
	int a = 0;
	int b = 1;

	//cout << integrate(a,b,2) << endl;
	//cout << integrate(a,b,3) << endl;
	//cout << integrate(a,b,4) << endl;

	run(a,b, 10e-6, 2);
	run(a,b, 10e-6, 3);
	run(a,b, 10e-6, 4);


	return 0;
}