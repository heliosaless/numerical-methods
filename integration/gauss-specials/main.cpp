#include <iostream>
#include <math.h> 
using namespace std;


double f(double x){
	return pow(sin(2*x) + 4*x*x + 3*x, 2);
}

double gh2point(double xi, double xf){
	double fx1 = f(-0.70710678118);
	double fx2 = f(0.70710678118);

	double w1 = 0.88622692545;
	double w2 = 0.88622692545;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;

	return i;
}

double gh3point(double xi, double xf){
	double fx1 = f(-1.22474487139);
	double fx2 = f(0);
	double fx3 = f(1.22474487139);

	double w1 = 0.29540897515;
	double w2 = 1.1816359006;
	double w3 = 0.29540897515;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;

	return i;
}

double gh4point(double xi, double xf){
	//Poly 4points: -2(-8x^4+24x^2-6)
	//wk at root x: (8*4!*sqrt(pi))/((4^2)*((2(4x^3-6x))^2))
	double fx1 = f(-1.6507);
	double fx2 = f(-0.52465);
	double fx3 = f(0.52465);
	double fx4 = f(1.6507);

	double w1 = 0.0813022;
	double w2 = 0.80491;
	double w3 = 0.80491;
	double w4 = 0.0813022;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;
	i += fx4*w4;

	return i;
}

double gl2point(double xi, double xf){
	double fx1 = f(-0.70710678118);
	double fx2 = f(3.41421356237);

	double w1 = 0.85355339059;
	double w2 = 0.1464466094;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;

	return i;
}

double gl3point(double xi, double xf){
	double fx1 = f(0.4157745568);
	double fx2 = f(2.2942803603);
	double fx3 = f(6.2899450829);

	double w1 = 0.7110930099;
	double w2 = 0.2785177336;
	double w3 = 0.0103892565;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;

	return i;
}

double gl4point(double xi, double xf){
	//x^4-16x^3+72x^2-96x+24/4!
	double fx1 = f(0.32255);
	double fx2 = f(1.7458);
	double fx3 = f(4.5366);
	double fx4 = f(9.3951);

	double w1 = 0.603115;
	double w2 = 0.357347;
	double w3 = 0.0388895;
	double w4 = 0.00053928;


	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;
	i += fx4*w4;

	return i;
}

double gt2point(double xi, double xf){
	double fx1 = f(-0.70710678118);
	double fx2 = f(0.70710678118);

	double w1 = 1.57079632679;
	double w2 = 1.57079632679;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;

	return i;
}

double gt3point(double xi, double xf){
	double fx1 = f(-0.86602540378);
	double fx2 = f(0);
	double fx3 = f(0.86602540378);

	double w1 = 1.0471975512;
	double w2 = 1.0471975512;
	double w3 = 1.0471975512;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;

	return i;
}

double gt4point(double xi, double xf){
	double fx1 = f(-0.92387953251);
	double fx2 = f(-0.38268343236);
	double fx3 = f(0.38268343236);
	double fx4 = f(0.92387953251);

	double w1 = 0.78539816339;
	double w2 = 0.78539816339;
	double w3 = 0.78539816339;
	double w4 = 0.78539816339;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;
	i += fx4*w4;


	return i;
}



double integrateH(double xi, double xf, int points){
	if(points == 2) return gh2point(xi,xf);
	if(points == 3) return gh3point(xi,xf);
	if(points == 4) return gh4point(xi,xf);
}

double integrateL(double xi, double xf, int points){
	if(points == 2) return gl2point(xi,xf);
	if(points == 3) return gl3point(xi,xf);
	if(points == 4) return gl4point(xi,xf);
}

double integrateT(double xi, double xf, int points){
	if(points == 2) return gt2point(xi,xf);
	if(points == 3) return gt3point(xi,xf);
	if(points == 4) return gt4point(xi,xf);
}
 

int main(int argc, char const *argv[])
{
	int inf = 1000;
	int minf = -1000;
	int zero = 0;
	int a = -1;
	int b = 1;

	//integrate e^(-x^2)*(sin(2*x) + 4*x^2 + 3*x)^2 from -inf to inf 
	cout << "Wolfram: " << 34.0278 << "----" << endl;
	cout << "2point Gauss-Hermit: " << integrateH(minf, inf,2) << endl;
	cout << "3point Gauss-Hermit: " << integrateH(minf, inf,3) << endl;
	cout << "4point Gauss-Hermit: " << integrateH(minf, inf,4) << endl;


	//integrate e^(-x)*(sin(2*x) + 4*x^2 + 3*x)^2 from 0 to inf
	cout << "Wolfram: " << 547.17 <<  "----" <<endl;
	cout << "2point Gauss-Laguerre: " << integrateL(zero, inf, 2) << endl;
	cout << "3point Gauss-Laguerre: " << integrateL(zero, inf, 3) << endl;
	cout << "4point Gauss-Laguerre: " << integrateL(zero, inf, 4) << endl;

	//integrate ((1)/(sqrt(1-x^2)))*(sin(2*x) + 4*x^2 + 3*x)^2 from -1 to 1
	cout << "Wolfram: " << 46.0524 << "----" << endl;
	cout << "2point Gauss-Chebyshev: " << integrateT(a,b, 2) << endl;
	cout << "3point Gauss-Chebyshev: " << integrateT(a,b, 3) << endl;
	cout << "4point Gauss-Chebyshev: " << integrateT(a,b, 4) << endl;


	return 0;
}