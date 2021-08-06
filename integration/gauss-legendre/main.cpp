#include <iostream>
#include <math.h> 
#define pi 	3.14159265358979323846

using namespace std;


double xExp(double s, double a, double b){
	return ((a+b) + (b-a)*tanh(s))/2;
}

double xExpExp(double s, double a, double b){
	return ((a+b) + (b-a)*tanh((pi*sinh(s))/2))/2;
}

double dxExp(double s, double a, double b){
	return (b-a)/(2*pow(cosh(s), 2));
}


double dxExpExp(double s, double a, double b){
	return ((b-a)*pi*cosh(s))/(4*pow(cosh((pi*sinh(s))/2),2));
}

double function(double x){
	//return pow(sin(2*x) + 4*x*x + 3*x, 2);
	return 1/sqrt(x);
	//return pow(x,-0.66666666666);
	//return 1/(sqrt(4-pow(x,2)));
}

double f(double x, int exp, double a, double b){
	if(exp == 0) return function(x);
	if(exp == 1) return function(xExp(x, a, b))*dxExp(x, a, b);
	if(exp == 2) return function(xExpExp(x, a, b))*dxExpExp(x, a, b);
}

double x(double alpha, double xi, double xf){
	return ((xi+xf)/2) + ((xf-xi)*alpha/2);
}

double gl2point(double xi, double xf, int exp, double a, double b){
	double alpha1 = -0.57735026919;
	double alpha2 = 0.57735026919;
	double fx1 = f(x(alpha1, xi, xf), exp, a, b );
	double fx2 = f(x(alpha2, xi, xf), exp, a, b );

	double w1 = 1;
	double w2 = 1;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	return ((xf-xi)/2) * i;
}

double gl3point(double xi, double xf, int exp, double a, double b){
	double alpha1 = -0.77459666924;
	double alpha2 = 0;
	double alpha3 = 0.77459666924;
	double fx1 = f(x(alpha1, xi, xf), exp, a, b );
	double fx2 = f(x(alpha2, xi, xf), exp, a, b );
	double fx3 = f(x(alpha3, xi, xf), exp, a, b );
	double w1 = 0.55555555555;
	double w2 = 0.88888888888;
	double w3 = 0.55555555555;

	double i = 0;

	i += fx1*w1;
	i += fx2*w2;
	i += fx3*w3;

	return ((xf-xi)/2) * i;
}

double gl4point(double xi, double xf, int exp, double a, double b){
	double alpha1 = -0.86113631159;
	double alpha2 = -0.33998104358;
	double alpha3 = 0.33998104358;
	double alpha4 = 0.86113631159;
	double fx1 = f(x(alpha1, xi, xf), exp, a, b);
	double fx2 = f(x(alpha2, xi, xf), exp, a, b);
	double fx3 = f(x(alpha3, xi, xf), exp, a, b);
	double fx4 = f(x(alpha3, xi, xf), exp, a, b);
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

double integrate_(double xi, double xf, int points, int exp, double a, double b){
	if(points == 2) return gl2point(xi,xf, exp, a, b);
	if(points == 3) return gl3point(xi,xf, exp, a, b);
	if(points == 4) return gl4point(xi,xf, exp, a, b);
}

struct info{
	double it;
	double value;
};
 
struct info integrate(double a, double b, int points, double eplison, int exp=0, double a_ = 0, double b_ = 0){
	double old_ = 0, new_ = 0, error = 10000;
	int n = 0;

	while(error >= eplison){
		n = n + 1;
		double delx = (b-a)/(double)n;
		new_ = 0;

		for (int k = 0; k < n; k++){
			double xi = a + k*delx;
			double xf = xi + delx;
			new_ += integrate_(xi, xf, points, exp, a_, b_);
		}

		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}
	struct info r;
	r.it = n;
	r.value = new_; 
	return r;
}


struct info runExp(double a, double b,  int points, int exp, double eplison1, double eplison2){
	double c = 1;
	double old_ = 0, new_ = 0;
	double error = 10000;
	while(error >= eplison1){
		new_ = integrate(-c, c, points,eplison2, exp, a, b).value;

		c = c + 0.1;
		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}

	struct info r;
	r.it = c;
	r.value = new_; 
	return r;
}

void print(double a, double b, int points, double eplison){
	struct info r = integrate(a, b, points, eplison);
	cout << "Integrate from " << a << " to " << b << " with " << points << " points: ";
	cout << r.value << " with " << r.it << " iterations." << endl;
}

void printExp(double a, double b, int points, int exp, double eplison1, double eplison2){
	struct info r = runExp(a, b, points, exp, eplison1, eplison2);
	cout << "Integrate from " << a << " to " << b << " with " << points << " points, using";
	if(exp == 1) cout << " Simple Exponential";
	else{
		if(exp == 2) cout << " Double Exponential";
		else {cout << "\n\n\n Use print instead.! \n";}
	} 
	cout << ": " <<  r.value << " with c = " << r.it << endl;
}

int main(int argc, char const *argv[])
{

	double a = 0;
	double b = 1;
	int points = 2;
	double eplison1 = 10e-6;
	double eplison2 = 10e-6;
	int exp = 1;


	print(a,b, points, eplison1);
	print(a,b, points+1, eplison1);
	print(a,b, points+2, eplison1);
	printExp(a,b, points, exp, eplison1, eplison2);
	printExp(a,b, points+1, exp, eplison1, eplison2);
	printExp(a,b, points+2, exp, eplison1, eplison2);
	printExp(a,b, points, exp+1, eplison1, eplison2);
	printExp(a,b, points+1, exp+1, eplison1, eplison2);
	printExp(a,b, points+2, exp+1, eplison1, eplison2);


	return 0;
}