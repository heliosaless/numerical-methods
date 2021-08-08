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
	//return 1/sqrt(x);
	return pow(x,-0.66666666666);
	//return 1/(sqrt(4-pow(x,2)));
}

double f(double x, int exp, double a, double b){
	if(exp == 0) return function(x);
	if(exp == 1) return function(xExp(x, a, b))*dxExp(x, a, b);
	if(exp == 2) return function(xExpExp(x, a, b))*dxExpExp(x, a, b);
}	

double closed_first(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx;
	return h*(f(xi, exp, a, b) + f(xf, exp, a, b))/2;
}

double closed_second(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/2;
	return h*(f(xi, exp, a, b) + 4*f(xi+h, exp, a, b) + f(xi+2*h, exp, a, b))/3;
}


double closed_third(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/3;
	return 3*h*(f(xi, exp, a, b) + 3*f(xi+h, exp, a, b) + 3*f(xi+2*h, exp, a, b) + f(xf, exp, a, b))/8;
}

double closed_fourth(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/4;

	//return h*( 896*f(xi)/120.0 - 1472*f(xi+h)/120.0
	//		+ 2528*f(xi+2*h)/120.0 -832*f(xi+3*h)/360.0
	//		+ 448*f(xi+4*h)/60.0  )  / 24.0;
	return 2*h*(7*f(xi, exp, a, b) + 32*f(xi+h, exp, a, b)+ 12*f(xi+2*h, exp, a, b)+ 32*f(xi+3*h, exp, a, b)+ 7*f(xi+4*h, exp, a, b))/45;
}

double open_first(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/3;
	return 3*h*(f(xi+h, exp, a, b) + f(xf-h, exp, a, b))/2;
}

double open_second(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/4;
	return 4*h*(2*f(xi+h, exp, a, b) - f(xi+2*h, exp, a, b) + 2*f(xf-h, exp, a, b))/3;
}

double open_third(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/5;
	return 5*h*(11*f(xi+h, exp, a, b) + f(xi+2*h, exp, a, b) + f(xf-2*h, exp, a, b) + 11*f(xf-h, exp, a, b))/24;
}

double open_fourth(double xi, double xf, double delx, int exp, double a, double b){
	double h = delx/6;
	//return h*( 9504*f(xi+h)/120.0 - 17853*f(xi+2*h)/120.0
			//+ 26961*f(xi+3*h)/120.0 -15693*f(xi+4*h)/360.0
			//+ 4734*f(xf-h)/60.0  )  / 24.0;

	return 6*h*(11*f(xi+h, exp, a, b) - 14*f(xi+2*h, exp, a, b) + 26*f(xi+3*h, exp, a, b) - 14*f(xi+4*h, exp, a, b) + 11*f(xi+5*h, exp, a, b))/20;
}

struct info{
	double it;
	double value;
};

double integrate_(double xi, double xf, double delx, bool closed, int points, int exp, double a, double b){
	if(closed){
		if(points == 2) return closed_first(xi, xf, delx, exp, a, b);
		else if(points == 3) return closed_second(xi, xf, delx, exp, a, b);
		else if(points == 4) return closed_third(xi, xf, delx, exp, a, b);
		else if(points == 5) return closed_fourth(xi, xf, delx, exp, a, b);
		else {return 0.;}
	}else{
		if(points == 2) return open_first(xi, xf, delx, exp, a, b);
		else if(points == 3) return open_second(xi, xf, delx, exp, a, b);
		else if(points == 4) return open_third(xi, xf, delx, exp, a, b);
		else if(points == 5) return open_fourth(xi, xf, delx, exp, a, b);
		else {return 0.;}
	}
}


 
struct info integrate(double a, double b, int points, bool closed, double eplison, int exp=0, double a_ = 0, double b_ = 0){
	double old_ = 0, new_ = 0, error = 10000;
	int n = 0;

	while(error >= eplison){
		n = n + 1;
		double delx = (b-a)/(double)n;
		new_ = 0;

		for (int k = 0; k < n; k++){
			double xi = a + k*delx;
			double xf = xi + delx;
			new_ += integrate_(xi, xf, delx, closed, points, exp, a_, b_);
		}

		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}

	struct info r;
	r.it = n;
	r.value = new_;
	return r;
}

struct info runExp(double a, double b,  int points, bool closed, int exp, double eplison1, double eplison2){
	double c = 1;
	double old_ = 0, new_ = 0;
	double error = 10000;
	while(error >= eplison1){
		new_ = integrate(-c, c, points, closed, eplison2, exp, a, b).value;

		c = c + 0.01;
		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}

	struct info r;
	r.it = c;
	r.value = new_; 
	return r;
}


void print(double a, double b, int points, bool closed, double eplison){
	struct info r = integrate(a, b, points, closed, eplison);
	cout << "Integrate from " << a << " to " << b << " with " << points << " points, using";
	if (closed) cout <<  " closed";
	else cout << " open";
	cout << " formula: " <<  r.value << " with " << r.it << " iterations." << endl;
}

void printExp(double a, double b, int points, bool closed, int exp, double eplison1, double eplison2){
	struct info r = runExp(a, b, points, closed, exp, eplison1, eplison2);
	cout << "Integrate from " << a << " to " << b << " with " << points << " points, using";
	if (closed) cout <<  " closed";
	else cout << " open";
	cout << " formula and";
	if(exp == 1) cout << " Simple Exponential";
	else{
		if(exp == 2) cout << " Double Exponential";
		else {cout << "\n\n\n Use print instead.! \n";}
	} 
	cout << ": " <<  r.value << " with c = " << r.it << endl;
}

int main(int argc, char const *argv[])
{
	bool closed = true;
	int points = 2;
	double a = 0;
	double b = 1;
	double eplison = 10e-6;

	print(a,b, points, closed, eplison);
	print(a,b, points, !closed, eplison);
	print(a,b, points+1, closed, eplison);
	print(a,b, points+1, !closed, eplison);
	print(a,b, points+2, closed, eplison);
	print(a,b, points+2, !closed, eplison);
	print(a,b, points+3, closed, eplison);
	print(a,b, points+3, !closed, eplison);


	printExp(a,b, points, closed, 1, eplison, eplison);
	printExp(a,b, points, closed, 2, eplison, eplison);
	printExp(a,b, points, !closed, 1, eplison, eplison);
	printExp(a,b, points, !closed, 2, eplison, eplison);

	printExp(a,b, points+3, closed, 1, eplison, eplison);
	printExp(a,b, points+3, closed, 2, eplison, eplison);
	printExp(a,b, points+3, !closed, 1, eplison, eplison);
	printExp(a,b, points+3, !closed, 2, eplison, eplison);



	return 0;
}