#include <iostream>
#include <math.h> 
using namespace std;

/*
*	Nome: Helio Matheus Sales Silva
*	Matr : 400800
*	
	
	Acho que o as 4derivadas não estão funcionando muito bem :/
*/

double function(double x){
	return pow(sin(2*x) + 4*x*x + 3*x, 2);
	//return 1/(sqrt(x));
}

double closed_first(double xi, double xf, double delx){
	double h = delx;
	return h*(function(xi) + function(xf))/2;
}

double closed_second(double xi, double xf, double delx){
	double h = delx/2;
	return h*(function(xi) + 4*function(xi+h) + function(xi+2*h))/3;
}


double closed_third(double xi, double xf, double delx){
	double h = delx/3;
	return 3*h*(function(xi) + 3*function(xi+h) + 3*function(xi+2*h) + function(xf))/8;
}

double closed_fourth(double xi, double xf, double delx){
	double h = delx/4;

	//return h*( 896*function(xi)/120.0 - 1472*function(xi+h)/120.0
	//		+ 2528*function(xi+2*h)/120.0 -832*function(xi+3*h)/360.0
	//		+ 448*function(xi+4*h)/60.0  )  / 24.0;
	return 2*h*(7*function(xi) + 32*function(xi+h)+ 12*function(xi+2*h)+ 32*function(xi+3*h)+ 7*function(xi+4*h))/45;
}

double open_first(double xi, double xf, double delx){
	double h = delx/3;
	return 3*h*(function(xi+h) + function(xf-h))/2;
}

double open_second(double xi, double xf, double delx){
	double h = delx/4;
	return 4*h*(2*function(xi+h) - function(xi+2*h) + 2*function(xf-h))/3;
}

double open_third(double xi, double xf, double delx){
	double h = delx/5;
	return 5*h*(11*function(xi+h) + function(xi+2*h) + function(xf-2*h) + 11*function(xf-h))/24;
}

double open_fourth(double xi, double xf, double delx){
	double h = delx/6;
	//return h*( 9504*function(xi+h)/120.0 - 17853*function(xi+2*h)/120.0
			//+ 26961*function(xi+3*h)/120.0 -15693*function(xi+4*h)/360.0
			//+ 4734*function(xf-h)/60.0  )  / 24.0;

	return 6*h*(11*function(xi+h) - 14*function(xi+2*h) + 26*function(xi+3*h) - 14*function(xi+4*h) + 11*function(xi+5*h))/20;
}

double integrate_(double xi, double xf, double delx, bool closed, int points){
	if(closed){
		if(points == 2) return closed_first(xi, xf, delx);
		else if(points == 3) return closed_second(xi, xf, delx);
		else if(points == 4) return closed_third(xi, xf, delx);
		else if(points == 5) return closed_fourth(xi, xf, delx);
		else {return 0.;}
	}else{
		if(points == 2) return open_first(xi, xf, delx);
		else if(points == 3) return open_second(xi, xf, delx);
		else if(points == 4) return open_third(xi, xf, delx);
		else if(points == 5) return open_fourth(xi, xf, delx);
		else {return 0.;}
	}
}

struct info{
	double it;
	double value;
};
 
struct info integrate(double a, double b, double eplison, bool closed, int points){
	double old_ = 0, new_ = 0, error = 10000;
	int n = 0;

	while(error >= eplison){
		n = n + 1;
		double delx = (b-a)/(double)n;
		new_ = 0;

		for (int k = 0; k < n; k++){
			double xi = a + k*delx;
			double xf = xi + delx;
			new_ += integrate_(xi, xf, delx, closed, points);
		}

		if(old_ != 0) error = fabs(new_ - old_)/old_;
		old_ = new_;
	}

	struct info r;
	r.it = n;
	r.value = new_;
	return r;
}

void print(double a, double b, int points, bool closed, double eplison){
	struct info r = integrate(a, b, eplison, closed, points);
	cout << "Integrate from " << a << " to " << b << " with " << points << " points, using";
	if (closed) cout <<  " closed";
	else cout << " open";
	cout << " formula: " <<  r.value << " with " << r.it << " iterations." << endl;
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


	return 0;
}