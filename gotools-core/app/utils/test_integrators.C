//===========================================================================
//                                                                           
// File: test_integrators.C                                                  
//                                                                           
// Created: Mon Nov 21 14:59:30 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_integrators.C,v 1.2 2006-04-19 09:22:36 jbt Exp $
//                                                                           
//===========================================================================

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

struct func { double operator()(double x) { return sin(x); } };

// Midpoint
// double interval(func f, double a, double b)
// {
//     return (b - a)*f((a + b)/2.0);
// }

// Simpson
// double interval(func f, double a, double b)
// {
//     return (b - a)*(f(a) + 4.0*f((a + b)/2.0) + f(b))/6.0;
// }

// 2-pt Gauss
// double interval(func f, double a, double b)
// {
//     double gamma = sqrt(3.0)/3.0;
//     double mid = (a + b)/2.0;
//     double diff = gamma*(b-a)/2.0;
//     return (b - a)*(f(mid-diff) + f(mid+diff))/2.0;
// }

// Weirdpoint
double interval(func f, double a, double b)
{
    return (b - a)*f(0.3*a + 0.7*b);
}

double composite(func f, double a, double b, int n)
{
    double h = (b - a)/double(n);
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
	total += interval(f, a + i*h, a + (i+1)*h);
    }
    return total;
}

double richardson(func f, double a, double b, int n)
{
    double coarse = composite(f, a, b, n);
    double fine = composite(f, a, b, 2*n);
    // for order 1
    return (2.0*fine - coarse)/1.0;
}

int main(int argc, char** argv)
{
    func f;
    double a = 0.0;
    double b = 1.0;
    double truth = 1 - cos(1.0);
    
    int n = 1;
    double lasterr = 0.0;
    for (int i = 0; i < 15; ++i) {
	double answer = richardson(f, a, b, n);
 	//double answer = composite(f, a, b, n);
	double err = fabs(truth - answer);
	if (i > 0) {
	    double slope = (log(err) - log(lasterr))/log(2.0);
	    cout << slope << endl;
	}
	n *= 2;
	lasterr = err;
    }
}
