/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
