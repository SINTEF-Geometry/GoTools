//==========================================================================
//                                                                          
// File: integration_test.C                                                  
//                                                                          
// Created: Thu Mar 21 13:24:12 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: integration_test.C,v 1.5 2006-02-15 12:38:54 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/utils/Integration.h"
#include <iostream>


using namespace std;
using namespace Go;


class Expo {
public:
    typedef double first_argument_type;
    typedef double second_argument_type;
    typedef double result_type;

    Expo() {}

    double operator() (double t) const
    {
	return exp(-t*t);
    }

    double operator() (double u, double v) const
    {
	return exp(-u*u-v*v);
    }
};



int main()
{
    Expo e;
    double xa = -2.0;
    double xb = 2.0;
//     double ya = -2.0;
//     double yb = 2.0;

    const int n = 1000;

    double s;
    for (int i = 0; i < n; ++i)
 	s = simpsons_rule(e, xa, xb);
//  	s = gaussian_quadrature(e, xa, xb);
//  	s = simpsons_rule2D(e, xa, xb, ya, yb);
// 	s = gaussian_quadrature2D(e, xa, xb, ya, yb);

    cout << "Integral = " << s << endl;

    return 0;
}
