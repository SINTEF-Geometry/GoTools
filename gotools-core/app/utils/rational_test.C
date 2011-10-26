//===========================================================================
//                                                                           
// File: rational_test.C                                                     
//                                                                           
// Created: Wed Dec 11 09:05:34 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: rational_test.C,v 1.1 2006-01-03 11:23:07 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/utils/Rational.h"

using namespace Go;
using namespace std;

int main()
{
    Rational r1(-20160,40320);
    cout << r1 << endl;
    r1.simplify();
    cout << r1 << endl;
    Rational r2(16,6);
    r2.simplify();
    cout << r2 << ' ' << r2*r1 << ' ' << r2/r1
	 << ' ' << r2+r1 << ' ' << r2-r1 << endl;
}
