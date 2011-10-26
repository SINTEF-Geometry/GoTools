//===========================================================================
//                                                                           
// File: example_BernsteinPoly.C
//                                                                           
// Created: Fri May  4 14:04:20 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: example_BernsteinPoly.C,v 1.6 2007-04-17 13:18:16 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/implicitization/BernsteinPoly.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    cout << "*** BernsteinPoly ***" << endl;
    cout << endl;

    // Construct from vector<double>
    double arr[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    int num = sizeof(arr)/sizeof(double);
    vector<double> coefs(arr, arr+num);
    BernsteinPoly b(coefs);

    cout << "*** Access functions ***" << endl;
    // degree()
    int deg = b.degree();
    cout << deg << endl;
    // operator[]
    double co = b[4];
    cout << co << endl;
    const BernsteinPoly c = b;
    double cco = c[3];
    cout << cco << endl;
    // coefsBegin() and coefsEnd()
    typedef vector<double>::iterator iter;
    typedef vector<double>::const_iterator const_iter;
    iter it = b.coefsBegin();
    iter jt = b.coefsEnd();
    const_iter cit = c.coefsBegin();
    const_iter cjt = c.coefsEnd();
    cout << endl;

    // operator()
    cout << "*** operator() ***" << endl;
    // Construct useful unit polynomial t
    double coefst[] = { 0.0, 1.0 };
    BernsteinPoly t(coefst, coefst+2);
    BernsteinPoly g = (t-0.5)*(t-0.5);
    cout << "g = (t-0.5)^2" << endl;
    double x = g(0.3);
    cout << "g(0.3) = " << x << endl;
    cout << endl;

    // isZero() and isStrictlyPositive()
    cout << "*** is Zero() and isStrictlyPositive() ***" << endl;
    vector<double> zero(5, 0.0);
    BernsteinPoly z(zero);
    cout << "z.isZero() = " << z.isZero() << endl
	 << "z.isStrictlyPositive() = " << z.isStrictlyPositive() << endl;
    double del = 1.0e-14;
    z[3] = del;
    cout << "z.isZero() = " << z.isZero() << endl;
    vector<double> pos(5, del);
    BernsteinPoly p(pos);
    p[3] = 0.5 * del;
    cout << "p.isStrictlyPositive() = " << p.isStrictlyPositive(del) << endl;
    cout << endl;

    // norm() and normalize()
    cout << "*** norm() and normalize() ***" << endl;
    BernsteinPoly h = g;
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    h.normalize();
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    cout << endl;

     // deriv()
    cout << "*** deriv() ***" << endl;
    cout << "g':\n" << g.deriv(1) << endl;
    cout << "g'':\n" << g.deriv(2) << endl;
    cout << endl;

    // integral()
    cout << "*** integral() ***" << endl;
    cout << "g.integral() = " << g.integral() << endl
	 << "g.integral(0.0, 0.5) = " << g.integral(0.0, 0.5) << endl;
    cout << endl;

    // blossom()
    cout << "*** blossom() ***" << endl;
    BernsteinPoly f = t*(t-1.0);
    deg = f.degree();
    vector<double> vec(deg, 0.0);
    for (int i = deg; i >= 0; --i) {
	cout << f.blossom(vec) << endl;
	if (i != 0) {
	    vec[i-1] = 1.0;
	}
    }
    vector<double> tvec(g.degree(), 0.3);
    x = g.blossom(tvec);
    cout << "g(0.3) = " << x << endl;
    cout << endl;

    // pickInterval()
    cout << "*** pickInterval() ***" << endl;
    cout << "g: [0, 1]\n" << g << endl;
    cout << "g: [0.25, 0.75]\n" << g.pickInterval(0.25, 0.75) << endl;
    cout << endl;

    // degreeElevate()
    cout << "*** degreeElevate() ***" << endl;
    g.degreeElevate(1);
    cout << g << endl;
    cout << endl;

    // read() and write()
    cout << "*** read() and write() ***" << endl;
    ifstream in("data/bernstein_poly.dat");
    in >> f;
    cout << f << endl;
    cout << endl;

    return 0;
}




