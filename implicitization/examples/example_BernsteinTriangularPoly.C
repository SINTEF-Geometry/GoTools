//==========================================================================
//                                                                          
// File: example_BernsteinTriangularPoly.C
//                                                                          
// Created: Fri Nov 29 17:34:39 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_BernsteinTriangularPoly.C,v 1.17 2007-05-09 17:57:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include <vector>
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    cout << "*** BernsteinTriangularPoly ***" << endl;
    cout << endl;

    // Construct from vector<double>
    double arr[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    int num = sizeof(arr)/sizeof(double);
    vector<double> coefs(arr, arr+num);
    BernsteinTriangularPoly b(2, coefs);

    cout << "*** Access functions ***" << endl;
    // degree()
    int deg = b.degree();
    cout << deg << endl;
    // operator[]
    double co = b[4];
    cout << co << endl;
    const BernsteinTriangularPoly c = b;
    double cco = c[3];
    cout << cco << endl;
    cout << endl;

    // operator()
    cout << "*** operator() ***" << endl;
    // Construct a linear polynomial
    int deg1 = 1;
    double c1[] = { 0.0, 1.0, 1.0 };
    BernsteinTriangularPoly p1(deg1, c1, c1+3);
    cout << "p1:\n" << p1 << endl;
    Vector3D u(0.3333, 0.3333, 0.3333);
    cout << u << endl << p1(u) << endl;
    cout << endl;

    // norm() and normalize()
    cout << "*** norm() and normalize() ***" << endl;
    BernsteinTriangularPoly h = b;
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    h.normalize();
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    cout << endl;

     // deriv()
    cout << "*** deriv() ***" << endl;
    Vector3D d(-1.0, 0.5, 0.5);
    BernsteinTriangularPoly p1der;
    p1.deriv(1, d, p1der);
    cout << "p1':\n" << p1der << endl;
    cout << endl;

    // blossom()
    cout << "*** blossom() ***" << endl;
    int degh = h.degree();
    Vector3D e1(1.0, 0.0, 0.0);
    Vector3D e2(0.0, 1.0, 0.0);
    Vector3D e3(0.0, 0.0, 1.0);
    for (int i = degh; i >= 0; --i) {
	for (int j = degh - i; j >= 0; --j) {
	    int k = degh - i - j;
	    vector<Vector3D> uvec(i, e1);
	    vector<Vector3D> v2(j, e2);
	    uvec.insert(uvec.end(), v2.begin(), v2.end());
	    vector<Vector3D> v3(k, e3);
	    uvec.insert(uvec.end(), v3.begin(), v3.end());
	    cout << h.blossom(uvec) << endl;
	}
    }
    Vector3D mid(0.3333, 0.3333, 0.3333);
    vector<Vector3D> uvec(p1.degree(), mid);
    double x = p1.blossom(uvec);
    cout << "p1(0.3333, 0.3333, 0.3333, 0.3333) = " << x << endl;
    cout << endl;

    // operator*=
    cout << "*** operator*= ***" << endl;
    // b2 is the square of the polynomial b.
    BernsteinTriangularPoly b2 = b * b;
    cout << "b2:\n" << b2 << endl;
    // Check multiplication for a number of points in the triangle
    Array<double, 3> uvw(0.1, 0.1, 0.8);
    cout << "(b*b)(u)" << "\t" << "b(u)*b(u)" << endl;
    cout << b2(uvw) << "\t" << b(uvw) * b(uvw) << endl;
    uvw = Array<double, 3>(0.5, 0.5, 0.0);
    cout << b2(uvw) << "\t" << b(uvw) * b(uvw) << endl;
    uvw = Array<double, 3>(0.33, 0.33, 0.0);
    cout << b2(uvw) << "\t" << b(uvw) * b(uvw) << endl;
    cout << endl;

    // read() and write()
    cout << "*** read() and write() ***" << endl;
    ifstream in("data/bernstein_triangular_poly.dat");
    BernsteinTriangularPoly f;
    in >> f;
    cout << f << endl;
    cout << endl;

    return 0;
}
