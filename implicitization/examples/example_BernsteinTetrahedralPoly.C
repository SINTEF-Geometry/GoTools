//==========================================================================
//                                                                          
// File: example_BernsteinTetrahedralPoly.C
//                                                                          
// Created: Tue Jan 21 10:15:32 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_BernsteinTetrahedralPoly.C,v 1.8 2007-05-09 17:57:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include <vector>
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    cout << "*** BernsteinTetrahedralPoly ***" << endl;
    cout << endl;

    // Construct from vector<double>
    double arr[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
    int num = sizeof(arr)/sizeof(double);
    vector<double> coefs(arr, arr+num);
    BernsteinTetrahedralPoly b(2, coefs);

    cout << "*** Access functions ***" << endl;
    // degree()
    int deg = b.degree();
    cout << deg << endl;
    // operator[]
    double co = b[4];
    cout << co << endl;
    const BernsteinTetrahedralPoly c = b;
    double cco = c[3];
    cout << cco << endl;
    cout << endl;

    // operator()
    cout << "*** operator() ***" << endl;
    // Construct a linear polynomial
    int deg1 = 1;
    double c1[] = { 0.0, 1.0, 1.0, 1.0 };
    BernsteinTetrahedralPoly p1(deg1, c1, c1+4);
    cout << "p1:\n" << p1 << endl;
    Vector4D u(0.25, 0.25, 0.25, 0.25);
    cout << u << endl << p1(u) << endl;
    cout << endl;

    // norm() and normalize()
    cout << "*** norm() and normalize() ***" << endl;
    BernsteinTetrahedralPoly h = b;
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    h.normalize();
    cout << "h:\n" << h << endl
	 << "h.norm() = " << h.norm() << endl;
    cout << endl;

     // deriv()
    cout << "*** deriv() ***" << endl;
    Vector4D d(-1.0, 0.3333, 0.3333, 0.3333);
    BernsteinTetrahedralPoly p1der;
    p1.deriv(1, d, p1der);
    cout << "p1':\n" << p1der << endl;
    cout << endl;

    // blossom()
    cout << "*** blossom() ***" << endl;
    int degh = h.degree();
    Vector4D e1(1.0, 0.0, 0.0, 0.0);
    Vector4D e2(0.0, 1.0, 0.0, 0.0);
    Vector4D e3(0.0, 0.0, 1.0, 0.0);
    Vector4D e4(0.0, 0.0, 0.0, 1.0);
    for (int i = degh; i >= 0; --i) {
	for (int j = degh - i; j >= 0; --j) {
	    for (int k = degh - i - j; k >= 0; --k) {
		int l = degh - i - j - k;
		vector<Vector4D> uvec(i, e1);
		vector<Vector4D> v2(j, e2);
		uvec.insert(uvec.end(), v2.begin(), v2.end());
		vector<Vector4D> v3(k, e3);
		uvec.insert(uvec.end(), v3.begin(), v3.end());
		vector<Vector4D> v4(l, e4);
		uvec.insert(uvec.end(), v4.begin(), v4.end());
		cout << h.blossom(uvec) << endl;
	    }
	}
    }
    Vector4D mid(0.25, 0.25, 0.25, 0.25);
    vector<Vector4D> uvec(p1.degree(), mid);
    double x = p1.blossom(uvec);
    cout << "p1(0.25, 0.25, 0.25, 0.25) = " << x << endl;
    cout << endl;

    // pickLine()
    cout << "*** pickLine() ***" << endl;
    Array<double, 4> aa(1.0, 0.0, 0.0, 0.0);
    Array<double, 4> bb(0.0, 1.0, 0.0, 0.0);
    BernsteinPoly line = h.pickLine(aa, bb);
    cout << "line: (0.0, 0.0, 0.0, 0.0) - (0.0, 1.0, 0.0, 0.0)\n"
	 << line << endl;
    cout << endl;

    // operator*=
    cout << "*** operator*= ***" << endl;
    // b2 is the square of the polynomial b.
    BernsteinTetrahedralPoly b2 = b * b;
    cout << "b2:\n" << b2 << endl;
    // Check multiplication for a number of points in the tetrahedron
    Array<double, 4> uvwz(0.1, 0.1, 0.1, 0.7);
    cout << "(b*b)(u)" << "\t" << "b(u)*b(u)" << endl;
    cout << b2(uvwz) << "\t" << b(uvwz) * b(uvwz) << endl;
    uvwz = Array<double, 4>(0.33, 0.33, 0.33, 0.01);
    cout << b2(uvwz) << "\t" << b(uvwz) * b(uvwz) << endl;
    uvwz = Array<double, 4>(0.25, 0.25, 0.25, 0.25);
    cout << b2(uvwz) << "\t" << b(uvwz) * b(uvwz) << endl;
    cout << endl;

    // read() and write()
    cout << "*** read() and write() ***" << endl;
    ifstream in("data/bernstein_tetrahedral_poly.dat");
    BernsteinTetrahedralPoly f;
    in >> f;
    cout << f << endl;
    cout << endl;

    return 0;
}
