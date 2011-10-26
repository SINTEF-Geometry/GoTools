//==========================================================================
//                                                                          
// File: example_BernsteinMulti.C
//                                                                          
// Created: Thu Jun 13 13:34:42 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: example_BernsteinMulti.C,v 1.6 2007-04-17 13:18:16 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/BernsteinMulti.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/utils/Array.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    cout << "*** BernsteinMulti ***" << endl;
    cout << endl;

    // Construct from vector
    double arr[] = { -0.2, -1.2, -5.2, -1.2, -0.2,
		    -0.2, 2.8, 2.8, 2.8, -0.2,
		    -0.2, -1.2, -5.2, -1.2, -0.2 };
    int num = sizeof(arr)/sizeof(double);
    vector<double> coefs(arr, arr+num);
    BernsteinMulti f(4, 2, coefs);
    cout << "f:\n" << f << endl;
 
    cout << "*** Access functions ***" << endl;
    // degreeU() and degreeV()
    int degu = f.degreeU();
    int degv = f.degreeV();
    cout << degu << '\t' << degv << endl;
    // operator[]
    double co = f[4]; // co = -0.2
    cout << co << endl;
    const BernsteinMulti g = f;
    double cco = f[3]; // cco = -1.2
    cout << cco << endl;
    // coefsBegin() and coefsEnd()
    typedef vector<double>::iterator iter;
    typedef vector<double>::const_iterator const_iter;
    iter it = f.coefsBegin();
    iter jt = f.coefsEnd();
    const_iter cit = g.coefsBegin();
    const_iter cjt = g.coefsEnd();
    cout << endl;

    // operator()
    cout << "*** operator() ***" << endl;
    // Construct useful unit polynomials u and v
    vector<double> coefsuv(2, 0.0);
    coefsuv[1] = 1.0;
    BernsteinMulti u(1, 0, coefsuv);
    BernsteinMulti v(0, 1, coefsuv);
    BernsteinMulti h = (u-0.25)*(u-0.25) + (v-0.5)*(v-0.5) - 0.1;
    cout << "h = (u-0.25)*(u-0.25) + (v-0.5)*(v-0.5) - 0.1" << endl;
    cout << "h:\n" << h << endl;
    double x = h(0.5, 0.5); // x = -0.0375
    cout << "h(0.5, 0.5) = " << x << endl;
    cout << endl;

    // isZero() etc.
    cout << "*** is Zero() etc. ***" << endl;
    vector<double> zero(36, 0.0);
    BernsteinMulti z(5, 5, zero);
    cout << "z.isZero() = " << z.isZero() << endl;
    double del = 1.0e-14;
    vector<double> pos(36, del);
    BernsteinMulti p(5, 5, pos);
    cout << "p.isStrictlyPositive() = "
	 << p.isStrictlyPositive(0.5 * del) << endl;
    BernsteinMulti n = -1.0 * p;
    cout << "n.isStrictlyNegative() = "
	 << n.isStrictlyNegative(0.5 * del) << endl;
    BernsteinMulti nn = -0.5 * p;
    cout << "nn.isNonNegative() = "
	 << nn.isNonNegative(del) << endl;
    cout << endl;

    // norm(), normalize(), and mean()
    cout << "*** norm(), normalize(), and mean() ***" << endl;
    BernsteinMulti k = g;
    cout << "k:\n" << k << endl
	 << "k.norm() = " << k.norm() << endl;
    k.normalize();
    cout << "k:\n" << k << endl
	 << "k.norm() = " << k.norm() << endl;
    cout << "k.mean() = " << k.mean() << endl;
    cout << endl;

    // deriv()
    cout << "*** deriv() ***" << endl;
    cout << "g:\n" << g << endl;
    cout << "g_u:\n" << g.deriv(1, 0) << endl;
    cout << "g_v:\n" << g.deriv(0, 1) << endl;
    cout << "g_uv:\n" << g.deriv(1, 1) << endl;
    cout << endl;

    // detHess() and traceHess()
    cout << "*** detHess() and traceHess() ***" << endl;
    cout << "g:\n" << g << endl;
    cout << "g.detHess():\n" << g.detHess() << endl;
    cout << "g.traceHess():\n" << g.traceHess() << endl;
    cout << endl;

    // blossom()
    cout << "*** blossom() ***" << endl;
    int du = h.degreeU();
    int dv = h.degreeV();
    vector<double> vvec(du, 0.0);
    for (int j = dv; j >= 0; --j) {
	vector<double> uvec(dv, 0.0);
	for (int i = du; i >= 0; --i) {
	    cout << h.blossom(uvec, vvec) << endl;
	    if (i != 0) {
		uvec[i-1] = 1.0;
	    }
	}
	if (j != 0) {
	    vvec[j-1] = 1.0;
	}
    }
    vector<double> uvec(h.degreeU(), 0.5);
    vvec = vector<double>(h.degreeV(), 0.5);
    x = h.blossom(uvec, vvec);
    cout << "h(0.5, 0.5) = " << x << endl;
    cout << endl;

    // pickDomain()
    cout << "*** pickDomain() ***" << endl;
    BernsteinMulti fll = f.pickDomain(0.0, 0.5, 0.0, 0.5);
    cout << "f: [0.0, 0.5] x [0.0, 0.5]\n" << fll << endl;
    cout << endl;

    // pickLine()
    cout << "*** pickLine() ***" << endl;
    Array<double, 2> a(0.0, 0.0);
    Array<double, 2> b(0.5, 0.5);
    BernsteinPoly line = f.pickLine(a, b);
    cout << "line: (0.0, 0.0) - (0.5, 0.5)\n" << line << endl;
    cout << endl;

    // degreeElevate()
    cout << "*** degreeElevate() ***" << endl;
    k.degreeElevate(1, 1);
    cout << "k:\n" << k << endl;
    cout << endl;

    // bindU() and bindV()
    cout << "*** bindU() and bindV() ***" << endl;
    BernsteinPoly bu = k.bindU(0.0);
    cout << "bu:\n" << bu << endl;
    BernsteinPoly bv = k.bindV(1.0);
    cout << "bv:\n" << bv << endl;
    cout << endl;

    // read() and write()
    cout << "*** read() and write() ***" << endl;
    ifstream in("data/bernstein_multi.dat");
    in >> f;
    cout << f << endl;
    cout << endl;

    return 0;
}




