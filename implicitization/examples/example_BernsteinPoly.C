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




