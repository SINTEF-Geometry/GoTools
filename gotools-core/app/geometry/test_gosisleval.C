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

#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sisl.h"

using namespace Go;
using namespace std;

void sisleval(const SplineSurface& sf, std::vector<Point>& p,
	      double u, double v, int derivs)
{
    SISLSurf* ss = GoSurf2SISL(sf, false);
    double epar[2];
    epar[0] = u;
    epar[1] = v;
    int ilfs = 0;
    int ilft = 0;
    vector<double> eder(sf.dimension()*((derivs+1)*(derivs+2)/2));
    double norm[3];
    int stat = 0;
    s1421(ss, derivs, epar, &ilfs, &ilft, &eder[0], norm, &stat);
    int cur = 0;
    for (int i = 0; i <= derivs; ++i) {
	for (int j = 0; j <= i; ++j) {
	    Point pt(sf.dimension());
	    pt.setValue(&eder[0] + cur*sf.dimension());
	    p[cur] = pt;
	    ++cur;
	}
    }
}


int main(int argc, char* argv[] )
{
    DEBUG_ERROR_IF(argc < 3, "Usage: " << argv[0]
		<< " inputsurf inputpoints" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    DEBUG_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read surface from file
    ObjectHeader head;
    SplineSurface sf;
    is >> head;
    ASSERT(head.classType() == SplineSurface::classType());
    is >> sf;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    std::vector<Point> p(3, Point(sf.dimension()));
//     for (int i = 0; i < 10000; ++i) {
	for (int j = 0; j < n; ++j) {
 	    sf.point(p, pt[2*j], pt[2*j+1], 1);
// 	    sisleval(sf, p, pt[2*j], pt[2*j+1], 1);
	    cout << p[0] << p[1] << p[2] << (p[1] % p[2]);
	}
//     }
}





