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

#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>


using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::cin;
using namespace Go;


int main()
{
    shared_ptr<SplineCurve> sc(new SplineCurve);
    shared_ptr<SplineCurve> pc(new SplineCurve);
    shared_ptr<SplineSurface> sf(new SplineSurface);
    ObjectHeader header;
    cin >> header >> (*sc);
    cin >> header >> (*pc);
    cin >> header >> (*sf);

    CurveOnSurface cv(sf, pc, sc, true);
    const int N = 100;
    Point pt(3);
    Point ppt(2);
    Point pt2(3);
    double t0 = cv.startparam();
    double t1 = cv.endparam();
    cout << "400 1 0 4 255 255 0 255\n" << N << endl;
    double mdiff = 1e100;
    for (int i = 0; i < N; ++i) {
	double fact = i/double(N-1);
	cv.point(pt, t0*(1.0-fact) + t1*fact);
	cout << pt;
	pc->point(ppt, t0*(1.0-fact) + t1*fact);
	sf->point(pt2, ppt[0], ppt[1]);
	mdiff = min(mdiff, pt.dist(pt2));
	//cout << pt2;
    }
    cerr << mdiff << endl;
}
