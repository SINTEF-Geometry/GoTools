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

#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineApproximator.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>

using namespace std;

int main()
{
    // Read data from standard input
    int dim, numpt;
    cin >> dim >> numpt;
    std::vector<double> param(numpt);
    std::vector<double> data(numpt * dim);
    for (int i = 0; i < numpt; ++i) {
	param[i] = static_cast<double>(i);
	for (int dd = 0; dd < dim; ++dd) {
	    cin >> data[i*dim + dd];
	}
    }

    // Interpolate and print out
     Go::SplineInterpolator interpol;
     interpol.setFreeConditions();
//      interpol.setNaturalConditions();
//      interpol.setHermiteConditions(Go::Point(0.0, 0.0, 1.0),
//  				  Go::Point(0.0, 0.0, -1.0));

//     Go::SplineApproximator interpol;
//     interpol.setNumCoefs(static_cast<int>(floor(max(4.0, numpt*0.2))));

    Go::SplineCurve cv;
    cv.interpolate(interpol, numpt, dim, &param[0], &data[0]);
//      std::vector<double> coefs;
//      interpol.interpolate(numpt, dim, param.begin(), data.begin(), coefs);

//      Go::SplineCurve cv(interpol.basis().numCoefs(),
//  			 interpol.basis().order(),
//  			 interpol.basis().begin(),
//  			 coefs.begin(),
//  			 dim);

    cout.precision(15);
    cv.writeStandardHeader(cout);
    cout << cv;

    vector<Go::Point> p(3, Go::Point(1));
    for (int i = 0; i < numpt; ++i) {
	cv.point(p, param[i], 2);
	cout << "Pts: " << p[0] << '\n' << p[1] << '\n' << p[2] << endl;
    }
}
