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
#include "GoTools/geometry/SplineSurface.h"
#include <vector>

using namespace std;

int main()
{
    // Read data from standard input
    int dim, numpt1, numpt2;
    cin >> dim >> numpt1 >> numpt2;
    std::vector<double> param1(numpt1);
    std::vector<double> param2(numpt2);
    std::vector<double> data(numpt1*numpt2 * dim);
    for (int i = 0; i < numpt1; ++i)
	cin >> param1[i];
     for (int i = 0; i < numpt2; ++i)
	cin >> param2[i];
    for (int i = 0; i < numpt1*numpt2; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    cin >> data[i*dim + dd];
	}
    }

    // Interpolate and print out
    Go::SplineInterpolator interpol1;
    interpol1.setNaturalConditions();
    Go::SplineInterpolator interpol2;
    interpol2.setNaturalConditions();
//      Go::SplineApproximator interpol1;
//      interpol1.setNumCoefs(4);
//      Go::SplineApproximator interpol2;
//      interpol2.setNumCoefs(4);

    Go::SplineSurface sf;
    sf.interpolate(interpol1, interpol2,
		   numpt1, numpt2, dim,
		   &param1[0], &param2[0], &data[0]);

    cout.precision(15);
    sf.writeStandardHeader(cout);
    cout << sf;
}
