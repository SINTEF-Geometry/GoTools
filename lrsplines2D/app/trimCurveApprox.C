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


#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include "GoTools/geometry/SplineCurve.h"

#include <iostream>
#include <fstream>
#include <vector>


using namespace Go;
using std::vector;
using std::cout;
using std::endl;


int main(int argc, char* argv[])
{
    if (argc != 4)
    {
	cout << "Usage: " << argv[0] << " pt_set.g2 max_iter app_cv_2d.g2" << endl;// translate_vec.g2" << endl;
	return -1;
    }

    std::ifstream filein(argv[1]);
    int max_iter(atoi(argv[2]));
    std::ofstream fileout_app_cv_2d(argv[3]);
//    std::ofstream fileout_translate(argv[4]);

    if (max_iter > 10)
    {
	cout << "Not allowing max_iter > 10." << endl;
	max_iter = 10;
    }

    const double epsgeo = 10;//5;//1e-04;
    const double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
    Point translate_vec;
    vector<double> pts_2d = TrimCrvUtils::readTrimPoints(filein, translate_vec);
    vector<vector<double> > split_pts_2d = 
      TrimCrvUtils::splitTrimPoints(pts_2d, epsgeo, kink_tol);

    // @@sbr201410 The epsgeo is very dependent on the input data.
    // For instance data with a staircase pattern needs a high tolerance, which should reflect the length of the steps.
    const int par_dim = 2;
    for (size_t ki = 0; ki < split_pts_2d.size(); ++ki)
    {
	shared_ptr<SplineCurve> spline_cv_appr_2d
	    (TrimCrvUtils::approximateTrimPts(split_pts_2d[ki], par_dim, epsgeo, max_iter));

	spline_cv_appr_2d->writeStandardHeader(fileout_app_cv_2d);
	spline_cv_appr_2d->write(fileout_app_cv_2d);
    }

}
