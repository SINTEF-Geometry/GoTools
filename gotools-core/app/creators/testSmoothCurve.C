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

#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"

#include <fstream>


using namespace Go;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	MESSAGE("Usage: inputfile outputfile.");
	return 0;
    }

    // Read input arguments
    std::ifstream filein(argv[1]);
    ALWAYS_ERROR_IF(filein.bad(), "Input file not found or file corrupt");

    std::ofstream fileout(argv[2]);

    // Input surface should be a Go::SplineCurve and a PointCloud.
    ObjectHeader header;
    header.read(filein);
    shared_ptr<SplineCurve> cv(new SplineCurve());
    cv->read(filein);
    // We then read the PointCloud. Dim is 4: par pos_x pos_y pos_z.
    header.read(filein);
    PointCloud4D pt_cl;
    pt_cl.read(filein);
    double* raw_data = pt_cl.rawData();
    int num_pts = pt_cl.numPoints();
    vector<double> pts(num_pts*3);
    vector<double> params(num_pts);
    for (int ki = 0; ki < num_pts; ++ki)
    {
	params[ki] = raw_data[ki*4];
	for (int kj = 0; kj < 3; ++kj)
	{
	    pts[ki*3+kj] = raw_data[ki*4+1+kj];
	}
    }

    const int dim = cv->dimension();
    SmoothCurve smooth_cv(dim);
    vector<int> coef_known(cv->numCoefs(), 0);
    coef_known[0] = coef_known[coef_known.size() - 1] = 1;
    for (int ki = 0; ki < dim; ++ki)
    {
	cv->coefs_begin()[ki] = pts[ki];
	cv->coefs_end()[-dim+ki] = pts[(num_pts-1)*3+ki];
    }
    smooth_cv.attach(cv, &coef_known[0]);
    double wgts[3];
    wgts[0] = 0.25;
    wgts[1] = 0.5;
    wgts[2] = 0.25;
    smooth_cv.setOptim(wgts[0], wgts[1], wgts[2]); // Typically start with 0.25, 0.5, 0.25, experiment.
    vector<double> pt_wgts(num_pts, 1.0);
    double appr_wgt = 1.0;
    smooth_cv.setLeastSquares(pts, params, pt_wgts, appr_wgt);
//    smooth_cv.setSideConstraints(); // May be used to express fixed end conditions, tangents, etc.
    shared_ptr<SplineCurve> res_cv;
    smooth_cv.equationSolve(res_cv);

    res_cv->writeStandardHeader(fileout);
    res_cv->write(fileout);
}

