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

#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Utils.h"

#include <fstream>


using namespace Go;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 6) {
	MESSAGE("Usage: input point file,  approximation tolerance, degree, initial no of coefs, output curve file.");
	return 0;
    }

    // Read input arguments
    std::ifstream filein(argv[1]);
    ALWAYS_ERROR_IF(filein.bad(), "Input file not found or file corrupt");
    double tol = atof(argv[2]);
    int degree = atoi(argv[3]);
    int num_coef = atoi(argv[4]);
    std::ofstream fileout(argv[5]);
    int dim = 3;

    // Input surface should be a a PointCloud.
    ObjectHeader header;
    header.read(filein);
    PointCloud3D pointseq;
    pointseq.read(filein);
    double* raw_data = pointseq.rawData();
    int num_pts = pointseq.numPoints();
    vector<double> points(raw_data, raw_data+3*num_pts);

    // Parameterize
    vector<double> param(num_pts);
    param[0] = 0.0;
    for (int ka=1; ka<num_pts; ++ka)
      param[ka] = param[ka-1] + sqrt(Utils::distance_squared(&points[3*(ka-1)], &points[3*ka], &points[3*ka]));

    ApproxCurve approx(points, param, 3, tol, num_coef, degree+1);
    double maxdist, avdist;
    int max_iter = 8;
    shared_ptr<SplineCurve> crv = approx.getApproxCurve(maxdist, avdist, max_iter);
    std::cout << "Maxdist: " << maxdist << ", avdist: " << avdist << std::endl;

    crv->writeStandardHeader(fileout);
    crv->write(fileout);
}

