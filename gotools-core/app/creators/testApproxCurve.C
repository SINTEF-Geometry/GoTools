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

#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <fstream>


using namespace Go;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 4) {
	MESSAGE("Usage: inputfile tolerance outputfile.");
	return 0;
    }

    // Read input arguments
    std::ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Input file not found or file corrupt");

    double epsge = atof(argv[2]);
    std::ofstream outfile(argv[3]);

    // Input surface is to be a GoSplineCurve.
    ObjectHeader header;
    header.read(infile);
    shared_ptr<SplineCurve> crv(new SplineCurve());
    crv->read(infile);

    int max_iter = 20;
    vector<shared_ptr<ParamCurve> > crvs;
    crvs.push_back(crv);
    vector<Point> start_pt, end_pt;
    try {
	double max_dist;
	crv = shared_ptr<SplineCurve>
	    (CurveCreators::approxCurves(&crvs[0], &crvs[1],
					 start_pt, end_pt, epsge, max_dist, max_iter));
	if (max_dist > epsge) {
	    MESSAGE("Failed approximating within tolerance (" << epsge <<
		       "), using cv anyway. Dist: " << max_dist);
	}
    } catch (...) {
	MESSAGE("Failed approximating input curve, returning input curve.");
    }
    crv->writeStandardHeader(outfile);
    crv->write(outfile);
}

