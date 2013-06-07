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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>


using namespace Go;
using std::ifstream;
using std::ofstream;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 3) {
	std::cout << "Usage:  input file, output file"
		  << std::endl;
	return -1;
    }

    SplineCurve cv_in;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> cv_in;
    ofstream outfile(argv[2]);

    // Fetch Greville parameters
    BsplineBasis basis = cv_in.basis();
    int nmb = basis.numCoefs();
    vector<double> par(nmb);
    int ki;
    for (ki=0; ki<nmb; ++ki)
      par[ki] = basis.grevilleParameter(ki);
    
    // Evaluate the curve in a regular grid
    vector<double> points;
    cv_in.gridEvaluator(points, par);
   
    // Fetch weights
    vector<double> weights;
    if (cv_in.rational())
      cv_in.getWeights(weights);

    // Make surface
    shared_ptr<SplineCurve> cv_out = 
      shared_ptr<SplineCurve>(CurveInterpolator::regularInterpolation(basis,
								      par,
								      points,
								      cv_in.dimension(),
								      cv_in.rational(),
								      weights));

    cv_out->writeStandardHeader(outfile);
    cv_out->write(outfile);
}
