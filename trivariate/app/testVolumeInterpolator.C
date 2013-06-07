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

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
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

    SplineVolume vol_in;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> vol_in;
    ofstream outfile(argv[2]);

    // Fetch Greville parameters
    BsplineBasis basis_u = vol_in.basis(0);
    int nmb_u = basis_u.numCoefs();
    BsplineBasis basis_v = vol_in.basis(1);
    int nmb_v = basis_v.numCoefs();
    BsplineBasis basis_w = vol_in.basis(2);
    int nmb_w = basis_w.numCoefs();
    vector<double> par_u(nmb_u);
    vector<double> par_v(nmb_v);
    vector<double> par_w(nmb_w);
    int ki;
    for (ki=0; ki<nmb_u; ++ki)
      par_u[ki] = basis_u.grevilleParameter(ki);
    for (ki=0; ki<nmb_v; ++ki)
      par_v[ki] = basis_v.grevilleParameter(ki);
    for (ki=0; ki<nmb_w; ++ki)
      par_w[ki] = basis_w.grevilleParameter(ki);
    
    // Evaluate the volume in a regular grid
    vector<double> points;
    vol_in.gridEvaluator(par_u, par_v, par_w, points);
   
    // Fetch weights
    vector<double> weights;
    if (vol_in.rational())
      vol_in.getWeights(weights);

    // Make volume
    shared_ptr<SplineVolume> vol_out = 
      shared_ptr<SplineVolume>(VolumeInterpolator::regularInterpolation(basis_u,
									basis_v,
									basis_w,
									par_u,
									par_v,
									par_w,
									points,
									vol_in.dimension(),
									vol_in.rational(),
									weights));

    vol_out->writeStandardHeader(outfile);
    vol_out->write(outfile);
}
