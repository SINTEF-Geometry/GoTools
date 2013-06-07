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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;


int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 5, "Usage: " << argv[0]
		    << " volumeinfile numpts_u numpts_v numpts_w" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read volume from file
    SplineVolume vol;
    is >> vol;

    vector<vector<double> > pars(3);

    for (int i = 0; i < 3; ++i)
      {
	double current_par = vol.startparam(i);
	int nmb = atoi(argv[2+i]);
	double step = (vol.endparam(i) - current_par)/(double)(nmb-1);
	for (int j = 0; j < nmb; ++j, current_par += step) pars[i].push_back(step);
      }

    vector<vector<double> > pts_only, pts_with_derivs, derivs_u, derivs_v, derivs_w;

    vol.computeBasisGrid(pars[0], pars[1], pars[2],
			 pts_only);
    vol.computeBasisGrid(pars[0], pars[1], pars[2],
			 pts_with_derivs, derivs_u, derivs_v, derivs_w);

    if (pts_only.size() != pts_with_derivs.size())
      {
	cout << "Test 1: Wrong base size" << endl;
	exit(-1);
      }
    for (size_t i = 0; i < pts_only.size(); ++i)
      {
	if (pts_only[i].size() != pts_with_derivs[i].size())
	  {
	    cout << "Test 1: Wrong base size at position " << i << endl;
	    exit(-1);
	  }
	for (size_t j = 0; j < pts_only[i].size(); ++j)
	  if (pts_only[i][j] != pts_with_derivs[i][j])
	    {
	      cout << "Test 1: Different entries at (" << i << "," << j << "): " << pts_only[i][j] << " and " << pts_with_derivs[i][j] << endl;
	      exit(-1);
	    }
      }

    cout << "Test 1 complete" << endl;

    vector<BasisPts> bp_pts_only;
    vector<BasisDerivs> bp_pts_with_derivs;
    vol.computeBasisGrid(pars[0], pars[1], pars[2],
			 bp_pts_only);
    vol.computeBasisGrid(pars[0], pars[1], pars[2],
			 bp_pts_with_derivs);

    if (bp_pts_only.size() != bp_pts_with_derivs.size())
      {
	cout << "Test 2: Wrong base size" << endl;
	exit(-1);
      }

    for (size_t i = 0; i < bp_pts_only.size(); ++i)
      {
	BasisPts pt1 = bp_pts_only[i];
	BasisDerivs pt2 = bp_pts_with_derivs[i];

	for (int j = 0; j < 3; ++j)
	  if (pt1.param[j] != pt2.param[j])
	    {
	      cout << "Test 2: Wrong param value at position " << i << "," << j << endl;
	      exit(-1);
	    }

	for (int j = 0; j < 3; ++j)
	  if (pt1.left_idx[j] != pt2.left_idx[j])
	    {
	      cout << "Test 2: Wrong left_idx value at position " << i << "," << j << endl;
	      exit(-1);
	    }

	if (pt1.basisValues.size() != pt2.basisValues.size())
	  {
	    cout << "Test 2: Wrong base size at position " << i << endl;
	    exit(-1);
	  }

	for (size_t j = 0; j < pt1.basisValues.size(); ++j)
	  if (pt1.basisValues[j] != pt2.basisValues[j])
	    {
	      cout << "Test 2: Different entries at (" << i << "," << j << "): " << pt1.basisValues[j] << " and " << pt2.basisValues[j] << endl;
	      exit(-1);
	    }
      }

    cout << "Test 2 complete" << endl;

}
