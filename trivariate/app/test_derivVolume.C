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

    ALWAYS_ERROR_IF(argc != 8 && argc != 10, "Usage: " << argv[0]
		    << " volumeinfile paru parv parw deru derv derw (swap pardir1, swap pardir2)" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read volume from file
    SplineVolume vol;
    ObjectHeader head;
    is >> head >> vol;

    double pu = atof(argv[2]);
    double pv = atof(argv[3]);
    double pw = atof(argv[4]);
    int deru = atoi(argv[5]);
    int derv = atoi(argv[6]);
    int derw = atoi(argv[7]);

    int dtot = deru + derv + derw;
    int dvwtot = derv + derw;

    int dir1=-1, dir2=-1;
    if (argc > 8)
    {
	dir1 = atoi(argv[8]);
	dir2 = atoi(argv[9]);
    }

    vector<Point> pts;
    pts.resize(((dtot+3)*(dtot+2)*(dtot+1))/6);
    vol.point(pts, pu, pv, pw, dtot);

    cout << "From point evaluation in original surface: " << pts[((dtot+2)*(dtot+1)*dtot)/6 + ((dvwtot+1)*dvwtot)/2 + derw] << endl;

  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      vol.swapParameterDirection(dir1, dir2);

    SplineVolume* dervol = vol.derivVolume(deru, derv, derw);
  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      dervol->swapParameterDirection(dir2, dir1);

    Point pderiv;
    dervol->point(pderiv, pu, pv, pw);

    cout << "From point evaluation in derivative surface: " << pderiv << endl;

    ofstream os("data/derivOut.g2");
    dervol->writeStandardHeader(os);
    dervol->write(os);

}
