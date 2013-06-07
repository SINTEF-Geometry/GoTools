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

#include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
void GeometryTools::findDominant(const SplineSurface& surface,
		  Vector3D& dominant_u, Vector3D& dominant_v)
//==========================================================================
{
    int nu = surface.numCoefs_u();
    int nv = surface.numCoefs_v();
    vector<double>::const_iterator start = surface.coefs_begin();
    Vector3D temp;
    // Dominant in u-direction
    dominant_u = Vector3D(0.0, 0.0, 0.0);
    for (int j = 0; j < nv; ++j) {
	for (int dd = 0; dd < 3; ++dd) {
	    temp[dd] = *(start + 3*(nu*j + (nu-1)) + dd)
		- *(start + 3*(nu*j) + dd);
	}
	dominant_u += temp;
    }
    // Dominant in v-direction
    dominant_v = Vector3D(0.0, 0.0, 0.0);
    for (int i = 0; i < nu; ++i) {
	for (int dd = 0; dd < 3; ++dd) {
	    temp[dd] = *(start + 3*(nu*(nv-1) + i) + dd)
		- *(start + 3*i + dd);
	}
	dominant_v += temp;
    }

    return;
}


//==========================================================================
bool GeometryTools::negativeProj(const SplineSurface& surface,
		  const Array<Vector3D, 2>& refvector,
		  const double eps)
//==========================================================================
{
    int num_u = surface.numCoefs_u();
    int num_v = surface.numCoefs_v();
    Vector3D temp;
    int i = 0, j = 0;
    while (i < num_u-1) {
	j = 0;
	while (j < num_v) {
	    temp[0] = *(surface.coefs_begin() + 3*(num_u*j + i+1))
		- *(surface.coefs_begin() + 3*(num_u*j + i));
	    temp[1] = *(surface.coefs_begin() + 3*(num_u*j + i+1) + 1)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 1);
	    temp[2] = *(surface.coefs_begin() + 3*(num_u*j + i+1) + 2)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 2);
	    // Positive tolerance means that there must be a small
	    // _nonzero_ negative projection before it is reported as
	    // negative!
	    if (temp * refvector[0] < -eps)
		return true;
	    ++j;
	}
	++i;
    }
    i = 0;
    while (i < num_u) {
	j = 0;
	while (j < num_v-1) {
	    temp[0] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i))
		- *(surface.coefs_begin() + 3*(num_u*j + i));
	    temp[1] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i) + 1)
		- *(surface.coefs_begin() + 3*(num_u*j + i) +1);
	    temp[2] = *(surface.coefs_begin() + 3*(num_u*(j+1) + i) + 2)
		- *(surface.coefs_begin() + 3*(num_u*j + i) + 2);
	    // Positive tolerance means that there must be a small
	    // _nonzero_ negative projection before it is reported as
	    // negative!
	    if (temp * refvector[1] < -eps)
		return true;
	    ++j;
	}
	++i;
    }

    return false;
}



//==========================================================================


} // namespace Go
