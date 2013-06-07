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

#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "newmatio.h"
#include "newmat.h"


using namespace std;


namespace Go {


//==========================================================================
void ImplicitizeSurfaceAlgo::perform()
//==========================================================================
{
    // Create barycentric coordinate system
    create_bary_coord_system3D(surf_, bc_);

    // Convert spline curve to barycentric coordinates
    SplineSurface surf_bc;
    cart_to_bary(surf_, bc_, surf_bc);

    // Check if the surface has a single patch
    bool single_patch = (surf_.order_u() == surf_.numCoefs_u()
			 && surf_.order_v() == surf_.numCoefs_v());

    // Make the matrix of numerical coefficients (the D-matrix). Any
    // vector in the nullspace of this matrix will be a solution.
    vector<vector<double> > mat;
    if (single_patch) {
	make_matrix(surf_bc, deg_, mat);
    } else {
	// The matrices from all the patches are stacked on top of each
	// other
	vector<SplineSurface> patches;
	GeometryTools::splitSurfaceIntoPatches(surf_bc, patches);
	int num = (int)patches.size();
	make_matrix(patches[0], deg_, mat);
	vector<vector<double> > tmp;
	for (int i = 1; i < num; ++i) {
	    make_matrix(patches[i], deg_, tmp);
	    mat.insert(mat.end(), tmp.begin(), tmp.end());
	}
    }

    // Find the nullspace and construct the implicit function.
    vector<double> b;
    make_implicit_svd(mat, b, sigma_min_);

    // Set the coefficients
    implicit_ = BernsteinTetrahedralPoly(deg_, b);

    return;
}


//==========================================================================


} // namespace Go
