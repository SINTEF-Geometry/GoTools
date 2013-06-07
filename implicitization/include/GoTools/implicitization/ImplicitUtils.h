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

#ifndef _IMPLICITUTILS_H
#define _IMPLICITUTILS_H


#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/geometry/PointCloud.h"
#include <vector>


namespace Go {


class BernsteinPoly;
class BernsteinMulti;
class BernsteinTriangularPoly;
class BernsteinTetrahedralPoly;
class SplineCurve;
class SplineSurface;
class BoundingBox;


/// Creates a barycentric coordinate system from a given spline curve
void create_bary_coord_system2D(const SplineCurve& curve,
				BaryCoordSystem2D& bc);

/// Creates a barycentric coordinate system from a given 3D spline curve
void create_bary_coord_system3D(const SplineCurve& curve,
				BaryCoordSystem3D& bc);

/// Creates a barycentric coordinate system from a given spline surface
void create_bary_coord_system3D(const SplineSurface& surface,
				BaryCoordSystem3D& bc);

/// Creates a barycentric coordinate system from a point cloud
void create_bary_coord_system3D(const PointCloud3D& cloud,
				BaryCoordSystem3D& bc);

/// Creates a barycentric coordinate system from a bounding box
void create_bary_coord_system3D(const BoundingBox& box,
				BaryCoordSystem3D& bc);

/// Creates a new curve with control points in 3 barycentric coordinates
/// from the 2D input curve.
void cart_to_bary(const SplineCurve& cv, const BaryCoordSystem2D& bc,
		  SplineCurve& cv_bc);

/// Creates a new curve with control points in 4 barycentric coordinates
/// from the 3D input curve.
void cart_to_bary(const SplineCurve& cv, const BaryCoordSystem3D& bc,
		  SplineCurve& cv_bc);

/// Creates a new surface with control points in 4 barycentric coordinates
/// from the 3D input surface.
void cart_to_bary(const SplineSurface& sf, const BaryCoordSystem3D& bc,
		  SplineSurface& sf_bc);

/// Creates a new point cloud in 4 barycentric coordinates
/// from the 3D input cloud.
void cart_to_bary(const PointCloud3D& cloud, const BaryCoordSystem3D& bc,
		  PointCloud4D& cloud_bc);

/// Make the matrix D
void make_matrix(const SplineCurve& curve, int deg,
		 std::vector<std::vector<double> >& mat);

/// Make the matrix D
void make_matrix(const SplineSurface& surf, int deg,
		 std::vector<std::vector<double> >& mat);

/// Make the matrix D
void make_matrix(const PointCloud4D& cloud, int deg,
		 std::vector<std::vector<double> >& mat);

/// Performs implicitization using SVD. This method is suitable when
/// the implicitization is approximate. If the implicitization is
/// exact, make_implicit_gauss() is better. Based on the function
/// SVD() from the newmat matrix library.
void make_implicit_svd(std::vector<std::vector<double> >& mat, 
		       std::vector<double>& b, double& sigma_min);

/// Performs implicitization using Gaussian elimination. This method
/// is suitable when the implicitization is exact. If the
/// implicitization is approximate, make_implicit_svd() is better.
void make_implicit_gauss(std::vector<std::vector<double> >& mat,
			 std::vector<double>& b);


} // namespace Go


#endif // _IMPLICITUTILS_H

