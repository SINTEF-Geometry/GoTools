//===========================================================================
//                                                                           
// File: ImplicitUtils.h                                                
//                                                                           
// Created: Tue Apr 24 15:08:42 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ImplicitUtils.h,v 1.44 2006-03-31 09:09:05 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

