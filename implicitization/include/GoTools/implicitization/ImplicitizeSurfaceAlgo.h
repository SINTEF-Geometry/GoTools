//==========================================================================
//                                                                          
// File: ImplicitizeSurfaceAlgo.h                                            
//                                                                          
// Created: Wed Feb  5 16:16:07 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ImplicitizeSurfaceAlgo.h,v 1.8 2006-03-31 09:09:06 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#ifndef _IMPLICITIZESURFACEALGO_H
#define _IMPLICITIZESURFACEALGO_H


#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BaryCoordSystem.h"


namespace Go {


/**
 * Class that implements an implicitization algorithm for spline
 * surfaces
 *
 * Input: A surface of type SplineSurface, and the degree that is
 * chosen for the implicitization.
 *
 * Output: A barycentric coordinate system, and the Bernstein
 * polynomial that represents the implicit surface.
 *
 * The algorithm first defines a barycentric coordinate system in
 * terms of a tetrahedron that is slightly larger than the bounding
 * box of the spline surface. It then computes the Bernstein
 * polynomial representing the implicit surface. If the chosen degree
 * is too low to produce the exact result, an approximate implicit
 * surface is produced. If the degree is too high, the resulting
 * polynomial will be reducible and contain the exact implicitization
 * as a factor. If the input surface consists of multiple patches, the
 * output polynomial will be a product of the separate implicit
 * surfaces, or an approximation of it.
 */

class ImplicitizeSurfaceAlgo {
public:
    /// Default constructor
    ImplicitizeSurfaceAlgo() : tol_(3.0e-15) { }
    /// Constructor.
    /// \param deg degree of the implicit representation
    explicit ImplicitizeSurfaceAlgo(int deg) : deg_(deg), tol_(3.0e-15) { }
    /// Constructor.
    /// \param surf spline surface to be implicitized
    /// \param deg degree of the implicit representation
    ImplicitizeSurfaceAlgo(const SplineSurface& surf, int deg)
	: surf_(surf), deg_(deg), tol_(3.0e-15) { }

    /// Load the spline surface to be implicitized.
    /// \param surf surfaceon SplineSurface form
    void useSplineSurface(const SplineSurface& surf)
    { surf_ = surf; }

    /// Choose the degree of the implicit representation.
    /// \param deg degree
    void setDegree(int deg)
    { deg_ = deg; }

    /// Set the tolerance.
    /// \param tol tolerance; default value is 3.0e-15
    void setTolerance(double tol)
    { tol_ = tol; }

    /// Perform the implicitization.
    /// This function runs the implicitization algorithm.
    void perform();

    /// Get the result of the implicitization.
    /// \retval implicit a BernsteinTetrahedralPoly representing the
    /// implicit surface
    /// \retval bc the barycentric coordinate system in which the
    /// BernsteinTetrahedralPoly is defined
    void getResultData(BernsteinTetrahedralPoly& implicit,
		       BaryCoordSystem3D& bc, double& sigma_min)
    {
	implicit = implicit_;
	bc = bc_;
	sigma_min = sigma_min_;
    }

private:
    SplineSurface surf_;
    BernsteinTetrahedralPoly implicit_;
    BaryCoordSystem3D bc_;
    int deg_;
    double tol_;
    double sigma_min_;

};


} // namespace Go


#endif // _IMPLICITIZESURFACEALGO_H

