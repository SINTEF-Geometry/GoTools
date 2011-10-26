//==========================================================================
//                                                                          
// File: ImplicitizeCurveAndVectorAlgo.h                                     
//                                                                          
// Created: Tue Jun  7 13:28:00 2005                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: 
// $Id: ImplicitizeCurveAndVectorAlgo.h,v 1.5 2006-03-31 09:09:06 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _IMPLICITIZECURVEANDVECTORALGO_H
#define _IMPLICITIZECURVEANDVECTORALGO_H


#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/BaryCoordSystem.h"


namespace Go {


/**
 * Class that implements an algorithm for finding an implicit ruled surface
 * defined by a curve and a vector.
 *
 * Input: A curve of type SplineCurve, a vector of type Point, and the
 * degree that is chosen for the implicit surface.
 *
 * Output: A barycentric coordinate system, and the Bernstein
 * polynomial that represents the implicit surface.
 *
 * The algorithm first defines a barycentric coordinate system in
 * terms of a tetrahedron that is slightly larger than the bounding
 * box of the spline curve. It then computes the Bernstein polynomial
 * representing the implicit ruled surface defined by the spline curve
 * and vector. If the chosen degree is too low to produce the exact
 * result, an approximate implicit surface is produced. If the degree
 * is too high, the resulting polynomial will be reducible and contain
 * the exact implicitization as a factor. If the input curve consists
 * of multiple segments, the output polynomial will be a product of
 * the separate implicit surfaces, or an approximation of it.
 */

class ImplicitizeCurveAndVectorAlgo {
public:
    /// Default constructor
    ImplicitizeCurveAndVectorAlgo() : tol_(3.0e-15) { }
    /// Constructor
    /// \param deg degree of the implicit representation
    explicit ImplicitizeCurveAndVectorAlgo(int deg) 
	: deg_(deg), tol_(3.0e-15) { }
    /// Constructor
    /// \param crv 3D spline curve defining the implicit surface
    /// \param pt 3D vector defining the implicit surface
    /// \param deg degree of the implicit surface
    ImplicitizeCurveAndVectorAlgo(const SplineCurve& crv,
				  const Point& pt, int deg)
	: crv_(crv), pt_(pt), deg_(deg), tol_(3.0e-15) { }

    /// Load the spline curve for the implicitization
    /// \param crv curve on SplineCurve form
    void useSplineCurve(const SplineCurve& crv)
    { crv_ = crv; }

    /// Load the vector for the implicitization.
    /// \param pt vector of type Point
    void useVector(const Point& pt)
    { pt_ = pt; }

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
    int perform();

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
    SplineCurve crv_;
    Point pt_;
    BernsteinTetrahedralPoly implicit_;
    BaryCoordSystem3D bc_;
    int deg_;
    double tol_;
    double sigma_min_;

};


} // namespace Go


#endif // _IMPLICITIZECURVEANDVECTORALGO_H

