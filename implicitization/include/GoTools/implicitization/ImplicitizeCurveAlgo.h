//==========================================================================
//                                                                          
// File: ImplicitizeCurveAlgo.h                                              
//                                                                          
// Created: Wed Feb  5 15:14:20 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ImplicitizeCurveAlgo.h,v 1.10 2006-03-31 09:09:05 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _IMPLICITIZECURVEALGO_H
#define _IMPLICITIZECURVEALGO_H


#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/BaryCoordSystem.h"


namespace Go {


/**
 * Class that implements an implicitization algorithm for spline curves.
 *
 * Input: A curve of type SplineCurve, and the degree that is chosen
 * for the implicit representation.
 *
 * Output: A barycentric coordinate system, and the Bernstein
 * polynomial that represents the implicit curve.
 *
 * The algorithm first defines a barycentric coordinate system in
 * terms of a triangle that is slightly larger than the bounding box
 * of the spline curve. It then computes the Bernstein polynomial
 * representing the implicit curve. If the chosen degree is too low to
 * produce the exact result, an approximate implicit curve is
 * produced. If the degree is too high, the resulting polynomial will
 * be reducible and contain the exact implicitization as a factor.  If
 * the input curve consists of multiple segments, the output
 * polynomial will be a product of the separate implicit curves, or an
 * approximation of it.
 */

class ImplicitizeCurveAlgo {
public:
    /// Default constructor
    ImplicitizeCurveAlgo() : tol_(3.0e-15) { }
    /// Constructor
    /// \param deg degree of the implicit representation
    explicit ImplicitizeCurveAlgo(int deg) : deg_(deg), tol_(3.0e-15) { }
    /// Constructor
    /// \param curve spline curve to be implicitized
    /// \param deg degree of the implicit representation
    ImplicitizeCurveAlgo(const SplineCurve& curve, int deg)
	: curve_(curve), deg_(deg), tol_(3.0e-15) { }

    /// Load the spline curve to be implicitized
    /// \param curve curve on SplineCurve form
    void useSplineCurve(const SplineCurve& curve)
    { curve_ = curve; }

    /// Choose the degree of the implicit representation
    /// \param deg degree
    void setDegree(int deg)
    { deg_ = deg; }

    /// Set the tolerance
    /// \param tol tolerance; default value is 3.0e-15
    void setTolerance(double tol)
    { tol_ = tol; }

    /// Perform the implicitization.
    /// This function runs the implicitization algorithm.
    void perform();

    /// Get the result of the implicitization
    /// \retval implicit a BernsteinTriangularPoly representing the
    /// implicit curve
    /// \retval bc the barycentric coordinate system in which the
    /// BernsteinTriangularPoly is defined
    void getResultData(BernsteinTriangularPoly& implicit,
		       BaryCoordSystem2D& bc, double& sigma_min)
    {
	implicit = implicit_;
	bc = bc_;
	sigma_min = sigma_min_;
    }

private:
    SplineCurve curve_;
    BernsteinTriangularPoly implicit_;
    BaryCoordSystem2D bc_;
    int deg_;
    double tol_;
    double sigma_min_;

};


} // namespace Go


#endif // _IMPLICITIZECURVEALGO_H

