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

