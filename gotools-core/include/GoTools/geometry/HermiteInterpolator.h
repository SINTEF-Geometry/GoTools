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

#ifndef _HERMITEINTERPOLATOR_H
#define _HERMITEINTERPOLATOR_H


#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{

/** An Interpolator that generates a hermite spline
 *  curve through given points and tangents.
 */

class GO_API HermiteInterpolator : public Interpolator
{
 public:
    /// Constructor takes no arguments
    HermiteInterpolator() {}
 
    /// Virtual destructor enables safe inheritance
    virtual ~HermiteInterpolator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// Hermite interpolation of a sequence of points with
    /// associated tangents and parameter values
    /// \param num_points number of points to interpolate
    /// \param dimension dimension of points to interpolate (2D, 3D, etc..)
    /// \param param_start pointer to the start of the array where
    ///                    the parameter values of the points are stored.
    ///                    This should be a strictly increasing sequence
    ///                    of 'num_points' values.
    /// \param data_start pointer to the start of the array where the
    ///                   points and tangents to be interpolated are stored.
    ///                   Each point and tangent consist of 'dimension'
    ///                   coordinates, and each tangent is stored immediately
    ///                   after its corresponding point.
    /// \retval coefs The control points of the computed hermite interpolation
    ///               curve will be returned in this vector. (Use the basis() 
    ///               function to get the associated b-spline basis).
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);

    /// Hermite interpolation of a sequence of points with 
    /// associated tangents and parameter values.
    /// \param data This vector contains the set of points and 
    ///             tangents to be interpolated.  The size of the
    ///             vector is twice the total number of points.
    ///             Each entry on the form [2i] represents a point,
    ///             and the entry [2i+1] represents the associated
    ///             tangent.
    /// \param param This vector represent the parameterization of the
    ///              points, and should have one entry per point.
    ///              (Making it half the size of the 'data' vector).
    ///              The parameter sequence must be strictly increaing.
    /// \retval coefs The control points of the computed hermite 
    ///               interpolation curve will be returned in this 
    ///               vector.  (Use the basis() function to get the 
    ///               associated b-spline basis.
    void interpolate(const std::vector<Point>& data,
		     const std::vector<double>& param,
		     std::vector<double>& coefs);

 private:
    BsplineBasis basis_;
}; 

} // namespace Go

#endif // _HERMITEINTERPOLATOR_H


 
