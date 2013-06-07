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

#ifndef _SPLINEAPPROXIMATOR_H
#define _SPLINEAPPROXIMATOR_H


#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/BsplineBasis.h"

namespace Go
{


    /** An Interpolator that generates a spline curve approximating
     *  the given dataset in the least squares sense.
     */
class SplineApproximator : public Interpolator
{
public:
    /// Constructor takes no arguments
    SplineApproximator()
	: num_coefs_(0)
    {
    }

    /// Virtual destructor ensures safe inheritance
    virtual ~SplineApproximator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// The interpolating function, as inherited by \ref Interpolator.  
    /// Prior to calling this function, the user must have specified:
    /// - \em either the number of control points to use in the 
    ///   approximating curve by \ref setNumCoefs().
    /// - \em or/and directly specified the BsplineBasis to use
    ///   by \ref setSplineSpace().
    /// A default BsplineBasis will be generated if the user has only set
    /// the number of control points prior to calling interpolate().
    /// For parameter list, see \ref Interpolator.
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);
    
    /// Specify the number of basis functions / control points to use
    /// in the approximating curve.
    void setNumCoefs(int num) {
	num_coefs_ = num;
    }

    /// Directly specify the spline space in which to search for the 
    /// approximating function.
    void setSplineSpace(const BsplineBasis& basis) {
	basis_ = basis;
    }


private:
    int num_coefs_;
    BsplineBasis basis_;
};


} // namespace Go


#endif // _SPLINEAPPROXIMATOR_H

