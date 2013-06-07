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

#ifndef _EVALCURVESET_H_
#define _EVALCURVESET_H_


#include "GoTools/utils/Point.h"

#include <vector>

namespace Go
{

/// This abstract class provides an interface to a set of curves that can
/// be evaluated.
/// Representing the actual geometry, typically used when iteratively
/// approximating the set of curves on the same basis.
class EvalCurveSet
{
public:
    /// virtual destructor ensures safe inheritance
    virtual ~EvalCurveSet()
    {}
  
    /// Evaluate the curves.
    /// \param t parameter in which to evaluate.
    /// \return the evaluated points for the curve set.
    virtual std::vector<Point> eval(double t)=0;

    /// Evaluate the curve derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n number of derivatives to compute.
    /// \param der the evaluated points up to the n'th derivative for the curve set.
    virtual void eval(double t, int n, std::vector<std::vector<Point> >& der)=0; // n = order of diff

    /// Start parameter of domain.
    /// \return start parameter of the spline space.
    virtual double start()=0;

    /// End parameter of domain.
    /// \return end parameter of the spline space.
    virtual double end()=0;

    /// The geometric dimension of the spline curves.
    /// \return geometric dimension of the space.
    virtual int dim() = 0;

    /// Whether the approximation is within tolerances in input parameter.
    /// \param par parameter in which to evaluate.
    /// \param approxpos whether the input points are within tolerance from the
    ///                  evaluated points (as given by eval()).
    /// \param tol1 tolerance used to decide approximation accuracy.
    /// \param tol2 tolerance used to decide approximation accuracy.
    /// \return whether the approximation is within tolerances in input parameter.
    virtual bool approximationOK(double par, const std::vector<Point>& approxpos,
				 double tol1, double tol2)=0;

    /// The number of curves in the curve set.
    /// \return the number of curves in the curve set.
    virtual int nmbCvs() = 0;

    /// Reset intermediate error. 
    virtual void resetErr()
    {
      ; // Nothing to do. Overruled when required
    }

};

} // namespace Go

#endif // _EVALCURVESET_H_
