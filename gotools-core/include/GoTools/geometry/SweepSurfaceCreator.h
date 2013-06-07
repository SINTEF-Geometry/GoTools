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

#ifndef __SWEEPSURFACECREATOR_H
#define __SWEEPSURFACECREATOR_H

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"

namespace Go
{

  /// Functionality to create a spline surface by linear or rotational sweep
  class SweepSurfaceCreator
  {

  public:

    // Destructor
    virtual ~SweepSurfaceCreator() { }

    /// Create a linearly swept surface as a tensor product of two curves. If \f$P\f$ is the point to
    /// be swept, and the curves have B-spline coefficients \f$C_i\f$ and \f$D_j\f$, then the surface has
    /// B-spline coefficients \f$C_i+D_j-P\f$. If the point lies on the first curve,
    /// the first curve will be swept along the second curve. If the point lies on the second curve,
    /// the second curve will be swept along the first curve.
    /// \param curv1 the first curve
    /// \param curv2 the second curve
    /// \param pt the point to be swept
    /// \return a pointer to the swept surface
    static SplineSurface* linearSweptSurface(const SplineCurve &curv1,
					     const SplineCurve &curv2,
					     const Point &pt);

    /// Create a surface as the rotation of a curve around an axis
    /// \param curve the curve
    /// \param angle the rotation angle (in radians). If less than \f$2\pi\f$, it is truncated down to \f$2\pi\f$.
    ///              If more than \f$-2\pi\f$, it is truncated up to \f$-2\pi\f$.
    /// \param pt a point on the rotation axis
    /// \param axis a vector describing the direction of the rotation axis
    /// \return a pointer to the rotated surface
    static SplineSurface* rotationalSweptSurface(const SplineCurve &curve,
						 double angle,
						 const Point &pt,
						 const Point &axis);


  };    // Class SweepSurfaceCreator


} // namespace Go


#endif    // #ifndef __SWEEPSURFACECREATOR_H
