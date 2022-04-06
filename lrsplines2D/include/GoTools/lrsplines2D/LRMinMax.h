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

#ifndef _LRMINMAX_H
#define _LRMINMAX_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/Point.h"

namespace Go {

  class CurveOnSurface;

/// Computation of extremal points on LR B-spline surface
namespace LRMinMax {

  /// Compute minimu and maximum points on a 1D surface
  /// \param surface the given surface
  /// \param contour_curves only extremal points found within an innermost closed
  /// contour curve (in the parameter domain of the surface) are considered to avoid
  /// reporting all local minima and maxima
  /// \param tol tolerance for approximation and domain extent
  /// \param epsge geometry tolerance 
void computeMinMaxPoints(shared_ptr<ParamSurface> surface,
			 std::vector<std::pair<shared_ptr<ParamCurve>, double> >& contour_crvs,
			 double tol, double epsge,
			 std::vector<std::pair<Point, Point> >& minpoints,
			 std::vector<std::pair<Point, Point> >& maxpoints);

  /// Compute global extremal (minimum or maximum) points on the possibly  trimmed 1D surface
  /// \param surface where to search for extremal points
  /// \param sng 1 = maximum point expected, -1 = minimum point
  /// \param epsge geomtry tolerance
  /// \param extpoints found extremal points
 int computeExtremalPoints(shared_ptr<ParamSurface> surface,
			    int sgn, double tol, double epsge,
			    std::vector<std::pair<Point, Point> >& extpoints);

} // End of namespace LRMinMax

} // End of namespace Go


#endif

