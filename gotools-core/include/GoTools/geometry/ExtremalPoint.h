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

#ifndef _EXTREMALPOINT_H
#define _EXTREMALPOINT_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/Point.h"


namespace Go {

/// Computation of extremal points
namespace ExtremalPoint {

  /// Most extreme points for a set of surfaces, either tensor-produc spline 
  /// surfaces, surfaces that can be converted to splines or trimmed versions 
  /// therof.
  /// Return value entities: pair of points and parameter pairs in that sequence
  void computeExtremalPoints(std::vector<shared_ptr<ParamSurface> >& surfs,
			     const Point& dir, double tol,
			     std::vector<std::pair<Point,Point> >& extremal_points);

  /// Most extreme points for one surface, either a tensor-produc spline 
  /// surface, a surface that can be converted to splines or a trimmed versions 
  /// therof.
  /// A maximum value for the constructed 1D surface is updated during computations
  /// Return value entities: pair of points and parameter pairs in that sequence
  void extremalPoints(shared_ptr<ParamSurface>& surf,
		      const Point& dir, double tol, double& maxval,
		      std::vector<std::pair<Point,Point> >& extremal_points);

} // End of namespace ExtremalPoint

} // End of namespace Go


#endif

