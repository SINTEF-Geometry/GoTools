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

#ifndef _MODIFYSURF_H
#define _MODIFYSURF_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

namespace Go {

  /// Functionality used to modify one or more spline surfaces with
  /// respect to given conditions
namespace ModifySurf
{
  /// Smoothing of SplineSurface. Keep a number of coeffients
  /// along the boundary fixed (1 <= fix <= 2)
  void smoothSurface(shared_ptr<SplineSurface>& surf,
		     int nmb_bd_keep);

  /// Replace one boundary curve of a SplineSurface
  /// Approximate the initial surface
  /// bd_idx = 0: umin
  ///          1: umax
  ///          2: vmin
  // /         3: vmax
  void replaceBoundary(shared_ptr<SplineSurface> surf,
		       shared_ptr<SplineCurve> curve,
		       int bd_idx, double tol);

  /// Enforce colinearity between coeffients at the common boundary between
  /// two spline surfaces. The two outer rows of coefficients are involved
  /// in the linearity constraints, but a few additional rows are modified
  /// for reasons of smoothness. All surface boundaries, except the affected
  /// one, are kept fixed.
  bool enforceCoefCoLinearity(shared_ptr<SplineSurface> sf1, int bd1, 
			      shared_ptr<SplineSurface> sf2, int bd2, 
			      double tol, 
			      std::vector<std::vector<int> >& enumeration);

  /// Enforce colinearity at vertices where 4 edge pairs meet. The two outer rows of 
  /// coefficients around the vertex are involved in the linearity constraints, 
  /// but a few additional rows are modified for reasons of smoothness.
  bool enforceVxCoefCoLinearity(std::vector<shared_ptr<SplineSurface> >& sfs, 
				std::vector<int>& vx_enum, 
				std::vector<std::pair<std::vector<int>, std::pair<int,int> > >& coef_cond,
				double tol);

} // of namespace ModifySurf

}; // end namespace Go

#endif // _MODIFYSURF_H
