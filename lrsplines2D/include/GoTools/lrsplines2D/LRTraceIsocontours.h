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

#ifndef _LRTRACEISOCONTOURS_H
#define _LRTRACEISOCONTOURS_H


#include <vector>
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/TraceContoursTypedefs.h"
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{
  /// Compute the level-set curves of the LR spline function 'lrs'.
  /// The input surface is split into a number of tensor-product spline patches,
  /// using the function LRSplineSurface::subdivideIntoSimpler,
  /// level curves are computed for every patch and the resulting
  /// curves are merged across patch boundaries.
  // In other words, compute the intersections between a set of horizontal planes and
  // 'lrs', interpreted as a 2 1/2 D surface. The height of the planes
  // (level-set values) are provided in the vector 'isovals'.  The output is a
  // vector with one entry per entry in 'isovals', containing the curves that
  // represents this level-set.  These are 2D curves in the parameter plane of
  // 'lrs'.  If 'include_3D_curves' is set to 'true', the corresponding 3D
  // curves will also be returned.  These will have constant z-value (equal to
  // the corresponding entry in 'isovals'), but are useful for visualization
  // purposes in a 3D viewer.  By default, the function employs the function
  // 'traceIsovals' to march out the curves.  If desired, SISL routine s1314 can
  // be used instead by setting 'use_sisl_marching' to true.  The SISL routine
  // is slower.
  std::vector<CurveVec> LRTraceIsocontours(const shared_ptr<ParamSurface>& surf,
					   const std::vector<double>& isovals,
					   const int threshold_missing,
					   const double tol); 

  std::vector<CurveVec> LRTraceIsocontours(const LRSplineSurface& lrs,
					   const std::vector<double>& isovals,
					   const int threshold_missing,
					   const double tol = 1e-6, 
					   const bool include_3D_curves = false,
					   const bool use_sisl_marching = false,
					   const CurveBoundedDomain* domain = NULL);
}; // end namespace Go;

#endif
