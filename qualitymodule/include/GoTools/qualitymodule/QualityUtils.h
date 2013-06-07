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

#ifndef __CMUTILS_H
#define __CMUTILS_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveLoop.h"


namespace Go
{

namespace qualityUtils
{

  // Determine whether a surface is a sliver face or not
  // Uses different method for SplineSurface and BoundedSurface
  // A surface is a sliver suface if its maximum length m1 in one
  // parameter direction is less than a thickness given as a parameter,
  // while its minimum length in the other parameter direction is
  // more than m1 * a factor, typically 2.0
  bool isSliverFace(shared_ptr<ParamSurface>,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace(const SplineSurface& sf,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace(const BoundedSurface& sf,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace2(const BoundedSurface& sf,
		     double thickness,
		     double factor = 2.0);

  bool hasIndistinctKnots(shared_ptr<ParamSurface> surf, double tol,
			  std::vector<shared_ptr<ParamCurve> >& trim_cv_knots);

  double estimateArea(shared_ptr<ParamSurface> surf);

  double estimateLoopArea(shared_ptr<CurveLoop> loop);

}   // namespace qualityUtils

}   // namespace Go


#endif    // #ifndef __CMUTILS_H
