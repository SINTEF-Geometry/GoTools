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

#ifndef LINDEP_UTILS_H
#define LINDEP_UTILS_H

#include "GoTools/lrsplines2D/LRSplineSurface.h"

//==============================================================================
namespace Go
//==============================================================================
{
  /// Utilities for checking an LR B-spline surface for potential linear
  /// dependency
  namespace LinDepUtils
  {
  //============================================================================
  /// Tests for potential linear dependence among the LR B-splines of a
  /// given LR spline. To be more precise: Tests whether a given LR
  /// spline is peelable. We say an LR spline is peelable if the
  /// incidence matrix contains NO non-zero entries after peeling it,
  /// i.e., if ALL the overloaded LR B-splines are peelable,
  /// cf. [Dokken, Lyche & Pettersen, 2012]. The LR spline PEELABILITY
  /// is a SUFFICIENT condition for linear INDEPENDENCE of the LR
  /// B-splines, or equivalently, the LR spline UNPEELABILITY is a
  /// NECESSARY condition for linear DEPENDENCE of the LR B-splines. If
  /// an LR spline IS peelable, the LR B-splines are linearly
  /// INdependent. If the LR spline is NOT peelable, the LR B-splines
  /// MAY be linearly dependent, and further investigations are
  /// required to determine whether some of the LR B-splines are
  /// ACTUALLY part of a linear dependence relation.
  //============================================================================
  bool isPeelable( const LRSplineSurface& );

  //============================================================================
  /// Finds unpeelable overloaded LR B-splines of a given LR spline. 
  // An LR B-spline is unpeelable if it cannot be peeled as described in
  // [Dokken, Lyche & Pettersen, 2012]. These unpeeplable LR B-splines
  // corresponds to rows in the incidence matrix with non-zero
  // entries.
  //============================================================================
  std::vector<LRBSpline2D*> unpeelableBasisFunctions ( const LRSplineSurface& );

  } // end namespace LinDepUtils

} // end namespace Go

#endif
