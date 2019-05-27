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

#ifndef _MESHLR_H
#define _MESHLR_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "GoTools/geometry/Streamable.h"
#include "GoTools/lrsplines2D/IndexMesh2DIterator.h"

namespace Go
{

// =============================================================================
class MeshLR : public Streamable 
// =============================================================================
{
public:

  // Get a pointer to the start of the knot vector in the given direction.
  virtual const double* const knotsBegin(int pardir) const = 0;
  
  // Get a pointer to the one-past-end of the knot vector in the given direction.
  virtual const double* const knotsEnd  (int pardir) const = 0;

  // Get the number of distinct knot valuess in a given direction (rows: 2, columns: 1).
  // Note that this is the number of _distinct_ knots, so multiplicities are not taken into
  // account.
  virtual int numDistinctKnots(int pardir) const = 0;

  // Return the knot value for the knot with index 'ix' along direction 'd'
  // (rows: 2, columns: 1).
  virtual double kval(int pardir, int ix) const = 0;


   
}; // end class MeshLR


}; // end namespace Go

#endif 
