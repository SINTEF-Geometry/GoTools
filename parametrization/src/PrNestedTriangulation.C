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

#include "GoTools/parametrization/PrNestedTriangulation.h"

//----------------------------------------------------------------------------
void PrNestedTriangulation::getParents(int i, int jlev, int& p1, int& p2)
//-----------------------------------------------------------------------------
// Assuming that node i belongs to V^j\V^{j-1} and that
// jlev >= 1, find i's two parent nodes in V^{j-1} p1 and p2.
{
  vector<int> neighbours;
  getNeighbours(i,jlev,neighbours);
  if(isBoundary(i))
  {
    // Then there are exactly four neighbours. We want the
    // first and last ones.
    p1 = neighbours[0];
    p2 = neighbours[3];
    return;
  }
  else
  {
    // Then there are exactly six neighbours. Find the first
    // coarse one.
    if(neighbours[0] < getNumNodes(jlev-1))
    {
      p1 = neighbours[0];
      p2 = neighbours[3];
      return;
    }
    else if(neighbours[1] < getNumNodes(jlev-1))
    {
      p1 = neighbours[1];
      p2 = neighbours[4];
      return;
    }
    else
    {
      p1 = neighbours[2];
      p2 = neighbours[5];
      return;
    }
  }
}

