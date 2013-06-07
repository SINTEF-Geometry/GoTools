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

#ifndef __BLOCKPOINTBDCOND_H
#define __BLOCKPOINTBDCOND_H


#include <vector>
#include "GoTools/utils/Point.h"

namespace Go
{

  // Abstract top class for a pointwise boundary condition of Dirichlet type
  class BlockPointBdCond
  {
  public:
    // Constructor
    BlockPointBdCond();

    // Destructor
    virtual ~BlockPointBdCond();

    // Get the index of the associated solution space
    virtual int getSolutionSpaceIdx() const;

    // Get the enumeration of affected surface coefficients (in the solution space)
    // NB! The number of a surface coefficient is not equal to the number of the degrees 
    // of freedom as a coefficient may have several components!
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) const = 0;

    // Get the value of the condition
    virtual Point getConditionValue() const = 0;

    // Get enumeration of affected surface coefficients and the factors in the
    // equation for requiring the current boundary to interpolate the condition
    virtual void 
      getInterpolationFactors(std::vector<std::pair<int,double> >& factors) const = 0;

  protected:

  };   // end class BlockPointBdCond
 
} // end namespace Go


#endif    // #ifndef __BLOCKPOINTBDCOND_H
