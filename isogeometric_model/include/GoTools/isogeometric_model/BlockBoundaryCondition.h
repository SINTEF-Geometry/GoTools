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

#ifndef __BLOCKBOUNDARYCONDITION_H
#define __BLOCKBOUNDARYCONDITION_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BdConditionType.h"
#include "GoTools/topology/tpTopologyTable.h"


namespace Go
{

  class BdCondFunctor;

  // Abstract top class representing the boundary conditions of one block
  // in a block structured model

  class BlockBoundaryCondition
  {
  public:
    // Constructor
    BlockBoundaryCondition(BdConditionType type);

    // Destructor
    virtual
      ~BlockBoundaryCondition();

    // Get the enumeration of affected surface coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) = 0;

    // We also include the inner coefficients next to bd.
    virtual void 
    getCoefficientsEnumeration(std::vector<int>& local_enumeration_bd,
			       std::vector<int>& local_enumeration_bd2) = 0;

    // Get the type
    BdConditionType getBdConditionType() const;

    bool isDirichlet() const;

    // Get the index of the associated solution space
    int getSolutionSpaceIdx() const;

    // Get the coefficients if Dirichlet
    // Represented by the local enumeration and the coefficient itself
    // If the type is not Dirichlet, no spline approximation of the
    // boundary condition exists, and no output is given
    virtual void 
      getBdCoefficients(std::vector<std::pair<int, Point> >& coefs) = 0;

    // We also include the next row of coefficients along the edge,
    // giving a C1 continuity interface.
    virtual void 
    getBdCoefficients(std::vector<std::pair<int, Point> >& coefs_bd,
		      std::vector<std::pair<int, Point> >& coefs_inner) = 0;

    // Update spline approximation if Dirichlet. If not Dirichlet, nothing is done
    virtual void update() = 0;

    // Update the boundary condition with a new functor (for FSI use)
    virtual void updateBoundaryValue(BdCondFunctor* fbd) = 0;

    // Get tolerances
    virtual tpTolerances getTolerances() const = 0;

  protected:
    // Type of boundary conditions
    BdConditionType type_;

  };   // end class BlockBoundaryCondition

} // end namespace Go


#endif    // #ifndef __BLOCKBOUNDARYCONDITION_H
