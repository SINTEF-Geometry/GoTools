//===========================================================================
//
// File : BlockBoundaryCondition.h
//
// Created: November, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================


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
