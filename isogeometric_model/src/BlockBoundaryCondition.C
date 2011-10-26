//===========================================================================
//
// File : BlockBoundaryCondition.C
//
// Created: Mon Mar  1 10:57:48 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================





#include "GoTools/isogeometric_model/BlockBoundaryCondition.h"


namespace Go
{


  //===========================================================================
  BlockBoundaryCondition::BlockBoundaryCondition(BdConditionType type):
    type_(type)
  //===========================================================================
  {
  }


  //===========================================================================
  BlockBoundaryCondition::~BlockBoundaryCondition()
  //===========================================================================
  {
  }


  //===========================================================================
  BdConditionType BlockBoundaryCondition::getBdConditionType() const
  //===========================================================================
  {
    return type_;
  }

  //===========================================================================
  bool BlockBoundaryCondition::isDirichlet() const
  //===========================================================================
  {
    return type_ == ZERO_DIRICHLET || type_ == CONSTANT_DIRICHLET || type_ == DIRICHLET;
  }
} // end namespace Go
