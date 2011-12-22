//===========================================================================
//                                                                           
// File: VolBoundaryCondition.C                                              
//                                                                           
// Created: Tue Nov  1 10:17:18 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#if 0

#include "GoTools/isogeometric_model/VolBoundaryCondition.h"


namespace Go
{

  //===========================================================================
  VolBoundaryCondition::VolBoundaryCondition(int face_nmb, BdConditionType type,
					     BdCondFunctor *fbd,
					     std::vector<std::pair<double, double> >& domain)
    : BlockBoundaryCondition(type),
      fbd_(fbd),
      domain_(domain)
  //===========================================================================
  {
    MESSAGE("VolBoundaryCondition() under construction!");
  }

  //===========================================================================
  VolBoundaryCondition::~VolBoundaryCondition()
  //===========================================================================
  {
    MESSAGE("~() not implemented");
  }

  //===========================================================================
  void 
  VolBoundaryCondition::getCoefficientsEnumeration(std::vector<int>& local_enumeration)
  //===========================================================================
  {
    MESSAGE("getCoefficientsEnumeration() not implemented");
  }

  //===========================================================================
  void
  VolBoundaryCondition::getBdCoefficients(std::vector<std::pair<int, Point> >& coefs)
  //===========================================================================
  {
    MESSAGE("getBdCoefficients() not implemented");
  }

  //===========================================================================
  void VolBoundaryCondition::getBasisFunctions(int index_of_Gauss_point1,
					       int index_of_Gauss_point2,
					       shared_ptr<BasisDerivs> result,
					       int solutionspace_idx) const
  //===========================================================================
  {
    MESSAGE("getBasisFunctions() not implemented");
  }

  //===========================================================================
  shared_ptr<SplineSurface> VolBoundaryCondition::getSplineApproximation() const
  //===========================================================================
  {
    MESSAGE("getSplineApproximation() not implemented");
    shared_ptr<SplineSurface> sf;
    return sf;
  }

  //===========================================================================
  void VolBoundaryCondition::update()
  //===========================================================================
  {
    MESSAGE("update() not implemented");
  }

  //===========================================================================
  void VolBoundaryCondition::updateBoundaryValue(BdCondFunctor* fbd)
  //===========================================================================
  {
    MESSAGE("updateBoundaryValue() not implemented");
  }

  //===========================================================================
  tpTolerances VolBoundaryCondition::getTolerances() const
  //===========================================================================
  {
    return parent_->getTolerances();
  }

}


#endif
