//===========================================================================
//
// File : VolBoundaryCondition.h
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



#ifndef __VOLBOUNDARYCONDITION_H
#define __VOLBOUNDARYCONDITION_H


#include <vector>
#include <memory>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/isogeometric_model/BlockBoundaryCondition.h"



namespace Go
{

  class VolSolution;

  // The VolBoundaryCondition class collects all information related to the
  // boundary conditions of one solution space corrsponding to the current
  // volume block.
  // One boundary surface of the volume block may correspond to more than one boundary
  // condition, but only if they are of different type (Dirichlet and Neumann,
  // not zero Dirichlet and Dirichlet equal to some other constant).
  // Non-matching boundary conditions of the same type is represented by the 
  // same functor, and in the Dirichlet case approximated as good as possible
  // in the current spline space.
  class VolBoundaryCondition : public BlockBoundaryCondition
  {
  public:
    // Constructor
    VolBoundaryCondition(int face_nmb, BdConditionType type, BdCondFunctor *fbd,
			 std::vector<std::pair<double, double> >& domain);

    // Destructor
    virtual ~VolBoundaryCondition();

    // Get the enumeration of affected volume coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration);

    // Get the coefficients if Dirichlet
    // Represented by the local enumeration and the coefficient itself
    // If the type is not Dirichlet, no spline approximation of the
    // boundary condition exists, and no output is given
    virtual void
      getBdCoefficients(std::vector<std::pair<int, Point> >& coefs);

    // Evaluation of the spline space related to the current boundary condition
    void getBasisFunctions(int index_of_Gauss_point1, int index_of_Gauss_point2,
			   std::shared_ptr<BasisDerivs> result,
			   int solutionspace_idx) const;

    // Get the spline approximation (if it exists)
    // The spline approximation is typically only created for Dirichlet boundary conditions
    std::shared_ptr<SplineSurface> getSplineApproximation() const;

    // Update spline approximation if Dirichlet. If not Dirichlet, nothing is done
    virtual void update();

    // Update the boundary condition with a new functor (for FSI use)
    virtual void updateBoundaryValue(BdCondFunctor* fbd);

    // Get tolerances
    virtual tpTolerances getTolerances() const;

  private:
    // The boundary surface it corresponds to
    // facenmb_ = 0: the boundary corresponding to the minimum parameter in the first parameter
    //                 direction (u_min)
    // facenmb_ = 1: the boundary corresponding to the maximum parameter in the first parameter
    //                 direction (u_max)
    // facenmb_ = 2: the boundary corresponding to the minimum parameter in the second parameter
    //                 direction (v_min)
    // facenmb_ = 3: the boundary corresponding to the maximum parameter in the second parameter
    //                 direction (v_max)
    // facenmb_ = 4: the boundary corresponding to the minimum parameter in the third parameter
    //                 direction (w_min)
    // facenmb_ = 5: the boundary corresponding to the maximum parameter in the third parameter
    //                 direction (w_max)
    int facenmb_;   
  
    // Functor that are able to evaluate the boundary condition
    BdCondFunctor *fbd_;

    // Approximating spline surface
    std::shared_ptr<SplineSurface> bdsrf_cond_;

    // Parameter domain associated with this condition
    // The domain is described by a polygon in the parameter domain
    std::vector<std::pair<double, double> > domain_;

    // Pointer to the block solution to which this boundary condition belongs
    VolSolution* parent_;

  };   // end class VolBoundaryCondition

} // end namespace Go


#endif    // #ifndef __VOLBOUNDARYCONDITION_H
