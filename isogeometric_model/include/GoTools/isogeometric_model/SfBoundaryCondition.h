//===========================================================================
//
// File : SfBoundaryCondition.h
//
// Created: November, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh Vu-Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================



#ifndef __SFBOUNDARYCONDITION_H
#define __SFBOUNDARYCONDITION_H


#include <vector>
#include <memory>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BlockBoundaryCondition.h"



namespace Go
{

  class SfSolution;

  // The SfBoundaryCondition class collects all information related to the
  // boundary conditions of one solution space corresponding to the current
  // surface block.
  // One edge of the surface block may correspond to more than one boundary
  // condition, but only if they are of different type (Dirichlet and Neumann,
  // not zero Dirichlet and Dirichlet equal to some other constant).
  // Non-matching boundary conditions of the same type is represented by the 
  // same functor, and in the Dirichlet case approximated as good as possible
  // in the current spline space.

  class SfBoundaryCondition : public BlockBoundaryCondition
  {
  public:
    // Constructors
    SfBoundaryCondition(int edge_nmb, BdConditionType type, BdCondFunctor *fbd,
			std::pair<double, double> end_par, SfSolution *solution);

    SfBoundaryCondition(int edge_nmb, BdConditionType type, const Point& const_val,
			std::pair<double, double> end_par, SfSolution *solution);

    // Destructor
    virtual ~SfBoundaryCondition();

    // Get the enumeration of affected surface coefficients. The coefficients come in
    // increasing or decreasing order, depending on whether the domain parameters are
    // in increasing or decreasing order.
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration);

    // We also include coefficient row number next to bd row.
    virtual void 
    getCoefficientsEnumeration(std::vector<int>& local_enumeration_bd,
			       std::vector<int>& local_enumeration_bd2);

    // Get the coefficients if Dirichlet
    // Represented by the local enumeration and the coefficient itself
    // If the type is not Dirichlet, no spline approximation of the
    // boundary condition exists, and no output is given
    /// \param coefs the boundary coefficients. Each pair consists
    ///        the coefficient and its index when considered a curve.
    virtual void 
      getBdCoefficients(std::vector<std::pair<int, Point> >& coefs);

    // We also include the next row of coefficients along the edge,
    // giving a C1 continuity interface.
    /// \param coefs_bd the boundary coefficients. Each pair consists
    ///        the coefficient and its index when considered a curve.
    /// \param coefs_bd2 the row of coefficients next to the boundary.
    ///        Each pair consists the coefficient and its index when
    ///        considered a curve.
    virtual void 
    getBdCoefficients(std::vector<std::pair<int, Point> >& coefs_bd,
		      std::vector<std::pair<int, Point> >& coefs_bd2);

    // Update spline approximation if Dirichlet. If not Dirichlet, nothing is done
    virtual void update();

    // Evalution of the spline space related to the current boundary condition
    // The evaluation takes place along the boundary curve, and the point and derivatives
    // are returned, as well as a boolean indicating whether the boundary curve parameter
    // direction is the u or v direction of the surface
    void getBasisFunctions(int index_of_Gauss_point,
			   bool& u_dir,
			   std::vector<double>& basisValues,
			   std::vector<double>& basisDerivs) const;

    // Get the spline approximation (if it exists)
    // The spline approximation is typically only created for Dirichlet boundary conditions
    shared_ptr<SplineCurve> getSplineApproximation() const;

    // Update the boundary condition with a new functor (for FSI use)
    virtual void updateBoundaryValue(BdCondFunctor* fbd);

    // Get edge number
    int edgeNumber() const;

    // Get tolerances
    virtual tpTolerances getTolerances() const;

  private:

    // Pointer to the block solution to which this boundary condition belongs
    SfSolution* parent_;

    // The edge it corresponds to
    // edgenmb_ = 0: the boundary corresponding to the minimum parameter in the first parameter
    //                 direction (u_min)
    // edgenmb_ = 1: the boundary corresponding to the maximum parameter in the first parameter
    //                 direction (u_max)
    // edgenmb_ = 2: the boundary corresponding to the minimum parameter in the second parameter
    //                 direction (v_min)
    // edgenmb_ = 3: the boundary corresponding to the maximum parameter in the second parameter
    //                 direction (v_max)
    int edgenmb_; 

    // Functor that are able to evaluate the boundary condition
    BdCondFunctor *fbd_;

    // Approximating spline curve
    shared_ptr<SplineCurve> bdcrv_cond_;

    // Constant value for constant Dirichlet boundary conditions
    Point const_val_;

    // Parameter domain associated with this condition. might be in decreasing order if the direction of the
    // condition curve is opposite of the parameter direction
    std::pair<double, double> domain_;

    // Approximation error at last spline curve approximation
    // A value of -1.0 means no approximation has occured yet
    double approx_err_;

  };   // class SfBoundaryCondition

} // end namespace Go


#endif    // #ifndef __SFBOUNDARYCONDITION_H
