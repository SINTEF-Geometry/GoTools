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
			   std::vector<std::pair<double, double> >& domain, VolSolution *solution);

      VolBoundaryCondition(int face_nmb, BdConditionType type,
			   const Point& const_val,
			   std::vector<std::pair<double, double> >& domain,
			   VolSolution *solution);

    // Destructor
    virtual ~VolBoundaryCondition();

    // Get the enumeration of affected volume coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration);

    // We also include the inner coefficients next to bd row.
    virtual void 
    getCoefficientsEnumeration(std::vector<int>& local_enumeration_bd,
			       std::vector<int>& local_enumeration_bd2);

    // Get the coefficients if Dirichlet
    // Represented by the local enumeration and the coefficient itself
    // If the type is not Dirichlet, no spline approximation of the
    // boundary condition exists, and no output is given
    virtual void
      getBdCoefficients(std::vector<std::pair<int, Point> >& coefs);

    // We also include the next row of coefficients along the edge,
    // giving a C1 continuity interface.
    virtual void
    getBdCoefficients(std::vector<std::pair<int, Point> >& coefs_bd,
		      std::vector<std::pair<int, Point> >& coefs_inner);

    // Evaluation of the spline space related to the current boundary condition
    void getBasisFunctions(int index_of_Gauss_point1, int index_of_Gauss_point2,
			   int& const_dir,
			   vector<double>& basisValues,
			   vector<double>& basisDerivs_u,
			   vector<double>& basisDerivs_v) const;
			   // shared_ptr<BasisDerivs> result,
			   // int solutionspace_idx) const;

    // Get the spline approximation (if it exists)
    // The spline approximation is typically only created for Dirichlet boundary conditions
    shared_ptr<SplineSurface> getSplineApproximation() const;

    // Update spline approximation if Dirichlet. If not Dirichlet, nothing is done
    virtual void update();

    // Update the boundary condition with a new functor (for FSI use)
    virtual void updateBoundaryValue(BdCondFunctor* fbd);

    // Get edge number
    int faceNumber() const;

    // Get tolerances
    virtual tpTolerances getTolerances() const;

  private:


    // Pointer to the block solution to which this boundary condition belongs
    VolSolution* parent_;

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
  
    // Functor that is able to evaluate the boundary condition
    BdCondFunctor *fbd_;

    // Approximating spline surface
    shared_ptr<SplineSurface> bdsrf_cond_;

    // Constant value for constant Dirichlet boundary conditions
    Point const_val_;

    // Parameter domain associated with this condition
    // The domain is described by a polygon in the parameter domain
    // @@sbr201209 Is the polygon closed, or could it describe a curve?
    std::vector<std::pair<double, double> > domain_;

    // Approximation error at last spline curve approximation
    // A value of -1.0 means no approximation has occured yet
    double approx_err_;

  };   // end class VolBoundaryCondition

} // end namespace Go


#endif    // #ifndef __VOLBOUNDARYCONDITION_H
