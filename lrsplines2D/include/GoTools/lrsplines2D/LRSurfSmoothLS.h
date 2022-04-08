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

#ifndef _LRSURFSMOOTHLS_H
#define _LRSURFSMOOTHLS_H

//===========================================================================
//
// File : LRSurfSmoothLS
//
// Created: May 2013
//
// Author: Vibeke Skytt
//
// Revision: $Id:$
//
// Description: Least squares approximation with smoothing for LR spline surface
//              Note that rational surfaces are currently not handled
//
//===========================================================================

#include <vector>
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"

namespace Go
{


  //   BsplineIndexMap: Maps each B-spline (represented by a pointer
  //     to a LRBSpline2D) to a linear Index.
  typedef std::map<LRBSpline2D*, size_t> BsplineIndexMap;

    /// This class modifies a LR spline surface with respect
    /// to conditions on smoothness, data points and boundary conditions.
  /// Method: Least squares approximation with a smoothing term.
class LRSurfSmoothLS
{
 public:
  /// Constructor.
  /// Note that the surface is modified
  /// \param coef_known indicates if the coefficient associated to the
  /// LR B-splines is known already
  LRSurfSmoothLS(shared_ptr<LRSplineSurface> surf, std::vector<int>& coef_known);
  /// Empty constructor
  LRSurfSmoothLS();

  /// Destructor
  ~LRSurfSmoothLS();

  /// Set initial surface
  /// Note that the surface is modified
  /// \param coef_known indicates if the coefficient associated to the
  /// LR B-splines is known already
  void setInitSf(shared_ptr<LRSplineSurface> surf, 
		 std::vector<int>& coef_known);

  /// Reset local arrays after changing LR B-spline surface
  void updateLocals();

  /// Check if data points already are available
  // (are stored in the elements)
  bool hasDataPoints() const;

  /// Add data points for approximation. The points are parameterized
  /// and stored as follows: parameter values for point 1 (2 doubles),
  /// position for point1 (dim doubles where dim is the dimension of the
  /// associated LR spline surface), parameter values for point 2 etc.
  void addDataPoints(std::vector<double>& points, 
		     LRSplineUtils::PointType type = LRSplineUtils::REGULAR_POINTS);

  /// Compute the smoothing part of the equation system.
  /// \param weight1 contribution weight with respect to the 1st derivative.
  /// \param weight2 contribution weight with respect to the 2nd derivative.
  /// \param weight3 contribution weight with respect to the 3rd derivative.
  void setOptimize(const double weight1, const double weight2,
		   const double weight3);

  /// Compute the boundary smoothing part of the equation system.
  /// \param weight1 contribution weight with respect to the 1st derivative.
  /// \param weight2 contribution weight with respect to the 2nd derivative.
  /// \param weight3 contribution weight with respect to the 3rd derivative.
  void smoothBoundary(const double weight1, const double weight2,
		      const double weight3);

  /// Compute matrices for least squares approximation.
  /// The data points are expected to be added already
  /// \param weight the contribution of the approximation of the pnts in the system.
  ///               weight should lie in the unit interval.
  void setLeastSquares(const double weight, const double significant_factor);

  /// OpenMP enabled version of the above function.
  void setLeastSquares_omp(const double weight, 
			   const double significant_factor);

  /// Compute matrices for least squares approximation.
  /// \param points Parameter values and point for each data point
  /// \param weight the contribution of the approximation of the pnts in the system.
  ///               weight should lie in the unit interval.
  void setLeastSquares(std::vector<double>& points, const double weight);

  /// Solve equation system, and produce output surface.
  /// If failing to solve the routine may throw an exception.
  /// \param surf the output surface.
  /// \return 0 = OK, negative = failed solving system.
  int equationSolve(shared_ptr<LRSplineSurface>& surf);

 private:
  shared_ptr<LRSplineSurface> srf_;  // Pointer to input surface.
  std::vector<int> coef_known_;
  int ncond_;                        // Number of unknown coefficients

  /// Storage of the equation system.
  std::vector<double> gmat_;         // Matrix at left side of equation system.  
  std::vector<double> gright_;       // Right side of equation system.      
 
  BsplineIndexMap BSmap_;   // Indices to all LR B-splines to associate
                            // a posistion in the stiffness matrix

  // Compute the least squares contributions to the stiffness matrix and
  // the right hand side for a specified set of B-splines
  void localLeastSquares(std::vector<double>& points, 
			 std::vector<double>& significant_points, 
			 std::vector<double>& ghost_points, int del,
			 const double significant_factor,
			 const std::vector<LRBSpline2D*>& bsplines,
			 double* mat, double* right, int ncond);

  void localLeastSquares_omp(std::vector<double>& points, 
			     std::vector<double>& significant_points, 
			     std::vector<double>& ghost_points, int del,
			     const double significant_factor,
			     const std::vector<LRBSpline2D*>& bsplines,
			     double* mat, double* right, int ncond);

  std::vector<double> getBasisValues(const std::vector<LRBSpline2D*>& bsplines,
				     double *par);

  void fetchBasisDerivs(const std::vector<LRBSpline2D*>& bsplines, 
			std::vector<double>& basis_derivs, 
			int der1, int der2, int der3, 
			double umin, double umax,
			double vmin, double vmax, int& nmbGauss);

  void evalAllBGridDer(const std::vector<LRBSpline2D*>& bsplines,
		       int nmb_der,
		       const std::vector<double>& par1, 
		       const std::vector<double>& par2, 
		       std::vector<double>& result);

  void fetchBasisLineDerivs(const std::vector<LRBSpline2D*>& bsplines, 
			    std::vector<double>& basis_derivs, 
			    int der1, int der2, int der3, 
			    Direction2D d, double tmin, 
			    double tmax, int& nmbGauss);

  void computeDer1Integrals(const std::vector<LRBSpline2D*>& bsplines, 
			    int nmbGauss, double* basis_derivs, double weight);
  void computeDer1LineIntegrals(const std::vector<LRBSpline2D*>& bsplines, 
				int nmbGauss, double* basis_derivs, double weight);

  void computeDer2Integrals(const std::vector<LRBSpline2D*>& bsplines, 
			    int nmbGauss, double* basis_derivs, double weight);
  void computeDer2LineIntegrals(const std::vector<LRBSpline2D*>& bsplines, 
				int nmbGauss, double* basis_derivs, double weight);

  void computeDer3Integrals(const std::vector<LRBSpline2D*>& bsplines, 
			    int nmbGauss, double* basis_derivs, double weight);
  void computeDer3LineIntegrals(const std::vector<LRBSpline2D*>& bsplines, 
				int nmbGauss, double* basis_derivs, double weight);

  std::vector<LRBSpline2D*> 
    bsplinesCoveringElement(std::vector<LRBSpline2D*>& cand, 
			    Direction2D d, double tmin, double tmax);  
}; // end of class LRSurfSmoothLS


} // end of namespace Go

#endif
