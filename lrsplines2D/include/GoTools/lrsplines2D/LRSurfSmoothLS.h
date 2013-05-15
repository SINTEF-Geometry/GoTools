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

namespace Go
{

    /// This class modifies a LR spline surface with respect
    /// to conditions on smoothness, data points and boundary conditions.

  //   BsplineIndexMap: Maps each B-spline (represented by a pointer
  //     to a LRBSpline2D) to a linear Index.
  typedef std::map<LRBSpline2D*, size_t> BsplineIndexMap;

class LRSurfSmoothLS
{
 public:
  /// Constructor
  /// Note that the surface is modified
  /// \param coef_known indicates if the coefficient associated to the
  /// LR B-splines is known already
  LRSurfSmoothLS(shared_ptr<LRSplineSurface> surf, std::vector<int>& coef_known);

  /// Destructor
  ~LRSurfSmoothLS();

  /// Reset local arrays after changing LR B-spline surface
  void updateLocals();

  /// Check if data points already are available
  // (are stored in the elements)
  bool hasDataPoints() const;

  /// Add data points for approximation. The points are parameterized
  /// and are store as follows: parameter values for point 1 (2 doubles),
  /// position for point1 (dim doubles where dim is the dimension of the
  /// associated LR spline surface), parameter values for point 2 etc.
  void addDataPoints(std::vector<double>& points);

  /// Compute the smoothing part of the equation system.
  /// \param weight1 contribution weight with respect to the 1st derivative.
  /// \param weight2 contribution weight with respect to the 2nd derivative.
  /// \param weight3 contribution weight with respect to the 3rd derivative.
  void setOptimize(const double weight1, const double weight2,
		   const double weight3);

  /// Compute matrices for least squares approximation.
  /// The data points are expected to be added already
  /// \param weight the contribution of the approximation of the pnts in the system.
  ///               weight should lie in the unit interval.
  void setLeastSquares(const double weight);

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
			 const std::vector<LRBSpline2D*>& bsplines,
			 double* mat, double* right, int ncond);

  std::vector<double> getBasisValues(const std::vector<LRBSpline2D*>& bsplines,
				     double *par);

}; // end of class LRSurfSmoothLS


} // end of namespace Go

#endif
