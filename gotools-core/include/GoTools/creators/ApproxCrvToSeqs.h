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

#ifndef _APPROXCRVTOSEQS_H_
#define _APPROXCRVTOSEQS_H_

//   -----------------------------------------------------------------------
//      Interface file for class ApproxCrvToSeqs
//   -----------------------------------------------------------------------
//
//       Approximate sequences of points by a B-spline curves in the same
//       spline space to satisfy a given accuracy
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. ApproxCrvToSeqs.C
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                           10-12
//   -----------------------------------------------------------------------

#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include "GoTools/creators/SmoothCurve.h"

namespace Go
{
    /// This class can generate B-spline curves that approximates
    /// sequences of points for a given accuracy.
class ApproxCrvToSeqs
{
public:
    /// Constructor where the user specifies a set of parameterized points
    /// and a tolerance to be used.  The generated curve will have a spline
    /// basis of order 4 (cubic), and a number of control points equal to 
    /// one-sixth of the number of input points.  The knotvector of the curve's
    /// basis will be set to uniform.
    /// \param points vector containing the sequences points to approximate.  
    ///               These are stored consecutively in "xyzxyzxyz-fashion".
    ///               The dimension of all points is expected to be the same, but
    ///               the number of points in the sequences may differ
    /// \param parvals the parameter values associated with the points to 
    ///                approximate. The parameter values must correspond for all
    ///                sequences
    /// \param dim the spatial dimension of the points (usually 3)
    /// \param aepsge the geometric tolerance to work with
  ApproxCrvToSeqs(const std::vector<std::vector<double> >& points, 
		  const std::vector<std::vector<double> >& parvals,
		  int dim, double aepsge);

    /// Constructor where the user specifies a set of parameterized points
    /// and a tolerance to be used, as well as the spline order and number of
    /// control points.  The knotvector of the curve's basis will be set to 
    /// uniform.
    /// \param points vector containing the sequences points to approximate.  
    ///               These are stored consecutively in "xyzxyzxyz-fashion".
    ///               The dimension of all points is expected to be the same, but
    ///               the number of points in the sequences may differ
    /// \param parvals the parameter values associated with the points to 
    ///                approximate. The parameter values must correspond for all
    ///                sequences
    /// \param dim the spatial dimension of the points (usually 3)
    /// \param aepsge the geometric tolerance to work with
    /// \param in the number of control points of the resulting spline curve
    /// \param ik the order of the resulting spline curve (pol. degree + 1)
    ApproxCrvToSeqs(const std::vector<std::vector<double> >& points, 
		    const std::vector<std::vector<double> >& parvals,
		    int dim, double aepsge, int in, int ik);

    /// Constructor where the user specifies a set of parameterized points,
    /// a tolerance and the spline space to use for constructing the curve.
    /// Constructor where the user specifies a set of parameterized points
    /// and a tolerance to be used, as well as the spline order and number of
    /// control points.  The knotvector of the curve's basis will be set to 
    /// uniform.
    /// \param points vector containing the sequences points to approximate.  
    ///               These are stored consecutively in "xyzxyzxyz-fashion".
    ///               The dimension of all points is expected to be the same, but
    ///               the number of points in the sequences may differ
    /// \param parvals the parameter values associated with the points to 
    ///                approximate. The parameter values must correspond for all
    ///                sequences
    /// \param dim the spatial dimension of the points (usually 3)
    /// \param aepsge the geometric tolerance to work with
    /// \param in the number of control points of the resulting spline curve
    /// \param ik the order of the resulting spline curve (pol. degree + 1)
     /// \param knots specifies the knotvector of the resulting spline curve.
    ApproxCrvToSeqs(const std::vector<std::vector<double> >& points, 
		    const std::vector<std::vector<double> >& parvals, int dim,
		    double aepsge, int in, int ik,
		    std::vector<double>& knots);

    /// Destructor
    ~ApproxCrvToSeqs();

    /// Set smoothing weight 0 < w < 1
    /// \param w the smoothing weight
    void setSmooth(double w);


    /// Unset smoothing weight (set it to zero)
    void unsetSmooth();

    /// The user may decide upon end tangents of approximation curve.
    /// If these are not set, the final end tangents result from smoothing equation.

    /// The user may specifically decide the start and end points and tangents of 
    /// the approximation curve (if these are not set, this information will result
    /// from the smoothing equation).  This function lets the user specify the start
    /// and end points and tangents. 
    /// \param start_point for each curve, this vector should contain up to two 
    ///                    elements: the start point and optionally the start 
    ///                    tangent of the curve.
    /// \param end_point for each curve, this vector should contain up to two 
    ///                  elements: the end point and optionally the end 
    ///                  tangent of the curve.
    void setEndPoints(const std::vector<std::vector<Point> >& start_point, 
		      const std::vector<std::vector<Point> >& end_point);

    /// Approximate C1 continuity with a given impartance 0 <= fac 
    void setC1Approx(double fac);

    /// When everything else is set, this function can be used to fetch the 
    /// approximating curves
    /// \retval maxdist report the maximum distance between the generated curves and 
    ///                 the data points
    /// \retval avdist report the average distance between the generated curves and
    ///                the datapoints
    /// \param max_iter specify the maximum number of iterations to use 
    /// \return a shared pointer to the generated SplineCurve, approximating the points
    ///         as specified.
    std::vector<shared_ptr<SplineCurve> > getApproxCurves(double& maxdist, 
								 double& avdist,
								 int max_iter = 5);
protected:
    /// Default constructor
    ApproxCrvToSeqs();

private:
  std::vector<shared_ptr<SplineCurve> > curr_crv_;
  double maxdist_;
  double avdist_;
  double aepsge_;
  double smoothweight_;
  double smoothfac_;
  double c1fac_;
  double a1_, a2_;

  int dim_;
  std::vector<std::vector<double> > points_;
  std::vector<std::vector<double> > parvals_;
  std::vector<std::vector<Point> > start_pt_; // Pt, der. May be empty.
  std::vector<std::vector<Point> > end_pt_; // Pt, der. May be empty.

  /// Generate an initial curve representing the spline space
  void makeInitCurves();

  void makeInitCurves(int in, int ik, std::vector<double> knots);

  void makeInitCurves(int in, int ik);

  /// Check distribution of knots and parameter values. If
  /// necessary increase the smoothing weight  
  void adjustSmoothWeight();

  /// Generate a smoothing curve
  void makeSmoothCurve(int idx);

  /// Check the accuracy of the current curve
  void checkAccuracy(std::vector<double>& newknots);

  /// Generate the approximating curve
  int doApprox(int max_iter);

  void getConstraints(int idx,
		      std::vector<sideConstraint>& pt_constraints,
		      std::vector<sideConstraint>& tangent_constraints);

};
}

#endif

