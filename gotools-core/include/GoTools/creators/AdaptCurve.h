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

#ifndef _ADAPTCURVE_H_
#define _ADAPTCURVE_H_

//   -----------------------------------------------------------------------
//      Interface file for class AdaptCurve
//   -----------------------------------------------------------------------
//
//       Approximate an evaluator based curve by a B-spline curve to
//       satisfy a given accuracy
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. AdaptCurve.C
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                           November 2009
//   -----------------------------------------------------------------------

#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/EvalCurve.h"

namespace Go
{
    /// This class can generate a B-spline curve that approximates
    /// an evaluator based curve within a given accuracy
class AdaptCurve
{
public:
  /// Constructor where the evaluator based curve and the approximation
  /// tolerance are specified.  The generated curve will have a spline
  /// basis of order 4 (cubic).  
  /// \param evalcrv curve to adapt to
  /// \param aepsge  geometric tolerance 
  AdaptCurve(const EvalCurve *evalcrv, double aepsge);

  /// Constructor where the evaluator based curve and the approximation
  /// tolerance are specified as well as the size of the initial spline space  
  /// \param evalcrv curve to adapt to
  /// \param aepsge  geometric tolerance 
  /// \param in the number of control points of the initial spline curve
  /// \param ik the order of the initial spline curve (pol. degree + 1)
  AdaptCurve(const EvalCurve *evalcrv, double aepsge, int in, int ik);

  /// Constructor where the evaluator based curve and the approximation
  /// tolerance are specified as well as the initial spline space  
  /// \param evalcrv curve to adapt to
  /// \param aepsge  geometric tolerance 
  /// \param in the number of control points of the initial spline curve
  /// \param ik the order of the initial spline curve (pol. degree + 1)
  /// \param knots specifies the knotvector of the resulting spline curve.
  AdaptCurve(const EvalCurve *evalcrv, double aepsge, int in, int ik,
	     std::vector<double>& knots);

  /// Constructor where the evaluator based curve and the approximation
  /// tolerance are specified as well as an initial spline curve
  /// \param evalcrv curve to adapt to
  /// \param aepsge  geometric tolerance 
  /// \param crv initial spline curve
  AdaptCurve(const EvalCurve *evalcrv, double aepsge, 
	     shared_ptr<SplineCurve> curve);

  /// Destructor
  ~AdaptCurve();

  /// Set smoothing weight 0 < w < 1
  /// \param w the smoothing weight
  void setSmooth(double w);

  /// Set initial number of sample points per polynomial segment
  void setNmbInitSamples(int nmb_init)
  {
    nmb_pr_seg_ = nmb_init;
  }

  /// Unset smoothing weight (set it to zero)
  void unsetSmooth();
    
  /// The user may decide upon end tangents of approximation curve.
  /// If these are not set, the final end tangents result from smoothing equation.

  /// The user may specifically decide the start and end points and tangents of 
  /// the approximation curve (if these are not set, this information will result
  /// from the smoothing equation).  This function lets the user specify the start
  /// and end points and tangents. 
  /// \param start_point this vector should contain one or two elements: the start
  ///                    point and optionally the start tangent of the curve.
  /// \param end_point this vector should contain one or two elements: the end point
  ///                  and optionally the end tangent of the curve.
  void setEndPoints(const std::vector<Point>& start_point, 
		    const std::vector<Point>& end_point);

  /// Perform the approximation
  /// \param max_iter specify the maximum number of iterations to use 
  int approximate(int max_iter = 5);

  /// Fetch the approximating curve
  /// \retval maxdist report the maximum distance between the generated curve and 
  ///                 the data points generated from the input curve
  /// \retval avdist report the average distance between the generated curve and
  ///                the datapoints
  /// \return a shared pointer to the generated SplineCurve, approximating the 
  ///  input evaluatorbased curve
  shared_ptr<SplineCurve> getAdaptCurve(double& maxdist, 
					       double& avdist,
					       int max_iter = 5);
 protected:
    /// Default constructor
    AdaptCurve();

private:
  shared_ptr<SplineCurve> prev_crv_;
  shared_ptr<SplineCurve> curr_crv_;
  double prev_maxdist_;
  double prev_avdist_;
  double maxdist_;
  double maxdist_inpoints_;
  double avdist_;
  double aepsge_;
  double smoothweight_;
  double smoothfac_;

  const EvalCurve *evalcrv_;

  int dim_;
  std::vector<double> points_;
  std::vector<double> parvals_;
  std::vector<double> pt_weight_;
  std::vector<Point> start_pt_; // Pt, der. May be empty.
  std::vector<Point> end_pt_; // Pt, der. May be empty.

  int cont_;     // Continuity of approximation
  int order_;    // Order of approximating curve

  int fix_[2];   // Number of coefficients to fix in the endpoints
  int init_sample_;  // Initial number of sample points
  int nmb_pr_seg_;   // Number of initial sample point for each polynomial segment
  double min_sample_par_;  // Parameter of first sample point
  double max_sample_par_;  // Parameter of last sample point

  /// Generate an initial curve representing the spline space
  void makeInitCurve();

  void makeInitCurve(int in, int ik, std::vector<double> knots);

  void makeInitCurve(int in, int ik);
  
  /// Initial sample points
  void initSamples();

  /// Check distribution of knots and parameter values. If
  /// necessary increase the smoothing weight  
  void adjustSmoothWeight();

  /// Generate a smoothing curve
  void makeSmoothCurve();

  /// Check the accuracy of the current curve
  void checkAccuracy(std::vector<double>& newknots, int uniform=1);

  void getConstraints(std::vector<sideConstraint>& pt_constraints,
		      std::vector<sideConstraint>& tangent_constraints);

};
}

#endif

