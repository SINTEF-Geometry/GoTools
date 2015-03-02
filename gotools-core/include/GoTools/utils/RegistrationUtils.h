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

#ifndef _REGISTRATIONUTILS_H
#define _REGISTRATIONUTILS_H


#include <vector>
#include "GoTools/utils/Point.h"


namespace Go
{
  /// Enumerator for registration error reasons
  /// RegistrationOK   No error found
  /// TooFewPoints     Less than three points are given
  /// PointSetSizeDiff The two point sets given have different length
  /// AreaTooSmall     Failed to find a tripple of points that are far enough from being almost colinear
  /// SolveFailed      Solving the linear system in the Newtons approach method failed.
  ///                  Error code is stored in registrationData.solve_code
  enum RegistrationReturnType { RegistrationOK, TooFewPoints, PointSetSizeDiff, AreaTooSmall, SolveFailed };

  /// Struct for input to registration, either raw, fine or combined
  struct RegistrationInput
  {
  public:
    RegistrationInput() : 
        area_tolerance_sq_(0.01),
        max_newton_iterations_(10),
        newton_tolerance_(1.0e-8),
        calculate_tolerance_weights_(false),
        tolerance_weight_rotation_(1.0),
        tolerance_weight_translation_(1.0),
        tolerance_weight_rescale_(1.0),
        max_solve_iterations_(150),
        solve_tolerance_(1.0e-8),
        multi_core_(true)
    {
    }


    /// The lower limit of the square of the sine value of the smallest angle in a triangle,
    /// to accept the triangle as good enough for a raw registration
    double area_tolerance_sq_;

    /// Maximum number of iterations in newton approach method
    int max_newton_iterations_;

    /// Upper limit for when weighed sum of squares of changes in Newton's method requires us to break out of the loop
    double newton_tolerance_;

    /// If true, the result of the first iteration of Newtons method will set the weights
    /// for the calculation of the sqaures of changes used to check if we have reached a satisfactory
    /// result and should break out of Newtons method. The weights will be set relatively according
    /// to the changes, with 1.0 as the weight for the translation vector
    bool calculate_tolerance_weights_;

    /// The weight of the square of the change in the rotation vector, used in the Newtons method
    double tolerance_weight_rotation_;

    /// The weight of the square of the change in the translation vector, used in the Newtons method
    double tolerance_weight_translation_;

    /// The weight of the square of the change in the rescaling, used in the Newtons method
    double tolerance_weight_rescale_;

    /// The maximum number of iterations used when solving the linear system in each iteration in Newtons method
    int max_solve_iterations_;

    /// The tolerance used when solving the linear system in each iteration in Newtons method
    double solve_tolerance_;

    /// If true, and if OPENMP is included, run the fine registration in multicore
    bool multi_core_;

    /// Set the tolerance weights used in Newtons method
    void setToleranceWeights(double w_rotation, double w_translation, double w_rescaling)
    {
      tolerance_weight_rotation_ = w_rotation;
      tolerance_weight_translation_ = w_translation;
      tolerance_weight_rescale_ = w_rescaling;
      calculate_tolerance_weights_ = false;
    }
  };

  /// Struct for result from registration process, either raw, fine or combined
  struct RegistrationResult
  {
  public:

    /// The result type of the registration process, either OK or an error code
    RegistrationReturnType result_type_;

    /// The rotation matrix of the registration
    std::vector<std::vector<double> > rotation_matrix_;

    /// The translation of the registration (to be performed after rescaling and rotation)
    Point translation_;

    /// The rescaling factor of the registration;
    double rescaling_;

    /// The iteration number (starting at 0) of the last newton iteration,
    /// or the maximum number of iterations if the weighed sum of squares of changes
    /// never became smaller than the tolerance limit.
    /// Only used for fine registration
    int last_newton_iteration_;

    /// The last weighed sum of squares of changes in Newton's method
    /// Only used for fine registration
    double last_change_;

    /// The return value of the last attempt to solve the linear system in the Newtons method iterations
    /// Only used for fine registration
    int solve_result_;

    /// Return wether the result of the registration was RegistrationOK
    bool ok()
    {
      return result_type_ == RegistrationOK;
    }
  };

  /// Given two sequences of points in 3D, get an approximate rotation, rescaling (optional) and translation that sends the second
  /// point set close to the first. The point sets should be close (but not necessarily identical) in shape, with matching
  /// points coming in the same order (i.e. point number i in the two sequences should be close to each other after performing
  /// the transformation on the second point set), but the original coordinate represenation of the two sets might be totally
  /// different.
  RegistrationResult rawRegistration(const std::vector<Point>& points_fixed, const std::vector<Point>& points_transform,
				     bool allow_rescaling, RegistrationInput params);

  /// Add the contribution from one single point to the linear system coefficients used during fine registration
  void addToLinearSystem(int pt_idx, const std::vector<Point>& points_fixed, const std::vector<Point>& points_transform, bool allow_rescaling,
			 const std::vector<std::vector<double> >& id, const Point& fine_R, const Point& fine_T, double fine_s,
			 const std::vector<std::vector<double> >& m_rot_R, double s2, double R2, bool zero_R,
			 const std::vector<std::vector<std::vector<double> > >& lhs_matrix,
			 const std::vector<std::vector<double> >& rhs_matrix);

  /// Given two sequences of points in 3D, get the rotation, rescaling (optional) and translation that sends the second point set
  /// as close as possible to the first (i.e. that minimizes the sum of the square distances). The sequences must be of
  /// same length (at least three), and the points will be matched in the order they come in the vectors, i.e.
  /// points_transform[i] should, after the transformation, be close to points_fixed[i].
  /// The method works best if the point sets already are close to each other, but not necessarily in optimal position.
  RegistrationResult fineRegistration(const std::vector<Point>& points_fixed, const std::vector<Point>& points_transform,
				      bool allow_rescaling, RegistrationInput params);

  /// Given two sequences of points in 3D, get the rotation, rescaling (optional) and translation that sends the second point set
  /// as close as possible to the first (i.e. that minimizes the sum of the square distances). The sequences must be of
  /// same length (at least three), and the points will be matched in the order they come in the vectors, i.e.
  /// points_transform[i] should, after the transformation, be close to points_fixed[i].
  RegistrationResult registration(const std::vector<Point>& points_fixed, const std::vector<Point>& points_transform,
				  bool allow_rescaling, RegistrationInput params);


} // namespace Go


#endif // _REGISTRATIONUTILS_H

