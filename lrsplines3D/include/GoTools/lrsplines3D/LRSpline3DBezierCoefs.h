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

#ifndef _LRSPLINE3DBEZIERCOEFS_H
#define _LRSPLINE3DBEZIERCOEFS_H

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/utils/BoundingBox.h"


#include <vector>
#include <string>
#include <fstream>

namespace Go 
{

  /// Prepare for visualization of LR B-spline volume by performing Bezier extraction.
  /// Coefficients of Bezier volumes are computed by interpolation of a set of
  /// sample points depending on the degrees of the LR B-spline volume.
// =============================================================================
class LRSpline3DBezierCoefs 
// =============================================================================
{
public:

  /// Empty constructor
  LRSpline3DBezierCoefs();

  /// Constructor given LRSplineVolume
  LRSpline3DBezierCoefs(LRSplineVolume& lr_spline);

  /// Fetch coefficients of Bezier volumes after extraction.
  void getBezierCoefs();

  /// Write Bezier coefficients and associated information to stream
  void writeToStream(std::ostream& os);

  /// Write Bezier coefficients and associated information to file
  void writeToFile(const std::string& filename);

private:
  
  void computeCoefsFromPts(const double *points, double *coefs);

  void calcMinMaxCoefficientsValue();

  LRSplineVolume lr_spline_; // TODO: remove the need for this
 
  int order_u_;
  int order_v_;
  int order_w_;
  int dim_;
  int num_elements_;
  
  Array<double,6> orig_dom_;
  std::vector<double> bezier_coefs_; 
  std::vector<double> min_coef_value_;
  std::vector<double> max_coef_value_;
  std::vector<Go::BoundingBox> boxes_;

  double min_box_diagonal_;
  double max_box_diagonal_;

  static constexpr double M3_[9] = // interpolation matrix for quadratic splines
  {
    1,   0,   0,
    -.5,   2, -.5,
    0,   0,   1
  };
  static constexpr double M4_[16] = // interpolation matrix for cubic splines
  {
         1,        0,      0,      0,
    -5/6.0,   18/6.0, -9/6.0,  2/6.0,
     2/6.0,   -9/6.0, 18/6.0, -5/6.0,
         0,        0,      0,      1
  };

};

} // end namespace Go

#endif // _LRSPLINE3DBEZIERCOEFS_H

