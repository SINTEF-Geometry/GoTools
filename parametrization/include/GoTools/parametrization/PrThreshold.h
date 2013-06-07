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

#ifndef PRTHRESHOLD_H
#define PRTHRESHOLD_H

#include "GoTools/parametrization/PrNestedTriangulation.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/config.h"
using Go::Vector3D;
#include <memory>

/*<PrThreshold-syntax: */

/** PrThreshold -  This class implements an algorithm for thresholding
 * a linear combination of coarse nodal functions and wavelets over several levels.
 */
class PrThreshold
{
protected:

  shared_ptr<PrNestedTriangulation> t_;
  double          threshold_; // for wavelet thresholding

public:
  /// Default constructor
  PrThreshold();
  /// Empty destructor
  virtual ~PrThreshold() {}

  virtual
  void        attach(shared_ptr<PrNestedTriangulation> t);
  void        setThreshold(double threshold = 0.0)
                    {threshold_ = threshold;}

  /** Threshold the wavelet coefficients by setting those
   * whose absolute value is less than the tolerance to 0.
   * Return the compression rate: ratio of no. of coarse hats plus
   * non-zero wavelets over no. of coarse hats plus all wavelets.
   * Return the wavelet compression rate: ratio of no. of non-zero
   * over no. of wavelets.
   *
   * If error_out = 1, do the opposite of thresholding,
   * i.e, retain only those wavelet coefficients which are smaller
   * in absolute value than the tolerance.
   * The resulting piecewise linear function represents the error after
   * thresholding and can be used to find the l2 error, max error, etc.
*/
  void       threshold(double& comp_rate, double& wavelet_comp_rate,
                       int error_out = 0);

  /** Find a threshold which will give the given compressiosn rate.
   * E.g. if wavelet_comp_rate = 0.05 and there are n wavelets then
   * we want a threshold so that m = int(0.05 * n) are retained.
   */
  double     findThreshold(double& wavelet_comp_rate);

  /// Threshold by given compression rate, e.g. 0.05.
  void       thresholdByCompRate(double& wavelet_comp_rate,
                       int error_out = 0);

  /// Set all coeffs at level jlev to zero (0 <= jlev <= getFinestLevel).
  void       truncateLevel(int jlev);

  /// Set all level zero coeffs to zero and all wavelet coeffs except those of level j.
  void       leaveLevel(int jlev);

  /// Calculate the max norm.
  double     maxNorm();

  /// Calculate the average absolute value.
  double     averageAbsValue();
  double     weightedL2Norm();
  double     triangleNorm(int i, int j, int k);
/* The nest two routines need 3D parameter points:
  double       L2Norm();
  double       triangleAreaNorm(int i, int j, int k);
*/
};

/*>PrThreshold-syntax: */

/*Class:PrThreshold

Name:              PrThreshold
Syntax:	           @PrThreshold-syntax
Keywords:
Description:       This class implements an algorithm for thresholding
                   a linear combination of coarse nodal functions and
                   wavelets over several levels.

Member functions:

Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Dec. 98
*/

#endif // PRTHRESHOLD_H
