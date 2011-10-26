/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRTHRESHOLD_H
#define PRTHRESHOLD_H

#include "GoTools/parametrization/PrNestedTriangulation.h"
#include "GoTools/utils/Array.h"
using Go::Vector3D;
#include <memory>

/*<PrThreshold-syntax: */

/** PrThreshold -  This class implements an algorithm for thresholding
 * a linear combination of coarse nodal functions and wavelets over several levels.
 */
class PrThreshold
{
protected:

  std::shared_ptr<PrNestedTriangulation> t_;
  double          threshold_; // for wavelet thresholding

public:
  /// Default constructor
  PrThreshold();
  /// Empty destructor
  virtual ~PrThreshold() {}

  virtual
  void        attach(std::shared_ptr<PrNestedTriangulation> t);
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
