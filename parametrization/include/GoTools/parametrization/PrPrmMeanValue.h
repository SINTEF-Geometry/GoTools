/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2002 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMMEANVALUE_H
#define PRPRMMEANVALUE_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmMeanValue-syntax: */

/** PrPrmMeanValue -  Implement the mean value parametrization
 * by implementing the virtual function makeWeights.
 * See preprint "Mean Value Coordinates", M.S. Floater.
 * This method is now THE recommended method.
 * The mapped points of the mean value parameterization
 * depend infinitely smoothly on the original data points,
 * unlike the shape-preserving parameterization.
 * However, in many examples that have been tested,
 * the two methods give similar results.
 * (MF. August 2002).
 */
class PrPrmMeanValue : public PrParametrizeInt
{
protected:

  virtual bool makeWeights(int i);

public:
  /// Default constructor
  PrPrmMeanValue();
  /// Empty destructor
  virtual ~PrPrmMeanValue();

};


/*>PrPrmMeanValue-syntax: */

/*Class:PrPrmMeanValue

Name:              PrPrmMeanValue
Syntax:	           @PrPrmMeanValue-syntax
Keywords:
Description:       Implement the mean value parametrization
                   by implementing the virtual function makeWeights.
                   See preprint "Mean Value Coordinates", M.S. Floater.
                   This method is now THE recommended method.
                   The mapped points of the mean value parameterization
                   depend infinitely smoothly on the original data points,
                   unlike the shape-preserving parameterization.
                   However, in many examples that have been tested,
                   the two methods give similar results.
                   (MF. August 2002).
Member functions:
                   "attach(PrOrganizedPoints& graph)" --\\
                   Set the planar graph.

Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Apr. 2002
*/

#endif // PRPRMMEANVALUE_H
