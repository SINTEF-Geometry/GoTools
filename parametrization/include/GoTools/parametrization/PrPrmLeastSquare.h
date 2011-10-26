/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMLEASTSQUARE_H
#define PRPRMLEASTSQUARE_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmLeastSquare-syntax: */

/** PrPrmLeastSquare - Implement the reciprocal least square parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmLeastSquare : public PrParametrizeInt
{
protected:

  virtual
  bool       makeWeights(int i);

public:
  /// Default constructor
  PrPrmLeastSquare();
  /// Empty destructor
  virtual ~PrPrmLeastSquare();

};


/*>PrPrmLeastSquare-syntax: */

/*Class:PrPrmLeastSquare

Name:              PrPrmLeastSquare
Syntax:	           @PrPrmLeastSquare-syntax
Keywords:
Description:       Implement the reciprocal least square parametrization
                   by implementing the virtual function makeWeights.
Member functions:
Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Mar. 97
*/

#endif // PRPRMLEASTSQUARE_H
