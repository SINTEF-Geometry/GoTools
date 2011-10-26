/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMUNIFORM_H
#define PRPRMUNIFORM_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmUniform-syntax: */

/** PrPrmUniform -  Implement the uniform parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmUniform : public PrParametrizeInt
{
private:

  virtual
  bool       makeWeights(int i);

public:
  /// Default constructor
  PrPrmUniform();
  /// Empty destructor
  virtual ~PrPrmUniform();

};


/*>PrPrmUniform-syntax: */

/*Class:PrPrmUniform

Name:              PrPrmUniform
Syntax:	           @PrPrmUniform-syntax
Keywords:
Description:       Implement the uniform parametrization
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

#endif // PRPRMUNIFORM_H
