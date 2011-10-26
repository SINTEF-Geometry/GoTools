/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMEDDHLS_H	
#define PRPRMEDDHLS_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmEDDHLS-syntax: */

/** PrPrmEDDHLS - Implement the EDDHLS parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmEDDHLS : public PrParametrizeInt
{
protected:
  virtual bool makeWeights(int i);

public:
  /// Default constructor 
  PrPrmEDDHLS();
  /// Empty destructor 
  virtual ~PrPrmEDDHLS();

};


/*>PrPrmEDDHLS-syntax: */

/*Class:PrPrmEDDHLS

Name:              PrPrmEDDHLS
Syntax:	           @PrPrmEDDHLS-syntax
Keywords:
Description:       Implement the EDDHLS parametrization
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

#endif // PRPRMEDDHLS_H
