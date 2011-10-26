/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMSHPPRES_H
#define PRPRMSHPPRES_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmShpPres-syntax: */

/** PrPrmShpPres - Implement the shape-preserving parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmShpPres : public PrParametrizeInt
{
protected:
  std::vector<double> u_;
  std::vector<double> v_;
  std::vector<double> alpha_;
  std::vector<double> len_;

  virtual bool makeWeights(int i);
  bool         localParam(int i);

public:
  /// Default constructor 
  PrPrmShpPres();
  /// Empty destructor
  virtual ~PrPrmShpPres();

};


/*>PrPrmShpPres-syntax: */

/*Class:PrPrmShpPres

Name:              PrPrmShpPres
Syntax:	           @PrPrmShpPres-syntax
Keywords:
Description:       Implement the shape-preserving parametrization
                   by implementing the virtual function makeWeights.
Member functions:
                   "attach(PrOrganizedPoints& graph)" --\\
                   Set the planar graph.

Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Mar. 97
*/

#endif // PRPRMSHPPRES_H
