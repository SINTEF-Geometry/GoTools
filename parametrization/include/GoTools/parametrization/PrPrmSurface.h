/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPRMSURFACE_H
#define PRPRMSURFACE_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmSurface-syntax: */

/** PrPrmSurface - Implement the shape-preserving parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmSurface : public PrParametrizeInt
{
private:
  std::vector<double> u_;
  std::vector<double> v_;
  std::vector<double> alpha_;
  std::vector<double> len_;

  std::vector<double> uold_;
  std::vector<double> vold_;
  std::vector<double> weightsold_;

  virtual
  bool       makeWeights(int i);
  bool       makeWeights(std::vector<double>& u,
			 std::vector<double>& v,
			 std::vector<double>& weights);
  bool       localParamXYZ(int i);
  bool       localParamUV(int i);

public:
  /// Default constructor
  PrPrmSurface();
  /// Empty destructor
  virtual ~PrPrmSurface();

};


/*>PrPrmSurface-syntax: */

/*Class:PrPrmSurface

Name:              PrPrmSurface
Syntax:	           @PrPrmSurface-syntax
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
Date:              Nov. 98
*/

#endif // PRPRMSURFACE_H
