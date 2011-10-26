//===========================================================================
//                                                                           
// File: PrPrmExperimental.h                                                 
//                                                                           
// Created: Wed Jul 16 16:46:42 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmExperimental.h,v 1.4 2007-03-02 16:18:19 jbt Exp $
//                                                                           
//===========================================================================

#ifndef _PRPRMEXPERIMENTAL_H
#define _PRPRMEXPERIMENTAL_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/** PrPrmExperimental -  Implement the experimental parametrization
 * by implementing the virtual function makeWeights  
 */
class PrPrmExperimental : public PrParametrizeInt
{
protected:

  virtual bool makeWeights(int i);

public:
  /// Default constructor 
  PrPrmExperimental();
  /// Empty destructor 
  virtual ~PrPrmExperimental();

};



#endif // _PRPRMEXPERIMENTAL_H

