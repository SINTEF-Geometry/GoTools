//===========================================================================
//
// File : BdCondFunctor.h
//
// Created: Fri Feb 26 11:07:45 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#ifndef __BDCONDFUNCTOR_H
#define __BDCONDFUNCTOR_H


#include "GoTools/utils/Point.h"


namespace Go
{

  class BdCondFunctor
  {

  public:

    // Constructor
    BdCondFunctor();

    // Destructor
    virtual ~BdCondFunctor();

    /// Functor evaluation in point
    virtual Point evaluate(const Point& geom_pos) = 0;

  private:

  };    // Class BdCondFunctor


} // namespace Go

#endif    // #ifndef __BDCONDFUNCTOR_H
