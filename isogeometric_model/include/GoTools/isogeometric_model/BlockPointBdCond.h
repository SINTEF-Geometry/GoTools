//===========================================================================
//
// File : BlockPointBdCond.h
//
// Created: November, 2009
//
// Author: Vibeke Skytt, SINTEF, and Anh-Vu Vuong, TU Munich
//
// Revision: 
//
// Description:
//
//===========================================================================


#ifndef __BLOCKPOINTBDCOND_H
#define __BLOCKPOINTBDCOND_H


#include <vector>
#include "GoTools/utils/Point.h"

namespace Go
{

  // Abstract top class for a pointwise boundary condition of Dirichlet type
  class BlockPointBdCond
  {
  public:
    // Constructor
    BlockPointBdCond();

    // Destructor
    virtual ~BlockPointBdCond();

    // Get the index of the associated solution space
    virtual int getSolutionSpaceIdx() const;

    // Get the enumeration of affected surface coefficients (in the solution space)
    // NB! The number of a surface coefficient is not equal to the number of the degrees 
    // of freedom as a coefficient may have several components!
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) const = 0;

    // Get the value of the condition
    virtual Point getConditionValue() const = 0;

    // Get enumeration of affected surface coefficients and the factors in the
    // equation for requiring the current boundary to interpolate the condition
    virtual void 
      getInterpolationFactors(std::vector<std::pair<int,double> >& factors) const = 0;

  protected:

  };   // end class BlockPointBdCond
 
} // end namespace Go


#endif    // #ifndef __BLOCKPOINTBDCOND_H
