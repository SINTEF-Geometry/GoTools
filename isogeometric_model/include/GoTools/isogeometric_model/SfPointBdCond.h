//===========================================================================
//
// File : SfPointBdCond.h
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


#ifndef __SFPOINTBDCOND_H
#define __SFPOINTBDCOND_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BlockPointBdCond.h"



namespace Go
{

  class SfSolution;

  // This class represents a pointwise boundary condition of Dirichlet type
  class SfPointBdCond : public BlockPointBdCond
  {
  public:
    // Constructor
    SfPointBdCond(int edge_nmb, double param, Point& condition_value);

    // Destructor
    virtual ~SfPointBdCond();

    // Get the enumeration of affected surface coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) const;

    // Get the value of the condition
    virtual Point getConditionValue() const;

    // Get the parameter corresponding to the condition value
    double getParam() const;

    // Get enumeration of affected surface coefficients and the factors in the
    // equation for requiring the current boundary to interpolate the condition
    virtual void 
      getInterpolationFactors(std::vector<std::pair<int,double> >& factors) const;

    // Get edge number
    int edgeNumber() const;

  private:
    // The edge it corresponds to
    // edgenmb_ = 0: the boundary corresponding to the minimum parameter in the first parameter
    //                 direction (u_min)
    // edgenmb_ = 1: the boundary corresponding to the maximum parameter in the first parameter
    //                 direction (u_max)
    // edgenmb_ = 2: the boundary corresponding to the minimum parameter in the second parameter
    //                 direction (v_min)
    // edgenmb_ = 3: the boundary corresponding to the maximum parameter in the second parameter
    //                 direction (v_max)
    int edgenmb_; 

    // The value of the boundary condition of the solution in the given point
    Point value_;

    // The parameter value corresponding to this pointwise condition
    double param_;

    // Pointer to the block solution to which this boundary condition belongs
    SfSolution* parent_;

  };  // end class SfPointBdCond

} // end namespace Go


#endif    // #ifndef __SFPOINTBDCOND_H
