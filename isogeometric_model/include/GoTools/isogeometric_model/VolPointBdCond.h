//===========================================================================
//
// File : VolPointBdCond.h
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


#ifndef __VOLPOINTBDCOND_H
#define __VOLPOINTBDCOND_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BlockPointBdCond.h"



namespace Go
{

  class VolSolution;

  // This class represents a pointwise boundary condition of Dirichlet type
  class VolPointBdCond : public BlockPointBdCond
  {
  public:
    // Constructor
    VolPointBdCond(int face_nmb, double param[], Point& condition_value);

    // Destructor
    ~VolPointBdCond();

    // Get the enumeration of affected surface coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) const;

    // Get the value of the condition
    virtual Point getConditionValue() const;

    // Get the parameter corresponding to the condition value
    double* getParam() const;

    // Get enumeration of affected surface coefficients and the factors in the
    // equation for requiring the current boundary to interpolate the condition
    virtual void 
      getInterpolationFactors(std::vector<std::pair<int,double> >& factors) const;

  private:
    // The boundary surface it corresponds to
    // facenmb_ = 0: the boundary corresponding to the minimum parameter in the first parameter
    //                 direction (u_min)
    // facenmb_ = 1: the boundary corresponding to the maximum parameter in the first parameter
    //                 direction (u_max)
    // facenmb_ = 2: the boundary corresponding to the minimum parameter in the second parameter
    //                 direction (v_min)
    // facenmb_ = 3: the boundary corresponding to the maximum parameter in the second parameter
    //                 direction (v_max)
    // facenmb_ = 4: the boundary corresponding to the minimum parameter in the third parameter
    //                 direction (w_min)
    // facenmb_ = 5: the boundary corresponding to the maximum parameter in the third parameter
    //                 direction (w_max)
    int facenmb_; 
  
    // The value of the boundary condition of the solution in the given point
    Point value_;

    // The parameter value corresponding to this pointwise condition
    double param_[2];

    // Pointer to the block solution to which this boundary condition belongs
    VolSolution* parent_;

  };   // end class VolPointBdCond

} // end namespace Go


#endif    // #ifndef __VOLPOINTBDCOND_H
