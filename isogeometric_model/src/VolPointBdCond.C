//===========================================================================
//                                                                           
// File: VolPointBdCond.C                                                    
//                                                                           
// Created: Mon Aug 13 17:23:40 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/isogeometric_model/VolPointBdCond.h"


using namespace std;


namespace Go
{

  //===========================================================================
  VolPointBdCond::VolPointBdCond(int face_nmb, double param[], Point& condition_value)
  //===========================================================================
  {
    MESSAGE("VolPointBdCond constructor() not implemented");
  }

  //===========================================================================
  VolPointBdCond::~VolPointBdCond()
  //===========================================================================
  {
    MESSAGE("VolPointBdCond destructor() not implemented");
  }

  //===========================================================================
  void VolPointBdCond::getCoefficientsEnumeration(vector<int>& local_enumeration) const
  //===========================================================================
  {
    MESSAGE("getCoefficientsEnumeration() not implemented");
  }

  //===========================================================================
  Point VolPointBdCond::getConditionValue() const
  //===========================================================================
  {
    MESSAGE("getConditionValue() not implemented");
    Point p(0.0, 0.0, 0.0);
    return p;
  }

  //===========================================================================
  double* VolPointBdCond::getParam() const
  //===========================================================================
  {
    MESSAGE("getParam() not implemented");
    return NULL;
  }

  //===========================================================================
  void VolPointBdCond::getInterpolationFactors(vector<pair<int,double> >& factors) const
  //===========================================================================
  {
    MESSAGE("getInterpolationFactors() not implemented");
  }

  int VolPointBdCond::faceNumber() const
  //===========================================================================
  {
    return facenmb_;
  }

}   // namespace Go
