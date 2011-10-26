//===========================================================================
//
// File : SfPointBdCond.C
//
// Created: Thu Aug  5 10:00:46 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================






#include "GoTools/isogeometric_model/SfPointBdCond.h"


using namespace std;


namespace Go
{

  //===========================================================================
  SfPointBdCond::SfPointBdCond(int edge_nmb, double param, Point& condition_value)
  //===========================================================================
  {
    MESSAGE("SfPointBdCond constructor() not implemented");
  }

  //===========================================================================
  SfPointBdCond::~SfPointBdCond()
  //===========================================================================
  {
    MESSAGE("SfPointBdCond destructor() not implemented");
  }

  //===========================================================================
  void SfPointBdCond::getCoefficientsEnumeration(vector<int>& local_enumeration) const
  //===========================================================================
  {
    MESSAGE("getCoefficientsEnumeration() not implemented");
  }

  //===========================================================================
  Point SfPointBdCond::getConditionValue() const
  //===========================================================================
  {
    MESSAGE("getConditionValue() not implemented");
    Point p(0.0, 0.0, 0.0);
    return p;
  }

  //===========================================================================
  double SfPointBdCond::getParam() const
  //===========================================================================
  {
    MESSAGE("getParam() not implemented");
    return 0.0;
  }

  //===========================================================================
  void SfPointBdCond::getInterpolationFactors(vector<pair<int,double> >& factors) const
  //===========================================================================
  {
    MESSAGE("getInterpolationFactors() not implemented");
  }

  int SfPointBdCond::edgeNumber() const
  //===========================================================================
  {
    return edgenmb_;
  }

}   // namespace Go
