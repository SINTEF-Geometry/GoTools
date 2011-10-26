//===========================================================================
//                                                                           
// File: PrPrmWachspress.C                                                   
//                                                                           
// Created: Mon Jun 30 15:19:28 2003                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PrPrmWachspress.C,v 1.4 2007-03-02 16:35:52 jbt Exp $
//                                                                           
//===========================================================================


#include "GoTools/parametrization/PrPrmWachspress.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmWachspress::makeWeights(int i)
//-----------------------------------------------------------------------------
{
  weights_.clear();
  int n = (int)neighbours_.size();

  int j;
  double tot_weight = 0.0, weight;
  for(j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    Vector3D d = g_->get3dNode(j2);

    weight = (cotangent(b, c, a) + cotangent(a, c, d))/a.dist2(c);
    weights_.push_back(weight);

    tot_weight += weight;
  }

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++)
  {
    weights_[j] *= ratio;
  }

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmWachspress::PrPrmWachspress()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmWachspress::~PrPrmWachspress()
//-----------------------------------------------------------------------------
{
}

