/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2002 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrPrmMeanValue.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Apr. 2002
 DESCRIPTION : Implementation of methods in the class PrPrmMeanValue.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPrmMeanValue.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmMeanValue::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate mean value weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is the parametrization in the article:
//  M. S. Floater, "Mean Value Coordinates", preprint 2002,
//  which may sometimes be visually smoother
//  than the "shape-preserving" parametrization.
{
  weights_.clear();
  int n = (int)neighbours_.size();

  //  std::cout << "Testing" << std::endl;


  int j;
  double tot_weight = 0.0, weight;
  for (j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    double t1 = tanThetaOverTwo(a,b,c);
    //std::cout << "t1 = " << t1 << endl;

    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D d = g_->get3dNode(j2);
    double t2 = tanThetaOverTwo(a,c,d);
    //std::cout << "t2 = " << t2 << endl;

    double len = a.dist(c);
    //std::cout << "len = " << len << endl;

    weight = (t1 + t2) / len;
    //std::cout << "weight = " << weight << endl;
    weights_.push_back(weight);
    tot_weight += weight;
  }

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++) weights_[j] *= ratio;

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmMeanValue::PrPrmMeanValue()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmMeanValue::~PrPrmMeanValue()
//-----------------------------------------------------------------------------
{
}

