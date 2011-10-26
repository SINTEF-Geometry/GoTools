/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrPrmLeastSquare.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Mar. 97
 DESCRIPTION : Implementation of methods in the class PrPrmLeastSquare.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPrmLeastSquare.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmLeastSquare::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate weights according to a weighted least squares method
//  for the interior node i of the graph.
//  Here we weight by the reciprocal of edge length.
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is the second parametrization in the article:
//  M. S. Floater, "Parametrization and smooth approximation of
//  surface triangulations", to appear in CAGD, 1997.
{
  weights_.clear();
  int n = (int)neighbours_.size();

  double weight_sum = 0.0, weight;

  int j;
  for (j=0; j<n; j++)
  {
    int j1 = neighbours_[j];
    weight = (g_->get3dNode(i)).dist(g_->get3dNode(j1));
    if(weight == 0.0) weight = 1.0;
    else weight = 1.0 / weight;
    weights_.push_back(weight);
    weight_sum += weight;
  }

  // Scale the weights so that they sum to 1.
  for(j=0; j<n; j++) weights_[j] /= weight_sum;

  return 1;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmLeastSquare::PrPrmLeastSquare()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmLeastSquare::~PrPrmLeastSquare()
//-----------------------------------------------------------------------------
{
}


