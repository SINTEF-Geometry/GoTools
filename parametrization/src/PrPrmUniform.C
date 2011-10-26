/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrPrmUniform.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Mar. 97
 DESCRIPTION : Implementation of methods in the class PrPrmUniform.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPrmUniform.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmUniform::makeWeights(int)
//-----------------------------------------------------------------------------
//  Calculate uniform weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is the first parametrization in the article:
//  M. S. Floater, "Parametrization and smooth approximation of
//  surface triangulations", CAGD 14 (1997), 231--250.
{
  weights_.clear();
  int n = (int)neighbours_.size();
  double oneOverDegree = 1.0 / (double)n;
  for(int j=1; j<=n; j++) weights_.push_back(oneOverDegree);
  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmUniform::PrPrmUniform()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmUniform::~PrPrmUniform()
//-----------------------------------------------------------------------------
{
}

