/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrNestedTriangulation.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Feb 1998
 DESCRIPTION : Implementation of methods in the class PrOrganizePoints.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrNestedTriangulation.h"

//----------------------------------------------------------------------------
void PrNestedTriangulation::getParents(int i, int jlev, int& p1, int& p2)
//-----------------------------------------------------------------------------
// Assuming that node i belongs to V^j\V^{j-1} and that
// jlev >= 1, find i's two parent nodes in V^{j-1} p1 and p2.
{
  vector<int> neighbours;
  getNeighbours(i,jlev,neighbours);
  if(isBoundary(i))
  {
    // Then there are exactly four neighbours. We want the
    // first and last ones.
    p1 = neighbours[0];
    p2 = neighbours[3];
    return;
  }
  else
  {
    // Then there are exactly six neighbours. Find the first
    // coarse one.
    if(neighbours[0] < getNumNodes(jlev-1))
    {
      p1 = neighbours[0];
      p2 = neighbours[3];
      return;
    }
    else if(neighbours[1] < getNumNodes(jlev-1))
    {
      p1 = neighbours[1];
      p2 = neighbours[4];
      return;
    }
    else
    {
      p1 = neighbours[2];
      p2 = neighbours[5];
      return;
    }
  }
}

