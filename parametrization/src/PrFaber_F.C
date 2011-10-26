/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrFaber_F.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Dec. 98
 DESCRIPTION : Implementation of methods in the class PrFaber_F.
*********************************************************************/

#include "GoTools/parametrization/PrFaber_F.h"

//-----------------------------------------------------------------------------
void
PrFaber_F::decompose(int jlev, int dim)
//-----------------------------------------------------------------------------
// The nodes in the triangulation belonging to V^j
// form a piecewise linear function f^j in the space S^j.
// This routines decomposes f^j
// into f^{j-1} + g^{j-1} where
// f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
{
  if(dim != 1 && dim != 3) return;
  if(jlev < 1 || jlev > t_->getFinestLevel()) return;

  int nprev = t_->getNumNodes(jlev-1);

  double newval;
  for(int i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    for(size_t j=0; j< neighbours_.size(); j++)
    {
      int jj = neighbours_[j];
      if(dim == 3)
      {
        newval = t_->getX(jj) - 0.5 * t_->getX(i);
        t_->setX(jj,newval);
        newval = t_->getY(jj) - 0.5 * t_->getY(i);
        t_->setY(jj,newval);
      }
      newval = t_->getZ(jj) - 0.5 * t_->getZ(i);
      t_->setZ(jj,newval);
    }
  }
}

//-----------------------------------------------------------------------------
void
PrFaber_F::compose(int jlev, int dim)
//-----------------------------------------------------------------------------
// The nodes in the triangulation belonging to V^j
// form a piecewise linear function f^j in the space S^j.
// This routines composes f^j
// from f^{j-1} + g^{j-1} where
// f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
{
  if(dim != 1 && dim != 3) return;
  if(jlev < 1 || jlev > t_->getFinestLevel()) return;

  int nprev = t_->getNumNodes(jlev-1);

  double newval;
  for(int i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    for(size_t j=0; j< neighbours_.size(); j++)
    {
      int jj = neighbours_[j];
      if(dim == 3)
      {
        newval = t_->getX(jj) + 0.5 * t_->getX(i);
        t_->setX(jj,newval);
        newval = t_->getY(jj) + 0.5 * t_->getY(i);
        t_->setY(jj,newval);
      }
      newval = t_->getZ(jj) + 0.5 * t_->getZ(i);
      t_->setZ(jj,newval);
    }
  }
}

