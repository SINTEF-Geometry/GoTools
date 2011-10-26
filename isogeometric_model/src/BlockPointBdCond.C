//===========================================================================
//
// File : BlockPointBdCond.C
//
// Created: Mon Mar  1 08:57:56 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/isogeometric_model/BlockPointBdCond.h"


namespace Go
{

  //===========================================================================
  BlockPointBdCond::BlockPointBdCond()
  //===========================================================================
  {
  }

  //===========================================================================
  BlockPointBdCond::~BlockPointBdCond()
  //===========================================================================
  {
  }

  //===========================================================================
  int BlockPointBdCond::getSolutionSpaceIdx() const
  //===========================================================================
  {
    MESSAGE("getSolutionSpaceIdx() not implemented");
    return 0;
  }

}   // namespace Go
