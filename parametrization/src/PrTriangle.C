/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrTriangle.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Aug. 97
 DESCRIPTION : Implementation of methods in the class PrTriangle.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrTriangle.h"
#include <iostream>

// MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
int PrTriangle::getOppositeTriangle(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node) return t1_;
  else if(n2_ == node) return t2_;
  else if(n3_ == node) return t3_;
  else return -1;
}

//-----------------------------------------------------------------------------
int PrTriangle::getLeftTriangle(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node) return t2_;
  else if(n2_ == node) return t3_;
  else if(n3_ == node) return t1_;
  else return -1;
}

//-----------------------------------------------------------------------------
int PrTriangle::getRightTriangle(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node) return t3_;
  else if(n2_ == node) return t1_;
  else if(n3_ == node) return t2_;
  else return -1;
}

//-----------------------------------------------------------------------------
int PrTriangle::getAnticlockwiseNode(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node) return n2_;
  else if(n2_ == node) return n3_;
  else if(n3_ == node) return n1_;
  else return -1;
}

//-----------------------------------------------------------------------------
int PrTriangle::getClockwiseNode(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node) return n3_;
  else if(n2_ == node) return n1_;
  else if(n3_ == node) return n2_;
  else return -1;
}

//-----------------------------------------------------------------------------
bool PrTriangle::isVertex(int node) const
//-----------------------------------------------------------------------------
{
  if(n1_ == node || n2_ == node || n3_ == node) return true;
  else return false;
}

//-----------------------------------------------------------------------------
void PrTriangle::replaceNode(int n1, int n2)
//-----------------------------------------------------------------------------
{
  if(n1_ == n1) n1_ = n2;
  else if(n2_ == n1) n2_ = n2;
  else if(n3_ == n1) n3_ = n2;
}

//-----------------------------------------------------------------------------
void PrTriangle::replaceTriangle(int t1, int t2)
//-----------------------------------------------------------------------------
{
  if(t1_ == t1) t1_ = t2;
  else if(t2_ == t1) t2_ = t2;
  else if(t3_ == t1) t3_ = t2;
}

//-----------------------------------------------------------------------------
int PrTriangle::getEdge(int triangle, int &n1, int &n2) const
//-----------------------------------------------------------------------------
{
   if (t1_==triangle)
   {
      n1=n2_;
      n2=n3_;
   } else if (t2_==triangle)
   {
      n1=n3_;
      n2=n1_;
   } else if (t3_==triangle)
   {
      n1=n1_;
      n2=n2_;
   } else
      return false;
   return true;
}



//-----------------------------------------------------------------------------
void PrTriangle::print(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  os << n1_ << ' ' << n2_ << ' ' << n3_ << "    "
     << t1_ << ' ' << t2_ << ' ' << t3_ << std::endl;
}
//-----------------------------------------------------------------------------

void PrTriangle::scan(std::istream& is)
//-----------------------------------------------------------------------------
{
  is >> n1_ >> n2_ >> n3_ >> t1_ >> t2_ >> t3_;
}
