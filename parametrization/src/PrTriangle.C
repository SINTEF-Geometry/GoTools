/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
