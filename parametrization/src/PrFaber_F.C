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

