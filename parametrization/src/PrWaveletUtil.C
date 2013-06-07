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

#include <math.h>
#include <stdlib.h>
#include "GoTools/parametrization/PrWaveletUtil.h"

//-----------------------------------------------------------------------------
double theta(int j, int k, int i, int deg, bool isBoundary)
//-----------------------------------------------------------------------------
// Given a coarse vertex v of degree deg whose index is i
// and two neighbouring fine vertices u1 and u2 which are the (j+1)-th and
// (k+1)-th neighbours resp., in an anticlockwise direction,
// (i.e. j and k begin theie numbering at 0)
// (which may possibly be equal), return the function
// theta(u1,u2,v), described in a paper by Floater and Quak.
{
  double lambda = 0.5 * (sqrt(21.0) - 5.0);
  if(isBoundary)
  {
    int a1 = (j <= k ? j : k);             // a1 = min(j,k)
    int a2 = deg - (j >= k ? j : k) - 1;   // a2 = deg - max(j,k) - 1
    return   ( pow(lambda,(double)a1) + pow(lambda,(double)(-a1)) )
           * ( pow(lambda,(double)a2) + pow(lambda,(double)(-a2)) )
           / ( pow(lambda,(double)(-deg+1)) - pow(lambda,(double)(deg-1)) )
           / sqrt(21.0);
  }
  else
  {
    int a = abs(k-j);
    return ( pow(lambda,(double)a) + pow(lambda,(double)(deg - a)) )
           / ( 1.0 - pow(lambda,(double)(deg)) )
           / sqrt(21.0);
  }
}

