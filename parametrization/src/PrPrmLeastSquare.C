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


