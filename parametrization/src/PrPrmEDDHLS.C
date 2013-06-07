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

#include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmEDDHLS::makeWeights(int i)
//-----------------------------------------------------------------------------
//   Calculate Eck, DeRose, Duchamp, Hoppe, Lounsbery, Stuetzle weights
//   for the interior node i of the graph.
//   The formula for the weight of edge (i,j) is
//   w_jk = (L_{i,k_1}^2 + L_{j,k_1}^2 - L_{i,j}^2) / Area_{i,j,k_1}
//             +
//          (L_{i,k_2}^2 + L_{j,k_2}^2 - L_{i,j}^2) / Area_{i,j,k_2}.
//
//   These weights are not necessarily positive!!!
//   So it's better to use the class PrPrmShpPres.
//   M.F. Mar. 97.
{
  weights_.clear();
  int n = (int)neighbours_.size();

  int j;
  double tot_weight = 0.0, weight;
  for(j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    double area1 = area(a,b,c);

    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D d = g_->get3dNode(j2);
    double area2 = area(a,c,d);

    double l1 = a.dist2(c);
    double l2 = a.dist2(b);
    double l3 = c.dist2(b);
    double l4 = a.dist2(d);
    double l5 = c.dist2(d);

    weight = (l2+l3-l1) / area1 + (l4+l5-l1) / area2;
    weights_.push_back(weight);
#ifdef PRDEBUG
    if(weight < 0.0)
    {
      cout << "negative weight=" << weight;
      cout << " for edge(" << i << "," << jj << ")" << "\n";
    }
#endif
    tot_weight += weight;
  }

#ifdef PRDEBUG
    if(tot_weight < 0.0)
    {
      cout << "negative sum of weight=" << tot_weight;
      cout << " for node " << i << "\n";
    }
#endif

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++)
  {
    weights_[j] *= ratio;
#ifdef PRDEBUG
    if(weights_[j] < 0.0)
    {
      cout << "negative weight=" << weights_[j];
      cout << " for edge(" << i << "," << neighbours_[j] << ")" << "\n";
    }
#endif
  }

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmEDDHLS::PrPrmEDDHLS()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmEDDHLS::~PrPrmEDDHLS()
//-----------------------------------------------------------------------------
{
}

