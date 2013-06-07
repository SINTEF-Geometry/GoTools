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

#include "GoTools/parametrization/PrPrmSurface.h"
#include "GoTools/parametrization/PrParamUtil.h"
#include "GoTools/utils/Values.h"
#ifdef _WIN32
#define M_PI 3.14159265358979
#endif

#ifdef __BORLANDC__
using std::cos;
using std::sin;
#endif

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmSurface::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate shape-preserving weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is the third parametrization in the article:
//  M. S. Floater, "Parametrization and smooth approximation of
//  surface triangulations", to appear in CAGD, 1997.
{
  weights_.clear();
  weightsold_.clear();
  int n = (int)neighbours_.size();

  weights_.resize(n, 0.0);

  // Find local u,v's
  localParamXYZ(i);
  makeWeights(u_,v_,weights_);

  localParamUV(i);
  makeWeights(uold_,vold_,weightsold_);

  double sum = 0.0;
  int j;
  for(j=0; j<n; j++)
  {
    //weights_[j] = weightsold_[j] / sqrt(weights_[j]);
    weights_[j] = weightsold_[j] / weights_[j];
    sum += weights_[j];
  }

  // Scale the weights so that they sum to 1.

  for(j=0; j<n; j++) weights_[j] /= sum;

  return true;
}

//-----------------------------------------------------------------------------
bool
PrPrmSurface::makeWeights(vector<double>& u,
			  vector<double>& v,
			  vector<double>& weights)
//-----------------------------------------------------------------------------
{
  int n = (int)u.size();
  weights.clear();
  weights.resize(n, 0.0);

  int j;
  for(j=0; j<n; j++)
  {
    /* Given the j-th neighbour of node i,
       find the two neighbours by intersecting the
       line through nodes i and j with all segments of the polygon
       made by the neighbours. Take the two neighbours on
       either side. Only one segment intersects this line. */

    int k;
    for(k=0; k<n; k++)
    {
      int kk = (k == (n-1) ? 0 : k+1);
      if(k == j || kk == j) continue;

      double cross1 = det(u[j],v[j],u[k],v[k]);
      double cross2 = det(u[j],v[j],u[kk],v[kk]);

      if(cross1 * cross2 <= 0.0)
      {
        double tau0,tau1,tau2;
        baryCoords0(u[j],v[j],u[k],v[k],u[kk],v[kk],tau0,tau1,tau2);
        weights[j]  += tau0;
        weights[k]  += tau1;
        weights[kk] += tau2;
        break;
      }
    }
  }

  // Scale the weights so that they sum to 1.
  // Before this, they sum to n

  double ratio = 1.0 / (double)n;
  for(j=0; j<n; j++) weights[j] *= ratio;

  // Temporary code.
  /*
  ratio = 0.0;
  for(j=1; j<=n; j++)
  {
    weights[j]  = 1.0 / weights[j];
    ratio += weights[j];
  }
  for(j=1; j<=n; j++) weights[j] /= ratio;

  double tmp;

  tmp = weights(1);
  weights(1) = weights(3);
  weights(3) = tmp;

  tmp = weights(2);
  weights(2) = weights(4);
  weights(4) = tmp;
  */
  // End of temporary code.

  return true;
}

//-----------------------------------------------------------------------------
bool
PrPrmSurface::localParamXYZ(int i)
//-----------------------------------------------------------------------------
//  Make a local parametrization for the interior
//  node i of the given graph.
//  This is based on a discretization of the geodesic polar map.
//  The distances of the neighbours from i are preserved
//  and the ratios of any two interior angles are also preserved
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  M.F. Mar. 97
{
  alpha_.clear();
  len_.clear();
  u_.clear();
  v_.clear();
  int n = (int)neighbours_.size();

  double alpha_sum = 0.0;
  double alpha;

  int j;
  for (j=0; j<n; j++)
  {
    int j1 = neighbours_[j];
    Vector3D v1 = g_->get3dNode(j1) - g_->get3dNode(i);

    // get previous node
    int jprev = (j == 0 ? n-1 : j-1);
    int j2 = neighbours_[jprev];
    Vector3D v2 = g_->get3dNode(j2) - g_->get3dNode(i);

    len_.push_back(v1.length());
    //len_[j] = sqrt(v1.norm()); // alternative
    //len_[j] = 1.0; // alternative
    alpha = v1.angle(v2);
    alpha_.push_back(alpha);
    //alpha_[j] = sqrt(v1.angle(v2)); // alternative
    //alpha_[j] = 1.0; // alternative
    alpha_sum += alpha;
  }

    double factor = 2.0 * M_PI / alpha_sum;

    for(j=0; j<n; j++) alpha_[j] *= factor;

    alpha_[0] = 0.0;
    for(j=1; j<n; j++) alpha_[j] += alpha_[j-1];

    for(j=0; j<n; j++)
    {
    //v_[j].init(len_[j] * cos(alpha_[j]), len_[j] * sin(alpha_[j]));
    //v_.push_back(CgVector2d(len_[j] * cos(alpha_[j]), len_[j] * sin(alpha_[j])));
      u_.push_back(len_[j] * cos(alpha_[j]));
      v_.push_back(len_[j] * sin(alpha_[j]));
    }

  return true;
}

//-----------------------------------------------------------------------------
bool
PrPrmSurface::localParamUV(int i)
//-----------------------------------------------------------------------------
//  Make a local parametrization for the interior
//  node i of the given graph.
//  This is based on a discretization of the geodesic polar map.
//  The distances of the neighbours from i are preserved
//  and the ratios of any two interior angles are also preserved
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  M.F. Mar. 97
{
  uold_.clear();
  vold_.clear();
  int n = (int)neighbours_.size();

  for (int j=0; j<n; j++)
  {
    int j1 = neighbours_[j];
    uold_.push_back(g_->getU(j1) - g_->getU(i));
    vold_.push_back(g_->getV(j1) - g_->getV(i));
  }

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmSurface::PrPrmSurface()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmSurface::~PrPrmSurface()
//-----------------------------------------------------------------------------
{
}

