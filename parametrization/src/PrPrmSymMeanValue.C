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

#include "GoTools/parametrization/PrPrmSymMeanValue.h"
#include "GoTools/parametrization/PrParamUtil.h"

// PRIVATE METHODS

//-----------------------------------------------------------------------------
bool
PrPrmSymMeanValue::makeWeights(int i)
//-----------------------------------------------------------------------------
//  Calculate mean value weights for the
//  interior node i of the graph
//  It is assumed here that the indices of the neighbours of i
//  have already been stored in neighbours_.
//  This is an attempt to improve the mean value method of
//  PrPrmMeanValue by making it symmetric.
{
    THROW("This method is bad, and suitable for experimentation only.");

  weights_.clear();
  int n = (int)neighbours_.size();

  //  std::cout << "Testing" << std::endl;


  int j;
  double tot_weight = 0.0, weight;
  for (j=0; j<n; j++)
  {
    int jj = neighbours_[j];
    int jprev = (j == 0 ? n-1 : j-1);
    int j1 = neighbours_[jprev];
    int jnext = (j == n-1 ? 0 : j+1);
    int j2 = neighbours_[jnext];
    Vector3D d = g_->get3dNode(j2);
    Vector3D a = g_->get3dNode(i);
    Vector3D b = g_->get3dNode(j1);
    Vector3D c = g_->get3dNode(jj);
    double t1 = tanThetaOverTwo(a,b,c);
    //std::cout << "t1 = " << t1 << endl;

    double t2 = tanThetaOverTwo(a,c,d);
    //std::cout << "t2 = " << t2 << endl;

    double t3 = tanThetaOverTwo(c,d,a);
    //std::cout << "t3 = " << t2 << endl;

    double t4 = tanThetaOverTwo(c,a,b);
    //std::cout << "t4 = " << t2 << endl;

    double len = a.dist(c);
    //std::cout << "len = " << len << endl;

    weight = (t1 + t2 + t3 + t4) / len;
    //std::cout << "weight = " << weight << endl;
    weights_.push_back(weight);
    tot_weight += weight;
  }

  // Scale the weights so that they sum to 1.

  double ratio = 1.0 / tot_weight;
  for(j=0; j<n; j++) weights_[j] *= ratio;

  return true;
}

//-----------------------------------------------------------------------------
double
PrPrmSymMeanValue::tanThetaOverTwo(Vector3D& a, Vector3D& b, Vector3D& c)
//-----------------------------------------------------------------------------
// Return tangent of half the angle between vectors b-a and c-a
// without using trig functions.
// Use fact that tan(alpha/2) = (1-cos(alpha)) / sin(alpha).
// and use scalar and dot products to get cos(alpha) and sin(alpha).
//  M.F. Apr. 2002.
{
  Vector3D bb = b - a;
  Vector3D cc = c - a;

  double cp = (bb.cross(cc)).length();
       // length of cross product of bb and cc
  double dp = bb * cc; // scalar product of bb and cc
  double bc = bb.length() * cc.length();
  return (bc - dp) / cp;

  return true;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrPrmSymMeanValue::PrPrmSymMeanValue()
: PrParametrizeInt()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
PrPrmSymMeanValue::~PrPrmSymMeanValue()
//-----------------------------------------------------------------------------
{
}

