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


#include "GoTools/creators/HermiteGrid2D.h"
#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalSurface.h"

using namespace std;

namespace Go
{


HermiteGrid2D::HermiteGrid2D(const EvalSurface& sf,
                             double u1, double u2, double v1, double v2)
    : dim_(sf.dim()), MM_(2), NN_(2), elem_size_(4), index_u_(0), index_v_(0)
{
  // Copy knot vector

  knots_u_.reserve(MM_);
  knots_u_.push_back(u1);
  knots_u_.push_back(u2);

  knots_v_.reserve(NN_);
  knots_v_.push_back(v1);
  knots_v_.push_back(v2);

  // Calculate the curve values at the parameter grid

  array_.reserve(elem_size_*MM_*NN_);
  Point derive[3];
  for (int kj = 0; kj < NN_; ++kj)
  {
      for (int ki = 0; ki < MM_; ++ki)
      {
          sf.eval(knots_u_[ki], knots_v_[kj], 1, derive);
          for (int kk = 0; kk < elem_size_; ++kk)
          {
              array_.push_back(derive[kk]);
          }
      }
  }
}

HermiteGrid2D::HermiteGrid2D(const EvalSurface& sf,
                             double param_u[], double param_v[], int mm, int nn)
    : dim_(sf.dim()), MM_(mm), NN_(nn), elem_size_(2), index_u_(mm/2), index_v_(nn/2)
//--------------------------------------------------------
//  Constructor
//
// INPUT:
//      param  - Array of strictly increasing parameters in the parameter
//               interval of "sf". First and last point defines the interval.
//      n      - Number of elements in "param"
//      sf    - Surface to sample
//--------------------------------------------------------
{
  // Check that knots are strictly increasing
  int i, j;
  for (i=1; i<mm; i++)
      if (param_u[i] <= param_u[i-1])
          THROW("Input grid illegal");

  for (i=1; i<nn; i++)
      if (param_v[i] <= param_v[i-1])
          THROW("Input grid illegal");


  // Calculate the curve values at the parameter grid

  array_.reserve(elem_size_*mm*nn);

  Point derive[2];
  for (j=0; j<nn; j++)
  {
      for (i=0; i<mm; i++)
      {
          sf.eval(param_u[i], param_v[j], 1, derive);
          for (j=0; j<elem_size_; j++)
          {
              array_.push_back(derive[j]);
          }
      }
  }

  // Copy knot vectors.

  knots_u_.reserve(mm);
  for (i=0; i<mm; i++)
    knots_u_.push_back(param_u[i]);

  knots_v_.reserve(nn);
  for (i=0; i<nn; i++)
    knots_v_.push_back(param_v[i]);
}

HermiteGrid2D::~HermiteGrid2D()
//--------------------------------------------------------
//  Destructor
//--------------------------------------------------------
{}


int HermiteGrid2D::addKnot(const EvalSurface& sf, double knot, bool dir_is_u)
//--------------------------------------------------------------------
// PURPOSE: Insert a new knot in the knotvector . Also the value and tangent 
//          of the curve at this knot is added to the Hermite grid
//
// INPUT:
//      sf	 - Surface to evaluate
//      knot	 - New knot
//      dir_is_u - True if we refine in the u-dir.
// OUTPUT:
//      addKnot() - The index of the new knot in the (sorted) knot vector
//                  after insertion. The first index for this knotvector is 0.
//--------------------------------------------------------------------
{
    // Evaluate the sf at the new knot and insert the values into
    // the Hermite grid.

    int index = getPosition(knot, dir_is_u);
    if (dir_is_u)
    {
        index_u_ = index;
    }
    else
    {
        index_v_ = index;
    }

    int num_knots_opp_dir = (dir_is_u) ? NN_ : MM_;
    // We run through all knot values in the opposite direction.
    for (size_t ki = 0; ki < num_knots_opp_dir; ++ki)
    {

        double knot_u = (dir_is_u) ? knot : knots_u_[ki];
        double knot_v = (dir_is_u) ? knots_v_[ki] : knot;
        Point derive[3];
        sf.eval(knot_u, knot_v, 1, derive);

        // Insert the new knot into the knot vector

        if (dir_is_u)
        {
            knots_u_.insert(knots_u_.begin() + index_u_ + 1, 1, knot);
            MM_++;

        }
        else
        {
            knots_v_.insert(knots_v_.begin() + index_v_ + 1, 1, knot);
            NN_++;
        }
      
        array_.insert(array_.begin()+elem_size_*(index+1), 1, derive[0]);
        array_.insert(array_.begin()+elem_size_*(index+1)+1, 1, derive[1]);
        array_.insert(array_.begin()+elem_size_*(index+1)+2, 1, derive[2]);
    }
    
    return index;
}


void HermiteGrid2D::getSegment(int left1, int right1,
                               int left2, int right2,
                               double& spar1, double& epar1,
                               double& spar2, double& epar2,
                               Point bezcoef[4])
//--------------------------------------------------------------------
// PURPOSE: Calculate Bezier coefficients of cubic curve interpolating the
//          Hermite values at grid nodes with indeces "left" and "right"
//
// INPUT:
//     left1 - indicating grid node for start of surface segment in u-dir
//     right1 - indicating grid node for end of surface segment in u-dir
//     left2 - indicating grid node for start of surface segment in v-dir
//     right2 - indicating grid node for end of surface segment in v-dir
// OUTPUT:
//      spar    - start parameter of segment
//      epar    - end parameter of segment
//      bezcoef - array of cubic Bezier coefficients
//--------------------------------------------------------------------
{
    spar1 = knots_u_[left1];
    epar1 = knots_u_[right1];
    spar2 = knots_v_[left2];
    epar2 = knots_v_[right2];

#if 1
    MESSAGE("Under construction!");
#else



  for (i=0; i<dim_; i++)
  {
    bezcoef[i] = sder1[i];
    bezcoef[dim_+i] = sder1[i]+sder1[dim_+i]*scale1;
    bezcoef[sdim2+i] = eder1[i]-eder1[dim_+i]*scale1;
    bezcoef[sdim3+i] = eder1[i];
  }
  for (i=0; i<dim_; i++)
  {
    bezcoef[sdim4+i] = sder1[i] + sder1[2*dim_+i]*scale2;
    bezcoef[sdim4+dim_+i] = bezcoef[dim_+i] + sder1[2*dim_+i]*scale2 +
      sder1[3*dim_+i]*scale1*scale2;
    bezcoef[sdim4+sdim2+i] = bezcoef[2*dim_+i] + eder1[2*dim_+i]*scale2 - 
      eder1[3*dim_+i]*scale1*scale2;
    bezcoef[sdim4+sdim3+i] = eder1[i] + eder1[2*dim_+i]*scale2;
  }
  for (i=0; i<dim_; i++)
  {
    bezcoef[3*sdim4+i] = sder2[i];
    bezcoef[3*sdim4+dim_+i] = sder2[i]+sder2[dim_+i]*scale1;
    bezcoef[3*sdim4+sdim2+i] = eder2[i]-eder2[dim_+i]*scale1;
    bezcoef[3*sdim4+sdim3+i] = eder2[i];
  }
  for (i=0; i<dim_; i++)
  {
    bezcoef[2*sdim4+i] = sder2[i] - sder2[2*dim_+i]*scale2;
    bezcoef[2*sdim4+dim_+i] = bezcoef[3*sdim4+dim_+i] - sder2[2*dim_+i]*scale2 -
      sder2[3*dim_+i]*scale1*scale2;
    bezcoef[2*sdim4+sdim2+i] = bezcoef[3*sdim4+sdim2+i] - eder2[2*dim_+i]*scale2 + 
      eder2[3*dim_+i]*scale1*scale2;
    bezcoef[2*sdim4+sdim3+i] = eder2[i] - eder2[2*dim_+i]*scale2;
  }

  
    double scale = (epar - spar)/3.0;
    bezcoef[0] = array_[elem_size_*left];
    bezcoef[3] = array_[elem_size_*right];
    bezcoef[1] = bezcoef[0] + array_[elem_size_*left+1]*scale;
    bezcoef[2] = bezcoef[3] - array_[elem_size_*right+1]*scale;

#endif
    
    return;
}


int HermiteGrid2D::getPosition(double knot, bool dir_is_u)
//---------------------------------------------------------
// PURPOSE: Find the index into the knot vector of the parameter knot
//----------------------------------------------------------
{
    const vector<double>& knots = (dir_is_u) ? knots_v_ : knots_v_;
    int index = (dir_is_u) ? index_u_ : index_v_;

    int id1=0, id2=(int)knots.size()-1;
    if ((index == id2 && knot >= knots[index]) ||
        (knots[index] <= knot && knot < knots[index]))
        return index;

    // Search
    index = id1; // We must reset index_ before we start searching.
    while (id2 > id1+1)
    {
        if (knots[index] > knot)
            id2 = index;
        if (knots[index] <= knot)
            id1 = index;
        index = (id1 + id2)/2;
    }

    if (dir_is_u)
    {
        index_u_ = index;
    }
    else
    {
        index_v_ = index;
    }
      
  
    return index;
}


}
