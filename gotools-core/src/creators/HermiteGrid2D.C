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
  Point derive[4]; // pos, 2*der, twist.
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

  no_split_status_.resize(MM_*NN_, 0);
  
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
// PURPOSE: Insert a new knot in the knot vector . Also the value and tangent 
//          of the curve at this knot is added to the Hermite grid
//
// INPUT:
//      sf	 - Surface to evaluate
//      knot	 - New knot
//      dir_is_u - True if we refine in the u-dir.
// OUTPUT:
//      addKnot() - The index of the new knot in the (sorted) knot vector
//                  after insertion. The first index for this knot vector is 0.
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
    for (int ki = 0; ki < num_knots_opp_dir; ++ki)
    {

        double knot_u = (dir_is_u) ? knot : knots_u_[ki];
        double knot_v = (dir_is_u) ? knots_v_[ki] : knot;
        Point derive[4];
        sf.eval(knot_u, knot_v, 1, derive);
        int index_2d = (dir_is_u) ? ki*(MM_ + 1) + index + 1: (index + 1)*MM_ + ki;
        array_.insert(array_.begin() + elem_size_*index_2d, 1, derive[0]);
        array_.insert(array_.begin() + elem_size_*index_2d + 1, 1, derive[1]);
        array_.insert(array_.begin() + elem_size_*index_2d + 2, 1, derive[2]);
        array_.insert(array_.begin() + elem_size_*index_2d + 3, 1, derive[3]);

        int index_2d_nb = (dir_is_u) ? ki*(MM_ + 1) + index : index*MM_ + ki;
        int no_split_status = no_split_status_[index_2d_nb];
        no_split_status_.insert(no_split_status_.begin() + index_2d, no_split_status);
    }

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

    return index;
}


void HermiteGrid2D::getSegment(int left1, int right1,
                               int left2, int right2,
                               double& spar1, double& epar1,
                               double& spar2, double& epar2,
                               Point bezcoef[16])
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

    const double scale1 = (epar1 - spar1)/3.0;
    const double scale2 = (epar2 - spar2)/3.0;

    // const int sdim4 = dim_*4;
    // const int sdim3 = dim_*3;
    // const int sdim2 = dim_*2;
//    std::cout << "array_.size(): " << array_.size() << std::endl;
    Point* sder1 = &array_[elem_size_*(left2*MM_+left1)]; // sder1 & eder1 contain values along vmin in Bezier patch.
    Point* eder1 = &array_[elem_size_*(left2*MM_+right1)];
    Point* sder2 = &array_[elem_size_*(right2*MM_+left1)]; // sder2 & eder2 contain values along vmax in Bezier patch.
    Point* eder2 = &array_[elem_size_*(right2*MM_+right1)];
    bezcoef[0] = sder1[0];
    bezcoef[1] = sder1[0]+sder1[1]*scale1;
    bezcoef[2] = eder1[0]-eder1[1]*scale1;
    bezcoef[3] = eder1[0];

    bezcoef[4] = sder1[0] + sder1[2]*scale2;
    bezcoef[5] = bezcoef[1] + sder1[2]*scale2 +
        sder1[3]*scale1*scale2;
    bezcoef[6] = bezcoef[2] + eder1[2]*scale2 - 
        eder1[3]*scale1*scale2;
    bezcoef[7] = eder1[0] + eder1[2]*scale2;

    bezcoef[12] = sder2[0];
    bezcoef[13] = sder2[0]+sder2[1]*scale1;
    bezcoef[14] = eder2[0]-eder2[1]*scale1;
    bezcoef[15] = eder2[0];

    bezcoef[8] = sder2[0] - sder2[2]*scale2;
    bezcoef[9] = bezcoef[13] - sder2[2]*scale2 -
        sder2[3]*scale1*scale2;
    bezcoef[10] = bezcoef[14] - eder2[2]*scale2 + 
        eder2[3]*scale1*scale2;
    bezcoef[11] = eder2[0] - eder2[2]*scale2;



  
    // double scale = (epar - spar)/3.0;
    // bezcoef[0] = array_[elem_size_*left];
    // bezcoef[3] = array_[elem_size_*right];
    // bezcoef[1] = bezcoef[0] + array_[elem_size_*left+1]*scale;
    // bezcoef[2] = bezcoef[3] - array_[elem_size_*right+1]*scale;

    
    return;
}


void HermiteGrid2D::removeGridLines(const std::vector<int>& grid_lines_u,
                                    const std::vector<int>& grid_lines_v)
//---------------------------------------------------------
// PURPOSE: Mark the grid lines as not to be used.
//----------------------------------------------------------
{
    MESSAGE("Function to be removed! Do not call!");
    
    removed_grid_u_ = grid_lines_u;
    removed_grid_v_ = grid_lines_v;
}

    
int HermiteGrid2D::getNoSplitStatus(int ind_u, int ind_v)
{
    int ind = ind_v*MM_ + ind_u;
    return no_split_status_[ind];
}

void HermiteGrid2D::setNoSplitStatus(int ind_u, int ind_v, int no_split_status)
{
    int ind = ind_v*MM_ + ind_u;
    no_split_status_[ind] = no_split_status;

}


int HermiteGrid2D::getPosition(double knot, bool dir_is_u)
//---------------------------------------------------------
// PURPOSE: Find the index into the knot vector of the parameter knot
//----------------------------------------------------------
{
    const vector<double>& knots = (dir_is_u) ? knots_u_ : knots_v_;
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
