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

#include "GoTools/creators/HermiteGrid1D.h"
#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"

using namespace std;


namespace Go
{

HermiteGrid1D::HermiteGrid1D(const EvalCurve& crv, double t1, double t2)
  : dim_(crv.dim()), MM_(2), elem_size_(2), index_(0)
{
  // Copy knot vector

  knots_.reserve(MM_);
  knots_.push_back(t1);
  knots_.push_back(t2);

  // Calculate the curve values at the parameter grid

  array_.reserve(elem_size_*MM_);
  Point derive[2];
  int i, j;
  for (i=0; i<MM_; i++)
  {
    crv.eval(knots_[i],1,derive);
    for (j=0; j<elem_size_; j++)
      array_.push_back(derive[j]);
  }
}

HermiteGrid1D::HermiteGrid1D(const EvalCurve& crv, double param[], int n)
  : dim_(crv.dim()), MM_(n), elem_size_(2), index_(n/2)
//--------------------------------------------------------
//  Constructor
//
// INPUT:
//      param  - Array of strictly increasing parameters in the parameter
//               interval of "crv". First and last point defines the interval.
//      n      - Number of elements in "param"
//      crv    - Curve to sample
//--------------------------------------------------------
{
  // Check that knots are strictly increasing
  int i, j;
  for (i=1; i<n; i++)
    if (param[i] <= param[i-1])
      THROW("Input grid illegal");
      

  // Calculate the curve values at the parameter grid

  array_.reserve(elem_size_*n);

  Point derive[2];
  for (i=0; i<n; i++)
  {
    crv.eval(param[i],1,derive);
    for (j=0; j<elem_size_; j++)
      array_.push_back(derive[j]);
  }

  // Copy knot vector

  knots_.reserve(n);
  for (i=0; i<n; i++)
    knots_.push_back(param[i]);
}

HermiteGrid1D::~HermiteGrid1D()
//--------------------------------------------------------
//  Destructor
//--------------------------------------------------------
{}


int HermiteGrid1D::addKnot(const EvalCurve& crv, double knot)
//--------------------------------------------------------------------
// PURPOSE: Insert a new knot in the knotvector . Also the value and tangent 
//          of the curve at this knot is added to the Hermite grid
//
// INPUT:
//      crv	- Curve to evaluate
//      param	- New knot
// OUTPUT:
//      addKnot() - The index of the new knot in the (sorted) knot vector
//                  after insertion. The first index for this knotvector is 0.
//--------------------------------------------------------------------
{
  // Evaluate the curve at the new knot and insert the values into
  // the Hermite grid.

  Point derive[2];
  crv.eval(knot,1,derive);

  // Insert the new knot into the knot vector

  index_ = getPosition(knot);
  knots_.insert(knots_.begin()+index_+1,1,knot);
  array_.insert(array_.begin()+elem_size_*(index_+1), 1, derive[0]);
  array_.insert(array_.begin()+elem_size_*(index_+1)+1, 1, derive[1]);

  MM_++;

  return index_;
}

int HermiteGrid1D::getPosition(double knot)
//---------------------------------------------------------
// PURPOSE: Find the index into the knot vector of the parameter knot
//----------------------------------------------------------
{
    int id1=0, id2=(int)knots_.size()-1;
  if ((index_ == id2 && knot >= knots_[index_]) ||
      (knots_[index_] <= knot && knot < knots_[index_]))
    return index_;

  // Search
  index_ = id1; // We must reset index_ before we start searching.
  while (id2 > id1+1)
    {
      if (knots_[index_] > knot)
	id2 = index_;
      if (knots_[index_] <= knot)
	id1 = index_;
      index_ = (id1 + id2)/2;
    }
  
  return index_;
}


void HermiteGrid1D::getSegment(int left, int right,
				 double& spar, double& epar,
				 Point bezcoef[4])
//--------------------------------------------------------------------
// PURPOSE: Calculate Bezier coefficients of cubic curve interpolating the
//          Hermite values at grid nodes with indeces "left" and "right"
//
// INPUT:
//     left - indicating grid node for start of curve segment
//     right - indicating grid node for end of curve segment
// OUTPUT:
//      spar    - start parameter of segment
//      epar    - end parameter of segment
//      bezcoef - array of cubic Bezier coefficients
//--------------------------------------------------------------------
{
  spar = knots_[left];
  epar = knots_[right];

  double scale = (epar - spar)/3.0;
  bezcoef[0] = array_[elem_size_*left];
  bezcoef[3] = array_[elem_size_*right];
  bezcoef[1] = bezcoef[0] + array_[elem_size_*left+1]*scale;
  bezcoef[2] = bezcoef[3] - array_[elem_size_*right+1]*scale;

  return;
}

}
