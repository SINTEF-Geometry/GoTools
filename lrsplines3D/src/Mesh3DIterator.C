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

#include <algorithm>
#include <stdexcept>
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/Mesh3DIterator.h"
#include "GoTools/lrsplines3D/Mesh3DUtils.h"

using namespace Go;
using std::vector;
using std::fill;

namespace {

  inline bool is_corner(const Mesh3D&m, int u, int v, int w) {
    return ((m.nu(XDIR, u, v, v+1, w, w+1) != 0) &&
	    (m.nu(YDIR, v, w, w+1, u, u+1) != 0) &&
	    (m.nu(ZDIR, w, u, u+1, v, v+1) != 0));
}

}; // end anonymous namespace


//==============================================================================
Mesh3DIterator::Mesh3DIterator(const Mesh3D& m, int u_ix, int v_ix, int w_ix)
//==============================================================================
  : m_(&m)
{
  if (u_ix < 0 || v_ix < 0 || w_ix < 0) 
    THROW("Mesh3DIterator constructor called with negative indices.");
  if ((u_ix >= m_->numDistinctKnots(XDIR) - 1) ||
      (v_ix >= m_->numDistinctKnots(YDIR) - 1) ||
      (w_ix >= m_->numDistinctKnots(ZDIR) - 1)) {

    fill(corners_.begin(), corners_.end(), -1);  // this iterator represents past-end of mesh

  } else {
    // locating lower-left corner
    corners_[0] = Mesh3DUtils::search_downwards_for_nonzero_multiplicity(*m_, XDIR, u_ix, v_ix, w_ix);
    corners_[1] = Mesh3DUtils::search_downwards_for_nonzero_multiplicity(*m_, YDIR, v_ix, w_ix, u_ix);
    corners_[2] = Mesh3DUtils::search_downwards_for_nonzero_multiplicity(*m_, ZDIR, w_ix, u_ix, v_ix);
    // locating upper-right corner
    ++u_ix; // we need to increase by one to ensure that the search_upwards routine below
    ++v_ix; // doesn't stop before it has started...
    ++w_ix; // doesn't stop before it has started...
    corners_[3] = Mesh3DUtils::search_upwards_for_nonzero_multiplicity(*m_, XDIR, u_ix, v_ix, w_ix);
    corners_[4] = Mesh3DUtils::search_upwards_for_nonzero_multiplicity(*m_, YDIR, v_ix, w_ix, u_ix);
    corners_[5] = Mesh3DUtils::search_upwards_for_nonzero_multiplicity(*m_, ZDIR, w_ix, u_ix, v_ix);
  }
}

//==============================================================================
Mesh3DIterator& Mesh3DIterator::operator++()
//==============================================================================
{
  // The boxes are ordered by ll[2] first, then by ll[1], finally by ll[0].
  const int start_w = corners_[2];
  const int end_w = m_->numDistinctKnots(ZDIR) - 1;
  const int start_v = corners_[1];//[4];//1];
  const int end_v = m_->numDistinctKnots(YDIR) - 1;
  const int start_u = corners_[3];
  const int end_u   = m_->numDistinctKnots(XDIR) - 1;
  for (int w = start_w; w != end_w; ++w)
    {
//      for (int v = start_v; v != end_v; ++v)
      for (int v = (w == start_w ? start_v : 0); v != end_v; ++v)
	{
	  for (int u = ((w == start_w && v == start_v) ? start_u : 0); u != end_u; ++u)
	    {
	      if (is_corner(*m_, u, v, w))
		{
		  Mesh3DIterator tmp(*m_, u, v, w);
		  this->swap(tmp);
		  return *this;
		}
	    }
	}
    }

  fill(corners_.begin(), corners_.end(), -1); // we have arrived past end of mesh
  return *this;
}


