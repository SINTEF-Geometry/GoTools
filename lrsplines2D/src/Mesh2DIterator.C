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
#include "GoTools/lrsplines2D/Mesh2DIterator.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"

using namespace Go;
using std::vector;
using std::fill;

namespace {

inline bool is_corner(const Mesh2D&m, int u, int v) {
  return ((m.nu(XFIXED, u, v, v+1) != 0) &&
	  (m.nu(YFIXED, v, u, u+1) != 0));
}

}; // end anonymous namespace


//==============================================================================
Mesh2DIterator::Mesh2DIterator(const Mesh2D& m, int u_ix, int v_ix)
//==============================================================================
  : m_(&m)
{
  if (u_ix < 0 || v_ix < 0) 
    THROW("Mesh2DIterator constructor called with negative indices.");
  if ((u_ix >= m_->numDistinctKnots(XFIXED) - 1) ||
      (v_ix >= m_->numDistinctKnots(YFIXED) - 1)) {

    fill(corners_.begin(), corners_.end(), -1);  // this iterator represents past-end of mesh

  } else {
    // locating lower-left corner
    corners_[0] = Mesh2DUtils::search_downwards_for_nonzero_multiplicity(*m_, XFIXED, u_ix, v_ix);
    corners_[1] = Mesh2DUtils::search_downwards_for_nonzero_multiplicity(*m_, YFIXED, v_ix, u_ix);
    // locating upper-right corner
    ++u_ix; // we need to increase by one to ensure that the search_upwards routine below
    ++v_ix; // doesn't stop before it has started...
    corners_[2] = Mesh2DUtils::search_upwards_for_nonzero_multiplicity(*m_, XFIXED, u_ix, v_ix);
    corners_[3] = Mesh2DUtils::search_upwards_for_nonzero_multiplicity(*m_, YFIXED, v_ix, u_ix);
  }
}

//==============================================================================
Mesh2DIterator& Mesh2DIterator::operator++()
//==============================================================================
{
  // looking for new corner
  const int start_v = corners_[1];
  const int end_v   = m_->numDistinctKnots(YFIXED) - 1;
  const int start_u = corners_[2]; // We start the search at the end of the u index.
  const int end_u   = m_->numDistinctKnots(XFIXED) - 1;
  for (int v = start_v; v != end_v; ++v) {
    for (int u = (v == start_v ? start_u : 0); u != end_u; ++u) {
      if (is_corner(*m_, u, v)) {
	Mesh2DIterator tmp(*m_, u, v);
	this->swap(tmp);
	return *this;
      }
    }
  }
  fill(corners_.begin(), corners_.end(), -1); // we have arrived past end of mesh
  return *this;
}


