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

#ifndef _MESH3DITERATOR_H
#define _MESH3DITERATOR_H


#include <array>



namespace Go
{
class Mesh3D; // forward declaration of Mesh3D 

// =============================================================================
class Mesh3DIterator
// =============================================================================
{
 public:

  // Constructor for empty mesh iterator
  Mesh3DIterator() {}

  // Constructor, taking a mesh and the column/row indices of the lower-left
  // corner of the element in the mesh that the iterator should refer to.  
  // (NB: if there is no element with a lower-left corner at (u_ix, v_ix), the 
  // constructor will search downwards and to the left until it finds a corner
  // for which there is an element).
    Mesh3DIterator(const Mesh3D& m, int u_ix, int v_ix, int w_ix);

  // Increment iterator - move it to refer to the next element in the mesh.
  Mesh3DIterator& operator++();
  
  // Test for equality with another Mesh3DIterator.  Two Mesh3DIterators are considered
  // equal if they refer to the same element in the same mesh.
  bool operator==(const Mesh3DIterator& rhs) const {
    return (m_ == rhs.m_) && (corners_ == rhs.corners_);
  }

  // Test for inequality with another Mesh3DIterator.
  bool operator!=(const Mesh3DIterator& rhs) const {
    return ! (rhs == *this);
  }

  // Swap two Mesh3DIterators
  void swap(Mesh3DIterator& rhs) {
    std::swap(m_, rhs.m_);
    std::swap(corners_, rhs.corners_);
  }

  // Dereference operator - return the defining information of the mesh element
  // that this iterator refers to.  
  // This is returned as an array of six indices.  The first three integers define
  // the indices of the knots of the lower corner of the element, whereas the
  // last three integers refer to the upper corner.
  const std::array<int, 6>& operator*() const { return corners_;}

 private:
  const Mesh3D* m_;
  std::array<int, 6> corners_; // start_u, start_v, start_w, end_u, end_v, end_w
};


} // end namespace Go





#endif // _MESH3DITERATOR_H

