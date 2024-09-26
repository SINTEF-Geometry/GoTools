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

#ifndef _INDEXMESH3DITERATOR_H
#define _INDEXMESH3DITERATOR_H


#include <array>

namespace Go 
{
class Mesh3D; // forward declaration of Mesh3D 

// =============================================================================
class IndexMesh3DIterator
// =============================================================================
{
 public:

  // Constructor: takes a mesh and the lower-left corner of index knot
  // positions (specified by knot indices and multiplicities, first in
  // u-direction then in v-direction) as input, and sets the iterator
  // to point to that element in the index mesh. The default knot
  // indices and multiplicities make the iterator refer to the
  // very first index element, i.e. the lower left element in the
  // index mesh.
  IndexMesh3DIterator(const Mesh3D& m,
		      int kind_u=0, int kmul_u=1,
		      int kind_v=0, int kmul_v=1,
		      int kind_w=0, int kmul_w=1);

  // Increment iterator: moves the iterator to make it refer to the
  // next element in the index mesh.
  IndexMesh3DIterator& operator++();
  
  // Equality operator: tests if two IndexMesh3DIterators are equal. Two
  // IndexMesh3DIterators are considered equal if they refer to the same
  // element in the same index mesh.
  bool operator==(const IndexMesh3DIterator& rhs) const {
    return (m_ == rhs.m_) && (corners_ == rhs.corners_);
  }

  // Inequality operators: negates the equality operator for two
  // IndexMesh3DIterators.
  bool operator!=(const IndexMesh3DIterator& rhs) const {
    return ! (rhs == *this);
  }

  // Dereference operator: returns information of the index mesh
  // element that the iterator refers to. This is returned as an array
  // of 4 IndexKnotPositions, where each IndexKnotPosition is a pair
  // of knot index and multiplicity. The first two IndexKnotPositions
  // refer to the lower-left corner of the element, first in the
  // u-direction and then in the v-direction, while the last two refer
  // to the upper-right corner, first in the u-direction and then in
  // the v-direction. 
  // @@ This dereference method is far from optimal. We could consider
  // to return 4 {index,multiplicity}-structures explicitly, rather
  // than 4x2 integers, to make it more obvious. Also, for examining
  // which basis functions have support on a given index mesh element,
  // we could consider to store multiplicities not only from below
  // (as is done now) but also from above, and hence represent an
  // index knot position as
  // {index,multiplicity_from_below,multiplicity_from_above},
  // resulting in a total of 4x3=12 integers.
  const std::array<int,8>& operator*() const { return corners_;}

 private:

  // Data storage
  const Mesh3D* m_;           // the underlying mesh
  std::array<int,8> corners_; // identifier for an element in the index mesh
};

}; // end namespace Go


#endif // _INDEXMESH3DITERATOR_H

