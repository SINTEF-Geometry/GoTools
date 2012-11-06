#ifndef INDEX_MESH2D_ITERATOR_H_
#define INDEX_MESH2D_ITERATOR_H_

#include <array>

namespace Go 
{
class Mesh2D; // forward declaration of Mesh2D 

// =============================================================================
class IndexMesh2DIterator
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
  IndexMesh2DIterator(const Mesh2D& m, int kind_u=0, int kmul_u=1, int kind_v=0, int kmul_v=1);

  // Increment iterator: moves the iterator to make it refer to the
  // next element in the index mesh.
  IndexMesh2DIterator& operator++();
  
  // Equality operator: tests if two IndexMesh2DIterators are equal. Two
  // IndexMesh2DIterators are considered equal if they refer to the same
  // element in the same index mesh.
  bool operator==(const IndexMesh2DIterator& rhs) const {
    return (m_ == rhs.m_) && (corners_ == rhs.corners_);
  }

  // Inequality operators: negates the equality operator for two
  // IndexMesh2DIterators.
  bool operator!=(const IndexMesh2DIterator& rhs) const {
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
  const Mesh2D* m_;           // the underlying mesh
  std::array<int,8> corners_; // identifier for an element in the index mesh
};

}; // end namespace Go

#endif
