#ifndef MESH2D_ITERATOR_H_
#define MESH2D_ITERATOR_H_

#include <array>



namespace Go
{
class Mesh2D; // forward declaration of Mesh2D 

// =============================================================================
class Mesh2DIterator
// =============================================================================
{
 public:

  // Constructor for empty mesh iterator
  Mesh2DIterator() {};

  // Constructor, taking a mesh and the column/row indices of the lower-left
  // corner of the element in the mesh that the iterator should refer to.  
  // (NB: if there is no element with a lower-left corner at (u_ix, v_ix), the 
  // constructor will search downwards and to the left until it finds a corner
  // for which there is an element).
  Mesh2DIterator(const Mesh2D& m, int u_ix, int v_ix);

  // Increment iterator - move it to refer to the next element in the mesh.
  Mesh2DIterator& operator++();
  
  // Test for equality with another Mesh2DIterator.  Two Mesh2DIterators are considered
  // equal if they refer to the same element in the same mesh.
  bool operator==(const Mesh2DIterator& rhs) const {
    return (m_ == rhs.m_) && (corners_ == rhs.corners_);
  }

  // Test for inequality with another Mesh2DIterator.
  bool operator!=(const Mesh2DIterator& rhs) const {
    return ! (rhs == *this);
  }

  // Swap two Mesh2DIterators
  void swap(Mesh2DIterator& rhs) {
    std::swap(m_, rhs.m_);
    std::swap(corners_, rhs.corners_);
  }

  // Dereference operator - return the defining information of the mesh element
  // that this iterator refers to.  
  // This is returned as an array of four indices.  The first two integers define
  // the indices of the knots of the lower-left corner of the element, whereas the
  // last two integers refer to the upper-right corner.
  const std::array<int, 4>& operator*() const { return corners_;}

 private:
  const Mesh2D* m_;
  std::array<int, 4> corners_;
};


}; // end namespace Go



#endif
