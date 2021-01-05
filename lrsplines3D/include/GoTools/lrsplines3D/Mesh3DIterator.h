//===========================================================================
//                                                                           
// File: Mesh3DIterator.h                                                    
//                                                                           
// Created: Mon Feb 25 11:08:01 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

