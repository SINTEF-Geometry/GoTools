#include <algorithm>
#include <stdexcept>
#include <utility>
#include <array>
#include "GoTools/lrsplines2D/IndexMesh2DIterator.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/lrsplines2D/Direction2D.h"

using namespace std;

namespace Go {
namespace { // anonymous namespace

// =============================================================================
// Below is: Header for IndexKnotPosition
// =============================================================================
class IndexKnotPosition
{ 
 public:

  // Constructor: empty.
  IndexKnotPosition();

  // Constructor: takes a mesh, a parametric direction, a knot index
  // and a knot multiplicity.
  IndexKnotPosition(const Mesh2D& m, const Direction2D dir, int kind=0, int kmul=1);

  // Operator increment: move the iterator to make it refer to the
  // next position in the index knot.
  IndexKnotPosition& operator++();
  
  // Dereference operator: returns the defining information of this
  // index knot position. The information is returned as an array of
  // two integers. The first is the index of knot, referring to the
  // global unique knot vector, and the second is the multiplicity of
  // the knot.
  const std::array<int, 2>& operator*() const { return knot_ind_mul_;}

  // Method: returns true if current IndexKnotsPosition is a valid
  // index knot position, i.e. if the knot index and multiplicity are
  // within their limits, and returns false otherwise.
  bool isValid() const;

  // Method: returns true if current IndexKnotsPosition is the last
  // index knot position.
  bool isLast() const;

 private:

  // Data
  std::vector<int> knots_mult_;    // highest multiplicity of all unique knots in the specified direction
  std::array<int,2> knot_ind_mul_; // identifier for index knot position: {knot index, knot multiplicity}
};

// =============================================================================
// Below is: Implementation for IndexKnotPosition
// =============================================================================

//==============================================================================
IndexKnotPosition::IndexKnotPosition(const Mesh2D& m, const Direction2D dir, int kind, int kmul)
//==============================================================================
  : knots_mult_(m.numDistinctKnots(dir))
{
  const int nknots = knots_mult_.size();
  for (int it=0; it<nknots; it++)                   // loop over unique knots
    knots_mult_[it] = m.largestMultInLine(dir, it); // store maximum knot multiplicity along mesh line
  knot_ind_mul_[0] = kind; // knot index
  knot_ind_mul_[1] = kmul; // knot multiplicity
  if (!(isValid())) // make sure this is valid IndexKnotPosition
    throw runtime_error("Constructor for IndexKnotPosition called with invalid knot index and multiplicity.");
}

//==============================================================================
IndexKnotPosition& IndexKnotPosition::operator++()
//==============================================================================
{
  int& kind = knot_ind_mul_[0]; // knot index
  int& kmul = knot_ind_mul_[1]; // knot multiplicity
  kmul++; // Assume the index is the same but with multiplicity increased by one.
  if (kmul > knots_mult_[kind]) { // If this was wrong then ...
    kind++;                       // ... increment the index ...
    kmul = 1;                     // ... and set the multiplicity to one.
  }
  return *this;
}

//==============================================================================
bool IndexKnotPosition::isValid() const
//==============================================================================
{
  const int kind = knot_ind_mul_[0]; // knot index
  const int kmul = knot_ind_mul_[1]; // knot multiplicity
  return ( (kind >= 0 && kind < knots_mult_.size()) &&
           (kmul >  0 && kmul <= knots_mult_[kind]) );
}

//==============================================================================
bool IndexKnotPosition::isLast() const
//==============================================================================
{
  const int kind = knot_ind_mul_[0]; // knot index
  const int kmul = knot_ind_mul_[1]; // knot multiplicity
  return ( ((kind==(knots_mult_.size()-1)) && (kmul==knots_mult_[kind])) );
}

// =============================================================================
// Below is: anonymous methods for IndexMesh2DIterator
// =============================================================================

// Returns the nearest valid index knot position, i.e. non-vanishing
// knot index and multiplicity in the mesh, from a given knot index
// position in one direction, searching upwards from that position in
// that direction along a specified knot index in the other direction.
IndexKnotPosition search_upwards_for_index_knot_position(const Mesh2D& m, Direction2D dir_this, IndexKnotPosition ikp_this, IndexKnotPosition ikp_other)
{
  // Initialize information for the OTHER knot index position
  Direction2D dir_other = (dir_this==YFIXED) ? XFIXED : YFIXED;         // Determine the direction and knot index ...
  int i_other = (*ikp_other)[0];                                      // ... of the OTHER knot index position and ...
  i_other += (i_other == (m.numDistinctKnots(dir_other)-1)) ? -1 : 0; // ... adjust if it is on the upper boundary.
  // Loop forward over THIS index knot position, and start by incrementing it by one
  for (++ikp_this; ikp_this.isValid(); ++ikp_this) {
    int i_this = (*ikp_this)[0];      // THIS knot index
    int m_this = (*ikp_this)[1];      // THIS knot multiplicity
    if ( m_this <= m.nu(dir_this,i_this,i_other,i_other+1) )
      break; // 1) THIS is a valid index knot position
  }
  // Reaching here, THIS index knot position is either 1) the next
  // valid index knot position or 2) an invalid index knot position
  // (i.e. the end).
  return ikp_this;
}

// Returns true if the specified knot indices and knot multiplicities
// together represent a valid corner in the index mesh.
bool is_index_mesh_corner(const Mesh2D& m, IndexKnotPosition ikp_u, IndexKnotPosition ikp_v)
{
  int kind_u = (*ikp_u)[0];
  int kmul_u = (*ikp_u)[1];
  int kind_v = (*ikp_v)[0];
  int kmul_v = (*ikp_v)[1];
  int kind_u0 = kind_u - ((kind_u==(m.numDistinctKnots(XFIXED)-1)) ? 1 : 0);
  int kind_v0 = kind_v - ((kind_v==(m.numDistinctKnots(YFIXED)-1)) ? 1 : 0);
  return ( (kmul_u <= m.nu(XFIXED,kind_u,kind_v0,kind_v0+1)) && 
           (kmul_v <= m.nu(YFIXED,kind_v,kind_u0,kind_u0+1)) );
}

}; // end anonymous namespace

// =============================================================================
// Below is: Implementation for IndexMesh2DIterator
// =============================================================================


//==============================================================================
IndexMesh2DIterator::IndexMesh2DIterator(const Mesh2D& m, int kind_u, int kmul_u, int kind_v, int kmul_v)
//==============================================================================
   : m_(&m)
{
  // If the input knot indices and multiplicities are all -1's, the
  // IndexMesh2DIterator should point to a past-end-element of the index
  // mesh, which is represented by filling the corner information with
  // -1's.
  if ( (kind_u==-1 && kmul_u==-1) && (kind_v==-1 && kmul_v==-1) ) {
    fill(corners_.begin(), corners_.end(), -1); 
    return;
  }
  // Else we must examine the specified knot indices and
  // multiplicities for the lower left corner, and set the index mesh
  // element pointer appropriately.
  IndexKnotPosition ikp_u_ll = IndexKnotPosition (*m_, XFIXED, kind_u, kmul_u);
  IndexKnotPosition ikp_v_ll = IndexKnotPosition (*m_, YFIXED, kind_v, kmul_v);
  if ( (ikp_u_ll.isLast() || ikp_v_ll.isLast()) ) {// this is a past-end element of index mesh
    fill(corners_.begin(), corners_.end(), -1);
  } else {
    if (is_index_mesh_corner(*m_, ikp_u_ll, ikp_v_ll)) {// this is a valid index mesh element
      // locate upper-right corner knot indices and knot multiplicities
      IndexKnotPosition ikp_u_ur = search_upwards_for_index_knot_position(*m_, XFIXED, ikp_u_ll, ikp_v_ll);
      IndexKnotPosition ikp_v_ur = search_upwards_for_index_knot_position(*m_, YFIXED, ikp_v_ll, ikp_u_ll);
      // store all corner knot indices and knot multiplicities
      corners_[0] = (*ikp_u_ll)[0];
      corners_[1] = (*ikp_u_ll)[1];
      corners_[2] = (*ikp_v_ll)[0];
      corners_[3] = (*ikp_v_ll)[1];
      corners_[4] = (*ikp_u_ur)[0];
      corners_[5] = (*ikp_u_ur)[1];
      corners_[6] = (*ikp_v_ur)[0];
      corners_[7] = (*ikp_v_ur)[1];
    } else {// this is an invalid index mesh element
      throw runtime_error("Constructor for IndexMesh2DIterator called with an invalid lower left corner (knot indices and multiplicities).");
    }
  }    
}

//==============================================================================
IndexMesh2DIterator& IndexMesh2DIterator::operator++()
//==============================================================================
{
  const int kind_u = corners_[0]; // knot index        in u-direction
  const int kmul_u = corners_[1]; // knot multiplicity in u-direction
  const int kind_v = corners_[2]; // knot index        in v-direction
  const int kmul_v = corners_[3]; // knot multiplicity in v-direction
  bool is_first_v_loop = true;
  // Loop over knot index positions in v-direction
  for ( IndexKnotPosition ikp_v_ll = IndexKnotPosition (*m_, YFIXED, kind_v, kmul_v); 
        (ikp_v_ll.isValid() && !ikp_v_ll.isLast());
        ++ikp_v_ll, is_first_v_loop = false ) {
    // Loop over knot index position in u-direction
    for ( IndexKnotPosition ikp_u_ll = (is_first_v_loop)
            ? ++(IndexKnotPosition (*m_, XFIXED, kind_u, kmul_u))  // first iteration in v-direction
            :   (IndexKnotPosition (*m_, XFIXED) );                // second or higher iteration in v-direction
	 (ikp_u_ll.isValid() && !ikp_u_ll.isLast());
         ++ikp_u_ll) {
      // Check if the current knot indices and multiplicities
      // represent the lower-left corner of a valid element in index
      // mesh.
      if (is_index_mesh_corner(*m_, ikp_u_ll, ikp_v_ll)) {
        // if so then locate the upper-right corner, store the
        // corner information and return it.
        IndexKnotPosition ikp_u_ur = search_upwards_for_index_knot_position(*m_, XFIXED, ikp_u_ll, ikp_v_ll);
        IndexKnotPosition ikp_v_ur = search_upwards_for_index_knot_position(*m_, YFIXED, ikp_v_ll, ikp_u_ll);
        corners_[0] = (*ikp_u_ll)[0];
        corners_[1] = (*ikp_u_ll)[1];
        corners_[2] = (*ikp_v_ll)[0];
        corners_[3] = (*ikp_v_ll)[1];
        corners_[4] = (*ikp_u_ur)[0];
        corners_[5] = (*ikp_u_ur)[1];
        corners_[6] = (*ikp_v_ur)[0];
        corners_[7] = (*ikp_v_ur)[1];
        return *this;
      } // end if
    } // end u-loop
  } // end v-loop
  // If we made it here, we have reached past-end of the index mesh,
  // and we represent this by filling the corner information with -1's
  // and return it.
  fill(corners_.begin(), corners_.end(), -1);
  return *this;
}

}; // end namespace Go
