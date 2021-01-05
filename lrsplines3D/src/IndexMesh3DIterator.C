//===========================================================================
//                                                                           
// File: IndexMesh3DIterator.C                                               
//                                                                           
// Created: Mon Mar 11 15:07:21 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include <algorithm>
#include <stdexcept>
#include <utility>
#include <array>
#include "GoTools/lrsplines3D/IndexMesh3DIterator.h"
#include "GoTools/lrsplines3D/Mesh3DUtils.h"
#include "GoTools/lrsplines3D/Direction3D.h"


namespace Go 
{
class Mesh3D; // forward declaration of Mesh3D 


//==============================================================================
IndexMesh3DIterator::IndexMesh3DIterator(const Mesh3D& m,
					 int kind_u, int kmul_u,
					 int kind_v, int kmul_v,
					 int kind_w, int kmul_w)
//==============================================================================
   : m_(&m)
{
    MESSAGE("Not implemented yet!");

#if 0
  // If the input knot indices and multiplicities are all -1's, the
  // IndexMesh3DIterator should point to a past-end-element of the index
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
      throw runtime_error("Constructor for IndexMesh3DIterator called with an invalid lower left corner (knot indices and multiplicities).");
    }
  }    
#endif
}


//==============================================================================
IndexMesh3DIterator& IndexMesh3DIterator::operator++()
//==============================================================================
{
    MESSAGE("Not implemented yet!");

#if 0
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
#endif

  return *this;
}


}; // end namespace Go
