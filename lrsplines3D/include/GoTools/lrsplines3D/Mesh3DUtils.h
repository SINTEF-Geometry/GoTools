//===========================================================================
//                                                                           
// File: Mesh3DUtils.h                                                       
//                                                                           
// Created: Mon Feb 25 11:08:10 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _MESH3DUTILS_H
#define _MESH3DUTILS_H


#include "GoTools/lrsplines3D/Mesh3D.h"

namespace Go
{
  /// Utility functions for operating on an LR Mesh (Mesh3D)
  namespace Mesh3DUtils
  {
  // Finds the largest [smallest] index of the knotvalue in the mesh 'm' along direction 'd' that is 
  // smaller or equal to [strictly larger than] parameter value 'par'.
  int last_nonlarger_knotvalue_ix(const Mesh3D&m, Direction3D d, double par);
  int first_larger_knotvalue_ix(const Mesh3D& m, Direction3D d, double par);

// Determine x_ix and y_ix so that the mesh point (x_ix, y_ix) represents the lower-left (upper-right)
// corner of the smallest rectangle in the mesh that contains the coordinate pair (u, v).
// Rectangles in the mesh are considered to be "closed" downwards and "open" upwards - this determines
// which rectangle is chosen when u or v is situated directly _on_ an internal meshline.  However,
// rectangles that borders the boundary of the mesh domain will always be considered closed along that
// border.  This ensures that all points of the mesh domain are included in exactly one mesh rectangle.
// Return 'true' if such a pair is found (which is <=> to (u, v) being inside the domain of 'm')
bool identify_patch_lower_left (const Mesh3D&m,
				double u, double v, double w,
				int& x_ix, int& y_ix, int& z_ix);
bool identify_patch_upper_right(const Mesh3D&m,
				double u, double v, double w,
				int& x_ix, int& y_ix, int& z_ix); 

// Looks for the last line in 'm' with index smaller than 'start_ix' in the specified direction 
// and that has nonzero multiplicity for the segment [other_ix, other-ix+1]
int search_downwards_for_nonzero_multiplicity(const Mesh3D& m, Direction3D d, 
					      int start_ix, int other1_ix, int other2_ix);
int search_upwards_for_nonzero_multiplicity(const Mesh3D& m, Direction3D d, 
					    int start_ix, int other1_ix, int other2_ix);
int search_downwards_for_nonzero_multiplicity(const Mesh3D& m, Direction3D d, 
					      int start_ix,
					      int other1_ix1, int other1_ix2,
					      int other2_ix1, int other2_ix2);
int search_upwards_for_nonzero_multiplicity(const Mesh3D& m, Direction3D d, 
					    int start_ix,
					    int other1_ix1, int other1_ix2,
					    int other2_ix1, int other2_ix2);

 void sectionKnots(const Mesh3D& mesh, Direction3D d, int ix,
		   std::vector<std::vector<int> >& knots1,
		   std::vector<std::vector<int> >& knots2);
}; // end namespace Mesh3DUtils

}; // end namespace Go


#endif // _MESH3DUTILS_H

