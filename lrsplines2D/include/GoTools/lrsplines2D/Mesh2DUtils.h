#ifndef MESH2D_UTILS_H_
#define MESH2D_UTILS_H_

#include "GoTools/lrsplines2D/Mesh2D.h"

namespace Go
{
  namespace Mesh2DUtils
  {
  // Finds the largest [smallest] index of the knotvalue in the mesh 'm' along direction 'd' that is 
  // smaller or equal to [strictly larger than] parameter value 'par'.
  int last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, double par);
  int first_larger_knotvalue_ix(const Mesh2D& m, Direction2D d, double par);

// Determine x_ix and y_ix so that the mesh point (x_ix, y_ix) represents the lower-left (upper-right)
// corner of the smallest rectangle in the mesh that contains the coordinate pair (u, v).
// Rectangles in the mesh are considered to be "closed" downwards and "open" upwards - this determines
// which rectangle is chosen when u or v is situated directly _on_ an internal meshline.  However,
// rectangles that borders the boundary of the mesh domain will always be considered closed along that
// border.  This ensures that all points of the mesh domain are included in exactly one mesh rectangle.
// Return 'true' if such a pair is found (which is <=> to (u, v) being inside the domain of 'm')
bool identify_patch_lower_left (const Mesh2D&m, double u, double v, int& x_ix, int& y_ix);
bool identify_patch_upper_right(const Mesh2D&m, double u, double v, int& x_ix, int& y_ix); 

// Looks for the last line in 'm' with index smaller than 'start_ix' in the specified direction 
// and that has nonzero multiplicity for the segment [other_ix, other-ix+1]
int search_downwards_for_nonzero_multiplicity(const Mesh2D& m, Direction2D d, int start_ix, int other_ix);
int search_upwards_for_nonzero_multiplicity(const Mesh2D& m, Direction2D d, int start_ix, int other_ix);

}; // end namespace Mesh2DUtils

}; // end namespace Go

#endif
