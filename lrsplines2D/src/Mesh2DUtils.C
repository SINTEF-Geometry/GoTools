#include "GoTools/lrsplines2D/Mesh2DUtils.h"

#include <map>
#include <array>
#include <stdexcept> // @@ for debug purposes only 

using namespace std;

namespace Go {
namespace { // anonymous namespace

  // Helper function for 'identify patch_lower_left' ['identify_patch_upper_right']
  // Finds the largest [smallest] index of the knotvalue in the mesh 'm' along direction 'd' that is 
  // smaller or equal to [strictly larger than] parameter value 'par'.
  int last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, double par);
  int first_larger_knotvalue_ix(const Mesh2D& m, Direction2D d, double par);

} // end anonymous namespace


// =============================================================================
  bool Mesh2DUtils::identify_patch_lower_left(const Mesh2D&m, double u, 
					      double v, int& x_ix, int& y_ix)
// =============================================================================
{
  x_ix = last_nonlarger_knotvalue_ix(m, XFIXED, u);
  y_ix = last_nonlarger_knotvalue_ix(m, YFIXED, v);

  // adjustment of index if positioned _exactly_ at upper bound of grid
  if (x_ix == m.numDistinctKnots(XFIXED) - 1 && u == m.maxParam(XFIXED)) --x_ix;
  if (y_ix == m.numDistinctKnots(YFIXED) - 1 && v == m.maxParam(YFIXED)) --y_ix;
  
  // checking if a valid corner was found
  if (x_ix < 0 || x_ix >= m.numDistinctKnots(XFIXED) - 1) return false; // u outside domain
  if (y_ix < 0 || y_ix >= m.numDistinctKnots(YFIXED) - 1) return false; // v outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  x_ix = search_downwards_for_nonzero_multiplicity(m, XFIXED, x_ix, y_ix);
  y_ix = search_downwards_for_nonzero_multiplicity(m, YFIXED, y_ix, x_ix);

  return true;
}

// =============================================================================
  bool Mesh2DUtils::identify_patch_upper_right(const Mesh2D&m, double u, 
					       double v, int& x_ix, int& y_ix)
// =============================================================================
{
  x_ix = first_larger_knotvalue_ix(m, XFIXED, u);
  y_ix = first_larger_knotvalue_ix(m, YFIXED, v);
  
  // adjustment of index if positioned _exactly_ at upper bound of grid
  if (x_ix == m.numDistinctKnots(XFIXED) && u == m.maxParam(XFIXED)) --x_ix;
  if (y_ix == m.numDistinctKnots(YFIXED) && v == m.maxParam(YFIXED)) --y_ix;

  // checking if a valid corner was found
  if (x_ix == 0 || x_ix >= m.numDistinctKnots(XFIXED)) return false; // u outside domain
  if (y_ix == 0 || y_ix >= m.numDistinctKnots(YFIXED)) return false; // v outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  x_ix = search_upwards_for_nonzero_multiplicity(m, XFIXED, x_ix, y_ix);
  y_ix = search_upwards_for_nonzero_multiplicity(m, YFIXED, y_ix, x_ix);

  return true;
}

// =============================================================================
int Mesh2DUtils::search_downwards_for_nonzero_multiplicity(const Mesh2D& m, 
							   Direction2D d, 
							   int start_ix, int other_ix)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix, other_ix, other_ix + 1) == 0; --ix);
  return ix;
}

//==============================================================================
int Mesh2DUtils::search_upwards_for_nonzero_multiplicity(const Mesh2D& m, 
							 Direction2D d, 
							 int start_ix, int other_ix)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix, other_ix - 1, other_ix) == 0; ++ix);
  return ix;  // @@ not yet tested!
}




namespace { // anonymous namespace 


// =============================================================================
int last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, double par)
// =============================================================================
{
  const double* a = m.knotsBegin(d);
  const double* b = m.knotsEnd(d);

  // searching for last nonlarger knotvalue using bisection
  for (int diff = (b-a)/2; diff != 0; diff = (b-a)/2) {
    const double* mid = a + diff;
    ( (*mid > par) ? b : a) = mid;
  }

  return (a - m.knotsBegin(d)); // if index becomes negative, it signalizes that 'par' 
                                // is smaller than even the first knot value
}



// // =============================================================================
// int last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, double par)
// // =============================================================================
// {
//   const int IMAX = m.numDistinctKnots(d);
//   int ix;
//   for (ix = 0; ix != IMAX && m.kval(d, ix) <= par; ++ix); // the first one that is larger

  
//   --ix; // by decrementing index by one, we get to the last nonlarger knot value.  If index becomes
//         // negative, it signalizes that 'par' is smaller than even the first knot value.
  

//   return ix;
// }

// =============================================================================
  int first_larger_knotvalue_ix(const Mesh2D& m, Direction2D d, double par)
// =============================================================================
{
  return last_nonlarger_knotvalue_ix(m, d, par) + 1; // @@ untested, but should be ok???
}

}; // end anonymous namespace
  
// =============================================================================
}; // end namespace Go
// =============================================================================




