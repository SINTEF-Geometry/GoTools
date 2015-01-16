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

#include "GoTools/lrsplines2D/Mesh2DUtils.h"

#include <map>
#include <array>
#include <stdexcept> // @@ for debug purposes only 

using namespace std;

namespace Go {


// =============================================================================
  bool Mesh2DUtils::identify_patch_lower_left(const Mesh2D&m, double u, 
					      double v, int& x_ix, int& y_ix)
// =============================================================================
{
  double tol = 1.0e-9;

  x_ix = last_nonlarger_knotvalue_ix(m, XFIXED, u);
  y_ix = last_nonlarger_knotvalue_ix(m, YFIXED, v);

  // adjustment of index if positioned _exactly_ at upper bound of grid
  if (x_ix == m.numDistinctKnots(XFIXED) - 1 && fabs(u-m.maxParam(XFIXED)) < tol) 
    --x_ix;
  if (y_ix == m.numDistinctKnots(YFIXED) - 1 && fabs(v-m.maxParam(YFIXED)) < tol)
    --y_ix;
  
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
  double tol = 1.0e-9;

  x_ix = first_larger_knotvalue_ix(m, XFIXED, u);
  y_ix = first_larger_knotvalue_ix(m, YFIXED, v);
  
  // adjustment of index if positioned _exactly_ at upper bound of grid
  if (x_ix == m.numDistinctKnots(XFIXED) - 1 && fabs(u-m.maxParam(XFIXED)) < tol) 
    --x_ix;
  if (y_ix == m.numDistinctKnots(YFIXED) - 1 && fabs(v-m.maxParam(YFIXED)) < tol)
    --y_ix;

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

// =============================================================================
int Mesh2DUtils::search_downwards_for_nonzero_multiplicity(const Mesh2D& m, 
							   Direction2D d, 
							   int start_ix, 
							   int ix1, int ix2)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix, ix1, ix2) == 0; --ix);
  return ix;
}

//==============================================================================
int Mesh2DUtils::search_upwards_for_nonzero_multiplicity(const Mesh2D& m, 
							 Direction2D d, 
							 int start_ix, 
							 int ix1, int ix2)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix, ix1, ix2) == 0; ++ix);
  return ix;  // @@ not yet tested!
}





// =============================================================================
int Mesh2DUtils::last_nonlarger_knotvalue_ix(const Mesh2D&m, Direction2D d, 
					     double par)
// =============================================================================
{
  const double* a = m.knotsBegin(d);
  const double* b = m.knotsEnd(d);
  
  double tol = 1.0e-9;
  // if (par < a[0] && par >= a[0]-tol)
  //   par = a[0];
  // if (par > b[-1] && par <= b[-1]+tol)
  //   par = b[-1];

  // searching for last nonlarger knotvalue using bisection
  for (int diff = (b-a)/2; diff != 0; diff = (b-a)/2) {
    if (par < a[0]+tol && par >= a[0]-tol)
      par = a[0];
    if (par > b[-1]-tol && par <= b[-1]+tol)
      par = b[-1];
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
  int Mesh2DUtils::first_larger_knotvalue_ix(const Mesh2D& m, Direction2D d, 
					     double par)
// =============================================================================
{
  int ix = last_nonlarger_knotvalue_ix(m, d, par);
  if (ix < m.numDistinctKnots(d)-1)
    ix++;
  return  ix; 
}

  
// =============================================================================
}; // end namespace Go
// =============================================================================




