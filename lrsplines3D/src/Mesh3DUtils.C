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

#include "GoTools/lrsplines3D/Mesh3DUtils.h"
#include "GoTools/lrsplines3D/LRBSpline3DUtils.h"

#include <map>
#include <array>
#include <stdexcept> // @@ for debug purposes only 

using namespace std;

namespace Go {


// =============================================================================
  bool Mesh3DUtils::identify_patch_lower_left(const Mesh3D&m,
					      double u, double v, double w,
					      int& x_ix, int& y_ix, int& z_ix)
// =============================================================================
{
  x_ix = last_nonlarger_knotvalue_ix(m, XDIR, u);
  y_ix = last_nonlarger_knotvalue_ix(m, YDIR, v);
  z_ix = last_nonlarger_knotvalue_ix(m, ZDIR, w);

  // adjustment of index if positioned _exactly_ at upper bound of grid
  if (x_ix == m.numDistinctKnots(XDIR) - 1 && u == m.maxParam(XDIR)) --x_ix;
  if (y_ix == m.numDistinctKnots(YDIR) - 1 && v == m.maxParam(YDIR)) --y_ix;
  if (z_ix == m.numDistinctKnots(ZDIR) - 1 && w == m.maxParam(ZDIR)) --z_ix;
  
  // checking if a valid corner was found
  if (x_ix < 0 || x_ix >= m.numDistinctKnots(XDIR) - 1) return false; // u outside domain
  if (y_ix < 0 || y_ix >= m.numDistinctKnots(YDIR) - 1) return false; // v outside domain
  if (z_ix < 0 || z_ix >= m.numDistinctKnots(ZDIR) - 1) return false; // v outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  while (true)
    {
      int x2_ix = search_downwards_for_nonzero_multiplicity(m, XDIR, x_ix, 
							    y_ix, z_ix);
      int y2_ix = search_downwards_for_nonzero_multiplicity(m, YDIR, y_ix, 
							    z_ix, x_ix);
      int z2_ix = search_downwards_for_nonzero_multiplicity(m, ZDIR, z_ix, 
							    x_ix, y_ix);
      if (x2_ix == x_ix && y2_ix == y_ix && z2_ix == z_ix)
	break;
      x_ix = x2_ix;
      y_ix = y2_ix;
      z_ix = z2_ix;
    }
  return true;
}

// =============================================================================
  bool Mesh3DUtils::identify_patch_upper_right(const Mesh3D&m,
					       double u, double v, double w,
					       int& x_ix, int& y_ix, int& z_ix)
// =============================================================================
{
  x_ix = first_larger_knotvalue_ix(m, XDIR, u);
  y_ix = first_larger_knotvalue_ix(m, YDIR, v);
  z_ix = first_larger_knotvalue_ix(m, ZDIR, w);
  
  // We do not need to adjust for a position _exactly_ at lower bound of grid as the value given by the
  // indices are guaranteed to be larger (or equal if at the end).

  // checking if a valid corner was found
  if (x_ix == 0 || x_ix >= m.numDistinctKnots(XDIR)) return false; // u outside domain
  if (y_ix == 0 || y_ix >= m.numDistinctKnots(YDIR)) return false; // v outside domain
  if (z_ix == 0 || z_ix >= m.numDistinctKnots(ZDIR)) return false; // w outside domain

  // We have now found the largest smaller knot in the u and v direction.  From here we
  // can search downwards to the lower-left corner of the containing mesh element, which
  // defines the sought-for "patch" of the surface.

  x_ix = search_upwards_for_nonzero_multiplicity(m, XDIR, x_ix, y_ix, z_ix);
  y_ix = search_upwards_for_nonzero_multiplicity(m, YDIR, y_ix, z_ix, x_ix);
  z_ix = search_upwards_for_nonzero_multiplicity(m, ZDIR, z_ix, x_ix, y_ix);
  
  return true;
}

// =============================================================================
int Mesh3DUtils::search_downwards_for_nonzero_multiplicity(const Mesh3D& m, 
							   Direction3D d, 
							   int start_ix,
							   int other1_ix,
							   int other2_ix)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix,
				      other1_ix, other1_ix + 1,
				      other2_ix, other2_ix + 1) == 0; --ix);

  return ix;
}

//==============================================================================
int Mesh3DUtils::search_upwards_for_nonzero_multiplicity(const Mesh3D& m, 
							 Direction3D d, 
							 int start_ix,
							 int other1_ix,
							 int other2_ix)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix,
							  other1_ix - 1, other1_ix,
							  other2_ix - 1, other2_ix) == 0; ++ix);

  return ix;  // @@ not yet tested!
}

// =============================================================================
int Mesh3DUtils::search_downwards_for_nonzero_multiplicity(const Mesh3D& m, 
							   Direction3D d, 
							   int start_ix, 
							   int other1_ix1, int other1_ix2, 
							   int other2_ix1, int other2_ix2)
// =============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != 0 && m.nu(d, ix,
				      other1_ix1, other1_ix2,
				      other2_ix1, other2_ix2) == 0; --ix);

  return ix;
}

//==============================================================================
int Mesh3DUtils::search_upwards_for_nonzero_multiplicity(const Mesh3D& m, 
							 Direction3D d, 
							 int start_ix, 
							 int other1_ix1, int other1_ix2, 
							 int other2_ix1, int other2_ix2)
//==============================================================================
{
  // provided that the mesh is a valid LR mesh, a valid index should always be found.
  int ix;
  for (ix = start_ix; ix != m.numDistinctKnots(d) && m.nu(d, ix,
							  other1_ix1, other1_ix2,
							  other2_ix1, other2_ix2) == 0; ++ix);

  return ix;  // @@ not yet tested!
}


// =============================================================================
int Mesh3DUtils::last_nonlarger_knotvalue_ix(const Mesh3D&m, Direction3D d, 
					     double par)
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
// int last_nonlarger_knotvalue_ix(const Mesh3D&m, Direction3D d, double par)
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
  int Mesh3DUtils::first_larger_knotvalue_ix(const Mesh3D& m, Direction3D d, 
					     double par)
// =============================================================================
{
  int ix = last_nonlarger_knotvalue_ix(m, d, par); 
  if (ix < m.numDistinctKnots(d)-1) ix++;
  return ix;
}

  
// =============================================================================
void Mesh3DUtils::sectionKnots(const Mesh3D& mesh, Direction3D d, int ix,
			       vector<vector<int> >& knots1,
			       vector<vector<int> >& knots2)
// =============================================================================
{
  Direction3D d1 = (d == XDIR) ? YDIR : ((d == YDIR) ? ZDIR : XDIR);
  Direction3D d2 = (d == XDIR) ? ZDIR : ((d == YDIR) ? XDIR : YDIR);
  bool atstart = (ix-mesh.firstMeshVecIx(d) < mesh.lastMeshVecIx(d)-ix);
  int beg1 = mesh.firstMeshVecIx(d1);
  int end1 = mesh.lastMeshVecIx(d1);
  for (int ix1=beg1; ix1<end1; ++ix1)
    {
      vector<int> curr_knots = LRBSpline3DUtils::derive_knots(mesh, d2, 
							      mesh.firstMeshVecIx(d2),
							      mesh.lastMeshVecIx(d2),
							      atstart ? ix : ix-1,
							      atstart ? ix+1 : ix,
							      ix1, ix1+1);
      knots2.push_back(curr_knots);
    }
  
  int beg2 = mesh.firstMeshVecIx(d2);
  int end2 = mesh.lastMeshVecIx(d2);
  for (int ix2=beg2; ix2<end2; ++ix2)
    {
      vector<int> curr_knots = LRBSpline3DUtils::derive_knots(mesh, d1, 
							      mesh.firstMeshVecIx(d1),
							      mesh.lastMeshVecIx(d1),
							      ix2, ix2+1,
							      atstart ? ix : ix-1,
							      atstart ? ix+1 : ix);
      knots1.push_back(curr_knots);
    }
  int stop_break = 1;
}

  
// =============================================================================
} // end namespace Go
// =============================================================================




