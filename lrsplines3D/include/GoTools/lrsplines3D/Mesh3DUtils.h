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

