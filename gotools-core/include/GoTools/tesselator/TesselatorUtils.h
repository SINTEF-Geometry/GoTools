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

#ifndef _TESSELATORUTILS_H
#define _TESSELATORUTILS_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/LineCloud.h"

namespace Go {

/// Related to the relative resolution of tesselation
namespace TesselatorUtils
{
  /// Get the mesh size of a surface in the two parameter directions
  /// given the total number of nodes in the mesh. NB! u_nmb*v_nmb
  /// are probably not exactly equal to uv_nmb
  void getResolution(const ParamSurface *surf, 
		     int& u_nmb, int& v_nmb, int uv_nmb = 400);

  /// Fetch the control polygon of some geometric entity
  shared_ptr<LineCloud> getCtrPol(GeomObject* obj);

}  // of namespace TesselatorUtils
}; // end namespace Go
#endif // _TESSELATORUTILS_H
