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

#ifndef _RECTANGULARVOLUMETESSELATOR_H
#define _RECTANGULARVOLUMETESSELATOR_H


#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/RegularVolMesh.h"
#include "GoTools/utils/config.h"

namespace Go
{

/** RectangularSurfaceTesselator: create a mesh for a boundary trimmed parametric
    surface. Visualization purposes.
*/

class GO_API RectangularVolumeTesselator : public Tesselator
{
public:
  /// Constructor. Volume and mesh size are given. The tesselator can be set
  /// to compute also iso parametric curves with a specified mesh size.
 RectangularVolumeTesselator(const ParamVolume& vol,
			     int res = 50)//,
			     // int vres = 20,
			     // int wres = 20)// ,
				   // bool iso = false,
				   // int uiso = 15,
				   // int viso = 15,
				   // int isores = 300)
	: vol_(vol)// ,
	// isolines_(iso), uiso_(uiso), viso_(viso), isores_(isores)
    {
	mesh_ = shared_ptr<RegularVolMesh>(new RegularVolMesh(res, true, true));
    }

  /// Destructor
    virtual ~RectangularVolumeTesselator();
  
    // @@sbr Note: Edge vertices are constructed for each face,
    // meaning that the same vertice may occur more than once in the
    // triangulation, possibly resulting in some small artifacts along
    // edges.
    virtual void tesselate();

    /// Fetch the computed mesh
    shared_ptr<RegularVolMesh> getMesh()
    {
	return mesh_;
    }

    //
    // 010430: I'm not sure if we're going to need these functions, an
    //         alternative is to trigger these actions when asking for
    //         pointers to the discretizations, which will be done by
    //         the 'painter' when that one is requested to redraw the scene...
    //         (jon)
    //

    /// Change mesh size
    void changeRes(int m)
    {
	mesh_->resize(m);
	tesselateVolume();
    }

    /// Fetch mesh size
    void getRes(int& m)
    {
	m = (mesh_->numStrips()/6) + 1;
	// m = mesh_->numVertices()/n;
    }

private:
    // We store triangle strips for each of the 6 sides.
    void tesselateVolume();

    const ParamVolume& vol_;
    shared_ptr<RegularVolMesh> mesh_;
};

} // namespace Go


#endif // _RECTANGULARVOLUMETESSELATOR_H

