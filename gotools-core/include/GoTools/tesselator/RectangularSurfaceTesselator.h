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

#ifndef RECTANGULARSURFACETESSELATOR_H
#define RECTANGULARSURFACETESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/config.h"

namespace Go
{

/** RectangularSurfaceTesselator: create a mesh for a boundary trimmed parametric
    surface. Visualization purposes.
*/

class GO_API RectangularSurfaceTesselator : public Tesselator
{
public:
  /// Constructor. Surface and mesh size are given. The tesselator can be set
  /// to compute also iso parametric curves with a specified mesh size.
 RectangularSurfaceTesselator(const ParamSurface& surf,
			      int ures = 20,
			      int vres = 20,
			      bool iso = false,
			      int uiso = 15,
			      int viso = 15,
			      int isores = 300)
	: surf_(surf),
	isolines_(iso), uiso_(uiso), viso_(viso), isores_(isores)
    {
	mesh_ = shared_ptr<RegularMesh>(new RegularMesh(ures, vres, true, true));
    }

  /// Destructor
    virtual ~RectangularSurfaceTesselator();
  
    virtual void tesselate();

    /// Fetch the computed mesh
    shared_ptr<RegularMesh> getMesh()
    {
	return mesh_;
    }

    /// Fetch tesselation of iso parametric curves
    std::vector<LineStrip>& getIsolineStrips()
    {
	return isolinestrips_;
    }

    //
    // 010430: I'm not sure if we're going to need these functions, an
    //         alternative is to trigger these actions when asking for
    //         pointers to the discretizations, which will be done by
    //         the 'painter' when that one is requested to redraw the scene...
    //         (jon)
    //

    /// Change mesh size
    void changeRes(int m, int n)
    {
	mesh_->resize(m, n);
	tesselateSurface();
    }

    /// Fetch mesh size
    void getRes(int& m, int& n)
    {
	n = mesh_->numStrips() + 1;
	m = mesh_->numVertices()/n;
    }

    /// Retesselate iso parametric curves
    void changeIsolines(bool isolines)
    {
	isolines_ = isolines;
	if (isolines_) {
	    tesselateIsolines();
	}
    }

    /// Fetch mesh of iso parametric curves
    void getIsolines(bool& isolines)
    {
	isolines = isolines_;
    }

    /// Change the number of iso parametric curves to compute
    void changeIsolineNumRes(int m, int n, int res)
    {
	uiso_ = m;
	viso_ = n;
	isores_ = res;
	tesselateIsolines();
    }
    /// Fetch the number of iso parametric curves to compute
    void getIsolineNumRes(int& m, int& n, int& res)
    {
	n = uiso_;
	m = viso_;
	res = isores_;
    }

private:
    void tesselateSurface();
    void tesselateIsolines();

    const ParamSurface& surf_;
    shared_ptr<RegularMesh> mesh_;
    std::vector<LineStrip> isolinestrips_;
    bool isolines_;
    int uiso_;
    int viso_;
    int isores_;
};

} // namespace Go


#endif // RECTANGULARSURFACETESSELATOR_H

