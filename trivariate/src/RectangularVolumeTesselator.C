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

#include "GoTools/trivariate/RectangularVolumeTesselator.h"
#include "GoTools/trivariate/SplineVolume.h"

#include <assert.h>
#include <fstream>

namespace Go
{

//===========================================================================
RectangularVolumeTesselator::~RectangularVolumeTesselator()
//===========================================================================
{
}

//===========================================================================
void RectangularVolumeTesselator::tesselate()
//===========================================================================
{
    tesselateVolume();
//    tesselateIsolines();
}

//===========================================================================
void RectangularVolumeTesselator::tesselateVolume()
//===========================================================================
{
//    int num_vert = mesh_->numVertices();
    int m = (mesh_->numStrips()/6) + 1;
//    int n = num_vert/(6*m);

    /// @@@ We can only tesselate properly rectangular-domain surfaces.
    // We extract the boundary surfaces for the volume.
    std::vector<shared_ptr<ParamSurface> >  bd_sfs = vol_.getAllBoundarySurfaces();
    for (size_t ki = 0; ki < bd_sfs.size(); ++ki)
    {
	shared_ptr<ParamSurface> bd_sf = bd_sfs[ki];
	RectDomain dom = bd_sf->containingDomain();
	int dim = bd_sf->dimension();
	Point pt(dim);
	for (int iu = 0; iu < m; ++iu) {
	    double ru = double(iu)/double(m-1);
	    double u = dom.umin()*(1.0-ru) + ru*dom.umax();
	    for (int iv = 0; iv < m; ++iv) {
		double rv = double(iv)/double(m-1);
		double v = dom.vmin()*(1.0-rv) + rv*dom.vmax();
		bd_sf->point(pt, u, v);
		//	    std::cout << pt << std::endl;
		int j;
		for (j=0; j<dim; ++j)
		    mesh_->vertexArray()[ki*3*m*m+(iv*m + iu)*3+j] = pt[j];
		for (; j<3; ++j)
		    mesh_->vertexArray()[ki*3*m*m+(iv*m + iu)*3 + j] = 0.0;
		mesh_->paramArray()[ki*2*m*m+(iv*m + iu)*2] = u;
		mesh_->paramArray()[ki*2*m*m+(iv*m + iu)*2+1] = v;
		if (mesh_->useNormals()) {
		    bd_sf->normal(pt, u, v);
		    mesh_->normalArray()[ki*3*m*m+(iv*m + iu)*3] = pt[0];
		    mesh_->normalArray()[ki*3*m*m+(iv*m + iu)*3 + 1] = pt[1];
		    mesh_->normalArray()[ki*3*m*m+(iv*m + iu)*3 + 2] = pt[2];
		}
		if (mesh_->useTexCoords()) {
		    mesh_->texcoordArray()[ki*2*m*m+(iv*m + iu)*2] = ru;
		    mesh_->texcoordArray()[ki*2*m*m+(iv*m + iu)*2+1] = rv;
		}
	    }
	}
    }
}


} // namespace Go
