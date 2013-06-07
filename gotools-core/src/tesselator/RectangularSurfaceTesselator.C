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

#include "GoTools/tesselator/RectangularSurfaceTesselator.h"


namespace Go
{

//===========================================================================
RectangularSurfaceTesselator::~RectangularSurfaceTesselator()
//===========================================================================
{
}

//===========================================================================
void RectangularSurfaceTesselator::tesselate()
//===========================================================================
{
    tesselateSurface();
    tesselateIsolines();
}

//===========================================================================
void RectangularSurfaceTesselator::tesselateSurface()
//===========================================================================
{
    int m = mesh_->numStrips() + 1;
    int n = mesh_->numVertices()/m;
    /// @@@ We can only tesselate properly rectangular-domain surfaces.
    RectDomain dom = surf_.containingDomain();
    int dim = surf_.dimension();
    Point pt(dim);
    for (int iu = 0; iu < n; ++iu) {
	double ru = double(iu)/double(n-1);
	double u = dom.umin()*(1.0-ru) + ru*dom.umax();
	for (int iv = 0; iv < m; ++iv) {
	    double rv = double(iv)/double(m-1);
	    double v = dom.vmin()*(1.0-rv) + rv*dom.vmax();
	    surf_.point(pt, u, v);
	    //	    std::cout << pt << std::endl;
//	    std::cout << "iu: " << iu << ", iv: " << iv << ", pt: " << pt << std::endl;
	    int j;
	    for (j=0; j<dim; ++j)
	      mesh_->vertexArray()[(iv*n + iu)*3+j] = pt[j];
	    for (; j<3; ++j)
	      mesh_->vertexArray()[(iv*n + iu)*3 + j] = 0.0;
	    mesh_->paramArray()[(iv*n + iu)*2] = u;
	    mesh_->paramArray()[(iv*n + iu)*2+1] = v;
	    if (mesh_->useNormals()) {
		surf_.normal(pt, u, v);
		mesh_->normalArray()[(iv*n + iu)*3] = pt[0];
		mesh_->normalArray()[(iv*n + iu)*3 + 1] = pt[1];
		mesh_->normalArray()[(iv*n + iu)*3 + 2] = pt[2];
	    }
	    if (mesh_->useTexCoords()) {
		mesh_->texcoordArray()[(iv*n + iu)*2] = ru;
		mesh_->texcoordArray()[(iv*n + iu)*2+1] = rv;
	    }
	}
    }
}



//===========================================================================
void RectangularSurfaceTesselator::tesselateIsolines()
//===========================================================================
{
    if (!isolines_) return;
    isolinestrips_.resize(uiso_ + viso_);

    RectDomain dom = surf_.containingDomain();
    Point pt(3);
    for (int i = 0; i < uiso_; ++i) {
	isolinestrips_[i].resize(isores_);
	double ru = double(i)/double(uiso_-1);
	double u = dom.umin()*(1.0-ru) + ru*dom.umax();
	for (int j = 0; j < isores_; ++j) {
	    double rv = double(j)/double(isores_-1);
	    double v = dom.vmin()*(1.0-rv) + rv*dom.vmax();
	    surf_.point(pt, u, v);
	    isolinestrips_[i].vertexArray()[j*3] = pt[0];
	    isolinestrips_[i].vertexArray()[j*3 + 1] = pt[1];
	    isolinestrips_[i].vertexArray()[j*3 + 2] = pt[2];
	    isolinestrips_[i].paramArray()[j] = v;
	}
    }
    for (int i = 0; i < viso_; ++i) {
	isolinestrips_[uiso_ + i].resize(isores_);
	double rv = double(i)/double(viso_-1);
	double v = dom.vmin()*(1.0-rv) + rv*dom.vmax();
	for (int j = 0; j < isores_; ++j) {
	    double ru = double(j)/double(isores_-1);
	    double u = dom.umin()*(1.0-ru) + ru*dom.umax();
	    surf_.point(pt, u, v);
	    isolinestrips_[uiso_ + i].vertexArray()[j*3] = pt[0];
	    isolinestrips_[uiso_ + i].vertexArray()[j*3 + 1] = pt[1];
	    isolinestrips_[uiso_ + i].vertexArray()[j*3 + 2] = pt[2];
	    isolinestrips_[uiso_ + i].paramArray()[j] = u;
	}
    }
}

} // namespace Go
