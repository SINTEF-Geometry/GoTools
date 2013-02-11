//===========================================================================
//                                                                           
// File: RectangularSurfaceTesselator.C                                               
//                                                                           
// Created: Wed Nov 28 17:15:55 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectangularSurfaceTesselator.C,v 1.2 2009-01-15 13:18:02 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================
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
