//===========================================================================
//                                                                           
// File: RectangularVolumeTesselator.C                                       
//                                                                           
// Created: Thu Jul  5 15:57:47 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
    const SplineVolume* spline_vol =
	dynamic_cast<const SplineVolume*>(&vol_);
    assert(spline_vol != NULL);

    std::vector<shared_ptr<SplineSurface> > bd_sfs =
	spline_vol->getBoundarySurfaces();

    for (size_t ki = 0; ki < bd_sfs.size(); ++ki)
    {
	shared_ptr<SplineSurface> bd_sf = bd_sfs[ki];
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
