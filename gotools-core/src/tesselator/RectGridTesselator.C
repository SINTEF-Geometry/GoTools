//===========================================================================
//                                                                           
// File: RectGridTesselator.C                                              
//                                                                           
// Created: Wed Jan 19 13:16:33 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectGridTesselator.C,v 1.1 2007-04-17 12:25:56 sbr Exp $
//                                                                           
//===========================================================================



#include "GoTools/tesselator/RectGridTesselator.h"

using std::shared_ptr;
namespace Go
{


//===========================================================================
RectGridTesselator::RectGridTesselator(const RectGrid& rg)
//===========================================================================
  : rectgrid_(rg)
{
  quadmesh_ = 
    shared_ptr<QuadMesh>(new QuadMesh(rg.numCoefs_u()*rg.numCoefs_v(), 
				      (rg.numCoefs_u() - 1) * (rg.numCoefs_v() - 1),
				      false, false));
}

//===========================================================================
RectGridTesselator::~RectGridTesselator()
//===========================================================================
{
}

//===========================================================================
void RectGridTesselator::tesselate()
//===========================================================================
{
    DEBUG_ERROR_IF(rectgrid_.dimension() != 3, "Dimension must be 3.");
    int numu = rectgrid_.numCoefs_u();
    int numv = rectgrid_.numCoefs_v();
    std::copy(rectgrid_.rawData(),
	      rectgrid_.rawData() + numu*numv*3,
	      quadmesh_->vertexArray());
    unsigned int* q = quadmesh_->quadIndexArray();
    int c = 0;
    for (int i = 1; i < numv; ++i) {
	for (int j = 1; j < numu; ++j) {
	    q[c++] = (i-1)*numu + (j-1);
	    q[c++] = (i)*numu + (j-1);
	    q[c++] = (i)*numu + (j);
	    q[c++] = (i-1)*numu + (j);
	}
    }
}

// //===========================================================================
// void RectGridTesselator::tesselate()
// //===========================================================================
// {
//     DEBUG_ERROR_IF(rectgrid_.dimension() != 3, "Dimension must be 3.");
//     int numu = rectgrid_.numCoefs_u();
//     int numv = rectgrid_.numCoefs_v();
//     int numl = 2*numu*numv - numu - numv;
//     // Build the line arrays. Two points define one line.
//     std::vector<Vector3D> p(2*numl);
//     int pct = 0;
//     const double* rp = rectgrid_.rawData();
//     // First the const-u-lines:
//     for (int i = 0; i < numv-1; ++i) {
// 	for (int j = 0; j < numu; ++j) {
// 	    p[pct][0] = rp[(i*numu + j)*3];
// 	    p[pct][1] = rp[(i*numu + j)*3 + 1];
// 	    p[pct][2] = rp[(i*numu + j)*3 + 2];
// 	    ++pct;
// 	    p[pct][0] = rp[((i+1)*numu + j)*3];
// 	    p[pct][1] = rp[((i+1)*numu + j)*3 + 1];
// 	    p[pct][2] = rp[((i+1)*numu + j)*3 + 2];
// 	    ++pct;
// 	}
//     }
//     // Then the const-v-lines:
//     for (int i = 0; i < numv; ++i) {
// 	for (int j = 0; j < numu-1; ++j) {
// 	    p[pct][0] = rp[(i*numu + j)*3];
// 	    p[pct][1] = rp[(i*numu + j)*3 + 1];
// 	    p[pct][2] = rp[(i*numu + j)*3 + 2];
// 	    ++pct;
// 	    p[pct][0] = rp[(i*numu + j+1)*3];
// 	    p[pct][1] = rp[(i*numu + j+1)*3 + 1];
// 	    p[pct][2] = rp[(i*numu + j+1)*3 + 2];
// 	    ++pct;
// 	}
//     }
//     ASSERT(pct == 2*numl);
//     render_cloud_.setCloud(p[0].begin(), numl);
// }

} // namespace Go
