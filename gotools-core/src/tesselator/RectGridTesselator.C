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

#include "GoTools/tesselator/RectGridTesselator.h"

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
