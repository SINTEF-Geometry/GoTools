//===========================================================================
//                                                                           
// File: LineCloudTesselator.C                                             
//                                                                           
// Created: Tue Oct 29 09:54:22 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/tesselator/LineCloudTesselator.h"

namespace Go
{


//===========================================================================
LineCloudTesselator::~LineCloudTesselator()
//===========================================================================
{
}

//===========================================================================
void LineCloudTesselator::tesselate()
//===========================================================================
{
    int numl = orig_cloud_.numLines();
    for (int i = 0; i < numl; ++i) {
	Vector3D raydir = orig_cloud_.point(2*i+1) - orig_cloud_.point(2*i);
	render_cloud_.point(2*i+1) = orig_cloud_.point(2*i) + raydir*scale_;
    }
}

} // namespace Go
