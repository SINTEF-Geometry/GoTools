//===========================================================================
//                                                                           
// File: TesselatorUtils.h
//                                                                           
// Created: Mars 2009                                         
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
