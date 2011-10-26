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
  void getResolution(const ParamSurface *surf, 
		     int& u_nmb, int& v_nmb, int uv_nmb = 400);

  std::shared_ptr<LineCloud> getCtrPol(GeomObject* obj);

}  // of namespace TesselatorUtils
}; // end namespace Go
#endif // _TESSELATORUTILS_H
