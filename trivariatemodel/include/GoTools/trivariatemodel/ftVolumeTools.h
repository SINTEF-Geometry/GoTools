//===========================================================================
//
// File : ftVolumeTools.h
//
// Created: Sept. 2010
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


#ifndef _FTVOLUMETOOLS_H
#define _FTVOLUMETOOLS_H

#include <vector>
#include <memory>
#include "GoTools/utils/config.h"

namespace Go
{
  class ftVolume;
  class ftSurface;
  class ftEdge;
  class SurfaceModel;

  /// This namespace contains a function for splitting of volumes
  namespace ftVolumeTools
  {
    /// Split two volumes with regard to the intersections between 
    /// the boundary surfaces corresponding to these two volumes
    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol1, 
		   shared_ptr<ftVolume>& vol2, double eps);

    /// Split one volume according to intersections with a given face
    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol, 
		   shared_ptr<ftSurface>& face, double eps);

    /// Specific functionality. Used from ftVolume::generateMissingBdSurf
    void updateWithSplitFaces(shared_ptr<SurfaceModel> shell,
			      shared_ptr<ftSurface>& face1,
			      shared_ptr<ftSurface>& face2,
			      std::vector<std::pair<ftEdge*, ftEdge*> >& replaced_wires);

 
  }  // namespace ftVolumeTools
} // namespace Go

#endif // _FTVOLUMETOOLS_H
