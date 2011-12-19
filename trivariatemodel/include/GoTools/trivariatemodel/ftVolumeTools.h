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

  /// This namespace contains a function for splitting two surfaces  
  namespace ftVolumeTools
  {
    /// Split two volumes with regard to the intersections between 
    /// the boundary surfaces corresponding to these two volumes
    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol1, 
		   shared_ptr<ftVolume>& vol2, double eps);

    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol, 
		   shared_ptr<ftSurface>& face, double eps);

 
  }  // namespace ftVolumeTools
} // namespace Go

#endif // _FTVOLUMETOOLS_H
