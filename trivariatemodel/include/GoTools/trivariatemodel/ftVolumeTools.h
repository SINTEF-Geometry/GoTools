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

namespace Go
{
  class ftVolume;
  class ftSurface;

  /// This namespace contains a function for splitting two surfaces  
  namespace ftVolumeTools
  {
    /// Split two volumes with regard to the intersections between 
    /// the boundary surfaces corresponding to these two volumes
    std::vector<std::shared_ptr<ftVolume> >
      splitVolumes(std::shared_ptr<ftVolume>& vol1, 
		   std::shared_ptr<ftVolume>& vol2, double eps);

    std::vector<std::shared_ptr<ftVolume> >
      splitVolumes(std::shared_ptr<ftVolume>& vol, 
		   std::shared_ptr<ftSurface>& face, double eps);

 
  }  // namespace ftVolumeTools
} // namespace Go

#endif // _FTVOLUMETOOLS_H
