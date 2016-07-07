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

  /// This namespace contains service functions related to ftVolume
  namespace ftVolumeTools
  {
    /// Split two volumes with regard to the intersections between 
    /// the boundary surfaces corresponding to these two volumes
    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol1, 
		   shared_ptr<ftVolume>& vol2, double eps,
		   std::vector<int>& config);

    /// Split one volume according to intersections with a given face
    std::vector<shared_ptr<ftVolume> >
      splitVolumes(shared_ptr<ftVolume>& vol, 
		   shared_ptr<ftSurface>& face, double eps);

    /// Specific functionality. Used from ftVolume::generateMissingBdSurf
    void updateWithSplitFaces(shared_ptr<SurfaceModel> shell,
			      shared_ptr<ftSurface>& face1,
			      shared_ptr<ftSurface>& face2,
			      std::vector<std::pair<ftEdge*, ftEdge*> >& replaced_wires);

    /// Given a boundary or trimming surface related to an ftVolume, check the status and
    /// update stored information if any new boundary status information is computed
    /// \param vol the volume from which the boundary status is to be checked
    /// \param bd_face the boundary face where the requested information should be fetched/computed
    /// \param tol tolerance to check for surface coincidence
    /// \return -1 = not a boundary surface, 0 = boundary surface corresponding to umin,
    /// 1 = boundary surface corresponding to umax, 2 = boundary surface corresponding to vmin,
    /// 3 = boundary surface corresponding to vmax, 4 = boundary surface corresponding to wmin,
    /// 5 = boundary surface corresponding to wmax
    int boundaryStatus(ftVolume* vol, shared_ptr<ftSurface>& bd_face,
		       double tol);
 
  }  // namespace ftVolumeTools
} // namespace Go

#endif // _FTVOLUMETOOLS_H
