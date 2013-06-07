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

#ifndef _VOLUMEADJACENCY_H
#define _VOLUMEADJACENCY_H

#include "GoTools/compositemodel/Body.h"

namespace Go
{
  class ftSurface;
  class ParamSurface;

  /// \brief Adjacency analysis of volume models

  class VolumeAdjacency
  {
  public:
    /// Constructor giving tolerances for adjacency analysis
    VolumeAdjacency(double gap, double neighbour);

    /// Destructor
    ~VolumeAdjacency();

    /// Perform topology analysis on a set of bodies
    void setAdjacency(std::vector<shared_ptr<Body> >& solids);

    /// Perform topology analysis on a set of bodies where it is assumed
    /// the the topological relationship between the first new_solid_pos
    /// bodies are known already
    void setAdjacency(std::vector<shared_ptr<Body> >& solids, 
		      int new_solid_pos);

    private:
    /// Gap between volumes
    double gap_;

    /// Tolerance to check if two volumes can be adjacent at all, 
    /// neighbour_ > gap_
    double neighbour_;

    /// Check for adjacency between two solids
    void setAdjacency(shared_ptr<Body> solid1, shared_ptr<Body> solid2);

    /// Check for adjacency between two boundary faces, split faces in case
    /// of partial adjacency (one face embedded in the other, general partial
    /// coincidence is not handled)
    int faceAdjacency(shared_ptr<ftSurface> face1, 
		      shared_ptr<ftSurface> face2,
		      std::vector<shared_ptr<ftSurface> >& new_faces1,
		      std::vector<shared_ptr<ftSurface> >& new_faces2);

    /// Handle embedded boundary surfaces
    void splitSurface(shared_ptr<ParamSurface> srf1,
		      shared_ptr<ParamSurface> srf2,
		      std::vector<shared_ptr<ParamSurface> >& new_sfs);
  };

} // namespace Go


#endif // _VOLUMEADJACENCY_H
