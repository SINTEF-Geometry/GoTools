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

#ifndef _TOPFACE_H
#define _TOPFACE_H

#include "GoTools/utils/Values.h"
#include "GoTools/utils/BoundingBox.h"
#include <memory>
#include <vector>

namespace Go
{

class tpEdge;


//===========================================================================
/** tpFace -  
 * Minimal structure of template faceType when using FaceAdjacency.
 * The interface is limited to the functionality required to compute
 * adjacency information in a set of faces and to remove this face
 * from a face set, i.e. remove all topology pointers related to this
 * face.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \see Class1
 */
//===========================================================================

class tpFace
{
public:
    // These first members functions are needed by template class faceType when
    // using the tpTopologyTable.

    /// Destructor
    virtual ~tpFace();
    /// Compute the edges associated to this face or fetch already existing
    /// edges
    virtual std::vector<shared_ptr<tpEdge> > 
      createInitialEdges(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    /// Return pointers to first part of all bd cvs.
    virtual std::vector<shared_ptr<tpEdge> > startEdges() = 0;
    /// Evaluate point on face
    virtual Point point(double u, double v) const = 0;
    /// Evaluate surface normal
    virtual Point normal(double u, double v) const = 0;
    /// The bounding box corresponding to this face
    virtual BoundingBox boundingBox() = 0;
    /// Return id, default id is -1
    virtual int getId() = 0;
    /// Remove all adjacency information related to this face
    virtual void isolateFace();
    //virtual std::vector<shared_ptr<tpEdge> > 
    //setOrientation(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    //void turnFace(std::vector<tpFace*>& turned);

    // The following member functions are not needed for the tpTopologyTable
    // to work, but are natural in this setting.
    //virtual void turnOrientation() = 0; // Called from our version of turnFace.

};

} // namespace Go

#endif // _TOPFACE_H

