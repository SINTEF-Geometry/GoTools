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

#ifndef _FTFACEBASE_H
#define _FTFACEBASE_H


#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/ftMessage.h"
#include "GoTools/igeslib/ftTangPriority.h"

namespace Go
{


// Forward declarations
class ftEdgeBase;
class ftSurface;

/** ftFaceBase -  An abstract interface to a topological face
 * 
 */
class GO_API ftFaceBase
{
public:
    /// Constructor
    ftFaceBase();
    /// Constructor
    ftFaceBase(int id/*, bool is_turned = false*/);
    /// Destructor
    virtual ~ftFaceBase();

    // Return as type ftSurface
    virtual ftSurface* asFtSurface();

    /// Reset loop information
    virtual void clearInitialEdges()
	{}  // Overriden when required

    // Evaluation and interrogation.
    /// Compute the edges associated to this face or fetch already existing
    /// edges
    virtual std::vector<shared_ptr<ftEdgeBase> > 
      createInitialEdges(double degenerate_epsilon = DEFAULT_SPACE_EPSILON,
			 double kink = 0.00015, bool no_split = false) = 0;
    /// Return pointers to first part of all bd cvs.
    virtual std::vector<shared_ptr<ftEdgeBase> > startEdges() = 0;
    /// Evaluate point on face
    virtual Point point(double u, double v) const = 0;
    /// Evaluate surface normal
    virtual Point normal(double u, double v) const = 0;
    /// The bounding box corresponding to this face
    virtual BoundingBox boundingBox() = 0;
    /// Set id for this face
    virtual void setId(int id);
    /// Return id, default id is -1
    virtual int getId();
    //virtual void turnOrientation() = 0;
    //virtual bool getOrientation() = 0;
    //virtual std::vector<shared_ptr<ftEdgeBase> > 
    //  setOrientation(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    //void turnFace(std::vector<ftFaceBase*>& turned);

    /// Fetch geometric surface
    virtual shared_ptr<ParamSurface> surface() = 0;
    /// Make sure that a geometric surface exists. Not applicable in normal
    /// configurations
    virtual ftMessage createSurf(double& max_error, double& mean_error) = 0;
    /// Fetch the error due to approximations in createSurf. Not applicable in normal
    /// configurations
    virtual void getError(double& max_error, double& mean_error) = 0;
    /// Priority of current face.  Not applicable in normal configurations
    virtual ftTangPriority getPrioType() const = 0;
    /// Update the boundary loops corresponding to this face with a new edge,
    /// i.e. one edge is split and all related topology information must be
    /// updated accordingly
    virtual void updateBoundaryLoops(shared_ptr<ftEdgeBase> new_edge);
    /// Remove all adjacency information related to this face
    virtual void isolateFace()
    {
      // Default no action
      ;
    }

#if 0
  // Functionality not through quality assurance. There are still problems
    /// Update topology based on changes in edge connectivity
    virtual void updateTopology(std::vector<ftEdgeBase*> removed_edgs)
    {
      // Default no action
      ;
    }
#endif

    /// Close gap between adjacent faces
    virtual ftMessage removeGap(ftEdgeBase* e1, ftEdgeBase* e2, ftFaceBase *other)
	{ return FT_NOT_SUPPORTED; }

    /// Closest point between this face and a point
    virtual void closestPoint(const Point& pt,
		      double&  clo_u,
		      double&  clo_v, 
		      Point& clo_pt,
		      double&  clo_dist,
		      double   epsilon) const = 0;


 
protected:
    int id_;
    //bool is_turned_;
};

} // namespace Go


#endif // _FTFACEBASE_H

