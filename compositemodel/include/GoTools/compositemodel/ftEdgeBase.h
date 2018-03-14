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

#ifndef _FTEDGEBASE_H
#define _FTEDGEBASE_H

#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/topology/tpJointType.h"
#include "GoTools/topology/FaceConnectivity.h"
#include <vector>
#include <memory>

namespace Go
{

class ftEdge;
class ftFaceBase;

/** \brief Base class for edges. Defines the interface used in topology analysis.
 */
class ftEdgeBase //: public ftEdgeBase
{
public:
    /// Constructor
    ftEdgeBase();
    /// Destructor
    virtual ~ftEdgeBase();

    /// Next edge in the edge loop
    ftEdgeBase* next();
    /// Previous edge in the edge loop
    ftEdgeBase* prev();
    /// The corresponding edge on the neighbouring face. If no such
    /// neighbour exists, the function returns null.
    ftEdgeBase* twin();

    /// Start parameter value of edge
    virtual double tMin() const = 0;
    /// End parameter value of edge
    virtual double tMax() const = 0;
    //virtual void turnOrientation() = 0;
    //virtual void setOrientation() = 0;
    //virtual bool isTurned() = 0;

    /// Set a flag to indicate that the orientation of the edge is
    /// reversed with respect to the parametrization of the underlying
    /// ParamCurve.
    /// \param is_reversed value of the orientation flag
    virtual void setReversed(bool is_reversed) = 0;
    /// Get the value of the flag indicating that the orientation of
    /// the edge is reversed with respect to the parametrization of
    /// the underlying ParamCurve.
    /// \return \a true if orientation is reversed, \a false otherwise
    virtual bool isReversed() = 0;

    /// The face corresponding to this edge
    virtual ftFaceBase* face() = 0;
    /// Set face pointer
    virtual void setFace(ftFaceBase* face) = 0;

    /// Coordinate box surrounding this edge
    virtual Go::BoundingBox boundingBox() = 0;
    /// Id of edge. Default set to -1
    virtual int entryId() = 0;
    /// Set id of edge
    virtual void setEntryId(int id) = 0;

    /// Split the edge at the given parameter t. This edge then represent the
    /// first part of the original edge, the second part is returned to the
    /// caller.
    virtual ftEdgeBase* split(double t) = 0;

    /// Evaluate position
    virtual Go::Point point(double t) const = 0;
    /// Evaluate tangent
    virtual Go::Point tangent(double t) const = 0;
    /// Evaluate normal of associated face
    virtual Go::Point normal(double t) const = 0;
    /// Evaluate normal of associated face given a seed to the search for
    /// the corresponding face parameter
    virtual Go::Point normal(double t, Go::Point& face_par_pt, double* face_seed) const = 0;
    /// Closest point on this edge to a given point
    virtual void closestPoint(const Go::Point& pt, double& clo_t,
			      Go::Point& clo_pt, double& clo_dist,
			      double const *seed = 0) const = 0;

    /// Add this edge to the edge loop after "edge".
    virtual void connectAfter(ftEdgeBase* edge);
    /// Close edge loop by setting the appropriate pointers between this
    /// edge and last edge
    virtual void closeLoop(ftEdgeBase* last);
    /// Remove this edge from the edge loop it belongs to
    virtual void disconnectThis();
    /// Set adjacency between the face connected to this edge and the face
    /// connected to twin by setting twin pointers
    virtual void connectTwin(ftEdgeBase* twin, int& status);
    /// Disconnect this edge from its twin, i.e. remove adjacency information
    /// between the face connected to this edge and the face connected to
    /// the twin edge
    virtual void disconnectTwin();
    /// Check if this edge lies on the boundary of the associated face set, i.e.
    /// it has no neighbour
    bool onBoundary() { return (twin_) ? false : true; }
    /// Consider the 'node' or 'corner' implied by the startpoint
    /// (if at_start_of_edge is true) or endpoint (if it is false).
    /// The vector adjacent is filled with all edges with that 'node'
    /// as a startpoint (indicated by a true at the corresponding place
    /// in the at_start vector) or endpoint (false in at_start).
    /// The edge originally asked for is not included.
    void adjacentEdges(bool at_start_of_edge,
		       std::vector<ftEdgeBase*>& adjacent,
		       std::vector<bool>& at_start);

    /// Compute the continuity between this edge and nextedge given the
    /// necessary continuity parameters.
    /// 0 = g1 continuity, 1 = not quite g1, 2 = g0 continuity, 3 = gap,
    /// 4 = discontinuous, 5 = last segment
    tpJointType checkContinuity(ftEdgeBase* nextedge,
				double neighbour, double gap,
				double bend, double kink) const;

    /// Compare orientation of curve and edge
    virtual bool orientationOK() const;

    /// Return edge pointer
    virtual ftEdge* geomEdge() = 0; // The only new functionalty.

    /// Compute an eventual overlap between this edge and the edge other
    bool checkOverlap(ftEdgeBase *other, double tol, int nmbsample,
		      double& t1, double& t2, double& t3, 
		      double& t4, bool& same_dir,
		      bool no_snap=true) const;

    /// Check if this edge has access to continuity information regarding
    /// the joint between the face corresponding to this edge and the adjacent
    /// face along this edge
    bool hasConnectivityInfo()
    {
	return connectivity_info_.get() ? true : false;
    }

    /// Fetch the continuity information regarding the joint between the face 
    /// corresponding to this edge and the adjacent face along this edge
     shared_ptr<FaceConnectivity<ftEdgeBase> > getConnectivityInfo()
      {
	return connectivity_info_;
      }

    /// Store the continuity information regarding the joint between the face 
    /// corresponding to this edge and the adjacent face along this edge
    void setConnectivityInfo(shared_ptr<FaceConnectivity<ftEdgeBase> > info)
    {
      connectivity_info_ = info;
    }

    /// Mark the connectivity information associated to this edge as outdated
    void resetConnectivityInfo()
    {
      connectivity_info_.reset();
    }


protected:
    ftEdgeBase* next_;
    ftEdgeBase* prev_;
    ftEdgeBase* twin_;
    shared_ptr<FaceConnectivity<ftEdgeBase> > connectivity_info_;

};

} // namespace Go

#endif // _FTEDGEBASE_H

