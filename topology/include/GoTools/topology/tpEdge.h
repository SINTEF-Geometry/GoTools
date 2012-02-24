//===========================================================================
//                                                                           
// File: tpEdge.h                                                            
//                                                                           
// Created: Tue Mar 21 15:24:20 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: tpEdge.h,v 1.27 2009-01-30 09:53:07 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _TOPEDGE_H
#define _TOPEDGE_H

#include "GoTools/topology/tpFace.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/topology/tpJointType.h"
#include "GoTools/topology/FaceConnectivity.h"

#include <memory>

namespace Go
{

//===========================================================================
/** tpEdge is an abstract interface for a half-edge topology structure.
 * Suggested (and minimal) structure of template edgeType when using 
 * FaceAdjacency.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \see faceType
 */
//===========================================================================

class tpEdge
{
public:
    // These first members functions are needed by template class edgeType when
    // using the tpTopologyTable.

    /// Empty destructor.
    virtual ~tpEdge();
    /// Next edge in the edge loop
    tpEdge* next();
    /// The corresponding edge on the neighbouring face. If no such
    /// neighbour exists, the function returns null.
    tpEdge* twin();
    /// Start parameter value of edge
    virtual double tMin() const = 0;
    /// End parameter value of edge
    virtual double tMax() const = 0;
    /*    virtual void setOrientation() = 0; */
/*     virtual bool isTurned() = 0; */
    /// The face corresponding to this edge
    virtual tpFace* face() = 0;
    /// Coordinate box surrounding this edge
    virtual BoundingBox boundingBox() = 0;
    /// Id of edge. Default set to -1
    virtual void setEntryId(int id) = 0;
    /// Split should be implemented so that the new edge
    /// (returned by the function), is fully connected.
    /// It should also be created by operator new, so
    /// callers should keep in mind that they are
    /// responsible for managing the lifetime of the
    /// new edge.
    virtual tpEdge* split(double t) = 0;
    /// Evaluate position
    virtual Point point(double t) const = 0;
    /// Evaluate normal of associated face
    virtual Point normal(double t) const = 0;
    /// Closest point on this edge to a given point
    virtual void closestpoint(const Point& pt, double& clo_t,
			      Point& clo_pt, double& clo_dist,
			      double const *seed = 0) const = 0;

    // The following member functions are not needed for the tpTopologyTable
    // to work, but are natural in this setting.
    /// Default constructor
    tpEdge();
    /// Previous edge in the edge loop
    tpEdge* prev();
    /*    virtual void turnOrientation() = 0; */
   /// Get id of edge
     virtual int entryId() = 0;
    /// Add the edge edge to the edge loop after this edge
    virtual void connectAfter(tpEdge* edge);
    /// Close edge loop by setting the appropriate pointers between this
    /// edge and last edge
    virtual void closeLoop(tpEdge* last);
    /// Remove this edge from the edge loop it belongs to
     virtual void disconnectThis();
    /// Set adjacency between the face connected to this edge and the face
    /// connected to twin by setting twin pointers
    virtual void connectTwin(tpEdge* twin, int& status);
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
		       std::vector<tpEdge*>& adjacent,
		       std::vector<bool>& at_start);
    // The tangent is only of importance when calling tpUtils::checkContinuity().
    /// Evaluate tangent
    virtual Point tangent(double t) const = 0;
    /// Compute the continuity between this edge and nextedge given the
    /// necessary continuity parameters.
    /// 0 = g1 continuity, 1 = not quite g1, 2 = g0 continuity, 3 = gap,
    /// 4 = discontinuous, 5 = last segment
    tpJointType checkContinuity(tpEdge* nextedge,
				double neighbour, double gap,
				double bend, double kink) const;

    /// Check if this edge has access to continuity information regarding
    /// the joint between the face corresponding to this edge and the adjacent
    /// face along this edge
     bool hasConnectivityInfo()
    {
      return (connectivity_info_.get() ? true : false);
    }

    /// Fetch the continuity information regarding the joint between the face 
    /// corresponding to this edge and the adjacent face along this edge
    shared_ptr<FaceConnectivity<tpEdge> > getConnectivityInfo()
      {
	return connectivity_info_;
      }

     /// Store the continuity information regarding the joint between the face 
    /// corresponding to this edge and the adjacent face along this edge
   void setConnectivityInfo(shared_ptr<FaceConnectivity<tpEdge> > info)
    {
      connectivity_info_ = info;
    }

     /// Mark the connectivity information associated to this edge as outdated
   void resetConnectivityInfo()
    {
      connectivity_info_.reset();
    }

    // The class should have some notion describing neighbourhood of an edge.
protected:
    tpEdge* next_;
    tpEdge* prev_;
    tpEdge* twin_;
    shared_ptr<FaceConnectivity<tpEdge> > connectivity_info_;

};

} // namespace Go
#endif // _TOPEDGE_H


