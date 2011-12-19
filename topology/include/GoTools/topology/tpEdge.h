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

//===========================================================================
/** tpEdge is an abstract interface for a half-edge topology structure.
 * Suggested (and minimal) structure of template edgeType when using 
 * tpTopologyTable.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \bug Not tested yet.
 * \see faceType
 */
//===========================================================================

#include "GoTools/topology/tpFace.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/topology/tpJointType.h"
#include "GoTools/topology/FaceConnectivity.h"

#include <memory>

class tpEdge
{
public:
    // These first members functions are needed by template class edgeType when
    // using the tpTopologyTable.

    /// Empty destructor.
    virtual ~tpEdge();
    tpEdge* next();
    tpEdge* twin();
    virtual double tMin() const = 0;
    virtual double tMax() const = 0;
    virtual void setOrientation() = 0;
/*     virtual bool isTurned() = 0; */
    virtual tpFace* face() = 0;
    virtual BoundingBox boundingBox() = 0;
    virtual void setEntryId(int id) = 0;
    /// Split should be implemented so that the new edge
    /// (returned by the function), is fully connected.
    /// It should also be created by operator new, so
    /// callers should keep in mind that they are
    /// responsible for managing the lifetime of the
    /// new edge.
    virtual tpEdge* split(double t) = 0;
    virtual Point point(double t) const = 0;
    virtual Point normal(double t) const = 0;
    virtual void closestpoint(const Point& pt, double& clo_t,
			      Point& clo_pt, double& clo_dist,
			      double const *seed = 0) const = 0;

    // The following member functions are not needed for the tpTopologyTable
    // to work, but are natural in this setting.
    tpEdge(); // Default constructor
    tpEdge* prev();
    virtual void turnOrientation() = 0;
    virtual int entryId() = 0;
    virtual void connectAfter(tpEdge* edge);
    virtual void closeLoop(tpEdge* last);
    virtual void disconnectThis();
    virtual void connectTwin(tpEdge* twin, int& status);
    virtual void disconnectTwin();
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
    virtual Point tangent(double t) const = 0;
    tpJointType checkContinuity(tpEdge* nextedge,
				double neighbour, double gap,
				double bend, double kink) const;

    bool hasConnectivityInfo()
    {
      return (connectivity_info_.get() ? true : false);
    }

    shared_ptr<FaceConnectivity<tpEdge> > getConnectivityInfo()
      {
	return connectivity_info_;
      }

    void setConnectivityInfo(shared_ptr<FaceConnectivity<tpEdge> > info)
    {
      connectivity_info_ = info;
    }

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


#endif // _TOPEDGE_H


