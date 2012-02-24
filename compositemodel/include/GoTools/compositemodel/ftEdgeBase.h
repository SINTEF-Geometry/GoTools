//===========================================================================
//                                                                           
// File: ftEdgeBase.h                                                        
//                                                                           
// Created: Mon Jul  8 15:19:56 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftEdgeBase.h,v 1.8 2009-01-30 09:54:51 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

    /// Split the edge at the given parameter t. This edge will then represent
    /// the first part of the original edge, the second part is returned to the
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

    /// Add the edge edge to the edge loop after this edge
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

