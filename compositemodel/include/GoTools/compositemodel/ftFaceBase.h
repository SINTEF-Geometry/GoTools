//===========================================================================
//                                                                           
// File: ftFaceBase.h                                                        
//                                                                           
// Created: Mon Jul  8 15:20:08 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftFaceBase.h,v 1.8 2009-05-13 07:30:01 vsk Exp $
//                                                                           
// Description: Base class for all faces in Fantastic, to be used with topology.
//              
//              
//===========================================================================

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

