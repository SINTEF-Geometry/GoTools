//===========================================================================
//                                                                           
// File: ftEdge.h                                                            
//                                                                           
// Created: Fri Mar  9 08:47:53 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftEdge.h,v 1.13 2009-06-12 08:55:14 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FTEDGE_H
#define _FTEDGE_H

#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftFaceBase.h"

namespace Go
{
    class SplineCurve;
    class Vertex;
    class EdgeVertex;

//===========================================================================
/** ftEdge - topological edge for Fantastic
 *
 * The ftEdge is a half-edge implementation of a topological data structure.
 * It implements the ftEdgeBase interface by using the Go geometry library.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \bug Not tested
 * \see ftEdgeBase
 */
//===========================================================================

///  \brief The ftEdge is a half-edge implementation of a topological 
/// data structure.

class ftEdge : public ftEdgeBase
{
public:

    /** Constructor.
     * Detailed description.
     */
    ftEdge(ftFaceBase* face, shared_ptr<ParamCurve> cv, 
	   double tmin, double tmax, int entry_id = -1);

    /// Constructor. Input is the curve representing the geometry of
    /// the edge and parameter values for the two vertices. There is
    /// no reference to a face object.
    ftEdge(shared_ptr<ParamCurve> cv, 
	   double tmin, double tmax, int entry_id = -1);

    /// Constructor.
    ftEdge(ftFaceBase* face, shared_ptr<ParamCurve> cv, 
	   shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, 
	   int entry_id = -1);

    /// Constructor. Input is the curve representing the geometry of
    /// the edge and two vertices. There is no reference to a face
    /// object.
    // @@@jbt - In STEP this corresponds to the entity
    // 'oriented_edge'.
    ftEdge(shared_ptr<ParamCurve> cv, 
	   shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, 
	   int entry_id = -1);

    /// Destructor
    ~ftEdge();

    /// Minimum parameter of curve restricted to edge
    virtual double tMin() const
    {
	return low_param_;
    }
    /// Maximum parameter of curve restricted to edge
    virtual double tMax() const
    {
	return high_param_;
    }
    //virtual void turnOrientation();
    //virtual void setOrientation();
    //virtual bool isTurned();

    virtual void setReversed(bool is_reversed);
    virtual bool isReversed();

    /// Reverse the geometric curve. A call to this function will
    /// reverse the direction of the geometry curve, swap the
    /// vertices, and change the 'is_reversed' flag. The geometry
    /// curve will be cloned.
    void reverseGeomCurve();

    /// The face associated to this edge
    virtual ftFaceBase* face();

    /// Set pointer to the face associated to this edge
    virtual void setFace(ftFaceBase* face)
	{
	    face_ = face;
	}

    /// Bounding box surrounding this edge. The box may be too large
    // The bounding box is not exact, it is much too large...
    // (in some cases). It is implemented as the bounding box
    // of the WHOLE edgecurve instead of only the piece covered
    // by this halfedge.
    virtual BoundingBox boundingBox();

    /// Split edge in a given parameter
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    virtual ftEdgeBase* split(double t);
#else
    virtual ftEdge* split(double t);
#endif

    shared_ptr<ftEdge> split2(double t);

    /// Split according to an already existing vertex
    shared_ptr<ftEdge> splitAtVertex(shared_ptr<Vertex> vx);

    /// Fetch Id corresponding to this edge. It is not necessarily uniquely set
    virtual int entryId() { return entry_id_; }

    /// Set Id corresponding to this edge
    virtual void setEntryId(int id) { entry_id_ = id; }

    /// Evaluate position given edge parameter
    virtual Point point(double t) const;

    /// Evaluate tangent of edge given edge parameter
    virtual Point tangent(double t) const;

    /// Evaluate normal of associated face given edge parameter
    virtual Point normal(double t) const;

    // Value of face_seed is altered par val of to closest_pt in sf.
    virtual Point normal(double t, Point& face_par_pt, double* face_seed) const;
    virtual void closestPoint(const Point& pt, double& clo_t,
			      Point& clo_pt, double& clo_dist,
			      double const *seed = 0) const;

    /// Connect this edge to the given one
    virtual void connectAfter(ftEdgeBase* edge);

    /// Closing of loop
    virtual void closeLoop(ftEdgeBase* last);

    /// Disconnect this edge from its neighbours in the loop of edges
    /// limiting a face
    virtual void disconnectThis();

    /// Represent adjacency by connecting this edge to its twin edge
    /// on the adjace face
    virtual void connectTwin(ftEdgeBase* twin, int& status);

    /// Remove twin information in this edge
    virtual void disconnectTwin();

    void joinVertices(ftEdgeBase* twin);

    /// Return edge pointer
    virtual ftEdge* geomEdge();

    /// Compare orientation of curve and edge
    virtual bool orientationOK() const;

    // New member

    /// Evaluate this edge at the parameter t, 
    /// including derivatives up to der
    void point(double t, int der, std::vector<Point>& derivs) const;

    /// Fetch the geometry curve associated with this edge
    shared_ptr<ParamCurve> geomCurve()
    {
	return geom_curve_;
    }

    /// Estimate the length of the curve corresponding to this edge
    double estimatedCurveLength();

    /// Estimate the length of the curve corresponding to this edge
    /// limited by edge parameters
    double estimatedCurveLength(double min_par, double max_par);

    /// Fetch the parameter in the associated face, given the edge
    /// parameter
    // Optional seed to closestBondaryPoint of face_.
    Point faceParameter(double t, double* seed = NULL) const;

    /// Convenience function to set or reset the geometry
    /// representation of the edge.
    /// \param geom_curve The curve representing the geometry of the
    /// edge. This will frequently be a CurveOnSurface.
    void setGeomCurve(shared_ptr<ParamCurve> geom_curve);

    virtual void updateGeomCurve(double tol);

    /// Update geometry info if possible
    bool updateEdgeInfo(double tol);

    /// Fetch index of geometry curve with respect to associated surface,
    /// if this information exists:
    /// -1 : No info or not a rectangular surface
    ///  0 : umin 
    ///  1 : umax 
    ///  2 : vmin 
    ///  3 : vmax 
    int getCurveIndex() const;

    /// Access function for vertices. The function takes into account
    /// whether or not the orientation is reversed.
    shared_ptr<Vertex> getVertex(bool at_start);
 
    // Access function for vertices. This function does not take into
    // account the orientation of the edge.
    void getVertices(shared_ptr<Vertex>& v1, 
		     shared_ptr<Vertex>& v2)
    {
      v1 = v1_;
      v2 = v2_;
    }

    /// Given one of the vertices belonging to this edge, fetch
    /// the other
    shared_ptr<Vertex> getOtherVertex(Vertex* vx)
      {
	shared_ptr<Vertex> dummy;
	if (vx == v1_.get())
	  return v2_;
	else if (vx == v2_.get())
	  return v1_;
	else 
	  return dummy;
      }
    
    shared_ptr<Vertex> getCommonVertex(ftEdge* other)
      {
	shared_ptr<Vertex> dummy;
	if (v1_.get() == other->v1_.get() || v1_.get() == other->v2_.get())
	  return v1_;
	else if (v2_.get() == other->v1_.get() || v2_.get() == other->v2_.get())
	  return v2_;
	else
	  return dummy;
      }
	
    bool hasVertex(Vertex* vx)
    {
      return (vx == v1_.get() || vx == v2_.get());
    }

    /// Assign verteices to this edge
    void setVertices(shared_ptr<Vertex> v1, 
		     shared_ptr<Vertex> v2);

    /// Fetch the parameter of the curve associated to this edge at a
    /// given vertex
    double parAtVertex(const Vertex* vx) const;

    /// Interchange one of this edge's vertex. Used in topology build.
    void replaceVertex(shared_ptr<Vertex>& this_vertex, 
		       shared_ptr<Vertex>& other_vertex);


    /// Functionality related to non-manifold models
    void addEdgeMultiplicityInstance(ftEdge *other);


    /// Check if this edge belongs to a radial edge (EdgeVertex)
    bool hasEdgeMultiplicity()
    {
      return (all_edges_.get() != 0);
    }

    /// Fetch the radial edge to which this edge belongs
    shared_ptr<EdgeVertex> getEdgeMultiplicityInstance()
      {
	return all_edges_;
      }

    /// Assign a radial edge to this edge. Used in topology build for
    /// volumes
    void setEdgeVertex(shared_ptr<EdgeVertex> radial_edge)
    {
      if (all_edges_.get() != radial_edge.get())
	all_edges_ = radial_edge;
    }
			 

    /// Connect to adjacent radial edges. Used in topology build for
    /// volumes
    void joinEdgeVertex(shared_ptr<EdgeVertex> radial_edge);

    /// Remove radial edge instance
    void removeEdgeVertex();

    /// Connect to adjacent vertices including updating edge information
    void joinVertex(shared_ptr<Vertex> this_vertex,
		    shared_ptr<Vertex> other_vertex); 
		    
    /// Get all (up to 2) faces meeting in this edge
    std::vector<ftSurface*> getAdjacentFaces() const;

    /// Get all faces meeting in this edge, including non-manifold cases
    std::vector<ftSurface*> getAllAdjacentFaces() const;

    /// Check for common vertices between two edges
    bool hasCommonVertices(ftEdge* other) const
    {
      if ((v1_.get() == other->v1_.get() && v2_.get() == other->v2_.get()) ||
	  (v1_.get() == other->v2_.get() && v2_.get() == other->v1_.get()))
	return true;
      else
	return false;
    }

    bool hasCommonRadialEdge(ftEdge* other) const
    {
      if (all_edges_.get() && other->all_edges_.get() &&
	  all_edges_.get() == other->all_edges_.get())
	return true;
      else
	return false;
    }

    bool checkEdgeTopology();

private:
    /// The face associated this edge
    ftFaceBase* face_;

    /// The geometrical representation of the curve associated this edge
    shared_ptr<ParamCurve> geom_curve_;

    /// Restriction of curve with regard to this edge, lower limit 
    /// parameter of the curve.
    double low_param_ ;

    /// Restriction of curve with regard to this edge, upper limit 
    /// parameter of the curve.
    double high_param_;

    /// Vertex corresponding to lower limit parameter
    shared_ptr<Vertex> v1_;

    /// Vertex corresponding to upper limit parameter
    shared_ptr<Vertex> v2_;

    int entry_id_;
    //int is_turned_;

    bool is_reversed_;

    /// This entity is intended for volume models or other non-manifold
    /// models. It is a radial edge collecting all half edges 
    /// meeting in an edge
    shared_ptr<EdgeVertex> all_edges_;

    ftEdge(ftFaceBase* face, shared_ptr<ParamCurve> cv, double t1,
	   shared_ptr<Vertex> v1, double t2, shared_ptr<Vertex> v2, 
	   int entry_id = -1);

};

} // namespace Go


#endif // _FTEDGE_H

