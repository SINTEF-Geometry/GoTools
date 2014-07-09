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

#ifndef _FTSURFACE_H
#define _FTSURFACE_H


#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/SurfaceModel.h"

// class ftFaceBase;

namespace Go
{
    class SplineCurve;
    class Loop;
    class PointOnEdge;
    class Body;
    class EdgeVertex;


    /// Struct to store information about adjacency relations between two faces
    struct AdjacencyInfo
    {
    public:
      /// True if adjacency was found, false if not
      bool adjacency_found_;
      /// Boundary index of first surface. 0 = umin, 1 = umax, 2 = vmin, 3 = vmax
      int bd_idx_1_;
      /// Boundary index of second surface. 0 = umin, 1 = umax, 2 = vmin, 3 = vmax
      int bd_idx_2_;
      /// True if the surfaces are equally oriented along the common boundary
      bool same_orient_;
      /// True if adjacency is found, user wanted to test if the surfaces meet
      /// in corner-to-corner configuration, and the test failed. Otherwise false.
      bool corner_failed_;

      AdjacencyInfo() :
      corner_failed_(false)
      { }
    };  // End struct AdjacencyInfo


//===========================================================================
/** ftSurface -  A topological face. The class also provides an interface for operations on the assosiated surface.
 * Detailed description.
 *
 * \author Vibeke Skytt
 * \see ftFaceBase
*/
//===========================================================================

class GO_API ftSurface : public ftFaceBase
{
public:
    /// Constructor.  Typically, the ParamSurface in this constructor
    /// will be a BoundedSurface. A loop will be built from this, but
    /// if we already have a loop and the edges in the loop has
    /// twins, then we must explicitly set this loop with
    /// addBoundaryLoops(). Note: In STEP, ftSurface
    /// corresponds to the entity 'advanced_face'.
    ftSurface(shared_ptr<ParamSurface> sf, int id);

    /// Constructor. Input is the geometrical surface representing the
    /// underlying face, and the topological loop that bounds the
    /// face. Typically, the ParamSurface in this constructor will be
    /// a SplineSurface or an ElementarySurface. Note: In STEP, ftSurface
    /// corresponds to the entity 'advanced_face'.
    ftSurface(shared_ptr<ParamSurface> sf,
	      shared_ptr<Loop> loop,
	      int id = -1);

    /// Empty destructor
    virtual ~ftSurface();

    /// Return as type ftSurface
    virtual ftSurface* asFtSurface();

    /// Reset loop information
    virtual void clearInitialEdges();

    // Evaluation and interrogation. Overridden from ftFaceBase.
    virtual std::vector<shared_ptr<ftEdgeBase> >
    createInitialEdges(double degenerate_epsilon = DEFAULT_SPACE_EPSILON,
		       double kink = 0.00015, bool no_split = false);
/*     virtual std::vector<shared_ptr<ftEdgeBase> > */
/*     setOrientation(double degenerate_epsilon = DEFAULT_SPACE_EPSILON); */
    virtual std::vector<shared_ptr<ftEdgeBase> > startEdges();
    /// Evaluate point on face
    virtual Point point(double u, double v) const;
    /// Evaluate surface normal
    virtual Point normal(double u, double v) const;
    /// The bounding box corresponding to this face
    virtual BoundingBox boundingBox();
    /// Set Id of this face
    virtual void setId(int id);
    /// Fetch Id of this face, may not be uniquely set
    virtual int getId();
    /// Turn orientation of the corresponding surface
    //virtual void turnOrientation();
    //virtual bool getOrientation();

    /// Get underlying surface
    virtual shared_ptr<ParamSurface> surface()
    { return surf_; }

    /// By inheritage, not relevant for this entity
    virtual ftMessage createSurf(double& max_error, double& mean_error);
    /// By inheritage, not relevant for this entity
    virtual void getError(double& max_error, double& mean_error);
    /// Priority type, not relevant for this entity
    virtual ftTangPriority getPrioType() const;

    /// Priority type, not relevant for this entity
    void setPrioType(ftTangPriority type)
    { prio_type_ = type; }

    /// Update the information stored in the boundary loops
    virtual void updateBoundaryLoops(shared_ptr<ftEdgeBase> new_edge);

    // Split vertices
    virtual void isolateFace();

    // Added in this class

    /// Set information about boundary loops. Intended for use when topology information
    /// exist prior to building a GoTools face set (i.e. SurfaceModel)
    /// Mark that the function will throw if the loop information is inconsistent with the
    /// already existing information in surf_ or the given loops are not consistent with the rules
    /// If more than one loop is given, the first loop is the outer one. Subsequent loops must lie
    /// inside the outer loop. The loops may not intersect.
    void addBoundaryLoops(std::vector<shared_ptr<Loop> >& bd_loops);

    /// Fetch the outer boundary loop
    void addOuterBoundaryLoop(shared_ptr<Loop> outer_loop);

    /// Number of loops, the first is the outer boundary loop, further loops
    /// represents holes
    int nmbBoundaryLoops()
	{
	    return (int)(boundary_loops_.size());
	}

    /// Get the specified boundary loops
    shared_ptr<Loop> getBoundaryLoop(int idx)
	{
	  shared_ptr<Loop> dummy;
	  return (idx >= 0 && idx < (int)(boundary_loops_.size())) ? 
	    boundary_loops_[idx] : dummy;
	}

    /// Number of edges in all loops
    int nmbEdges() const;

    /// Fetch all edges in all loops
    std::vector<shared_ptr<ftEdge> > getAllEdges() const;
 
    /// Fetch all edges in given loop
    std::vector<shared_ptr<ftEdge> > getAllEdges(int loop_idx) const;
 
    /// Fetch pointers to all edges in all loops
    std::vector<ftEdge*> getAllEdgePtrs() const;

    /// Fetch pointers to all edges in specified loop
    std::vector<ftEdge*> getAllEdgePtrs(int loop_idx) const;

    /// Check if this face contains any holes
    bool onlyOuterTrim() const
    {
      return (boundary_loops_.size() == 1);
    }

    /// Count the number of curves in the outer boundary loop with regard
    /// to how many corners there is in this loop, but without regard to
    /// how the loop is actually divided into curves
    int nmbOuterBdCrvs(double gap, double neighbour, double angtol) const;

    /// Approximate a regular face with a non-trimmed spline surface
    /// If the initial surface is not regular, no output is created
    shared_ptr<ParamSurface> getUntrimmed(double gap, double neighbour, 
					  double angtol, bool only_corner=false);

     /// Closest point between this face and a point
    virtual void closestPoint(const Point& pt,
			      double&  clo_u,
			      double&  clo_v, 
			      Point& clo_pt,
			      double&  clo_dist,
			      double   epsilon) const;

    /// Closest point between a domain of this face and a point.
    void closestPoint(const Point& pt,
		      double&  clo_u,
		      double&  clo_v, 
		      Point& clo_pt,
		      double&  clo_dist,
		      double   epsilon,
		      const RectDomain* domain_of_interest,
		      double   *seed) const;

    /// Closest point between a boundary on this face and a point
    ftEdgeBase* closestBoundaryPoint(const Point& pt,
				     const Point& in_vec,
				     double&  clo_u,
				     double&  clo_v, 
				     Point& clo_pt,
				     double&  clo_dist,
				     double& clo_par) const;

    /// Closest point between the outer boundary on this face and a point
    ftEdgeBase* closestOuterBoundaryPoint(const Point& pt,
					  const Point& in_vec,
					  double&  clo_u,
					  double&  clo_v, 
					  Point& clo_pt,
					  double&  clo_dist,
					  double& clo_par) const;

    /// Return the edge on this face closest to a point
    ftEdgeBase* edgeClosestToPoint(double u, double v);

    /// Write surface to stream
    void write(std::ostream& os);
 
    /// Smooth the face, the associated surface is required 
    /// to be of type SplineSurface.
    ftMessage smoothOutFace(int edge_cont, double approx_orig_tol, 
			    double deg_tol, double& maxerr, double& meanerr,
			    double approx_weight = 0.8);

    /// Return a part of the boundary of the current surface, specified by
    /// two points belonging to this boundary curve
    shared_ptr<SplineCurve> 
      getBoundaryPiece(Point& pt1, Point& pt2, double eps);


    /// Check for constant parameter line kinks
    bool getSurfaceKinks(double angtol, std::vector<double>& g1_disc_u,
			 std::vector<double>& g1_disc_v);

    /// Check for constant parameter line kinks, i.e. C1 discontinuities
    bool getSurfaceDisconts(double tol, std::vector<double>& disc_u,
			    std::vector<double>& disc_v);

    /// Split a surface along constant parameter line kinks. Currently
    /// only spline surface is implemented.
    std::vector<shared_ptr<ftSurface> > splitAlongKinks(double angtol);

    /// Close the gap between adjacent surfaces if the configuration allows it
    virtual ftMessage removeGap(ftEdgeBase* e1, ftEdgeBase* e2, ftFaceBase *other,
				double epsge);

    /// Test if a Vertex is close to the surface within some tolerance
    bool isClose(shared_ptr<Vertex> v, double tol) const;

    /// Get all vertices, duplicates removed
    std::vector<shared_ptr<Vertex> > vertices() const;

    /// Get non corner vertices
    std::vector<shared_ptr<Vertex> > getNonCornerVertices(double kink) const;

    /// Get non corner vertices restricted to a particular loop
    std::vector<shared_ptr<Vertex> > getNonCornerVertices(double kink,
								 int loop_idx) const;

    /// Get corner vertices
    std::vector<shared_ptr<Vertex> > getCornerVertices(double kink) const;
    /// Get corner vertices restricted to a particular loop
    std::vector<shared_ptr<Vertex> > getCornerVertices(double kink,
							      int loop_idx) const;

    /// Get all vertices common to this face and another face
    std::vector<shared_ptr<Vertex> > 
      getCommonVertices(ftSurface* other) const;

    /// Get all edges common to this face and another face represented by
    /// the half edges in this face
    std::vector<shared_ptr<ftEdge> > getCommonEdges(ftSurface *other) const;

    /// Get the vertex closest to a given point
    shared_ptr<Vertex> getClosestVertex(const Point& pnt) const;

    /// Collect all pairs of surface and vertex points where distance is greater than a tolerance
    void getBadDistance(std::vector<std::pair<ftSurface*, shared_ptr<Vertex> > >& badPairs,
			double tol);

    /// Collect all pairs of surface and edge where some part of the 
    /// edge has a distance to the surface greater than a tolerance
    void getBadDistance(std::vector<std::pair<ftSurface*, ftEdge*> >& badPairs,
			double tol) const;

    /// Collect all pairs of edge and vertex where the vertex has a 
    /// distance to the edge greater than a tolerance
    void getBadDistance(std::vector<std::pair<ftEdge*, shared_ptr<Vertex> > >& badPairs,
			double tol) const;

    /// Information about discontinuites, used in quality checing of a 
    /// surface model
    void getPosTangentSurfaceDiscont(std::vector<ftEdge*>& badPos,
				     std::vector<ftEdge*>& badTangent,
				     double tol, double kink, double bend, int leastSurfIndex,
				     shared_ptr<SurfaceModel> sm) const;

    /// Check consistency of boundary loops
    bool checkLoopOrientation(std::vector<shared_ptr<Loop> >& inconsistent_loops) const;

    /// Used in quality checing of a surface model
    bool hasAcuteAngle(ftEdge* along_edge, double angtol) const;

    ///  Get pairs of close boundary loop points
    void getNarrowRegion(double gap_tol, double tol, 
			 std::vector<std::pair<shared_ptr<PointOnEdge>, 
			 shared_ptr<PointOnEdge> > >& narrow_pt);

    /// Check and fix orientation of boundary loops of trimmed surfaces
    /// \return \c true if the orientation of a boundary loop was
    /// fixed, \c false otherwise
    bool checkAndFixBoundaries();

    /// Compute face area
    double area(double tol) const;

    /// Get neighbouring faces
    void getAdjacentFaces(std::vector<ftSurface*>& neighbours) const;

    /// Fetch tolerance used in check for degenerace
    double getCurrEps() const
	{
	    return degenerate_eps_;
	}


    /// Boundary model / volume model, set associated body
    void setBody(Body *body)
    {
      body_ = body;
    }

    /// Fetch associated body, if any
    Body* getBody() const
    {
      return body_;
    }

    /// Whether or not this face belongs to a solid
    bool hasBody() const
    {
      return (body_ != 0);
    }

    /// Non-manifold, set adjacency information between two bodies
    ftSurface *twin()
    {
      return twin_;
    }

    /// Get all (up to 2) adjacent bodies
    std::vector<Body*> getAdjacentBodies() const;

    /// Set adjacency between bodies
    /// Radial edge information is existed to exist
    void setTwin(ftSurface* newtwin);

    /// Set adjacency between bodies
    /// Radial edge information is set
    bool connectTwin(ftSurface* newtwin, double tol, bool no_snap=true);

    /// Set adjacency between bodies
    /// Correspondence of surface loops is given as input
    void connectTwin(ftSurface* newtwin, std::vector<ftEdgeBase*> first1, 
		     std::vector<ftEdgeBase*> first2, 
		     std::vector<int> forward);

    /// Disconnect twin info. NB radial edge information is NOT removed
    void disconnectTwin();
   
   /// Check if two faces are adjacent, and return information about
    /// the edge along which this acjacency ocurrs
    bool areNeighbours(ftSurface *other, shared_ptr<ftEdge>& edge1, 
		       shared_ptr<ftEdge>& edge2, int adj_idx = 0) const;

    /// Number of edges along which the two faces are adjacent
    int nmbAdjacencies(ftSurface *other) const;

    /// Number neighbouring faces common to the two given faces
    int nmbNextNeighbours(ftSurface *other) const;

    /// Check if the face is represented as a spline
    bool isSpline() const;

    /// Check if two bodies are neighbours, are splines and have common
    /// a spline space at the interface
    bool commonSplineSpace(ftSurface *other, double tol);

    /// Ensure that two neighbouring spline surfaces have a common
    /// spline space at the interface
    /// Averaging of coefficients at the common boundary is also performed
    void makeCommonSplineSpace(ftSurface *other);

    /// Check if two volumes are splines and meet in a corner-to-corner
    /// configuration
    bool isCornerToCorner(ftSurface* other, double tol, int adj_idx = 0);

    /// Ensure that two spline volumes meet in a corner to corner
    /// configuration
    void splitAtInternalCorner(ftSurface* other,
			       std::vector<shared_ptr<ftSurface> >& new_face1,
			       std::vector<shared_ptr<ftSurface> >& new_face2,
			       double tol=DEFAULT_SPACE_EPSILON);

    /// DEPRECATED METHOD. USE THE OTHER getAdjacencyInfo() instead
    /// Fetch info on adjacency between neighbouring faces
    /// \param other The other face
    /// \param tol Adjacency tolerance
    /// \param bd1 Boundary index of first surface 
    ///  0 = umin, 1 = umax, 2 = vmin,  3 = vmax
    /// \param bd2 Boundary index of second surface 
    ///  0 = umin, 1 = umax, 2 = vmin,  3 = vmax
    /// \param same_orient Indicates whether the surfaces are equally oriented along
    /// the common boundary
    /// \return \a true if adjacency is found and both faces are non-trimmed or
    /// boundary trimmed
    bool getAdjacencyInfo(ftSurface *other, double tol,
			  int& bd1, int& bd2, bool& same_orient);

    /// Fetch info on adjacency between neighbouring faces
    AdjacencyInfo getAdjacencyInfo(ftSurface *other, double tol,
				   int adj_idx = 0, bool test_corner = false);

    /// Fetch info on adjacency between neighbouring faces given information
    /// about the common edge
    AdjacencyInfo getAdjacencyInfo(ftEdge *edge, ftSurface *other, double tol);

/// Check if two degenerate surface boundaries meet in a vertex
    bool checkDegAdjacency(ftSurface *other, shared_ptr<Vertex> vx,
			   double tol,
			   shared_ptr<ParamCurve>& bdcv1, 
			   shared_ptr<ParamCurve>& bdcv2);

    /// Given two adjacent spline faces represented in the same spline space,
    /// return local enumeration of the coefficients. This surface is given 
    /// first and the other surface second
    bool getCorrCoefEnumeration(ftSurface *other, double tol,
				std::vector<std::pair<int,int> >& enumeration);

    /// Return the boundary number of outer boundary edges that follow surface
    /// boundaries
    /// return value = true: All free boundaries are associated with surface
    ///                      boundaries
    ///              = false:Some free boundaries are not associated with 
    ///                      surface boundaries
    bool getFreeBoundaryInfo(double tol, std::vector<int>& free_boundaries);

    /// For a spline surface, get the local enumeration of coefficients
    /// along a specified boundary
    /// return value: true if spline, false if not
    bool getBoundaryCoefEnumeration(int bd, std::vector<int>& enumeration);

    /// Get all instances of radial edges
    std::vector<shared_ptr<EdgeVertex> > getRadialEdges() const;

    /// Check for radial edges
    /// Existence
    bool hasRadialEdges() const;
    bool hasRealRadialEdges() const;

    /// All edges is connected a radial edge
    bool allRadialEdges() const;

    /// Fetch all faces that belongs to all radial edges corresponding
    /// to this face
    std::vector<ftSurface*> fetchCorrespondingFaces() const;

    /// Debug functionality
    bool checkFaceTopology();

 protected:
    void replaceSurf(shared_ptr<ParamSurface> sf)
	{ surf_ = sf;}

private:
    /// Geometric description of the surface associated to this face
    shared_ptr<ParamSurface> surf_;

    /// Related to inheritance 
    ftTangPriority prio_type_;

    /// The boundary loops limiting this face, the first loop represents
    /// the initial boundary
    std::vector<shared_ptr<Loop> > boundary_loops_;
    mutable double degenerate_eps_; /// Tolerance used in createInitialEdges
    mutable double kink_;  /// Tolerance used in createInitialEdges

    /// Represents adjacency between bodies
    ftSurface* twin_;  // Intended for use in volume models and other
                       // non-manifold models

    /// A solid associated to this faces, if such an associaction exists.
    /// May point to nothing.
    Body* body_;
//     int id_;
//     bool is_turned_;

    // Private functions
    void 
      getApproxCurves(std::vector<std::pair<shared_ptr<ParamCurve>,
		      shared_ptr<ParamCurve> > >::iterator cvs_in,
		      int nmb_cvs, 
		      std::vector<shared_ptr<SplineCurve> >& cvs_out,
		      double tol);

    void 
      getBoundaryCurves(double kink,
			std::vector<std::pair<shared_ptr<ParamCurve>,
			shared_ptr<ParamCurve> > >& cvs, 
			bool only_corner=false);

};

} // namespace Go



#endif // _FTSURFACE_H

