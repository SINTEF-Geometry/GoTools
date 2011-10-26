//===========================================================================
//                                                                           
// File: Loop.h                                                       
//                                                                           
// Created: August 26th 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LOOP_H
#define _LOOP_H

#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include <vector>

namespace Go
{

//===========================================================================
/** Loop. Boundary loop connected to a face (see ftSurface)
 * Detailed description.
 *
 * \author Vibeke Skytt
 * \see ftFaceBase
*/
//===========================================================================

    class ftEdgeBase;
    class ftFaceBase;
    class PointOnEdge;

    /// \brief Primarily a loop connected to a face in a boundary represented
    /// solid or face set. May also be used to represent a general closed
    /// sequence of edges

class Loop
    {
    public:
	/// Constructor
      Loop(ftFaceBase* face, CurveLoop& curve_loop, double kink,
	   bool split_in_kinks = true, bool no_split = false);

	/// This constructor takes an ordered sequence of edges as input
	/// Note that the function may throw
	Loop(ftFaceBase* face, std::vector<std::shared_ptr<ftEdgeBase> >& edges, 
	     double space_epsilon);

	/// Constructor that takes an ordered sequence of edges as
	/// input. There is no reference to a face object.
	// @@@jbt - In STEP this corresponds to the entity
	// 'face_outer_bound', which has no reference to "underlying
	// surfaces".
	Loop(std::vector<std::shared_ptr<ftEdgeBase> >& edges, 
	     double space_epsilon);
	
	/// Destructor
	~Loop();

	/// Number of edges in loop
	size_t size()
	    {
		return edges_.size();
	    }

	/// Get all edges oriented head to tail
	std::vector<std::shared_ptr<ftEdgeBase> >& getEdges()
	    {
		return edges_;
	    }

	/// Fetch an edge from the loop specified by the index in the
	/// sequence of edges
	std::shared_ptr<ftEdgeBase> getEdge(size_t idx)
	    {
		std::shared_ptr<ftEdgeBase> dummy;
		if (idx < edges_.size())
		    return edges_[idx];
		else
		    return dummy;
	    }

	/// Fetch all vertices between edges in the loop
	std::vector<std::shared_ptr<Vertex> > getVertices() const;

	/// Get tolerance
	double getTol()
	    {
		return eps_;
	    }

	/// Check consistency with regard to face
	bool isFaceConsistent();

	/// Set face pointer
	void setFace(ftFaceBase* face);


	/// The face associated to this loop
	const ftFaceBase* getFace()
	    {
		return face_;
	    }

	/// Make sure that the edges belonging to this loop is complete
	void updateLoop(std::shared_ptr<ftEdgeBase> new_edge);

	/// Test if an edge is close to the surface within some tolerance
	bool isClose(ftEdge* edge,
		     RectDomain* domain,
		     double tol) const;

	/// Collect all pairs of surface and edges where some part of the edge has a distance to the surface greater than a tolerance
	void getBadDistance(std::vector<std::pair<ftSurface*, ftEdge* > >& badPairs,
			    RectDomain* domain,
			    double tol) const;

	// Collect all pairs of edge and vertex where the vertex has a distance to the edge greater than a tolerance
	void getBadDistance(std::vector<std::pair<ftEdge*, std::shared_ptr<Vertex> > >& badPairs,
			    double tol) const;

	/// Fetch information on continuity 
	void getPosTangentSurfaceDiscont(std::vector<ftEdge*>& badPos,
					 std::vector<ftEdge*>& badTangent,
					 double tol, double kink, double bend, int leastSurfIndex,
					 std::shared_ptr<SurfaceModel> sm) const;

	/// Check if the edges of the loop are consistent with the corresponding curves
	/// with regard to orientation
	bool checkConsistency() const;

	/// Check for acute edges in boundary loop
	void getAcuteEdges(std::vector<std::pair<ftEdge*, ftEdge*> >& acute_edges, double angtol) const;

	/// Compute intersections between boundary loops
	void getLoopIntersections(std::shared_ptr<Loop> loop2, double tol, 
				  std::vector<std::pair<std::shared_ptr<PointOnEdge>, 
				  std::shared_ptr<PointOnEdge> > >& int_pt) const;

	/// Compute self intersections of boundary loop
	void getLoopSelfIntersections(double tol, 
				      std::vector<std::pair<std::shared_ptr<PointOnEdge>, 
				      std::shared_ptr<PointOnEdge> > >& int_pt) const;

	/// Non-manifold functionality
	/// Check for loop correspondance, split edges if required and
	/// return the start edges for corresponding loops
	bool correspondingEdges(std::shared_ptr<Loop> other, double tol,
				ftEdgeBase* &first1, ftEdgeBase* &first2,
				bool& same_dir, bool no_snap=true);
	/// Check for radial edges
	/// Existance
	bool hasRadialEdges() const;
	
	/// All edges is connected a radial edge
	bool allRadialEdges() const;

    private:
	/// The face which the loop belongs to. In cases where the loop
	/// is used to represent a closed sequence of edges in general, 
	/// this pointer is not set
	ftFaceBase* face_;

	/// The edges which make up the loop
	std::vector<std::shared_ptr<ftEdgeBase> > edges_;

	/// Positional tolerance regarding the distance between
	/// subsequent edges
	double eps_;
	

	void setEdges(CurveLoop& curve_loop, double kink, bool split_in_kinks,
		      bool no_split);

	void setEdges(std::vector<std::shared_ptr<ftEdgeBase> >& edges);

    };

} // namespace Go



#endif // _LOOP_H
