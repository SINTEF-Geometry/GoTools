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
	Loop(ftFaceBase* face, std::vector<shared_ptr<ftEdgeBase> >& edges, 
	     double space_epsilon);

	/// Constructor that takes an ordered sequence of edges as
	/// input. There is no reference to a face object.
	// @@@jbt - In STEP this corresponds to the entity
	// 'face_outer_bound', which has no reference to "underlying
	// surfaces".
	Loop(std::vector<shared_ptr<ftEdgeBase> >& edges, 
	     double space_epsilon);
	
	/// Destructor
	~Loop();

	/// Number of edges in loop
	size_t size() const
	    {
		return edges_.size();
	    }

	/// Get all edges oriented head to tail
	std::vector<shared_ptr<ftEdgeBase> >& getEdges()
	    {
		return edges_;
	    }

	/// Fetch an edge from the loop specified by the index in the
	/// sequence of edges
	shared_ptr<ftEdgeBase> getEdge(size_t idx)
	    {
		shared_ptr<ftEdgeBase> dummy;
		if (idx < edges_.size())
		    return edges_[idx];
		else
		    return dummy;
	    }

	/// Fetch all vertices between edges in the loop
	std::vector<shared_ptr<Vertex> > getVertices() const;

	/// Fetch all vertices between edges in the loop in sequence
	std::vector<shared_ptr<Vertex> > getSeqVertices() const;

	/// Get tolerance
	double getTol() const
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
	void updateLoop(shared_ptr<ftEdgeBase> new_edge);

	/// Split loop in a given parameter of an edge given by index
	void split(int ind, double par);

	/// Test if an edge is close to the surface within some tolerance
	bool isClose(ftEdge* edge,
		     RectDomain* domain,
		     double tol) const;

	/// Collect all pairs of surface and edges where some part of the edge has a distance to the surface greater than a tolerance
	void getBadDistance(std::vector<std::pair<ftSurface*, ftEdge* > >& badPairs,
			    RectDomain* domain,
			    double tol) const;

	// Collect all pairs of edge and vertex where the vertex has a distance to the edge greater than a tolerance
	void getBadDistance(std::vector<std::pair<ftEdge*, shared_ptr<Vertex> > >& badPairs,
			    double tol) const;

	/// Fetch information on continuity 
	void getPosTangentSurfaceDiscont(std::vector<ftEdge*>& badPos,
					 std::vector<ftEdge*>& badTangent,
					 double tol, double kink, double bend, int leastSurfIndex,
					 shared_ptr<SurfaceModel> sm) const;

	/// Check if the edges of the loop are consistent with the corresponding curves
	/// with regard to orientation
	bool checkConsistency() const;

	/// Check for acute edges in boundary loop
	void getAcuteEdges(std::vector<std::pair<ftEdge*, ftEdge*> >& acute_edges, double angtol) const;

	/// Compute intersections between boundary loops
	void getLoopIntersections(shared_ptr<Loop> loop2, double tol, 
				  std::vector<std::pair<shared_ptr<PointOnEdge>, 
				  shared_ptr<PointOnEdge> > >& int_pt) const;

	/// Compute self intersections of boundary loop
	void getLoopSelfIntersections(double tol, 
				      std::vector<std::pair<shared_ptr<PointOnEdge>, 
				      shared_ptr<PointOnEdge> > >& int_pt) const;

	/// Non-manifold functionality
	/// Check for loop correspondance, split edges if required and
	/// return the start edges for corresponding loops
	bool correspondingEdges(shared_ptr<Loop> other, double tol,
				ftEdgeBase* &first1, ftEdgeBase* &first2,
				bool& same_dir, bool no_snap=true);
	/// Check for radial edges
	/// Existance
	bool hasRadialEdges() const;
	
	/// All edges is connected a radial edge
	bool allRadialEdges() const;

    /// Find the closest point on the curve loop to a point specified 
    /// by the user.
    /// \param pt The point given by the user.  We want to determine the closest
    ///           point to this on the CurveLoop.
    /// \param clo_ind Upon return: the index of the curve segment on which the 
    ///                closest point was found.
    /// \param clo_par Upon return: the parameter of the detected closest point on 
    ///                the curve containing it.
    /// \param clo_pt  Upon return: the geometric position of the detected closest 
    ///                point
    /// \param clo_dist Upon return: the distance to the detected closest point.
    void closestPoint(const Point& pt, int& clo_ind, double& clo_par, 
		      Point& clo_pt, double& clo_dist) const;

    /// Check if a given edge is in this loop
    bool isInLoop(ftEdgeBase* edge)
    {
      for (size_t ki=0; ki<edges_.size(); ++ki)
	if (edges_[ki].get() == edge)
	  return true;
      return false;
    }

    /// Remove specified edge. Make sure that this does not violate the
    /// loop connectivity
    void removeEdge(ftEdgeBase* edge);

    /// Group edges that are smoothly joined together. The sequence of edges
    /// corresponds to the sequence in the Loop, i.e. head to tail connected
    void 
      groupSmoothEdges(double tol, double angtol,
		       std::vector<std::vector<shared_ptr<ftEdgeBase> > >& edge_groups);

    /// Return true if prev_ and next_ corresponds to the ordering of the edge vector.
    bool consistentEdgeTopology();

    private:
	/// The face which the loop belongs to. In cases where the loop
	/// is used to represent a closed sequence of edges in general, 
	/// this pointer is not set
	ftFaceBase* face_;

	/// The edges which make up the loop
	std::vector<shared_ptr<ftEdgeBase> > edges_;

	/// Positional tolerance regarding the distance between
	/// subsequent edges
	double eps_;
	

	void setEdges(CurveLoop& curve_loop, double kink, bool split_in_kinks,
		      bool no_split);

	void setEdges(std::vector<shared_ptr<ftEdgeBase> >& edges);

    };

} // namespace Go



#endif // _LOOP_H
