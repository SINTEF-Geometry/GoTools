//===========================================================================
//                                                                           
// File: tpUtils.h                                                           
//                                                                           
// Created: Wed Jul 10 09:49:17 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: tpUtils.h,v 1.5 2004-01-20 11:24:34 bsp Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _TPUTILS_H
#define _TPUTILS_H

// #include "edgeType.h"
#include "GoTools/topology/tpJointType.h"

namespace Go
{

namespace tpUtils
{

    /// Consider the 'node' or 'corner' implied by the startpoint
    /// (if at_start_of_edge is true) or endpoint (if it is false).
    /// The vector adjacent is filled with all edges with that 'node'
    /// as a startpoint (indicated by a true at the corresponding place
    /// in the at_start vector) or endpoint (false in at_start).
    /// The edge originally asked for is not included.
    //=======================================================================
    template <class edgeType>
    void adjacentEdges(edgeType* edge, bool at_start_of_edge,
		       std::vector<edgeType*>& adjacent,
		       std::vector<bool>& at_start)
    //===========================================================================
    {
	// Clear the input vectors
	adjacent.clear();
	at_start.clear();

	// Every edge in a valid structure must have a next and prev edge.
	// We use this to get an edge pointing out from the edge-point in
	// question.
	edgeType* orige = 0;
	if (!at_start_of_edge) {
	    orige = edge->next();
	} else {
	    orige = edge;
	}

	edgeType* e = orige;

	adjacent.push_back(e);
	at_start.push_back(true);

	// Iterate clockwise until we meet a boundary or we've come back to
	// this edge.
	bool finished = false;
	bool full_circle = false;
	while (!finished) {
	    edgeType* t = e->twin();
	    if (t == 0) {
		finished = true;
	    } else {
		adjacent.push_back(t);
		at_start.push_back(false);
		e = t->next();
		if (e == orige) {
		    finished = true;
		    full_circle = true;
		} else {
		    adjacent.push_back(e);
		    at_start.push_back(true);
		}
	    }
	}

	// If we met a boundary, we iterate counterclockwise until we meet
	// a boundary.
	e = orige;
	if (!full_circle) {
	    finished = false;
	    while (!finished) {
		edgeType* p = e->prev();
		adjacent.push_back(p);
		at_start.push_back(false);
		e = p->twin();
		if (e == 0) {
		    finished = true;
		} else {
		    adjacent.push_back(e);
		    at_start.push_back(true);
		}
	    }
	}

	// We have to remove this from the vector, because we promised
	// to do so in the doc...
	typename std::vector<edgeType*>::iterator it
	    = std::find(adjacent.begin(), adjacent.end(), edge);
	adjacent.erase(it);
	at_start.erase(at_start.begin() + (it-adjacent.begin()));
    }


    /// Return continuity between end of edge1 and start of edge2.
    //=======================================================================
    template <class edgeType>
    tpJointType checkContinuity(const edgeType* edge1, const edgeType* edge2,
				double neighbour, double gap,
				double bend, double kink)
    //=======================================================================
    {
	// Get endparameters of segments
	double tmax1 = edge1->tMax();
	double tmin2 = edge2->tMin();

	// Get endpoints
	Go::Point p1 = edge1->point(tmax1);
	Go::Point p1_0 = edge1->point(edge1->tMin());
	Go::Point p2 = edge2->point(tmin2);
	Go::Point p2_1 = edge2->point(edge2->tMax());
	//    std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
	//    std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
	//    std::cout << p1_0[0] << " " << p1_0[1] << " " << p1_0[2] << std::endl;
	//    std::cout << p2_1[0] << " " << p2_1[1] << " " << p2_1[2] << std::endl;

	// Compute distance between endpoints
	double dist = p1.dist(p2);
	//  std::cout << "Distance: " << dist << std::endl;
	if (dist > neighbour)
	    return JOINT_DISC;

	if (dist > gap)
	    return JOINT_GAP;

	// Compute tangents at endpoints
	Go::Point t1 = edge1->tangent(tmax1);
	Go::Point t2 = edge2->tangent(tmin2);
	double ang_dist = t1.angle(t2);

	if (ang_dist > bend)
	    return JOINT_G0;

	if (ang_dist > kink)
	    return JOINT_KINK;

	return JOINT_G1;
    }


}

} // namespace Go
#endif // _TPUTILS_H

