//===========================================================================
//                                                                           
// File: tpTopologyTable.h                                                   
//                                                                           
// Created: Tue Mar 21 15:25:28 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: tpTopologyTable.h,v 1.57 2009-03-23 10:25:06 vsk Exp $
//                                                                           
// Description: Topology container and some related classes/structs
//
// Implementation files: tpTopologyTable.C, topanalysis.C
//                                                                           
//===========================================================================

#ifndef _TPTOPOLOGYTABLE_H
#define _TPTOPOLOGYTABLE_H

#include "GoTools/topology/tpTolerances.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/CurveOnSurface.h"

#include <vector>
#include <set>
#include <memory>
#include <fstream>


using std::vector;
// #include "edgeType.h"
// class faceType;

/** tpTopologicalInfo -  Short description.
 * Detailed description.
 */
struct tpTopologicalInfo
{
public:

    // The status is:
    // 0 : edges join smoothly. G1.
    // 1 : edges join, but normals are slightly discontinous. A kink.
    // 2 : edges join, but the normals are discontinous. G0.
    // 3 : edges almost join. A gap.
    // 4 : edges are totally discontinous.
    // The minimal tpTopologicalInfo has a one-element status vector
    // and a two-element parameters vector
    std::vector<int> status_;
    std::vector< std::pair<double, double> > parameters_;

    int BestStatus() const
    {
	int s = 4;
	for (size_t i = 0; i < status_.size(); ++i)
	    if (status_[i] < s)
		s = status_[i];
	return s;
    }
    int WorstStatus() const
    {
	int s = 0;
	for (size_t i = 0; i < status_.size(); ++i)
	    if (status_[i] > s) s = status_[i];
	return s;
    }
};

/** tpTableEntry -  Short description.
 * Detailed description.
 */
template <class edgeType>
class tpTableEntry
{
public:
    edgeType* e1_;
    edgeType* e2_;
    tpTopologicalInfo info_;
    tpTableEntry(edgeType* e1,
		 edgeType* e2,
		 const tpTopologicalInfo& info)
	: e1_(e1), e2_(e2), info_(info)
    {}
    const tpTopologicalInfo& Info() { return info_; }
}; // End of tpTableEntry protected class


/// Helper class for marching
class tpMarchPoint
{
public:
    Go::Point pt;
    double par;
    double par2;
    double dist;
    double cos_ang;
    int status;
    tpMarchPoint(const Go::Point& p,  double pa, double pa2,
		 double d, double ca, int s)
	: pt(p), par(pa), par2(pa2),  dist(d), cos_ang(ca), status(s)
    {}
    bool operator < (const tpMarchPoint& other) const
    { return par < other.par; }
};

    

// ------------ The main class definition starts here ------------

//===========================================================================
/** tpTopologyTable - Encapsulates a topology structure.
 * Detailed description.
 *
 * This class encapsulates a topology structure. It points to external (to
 * this class) surface objects, which are gotten as input to the
 * PrepareTable() and ConstructTable() members. It is implemented in terms
 * of edgeType objects in a half-edge data structure, and many members
 * operate on and/or return pointers to the edgeType objects owned by the
 * tpTopologyTable object.
 *
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \see edgeType
 * \see faceType
 */
//===========================================================================

// edgeType and faceType encapsulate structures with certain given operations.
// As of now the structures edgeType and faceType fullfills these requirements.
template <class edgeType, class faceType>
class tpTopologyTable
{
protected:

    std::vector< tpTableEntry<edgeType> > table_;
/*  VSK, 25.09.08. The lookup table is not in use an makes it more complex to add
    and remove faces from the topology table. */
/*     std::vector< std::vector<int> > lookup_; */
    std::vector< shared_ptr<edgeType> > edges_;
    tpTolerances tol_;

    std::vector<std::pair<faceType*, faceType* > > orientation_inconsist_;

public:

    /** Constructor.
     * Detailed description.
     */
    tpTopologyTable(double tol_gap, double tol_neighbour,
		    double tol_kink, double tol_bend)
	: tol_(tpTolerances(tol_gap, tol_neighbour, tol_kink, tol_bend))
    {}


    /** Destructor.
     * Detailed description.
     */
    ~tpTopologyTable()
    {
    }


    /// Changes the tolerances used in constructing the table.
    //=======================================================================
    void setTolerances(const tpTolerances& tol)
    //=======================================================================
    {
	tol_ = tol;
    }


    tpTolerances getTolerances() 
    {
      return tol_;
    }

    // For manually adding entries to table
    //=======================================================================
    void prepareTable(const std::vector<shared_ptr<faceType> >& faces)
    //=======================================================================
    {
	initiateTable(faces);

    }


    //=======================================================================
    void constructTable(const std::vector<shared_ptr<faceType> >& faces,
			bool test_orientation = false)
    //=======================================================================
    {
	int i, j, k, l;
	int num_faces = faces.size();
	std::vector<Go::BoundingBox> boxes;
	boxes.reserve(num_faces);
	for (i = 0; i < num_faces; ++i) {
	    boxes.push_back(faces[i]->boundingBox());
	}

	std::vector<shared_ptr<edgeType> > startedges0, startedges1;
	for (i = 0; i < num_faces - 1; ++i) {
	    for (j = i + 1; j < num_faces; ++j) {
		// For every combination of faces, do a boxtest.
		if (boxes[i].overlaps(boxes[j], tol_.neighbour)) {
		    // We have some possible neighbourhood incidents.
		    // Now do a box test on every combination of edges
		    startedges0 = faces[i]->startEdges();
		    startedges1 = faces[j]->startEdges();
		    // Testing all loops in one surface against
		    // all loops in the other.
		    for (k = 0; k < int(startedges0.size()); ++k) {
			for (l = 0; l < int(startedges1.size()); ++l) {
			    edgeType* s0 = startedges0[k].get();
			    edgeType* s1 = startedges1[l].get();
			    if (s0 ==0 || s1 == 0) break;
			    edgeType* e[2];
			    e[0] = s0;
			    e[1] = s1;
			    edgeType* en[2];
			    bool finished = false;
			    while(!finished) {
				en[0] = e[0]->next();
				en[1] = e[1]->next();
#ifdef TOPOLOGY_DEBUG
				std::ofstream debug("data/debug.g2");
				for (int ki = 0; ki < 2; ++ki) {
				    e[ki]->face()->surface()->writeStandardHeader(debug);
				    e[ki]->face()->surface()->write(debug);
				    std::vector<double> pts(12);
				    Go::Point from = e[ki]->point(e[ki]->tMin());
				    double tmid = 0.5*(e[ki]->tMin() + e[ki]->tMax());
				    Go::Point mid = e[ki]->point(tmid);
				    Go::Point to = e[ki]->point(e[ki]->tMax());
				    std::copy(from.begin(), from.end(), pts.begin());
				    std::copy(mid.begin(), mid.end(), pts.begin() + 3);
				    std::copy(mid.begin(), mid.end(), pts.begin() + 6);
				    std::copy(to.begin(), to.end(), pts.begin() + 9);
				    Go::LineCloud lc(pts.begin(), 2);
				    lc.writeStandardHeader(debug);
				    lc.write(debug);
				}
#endif // TOPOLOGY_DEBUG
				if (e[0]->boundingBox().overlaps(e[1]->boundingBox(),
								 tol_.neighbour)) {
				    // We found an edge overlap. Possible incident.
				    int incident_occurred = 
					testEdges(e);
				    if (incident_occurred) {
					// We skip the rest of this subloop (looping
					// over edges e[1] in face faces[j]) by
					// making en[1] so that e[0] will be
					// incremented.
					// If e[0] was split w/t-value higher than
					// start value, do not forget first part of
					// edge.

					if (incident_occurred >= 2)
					{
					    // Inconsistence in face orientation
					    // Remember incident
					    // Check if it has occured before
					    size_t kr;
					    for (kr=0; kr<orientation_inconsist_.size(); ++kr)
						if ((orientation_inconsist_[kr].first == faces[i].get() &&
						     orientation_inconsist_[kr].second == faces[j].get()) ||
						    (orientation_inconsist_[kr].first == faces[j].get() &&
						     orientation_inconsist_[kr].second == faces[i].get()))
						    break;

					    if (kr == orientation_inconsist_.size())
					      orientation_inconsist_.push_back(std::make_pair(faces[i].get(),
											   faces[j].get()));
					}
				    }
				}
				// Pick next edges, check if we're done
				e[1] = en[1];
				if (e[1] == s1) {
				    e[0] = en[0];
				    if (e[0] == s0)
					finished = true;
				}
			    }
			}
		    }
		}
	    }
	}
	if (!test_orientation) {
	    Validate();
	}

    }


    /// In the AddEntry() function, the min? and max? arguments
    /// are the parameter values associated with the endpoints
    /// of the curve segments. min1 <= max1 and min2 <= max2.
    /// BUT the edges may very well run in different directions,
    /// and the point on the first edge corresponding to
    /// min1 may be adjacent to either the min2- or the max2-
    /// corresponding point on the second edge.
    //=======================================================================
    void addEntry(edgeType* e1, edgeType* e2,
		  double min1, double max1,
		  double min2, double max2,
		  const tpTopologicalInfo& info,
		  int status)
    //=======================================================================
    {
/*  VSK, 25.09.08. The lookup table is not in use an makes it more complex to add
    and remove faces from the topology table. */
/* 	faceType* f1 = e1->face(); */
/* 	faceType* f2 = e2->face(); */
	int index_in_table = table_.size();
/* 	int id_of_face_1 = f1->getId(); */
/* 	int id_of_face_2 = f2->getId(); */
/* 	// Here I assume IDs are from 0 up to a reasonable number (such as */
/* 	// the number of faces minus one...). */
/* 	int lookup_size_needed = std::max(id_of_face_1, id_of_face_2) + 1; */
/* 	lookup_.resize(lookup_size_needed); */
/* 	lookup_[id_of_face_1].push_back(index_in_table); */
/* 	lookup_[id_of_face_2].push_back(index_in_table); */

	// Now, we create the actual twin-edges. We assume that any necessary
	// tolerance tests have taken place already (i.e. outside this func.),
	// so that if the min?/max? arguments *should* be equal to high_param_
	// or low_param_ of the edges involved, then they *are* equal.

	std::vector<edgeType*> twins(2);
	twins = connectTwins(e1, e2, min1, max1, min2, max2, status);
	for (int i = 0; i < 2; ++i)
	    twins[i]->setEntryId(index_in_table);

	if (status != 2)
	    table_.push_back(tpTableEntry<edgeType>(twins[0], twins[1], info));
    }


    /// Returns pointer to (vector of) twin edges.
    /// Updates topological structure and, if necessary, splits edges.
    //=======================================================================
    std::vector<edgeType*> connectTwins(edgeType* e1, edgeType* e2,
				  double min1, double max1,
				  double min2, double max2,
				   int& status)
    //=======================================================================
    {
	double lp[2];
	lp[0] = e1->tMin(); lp[1] = e2->tMin();
	double hp[2];
	hp[0] = e1->tMax(); hp[1] = e2->tMax();
	edgeType* ed[2];
	ed[0] = e1; ed[1] = e2;
	double mins[2];
	mins[0] = min1; mins[1] = min2;
	double maxs[2];
	maxs[0] = max1; maxs[1] = max2;
	std::vector<edgeType*> twins(2);
	edgeType* newedge;
	//    int index_in_table = table_.size();

	for (int i = 0; i < 2; ++i) {
	    twins[i] = ed[i];
	    if (maxs[i] < hp[i]) {
		newedge = twins[i]->split(maxs[i]);
		shared_ptr<edgeType> newtwin = shared_ptr<edgeType>(newedge);
		newtwin->face()->updateBoundaryLoops(newtwin);
		edges_.push_back(newtwin);
		// ed[i] = twins[i]; // We need to update ed[i] with the new limits
	    }
	    // else {
	    // ...
	    //    twins[i] = ed[i];
	    //	}
	    if (mins[i] > lp[i]) {
		newedge = ed[i];
		twins[i] = newedge->split(mins[i]);
		shared_ptr<edgeType> newtwin = shared_ptr<edgeType>(twins[i]);
		newtwin->face()->updateBoundaryLoops(newtwin);
		edges_.push_back(shared_ptr<edgeType>(newtwin));
	    }
	    // else {
	    //    twins[i] = ed[i];
	    //	}
	    //	twins[i]->SetEntryId(index_in_table);
	}

	twins[0]->connectTwin(twins[1], status); // Connects both twins to each other

	//    e1 = ed[0]; e2 = ed[1];
	return twins;
    }


    /// Returns all boundary loops. Keep in mind that
    /// calling next() on elements on the boundary will
    /// not iterate around the boundary, but around the face
    /// the edge comes from.
    //=======================================================================
    void BoundaryLoops(std::vector< std::vector<edgeType*> > & loopvec)
    //=======================================================================
    {
	loopvec.clear();
	std::set<edgeType* > used_edges;

	//    std::ofstream dump("debugdump1.g2");
	int loop = 0;
	for (size_t i = 0; i < edges_.size(); ++i) {
	    edgeType* edge = edges_[i].get();
	    bool in_set = (used_edges.find(edge) != used_edges.end());
	    if (edges_[i]->twin() == 0 && !in_set) {
		loopvec.push_back(std::vector<edgeType*>());
		while (!in_set) {
		    //  	  ftEdge *geomedge = edge->geomEdge();
		    //            geomedge->SpaceCurve()->writeStandardHeader(dump);
		    //  	      geomedge->SpaceCurve()->write(dump);
		    used_edges.insert(edge);
		    loopvec[loop].push_back(edge);
		    edge = edge->next();
		    std::set<edgeType* > tmp_edges; // To avoid infinite loop in case
		    tmp_edges.insert(edge);
		    // of inconsistency
		    while (edge->twin())
		      {
			edge = edge->twin()->next();
			if (tmp_edges.find(edge) != tmp_edges.end())
			  break;
			tmp_edges.insert(edge);
		      }
		    in_set = (used_edges.find(edge) != used_edges.end());
		}
		++loop;
	    }
	}
    }


    //=======================================================================
    void BoundaryLoops(std::vector< std::vector<shared_ptr<edgeType> > > & loopvec)
    //=======================================================================
    {
	loopvec.clear();
	std::set<edgeType* > used_edges;
	edgeType* e2;

	//    std::ofstream dump("debugdump1.g2");
	int loop = 0;
	for (size_t i = 0; i < edges_.size(); ++i) {
	    shared_ptr<edgeType> edge = edges_[i];
	    bool in_set = (used_edges.find(edge.get()) != used_edges.end());
	    if (edges_[i]->twin() == 0 && !in_set) {
		loopvec.push_back(std::vector<shared_ptr<edgeType> >());
		while (!in_set) {
		    //  	  ftEdge *geomedge = edge->geomEdge();
		    //            geomedge->SpaceCurve()->writeStandardHeader(dump);
		    //  	      geomedge->SpaceCurve()->write(dump);
		  used_edges.insert(edge.get());
		    loopvec[loop].push_back(edge);
		    e2 = edge->next();
		    for (size_t kr=0; kr<edges_.size(); ++kr)
		      if (edges_[kr].get() == e2)
			{
			  edge = edges_[kr];
			  break;
			}
		    std::set<edgeType* > tmp_edges; // To avoid infinite loop in case
		    tmp_edges.insert(edge.get());
		    // of inconsistency
		    while (edge->twin())
		      {
			e2 = edge->twin()->next();
			for (size_t kr=0; kr<edges_.size(); ++kr)
			  if (edges_[kr].get() == e2)
			    {
			      edge = edges_[kr];
			      break;
			    }
			if (tmp_edges.find(e2) != tmp_edges.end())
			  break;
			tmp_edges.insert(e2);
		      }
		    in_set = (used_edges.find(e2) != used_edges.end());
		}
		++loop;
	    }
	}
    }


    // Routine separates objects which are not path connected.
    //=======================================================================
    void disjointObjects(std::vector<std::vector<faceType*> >& grouped_faces)
    //=======================================================================
    {
	grouped_faces.clear();
	std::vector<faceType*> temp_faces(2);
	size_t i, j, k;

	if (table_.size() != 0) {
	    // We search through table_.
	    temp_faces[0] = table_[0].e1_->face();
	    temp_faces[1] = table_[0].e2_->face();
	    grouped_faces.push_back(temp_faces);

	    for (i = 1; i < table_.size(); ++i) {
		temp_faces[0] = table_[i].e1_->face();
		temp_faces[1] = table_[i].e2_->face();
		int found1 = -1, found2 = -1;
		for (j = 0; j < grouped_faces.size(); ++j)
		    for (k = 0; k < grouped_faces[j].size(); ++k) {
			if (grouped_faces[j][k] == temp_faces[0])
			    found1 = j;
			if (grouped_faces[j][k] == temp_faces[1])
			    found2 = j;
		    }

		if ((found1 == -1) && (found2 == -1))
		    grouped_faces.push_back(temp_faces);
		else if ((found1 != -1) && (found2 != -1)) {
		    if (found1 != found2) {
			grouped_faces[found1].insert(grouped_faces[found1].end(),
						     grouped_faces[found2].begin(),
						     grouped_faces[found2].end());
			grouped_faces.erase(grouped_faces.begin() + found2,
					    grouped_faces.begin() + found2 + 1);
		    }
		} else if (found1 != -1)
		    grouped_faces[found1].push_back(temp_faces[1]);
		else
		    grouped_faces[found2].push_back(temp_faces[0]);
	    }
	}

	// We have pushed back faces for all connected edges, on to the twinless!
	temp_faces.resize(1);
	for (i = 0; i < edges_.size(); ++i)
	    if (edges_[i]->twin() == 0) {
		for (j = 0; j < grouped_faces.size(); ++j) {
		    for (k = 0; k < grouped_faces[j].size(); ++k)
			if (grouped_faces[j][k] == edges_[i]->face())
			    break;
		    if (k < grouped_faces[j].size())
			break;
		}
		if (j == grouped_faces.size()) {
		    temp_faces[0] = edges_[i]->face();
		    grouped_faces.push_back(temp_faces);
		}
	    }

	return;
    }


    /// Only gives the first of every pair of edges
    /// representing a cornering edge or kink edge.
    //=======================================================================
    void cornersAndKinks(std::vector<edgeType*>& vec)
    //=======================================================================
    {
	vec.clear();
	for (size_t i = 0; i < table_.size(); ++i) {
	    if (table_[i].Info().WorstStatus() > 0) {
		vec.push_back(table_[i].e1_); // Always e1_, e2_ is ignored...
	    }
	}
    }


    /// Returns the raw topological info element.
    //=======================================================================
    const tpTopologicalInfo& info(int n)
    //=======================================================================
    {
	DEBUG_ERROR_IF(n<0 || n>=int(table_.size()), "Info index out of bounds");
	return table_[n].Info();
    }


    /// Returns all edges bounding a specific face
    //=======================================================================
    std::vector<edgeType*> edgesBoundingFace(faceType* face) const
    //=======================================================================
    {
	std::vector<edgeType*> bounding_edges;
	edgeType* e = 0;
	for (size_t i = 0; i < edges_.size(); ++i) {
	    if (edges_[i]->face() == face) {
		e = edges_[i].get();
		break;
	    }
	}
	if (e == 0)
	    return bounding_edges;
	bounding_edges.push_back(e);
	edgeType* orig = e->next();
	while (e != orig) {
	    bounding_edges.push_back(e);
	    e = e->next();
	}
	return bounding_edges;
    }


    /// Returns a const reference to a vector representing twin edges.
    //=======================================================================
    const std::vector<tpTableEntry<edgeType> >& connectedEdges() const
    //=======================================================================
    {
	return table_;
    }

    /// Returns information about face pair with an orientation inconsistency
    //=======================================================================
    const void getInconsistentFaces(std::vector<std::pair<faceType*, faceType* > >& faces)
    //=======================================================================
    {
	faces = orientation_inconsist_;
    }

    // For debugging
    //=======================================================================
    void Validate()
    //=======================================================================
    {
	// First, check that all edges have a nonzero parameter interval
	int n = edges_.size();
	for (int i = 0; i < n; ++i)
	    if (edges_[i]->tMax() <= edges_[i]->tMin())
		std::cout << "Edge " << i << " has zero parameter interval." <<
		    std::endl;
	// Check all table entries for same
	int nt = table_.size();
	for (int i = 0; i < nt; ++i) {
	    const tpTopologicalInfo& info = table_[i].info_;
	    for (size_t j = 0; j < info.status_.size(); ++j)
		if (info.parameters_[j].first >= info.parameters_[j+1].first)
		    std::cout << "Error in entry " << i << std::endl;
	}
    }


    //=======================================================================
    // Fetch existing adjacency information between faces and build a 
    // topology table representing this adjacency.    
    void setTable(const std::vector<shared_ptr<faceType> >& faces)
//=======================================================================
    {
	initiateTable(faces);

	size_t ki, kj;
	edgeType *e1, *e2;
	for (ki=0; ki<edges_.size(); ki++)
	{
	    // Check if the edge has a twin
	    e1 = edges_[ki].get();
	    e2 = e1->twin();
	    if (e2)
	    {
		// Check if the entry is set already
		for (kj=0; kj<table_.size(); kj++)
		{
		    if ((table_[kj].e1_ == e1 && table_[kj].e2_ == e2) ||
			(table_[kj].e1_ == e2 && table_[kj].e2_ == e1))
			break;
		}

		if (kj < table_.size())
		    continue;  // Table entry already computed

		tpTopologicalInfo topinfo;  // Dummy information
		addEntry(e1, e2, e1->tMin(), e1->tMax(), e2->tMin(), e2->tMax(), topinfo, 1);

		// Update topological information
		updateTableEntry(e1, e2);
	    }
	}
    }

    
    //=======================================================================
    // Remove one face from the current table
    void removeFaceFromTable(shared_ptr<faceType> face)
    //=======================================================================
    {
      // Check if the face is listed as inconsistent. In that case remove 
      // table entry
      size_t ki;
      for (ki=0; ki<orientation_inconsist_.size();)
	{
	  if (orientation_inconsist_[ki].first == face.get() ||
	      orientation_inconsist_[ki].second == face.get())
	      orientation_inconsist_.erase(orientation_inconsist_.begin() + ki);
	  else
	    ki++;
	}

      // Fetch face edges
      std::vector<shared_ptr<edgeType> > tmp_edges;
      tmp_edges = face->createInitialEdges(tol_.neighbour, tol_.kink);

      // For each edge, remove edges from edge list and entries from the table
      for (ki=0; ki<tmp_edges.size(); ki++)
	{
	  edgeType *curr = tmp_edges[ki].get();
	  size_t kj;
	  for (kj=0; kj<edges_.size(); kj++)
	    if (edges_[kj].get() == curr)
	      {
		// Remove entry
		edges_.erase(edges_.begin() + kj);
		break;
	      }

	  for (kj=0; kj<table_.size(); kj++)
	    if (table_[kj].e1_ == curr || table_[kj].e2_ == curr)
	      {
		// Remove twin information
		//table_[kj].e1_->disconnectTwin(table_[kj].e2_);
		curr->disconnectTwin();
		table_.erase(table_.begin() + kj);
		kj--;
	      }
	}

      face->isolateFace();
      int stop_here;
      stop_here = 1;
    }


    
    //=======================================================================
// Add one face to the current table
    void addFaceToTable(shared_ptr<faceType> face)
    //=======================================================================
	{
	    // Fetch face edges
	    std::vector<shared_ptr<edgeType> > tmp_edges;
	    tmp_edges = face->createInitialEdges(tol_.neighbour);

	    // Add new edges to list
	    size_t nmb0 = edges_.size();
	    edges_.insert(edges_.end(), tmp_edges.begin(), tmp_edges.end());

	    // Store edge boxes to avoid multiple computation
	    size_t ki, kj, kr;
	    size_t nmb1 = edges_.size();
	    std::vector<Go::BoundingBox> boxes;
	    boxes.reserve((int)nmb1);
	    for (ki=0; ki<nmb1; ki++)
		boxes.push_back(edges_[ki]->boundingBox());

#ifdef DEBUG
	    std::ofstream of("top.txt");
#endif

	    edgeType* e[2];
	    bool split1 = false, split2 = false;
	    for (ki=nmb0; ki<edges_.size(); )
	      {
		split1 = false;
		for (kj=0; kj<edges_.size(); )
		  {
		    split2 = false;
		    if (ki == kj)
		      {
			kj++;
			continue;  // Same edge
		      }

		    if (boxes[kj].overlaps(boxes[ki], tol_.neighbour))
		      {
			// We have some possible neighbourhood incidents.
			e[0] = edges_[kj].get();
			e[1] = edges_[ki].get();

#ifdef DEBUG
			of << kj << "; " << e[0] << ": [" << e[0]->tMin() << ",";
			of << e[0]->tMax() << "]  ";
			of << e[0]->point(e[0]->tMin()) << ", ";
			of << e[0]->point(e[0]->tMax()) << std::endl;
			of << ki << "; " << e[1] << ": [" << e[1]->tMin() << ",";
			of << e[1]->tMax() << "]  ";
			of << e[1]->point(e[1]->tMin()) << ", ";
			of << e[1]->point(e[1]->tMax()) << std::endl;
#endif
			int status = testEdges(e);

#ifdef DEBUG
			of << "Status: " << status << std::endl;
			if (status > 0)
			  {
			    of << e[0] << ": [" << e[0]->tMin() << ",";
			    of << e[0]->tMax() << "]  ";
			    of << e[0]->point(e[0]->tMin()) << ", ";
			    of << e[0]->point(e[0]->tMax()) << std::endl;
			    of << e[1] << ": [" << e[1]->tMin() << ",";
			    of << e[1]->tMax() << "]  ";
			    of << e[1]->point(e[1]->tMin()) << ", ";
			    of << e[1]->point(e[1]->tMax()) << std::endl;
			  }
			of << std::endl;
#endif
			if (status >= 2)
			  {
			    // Inconsistence in face orientation
			    // Remember incident
			    // Check if it has occured before
			    size_t kr;
			    for (kr=0; kr<orientation_inconsist_.size(); ++kr)
			      if ((orientation_inconsist_[kr].first == e[0]->face() &&
				   orientation_inconsist_[kr].second == e[1]->face()) ||
				  (orientation_inconsist_[kr].first == e[1]->face() &&
				   orientation_inconsist_[kr].second == e[0]->face()))
				break;

			    if (kr == orientation_inconsist_.size())
			      orientation_inconsist_.push_back(make_pair(e[0]->face(),e[1]->face()));
			  }

			if (edges_.size() > nmb1)
			  {
			    // Some edge is split. Make box
			    for (kr=nmb1; kr<edges_.size(); kr++)
			      boxes.push_back(edges_[kr]->boundingBox());
			    nmb1 = edges_.size();

			    split1 = true;
			    split2 = true;
			  }
		      }
		    if (!split2)
		      kj++;
		  }
		if (!split1)
		  ki++;
	      }
	}


//=======================================================================
    // Update topological information in a given table entry
    void updateTableEntry(edgeType* e1, edgeType *e2)
    //=======================================================================
	{
	    for (size_t ki=0; ki<table_.size(); ki++)
	    {
		if (!((table_[ki].e1_ == e1 && table_[ki].e2_ == e2) ||
		      (table_[ki].e1_ == e2 && table_[ki].e2_ == e1)))
		    continue;   // Not this entry

		// Compute new topological information
		tpTopologicalInfo topinfo;
		edgeType* curr_edges[2];
		curr_edges[0] = table_[ki].e1_;
		curr_edges[1] = table_[ki].e2_;
		double param0[2], param1[2];
		param0[0] = curr_edges[0]->tMin();
		param0[1] = curr_edges[0]->tMax();
		Go::Point p1, p3, p4;
		p1 = curr_edges[0]->point(curr_edges[0]->tMin());
		p3 = curr_edges[1]->point(curr_edges[1]->tMin());
		p4 = curr_edges[1]->point(curr_edges[1]->tMax());
		if (p1.dist(p3) < p1.dist(p4))
		{
		    // Same orientation of edge
		    param1[0] = curr_edges[1]->tMin();
		    param1[1] = curr_edges[1]->tMax();
		}
		else
		{
		    // Opposite orientation
		    param1[0] = curr_edges[1]->tMax();
		    param1[1] = curr_edges[1]->tMin();
		}

		int status;
		status = march2(curr_edges, param0, param1, topinfo);
		table_[ki].info_ = topinfo;
	    }
		
	}

    
private:

    //=======================================================================
    void initiateTable(const std::vector<shared_ptr<faceType> >& faces)
    //=======================================================================
    {
	// Make the edge vector from the faces
	int n = faces.size();
	edges_.clear();
	//    edges_.reserve(n*4);
	int i;
	int n2 = 0;
	for (i=0; i<n; ++i) {
	    faces[i]->setId(i);
	    std::vector<shared_ptr<edgeType> > tmp_edges;
	    tmp_edges = faces[i]->createInitialEdges(tol_.neighbour);
	    n2 += tmp_edges.size();
	    edges_.insert(edges_.end(), tmp_edges.begin(), tmp_edges.end());
	}
	//      std::cout << "n*4 : " << n*4 << ", n2 : " << n2;
	//      std::cout << ", size : " << edges_.size() << std::endl;
	table_.clear();
	table_.reserve(n*2);
/*  VSK, 25.09.08. The lookup table is not in use an makes it more complex to add
    and remove faces from the topology table. */
/* 	lookup_.resize(n); */
/* 	for (i = 0; i<n ; ++i) */
/* 	    lookup_[i].reserve(4); */
	orientation_inconsist_.clear();  // No inconsistencies in orientation of
	                                 // faces yet

    }



    //=======================================================================
    int testEdges(edgeType* e[2])
    //=======================================================================
    {
	int k, l;
	// Find the endpoints
	Go::Point p[2][2];
	for (k = 0; k < 2; ++k) {
	    p[k][0] = e[k]->point(e[k]->tMin());
	    p[k][1] = e[k]->point(e[k]->tMax());
	}


	// Compare endpoints
	double dist;
	// The meaning of hithere[k] is that point p[0][k],
	// which is an endpoint of e[0], is:
	// -1: not close to e[1]
	// 0:  close to the tMin() e[1] endpoint
	// 1:  close to the tMax() e[1] endpoint
	// 2:  close to a point on edge e[1]
	// In cases 0-2, the corresponding param value on e[1] is in params[k]
	int hithere[2];
	int num_hits = 0;
	int pt_that_matched = -1;
	bool inner_hit = false;
	double params[2];
	double clo_t[2][2];
	Go::Point clo_pt[2][2];
	clo_pt[0][0] = clo_pt[0][1] = clo_pt[1][0] = clo_pt[1][1] = Go::Point(3);
	// For each endpoint of e[0]...
	for (k = 0; k < 2; ++k) {
	    hithere[k] = -1;
	    // ... we check against both endpoints of e[1]...
	    for (l = 0; l < 2; ++l) {
		dist = p[0][k].dist(p[1][l]);
		if (dist < tol_.neighbour) {
		    // We are close to an endpoint on e[1]!
		    hithere[k] = l;
		    pt_that_matched = k;
		    params[k] = (l == 0) ? e[1]->tMin() : e[1]->tMax();
		    break;
		}
	    }
	    // ... then, if that failed, we do a closest point
	    if (hithere[k] == -1) {
		try
		    {
			e[1]->closestPoint(p[0][k],
					   clo_t[0][k],
					   clo_pt[0][k],
					   dist);
		    }
		catch (...) {
		    dist = 2*tol_.neighbour;
		}
		if (dist < tol_.neighbour) {
		    // We are close to an point in the interior of e[1]!
		    hithere[k] = 2;
		    params[k] = clo_t[0][k];
		    inner_hit = true;
		}
	    }
	    if (hithere[k] != -1) ++num_hits;
	}

	// Check if we already have enough (2) hits.
	// If so, the whole of edge e[0] is marched
	if (num_hits == 2) {
	    double p0[2] = { e[0]->tMin(), e[0]->tMax() };
	    return march(e, p0, params);
	}

	// We have to check how the endpoints of e[1] line up against e[0].
	// If they match, it must be in the interior of e[0].
	bool onehit = (num_hits == 1);
	int from = 0;
	int to = 1;
	if (onehit && (pt_that_matched != -1)) { // For test to be valid we must have end_pt match.
	    if (hithere[pt_that_matched] == 0)
		from = 1;
	    else
		to = 0;
	}
	for (k = from; k <= to; ++k) {
	    try
		{
		    e[0]->closestPoint(p[1][k],
				       clo_t[1][k],
				       clo_pt[1][k],
				       dist);
		}
	    catch (...) {
		dist = 2*tol_.neighbour;
	    }
	    if (dist < tol_.neighbour) {
		++num_hits;
		// Now, we have two different cases.
		//    o We had one hit from the e[0] endpoints
		//    o We had no hits from the e[1] endpoints
		if (onehit) {
		    double p0[2];
		    double p1[2];
		    if (hithere[0] > -1) {
			p0[0] = e[0]->tMin();
			p0[1] = clo_t[1][k];
			p1[0] = params[0];
			p1[1] = (k==0) ? e[1]->tMin() : e[1]->tMax();
		    } else {      // That is, hithere[1] > -1
			p0[0] = clo_t[1][k];
			p0[1] = e[0]->tMax();
			p1[0] = (k==0) ? e[1]->tMin() : e[1]->tMax();
			p1[1] = params[1];
		    }
		    return march(e, p0, p1);
		}
	    }
	}

	// Finally, if we had 2 hits from endpoints of e[1] to interior points
	// of e[0]:
	if (num_hits == 2) {
	    double p1[2] = { e[1]->tMin(), e[1]->tMax() };
/* 	    if (((clo_t[1][0] > clo_t[1][1]) && !(e[0]->isTurned())) || */
/* 		((clo_t[1][0] < clo_t[1][1]) && (e[0]->isTurned()))) { */
	    if (clo_t[1][0] > clo_t[1][1]) {
		std::swap(clo_t[1][0], clo_t[1][1]);
		std::swap(p1[0], p1[1]);
	    }
	    return march(e, clo_t[1], p1);
	}

	// If we got here, we did not encounter any neighbourhood incidents
	return 0;
    }

// Assuming the edge given by (param0[0], param0[1]) has same direction as e[0].
    //=======================================================================
    int march(edgeType* e[2], double param0[2], double param1[2])
    //=======================================================================
    {
	// This function assumes that the point on edge 0 with parameter
	// param0[j] corresponds to the point on edge 1 with parameter
	// param1[j]. So param1[0] can be > param1[1]!

	if (std::min(fabs(param0[1]-param0[0]),
		     fabs(param1[1]-param1[0])) < 1e-10) {
	    // MESSAGE("Tiny incident ignored.");
	    return 0;
	}

	// @@sbr This should be replaced by a stable traversal routine.
	double num_tol = 1e-12;
	// We run through table_ to see if edges have already been connected.
	for (size_t ki = 0; ki < table_.size(); ++ki) {
	  if ((e[0] == table_[ki].e1_) &&
	      (fabs(e[0]->tMin() - std::min(param0[0], param0[1])) < num_tol)
	      &&
	      (fabs(e[0]->tMax() - std::max(param0[0], param0[1])) < num_tol)
	      &&
	      (e[1] == table_[ki].e2_) &&
	      (fabs(e[1]->tMin() - std::min(param1[0], param1[1])) < num_tol)
	      &&
	      (fabs(e[1]->tMax() - std::max(param1[0], param1[1])) < num_tol))
	    {
// 	      MESSAGE("Trying to connect already connected edges, "
// 		      "moving on to next edges.");
	      return 0;
	    }
	}

	tpTopologicalInfo topinfo;
	int status = march2(e, param0, param1, topinfo);
	if (status > 0)
	    addEntry(e[0], e[1],
		     std::min(param0[0], param0[1]),  std::max(param0[0], param0[1]),
		     std::min(param1[0], param1[1]),  std::max(param1[0], param1[1]),
		     topinfo, status);

	return status;
    }


    //=======================================================================
    int march2(edgeType* e[2], double param0[2], double param1[2],
	      tpTopologicalInfo& tinfo)
    //=======================================================================
	{
	Go::Point pnt[2][2];
	for (int k = 0; k < 2; ++k) {
	    pnt[0][k] = e[0]->point(param0[k]);
	    pnt[1][k] = e[1]->point(param1[k]);
	}
	double dist0 = pnt[0][0].dist(pnt[1][0]);
	double dist1 = pnt[0][1].dist(pnt[1][1]);

	double cos_ang0 = e[0]->normal(param0[0]).
	    cosAngle(e[1]->normal(param1[0]));
	double cos_ang1 = e[0]->normal(param0[1]).
	    cosAngle(e[1]->normal(param1[1]));

	// Let's march along the curve. We use an STL set to store our marched
	// points. The set keeps them sorted according to the
	// tpMarchPoint.par member (see tpMarchPoint::operator <).
	std::set<tpMarchPoint> mpset;
	mpset.insert(tpMarchPoint(pnt[0][0], param0[0], param1[0], dist0, cos_ang0,
				  status(dist0, cos_ang0)));
	mpset.insert(tpMarchPoint(pnt[0][1], param0[1], param1[1], dist1, cos_ang1,
				  status(dist1, cos_ang1)));

	// Make an initial subdivision to make closed-curve cases work
	double par = 0.5*param0[0] + 0.5*param0[1];
	Go::Point midpt0 = e[0]->point(par);
	double chordlen1 = pnt[0][0].dist(midpt0);
	double chordlen2 = pnt[0][1].dist(midpt0);
	if (std::max(chordlen1, chordlen2) < 0.5*tol_.neighbour) {
	    // MESSAGE("Trivial incident ignored.");
	    return 0;
	}
	// We must check whether closest point to other edge is within top gap.
	double dist, cos_ang, clo_par;
	std::vector<Go::Point> sf_seeds(2);
	getDistAndCosAngle(e, par, clo_par, dist, cos_ang, sf_seeds);
	// @@sbr An attempt to handle edges which correspond in end pts, but not in the middle.
	// Should already be treated by later tests, but it seems to fail if test_orientation == true.
	if (dist > tol_.neighbour) {
	    return 0;
	}
	mpset.insert(tpMarchPoint(e[0]->point(par), par, clo_par, dist, cos_ang,
				  status(dist, cos_ang)));

	// If the edges are oriented the same way, they will have normals
	// pointing in opposite directions. This condition will screw up
	// most functions trying to traverse the edge structure, such as
	// a function trying to find the outer edges. We test for that
	// condition and issue a warning.
    
/* 	MESSAGE_IF(param1[0] < param1[1], */
/* 		      "Marching edges with same orientation.\n" */
/* 		      << "This means the underlying surfaces have different " */
/* 		      << "orientation."); */

	// We subdivide the curve until the smallest distance from one point
	// to the next is less than max_step (but no more than 100 inserts are
	// done in any case).
	const double max_step = std::max(tol_.neighbour*1e2,
					 (chordlen1+chordlen2)/100.0);
    
	// Scan through set and compute distances to previous point.
	std::set<tpMarchPoint>::iterator prev = mpset.begin();
	std::set<tpMarchPoint>::iterator it = mpset.begin();
	++it;
	int inserted = 0;
	int max_inserts = 100; // @@sbr Hack!!! 10. Should be 100;
	sf_seeds[0].resize(0);
	sf_seeds[1].resize(0);
	for (; it != mpset.end(); ++it, ++prev) {
	    double d = it->pt.dist(prev->pt);
	    while (((prev->status != it->status) && d > tol_.neighbour)
		   || (d > max_step && ++inserted < max_inserts)) {
		// Distance between consecutive points too large
		par = 0.5*prev->par + 0.5*it->par;
		getDistAndCosAngle(e, par, clo_par, dist, cos_ang, sf_seeds);
		it = mpset
		    .insert(tpMarchPoint(e[0]->point(par), par, clo_par, dist,
					 cos_ang, status(dist, cos_ang)))
		    .first;
		d = it->pt.dist(prev->pt);
	    }
	}

	tinfo.parameters_.push_back(std::pair<double, double>(param0[0], param1[0]));
	// Add the entry (or entries)
	prev = mpset.begin();
	it = mpset.begin();
	++it;
	for (; it != mpset.end(); ++it) {
	    if (prev->status != it->status) {
		// In this AddEntry and the one below, the reason for using
		// the face-and-location instead of the specific edges
		// version of AddEntry is:
		// After adding an entry, the edges e[0] and e[1] have
		// been modified to be only a part of the original edge
		// and the parameter at which we want to split may no
		// longer be in that edge.
		tinfo.status_.push_back(prev->status);
		// MESSAGE_IF(prev->status >=4, "Edges being marched do"
		//	  << " not touch at some point(s).");
		tinfo.parameters_
		    .push_back(std::pair<double, double>(it->par, it->par2));
		prev = it;
	    }
	}

	// If we have not already come to the end, go to the end!
	if (tinfo.parameters_.back().first != param0[1]
	    && tinfo.parameters_.back().second != param1[1]) {
	    std::set<tpMarchPoint>::reverse_iterator rit = mpset.rbegin();
	    if (rit->status == prev->status) {
		tinfo.status_.push_back(prev->status);
		tinfo.parameters_
		    .push_back(std::pair<double, double>(param0[1], param1[1]));	    
	    }
	}

	return (param1[0] < param1[1]) ? 2 : 1;
    }


    //=======================================================================
    void getDistAndCosAngle(edgeType* e[2], double param,
			    double& other_par,
			    double& dist, double& cos_ang,
			    std::vector<Go::Point>& sf_seeds)
    //=======================================================================
    {
	Go::Point pt = e[0]->point(param);
	Go::Point clo_pt;
	e[1]->closestPoint(pt, other_par, clo_pt, dist);
	if (true) {
	    // Not that reliable, depending on isometric parametrization.
	   double seed; // We use some effort to come up with a seed.
	   // Typically e[0] and e[1] match in end pts, I assume
	   Go::Point from0 = e[0]->point(e[0]->tMin());
	   Go::Point to0 = e[0]->point(e[0]->tMax());
	   Go::Point from1 = e[1]->point(e[1]->tMin());
	   Go::Point to1 = e[1]->point(e[1]->tMax());
	   double frac =
	      (e[1]->tMax() - e[1]->tMin())*(param - e[0]->tMin())/(e[0]->tMax() - e[0]->tMin());
	   if (from0.dist(from1) < from0.dist(to1))
	      seed = e[1]->tMin() + frac;
	   else
	      seed = e[1]->tMax() - frac;
	   double seed_dist, seed_other_par;
	   e[1]->closestPoint(pt, seed_other_par, clo_pt, seed_dist, &seed);
	   if (seed_dist < dist) {
	      other_par = seed_other_par;
	      dist = seed_dist;
	   }
	}

	Go::Point normal1, normal2;
	// To speed things up we do not perform closest point calc if we don't have to.
	if (e[0]->geomEdge()->geomCurve()->instanceType() == Go::Class_CurveOnSurface) {
	    shared_ptr<Go::CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<Go::CurveOnSurface, Go::ParamCurve>
		(e[0]->geomEdge()->geomCurve());
	    // Face need not be created.
	    if (cv_on_sf->parPref() && (e[0]->face()->surface()).get() != 0) {
		Go::Point par_pt = cv_on_sf->parameterCurve()->point(param);
		try {
		    normal1 = e[0]->face()->normal(par_pt[0], par_pt[1]);
		} catch (...) {
		    MESSAGE("Failed evaluating normal!");
		}
	    }
	}
	if (normal1.size() == 0) {
	    double* sf_seed = (sf_seeds[0].size() != 2) ? NULL : sf_seeds[0].begin();
	    normal1 = e[0]->normal(param, sf_seeds[0], sf_seed);
	}
	if (e[1]->geomEdge()->geomCurve()->instanceType() == Go::Class_CurveOnSurface) {
	    shared_ptr<Go::CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<Go::CurveOnSurface, Go::ParamCurve>
		(e[1]->geomEdge()->geomCurve());
	    if (cv_on_sf->parPref() && (e[1]->face()->surface()).get() != 0) {
		Go::Point par_pt = cv_on_sf->parameterCurve()->point(other_par);
		try {
		    normal2 = e[1]->face()->normal(par_pt[0], par_pt[1]);
		} catch (...) {
		    MESSAGE("Failed evaluating normal!");
		}
	    }
	}
	if (normal2.size() == 0) {
	    double* sf_seed = (sf_seeds[1].size() != 2) ? NULL : sf_seeds[1].begin();
	    normal2 = e[1]->normal(other_par, sf_seeds[1], sf_seed);
	}

	cos_ang = normal1.cosAngle(normal2);
    }


    //=======================================================================
    int status(double dist, double cos_ang)
    //=======================================================================
    {
	if (dist > tol_.neighbour) return 4;
	if (dist > tol_.gap) return 3;
	if (cos_ang < cos(tol_.bend)) return 2;
	if (cos_ang < cos(tol_.kink)) return 1;
	return 0;
    }

};



#endif // _TPTOPOLOGYTABLE_H

