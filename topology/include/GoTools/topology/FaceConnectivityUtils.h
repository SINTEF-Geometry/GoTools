#ifndef _FACECONNECTIVITYUTILS_H
#define _FACECONNECTIVITYUTILS_H


#include "GoTools/topology/FaceConnectivity.h"

#include <vector>
#include <set>
#include <memory>
#include <fstream>


template <class edgeType, class faceType>
class FaceConnectivityUtils
{
 public:

  /// Returns all boundary loops. Keep in mind that
  /// calling next() on elements on the boundary will
  /// not iterate around the boundary, but around the face
  /// the edge comes from.
  //=======================================================================
  void BoundaryLoops(const std::vector<shared_ptr<faceType> >& faces,
		     std::vector< std::vector<edgeType*> > & loopvec)
  //=======================================================================
  {
    loopvec.clear();
    int loop = 0;
     std::set<edgeType* > used_edges;

     for (size_t i=0; i<faces.size(); ++i)
      {
	std::vector<shared_ptr<edgeType> > start_edges = faces[i]->startEdges();
	for (size_t j=0; j<start_edges.size(); ++j)
	  {
	    edgeType* edge = start_edges[j].get();
	    edgeType* orig = edge;
	    while (true) 
	      {
		bool in_set = (used_edges.find(edge) != used_edges.end());
		if (edge->twin() == 0 && !in_set)
		  {
		    loopvec.push_back(std::vector<edgeType*>());
		    while (!in_set) 
		      {
			used_edges.insert(edge);
			loopvec[loop].push_back(edge);
			edge = edge->next();
			std::set<edgeType* > tmp_edges; // To avoid infinite loop in case
			// of inconsistency
			tmp_edges.insert(edge);
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

		edge = edge->next();
		if (edge == orig)
		  break;
	      }
	  }
      }
  }


  //=======================================================================
  void BoundaryLoops(const std::vector<shared_ptr<faceType> >& faces,
		     std::vector< std::vector<shared_ptr<edgeType> > > & loopvec)
  //=======================================================================
  {
    loopvec.clear();
    std::set<edgeType* > used_edges;
    edgeType* e2;

    // Collect all edges
    // NB! This operation requires that createInitialEdges returns the existing
    // face edges if such edges exist, otherwise the boundary loops of all faces
    // will be collected
    std::vector<shared_ptr<edgeType> > edges;
    for (size_t i=0; i<faces.size(); ++i)
      {
	std::vector<shared_ptr<edgeType> > face_edges = faces[i]->createInitialEdges();
	edges.insert(edges.end(), face_edges.begin(), face_edges.end());
      }

    int loop = 0;
    for (size_t i = 0; i < edges.size(); ++i) {
      shared_ptr<edgeType> edge = edges[i];
      bool in_set = (used_edges.find(edge.get()) != used_edges.end());
      if (edges[i]->twin() == 0 && !in_set) {
	loopvec.push_back(std::vector<shared_ptr<edgeType> >());
	while (!in_set) {
	  used_edges.insert(edge.get());
	  loopvec[loop].push_back(edge);
	  e2 = edge->next();
	  for (size_t kr=0; kr<edges.size(); ++kr)
	    if (edges[kr].get() == e2)
	      {
		edge = edges[kr];
		break;
	      }
	  std::set<edgeType* > tmp_edges; // To avoid infinite loop in case
	  // of inconsistency
	  tmp_edges.insert(edge.get());
	  while (edge->twin())
	    {
	      e2 = edge->twin()->next();
	      for (size_t kr=0; kr<edges.size(); ++kr)
		if (edges[kr].get() == e2)
		  {
		    edge = edges[kr];
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
  void disjointObjects(const std::vector<shared_ptr<faceType> >& faces,
		       std::vector<std::vector<faceType*> >& grouped_faces)
  //=======================================================================
  {
    grouped_faces.clear();

    // Make a working copy of pointers to input faces
    std::vector<shared_ptr<faceType> > faces2;
    faces2.insert(faces2.end(), faces.begin(), faces.end());

    while (faces2.size() > 0)
      {
	// Store start face in connected set
	std::vector<faceType*> curr_group;
	shared_ptr<faceType> curr_face = faces2[0];      ;
	curr_group.push_back(curr_face.get());
	faces2.erase(faces2.begin());
	currDisjointObject(curr_face.get(), faces2, curr_group);
	grouped_faces.push_back(curr_group);
      }
  }

  /// Only gives the first of every pair of edges
  /// representing a cornering edge or kink edge.
  //=======================================================================
  void cornersAndKinks(const std::vector<shared_ptr<faceType> >& faces,
		       std::vector<edgeType*>& vec)
  //=======================================================================
  {
    vec.clear();
    for (size_t i=0; i<faces.size(); ++i)
      {
	std::vector<shared_ptr<edgeType> > start_edges = faces[i]->startEdges();
	for (size_t j=0; j<start_edges.size(); ++j)
	  {
	    edgeType* e = start_edges[j].get();
	    edgeType* orig = e;
	    while (true) 
	      {
		if (e->hasConnectivityInfo())
		  {
		    int status = e->getConnectivityInfo()->WorstStatus();
		    if (status > 0)
		      {
			// A continuity issue
			// Check if this instance is stored already
			size_t r;
			for (r=0; r<vec.size(); ++r)
			  {
			    if (vec[r] == e || vec[r] == e->twin())
			      break;
			  }
			if (r == vec.size())
			  vec.push_back(e);   // Only one edge in a twin pair is stored
		      }
		  }
		
		e = e->next();
		if (e == orig)
		  break;
	      }
	  }
      }
  }


  /// Returns all edges bounding a specific face
  //=======================================================================
  std::vector<edgeType*> edgesBoundingFace(faceType* face) 
    //=======================================================================
    {
      std::vector<edgeType*> bounding_edges;
      std::vector<shared_ptr<edgeType> > start_edges = face->startEdges();
      edgeType* e = 0;
      for (size_t i = 0; i < start_edges.size(); ++i) 
	{
	  e = start_edges[i].get();

	  edgeType* orig = e;
	  while (true) 
	    {
	      bounding_edges.push_back(e);
	      e = e->next();
	      if (e == orig)
		break;
	  }
	}
      return bounding_edges;
    }

 private:

  void currDisjointObject(faceType* curr,
			  std::vector<shared_ptr<faceType> >& faces,
			  std::vector<faceType*>& group)
  {
    // Store all neighbours
    std::vector<shared_ptr<edgeType> > start_edges = curr->startEdges();
    for (size_t kj=0; kj<start_edges.size(); ++kj)
      {
	edgeType* e1 = start_edges[kj].get();
	edgeType* orig = e1;
	while (true) 
	  {
	    edgeType* e2 = e1->twin();
	    if (e2)
	      {
		faceType* f2 = e2->face();
		
		// Check if this face is not grouped already
		size_t kr;
		for (kr=0; kr<faces.size(); ++kr)
		  if (faces[kr].get() == f2)
		    break;
		if (kr < faces.size())
		  {
		    group.push_back(f2);
		    faces.erase(faces.begin()+kr);
		    currDisjointObject(f2, faces, group);
		  }
	      }
	    e1 = e1->next();
	    if (e1 == orig)
	      break;
	  }
      }
  }




};

#endif // _FACECONNECTIVITYUTILS_H

