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

#ifndef _FACEADJACENCY_H
#define _FACEADJACENCY_H

//#define DEBUG

#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/topology/FaceConnectivity.h"
#include "GoTools/topology/tpTolerances.h"

#include <vector>
#include <set>
#include <memory>
#include <fstream>

namespace Go
{

/// Helper class for marching
class MarchPoint
{
public:
    Go::Point pt;
    double par;
    double par2;
    double dist;
    double cos_ang;
    int status;
    MarchPoint(const Go::Point& p,  double pa, double pa2,
		 double d, double ca, int s)
	: pt(p), par(pa), par2(pa2),  dist(d), cos_ang(ca), status(s)
    {}
    bool operator < (const MarchPoint& other) const
    { return par < other.par; }
};
    
// edgeType and faceType encapsulate structures with certain given operations.
// As of now the structures edgeType and faceType fullfills these requirements.

/// Compute adjacency information for a set of surfaces
template <class edgeType, class faceType>
class FaceAdjacency
{
protected:

    tpTolerances tol_;
    std::vector< shared_ptr<edgeType> > new_edges_;  // Intermediate storage of new edges

public:

    /** Constructor.
     * \param tol_gap Tolerance for when two faces are assumed to be C0
     * continuous
     * \param tol_neighbour Tolerance for when two faces are assumed to
     * be neighbours
     * \param tol_kink Tolerance for when two adjacent faces are assumed to be
     * C1 continous
     * \param tol_bend Tolerance for when two adjacent faces are assumed
     * to have an intentially smooth connection
     */
 FaceAdjacency(double tol_gap, double tol_neighbour,
	       double tol_kink, double tol_bend)
   : tol_(tpTolerances(tol_gap, tol_neighbour, tol_kink, tol_bend))
      {}

    /** Constructor.
     * \param tol Topological tolerances
     */
  FaceAdjacency(const tpTolerances& tol)
   : tol_(tol)
    {}


    /** Destructor.
     * Detailed description.
     */
    ~FaceAdjacency()
      {
      }


    /// Changes the tolerances used in adjacency analysis
    //=======================================================================
    void setTolerances(const tpTolerances& tol)
    //=======================================================================
    {
      tol_ = tol;
    }


    /// Fetch the topological tolerances used in the adjaceny analysis
    tpTolerances getTolerances() 
    {
      return tol_;
    }


    /// Compute the adjecency between the given faces 
    //=======================================================================
    void 
      computeAdjacency(const std::vector<shared_ptr<faceType> >& faces,
		       int first_idx)
    //=======================================================================
    {
      std::vector<std::pair<faceType*,faceType*> > orient_inconsist;

      computeAdjacency(faces, orient_inconsist, first_idx);
    }

    /// Compute the adjecency between the given faces 
    /// \param faces The face set on which to compute adjacency
    /// \param orient_inconsist Information of adjacent faces where the 
    /// direction of the face normal is inconsistent
     //=======================================================================
    void 
      computeAdjacency(const std::vector<shared_ptr<faceType> >& faces,
		       std::vector<std::pair<faceType*,faceType*> >& orient_inconsist,
		       int first_idx)
    //=======================================================================
    {
      int i, j, k, l;
      int num_faces = (int)faces.size();
      std::vector<Go::BoundingBox> boxes;
      boxes.reserve(num_faces);
      // Make sure that the faces are equipped with edges and compute face boxes
      for (i = 0; i < num_faces; ++i) {
	(void)faces[i]->createInitialEdges(tol_.neighbour);
	boxes.push_back(faces[i]->boundingBox());
      }

      orient_inconsist.clear();

      std::vector<shared_ptr<edgeType> > startedges0, startedges1;
      for (i = 0; i < num_faces - 1; ++i) {
	int first = std::max(first_idx, i+1);
	for (j = first; j < num_faces; ++j) {
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
		  if (e[0]->boundingBox().overlaps(e[1]->boundingBox(),
						   tol_.neighbour)) {
/* #ifdef DEBUG */
/* 		  std::ofstream debug("top_debug.g2"); */
/* 		  for (int ki = 0; ki < 2; ++ki) { */
/* 		    e[ki]->face()->surface()->writeStandardHeader(debug); */
/* 		    e[ki]->face()->surface()->write(debug); */
/* 		    std::vector<double> pts(12); */
/* 		    Go::Point from = e[ki]->point(e[ki]->tMin()); */
/* 		    double tmid = 0.5*(e[ki]->tMin() + e[ki]->tMax()); */
/* 		    Go::Point mid = e[ki]->point(tmid); */
/* 		    Go::Point to = e[ki]->point(e[ki]->tMax()); */
/* 		    std::copy(from.begin(), from.end(), pts.begin()); */
/* 		    std::copy(mid.begin(), mid.end(), pts.begin() + 3); */
/* 		    std::copy(mid.begin(), mid.end(), pts.begin() + 6); */
/* 		    std::copy(to.begin(), to.end(), pts.begin() + 9); */
/* 		    Go::LineCloud lc(pts.begin(), 2); */
/* 		    lc.writeStandardHeader(debug); */
/* 		    lc.write(debug); */
/* 		  } */
/* #endif  */
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
			  for (kr=0; kr<orient_inconsist.size(); ++kr)
			    if ((orient_inconsist[kr].first == faces[i].get() &&
				 orient_inconsist[kr].second == faces[j].get()) ||
				(orient_inconsist[kr].first == faces[j].get() &&
				 orient_inconsist[kr].second == faces[i].get()))
			      break;

			  if (kr == orient_inconsist.size())
			    orient_inconsist.push_back(std::make_pair(faces[i].get(),
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
    }

    //=======================================================================
    /// Fetch existing adjacency information between faces and set
    /// topological pointers representing this adjacency.    
    void setConnectivity(const std::vector<shared_ptr<faceType> >& faces)
    //=======================================================================
    {
      for (size_t ki=0; ki<faces.size(); ++ki)
	{
	  std::vector<shared_ptr<edgeType> > startedges = faces[ki]->startEdges();
	  for (size_t kj=0; kj<startedges.size(); ++kj)
	    {
	      edgeType *e1 = startedges[kj].get();
	      edgeType *orig = e1;
	      while (true)
		{
		  if (e1->twin() && !e1->hasConnectivityInfo())
		    {
		      // Compute missing connectivity information
		      updateConnectivity(e1, e1->twin());
		    }

		  e1 = e1->next();
		  if (e1 == orig)
		    break;
		}
	    }
	}
    }

    
    //=======================================================================
    /// Fetch existing adjacency information between faces and set
    /// topological pointers representing this adjacency.    
    void setConnectivity(const std::vector<faceType*>& faces)
    //=======================================================================
    {
      for (size_t ki=0; ki<faces.size(); ++ki)
	{
	  std::vector<shared_ptr<edgeType> > startedges = faces[ki]->startEdges();
	  for (size_t kj=0; kj<startedges.size(); ++kj)
	    {
	      edgeType *e1 = startedges[kj].get();
	      edgeType *orig = e1;
	      while (true)
		{
		  if (e1->twin() && !e1->hasConnectivityInfo())
		    {
		      // Compute missing connectivity information
		      updateConnectivity(e1, e1->twin());
		    }

		  e1 = e1->next();
		  if (e1 == orig)
		    break;
		}
	    }
	}
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

	for (int i = 0; i < 2; ++i) {
	  twins[i] = ed[i];
	  if (maxs[i] < hp[i]) {
	    newedge = twins[i]->split(maxs[i]);
	    if (!newedge->geomEdge())
	      {
	    	std::cout << "Connecttwins1: " << newedge << " " << newedge->geomEdge() << std::endl;
	    	std::cout << i << ": " << twins[i] << ", " << twins[i]->geomEdge() << std::endl;
	      }

	    shared_ptr<edgeType> newtwin = shared_ptr<edgeType>(newedge);
	    newtwin->face()->updateBoundaryLoops(newtwin);
	    new_edges_.push_back(newtwin);
	  }

	  if (mins[i] > lp[i]) {
	    newedge = ed[i];
	    twins[i] =  newedge->split(mins[i]);
	    if (!twins[i]->geomEdge())
	      {
	    	std::cout << i << ": " << twins[i] << ", " << twins[i]->geomEdge() << std::endl;
	    	std::cout << "Connecttwins2: " << newedge << " " << newedge->geomEdge() << std::endl;
	      }
	    shared_ptr<edgeType> newtwin = shared_ptr<edgeType>(twins[i]);
	    newtwin->face()->updateBoundaryLoops(newtwin);
	    new_edges_.push_back(newtwin);
	  }
	}

	twins[0]->connectTwin(twins[1], status); // Connects both twins to each other

	return twins;
      }


    
    //=======================================================================
    /// Remove one face from the topological structures of the face set
    /// it belongs to
    void releaseFaceAdjacency(shared_ptr<faceType> face)
    //=======================================================================
    {
      // Fetch face edges
      std::vector<shared_ptr<edgeType> >  startedges =
	face->startEdges();

      // For each edge, unset twin pointers
      for (size_t ki=0; ki<startedges.size(); ki++)
	{
	  edgeType* e1 = startedges[ki].get();
	  edgeType* orig = e1;
	  while (true)
	    {
	      e1->disconnectTwin();
	      e1 = e1->next();
	      if (e1 == orig)
		break;
	    }
	}

      face->isolateFace();
      int stop_here;
      stop_here = 1;

      /*  // Fetch face edges */
      /* std::vector<shared_ptr<edgeType> > tmp_edges; */
      /* tmp_edges = face->createInitialEdges(tol_.neighbour, tol_.kink); */

      /* // For each edge, unset twin pointers */
      /* for (size_t ki=0; ki<tmp_edges.size(); ki++) */
      /* 	{ */
      /* 	  edgeType *curr = tmp_edges[ki].get(); */
      /* 	  // Remove twin information */
      /* 	  curr->disconnectTwin(); */
      /* 	} */

      /* face->isolateFace(); */
      /* int stop_here; */
      /* stop_here = 1; */
    }


    
    //=======================================================================
    /// Add one face to the topological structures of a face set 
    /// \param faces The set of faces where all topological information is
    /// computed
    /// \param new_face The face to add to the face set
    void computeFaceAdjacency(std::vector<shared_ptr<faceType> > faces,
			      shared_ptr<faceType> new_face)
    //=======================================================================
    {
      std::vector<std::pair<faceType*,faceType*> > orient_inconsist;
      computeFaceAdjacency(faces, new_face, orient_inconsist);
    }

    void computeFaceAdjacency(std::vector<faceType*> faces, faceType* new_face)
    //=======================================================================
    {
      std::vector<std::pair<faceType*,faceType*> > orient_inconsist;
      computeFaceAdjacency(faces, new_face, orient_inconsist);
    }

     //=======================================================================
    /// Add one face to the topological structures of a face set 
    /// \param faces The set of faces where all topological information is
    /// computed
    /// \param new_face The face to add to the face set
    /// \param orient_inconsist Information of adjacent faces where the 
    /// direction of the face normal is inconsistent  
    void computeFaceAdjacency(std::vector<shared_ptr<faceType> > faces,
			      shared_ptr<faceType> new_face,
			      std::vector<std::pair<faceType*,faceType*> >& orient_inconsist)
    //=======================================================================
    {
      orient_inconsist.clear();
      new_edges_.clear();  // Prepare intermediate storage

      // Fetch existing edges
      size_t ki, kj, kr;
      std::vector<shared_ptr<edgeType> > edges;
      for (ki=0; ki<faces.size(); ++ki)
	{
	  std::vector<shared_ptr<edgeType> > tmp_edges;
	  tmp_edges = faces[ki]->createInitialEdges(tol_.neighbour);
	  edges.insert(edges.end(), tmp_edges.begin(), tmp_edges.end());
	}

      // Fetch face edges
      std::vector<shared_ptr<edgeType> > tmp_edges;
      tmp_edges = new_face->createInitialEdges(tol_.neighbour);

      // Add new edges to list
      size_t nmb0 = edges.size();
      edges.insert(edges.end(), tmp_edges.begin(), tmp_edges.end());

      // Store edge boxes to avoid multiple computation
      size_t nmb1 = edges.size();
      std::vector<Go::BoundingBox> boxes;
      boxes.reserve((int)nmb1);
      for (ki=0; ki<nmb1; ki++)
	boxes.push_back(edges[ki]->boundingBox());

#ifdef DEBUG
      std::ofstream of("top.txt");
#endif

      edgeType* e[2];
      bool split1 = false, split2 = false;
      for (ki=nmb0; ki<edges.size(); )
	{
	  split1 = false;
	  for (kj=0; kj<edges.size(); )
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
		  e[0] = edges[kj].get();
		  e[1] = edges[ki].get();

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
		      for (kr=0; kr<orient_inconsist.size(); ++kr)
			if ((orient_inconsist[kr].first == e[0]->face() &&
			     orient_inconsist[kr].second == e[1]->face()) ||
			    (orient_inconsist[kr].first == e[1]->face() &&
			     orient_inconsist[kr].second == e[0]->face()))
			  break;

		      if (kr == orient_inconsist.size())
			orient_inconsist.push_back(std::make_pair(e[0]->face(),
								  e[1]->face()));
		    }

		  if (new_edges_.size() > 0)
		    {
		      // Some edge is split. Store new edges and make box
		      edges.insert(edges.end(), new_edges_.begin(), new_edges_.end());
		      new_edges_.clear();
		      for (kr=nmb1; kr<edges.size(); kr++)
			boxes.push_back(edges[kr]->boundingBox());
		      nmb1 = edges.size();

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
    /// Add one face to the topological structures of a face set 
    /// \param faces The set of faces where all topological information is
    /// computed
    /// \param new_face The face to add to the face set
    /// \param orient_inconsist Information of adjacent faces where the 
    /// direction of the face normal is inconsistent  
    void computeFaceAdjacency(std::vector<faceType*> faces,
			      faceType* new_face,
			      std::vector<std::pair<faceType*,faceType*> >& orient_inconsist)
    //=======================================================================
    {
      orient_inconsist.clear();
      new_edges_.clear();  // Prepare intermediate storage

      // Fetch existing edges
      size_t ki, kj, kr;
      std::vector<shared_ptr<edgeType> > edges;
      for (ki=0; ki<faces.size(); ++ki)
	{
	  std::vector<shared_ptr<edgeType> > tmp_edges;
	  tmp_edges = faces[ki]->createInitialEdges(tol_.neighbour);
	  edges.insert(edges.end(), tmp_edges.begin(), tmp_edges.end());
	}

      // Fetch face edges
      std::vector<shared_ptr<edgeType> > tmp_edges;
      tmp_edges = new_face->createInitialEdges(tol_.neighbour);

      // Add new edges to list
      size_t nmb0 = edges.size();
      edges.insert(edges.end(), tmp_edges.begin(), tmp_edges.end());

      // Store edge boxes to avoid multiple computation
      size_t nmb1 = edges.size();
      std::vector<Go::BoundingBox> boxes;
      boxes.reserve((int)nmb1);
      for (ki=0; ki<nmb1; ki++)
	boxes.push_back(edges[ki]->boundingBox());

#ifdef DEBUG
      std::ofstream of("top.txt");
#endif

      edgeType* e[2];
      bool split1 = false, split2 = false;
      for (ki=nmb0; ki<edges.size(); )
	{
	  split1 = false;
	  for (kj=0; kj<edges.size(); )
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
		  e[0] = edges[kj].get();
		  e[1] = edges[ki].get();

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
		      for (kr=0; kr<orient_inconsist.size(); ++kr)
			if ((orient_inconsist[kr].first == e[0]->face() &&
			     orient_inconsist[kr].second == e[1]->face()) ||
			    (orient_inconsist[kr].first == e[1]->face() &&
			     orient_inconsist[kr].second == e[0]->face()))
			  break;

		      if (kr == orient_inconsist.size())
			orient_inconsist.push_back(std::make_pair(e[0]->face(),
								  e[1]->face()));
		    }

		  if (new_edges_.size() > 0)
		    {
		      // Some edge is split. Store new edges and make box
		      edges.insert(edges.end(), new_edges_.begin(), new_edges_.end());
		      new_edges_.clear();
		      for (kr=nmb1; kr<edges.size(); kr++)
			boxes.push_back(edges[kr]->boundingBox());
		      nmb1 = edges.size();

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
    /// Update topological information related to two neighbouring faces,
    /// i.e. recompute information regarding the type of connectivity between
    /// the faces
    void updateConnectivity(edgeType* e1, edgeType *e2)
    //=======================================================================
    {
      // Compute new topological information
      shared_ptr<FaceConnectivity<edgeType> > topinfo;
      double param0[2], param1[2];
      param0[0] = e1->tMin();
      param0[1] = e1->tMax();
      Go::Point p1, p3, p4;
      p1 = e1->point(e1->tMin());
      p3 = e2->point(e2->tMin());
      p4 = e2->point(e2->tMax());
      if (p1.dist(p3) < p1.dist(p4))
	{
	  // Same orientation of edge
	  param1[0] = e2->tMin();
	  param1[1] = e2->tMax();
	}
      else
	{
	  // Opposite orientation
	  param1[0] = e2->tMax();
	  param1[1] = e2->tMin();
	}

      int status;
      edgeType* curr_edges[2];
      curr_edges[0] = e1;
      curr_edges[1] = e2;
      status = march2(curr_edges, param0, param1, topinfo);

      e1->setConnectivityInfo(topinfo);
      e2->setConnectivityInfo(topinfo);
    }

    
 private:



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
      double dd[2];
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
	for (l = 0; l < 2; ++l) 
	  dd[l] = p[0][k].dist(p[1][l]);

	if (dd[0] < tol_.neighbour && dd[0] < dd[1])
	  {
	    // We are close to an endpoint on e[1]!
	    hithere[k] = 0;
	    pt_that_matched = k;
	    params[k] =  e[1]->tMin();
	  }
	else if (dd[1] < tol_.neighbour)
	  {
	    // We are close to an endpoint on e[1]!
	    hithere[k] = 1;
	    pt_that_matched = k;
	    params[k] =  e[1]->tMax();
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

      /* 	// @@sbr This should be replaced by a stable traversal routine. */
      /* 	double num_tol = 1e-12; */
      /* 	// We run through table_ to see if edges have already been connected. */
      /* 	for (size_t ki = 0; ki < table_.size(); ++ki) { */
      /* 	  if ((e[0] == table_[ki].e1_) && */
      /* 	      (fabs(e[0]->tMin() - std::min(param0[0], param0[1])) < num_tol) */
      /* 	      && */
      /* 	      (fabs(e[0]->tMax() - std::max(param0[0], param0[1])) < num_tol) */
      /* 	      && */
      /* 	      (e[1] == table_[ki].e2_) && */
      /* 	      (fabs(e[1]->tMin() - std::min(param1[0], param1[1])) < num_tol) */
      /* 	      && */
      /* 	      (fabs(e[1]->tMax() - std::max(param1[0], param1[1])) < num_tol)) */
      /* 	    { */
      /* // 	      MESSAGE("Trying to connect already connected edges, " */
      /* // 		      "moving on to next edges."); */
      /* 	      return 0; */
      /* 	    } */
      /* 	} */

      shared_ptr<FaceConnectivity<edgeType> > topinfo;
      int status = march2(e, param0, param1, topinfo);
      if (status > 0)
	{
	  std::vector<edgeType*> twins(2);
	  twins = connectTwins(e[0], e[1], 
			       std::min(param0[0], param0[1]),  std::max(param0[0], param0[1]),
			       std::min(param1[0], param1[1]),  std::max(param1[0], param1[1]),
			       status);
	  topinfo->setEdges(twins[0], twins[1]);
	  twins[0]->setConnectivityInfo(topinfo);
	  twins[1]->setConnectivityInfo(topinfo);
	}

      return status;
    }


    //=======================================================================
    int march2(edgeType* e[2], double param0[2], double param1[2],
	       shared_ptr<FaceConnectivity<edgeType> >& tinfo)
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
      // MarchPoint.par member (see MarchPoint::operator <).
      std::set<MarchPoint> mpset;
      mpset.insert(MarchPoint(pnt[0][0], param0[0], param1[0], dist0, cos_ang0,
			      status(dist0, cos_ang0)));
      mpset.insert(MarchPoint(pnt[0][1], param0[1], param1[1], dist1, cos_ang1,
			      status(dist1, cos_ang1)));

      // Make an initial subdivision to make closed-curve cases work
      double par = 0.5*param0[0] + 0.5*param0[1];
      Go::Point midpt0 = e[0]->point(par);
      double chordlen1 = pnt[0][0].dist(midpt0);
      double chordlen2 = pnt[0][1].dist(midpt0);

      if (chordlen1+chordlen2 < std::max(dist0, dist1))
	return 0;  // The curve lenght is less than the distance between the curves
      // in one of their endpoints. Likely that we march in the wrong direction

      /* if (std::max(chordlen1, chordlen2) < 0.5*tol_.neighbour) { */
      /* 	// MESSAGE("Trivial incident ignored."); */
      /* 	return 0; */
      /* } */
      // We must check whether closest point to other edge is within top gap.
      double dist, cos_ang, clo_par;
      std::vector<Go::Point> sf_seeds(2);
      getDistAndCosAngle(e, par, clo_par, dist, cos_ang, sf_seeds);
      // @@sbr An attempt to handle edges which correspond in end pts, but not in the middle.
      // Should already be treated by later tests, but it seems to fail if test_orientation == true.
      if (dist > tol_.neighbour) {
#ifdef DEBUG
	std::ofstream of("adjacency.g2");
	e[0]->face()->surface()->writeStandardHeader(of); 
	e[0]->face()->surface()->write(of); 
	e[1]->face()->surface()->writeStandardHeader(of); 
	e[1]->face()->surface()->write(of); 
	of << "400 1 0 4 100 150 0 255" << std::endl;
	of << " 1" << std::endl;
	of << midpt0 << std::endl;
#endif 

	return 0;
      }
      mpset.insert(MarchPoint(e[0]->point(par), par, clo_par, dist, cos_ang,
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
      std::set<MarchPoint>::iterator prev = mpset.begin();
      std::set<MarchPoint>::iterator it = mpset.begin();
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
	    .insert(MarchPoint(e[0]->point(par), par, clo_par, dist,
			       cos_ang, status(dist, cos_ang)))
	    .first;
	  d = it->pt.dist(prev->pt);
	}
      }

      tinfo = shared_ptr<FaceConnectivity<edgeType> >(new FaceConnectivity<edgeType>(e[0], e[1]));
      tinfo->parameters_.push_back(std::pair<double, double>(param0[0], param1[0]));
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
	  tinfo->status_.push_back(prev->status);
	  // MESSAGE_IF(prev->status >=4, "Edges being marched do"
	  //	  << " not touch at some point(s).");
	  tinfo->parameters_
	    .push_back(std::pair<double, double>(it->par, it->par2));
	  prev = it;
	}
      }

      // If we have not already come to the end, go to the end!
      if (tinfo->parameters_.back().first != param0[1]
	  && tinfo->parameters_.back().second != param1[1]) {
	std::set<MarchPoint>::reverse_iterator rit = mpset.rbegin();
	if (rit->status == prev->status) {
	  tinfo->status_.push_back(prev->status);
	  tinfo->parameters_
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
	if (cv_on_sf->parPref() && e[0]->face()->surface().get() != 0) {
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
	if (cv_on_sf->parPref() && e[1]->face()->surface().get() != 0) {
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
      //if (cos_ang <  cos(tol_.bend)) return 2;
      if (fabs(cos_ang) <  cos(tol_.bend)) return 2;
      if (cos_ang < cos(tol_.kink)) return 1;
      return 0;
    }

};

} // namespace Go

#endif // _FACEADJACENCY_H

