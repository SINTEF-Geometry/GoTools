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

#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/intersections/IntersectionInterface.h"
#include "GoTools/compositemodel/PointOnEdge.h"

#include <fstream>

using std::vector;
using std::make_pair;
using std::pair;
using std::min;

namespace Go
{

//===========================================================================
  Loop::Loop(ftFaceBase* face, CurveLoop& curve_loop, double kink, 
	     bool split_in_kinks, bool no_split)
// Constructor taking a loop of curves as input
// Edges are constructed as ftEdge
	: face_(face)
//===========================================================================
    {
      setEdges(curve_loop, kink, split_in_kinks, no_split);
    }

    
//===========================================================================
    // This constructor takes an ordered sequence of edges as input
    // Note that the function may throw
    Loop::Loop(ftFaceBase* face, vector<shared_ptr<ftEdgeBase> >& edges, 
	       double space_epsilon)
	: face_(face), eps_(space_epsilon)
//===========================================================================
    {
	setEdges(edges);
    }

//===========================================================================
    // This constructor takes an ordered sequence of edges as input
    // Note that the function may throw
    Loop::Loop(vector<shared_ptr<ftEdgeBase> >& edges, double space_epsilon)
	: face_(0), eps_(space_epsilon)
//===========================================================================
    {
	setEdges(edges);
    }

//===========================================================================
   Loop::~Loop()
// Destructor
//===========================================================================
   {
   }

    
//===========================================================================
  vector<shared_ptr<Vertex> > Loop::getVertices() const
//===========================================================================
    {
      std::set<shared_ptr<Vertex> > vertices;
      for (size_t ki=0; ki<edges_.size(); ++ki)
	{
	  ftEdge *curr = edges_[ki]->geomEdge();
	  if (curr)
	    {
	      shared_ptr<Vertex> v1, v2;
	      curr->getVertices(v1, v2);
	      vertices.insert(v1);
	      vertices.insert(v2);
	    }
	}
      vector<shared_ptr<Vertex> > vertices2;
      vertices2.insert(vertices2.end(), vertices.begin(), vertices.end());
      return vertices2;
    }

//===========================================================================
  vector<shared_ptr<Vertex> > Loop::getSeqVertices() const
//===========================================================================
  {
    vector<shared_ptr<Vertex> > result;
    for (size_t j = 0; j < edges_.size(); ++j)
      {
	shared_ptr<Vertex> vert = edges_[j]->geomEdge()->getVertex(false);
	result.push_back(vert);
      }

	   return result;
	 }

//===========================================================================
    /// Check consistency with regard to face
    bool Loop::isFaceConsistent()
//===========================================================================
    {
	for (size_t ki=0; ki<edges_.size(); ki++)
	    if (edges_[ki]->face() != face_)
		return false;

	return true;
    }

//===========================================================================
    void Loop::setFace(ftFaceBase* face)
//===========================================================================
    {
	face_ = face;

	// Update the face pointers of the edges
	for (size_t ki=0; ki<edges_.size(); ki++)
	{
	    ftFaceBase* edge_face = edges_[ki]->face();
	    ALWAYS_ERROR_IF(edge_face && edge_face != face, 
			    "Inconsistence in face pointers");

	    edges_[ki]->setFace(face);
	}

    }
//===========================================================================
    void Loop::setEdges(vector<shared_ptr<ftEdgeBase> >& edges)
//===========================================================================
    {
	ALWAYS_ERROR_IF(edges.size() == 0, "No edges in loop");

	// Check if the face pointer in the edges and the loop is consistent
	size_t ki;
	for (ki=0; ki<edges.size(); ki++)
	    ALWAYS_ERROR_IF(face_ && edges[ki]->face() != face_, "Face mismatch");

#ifndef NDEBUG
	{
	    std::ofstream fileout_debug("tmp/edges.g2");
	    for (ki = 0; ki < edges.size(); ++ki)
	    {
		shared_ptr<ParamCurve> geom_cv = edges[ki]->geomEdge()->geomCurve();
		double tmin = edges[ki]->tMin();
		double tmax = edges[ki]->tMax();
		try
		{
		    shared_ptr<ParamCurve> sub_cv(geom_cv->subCurve(tmin, tmax));
		    if (sub_cv.get() != NULL)
		    {
			sub_cv->writeStandardHeader(fileout_debug);
			sub_cv->write(fileout_debug);
		    }
		}
		catch (...)
		{
		    MESSAGE("Fail!");
		}
	    }
	}
#endif NDEBUG

	// A check on the consistence of the loop with respect to sequence and
	// orientation of edges should be implemented here. Use given tolerance.

	edges_.clear();
	edges_.reserve(edges.size());
	edges_.push_back(edges[0]);
	for (ki=1; ki<edges.size(); ki++)
	{
	    edges_.push_back(edges[ki]);
	    if (!edges_[ki-1]->next())
	      edges_[ki]->connectAfter(edges_[ki-1].get());
	}
	if (!edges_[edges_.size()-1]->next())
	    edges_[0]->closeLoop(edges_[edges_.size()-1].get());
    }

//===========================================================================
  void Loop::setEdges(CurveLoop& curve_loop, double kink, bool split_in_kinks,
		      bool no_split)
//===========================================================================
    {
	shared_ptr<ParamCurve> cv;
	ftEdgeBase* e;
	double degenerate_epsilon = curve_loop.getSpaceEpsilon();
	eps_ = degenerate_epsilon;
	edges_.clear();
	edges_.reserve(curve_loop.size());

	for (int kj = 0; kj < curve_loop.size(); ++kj) {
	    cv = curve_loop[kj];

	    double len = cv->estimatedCurveLength();
	    // if (len < 1.0e-2)
	    //   {
	    // 	std::cout << "Small curve length in Loop::setEdges. Length= " << len << std::endl;
	    //   }

	    vector<double> split_params;

	    // As class is used by tpTopologyTable, which only splits in endpoints, we
	    // treat the special cases of curves with inner G1-disc and loops repr as one curve.
	    if (cv->instanceType() == Class_SplineCurve ||
		cv->instanceType() == Class_CurveOnSurface) {
		shared_ptr<SplineCurve> spline_cv;
		if (cv->instanceType() == Class_SplineCurve)
		    spline_cv = dynamic_pointer_cast<SplineCurve, ParamCurve>(cv);
		else {
		    shared_ptr<CurveOnSurface> cv_on_sf =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		    if (cv_on_sf->parPref())
			spline_cv = dynamic_pointer_cast<SplineCurve, ParamCurve>
			    (cv_on_sf->parameterCurve());
		    else
			spline_cv = dynamic_pointer_cast<SplineCurve, ParamCurve>
			    (cv_on_sf->spaceCurve());
		}
		// We next must get G1 joints for spline curve.
		vector<double> cont(2);
		cont[0] = degenerate_epsilon; // @@ No reason to choose this I suppose...
		cont[1] = kink;
		if (split_in_kinks && spline_cv.get())
		  GeometryTools::getGnJoints(*spline_cv, cont, split_params);
		else
		  {
		    split_params.push_back(cv->startparam());
		    split_params.push_back(cv->endparam());
		  }
	    } else {
		MESSAGE("Unknown curve type. Curve may have inner tangent discontinuities.");
		split_params.push_back(cv->startparam());
		split_params.push_back(cv->endparam());
	    }

	    // We introduce another reason to split: if curve is a loop, and not degenerate.
	    // Otherwise two meeting loops (consisting of single curves) will not be connected.
	    // @@ This is a hack and not too elegant or natural!!!
	    // @@@ VSK, 0111. TESTING
	    //no_split = true;
	    if (split_params.size() == 2 && !no_split) {
		Point start_pt = cv->point(cv->startparam());
		Point end_pt = cv->point(cv->endparam());
		if (/*curve_loop.size() == 1 &&*/
		    (start_pt.dist(end_pt) < degenerate_epsilon) &&
		    (cv->estimatedCurveLength() > degenerate_epsilon)) {
		    int nmb_to_insert = 2; // Only one may yield troublesome outcome.
		    double tmin = split_params[0];
		    double tmax = split_params[1];
		    double tstep = (tmax - tmin)/(nmb_to_insert + 1);
		    for (int k = 1; k < 3; ++k)
			split_params.insert(split_params.begin() + split_params.size() - 1,
					    tmin + k*tstep);
		}
	    }

	    for (size_t kk = 0; kk < split_params.size() - 1; ++kk) {
		ASSERT(split_params[kk] < split_params[kk+1]);
		e = new ftEdge(face_, cv, split_params[kk], split_params[kk+1]);

		edges_.push_back(shared_ptr<ftEdgeBase> (e));
		if ((kk > 0) || (kj > 0 && kk == 0)) {
		    int ne = (int)edges_.size();
		    edges_[ne - 1]->connectAfter(edges_[ne - 2].get());
		}
	    }
	}
	int ne = (int)edges_.size();
	edges_[0]->closeLoop(edges_[ne - 1].get());
    }

//===========================================================================
void Loop::updateLoop(shared_ptr<ftEdgeBase> new_edge)
//===========================================================================
{
    // After a edge split, the loop may no longer be consistent. Remake the loop
    // from the first edge in the loop
    ftEdgeBase* first = edges_[0].get();
    size_t ki;
    for (ki=0; ki<edges_.size(); ++ki, first=first->next())
	if (edges_[ki].get() != first)
	{
	    // The position of the missing edge is found
	    edges_.insert(edges_.begin()+ki, new_edge);
	    break;
	}

    if (ki == edges_.size() && first == new_edge.get())
	edges_.insert(edges_.begin()+ki, new_edge);
}

//===========================================================================
void Loop::split(int ind, double par)
//===========================================================================
{
  // The edge split performs the task including splitting twin edges
  // and updating appropriate loops
  edges_[ind]->geomEdge()->split2(par);
}

//===========================================================================
bool Loop::isClose(ftEdge* edge,
		   RectDomain* domain,
		   double tol) const
//===========================================================================
{
  shared_ptr<ParamCurve> crv = edge->geomCurve();
  shared_ptr<SplineCurve> scurve = 
    shared_ptr<SplineCurve>(crv->geometryCurve());
  shared_ptr<CurveOnSurface> sf_cv = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(crv);
  shared_ptr<ParamSurface> srf = face_ -> asFtSurface() -> surface();

  
  Point param_last;
  bool p_last_found = false;

  vector<double> all_knots;
  scurve->basis().knotsSimple(all_knots);
  vector<double> knots;
  double t_min = edge->tMin();
  double t_max = edge->tMax();
  knots.push_back(t_min);
  for (size_t i = 0; i < all_knots.size(); ++i)
    if (t_min < all_knots[i] && all_knots[i] < t_max)
      knots.push_back(all_knots[i]);
  knots.push_back(t_max);

  int deg = scurve->order() - 1;

  for (size_t i = 0; i < knots.size()-1; ++i)
    {
      double step = (knots[i+1] - knots[i]) / double(deg);
      double par;
      int j;
      for (j=0, par=knots[i]; j <= deg; ++j, par+=step)
	{
	  if (j==deg && i < knots.size()-2) continue;
	  Point curve_p;
	  scurve -> point(curve_p, par);
	  double clo_u, clo_v, clo_dist;
	  // double epsilon = 1.0e-9;
	  Point p_now;

	  if (sf_cv.get())
	    {
	      Point sf_par = sf_cv->faceParameter(par);
	     srf->closestPoint(curve_p, clo_u, clo_v, p_now, clo_dist, 
			       tol /*epsilon*/, domain, 
			       sf_par.begin());
	    }
	  else if (p_last_found)
	     srf->closestPoint(curve_p, clo_u, clo_v, p_now, clo_dist, 
			       tol /*epsilon*/, domain, 
			       param_last.begin());
	  else
	    {
	      p_last_found = true;
	      srf->closestPoint(curve_p, clo_u, clo_v, p_now, clo_dist, 
				tol /*epsilon*/, domain);
	    }

	  param_last = Point(clo_u, clo_v);
	  if (clo_dist > tol) 
	    return false;
	}
    }

  return true;
}



//===========================================================================
void Loop::getBadDistance(vector<pair<ftSurface*, ftEdge*> >& badPairs,
			 RectDomain* domain,
			 double tol) const
//===========================================================================
{
  ftSurface* ft_surf = face_ -> asFtSurface();

  for (size_t i = 0; i < edges_.size(); ++i)
    {
      ftEdge* edge = edges_[i]->geomEdge();
      if (!isClose(edge, domain, tol))
	badPairs.push_back(make_pair(ft_surf, edge));
    }
}


//===========================================================================
void Loop::getBadDistance(vector<pair<ftEdge*, shared_ptr<Vertex> > >& badPairs,
			  double tol) const
//===========================================================================
{
  for (size_t i = 0; i < edges_.size(); ++i)
    {
      ftEdge* edge = edges_[i]->geomEdge();
      shared_ptr<Vertex> v0 = edge->getVertex(true);
      shared_ptr<Vertex> v1 = edge->getVertex(true);
      shared_ptr<ParamCurve> pcurve = edge->geomCurve();

      Point p0 = v0->getVertexPoint();
      Point p1 = v1->getVertexPoint();

      Point q0, q1;
      pcurve->point(q0, edge -> tMin());
      pcurve->point(q1, edge -> tMax());

      if (min(p0.dist2(q0), p0.dist2(q1)) > tol*tol)
	badPairs.push_back(make_pair(edge, v0));
      if (min(p1.dist2(q0), p1.dist2(q1)) > tol*tol)
	badPairs.push_back(make_pair(edge, v1));
    }
}


//===========================================================================
void Loop::getPosTangentSurfaceDiscont(vector<ftEdge*>& badPos,
				       vector<ftEdge*>& badTangent,
				       double tol, double kink, double bend, int leastSurfIndex,
				       shared_ptr<SurfaceModel> sm) const
//===========================================================================
{
  double cos2_kink = (cos(2.0*kink)+1)/2;
  double cos2_bend = (cos(2.0*bend)+1)/2;

  for (size_t i = 0; i < edges_.size(); ++i)
    {
      ftEdge* edge1 = edges_[i]->geomEdge();
      ftEdge* edge2;
      edge2 = (ftEdge*) edge1->twin();
      ftSurface* face1 = (ftSurface*) edge1->face();
      ftSurface* face2 = (ftSurface*) edge2->face();

      int twinIndex = sm->getIndex(face2);
      if (twinIndex < leastSurfIndex) continue;
      if (twinIndex == leastSurfIndex)
	{
	  size_t j = 0;
	  for (; j < i; ++j)
	    if (edges_[j].get() == edge2) break;
	  if (j < i) continue;
	}

      SplineCurve* scurve1 = edge1->geomCurve()->geometryCurve();

      Point p_last;
      bool p_last_found = false;

      vector<double> all_knots;
      scurve1->basis().knotsSimple(all_knots);
      vector<double> knots;
      double t_min = edge1->tMin();
      double t_max = edge1->tMax();
      knots.push_back(t_min);
      for (size_t i = 0; i < all_knots.size(); ++i)
	if (t_min < all_knots[i] && all_knots[i] < t_max)
	  knots.push_back(all_knots[i]);
      knots.push_back(t_max);

      int deg1 = scurve1->order() - 1;

      bool farPos = false;
      bool g1Discont = false;

      for (size_t j = 0; j < knots.size()-1; ++j)
	{
	  double step = (knots[j+1] - knots[j]) / double(deg1);
	  for (int k = 0; k <= deg1; ++k)
	    {
	      if (k==deg1 && j < knots.size()-2) continue;
	      Point p_curve_1;
	      scurve1 -> point(p_curve_1, knots[j] + step*double(k));
	      double clo_t, clo_dist;
	      Point p_curve_2;

	      if (p_last_found)
		edge2 -> closestPoint(p_curve_1, clo_t, p_curve_2, clo_dist, p_last.begin());
	      else
		{
		  p_last_found = true;
		  edge2 -> closestPoint(p_curve_1, clo_t, p_curve_2, clo_dist);
		}

	      p_last = p_curve_2;

	      if (!farPos && p_curve_2.dist2(p_curve_1) > tol*tol)
		farPos = true;

	      if (!g1Discont)
		{

		  Point p_surf_1, p_surf_2;
		  double clo_u_1, clo_v_1, clo_u_2, clo_v_2, clo_dist_1, clo_dist_2;
		  face1 -> closestPoint(p_curve_1, clo_u_1, clo_v_1, p_surf_1, clo_dist_1, 1.0e-9);
		  face2 -> closestPoint(p_curve_2, clo_u_2, clo_v_2, p_surf_2, clo_dist_2, 1.0e-9);

		  Point normal1 = face1 -> normal(clo_u_1, clo_v_1);
		  Point normal2 = face2 -> normal(clo_u_2, clo_v_2);

		  double norm2 = normal1.length2() * normal2.length2();
		  double normProd = normal1 * normal2;
		  if (norm2 * cos2_kink < normProd * normProd &&
		      normProd * normProd < norm2 * cos2_bend)
		    g1Discont = true;
		}
	      if (farPos && g1Discont) break;
	    }
	  if (farPos && g1Discont) break;
	}

      if (farPos) badPos.push_back(edge1);
      if (g1Discont) badTangent.push_back(edge1);

    }
}

//===========================================================================
// Check if the edges of the loop are consistent with the corresponding curves
// with regard to orientation
bool Loop::checkConsistency() const
//===========================================================================
{
    // First check orientation of edges compared to curves
    size_t ki;
    for (ki=0; ki<edges_.size(); ++ki)
    {
	bool isOK = edges_[ki]->orientationOK();
	if (!isOK)
	    return false;
    }
    return true;
}

//===========================================================================
// Check for acute edges in boundary loop
void Loop::getAcuteEdges(vector<pair<ftEdge*, ftEdge*> >& acute_edges, double angtol) const
//===========================================================================
{
    ftEdgeBase *first = edges_[0].get();
    ftEdgeBase *curr = first;
    ftEdgeBase *next = curr->next(); 

    // Traverse loop
    while (true)
    {
	// Evaluate tangents in common endpoint
	double t_currend = curr->isReversed() ? curr->tMin() : curr->tMax();
	double t_nextstart = next->isReversed() ? next->tMax() : next->tMin();
	Point tan1 = curr->tangent(t_currend);
	Point tan2 = next->tangent(t_nextstart);

	double ang = fabs(M_PI - tan1.angle(tan2));
	if (ang < angtol)
	    acute_edges.push_back(make_pair(curr->geomEdge(), next->geomEdge()));

	curr = next;
	next = curr->next();
	if (curr == first)
	    break;
    }

}

//===========================================================================
// Compute intersections between boundary loops
void Loop::getLoopIntersections(shared_ptr<Loop> loop2, double tol, 
				vector<pair<shared_ptr<PointOnEdge>, 
				shared_ptr<PointOnEdge> > >& int_pt) const
//===========================================================================
{
    if (this == loop2.get())
	return;   // Same loop

    // Compute intersections in all edge combinations
    size_t ki, kj, kr;
    for (ki=0; ki<edges_.size(); ++ki)
    {
	ftEdge  *e1 = edges_[ki]->geomEdge();
	if (!e1)
	    continue;
	shared_ptr<ParamCurve> crv1 = shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(),
										       e1->tMax()));
	for (kj=0; kj<loop2->size(); ++kj)
	{
	    ftEdge *e2 = loop2->getEdge(kj)->geomEdge();
	    if (!e2)
		continue;
	    shared_ptr<ParamCurve> crv2 = shared_ptr<ParamCurve>(e2->geomCurve()->subCurve(e2->tMin(),
											   e2->tMax()));
	    vector<pair<double, double> > intersections;
	    intersectCurves(crv1, crv2, tol, intersections);
	    for (kr=0; kr<intersections.size(); kr++)
	    {
		shared_ptr<PointOnEdge> pt_e1 = 
		    shared_ptr<PointOnEdge>(new PointOnEdge(e1, intersections[kr].first));
		shared_ptr<PointOnEdge> pt_e2 = 
		    shared_ptr<PointOnEdge>(new PointOnEdge(e2, intersections[kr].second));
		int_pt.push_back(make_pair(pt_e1, pt_e2));
	    }
	}
    }
	
}

//===========================================================================
// Compute intersections between boundary loops
void Loop::getLoopSelfIntersections(double tol, 
				    vector<pair<shared_ptr<PointOnEdge>, 
				    shared_ptr<PointOnEdge> > >& int_pt) const
//===========================================================================
{
    // Compute intersections in all edge combinations
    size_t ki, kj, kr;
    double epsilon = 1.0e-10;  // In removal of trivial intersections
    for (ki=0; ki<edges_.size(); ++ki)
    {
	ftEdge  *e1 = edges_[ki]->geomEdge();
	if (!e1)
	    continue;
	shared_ptr<ParamCurve> crv1 = shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(),
										       e1->tMax()));
	for (kj=ki+1; kj<edges_.size(); ++kj)
	{
	    ftEdge *e2 = edges_[kj]->geomEdge();
	    if (!e2)
		continue;
	    shared_ptr<ParamCurve> crv2 = shared_ptr<ParamCurve>(e2->geomCurve()->subCurve(e2->tMin(),
											   e2->tMax()));
	    std::vector<std::pair<double, double> > intersections;
	    intersectCurves(crv1, crv2, tol, intersections);

	    // Remove the trivial intersection between adjacent edges
	    if (e1->next() == e2)
	    {
		double t1 = e1->tMax();
		double t2 = e2->tMin();

		for (kr=0; kr<intersections.size(); ++kr)
		    if (fabs(intersections[kr].first-t1) < epsilon && 
			fabs(intersections[kr].second-t2) < epsilon)
			intersections.erase(intersections.begin()+kr);
	    }
	    if (e2->next() == e1)
	    {
		double t1 = e1->tMin();
		double t2 = e2->tMax();

		for (kr=0; kr<intersections.size(); ++kr)
		    if (fabs(intersections[kr].first-t1) < epsilon && 
			fabs(intersections[kr].second-t2) < epsilon)
			intersections.erase(intersections.begin()+kr);
	    }

	    // Store the remaining intersections
	    for (kr=0; kr<intersections.size(); kr++)
	    {
		shared_ptr<PointOnEdge> pt_e1 = 
		    shared_ptr<PointOnEdge>(new PointOnEdge(e1, intersections[kr].first));
		shared_ptr<PointOnEdge> pt_e2 = 
		    shared_ptr<PointOnEdge>(new PointOnEdge(e2, intersections[kr].second));
		int_pt.push_back(make_pair(pt_e1, pt_e2));
	    }
	}
    }
}


//===========================================================================
// Check for and if possible ensure loop correspondance
bool Loop::correspondingEdges(shared_ptr<Loop> other, double tol,
			      ftEdgeBase* &first1, ftEdgeBase* &first2,
			      bool& same_dir, bool no_snap)
//===========================================================================
{
  ftEdgeBase *e1, *e2;

  int nmbsample = 5;

  first1 = e1 = edges_[0].get();
  first2 = other->edges_[0].get();

  int ki, kj;
  bool same;
  double par1, par2, par3, par4;
  int idx;
  for (ki=0; ki<(int)other->edges_.size(); ++ki)
    {
       e2 = other->edges_[ki].get();

       bool overlaps = e1->checkOverlap(e2, tol, nmbsample, par1, par2, 
					par3, par4, same, no_snap);
       if (overlaps)
	 {
	   double par01, par02, par03, par04;
	   bool same0;
	   ftEdgeBase* e01 = (par2 < e1->tMax()) ? e1 : e1->next();
	   ftEdgeBase* e02;
	   if (same)
	     e02 = (par4 < e2->tMax()) ? e2 : e2->next();
	   else
	     e02 = (par3 > e2->tMin()) ? e2 : e2->prev();
	   bool overlaps2 = e01->geomEdge()->checkOverlap(e02, tol, 
							  nmbsample, 
							  par01,par02,  
							  par03, par04, 
							  same0, no_snap);
	   if (!overlaps2)
	     continue;  // Not correct edge in loop

	   idx = 0;
	   if (par1 > e1->tMin())
	     {
	       ftEdgeBase *e3 = e1->split(par1);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	       edges_.insert(edges_.begin()+1, tmp_edge);
	       e1 = e3;
	       idx = 1;
	     }
	   if (par2 < e1->tMax())
	     {
	       ftEdgeBase *e3 = e1->split(par2);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	       edges_.insert(edges_.begin()+idx+1, tmp_edge);
	     }
	   first1 = e1;

	   idx = 0;
	   if (par3 > e2->tMin())
	     {
	       ftEdgeBase *e3 = e2->split(par3);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	       if (same)
		 other->edges_.insert(other->edges_.begin()+ki+1, tmp_edge);
	       else
		   other->edges_.insert(other->edges_.begin()+ki, tmp_edge);
	       e2 = e3;
	       idx = 1;
	     }
	   if (par4 < e2->tMax())
	     {
	       ftEdgeBase *e3 = e2->split(par4);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	       if (same)
		 other->edges_.insert(other->edges_.begin()+ki+idx+1, tmp_edge);
	       else
		 other->edges_.insert(other->edges_.begin()+ki+idx, tmp_edge);
	     }
	   first2 = e2;

	   break;
	 }
    }

  if (ki == (int)other->edges_.size())
    {
      // No correspondance
      return false;
    }

  same_dir = same;
  e1 = first1->next();
  e2 = (same) ? first2->next() : first2->prev();
  kj = 0;
  while (kj < (int)edges_.size() && e1 != edges_[kj].get())
    kj++;
  ki = 0;
  while (ki < (int)other->edges_.size() && e2 != other->edges_[ki].get())
    ki++;
  while (e1 != first1)
    {
      if (kj == (int)edges_.size())
	kj = 0;
      if (same && ki == (int)other->edges_.size())
	ki = 0;
      else if (!same && ki<0)
	  ki = (int)other->edges_.size()-1;
      bool overlaps = e1->checkOverlap(e2, tol, nmbsample, par1, par2, 
				       par3, par4, same, no_snap);
      if (overlaps)
	{
	  if (par1 > e1->tMin())
	    {
	      return false;
	    }
	  if (par2 < e1->tMax())
	    {
	      ftEdgeBase *e3 = e1->split(par2);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	      edges_.insert(edges_.begin()+kj+1, tmp_edge);
	    }

	  idx = 0;
	  if (par3 > e2->tMin())
	    {
	      if (same)
		return false;

	      ftEdgeBase *e3 = e2->split(par3);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	      other->edges_.insert(other->edges_.begin()+ki+1, tmp_edge);
	      e2 = e3;
	      idx = 1;
	    }
	  if (par4 < e2->tMax())
	    {
	      if (!same)
		return false;
	      ftEdgeBase *e3 = e2->split(par4);
	      shared_ptr<ftEdgeBase> tmp_edge = shared_ptr<ftEdgeBase>(e3);
	      other->edges_.insert(other->edges_.begin()+ki+1+idx, tmp_edge);
	    }

	}
      ki = (same) ? ki+1 : ki-1;
      kj++;
      e1 = e1->next();
      e2 = (same) ? e2->next() : e2->prev();
    }
  
  return true;
}

//===========================================================================
bool Loop::hasRadialEdges() const
//===========================================================================
{
  for (size_t ki=0; ki<edges_.size(); ++ki)
    {
      ftEdge *curr = edges_[ki]->geomEdge();
      if (curr)
	{
	  if (curr->hasEdgeMultiplicity())
	    return true;
	}
    }
  return false;
}

//===========================================================================
bool Loop::allRadialEdges() const
//===========================================================================
{
  for (size_t ki=0; ki<edges_.size(); ++ki)
    {
      ftEdge *curr = edges_[ki]->geomEdge();
      if (!curr)
	return false;
      else
	{
	  if (!curr->hasEdgeMultiplicity())
	    return false;
	}
    }
  return true;
}

//===========================================================================
void Loop::closestPoint(const Point& pt, int& clo_ind, double& clo_par, 
		      Point& clo_pt, double& clo_dist) const
//===========================================================================
{
    clo_ind = 0;
    double tmp_par, tmp_dist;
    Point tmp_pt;
    edges_[0]->closestPoint(pt, clo_par, clo_pt, clo_dist);
    size_t ki;
    for (ki=1; ki < edges_.size(); ki++) {
	edges_[ki]->closestPoint(pt, tmp_par, tmp_pt, tmp_dist);
	if (tmp_dist < clo_dist) {
	    clo_dist = tmp_dist;
	    clo_pt = tmp_pt;
	    clo_par = tmp_par;
	    clo_ind = (int)ki;
	}
    }
}

//===========================================================================
void Loop::groupSmoothEdges(double tol, double angtol,
			    vector<vector<shared_ptr<ftEdgeBase> > >& edge_groups)
//===========================================================================
{
  int ki, kj;
  int nmbedge = (int)edges_.size();
  vector<shared_ptr<ftEdgeBase> > curr_group;
  for (ki=nmbedge-1, kj=0; kj<nmbedge; ki=(ki+1)%nmbedge, kj++)
    {
      Point p1 = edges_[ki]->point(edges_[ki]->tMax());
      Point p2 = edges_[kj]->point(edges_[kj]->tMin());
      Point tan1 = edges_[ki]->tangent(edges_[ki]->tMax());
      Point tan2 = edges_[kj]->tangent(edges_[kj]->tMin());
      
      double dist = p1.dist(p2);
      double ang = tan1.angle(tan2);
      if (curr_group.size() > 0 && (dist > tol || ang > angtol)) 
	{
	  edge_groups.push_back(curr_group);
	  curr_group.clear();
	}
      curr_group.push_back(edges_[kj]);
    }
  edge_groups.push_back(curr_group);

  // Check if the first and last edge group may be joined
  shared_ptr<ftEdgeBase> e1 = 
    edge_groups[edge_groups.size()-1][edge_groups[edge_groups.size()-1].size()-1];
  shared_ptr<ftEdgeBase> e2  = edge_groups[0][0];
  Point p1 = e1->point(e1->tMax());
  Point p2 = e2->point(e2->tMin());
  Point tan1 = e1->tangent(e1->tMax());
  Point tan2 = e2->tangent(e2->tMin());
      
  double dist = p1.dist(p2);
  double ang = tan1.angle(tan2);
  if (dist <= tol && ang <= angtol) 
    {
      edge_groups[edge_groups.size()-1].insert(edge_groups[edge_groups.size()-1].end(),
					       edge_groups[0].begin(), edge_groups[0].end());
      edge_groups.erase(edge_groups.begin());
    }	  
}

} // namespace Go
