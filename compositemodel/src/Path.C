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

#include "GoTools/compositemodel/Path.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include <fstream>

using std::vector;

namespace Go {

//==========================================================================
bool Path::estimateHoleInfo(const vector<ftEdge*>& edges, Point& centre, 
			    Point& axis, double& radius, double& angle)
//==========================================================================
{
  centre.resize(edges[0]->geomCurve()->dimension());
  centre.setValue(0.0);

  // Select four points
  // First make parameterization
  vector<double> parval(edges.size()+1);
  vector<double> edgelen(edges.size()+1);
  parval[0] = 0.0;
  edgelen[0] = 0.0;
  size_t ki;
  for (ki=0; ki<edges.size(); ++ki)
    {
      double tdel = edges[ki]->tMax() - edges[ki]->tMin();
      double len = edges[ki]->estimatedCurveLength();
      parval[ki+1] = parval[ki]+tdel;
      edgelen[ki+1] = edgelen[ki] + len;
      centre += edges[ki]->point(edges[ki]->tMin());
    }
  int nmb_vx = (int)edges.size();
  if (edges[edges.size()-1]->next() != edges[0])
    {
      centre += edges[edges.size()-1]->point(edges[edges.size()-1]->tMax());
      nmb_vx++;
    }
  centre /= nmb_vx;

  vector<Point> pnt(4);
  int kj;
  //double tdel = (parval[parval.size()-1] - parval[0])/(double)4;
  double del = (edgelen[edgelen.size()-1] - edgelen[0])/(double)4;
  double tpar, len;
  for (kj=0, len=edgelen[0]+0.5*del; kj<4; ++kj, len+=del)
    {
      for (ki=0; ki<edgelen.size()-1; ++ki)
	if (len < edgelen[ki+1])
	  break;

      double frac = (len - edgelen[ki])/(edgelen[ki+1] - edgelen[ki]);
      tpar = parval[ki] + frac*(parval[ki+1]-parval[ki]);

      pnt[kj] = edges[ki]->point(edges[ki]->tMin() + tpar - parval[ki]);
    }

  Point centre2 = 0.25*(pnt[0] + pnt[1] + pnt[2] + pnt[3]);

  // std::ofstream of("hole_pts.g2");
  // of << "400 1 0 4 0 155 100 255" << std::endl;
  // of << pnt.size() << std::endl;
  // for (ki=0; ki<pnt.size(); ++ki)
  //   of << pnt[ki] << std::endl;
  // of << "400 1 0 4 255 0 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << centre << std::endl;
  // of << "400 1 0 4 0 255 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << centre2 << std::endl;

  centre = centre2;
  //axis = (pnt[1] - pnt[0]).cross(pnt[3] - pnt[2]);
  axis = (pnt[2] - pnt[0]).cross(pnt[3] - pnt[1]);
  double axlen = axis.length();
  double lentol = 1.0e-10;
  if (axlen < lentol)
    return false;

  axis.normalize();
  
  Point x1 = 0.5*(pnt[0] + pnt[1]);
  Point x2 = 0.5*(pnt[1] + pnt[2]);

  Point d1 = (pnt[1]-pnt[0]).cross(axis);
  d1.normalize();
  Point d2 = (pnt[2]-pnt[1]).cross(axis);
  d2.normalize();

  Point tmp1 = x1 + d1;
  Point tmp2 = x2 + d2;

  //double d1d2 = d1*d2;
  //double tdiv = 1.0 - d1d2*d1d2;
  //double t1 = (-x1*d1 + x2*d1 + d1d2*x1*d2 - d1d2*x2*d2)/tdiv;

  //centre = x1 + t1*d1;
  radius = (centre - pnt[0]).length();

//   // Adjust radius relative to bounding box
//   vector<Point> pos;
//   pos.insert(pos.end(), pnt.begin(), pnt.end());
//   for (ki=0; ki<edges.size(); ++ki)
//     {
//       pos.push_back(edges[ki]->point(edges[ki]->tMin()));
//       pos.push_back(edges[ki]->point(edges[ki]->tMax()));
//     }

//   BoundingBox box;
//   box.setFromPoints(pos);
//   double len = box.high().dist(box.low());
//   radius = std::min(radius, 0.5*len);
  Point vec1 = edges[0]->tangent(edges[0]->tMin());
  Point vec2 = edges[edges.size()-1]->tangent(edges[edges.size()-1]->tMax());
  vec2 *= -1;
  angle = vec1.angle(vec2);
  if (vec1*vec2 < 0.0)
    angle = 2*M_PI - angle;

  return true;
}

//==========================================================================
  void Path::classifyCorners(const vector<ftEdge*>& edges, double tol,
			     vector<shared_ptr<Vertex> >& corner,
			     vector<shared_ptr<Vertex> >& non_corner)
			     
//==========================================================================
  {
    corner.clear();
    non_corner.clear();

    size_t ki;
    size_t nmb = edges.size();
    shared_ptr<Vertex> vx1 = edges[0]->getVertex(true);
    shared_ptr<Vertex> vx2 = edges[edges.size()-1]->getVertex(false);
    if (vx1.get() == vx2.get())
      {
	double t1 = edges[nmb-1]->parAtVertex(vx1.get());
	double t2 = edges[0]->parAtVertex(vx1.get());
	Point tan1 = edges[nmb-1]->tangent(t1);
	Point tan2 = edges[0]->tangent(t2);
	double ang = tan1.angle(tan2);
	if (ang < tol)
	  non_corner.push_back(vx1);
	else
	  corner.push_back(vx1);
      }
    else
      corner.push_back(vx1);
    for (ki=1; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> vx = edges[ki]->getVertex(true);
	double t1 = edges[ki-1]->parAtVertex(vx.get());
	double t2 = edges[ki]->parAtVertex(vx.get());
	Point tan1 = edges[ki-1]->tangent(t1);
	Point tan2 = edges[ki]->tangent(t2);
	double ang = tan1.angle(tan2);
	if (ang < tol)
	  non_corner.push_back(vx);
	else
	  corner.push_back(vx);
      }
    if (vx1.get() != vx2.get())
      corner.push_back(vx2);
  }

//==========================================================================
  vector<ftEdge*> Path::identifyLoop(vector<ftEdge*> edges, shared_ptr<Vertex> vx)
//==========================================================================
{
  vector<ftEdge*> loop;

  // Traverse the edge list and check if the given vertex can be found in
  // that list
  size_t ki, kj;
  shared_ptr<Vertex> curr_vx;
  for (ki=0; ki<edges.size()-1; ++ki)
    {
      shared_ptr<Vertex> tmp = edges[ki]->getVertex(true);
      if (edges[ki+1]->hasVertex(tmp.get()))
	curr_vx = edges[ki]->getVertex(false);
      else curr_vx = tmp;
	  
      if (curr_vx.get() == vx.get())
	break;
    }

  if (ki == edges.size())
    return loop;   // No loop is found, returning empty vector

  for (kj=ki+1; kj<edges.size(); ++kj)
    {
       shared_ptr<Vertex> tmp = edges[kj]->getVertex(true);
      if (edges[kj-1]->hasVertex(tmp.get()))
	curr_vx = edges[kj]->getVertex(false);
      else curr_vx = tmp;
	  
      if (curr_vx.get() == vx.get())
	break;
    }

  if (kj == edges.size())
    return loop;  // The vertex can be found only once in the list of edges

  // A loop is found, collect it
  loop.insert(loop.end(), edges.begin()+ki, edges.begin()+kj+1);
  return loop;
}

//==========================================================================
void Path::closestPoint(vector<ftEdge*> edges, const Point& pt, 
			int& clo_ind, double& clo_par, 
			Point& clo_pt, double& clo_dist)
//==========================================================================
{
    clo_ind = 0;
    double tmp_par, tmp_dist;
    Point tmp_pt;
    edges[0]->closestPoint(pt, clo_par, clo_pt, clo_dist);
    size_t ki;
    for (ki=1; ki < edges.size(); ki++) {
	edges[ki]->closestPoint(pt, tmp_par, tmp_pt, tmp_dist);
	if (tmp_dist < clo_dist) {
	    clo_dist = tmp_dist;
	    clo_pt = tmp_pt;
	    clo_par = tmp_par;
	    clo_ind = (int)ki;
	}
    }
}

//===========================================================================
void Path::getEdgeCurves(vector<ftEdge*>& loop, 
			 vector<shared_ptr<ParamCurve> >& space_cvs,
			 vector<Point>& joint_points,
			 double eps, double tol,
			 bool corner_in_Tjoint)
//===========================================================================
{
  // Modify the start point of the loop such that there is a corner
  // between the last and the first edge
  size_t ki, kr;
  for (ki=0; ki<loop.size(); ++ki)
    {
      kr = loop.size()-1;
      double ang = M_PI;
      shared_ptr<Vertex> common_vx = loop[kr]->getCommonVertex(loop[0]);
      if (common_vx->nmbUniqueEdges() == 2 || !corner_in_Tjoint)
	{
	  // Check for corner
	  double t1 = loop[kr]->parAtVertex(common_vx.get());
	  double t2 = loop[0]->parAtVertex(common_vx.get());
	  Point tan1 = loop[kr]->tangent(t1);
	  Point tan2 = loop[0]->tangent(t2);
	  ang = tan1.angle(tan2);
	}
      if (ang < tol)
	{
	  loop.push_back(loop[0]);
	  loop.erase(loop.begin());
	}
      else
	break;
    }

      
  // Sort edges into sequences between each corner between edges
  vector<vector<ftEdge*> > joined_loop;
  for (ki=0; ki<loop.size(); ++ki)
    {
      double ang = M_PI;
      if (ki > 0)
	{
	  // Check if there is a corner or T-joint between this curve and
	  // the previous
	  shared_ptr<Vertex> common_vx = loop[ki-1]->getCommonVertex(loop[ki]);
	  if (common_vx->nmbUniqueEdges() == 2 || !corner_in_Tjoint)
	    {
	      // Check for corner
	      double t1 = loop[ki-1]->parAtVertex(common_vx.get());
	      double t2 = loop[ki]->parAtVertex(common_vx.get());
	      Point tan1 = loop[ki-1]->tangent(t1);
	      Point tan2 = loop[ki]->tangent(t2);
	      ang = tan1.angle(tan2);
	    }
	}
      if (ang >= tol)
	{
	  vector<ftEdge*> curr_loop;
	  curr_loop.push_back(loop[ki]);
	  joined_loop.push_back(curr_loop);
	}
      else
	joined_loop[joined_loop.size()-1].push_back(loop[ki]);
    }
 
  if (joined_loop.size() < 4 && loop.size() >= 4)
    {
      // Ensure that there is enough curves to create a coons patch.
      // This split could be done in a more smart way, but don't 
      // expect this to be a probable case
      for (ki=0; ki<joined_loop.size(); ++ki)
	{
	  if (joined_loop[ki].size() > 1)
	    {
	      vector<ftEdge*> curr_loop(joined_loop[ki].begin()+1,
					joined_loop[ki].end());
	      joined_loop[ki].erase(joined_loop[ki].begin()+1, 
				    joined_loop[ki].end());
	      joined_loop.insert(joined_loop.begin()+ki, curr_loop);
	      if (joined_loop.size() == 4)
		break;
	    }
	}
    }
	      
  
  
  // Fetch curves, and make parameter curves corresponding to the volume
  // Join curves that should belong to the same boundary curve of the
  // missing surface, but are represented as several edges
  // Remember the positions at these joints
  space_cvs.resize(joined_loop.size());
  for (ki=0; ki<joined_loop.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = 
	shared_ptr<ParamCurve>(joined_loop[ki][0]->geomCurve()->subCurve(joined_loop[ki][0]->tMin(),
									 joined_loop[ki][0]->tMax()));
      for (kr=1; kr<joined_loop[ki].size(); ++kr)
	{
	  shared_ptr<ParamCurve> tmp_cv2 = 
	    shared_ptr<ParamCurve>(joined_loop[ki][kr]->geomCurve()->subCurve(joined_loop[ki][kr]->tMin(),
									      joined_loop[ki][kr]->tMax()));

	  // NOTE. In the curve-on-surface case, only the space curve is
	  // considered
	  shared_ptr<CurveOnSurface> sfcv1 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_cv);
	  if (sfcv1.get())
	    {
	      sfcv1->ensureSpaceCrvExistence(eps);
	      tmp_cv = sfcv1->spaceCurve();
	    }
	  shared_ptr<CurveOnSurface> sfcv2 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_cv2);
	  if (sfcv2.get())
	    {
	      sfcv2->ensureSpaceCrvExistence(eps);
	      tmp_cv2 = sfcv2->spaceCurve();
	    }

	  // Make sure that the curves are consistently oriented
	  Point pos1 = tmp_cv->point(tmp_cv->startparam());
	  Point pos2 = tmp_cv->point(tmp_cv->endparam());
	  Point pos3 = tmp_cv2->point(tmp_cv2->startparam());
	  Point pos4 = tmp_cv2->point(tmp_cv2->endparam());
	  double d1 = pos1.dist(pos3);
	  double d2 = pos1.dist(pos4);
	  double d3 = pos2.dist(pos3);
	  double d4 = pos2.dist(pos4);
	  if (std::min(d3, d4) > std::min(d1, d2))
	    {
	      tmp_cv->reverseParameterDirection();
	      std::swap(pos1, pos2);
	      std::swap(d1, d3);
	      std::swap(d2, d4);
	    }
	  if (d4 < d3)
	    tmp_cv2->reverseParameterDirection();
	  Point joint = tmp_cv2->point(tmp_cv2->startparam());
	  joint_points.push_back(joint);
	  double dist;
	  vector<Point> pts1(2);
	  tmp_cv->point(pts1, tmp_cv->endparam(), 1);
	  vector<Point> pts2(2);
	  tmp_cv2->point(pts2, tmp_cv2->startparam(), 1);
	  double fac = pts2[1].length()/pts1[1].length();

	  // TEST
	  fac = 1.0;

	  double s1 = tmp_cv->startparam();
	  double s2 = tmp_cv->endparam();
	  double t1 = tmp_cv2->startparam();
	  double t2 = tmp_cv2->endparam();
	  double len1 = tmp_cv->estimatedCurveLength();
	  double len2 = tmp_cv2->estimatedCurveLength();
	  fac = len2*(s2-s1)/(len1*(t2-t1));

	  tmp_cv2->setParameterInterval(t1, t1+fac*(t2-t1));
	  tmp_cv->appendCurve(tmp_cv2.get(), 0, dist, false);
	}
      space_cvs[ki] = tmp_cv;
    }
  }

//===========================================================================
vector<ftEdge*> Path::edgeChain(ftEdge *edg, double angtol, shared_ptr<Vertex>& v1,
				shared_ptr<Vertex>& v2)
//===========================================================================
{
  // Extract edge chain with no corners and no joints between more than two edges 
  // in the same underlying surface
  // Fetch surface
  shared_ptr<ParamSurface> surf = edg->face()->surface();

  vector<ftEdge*> edges;
  edges.push_back(edg);
  edg->getVertices(v1, v2);

  // Traverse backwards from start of given edge
  ftEdge* curr = edg;
  // Fetch edges in the underlying surface
  vector<ftEdge*> curr_edges = v1->uniqueEdges();
  int ki;
  for (ki=(int)curr_edges.size()-1; ki>=0; --ki)
    {
      shared_ptr<ParamSurface> surf1 = curr_edges[ki]->face()->surface();
      shared_ptr<ParamSurface> surf2;
      surf2 = (curr_edges[ki]->twin()) ? 
	curr_edges[ki]->twin()->geomEdge()->face()->surface() : surf1;
      if (surf1.get() != surf.get() && surf2.get() != surf.get())
	curr_edges.erase(curr_edges.begin()+ki);
    }
  while (curr_edges.size() == 2)
    {
      // Check if the vertex represents a corner
      ftEdge *other = (curr_edges[0] == curr || curr_edges[0]->twin() == curr) ?
	curr_edges[1] : curr_edges[0];
      double t1 = curr->parAtVertex(v1.get());
      double t2 = other->parAtVertex(v1.get());
      Point tan1 = curr->tangent(t1);
      Point tan2 = other->tangent(t2);
      double ang = tan1.angle(tan2);
      if (ang > angtol)
	break; // Corner

      edges.insert(edges.begin(), other);
      v1 = other->getOtherVertex(v1.get());
      curr = other;

      if (v1.get() == v2.get())
	break;

      curr_edges = v1->uniqueEdges();
      for (ki=(int)curr_edges.size()-1; ki>=0; --ki)
	{
	  shared_ptr<ParamSurface> surf1 = curr_edges[ki]->face()->surface();
	  shared_ptr<ParamSurface> surf2;
	  surf2 = (curr_edges[ki]->twin()) ? 
	    curr_edges[ki]->twin()->geomEdge()->face()->surface() : surf1;
	  if (surf1.get() != surf.get() && surf2.get() != surf.get())
	    curr_edges.erase(curr_edges.begin()+ki);
	}
    }
      
  // Traverse forwards
  curr = edg;
  curr_edges = v2->uniqueEdges();
  for (ki=(int)curr_edges.size()-1; ki>=0; --ki)
    {
      shared_ptr<ParamSurface> surf1 = curr_edges[ki]->face()->surface();
      shared_ptr<ParamSurface> surf2;
      surf2 = (curr_edges[ki]->twin()) ? 
	curr_edges[ki]->twin()->geomEdge()->face()->surface() : surf1;
      if (surf1.get() != surf.get() && surf2.get() != surf.get())
	curr_edges.erase(curr_edges.begin()+ki);
    }
  while (curr_edges.size() == 2)
    {
      // Check if the vertex represents a corner
      vector<ftEdge*> curr_edges = v2->uniqueEdges();
      ftEdge *other = (curr_edges[0] == curr || curr_edges[0]->twin() == curr) ?
	curr_edges[1] : curr_edges[0];
      double t1 = curr->parAtVertex(v2.get());
      double t2 = other->parAtVertex(v2.get());
      Point tan1 = curr->tangent(t1);
      Point tan2 = other->tangent(t2);
      double ang = tan1.angle(tan2);
      if (ang > angtol)
	break; // Corner

      edges.push_back(other);
      v2 = other->getOtherVertex(v2.get());
      curr = other;

      if (v2.get() == v1.get())
	break;

      curr_edges = v2->uniqueEdges();
      for (ki=(int)curr_edges.size()-1; ki>=0; --ki)
	{
	  shared_ptr<ParamSurface> surf1 = curr_edges[ki]->face()->surface();
	  shared_ptr<ParamSurface> surf2;
	  surf2 = (curr_edges[ki]->twin()) ? 
	    curr_edges[ki]->twin()->geomEdge()->face()->surface() : surf1;
	  if (surf1.get() != surf.get() && surf2.get() != surf.get())
	    curr_edges.erase(curr_edges.begin()+ki);
	}
    }
 
  return edges;
}

}   // namespace Go
