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
#include <fstream>

using std::vector;

namespace Go {

//==========================================================================
bool Path::estimateHoleInfo(vector<ftEdge*> edges, Point& centre, 
			    Point& axis, double& radius)
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
  return true;
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



}   // namespace Go
