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

#include <algorithm>
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"
#include <fstream>

#define DEBUG

namespace Go
{

//===========================================================================
ftSamplePoint::ftSamplePoint(Vector3D xyz, int bnd)
    : xyz_(xyz), uv_(Vector2D(0.0, 0.0)),
      dist_(-1.0), index_(0), at_boundary_(bnd)
//---------------------------------------------------------------------------
//
// Purpose: Constructor
//
//===========================================================================
{
}


//===========================================================================
void ftSamplePoint::addNeighbour(PointIter next)
//---------------------------------------------------------------------------
//
// Purpose: Add a new neighbouring point
//
//===========================================================================
{
    size_t ki;
    for (ki=0; ki<next_.size(); ++ki)
	if (next_[ki] == next)
	    break;
    if (ki == next_.size())
	next_.push_back(next);
}

//===========================================================================
void ftSamplePoint::removeNeighbour(PointIter neighbour)
//===========================================================================
{
    for (size_t ki = 0; ki < next_.size(); ++ki) {
	if (neighbour == next_[ki]) {
	    next_.erase(next_.begin() + ki);
	    break;
	}
    }
}

/// Function object for isOnBoundary test.
class notBoundaryPredicate
{
public:
    bool operator() (PointIter p)
    {
	return !(p->isOnBoundary());
    }
};                



//===========================================================================
void ftSamplePoint::orderNeighbours(ftSamplePoint* nextpoint, bool forward)
//===========================================================================
{
    int n = (int)next_.size();
    // Sort it so that the boundary elements are last
    std::vector<PointIter>::iterator new_last
	= std::partition(next_.begin(), next_.end(),
			 //mem_fun(&ftSamplePoint::isOnBoundary));
			 notBoundaryPredicate());
    ALWAYS_ERROR_IF(next_.end() - new_last != 2,
		    "Boundary nodes must have precisely two boundary-node "
		    "neighbours");

    // Now, swap the second-last and the first element
    std::swap(next_[n-2], next_[0]);
	// Now, we may have to swap the first and last elements in order to
	// satisfy our requirements:
    if (&(*next_[0]) != nextpoint  && forward) {
      std::swap(next_[n-1], next_[0]);
    } else if (&(*next_[0]) == nextpoint && !forward) {
      std::swap(next_[n-1], next_[0]);
    }
}

//===========================================================================
double ftSamplePoint::pntDist(ftSamplePoint* other) const
//===========================================================================
{
    return xyz_.dist(other->xyz_);
}

//===========================================================================
  void ftSamplePoint::getAttachedTriangles(vector<vector<int> >& triangles) const
//===========================================================================
{
  PointIter pnt;
  for (size_t ki=0; ki<next_.size(); ++ki)
    for (size_t kj=0; kj<next_[ki]->next_.size(); ++kj)
      {
	pnt = next_[ki]->next_[kj];
	size_t kh = 0;
	for (kh=0; kh<pnt->next_.size(); ++kh)
	  if (pnt->next_[kh] == this)
	    {
	      // A triangle is found. Arrange points after increasing
	      // index
	      vector<int> index(3);
	      index[0] = index_;
	      index[1] = next_[ki]->index_;
	      index[2] = pnt->index_;
	      std::sort(index.begin(), index.end());

	      triangles.push_back(index);
	      break;
	    }
	if (kh<pnt->next_.size())
	    break;
      }
	      
}

//===========================================================================
void ftSamplePoint::write2Dval(std::ostream& os) const
//===========================================================================
{
  uv_.write(os);
}

//===========================================================================
void ftPointSet::orderNeighbours()
//===========================================================================
{
    ALWAYS_ERROR_IF((first_ == 0) || (second_ == 0),
		"Direction pointers have not been set!");
    if (!first_->isOnBoundary()) {
	return;
    }
    // The neighbours must be in the correct order!
    // This implies that the first and last elements of neighbours
    // must be the neighbouring boundary elements, if this is a boundary
    // element. The rest may be random. The first must be the CCW boundary
    // element and the last the CW one.
    first_->orderNeighbours(&(*second_), true);

    // Now, I'll do the same to all neighbours around the whole boundary.
    // At least two nodes are assumed to exist.
    bool finished = false;
    ftSamplePoint* last = &(*first_);
    ftSamplePoint* current = &(*first_->getFirstNeighbour());
    while (!finished) {
      current->orderNeighbours(last, false);
      last = current;
	current = &(*current->getFirstNeighbour());
	if (current == &(*first_))
	    finished = true;
    }
    nb_ordered_ = true;
}



//===========================================================================
ftPointSet::ftPointSet()
//---------------------------------------------------------------------------
//
// Purpose: Constructor
//
//===========================================================================
    : nb_ordered_(false)
{ }

//===========================================================================
ftPointSet::~ftPointSet()
//---------------------------------------------------------------------------
//
// Purpose: Destructor
//
//===========================================================================
{ }

//===========================================================================
void ftPointSet::removePoint(PointIter point)
//===========================================================================
{
  int idx = point->getIndex();

  // First find point instances

  std::list<shared_ptr<ftSamplePoint> >::iterator pnt = points_.begin();
  std::list<shared_ptr<ftSamplePoint> >::iterator begin = points_.begin();
  std::list<shared_ptr<ftSamplePoint> >::iterator end = points_.end();
  
  if (idx < int(index_to_iter_.size()/2))
    {
      for (; pnt!=end; pnt++)
	{
	  if ((*pnt).get() == point)
	    break;
	}

      if (pnt == end)
	return;         // Point not found
    }
    else
      {
	for (pnt=end; pnt!=begin;)
	{
	  pnt--;
	  if ((*pnt).get() == point)
	    break;
	}

	if (pnt == begin && (*pnt).get() != point)
	  return;         // Point not found
      }

    // Remove neighbour links
    vector<PointIter> neighbours = point->getNeighbours();
    for (size_t kj=0; kj<neighbours.size(); ++kj)
	neighbours[kj]->removeNeighbour(point);

    points_.erase(pnt);
    index_to_iter_.erase(index_to_iter_.begin()+idx);

    int ki;
    for (ki=idx; ki<(int)index_to_iter_.size(); ++ki)
      index_to_iter_[ki]->setIndex(ki);
}

//===========================================================================
double ftPointSet::getMaxDist() const
//---------------------------------------------------------------------------
//
// Purpose: Return the maximum distance in the pointset
//
//===========================================================================
{
  double maxdist = -1.0;
  double dist;

  for (PointList::const_iterator it
	 =points_.begin();
       it != points_.end(); it++)
  {
    dist = (*it)->getDist();
    if (dist < 0.0)
      return -1.0;
    else
      maxdist = std::max(maxdist, dist);
  }

  return maxdist;
}

//===========================================================================
double ftPointSet::getMeanDist() const
//---------------------------------------------------------------------------
//
// Purpose: Return the medium distance in the pointset
//
//===========================================================================
{
  double dist = 0.0;
  for (PointList::const_iterator it=points_.begin();
       it != points_.end(); it++)
    dist += (*it)->getDist();
  dist /= (double)(index_to_iter_.size());

  return dist;
}

//===========================================================================
double ftPointSet::reparBdy(shared_ptr<ParamSurface> surf, bool use_seed)
//---------------------------------------------------------------------------
//
// Purpose: Compute the distances from the points in the point set
//          to the given surface and reparametrize.
//
//===========================================================================
{
  // Very small: We should not move very far on the surface,
  // so it should be OK.
  double eps = 1e-13;
  double clo_u, clo_v, clo_dist;
  Point pt(3);
  Point clo_pt(3);
  RectDomain *rd = 0;
  Vector2D par;
  double maxdist = 0.0;
  for (size_t ki=0; ki<index_to_iter_.size(); ki++) {
    // @@@ We should make sure the existing parameter values
    // are used as starting values for the closest point
    // iteration
      if (this->operator[]((int)ki)->isOnBoundary())
      {
	double *seed = NULL;
	double *tmp = NULL;
	par = this->operator[]((int)ki)->getPar();
	  tmp = &par[0];
	if (use_seed) {
	    par = this->operator[]((int)ki)->getPar();
	  seed = &par[0];
	}
	pt.setValue(this->operator[]((int)ki)->getPoint().begin());

	surf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist,
				   eps, rd, seed);
	maxdist = std::max(maxdist, clo_dist);

	this->operator[]((int)ki)->setPar(Vector2D(clo_u, clo_v));
	this->operator[]((int)ki)->setDist(clo_dist);
      }
  }
  return maxdist;
}


//===========================================================================
double ftPointSet::reparInnerPoints(shared_ptr<ParamSurface> surf, bool use_seed)
//---------------------------------------------------------------------------
//
// Purpose: Compute the distances from the points in the point set
//          to the given surface and reparametrize.
//
//===========================================================================
{
  // Very small: We should not move very far on the surface,
  // so it should be OK.
  double eps = 1e-13;
  double u, v, dist;
  Point pt(3);
  Point clo_pnt(3);
  double maxdist = 0.0;
  for (size_t ki=0; ki<index_to_iter_.size(); ki++) {
    // @@@ We should make sure the existing parameter values
    // are used as starting values for the closest point
    // iteration
      if (!(this->operator[]((int)ki)->isOnBoundary()))
	{
	    pt.setValue(this->operator[]((int)ki)->getPoint().begin());
	  // 	surf->closestPoint(pt, u, v, clo_pnt, dist, eps);

	  double *seed = NULL;
	  double *tmp = NULL;
	  Vector2D par = this->operator[]((int)ki)->getPar();
	  tmp = &par[0];
	  if (use_seed) {
	    seed = &par[0]; //u;
	  }
	  RectDomain *rd = 0;
	  surf->closestPoint(pt, u, v, clo_pnt, dist, eps, rd, seed);
	  maxdist = std::max(maxdist, dist);

	  // We check whether the new point was closer.
	  Point uv_pt = surf->ParamSurface::point(par[0], par[1]);
	  double curr_dist = pt.dist(uv_pt);
	  if (curr_dist < dist) {
	      this->operator[]((int)ki)->setDist(curr_dist);
	      std::cout << "Curr dist: " << this->operator[]((int)ki)->getDist() <<
	      ", proj dist: " << dist << std::endl;
	  } else {
	      this->operator[]((int)ki)->setPar(Vector2D(u, v));
	      this->operator[]((int)ki)->setDist(dist);
	  }
	}
  }
  return maxdist;
}

//===========================================================================
void ftPointSet::computeDist(shared_ptr<ParamSurface> surf)
//---------------------------------------------------------------------------
//
// Purpose: Compute the distances from the points in the point set
//          to the given surface.
//
//===========================================================================
{
  // Very small: We should not move very far on the surface,
  // so it should be OK.
  double eps = 1e-13;
  double u, v, dist;
  Point pt(3);
  Point clo_pnt(3);
  double seed[2];

  for (PointList::iterator it=points_.begin();
       it != points_.end(); it++)
    {
      Vector2D par = (*it)->getPar();
      pt.setValue((*it)->getPoint().begin());
      seed[0] = par[0];
      seed[1] = par[1];
      if ((*it)->isOnBoundary())
	  surf->closestBoundaryPoint(pt, u, v, clo_pnt, dist, eps, NULL, seed);
      else
	  surf->closestPoint(pt, u, v, clo_pnt, dist, eps, NULL, seed);

      (*it)->setDist(dist);
  }
}

//===========================================================================
void ftPointSet::computeParametricDist(shared_ptr<ParamSurface> surf)
//---------------------------------------------------------------------------
//
// Purpose: Compute the distances from the points in the point set
//          to the given surface.
//
//===========================================================================
{
  // Very small: We should not move very far on the surface,
  // so it should be OK.
//   double eps = 1e-13;
//   double u, v, dist;
//   Point clo_pnt(3);
//   double seed[2];
    Point pt(3);

    for (PointList::iterator it=points_.begin();
	 it != points_.end(); it++) {
	Vector2D par = (*it)->getPar();
	surf->point(pt, par[0], par[1]);
	Vector3D pos(pt.begin());
	(*it)->setDist((*it)->getPoint().dist(pos));
    }
}


//===========================================================================
void ftPointSet::computeDistAndRepar(shared_ptr<ParamSurface> surf)
//---------------------------------------------------------------------------
//
// Purpose: Compute the distances from the points in the point set
//          to the given surface and reparametrize.
//
//===========================================================================
{
  // Very small: We should not move very far on the surface,
  // so it should be OK.
  double eps = 1e-13;
  int ki;
  double u, v, dist;
  Point pt(3);
  Point clo_pnt(3);
  double seed[2];
  ki = 0;
  for (PointList::iterator it = points_.begin(); it != points_.end(); it++) {
      ki++;
      Vector2D par = (*it)->getPar();
      pt.setValue((*it)->getPoint().begin());
      u = par[0];
      v = par[1];
      surf->point(clo_pnt, u, v);
      double curr_dist = pt.dist(clo_pnt);
      seed[0] = par[0];
      seed[1] = par[1];
      try {
	if ((*it)->isOnBoundary())
	  surf->closestBoundaryPoint(pt, u, v, clo_pnt, dist, eps, NULL, seed);
	else
	  surf->closestPoint(pt, u, v, clo_pnt, dist, eps, NULL, seed);
      }
      catch(...)
	{
	  dist = curr_dist;
	}
      if (dist < curr_dist) {
	(*it)->setPar(Vector2D(u, v));
	(*it)->setDist(dist);
      } else // We do not reparametrize if found point is worse than the current.
	(*it)->setDist(curr_dist);

//       // debugging
//       if (dist > 1e-01) {
// 	std::cout << ki << " (" << u << "," << v << "), (";
// 	std::cout << pt[0] << "," << pt[1] << "," << pt[2] << "), ";
// 	std::cout << dist << " " << boundary_pt << std::endl;
//       }
//       // end of debugging
  }
}

//===========================================================================
void ftPointSet::append(shared_ptr<ftPointSet> triang)
//---------------------------------------------------------------------------
//
// Purpose: Extend the point set with points from another set. Avoid points with no
//          connectivity
//
//===========================================================================
{
    std::list<shared_ptr<ftSamplePoint> >::iterator pnt = triang->points_.begin();
    std::list<shared_ptr<ftSamplePoint> >::iterator end = triang->points_.end();
    for (; pnt!=end; pnt++)
    {
	int nmb = (*pnt)->getNmbNeighbour();
	if (nmb == 1)
	{
	    // Too low connectivity. Remove pointer to this point
	    vector<ftSamplePoint*> neighbours = (*pnt)->getNeighbours();
	    neighbours[0]->removeNeighbour((*pnt).get());
	}
	
	if (nmb <= 1)
	    continue;  // Point dismissed

	(void)addEntry(*pnt);
    }
	
}

  //===========================================================================
    bool compare_par(pair<ftSamplePoint*,double> f1, pair<ftSamplePoint*,double> f2)
    {
	return (f1.second < f2.second);
    }

//===========================================================================
void ftPointSet::cleanNodeIdentity(double tol)
//---------------------------------------------------------------------------
//
// Purpose: Remove identitical boundary nodes
//
//===========================================================================
{
    // Fetch boundary nodes
    vector<ftSamplePoint*> bd_nodes;
    size_t ki, kj;
    for (ki=0; ki<index_to_iter_.size(); ++ki)
	if (index_to_iter_[ki]->isOnSubSurfaceBoundary())
	    bd_nodes.push_back(index_to_iter_[ki]);

    // Compare distances. NB! Some sorting may speed up this process
    for (ki=0; ki<bd_nodes.size(); ki++)
	for (kj=ki+1; kj<bd_nodes.size(); kj++)
	    if (bd_nodes[ki]->pntDist(bd_nodes[kj]) < tol)
	    {
		ftSurfaceSetPoint* pnt1 = bd_nodes[ki]->asSurfaceSetPoint();
		ftSurfaceSetPoint* pnt2 = bd_nodes[kj]->asSurfaceSetPoint();
		if (!(pnt1 && pnt2))
		    continue;
		pnt1->addInfo(pnt2);
		removePoint(bd_nodes[kj]);
		bd_nodes.erase(bd_nodes.begin()+kj);
		kj--;
	    }
}

//===========================================================================
void ftPointSet::mergeBoundary(shared_ptr<ftFaceBase> face1, int range1_idx1, 
			       int range1_idx2, shared_ptr<ftFaceBase> face2,
			       int range2_idx1, int range2_idx2, double eps)

//===========================================================================
{
  // Identify twin edges related to the two faces
  vector<shared_ptr<ftEdgeBase> > edges1 = face1->createInitialEdges();
  vector<shared_ptr<ftEdgeBase> > edges2;
  size_t ki, kj;
  double close_dist = 1.0e-4; //1.0e-3; //0.05;   // Must be set appropriately
  for (ki=0; ki<edges1.size(); ++ki)
    if (edges1[ki]->twin() && edges1[ki]->twin()->face() == face2.get())
      edges2.push_back(edges1[ki]);

  // Join edges into smooth curves to avoid unnecesary fractioning
  vector<shared_ptr<ParamCurve> > bd_crvs;
  mergeBoundaryEdges(edges2, bd_crvs, eps);

  // Compute endpoints of merging curve(s)
  vector<Point> corners(2*bd_crvs.size());
  for (ki=0; ki<bd_crvs.size(); ++ki)
    {
      corners[2*ki] = bd_crvs[ki]->point(bd_crvs[ki]->startparam());
      corners[2*ki+1] = bd_crvs[ki]->point(bd_crvs[ki]->endparam());
    }

  // Dismiss identical points
  for (ki=0; ki<corners.size(); ++ki)
    for (kj=ki+1; kj<corners.size(); )
      {
	if (corners[ki].dist(corners[kj]) < close_dist)
	  corners.erase(corners.begin()+kj);
	else
	  kj++;
      }
      
  for (ki=0; ki<bd_crvs.size(); ++ki)
    {
      // For all relevant boundary points, check if they lie along this edge
      // and sort them with respect to the edge
      vector<pair<ftSamplePoint*,double> > points_on_edge;

      shared_ptr<ParamCurve> tmp_crv = bd_crvs[ki];
#ifdef DEBUG
      std::ofstream pt("debug_point.g2");
      std::ofstream out("debug_edge.g2");
      tmp_crv->geometryCurve()->writeStandardHeader(out);
      tmp_crv->geometryCurve()->write(out);
#endif

      int kr, kh, idx;
      int last_idx = std::min((int)index_to_iter_.size(), range2_idx2);
      for (kr=range1_idx1; kr<last_idx; ++kr)
	{
	  if (kr >= range1_idx2 && kr<range2_idx1)
	    continue;  // Not a relevant point

	  ftSamplePoint *curr = (*this)[kr];
	  if (!curr->isOnBoundary())
	    continue;  // Not a boundary point

	  Vector3D pos = curr->getPoint();
	  Point pos2(pos.begin(), pos.end());
	  Point close;
	  double par, dist;
	  tmp_crv->closestPoint(pos2, tmp_crv->startparam(), tmp_crv->endparam(),
				par, close, dist);
#ifdef DEBUG
	  pt << "400 1 0 4 255 0 0 255 " << std::endl;
	  pt << "1" << std::endl;
	  pt << pos2[0] << " " << pos2[1] << " " << pos2[2] << std::endl;
	  pt << "400 1 0 4 0 255 0 255 " << std::endl;
	  pt << "1" << std::endl;
	  pt << close[0] << " " << close[1] << " " << close[2] << std::endl;
#endif

	  if (dist < close_dist)
	    points_on_edge.push_back(std::make_pair(curr, par));
	}

      if (points_on_edge.size() < 2)
	return;  // Nothing to do

      // Sort along edge
      std::sort(points_on_edge.begin(), points_on_edge.end(), compare_par);

      // Check if the edge is complete
      for (idx=1; idx<(int)points_on_edge.size(); idx++)
	if ((points_on_edge[idx].first->containsFace(face1.get()) && 
	     points_on_edge[0].first->containsFace(face2.get())) ||
	    (points_on_edge[idx].first->containsFace(face2.get()) && 
	     points_on_edge[0].first->containsFace(face1.get())))
	  break;

      if (idx >= (int)points_on_edge.size())
	idx = -1;
      else if (points_on_edge[0].second > tmp_crv->startparam() + eps)
	idx = 0;
      else if (points_on_edge[idx].second > tmp_crv->startparam() + eps);
      else
	idx = -1;
	 
      if (idx >= 0)
	{
	  // Check if any neighbouring points lie on the edge
	  vector<PointIter> next = points_on_edge[idx].first->getNeighbours();
	  for (kj=0; kj<next.size(); ++kj)
	    {
	      if (!next[kj]->isOnBoundary())
		continue;  // Not a boundary point

	      for (kh=0; kh<(int)points_on_edge.size(); ++kh)
		if (next[kj] == points_on_edge[kh].first)
		  break;
	      if (kh < (int)points_on_edge.size())
		continue;   // Point exist already

	      Vector3D pos = next[kj]->getPoint();
	      Point pos2(pos.begin(), pos.end());
	      Point close;
	      double par, dist;
	      tmp_crv->closestPoint(pos2, tmp_crv->startparam(), tmp_crv->endparam(),
				    par, close, dist);
	      if (dist < close_dist)
		points_on_edge.insert(points_on_edge.begin(), std::make_pair(next[kj], par));
	    }
	}

      for (idx=(int)points_on_edge.size()-2; idx>=(int)0; idx--)
	if ((points_on_edge[idx].first->containsFace(face1.get()) && 
	     points_on_edge[(int)points_on_edge.size()-1].first->containsFace(face2.get())) ||
	    (points_on_edge[idx].first->containsFace(face2.get()) && 
	     points_on_edge[points_on_edge.size()-1].first->containsFace(face1.get())))
	  break;

      if (idx < 0)
	idx = -1;
      else if (points_on_edge[points_on_edge.size()-1].second < tmp_crv->endparam() - eps)
	idx = (int)points_on_edge.size()-1;
      else if (points_on_edge[idx].second < tmp_crv->endparam() - eps);
      else
	idx = -1;
	 
      if (idx >= 0)
	{
	  // Check if any neighbouring points lie on the edge
	  vector<PointIter> next = points_on_edge[idx].first->getNeighbours();
	  for (kj=0; kj<next.size(); ++kj)
	    {
	      if (!next[kj]->isOnBoundary())
		continue;  // Not a boundary point

	      for (kh=0; kh<(int)points_on_edge.size(); ++kh)
		if (next[kj] == points_on_edge[kh].first)
		  break;
	      if (kh < (int)points_on_edge.size())
		continue;   // Point exist already

	      Vector3D pos = next[kj]->getPoint();
	      Point pos2(pos.begin(), pos.end());
	      Point close;
	      double par, dist;
	      tmp_crv->closestPoint(pos2, tmp_crv->startparam(), tmp_crv->endparam(),
				    par, close, dist);
	      if (dist < close_dist)
		points_on_edge.push_back(std::make_pair(next[kj], par));
	    }
	}

      // Sort again
      std::sort(points_on_edge.begin(), points_on_edge.end(), compare_par);

      // Alter boundary classification of points being distant from the
      // endpoints of the edge to inner edge
      for (kr=0; kr<(int)points_on_edge.size(); ++kr)
	{
	  Vector3D tmp_pt1 = points_on_edge[kr].first->getPoint();
	  Point tmp_pt2(tmp_pt1.begin(), tmp_pt1.end());
	  points_on_edge[kr].first->setBoundary(2);
	  for (kj=0; kj<corners.size(); ++kj)
	    {
	      double dist = corners[kj].dist(tmp_pt2);
	      if (dist < close_dist)
		points_on_edge[kr].first->setBoundary(1);
	    }
	}

#ifdef DEBUG
	     std::ofstream debug0("debug0.g2");
	   this->write(debug0);
#endif
	
	   // Merge boundary information
	   ftSurfaceSetPoint *pnt0=0, *pnt1=0, *pnt2=0, *pnt3=0, *pnt4=0;
	   pnt1 = points_on_edge[0].first->asSurfaceSetPoint();
	   if (points_on_edge.size() > 1)
	     pnt2 = points_on_edge[1].first->asSurfaceSetPoint();
	   if (points_on_edge.size() > 2)
	     pnt3 = points_on_edge[2].first->asSurfaceSetPoint();
	   size_t index = 3;
	   size_t del;
	   while (pnt1)
	     {
#ifdef DEBUG
	       std::ofstream outp("debug_p.txt");
	       printPoints(outp);
#endif

	       pnt4 = (index < points_on_edge.size()) ? 
		 points_on_edge[index].first->asSurfaceSetPoint() : 0;

#ifdef DEBUG 
	       std::ofstream outp2("bd_pnts.g2");
	       outp2 << "400 1 0 4 100 0 155 255" << std::endl;
	       outp2 << "1" << std::endl;
	       outp2 << pnt1->getPoint() << std::endl;
	       if (pnt2)
		 {
		   outp2 << "400 1 0 4 100 0 155 255" << std::endl;
		   outp2 << "1" << std::endl;
		   outp2 << pnt2->getPoint() << std::endl;
		 }
	       if (pnt3)
		 {
		   outp2 << "400 1 0 4 100 0 155 255" << std::endl;
		   outp2 << "1" << std::endl;
		   outp2 << pnt3->getPoint() << std::endl;
		 }
	       if (pnt4)
		 {
		   outp2 << "400 1 0 4 100 0 155 255" << std::endl;
		   outp2 << "1" << std::endl;
		   outp2 << pnt4->getPoint() << std::endl;
		 }
#endif

	       if (!pnt2 || pnt1->pntDist(pnt2) >= eps)
		 {
		   if ((pnt0 == 0 || pnt2 == 0 || (pnt3 && pnt2->pntDist(pnt3) < eps) ||
			(pnt1->containsFace(face1.get()) && pnt2->containsFace(face1.get())) ||
			(pnt1->containsFace(face2.get()) && pnt2->containsFace(face2.get()))))
		     {
		       // Check if the point already contains all relevant information
		       if (!(pnt1->containsFace(face1.get()) && pnt1->containsFace(face2.get())))
			 {
			   // Point one must be kept, add information related to the other face
			   // Face information
			   shared_ptr<ftFaceBase> cf = (pnt1->containsFace(face1.get())) ? face2 : face1;
			   shared_ptr<ftFaceBase> cf2 = (pnt1->containsFace(face1.get())) ? face1 : face2;
			   pnt1->addFace(cf);
		
			   // Triangle information
			   // Find first point along the edge related to the other face
			   for (kr=(int)index-2; kr<(int)points_on_edge.size(); ++kr)
			     if (points_on_edge[kr].first->containsFace(cf.get()))
			       {
				 addConnectivityInfo(pnt1, points_on_edge[kr].first, 0);
				 break;
			       }
			   if (kr >= (int)points_on_edge.size() && pnt0)
			     addConnectivityInfo(pnt1, pnt0, cf2.get());
			 }
		       del = 1;
		     }
		   else if (pnt3 == 0)
		     {
		       // Move point 1 to point 2
		       // Find new parameter value
		       Vector3D pos = pnt2->getPoint();
		       int bd = (pnt2->isOnBoundary()) ? 1 :
			 ((pnt2->isOnSubSurfaceBoundary() ? 2 : 0));
		       pnt1->resetPosition(pos, bd);
		     }
		   else
		     {
		       // Position the new point between point 1 and point 2
		       double t1 = 0.5*(points_on_edge[index-3].second + points_on_edge[index-2].second);
		       Point pos = tmp_crv->point(t1);
		       Vector3D pos2(pos.begin());
		       int bd = 2;
		       for (kr=0; kr<(int)corners.size(); ++kr)
			 if (corners[kr].dist(pos) < close_dist)
			   bd = 1;
		       pnt1->resetPosition(pos2,bd);
		       pnt2->resetPosition(pos2,bd);
		     }
		 }	
					 
					
	       if (pnt2 && pnt1->pntDist(pnt2) < eps)
		 {
		   if (pnt1->containsFace(face1.get()))
		     {
		       // Keep point 1, remove point 2 after transferring information
		       pnt1->addInfo(pnt2);
		       removePoint(pnt2);
		       pnt0 = pnt1;
		     }
		   else
		     {
		       // Keep point 2, remove point 1 after transferring information
		       pnt2->addInfo(pnt1);
		       removePoint(pnt1);
		       pnt0 = pnt2;
		     }
		   del = 2;
		 }

#ifdef DEBUG
	       std::ofstream debug("debug.g2");
	       this->write(debug);
#endif
	
	       // Update pointers
	       if (del == 1)
		 {
		   pnt0 = pnt1;
		   pnt1 = pnt2;
		   pnt2 = pnt3;
		   pnt3 = pnt4;
		   index++;
		 }
	       else
		 {
		   pnt2 = pnt4;
		   pnt1 = pnt3;
		   pnt3 = (index < points_on_edge.size()-1) ? 
		     points_on_edge[index+1].first->asSurfaceSetPoint() : 0;
		   index += del;
		 }
	    
	     }
	     }
}

 //===========================================================================
void  ftPointSet::addConnectivityInfo(PointIter pnt, PointIter pnt2, ftFaceBase* other_face)
//===========================================================================
{
    // Select neighbouring point to pnt2
    vector<PointIter> neighbours = pnt2->getNeighbours();

    // Find the boundary neighbouring point
    PointIter pnt3 = 0;
    PointIter pnt4 = 0;
    size_t ki, kj;
    double max_ang = 0.0;
    int max_idx = -1;
    Vector3D xyz1 = pnt->getPoint();
    Vector3D xyz2 = pnt2->getPoint();
    Point pos(xyz1.begin(), xyz1.end());
    Point pos2(xyz2.begin(), xyz2.end());
    Point vec1 = pos2 - pos;
    for (kj=0; kj<neighbours.size(); ++kj)
    {
	if (other_face && neighbours[kj]->containsFace(other_face))
	    continue;

	if (neighbours[kj]->isOnBoundary() || 
	    neighbours[kj]->isOnSubSurfaceBoundary())
	{
	    Vector3D xyz3 = neighbours[kj]->getPoint();
	    Point pos3(xyz3.begin(), xyz3.end());
	    Point vec2 = pos3 - pos2;
	    double ang1 = vec1.angle(vec2);
	    if (ang1 > max_ang)
	    {
		max_ang = ang1;
		max_idx = (int)kj;
	    }
	}
    }
    pnt3 = (max_idx >= 0) ? neighbours[max_idx] : 0;

    if (pnt3)
    {
	if (other_face && pnt2->containsFace(other_face))
	    pnt4 = 0;
	else
	{
	    // Find the third point in the triangle containing pnt2 and pnt3
	    vector<PointIter> neighbours2 = pnt3->getNeighbours();
	    for (kj=0; kj<neighbours2.size(); ++kj)
	    {
		for (ki=0; ki<neighbours.size(); ++ki)
		{
		    if (neighbours2[kj] == neighbours[ki] && neighbours[ki] != pnt2 &&
			neighbours[ki] != pnt3)
		    {
			pnt4 = neighbours[ki];
			break;
		    }
		}
		if (pnt4)
		    break;
	    }
	}

	if (pnt4)
	{
	    // Split long triangle
	    pnt3->removeNeighbour(pnt2);
	    pnt2->removeNeighbour(pnt3);

	    // Make connections
	    pnt->addNeighbour(pnt3);
	    pnt3->addNeighbour(pnt);
	    pnt->addNeighbour(pnt4);
	    pnt4->addNeighbour(pnt);
	    pnt->addNeighbour(pnt2);
	    pnt2->addNeighbour(pnt);
	}
	else
	{
	    pnt->addNeighbour(pnt3);
	    pnt3->addNeighbour(pnt);
	}
    }
		
}

//===========================================================================
void ftPointSet::mergeBoundaryEdges(vector<shared_ptr<ftEdgeBase> >& edges,
				    vector<shared_ptr<ParamCurve> >& crvs,
				    double tol) const
//===========================================================================
{
  double ang_tol = 0.02; // 0.01;
    shared_ptr<ParamCurve> prev;
    shared_ptr<ParamCurve> curr;

    size_t ki, kj;
    for (ki=0, kj=0; ki<edges.size(); ++ki)
    {
	curr = edges[ki]->geomEdge()->geomCurve();
	if (prev.get() == 0)
	{
	    crvs.push_back(shared_ptr<ParamCurve>(curr->clone()));
	}
	else
	{
	    // Check continuity
	    vector<Point> der1(2);
	    vector<Point> der2(2);
	    crvs[kj]->point(der1, crvs[kj]->endparam(), 1);
	    curr->point(der2, curr->startparam(), 1);
	    if (der1[0].dist(der2[0]) < tol && der1[1].angle(der2[1]) < ang_tol)
	    {
		// Sufficient continuity
		crvs[kj]->appendCurve(curr->clone());
	    }
	    else
	    {
		// Check continuity in start of chain
		crvs[0]->point(der1, crvs[0]->startparam(), 1);
		curr->point(der2, curr->endparam(), 1);
		if (der1[0].dist(der2[0]) < tol && der1[1].angle(der2[1]) < ang_tol)
		{
		    // Sufficient continuity. Turn curves and append
		    crvs[0]->reverseParameterDirection();
		    ParamCurve *curr2 = curr->clone();
		    curr2->reverseParameterDirection();
		    crvs[0]->appendCurve(curr2);
		}
		else
		{
		    // Not sufficient continuity. Make new boundary curve
		    crvs.push_back(shared_ptr<ParamCurve>(curr->clone()));
		    kj++;
		}
	    }
	}
	prev = curr;
    }	
}

//===========================================================================
void ftPointSet::identifyBdPnts(vector<Point>& points, vector<int>& pnt_ix)
//===========================================================================
{
  pnt_ix.resize(points.size());

  // For each point
  size_t nmb = index_to_iter_.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      // Traverse all boundary points and select the closest
      double min_dist = HUGE;
      int min_idx = -1;

      for (size_t kj=0; kj<nmb; ++kj)
	{
	  ftSamplePoint *curr = (*this)[(int)kj];
	  if (!curr->isOnBoundary())
	    continue;  // Not a boundary point

	  Vector3D pos = curr->getPoint();
	  Point pos2(pos.begin(), pos.end());
	  double dist = pos2.dist(points[ki]);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      min_idx = (int)kj;
	    }
	}
      pnt_ix[ki] = min_idx;
    }
}

//===========================================================================
// Fetch all triangles in the connectivity graph
void ftPointSet::getTriangles(vector<vector<int> >& triangles) const
//===========================================================================
{
  // For each point, get all triangles and make sure to store them only
  // once
  std::set<vector<int> > all_triangles;
  for (size_t ki=0; ki<index_to_iter_.size(); ++ki)
    {
      vector<vector<int> > tri;
      index_to_iter_[ki]->getAttachedTriangles(tri);
      all_triangles.insert(tri.begin(), tri.end());
    }

  triangles.insert(triangles.end(), all_triangles.begin(), 
		   all_triangles.end());
}


//===========================================================================
// Get the position of all points
void ftPointSet::getPoints(std::vector<Vector3D>& positions) const
//===========================================================================
{
  for (size_t ki=0; ki<index_to_iter_.size(); ++ki)
    {
      positions.push_back(index_to_iter_[ki]->getPoint());
    }
//     std::ofstream out("debug_next.txt");
//     printPoints(out);
}

//===========================================================================
void ftPointSet::checkAndUpdateTriangCorners()
//===========================================================================
{
  for (size_t ki=0; ki<index_to_iter_.size(); ++ki)
    {
      if (!index_to_iter_[ki]->isOnBoundary())
	continue;  // Not a boundary point

      // Count the number of neighbouring boundary points
      vector<PointIter> neighbours = index_to_iter_[ki]->getNeighbours();
      int nmb = 0;
      size_t kk;
      for (kk=0; kk<neighbours.size(); ++kk)
	if (neighbours[kk]->isOnBoundary())
	  nmb++;

      if (nmb == 2 && neighbours.size() == 2 && 
	  neighbours[0]->isConnected(neighbours[1]))
	{
	  // Try to swap a triangle edge to avoid a corner configuration
	  // where one boundary point has more than two boundary neighbours

	  // Check if the two points have a common non-boundary
	  // neighbour
	  vector<PointIter> next1 = neighbours[0]->getNeighbours();
	  vector<PointIter> next2 = neighbours[1]->getNeighbours();
	  
	  size_t kr, kh;
	  for (kr=0; kr<next1.size(); ++kr)
	    {
	      for (kh=0; kh<next2.size(); ++kh)
		{
		  if (next1[kr] == next2[kh] && !next1[kr]->isOnBoundary())
		    {
		      neighbours[0]->removeNeighbour(neighbours[1]);
		      neighbours[1]->removeNeighbour(neighbours[0]);
		      next1[kr]->addNeighbour(index_to_iter_[ki]);
		      index_to_iter_[ki]->addNeighbour(next1[kr]);
		      break;
		    }
		}
	      if (kh < next2.size())
		break;
	    }
	}
    }

  if (false)
    {
  for (size_t ki=0; ki<index_to_iter_.size(); ++ki)
    {
      if (!index_to_iter_[ki]->isOnBoundary())
	continue;  // Not a boundary point

      // Count the number of neighbouring boundary points
      vector<PointIter> neighbours = index_to_iter_[ki]->getNeighbours();
      int nmb = 0;
      size_t kk;
      for (kk=0; kk<neighbours.size(); ++kk)
	if (neighbours[kk]->isOnBoundary())
	  nmb++;
      if (nmb == 3)
	{
	  // Check if one connection can be removed
	  size_t kr, kh;
	  for (kr=0; kr<neighbours.size(); ++kr)
	    {
	      if (!neighbours[kr]->isOnBoundary())
		continue;
	      if (index_to_iter_[ki] == neighbours[kr])
		{
		  index_to_iter_[ki]->removeNeighbour(neighbours[kr]);
		  break;
		}
	    }

	  if (kr == neighbours.size())
	    {
	      for (kr=0; kr<neighbours.size(); ++kr)
		{
		  for (kh=kr+1; kh<neighbours.size(); ++kh)
		    {
		      if (!neighbours[kh]->isOnBoundary())
			continue;
		      if (neighbours[kr]->isConnected(neighbours[kh]))
			{
			  double dist1 = index_to_iter_[ki]->pntDist(neighbours[kr]);
			  double dist2 = index_to_iter_[ki]->pntDist(neighbours[kh]);
			  if (dist1 < dist2)
			    index_to_iter_[ki]->removeNeighbour(neighbours[kh]);
			  else
			    index_to_iter_[ki]->removeNeighbour(neighbours[kr]);
			  break;
			}
		    }
		  if (kh < neighbours.size())
		    break;
		}
	    }
	}
    }
    }	      
}

//===========================================================================
void ftPointSet::getOrientedTriangles(vector<vector<int> >& triangles)
//===========================================================================
{
  // Fetch all triangles
  getTriangles(triangles);

  // For each pair of triangles, check that the edge orientation between
  // two nodes are opposite
  size_t ki, kj;
  for (ki=0; ki<triangles.size(); ++ki)
    {
      int not_swapped = -1;
      for (kj=ki+1; kj<triangles.size(); ++kj)
	{
	  // Check if the two triangles have the same two nodes
	  int ki1=-1, ki2=-1, kj1=-1, kj2=-1;
	  for (int ix1=0; ix1<3; ++ix1)
	    for (int ix2=0; ix2<3; ++ix2)
	      {
		if (triangles[ki][ix1] == triangles[kj][ix2])
		  {
		    if (ki1<0)
		      {
			ki1 = ix1;
			kj1 = ix2;
		      }
		    else
		      {
			ki2 = ix1;
			kj2 = ix2;
		      }
		  }
	      }

	  if (ki2 >= 0)
	    {
	      // Common edge
	      // @@@ VSK, 0214. Uncertain whether this test always will
	      // be correct
	      if ((ki2-ki1 == kj2-kj1 && (ki2-ki1)*(kj2-kj1) > 0) ||
		  (ki1==0 && ki2==2  && (ki2-ki1)*(kj2-kj1)<0) ||
		 (kj1==0 && kj2==2 && (ki2-ki1)*(kj2-kj1)<0))
		{
		  // Same orientation. Swap
		  std::swap(triangles[kj][kj1], triangles[kj][kj2]);

		  if (not_swapped >= 0)
		    {
		      // Reorganize triangles to avoid changing sequence back
		      std::swap(triangles[not_swapped], triangles[kj]);
		      not_swapped++;
		      kj--;
		    }
		}
	    }
	  else if (not_swapped == -1)
	    {
	      not_swapped = (int)kj;
	    }
	}
    }
}

//===========================================================================
void ftPointSet::printPoints(std::ostream& os) const
//===========================================================================
{
  for (size_t ki=0; ki<index_to_iter_.size(); ++ki)
    {
	os << ki << " " << index_to_iter_[ki]->getIndex() << ";   ";
	index_to_iter_[ki]->getPoint().write(os);
	os << std::endl;
	vector<PointIter> neighbours = index_to_iter_[ki]->getNeighbours();
	os << neighbours.size() << ": ";
	for (size_t kj=0; kj<neighbours.size(); ++kj)
	    os << neighbours[kj]->getIndex() << "  ";
	os << std::endl;
	os << std::endl;
    }
}
// From PrOrganizedPoints:
//===========================================================================
int         ftPointSet::getNumNodes() const
//===========================================================================
{
    return size();
}

//===========================================================================
Vector3D  ftPointSet::get3dNode(int i) const
//===========================================================================
{
    return this->operator[](i)->getPoint();
}


//===========================================================================
void        ftPointSet::set3dNode(int i, const Vector3D& p)
//===========================================================================
{
    this->operator[](i)->setPoint(p);
}


//===========================================================================
void        ftPointSet::getNeighbours(int i, vector<int>& neighbours) const
//===========================================================================
{
    if (!nb_ordered_) {
// 	THROW("Neighbours are not ordered yet.");
	// @@sbr This should be fixed, but not now.

      // We make sure that at least the neighbours of current pt are ordered.
      // suppose we could speed things up a bit by ordering local pt only.
	//	orderNeighbours();
//       nb_ordered_ = true;
    }

    neighbours.clear();
    const vector<PointIter>& nb = this->operator[](i)->getNeighbours();
    int n = (int)nb.size();
    for (int ki = 0; ki < n; ++ki) {
	neighbours.push_back(nb[ki]->getIndex());
    }
}


//===========================================================================
bool        ftPointSet::isBoundary(int i) const
//===========================================================================
{
    return this->operator[](i)->isOnBoundary();
}


//===========================================================================
double      ftPointSet::getU(int i) const
//===========================================================================
{
    Vector2D uv = this->operator[](i)->getPar();
    return uv[0];
}

//===========================================================================
double      ftPointSet::getV(int i) const
//===========================================================================
{
    Vector2D uv = this->operator[](i)->getPar();
    return uv[1];
}


//===========================================================================
void        ftPointSet::setU(int i, double u)
//===========================================================================
{
    Vector2D uv = this->operator[](i)->getPar();
    uv[0] = u;
    this->operator[](i)->setPar(uv);
}



//===========================================================================
void ftPointSet::setV(int i, double v)
//===========================================================================
{
    Vector2D uv = this->operator[](i)->getPar();
    uv[1] = v;
    this->operator[](i)->setPar(uv);
}


//-----------------------------------------------------------------------------
void ftPointSet::write(std::ostream& os) const
//-----------------------------------------------------------------------------
{
    for (size_t ix=0; ix<index_to_iter_.size(); ++ix)
    {
	Vector3D p1 = (*this)[(int)ix]->getPoint();
	vector<PointIter> ngh = (*this)[(int)ix]->getNeighbours();
	os << "400 1 0 4 255 0 0 255" << std::endl;
	os << "1" << std::endl;
	p1.write(os);
	os << std::endl;
	if (ngh.size() > 0)
	  {
	    os << "410 1 0 4 0 55 200 255" << std::endl;
	    os << ngh.size() << std::endl;
	  }
	for (size_t iy=0; iy<ngh.size(); ++iy)
	{
	    Vector3D p2 = ngh[iy]->getPoint();
	    p1.write(os);
	    os << "  ";
	    p2.write(os);
	    os << std::endl;
	}
    }
	    
}

//-----------------------------------------------------------------------------
void ftPointSet::write2D(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  for (size_t ix=0; ix<index_to_iter_.size(); ++ix)
    {
      os << "400 1 0 4 255 0 0 255" << std::endl;
      os << "1" << std::endl;
      Vector2D pnt = (*this)[(int)ix]->getPar();
      os << pnt;
      //(*this)[ix]->write2Dval(os);
      os << "  0 " << std::endl;
      vector<PointIter> ngh = (*this)[(int)ix]->getNeighbours();
      if (ngh.size() > 0)
	{
	  os << "410 1 0 4 0 55 200 255" << std::endl;
	  os << ngh.size() << std::endl;
	  for (size_t iy=0; iy<ngh.size(); ++iy)
	    {
	      os << pnt;
	      //(*this)[ix]->write2Dval(os);
	      os << " 0     ";
	      Vector2D pnt2 = ngh[iy]->getPar();
	      os << pnt2;
 	      //ngh[iy]->write2Dval(os);
	      os << "  0 " << std::endl;
	    }
	}
    }
}

} // namespace Go
