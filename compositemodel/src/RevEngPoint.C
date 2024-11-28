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

#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Point.h"
#include <vector>

using namespace Go;
using std::vector;
  
//===========================================================================
RevEngPoint::RevEngPoint()
  : ftSamplePoint()
//===========================================================================
{
  int dim = 3;
  eigen1_ = Point(dim);
  eigen2_ = Point(dim);
  eigen3_ = Point(dim);
  LocFuncnormal_ = Point(dim);
  kvecmin_ = Point(dim);
  kvecmax_ = Point(dim);
  lambda1_ = lambda2_ = lambda3_ = -1.0;
  kmin_ = kmax_ = 0.0;
  ptdist_ = avdist_ = 0.0;
  nmb_eigen_ = nmb_locfunc_ = 0;
  Point dummy(0.0, 0.0, 0.0);
  // normalcone_.setFromArray(dummy.begin(), dummy.end(), 3);
  normalcone_ = DirectionCone(dummy);
  avedglen_ = -1.0;
  region_ = 0;
  visited_ = 0;
  moved_ = 0;
  outlier_ = false;
  mark_ix_ = -1;
  sfdist_ = -1.0;
  sfang_ = -1.0;
  nmb_move_ = 0;
  meancurv0_ = meancurv_ = gausscurv0_ = gausscurv_ = 0.0;
  curvedness_ = 0.0;
  edge_[0] = PCA_EDGE_UNDEF;
  edge_[1] = C1_EDGE_UNDEF;
  edge_[2] = C2_EDGE_UNDEF;
  surf_[0] = PCA_UNDEF;
  surf_[1] = C1_UNDEF;
}

//===========================================================================
RevEngPoint::RevEngPoint(Vector3D xyz, int bnd)
  : ftSamplePoint(xyz, bnd)
//===========================================================================
{
  int dim = 3;
  eigen1_ = Point(dim);
  eigen2_ = Point(dim);
  eigen3_ = Point(dim);
  LocFuncnormal_ = Point(dim);
  kvecmin_ = Point(dim);
  kvecmax_ = Point(dim);
  lambda1_ = lambda2_ = lambda3_ = -1.0;
  kmin_ = kmax_ = 0.0;
  ptdist_ = avdist_ = 0.0;
  nmb_eigen_ = nmb_locfunc_ = 0;
  Point dummy(0.0, 0.0, 0.0);
  // normalcone_.setFromArray(dummy.begin(), dummy.end(), 3);
  normalcone_ = DirectionCone(dummy);
  avedglen_ = -1.0;
  region_ = 0;
  visited_ = 0;
  outlier_ = false;
  mark_ix_ = -1;
  sfdist_ = -1.0;
  sfang_ = -1.0;
  nmb_move_ = 0;
  meancurv0_ = meancurv_ = gausscurv0_ = gausscurv_ = 0.0;
  curvedness_ = 0.0;
  edge_[0] = PCA_EDGE_UNDEF;
  edge_[1] = C1_EDGE_UNDEF;
  edge_[2] = C2_EDGE_UNDEF;
  surf_[0] = PCA_UNDEF;
  surf_[1] = C1_UNDEF;
}

//===========================================================================
RevEngPoint::~RevEngPoint()
//===========================================================================
{
}

//===========================================================================
double RevEngPoint::getMeanEdgLen()
//===========================================================================
{
  if (avedglen_ < 0.0)
    {
      double len = 0.0;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  double currlen = xyz_.dist(next_[ki]->getPoint());
	  len += currlen;
	}
      if (next_.size() > 0)
	len /= (double)next_.size();
      avedglen_ = len;
    }
  return avedglen_;

}

//===========================================================================
double RevEngPoint::getMeanEdgLen(double maxlen)
//===========================================================================
{
  int nmb = 0;
  double len = 0.0;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      double currlen = xyz_.dist(next_[ki]->getPoint());
      if (currlen > maxlen)
	continue;
	  len += currlen;
      ++nmb;
    }
  if (nmb > 0)
    len /= (double)nmb;
  return len;

}

//===========================================================================
void RevEngPoint::computeTriangNormal(double lim)
//===========================================================================
{
  if (next_.size() == 0)
    return;
  double eps = 1.0e-10;
  size_t prev = next_.size()-1;
  Vector3D vec1 = next_[prev]->getPoint() - xyz_;

  size_t ki, kj;
  for (ki=0, kj=0; ki<next_.size(); prev=ki, ++ki)
    {
      Vector3D vec2 = next_[ki]->getPoint() - xyz_;
      Vector3D norm = vec1 % vec2;
      bool neighbour = next_[ki]->isNeighbour(next_[prev]);
  
      if (neighbour && vec1.length() <= lim && vec2.length() <= lim &&
	  norm.length() > eps)
	{
	  if (kj == 0)
	    normalcone_.setFromArray(norm.begin(), norm.end(), 3);
	  else
	    {
	      Point norm2(norm[0], norm[1], norm[2]);
	      normalcone_.addUnionWith(norm2);
	    }
	  ++kj;
	}
      else
	{
	  int stop_break = 1;
	}
      vec1 = vec2;
    }
  if (kj == 0)
    setOutlier();
}

//===========================================================================
int RevEngPoint::surfaceClassification(int classification_type) const
//===========================================================================
{
  if (classification_type == CLASSIFICATION_CURVATURE)
    return surf_[1];
  else if (classification_type == CLASSIFICATION_SHAPEINDEX)
    return surf_[2];
  else if (classification_type == CLASSIFICATION_POINTASSOCIATION)
    return surf_[3];
  else
    return CLASSIFICATION_UNDEF;
 }

//===========================================================================
Point RevEngPoint::fetchClosePoints(double radius, int min_nmb, int max_nmb,
				    vector<Point>& nearpts)
//===========================================================================
{
  int nmb_iter = 0;
  int max_iter = 5;
  // Debug
  DirectionCone pcone = normalcone_;
  while ((int)nearpts.size() < min_nmb)
    {
      setVisited();
      vector<RevEngPoint*> near;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
	  if (curr->visited())
	    continue;
	  if (xyz_.dist(curr->getPoint()) <= radius)
	    {
	      curr->setVisited();
	      near.push_back(curr);
	      pcone.addUnionWith(curr->normalcone_);
	      curr->getNearby(xyz_, radius, max_nmb, near);
	    }
	}
      
      unsetVisited();
      for (size_t ki=0; ki<near.size(); ++ki)
	{
	  Vector3D vx = near[ki]->getPoint();
	  nearpts.push_back(Point(vx[0], vx[1], vx[2]));
	  near[ki]->unsetVisited();
	  pcone.addUnionWith(near[ki]->normalcone_);
	}

      if (nmb_iter > max_iter)
	break;
      
      if (nearpts.size() < min_nmb)
	{
	  radius *= std::max(1.1, (double)min_nmb/(double)nearpts.size());
	  nearpts.clear();
	}
      ++nmb_iter;
    }

  return Point(xyz_[0], xyz_[1], xyz_[2]);
}

//===========================================================================
void RevEngPoint::fetchClosePoints2(double radius, int min_nmb, int max_nmb,
				    vector<RevEngPoint*>& nearpts,
				    RevEngRegion *region)
//===========================================================================
{
  int nmb_iter = 0;
  int max_iter = 5;
  size_t prev_nmb = nearpts.size();
  int nmb_same = 0;
  while ((int)nearpts.size() < min_nmb)
    {
      setVisited();
      vector<RevEngPoint*> near;
      for (size_t ki=0; ki<next_.size(); ++ki)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
	  if (curr->visited())
	    continue;
	  if (region && curr->region() != region)
	    continue;
	  if (xyz_.dist(curr->getPoint()) <= radius)
	    {
	      curr->setVisited();
	      near.push_back(curr);
	      curr->getNearby(xyz_, radius, max_nmb, near, region);
	    }
	}

      unsetVisited();
      for (size_t ki=0; ki<near.size(); ++ki)
	{
	  nearpts.push_back(near[ki]);
	  near[ki]->unsetVisited();
	}

      if (nmb_iter > max_iter)
	break;

      if (nearpts.size() == prev_nmb)
	++nmb_same;
      prev_nmb = nearpts.size();
      if (nearpts.size() < min_nmb && nmb_iter < max_iter)
	{
	  radius *= std::max(1.1, (double)min_nmb/(double)nearpts.size());
	  nearpts.clear();
	}
      ++nmb_iter;
    }
  // // Debug
  // DirectionCone pcone = normalcone_;
  // for (size_t ki=0; ki<nearpts.size(); ++ki)
  //   pcone.addUnionWith(nearpts[ki]->normalcone_);
    if (nearpts.size() < min_nmb && nmb_same > 0)
      {
	for (size_t ki=0; ki<nearpts.size(); ++ki)
	  nearpts[ki]->setOutlier();
	nearpts.clear();
    }
  int stop_break = 1;
}

//===========================================================================
void RevEngPoint::fetchConnected(RevEngRegion *region, int max_nmb,
				 vector<RevEngPoint*>& group)
//===========================================================================
{
  setVisited();
  group.push_back(this);
  for (size_t ki=0; ki<group.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (curr->visited())
	    continue;
	  if (curr->region() != region)
	    continue;
	  curr->setVisited();
	  group.push_back(curr);
	}
      if ((int)group.size() >= max_nmb)
	break;
    }
//   double radius = std::numeric_limits<double>::max();

//   setVisited();
//   vector<RevEngPoint*> connected;
//   for (size_t ki=0; ki<next_.size(); ++ki)
//     {
//       RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
//       if (curr->visited())
// 	continue;
//       if (curr->region() != region)
// 	continue;
//       curr->setVisited();
//       connected.push_back(curr);
//       curr->getNearby(xyz_, radius, max_nmb, connected, region);
//     }
//   group.push_back(this);
//   group.insert(group.end(), connected.begin(), connected.end());
  int stop_break = 1;
}

//===========================================================================
void RevEngPoint::fetchConnectedMarked(int mark,
				       vector<RevEngPoint*>& group)
//===========================================================================
{
  setVisited();
  group.push_back(this);
  for (size_t ki=0; ki<group.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (curr->visited())
	    continue;
	  if (curr->getMarkIx() != mark)
	    continue;
	  curr->setVisited();
	  group.push_back(curr);
	}
    }
}
  
///===========================================================================
void RevEngPoint::getNearby(Vector3D xyz, double radius, int max_nmb,
			    vector<RevEngPoint*>& near, RevEngRegion* region)
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (curr->visited())
	continue;
      if (region && curr->region() != region)
	continue;
      if (xyz.dist(curr->getPoint()) <= radius)
	{
	  curr->setVisited();
	  near.push_back(curr);
	  if (near.size() < max_nmb)
	    curr->getNearby(xyz, radius, max_nmb, near, region);
	}
    }
}
 
//===========================================================================
void
RevEngPoint::addCovarianceEigen(Point& eigen1, double lambda1, Point& eigen2,
				double lambda2, Point& eigen3, double lambda3)
//===========================================================================
{
  if (nmb_eigen_ == 0)
    {
      eigen1_ = eigen1;
      eigen2_ = eigen2;
      eigen3_ = eigen3;
      lambda1_ = lambda1;
      lambda2_ = lambda2;
      lambda3_ = lambda3;
    }
  else
    {
      eigen1_ = nmb_eigen_*eigen1_ + eigen1;
      eigen2_ = nmb_eigen_*eigen2_ + eigen2;
      eigen3_ = nmb_eigen_*eigen3_ + eigen3;
      lambda1_ = nmb_eigen_*lambda1_ + lambda1;
      lambda2_ = nmb_eigen_*lambda2_ + lambda2;
      lambda3_ = nmb_eigen_*lambda3_ + lambda3;
      eigen1_ /= (double)(nmb_eigen_+1);
      eigen2_ /= (double)(nmb_eigen_+1);
      eigen3_ /= (double)(nmb_eigen_+1);
      lambda1_ /= (double)(nmb_eigen_+1);
      lambda2_ /= (double)(nmb_eigen_+1);
      lambda3_ /= (double)(nmb_eigen_+1);
    }
  nmb_eigen_++;
}

//===========================================================================
void RevEngPoint::addLocFuncInfo(Point& norm, Point& mincvec, double minc, Point& maxcvec,
			       double maxc, double currdist, double avdist)
//===========================================================================
{
  if (nmb_locfunc_ == 0)
    {
      LocFuncnormal_ = norm;
      kvecmin_ = mincvec;
      kvecmax_ = maxcvec;
      kmin_ = minc;
      kmax_ = maxc;
      ptdist_ = currdist;
      avdist_ = avdist;
    }
  else
    {
      LocFuncnormal_ = nmb_locfunc_*LocFuncnormal_ + norm;
      kvecmin_ = nmb_locfunc_*kvecmin_ + mincvec;
      kvecmax_ = nmb_locfunc_*kvecmax_ + maxcvec;
      kmin_ = nmb_locfunc_*kmin_ + minc;
      kmax_ = nmb_locfunc_*kmax_ + maxc;
      ptdist_ = nmb_locfunc_*ptdist_ + currdist;
      avdist_ = nmb_locfunc_*avdist_ + avdist;
      LocFuncnormal_ /= (double)(nmb_locfunc_+1);
      kvecmin_ /= (double)(nmb_locfunc_+1);
      kvecmax_ /= (double)(nmb_locfunc_+1);
      kmin_ /= (double)(nmb_locfunc_+1);
      kmax_ /= (double)(nmb_locfunc_+1);
      ptdist_ /= (double)(nmb_locfunc_+1);
      avdist_ /= (double)(nmb_locfunc_+1);
    }
  nmb_locfunc_++;

  meancurv0_ = meancurv_ = 0.5*(kmin_ + kmax_);
  gausscurv0_ = gausscurv_ = kmin_*kmax_;
  curvedness_ = sqrt(0.5*(kmin_*kmin_ + kmax_*kmax_));

}


//===========================================================================
bool RevEngPoint::isolatedEdge(int edge_class_type, int nmb, bool close)
//===========================================================================
{
  if (notEdge(edge_class_type))
    return false;

  int nn = 0;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint* curr = dynamic_cast<RevEngPoint*>(next_[ki]);
      bool found = (close) ? curr->closeEdge(edge_class_type) :
	curr->isEdge(edge_class_type);
      if (found)
	++nn;
    }
  return (nn <= nmb);
}


//===========================================================================
void RevEngPoint::adjustWithTriangNorm(double anglim)
//===========================================================================
{
  double ang = getTriangAngle();
  if (ang < anglim)
    {
      if (edge_[0] == PCA_EDGE)
	edge_[0] = PCA_CLOSE_EDGE;
     if (edge_[1] == C1_EDGE)
	edge_[1] = C1_CLOSE_EDGE;
     if (edge_[2] == C2_EDGE)
	edge_[2] = C2_CLOSE_EDGE;
    }
}

//===========================================================================
void RevEngPoint::adjacentRegions(vector<RevEngRegion*>& adj) const
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (pt->region_ && pt->region_ != region_)
	{
	  size_t kj;
	  for (kj=0; kj<adj.size(); ++kj)
	    if (pt->region_ == adj[kj])
	      break;
	  if (kj == adj.size())
	    adj.push_back(pt->region_);
	}
    }
}

//===========================================================================
bool RevEngPoint::nextToRegion(RevEngRegion *reg)
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (pt->region_ && pt->region_ == reg)
	return true;
    }
  return false;
}


//===========================================================================
int RevEngPoint::nmbSameClassification(int classification_type) const
//===========================================================================
{
  int same = 0;
  int type = surfaceClassification(classification_type);
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (pt->surfaceClassification(classification_type) == type)
	same++;
    }
  return same;
 }

//===========================================================================
bool RevEngPoint::isNeighbour(RevEngRegion* reg) const
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      if (pt->region() == reg)
	return true;
    }
  return false;
}

//===========================================================================
bool RevEngPoint::isNeighbour(RevEngPoint* pt) const
//===========================================================================
{
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      if (pt == next_[ki])
	return true;
    }
  return false;
}



//===========================================================================
void RevEngPoint::getAdjInfo(double mean_edge_len, vector<RevEngRegion*>& adj_reg,
			     vector<RevEngPoint*>& adj_pt)
//===========================================================================
{
  double fac = 10.0;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      if (pntDist(next_[ki]) > fac*mean_edge_len)
	  continue;
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      RevEngRegion *curr = pt->region();
      if (curr)
	{
	  size_t kj;
	  for (kj=0; kj<adj_reg.size(); ++kj)
	    {
	      if (adj_reg[kj] == curr)
		{
		  if (pntDist(pt) < pntDist(adj_pt[kj]))
		    adj_pt[kj] = pt;
		  break;
		}
	    }
	  if (kj == adj_reg.size())
	    {
	      adj_reg.push_back(curr);
	      adj_pt.push_back(pt);
	    }
	}
    }
}

//===========================================================================
vector<RevEngRegion*> RevEngPoint::getAdjacentRegions() const
//===========================================================================
{
  vector<RevEngRegion*> adj_reg;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      RevEngRegion *curr = pt->region();
      size_t kj;
      for (kj=0; kj<adj_reg.size(); ++kj)
	if (adj_reg[kj] == curr)
	  break;
      if (kj == adj_reg.size())
	adj_reg.push_back(curr);
    }
		    
  return adj_reg;
}

//===========================================================================
vector<RevEngRegion*> RevEngPoint::adjacentRegsWithSurf() const
//===========================================================================
{
  vector<RevEngRegion*> adj_reg;
  for (size_t ki=0; ki<next_.size(); ++ki)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next_[ki]);
      RevEngRegion *curr = pt->region();
      if (!curr)
	continue;
      if (!curr->hasSurface())
	continue;
      size_t kj;
      for (kj=0; kj<adj_reg.size(); ++kj)
	if (adj_reg[kj] == curr)
	  break;
      if (kj == adj_reg.size())
	adj_reg.push_back(curr);
    }
		    
  return adj_reg;
}

//===========================================================================
int RevEngPoint::numAdjacentRegions() const
//===========================================================================
{
  vector<RevEngRegion*> adj_reg = getAdjacentRegions();
  return (int)adj_reg.size();
}

//===========================================================================
bool RevEngPoint::mergeWithAdjacent(double mean_edge_len)
//===========================================================================
{
  vector<RevEngRegion*> adj_reg;
  vector<RevEngPoint*> adj_pt;
  getAdjInfo(mean_edge_len, adj_reg, adj_pt);

  if (adj_reg.size() == 0)
    return false;

  double lentol = 2.0*mean_edge_len;
  double angtol = 0.1*M_PI;
  int minlen = std::numeric_limits<double>::max();
  int minix = -1;
  for (size_t ki=0; ki<adj_pt.size(); ++ki)
    {
      double len = pntDist(adj_pt[ki]);
      if (len > lentol)
	continue;
      Point monge = adj_pt[ki]->getLocFuncNormal();
      if (LocFuncnormal_*monge < 0.0 || LocFuncnormal_.angle(monge) > angtol)
	continue;
      int ka;
      for (ka=1; ka<4; ++ka)
	if (surf_[ka] == adj_pt[ki]->surf_[ka])
	  break;
      if (ka == 4)
	continue;
      if (len < minlen)
	{
	  minlen = len;
	  minix = (int)ki;
	}
    }

  if (minix >= 0)
    {
      adj_reg[minix]->addPoint(this);
      return true;
    }
  
  int stop_break = 1;
  return false;
}

//===========================================================================
void RevEngPoint::store(std::ostream& os) const
//===========================================================================
{
  os << index_ << " " << xyz_ << " " << uv_ << std::endl;
  os << next_.size();
  for (size_t ki=0; ki<next_.size(); ++ki)
    os << " " << next_[ki]->getIndex();
  os << std::endl;
  os << avedglen_ << " " << eigen1_ << " " << lambda1_ << " " << eigen2_;
  os << " " << lambda2_ << " " << eigen3_ << " " << lambda3_ << std::endl;
  os << LocFuncnormal_ << " " << kvecmin_ << " " << kmin_ << " " << kvecmax_;
  os << " " << kmax_ << std::endl;
  os << ptdist_ << " " << avdist_ << std::endl;
  normalcone_.write(os);
  for (int ka=0; ka<3; ++ka)
    os << " " << edge_[ka];
  for (int ka=0; ka<2; ++ka)
    os << " " << surf_[ka];
  os << " " << outlier_ << " " << sfdist_ << " " << sfang_ << std::endl;
  os << std::endl;
}

//===========================================================================
void RevEngPoint::read(std::istream& is, vector<int>& next_ix) 
//===========================================================================
{
  is >> index_ >> xyz_ >> uv_;
  int nmb_next;
  is >> nmb_next;
  if (nmb_next > 0)
    next_ix.resize(nmb_next);
  for (int ki=0; ki<nmb_next; ++ki)
    is >> next_ix[ki];
  is >> avedglen_ >> eigen1_ >> lambda1_ >> eigen2_ >> lambda2_;
  is >> eigen3_ >> lambda3_ >> LocFuncnormal_ >> kvecmin_ >> kmin_;
  is >> kvecmax_ >> kmax_ >> ptdist_ >> avdist_;
  Point dummy(3);
  //normalcone_ = DirectionCone(dummy);
  normalcone_.read(is);
    meancurv0_ = meancurv_ = 0.5*(kmin_ + kmax_);
  gausscurv0_ = gausscurv_ = kmin_*kmax_;
  curvedness_ = sqrt(0.5*(kmin_*kmin_ + kmax_*kmax_));
  for (int ka=0; ka<3; ++ka)
    is >> edge_[ka];
  for (int ka=0; ka<2; ++ka)
    is >> surf_[ka];
  is >> outlier_ >> sfdist_ >> sfang_;

}


