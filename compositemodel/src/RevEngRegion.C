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

#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/RevEngEdge.h"
#include "GoTools/compositemodel/RevEng.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/PointCloud.h"
//#include "GoTools/utils/Array.h"
//#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/Curvature.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/CurveCreators.h"
// #include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "newmat.h"
#include "newmatap.h"
#include <vector>
#include <fstream>

//#define DEBUG_JOIN
//#define DEBUG_CHECK
//#define DEBUG_UPDATE
//#define DEBUG_INTEGRATE
//#define DEBUG_TORUSCONTEXT
//#define DEBUG_CYLCONTEXT
//#define DEBUG
//#define DEBUG0
//#define DEBUG_EXTRACT
//#define DEBUG_CYL
//#define DEBUG_PLANAR
//#define DEBUG_SEGMENT
//#define DEBUG_REPAR
//#define DEBUG_ADJUST
//#define DEBUG_GROW
//#define DEBUG_MERGE
//#define DEBUG_ADJACENT
//#define DEBUG_VALIDATE
//#define DEBUG_COLLECT
//#define DEBUG_AXIS
//#define DEBUG_GROWNEIGHBOUR
//#define DEBUG_TRIM
//#define DEBUG_BLEND

using namespace Go;
using std::vector;
using std::set;
using std::pair;

//===========================================================================
RevEngRegion::RevEngRegion(int edge_class_type)
//===========================================================================
  : Id_(0), classification_type_(CLASSIFICATION_UNDEF),
    edge_class_type_(edge_class_type),
    associated_sf_(0), associated_blend_(0), blend_edge_(0), surfflag_(NOT_SET),
    surf_adaption_(INITIAL), mink1_(0.0), maxk1_(0.0), 
    mink2_(0.0), maxk2_(0.0), avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), frac_norm_in2_(0.0), maxdist_(0.0), avdist_(0.0), 
    num_inside_(0), num_inside2_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    visited_(false), to_be_removed_(false)
{
  domain_[0] = domain_[1] = domain_[2] = domain_[3] = 0.0;
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type, int edge_class_type)
//===========================================================================
  : Id_(0), classification_type_(classification_type),
    edge_class_type_(edge_class_type),
    associated_sf_(0), associated_blend_(0), blend_edge_(0), surfflag_(NOT_SET),
    surf_adaption_(INITIAL), mink1_(0.0), maxk1_(0.0), 
    mink2_(0.0), maxk2_(0.0), avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), frac_norm_in2_(0.0), maxdist_(0.0), avdist_(0.0), 
    num_inside_(0), num_inside2_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    visited_(false), to_be_removed_(false)
{
  domain_[0] = domain_[1] = domain_[2] = domain_[3] = 0.0;
  bbox_ = BoundingBox(3);
}

//===========================================================================
RevEngRegion::RevEngRegion(int classification_type,
			   int edge_class_type,
			   vector<RevEngPoint*> & points)
//===========================================================================
  : Id_(0), group_points_(points), classification_type_(classification_type),
    edge_class_type_(edge_class_type), associated_sf_(0), associated_blend_(0),
    blend_edge_(0), surfflag_(NOT_SET), surf_adaption_(INITIAL),
    avH_(0.0), avK_(0.0), MAH_(0.0), MAK_(0.0),
    frac_norm_in_(0.0), frac_norm_in2_(0.0), maxdist_(0.0), avdist_(0.0), 
    num_inside_(0), num_inside2_(0), 
    prev_region_(0), maxdist_base_(0.0), avdist_base_(0.0), num_in_base_(0),
    visited_(false), to_be_removed_(false)
{
  domain_[0] = domain_[1] = domain_[2] = domain_[3] = 0.0;

  for (size_t kj=0; kj<group_points_.size(); ++kj)
    group_points_[kj]->setRegion(this);
  
  // Bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  bbox_ = BoundingBox(3);
  double fac = 1.0/(double)group_points_.size();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      double H = group_points_[kj]->meanCurvature();
      double K = group_points_[kj]->GaussCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      avH_ += fac*H;
      avK_ += fac*K;
      MAH_ += fac*fabs(H);
      MAK_ += fac*fabs(K);
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }

  if (group_points_.size() > 0)
    normalcone_ = DirectionCone(group_points_[0]->getLocFuncNormal());
  avnorm_ = Point(0.0, 0.0, 0.0);
  if (group_points_.size() > 0)
    normalcone2_ = DirectionCone(group_points_[0]->getTriangNormal());
  avnorm2_ = Point(0.0, 0.0, 0.0);
  for  (size_t kj=1; kj<group_points_.size(); ++kj)
    {
      Point norm = group_points_[kj]->getLocFuncNormal();
      normalcone_.addUnionWith(norm);
      avnorm_ += fac*norm;
      Point norm2 = group_points_[kj]->getTriangNormal();
      normalcone2_.addUnionWith(norm2);
      avnorm2_ += fac*norm2;
    }
  // (void)avnorm_.normalize_checked();
  // (void)avnorm2_.normalize_checked();
}

//===========================================================================
RevEngRegion::~RevEngRegion()
//===========================================================================
{
  if (associated_blend_ != 0)
    associated_blend_->removeBlendReg(this);
  removeFromAdjacent();
  for (size_t ki=0; ki<rev_edges_.size(); ++ki)
    {
      rev_edges_[ki]->eraseAdjacent(this);
    }
  for (size_t ki=0; ki<associated_sf_.size(); ++ki)
    associated_sf_[ki]->removeRegion(this);
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::setAccuracy(double maxdist, double avdist, int num_inside,
			       int num_inside2)
//===========================================================================
{
  maxdist_ = maxdist;
  avdist_ = avdist;
  num_inside_ = num_inside;
  num_inside2_ = num_inside2;
 
}

//===========================================================================
void RevEngRegion::joinToCurrent(double tol, double angtol, int small_lim,
				vector<RevEngRegion*>& adapted_regions)
//===========================================================================
{
  if (!hasSurface())
    return;
  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  int classtype = surf->instanceType();
  bool cyllike = (classtype == Class_Cylinder || classtype == Class_Cone);
  double small_frac = 0.1;
  double tol_fac = 2.0;
  
#ifdef DEBUG_GROW
  std::ofstream of("master_reg.g2");
  writeRegionPoints(of);
  writeSurface(of);
#endif
  
  vector<RevEngRegion*> adj_reg(adjacent_regions_.begin(), adjacent_regions_.end());

  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      if (adj_reg[ki]->hasSurface())
	  continue;  // Not the correct grow method
      if (adj_reg[ki]->hasAssociatedBlend())
	continue;
      

#ifdef DEBUG_GROW
      std::ofstream ofn("added_region_cand.g2");
      adj_reg[ki]->writeRegionPoints(ofn);
#endif
      
     vector<RevEngPoint*> points = adj_reg[ki]->getPoints();
     int num_adj = (int)points.size();
     int num_curr = (int)group_points_.size();

      // Check if the points can be approximated by the current surface
       double maxd, avd;
       int num_in, num2_in;
       vector<RevEngPoint*> in, out;
       vector<pair<double, double> > dist_ang;
       vector<double> parvals;
       RevEngUtils::distToSurf(points.begin(), points.end(), surf, tol,
			       maxd, avd, num_in, num2_in,
			       in, out, parvals, dist_ang, angtol);
       
       int num_ang_in = 0;
       double avang = 0.0;
       double nfac = 1.0/(double)num_adj;
       for (size_t kh=0; kh<dist_ang.size(); ++kh)
	 {
	   avang += nfac*dist_ang[kh].second;
	   if (dist_ang[kh].second <= angtol)
	     num_ang_in++;
	 }
       
       int adj_sf_flag = adj_reg[ki]->defineSfFlag(0, tol, num_in, num2_in, 
						   avd, cyllike);

       // Combined numbers
       double maxd2 = std::max(maxdist_, maxd);
       double avd2 = ((double)num_curr*avdist_ + (double)num_adj*avd)/
	 (double)(num_curr+num_adj);
       int num_in2 = num_inside_ + num_in;
       int num2_in2 = num_inside2_ + num2_in;
       int sf_flag2 = defineSfFlag(num_curr+num_adj, 0, tol, num2_in2, num_in2,
				   avd2, cyllike);
       if (sf_flag2 <= surfflag_ &&
	   (adj_sf_flag < ACCURACY_POOR ||
	    (avd <= tol_fac*tol && num_adj < small_lim &&
	     (double)num_adj < small_frac*(double)num_curr)))
	 {
	   vector<RevEngRegion*> added_adjacent;
	   includeAdjacentRegion(adj_reg[ki], maxd, avd, num_in, 
				 num2_in, parvals, dist_ang,
				 added_adjacent);
	   adapted_regions.push_back(adj_reg[ki]);
	   setSurfaceFlag(sf_flag2);
	 }
    }
#ifdef DEBUG_GROW
  std::ofstream of2("updated_region.g2");
  writeRegionPoints(of2);
#endif
  int stop_break = 1;
}
 


//===========================================================================
void RevEngRegion::joinRegions(Point mainaxis[3], double approx_tol, double anglim,
				vector<RevEngRegion*>& adapted_regions)
//===========================================================================
{
  if (group_points_.size() < 5)
    return;

#ifdef DEBUG_JOIN
  std::ofstream of("region_grow.g2");
  writeRegionInfo(of);
#endif

  //double eps = 1.0e-6;
  shared_ptr<ParamSurface> surf1;
  double avdist1, maxdist1;
  int num_in1, num2_in1;
  if (basesf_.get())
    {
      surf1 = basesf_;
      getBaseDist(maxdist1, avdist1, num_in1, num2_in1);
    }
  else
    {
      surf1 = surfApprox(group_points_, bbox_);
      vector<RevEngPoint*> in1, out1;
      approximationAccuracy(group_points_, surf1, approx_tol, anglim, maxdist1, avdist1,
			    in1, out1);
      num_in1 = (int)in1.size();
    }
#ifdef DEBUG_JOIN
  std::ofstream of2("approx_sf.g2");
  surf1->writeStandardHeader(of2);
  surf1->write(of2);
#endif
  
  if (num_in1 < (int)group_points_.size()/2 || avdist1 > approx_tol)
    return;

  int kv=0;
  int numfac = 10;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end();)
    {
      ++kv;
      if ((*it)->hasSurface())
	{
	  ++it;
	  continue;  // Not the correct grow method
	}
      
      if ((*it)->visited())
	{
	  ++it;
	  continue;
	}
      (*it)->setVisited(true);

      int num = (int)group_points_.size();
      int num2 = (*it)->numPoints();
      if (num2 > numfac*num)
	{
	  ++it;
	  continue;
	}

      auto itnext = it;
      ++itnext;

#ifdef DEBUG_JOIN
      std::ofstream ofn("region_grow_cand.g2");
      (*it)->writeRegionInfo(ofn);
#endif
      
     vector<RevEngPoint*> points = (*it)->getPoints();

      // Check if the points can be approximated by the current surface
      double avdist2_0, maxdist2_0;
      vector<RevEngPoint*> in2_0, out2_0;
      approximationAccuracy(points, surf1, approx_tol, anglim, maxdist2_0, avdist2_0,
			    in2_0, out2_0);
      int num_in2_0 = (int)in2_0.size();

      // Compute overall numbers
      int tot_num = (int)group_points_.size()+(int)points.size();
      double frac = (double)(num_in1+num_in2_0)/(double)(tot_num);
      double avd = (avdist1*(double)group_points_.size() +
		    avdist2_0*(double)points.size())/(double)tot_num;
      //if (in2_0.size() > points.size()/2 && avdist2_0 < approx_tol)
      if (frac > 0.5 && avd < approx_tol)
	{
	  // Include region in current
	  adapted_regions.push_back((*it));
	  for (auto it2=(*it)->adjacent_regions_.begin();
	       it2!=(*it)->adjacent_regions_.end(); ++it2)
	    {
	      if (*it2 != this)
		{
		  addAdjacentRegion(*it2);
		  (*it2)->addAdjacentRegion(this);
		  (*it2)->removeAdjacentRegion(*it);
		}
	    }
	  maxdist1 = std::max(maxdist1, maxdist2_0);
	  avdist1 = (num*avdist1 + num2*avdist2_0)/(double)(num+num2);
	  for (auto it3=(*it)->pointsBegin(); it3!=(*it)->pointsEnd(); ++it3)
	    (*it3)->addMove();
	  vector<pair<double,double> > dummy;
	  addRegion((*it), dummy);
	  removeAdjacentRegion(*it);
	  it = itnext;
	  continue;
	}
			     

      bool always_stop = true;
      if (num2 < 5 || num2 > num || always_stop)
	{
	  it = itnext;
	  continue;
	}
      
      shared_ptr<ParamSurface> surf2;

      if (basesf_.get() &&
	  ((!(*it)->basesf_.get()) ||
	   ((*it)->basesf_.get() &&
	    basesf_->instanceType() == (*it)->basesf_->instanceType())))
	{
	  vector<RevEngPoint*> all_pts;
	  all_pts.insert(all_pts.end(), group_points_.begin(), group_points_.end());
	  all_pts.insert(all_pts.end(), points.begin(), points.end());
	  if (basesf_->instanceType() == Class_Plane)
	    surf2 = computePlane(all_pts, avnorm_, mainaxis);
	  else if (basesf_->instanceType() == Class_Cylinder)
	    {
	      surf2 = computeCylinder(all_pts, approx_tol);
	    }
	}
      else
	surf2 = surfApprox(points, (*it)->getBbox());

      if (!surf2.get())
	{
	  it = itnext;
	  continue;
	}
	
#ifdef DEBUG_JOIN
      std::ofstream of3("approx_sf2.g2");
      surf2->writeStandardHeader(of3);
      surf2->write(of3);
#endif
      
      double avdist2, maxdist2;
      vector<RevEngPoint*> in2, out2;
      approximationAccuracy(points, surf2, approx_tol, anglim,
			    maxdist2, avdist2, in2, out2);
      if (in2.size() < points.size()/2 && maxdist2 > approx_tol)
	{
	  it = itnext;
	  continue;
	}
      int num_in2 = (int)in2.size();

      vector<RevEngPoint*> mergept(group_points_.begin(), group_points_.end());
      mergept.insert(mergept.end(), (*it)->pointsBegin(),
		     (*it)->pointsEnd());
      BoundingBox bbox = bbox_;
      bbox.addUnionWith((*it)->getBbox());

      shared_ptr<SplineSurface> surf3 = surfApprox(mergept, bbox);
#ifdef DEBUG_JOIN
      std::ofstream of4("approx_sf3.g2");
      surf3->writeStandardHeader(of4);
      surf3->write(of4);
#endif
      
      double avdist3_0, maxdist3_0;
      vector<RevEngPoint*> in3_0, out3_0;
      approximationAccuracy(points, surf3, approx_tol, anglim,
			    maxdist3_0, avdist3_0, in3_0, out3_0);
      if (in3_0.size() < points.size()/2 || avdist3_0 > approx_tol)
	{
	  it = itnext;
	  continue;
	}

      
      double avdist3, maxdist3;
      vector<RevEngPoint*> in3, out3;
      approximationAccuracy(mergept, surf3, approx_tol, anglim,
			    maxdist3, avdist3, in3, out3);
      int num_in3 = (int)in3.size();
      
      if (/*maxdist3 <= 1.1*std::max(maxdist1, maxdist2) && */
	  avdist3 <= 1.1*std::max(avdist1, avdist2) && avdist3 < approx_tol &&
	  num_in3 > num_in1 && num_in3 > 3*(num_in1+num_in2)/4)
	{
	  // Include region in current
	  adapted_regions.push_back((*it));
	  for (auto it2=(*it)->adjacent_regions_.begin();
	       it2!=(*it)->adjacent_regions_.end(); ++it2)
	    {
	      if (*it2 != this)
		{
		  addAdjacentRegion(*it2);
		  (*it2)->addAdjacentRegion(this);
		  (*it2)->removeAdjacentRegion(*it);
		}
	    }
	  maxdist1 = maxdist3;
	  avdist1 = avdist3;
	  for (auto it3=(*it)->pointsBegin(); it3!=(*it)->pointsEnd(); ++it3)
	    (*it3)->addMove();
	  vector<pair<double,double> > dummy;
	  addRegion((*it), dummy);
	  removeAdjacentRegion(*it);
	  num_in1 = num_in3;
	  surf1 = surf3;
	}

#ifdef DEBUG_JOIN
      std::ofstream ofl("region_grow2.g2");
      writeRegionInfo(ofl);
#endif
      int stop_break = 1;
      it = itnext;
    }

  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    (*it)->setVisited(false);
      
}

  
//===========================================================================
void RevEngRegion::approximationAccuracy(vector<RevEngPoint*>& points,
					 shared_ptr<ParamSurface> surf,
					 double tol, double angtol,
					 double& maxd, double& avd,
					 vector<RevEngPoint*>& in,
					 vector<RevEngPoint*>& out)
//===========================================================================
{
  avd = maxd = 0.0;
  double eps = 1.0e-6;
  double wgt = 1.0/(double)points.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pos, upar, vpar, close, dist, eps);
      avd += wgt*dist;
			       maxd = std::max(maxd, dist);
      if (dist < tol)
	{
	  Point normal;
	  surf->normal(normal, upar, vpar);
	  double ang = normal.angle(points[ki]->getLocFuncNormal());
	  ang = std::min(ang, M_PI-ang);
	  if (ang < angtol)
	    in.push_back(points[ki]);
	  else
	    out.push_back(points[ki]);
	}
      else
	out.push_back(points[ki]);
    }
}

//===========================================================================
void RevEngRegion::splitPlanar(double lim_cone, int min_point_reg,
			       vector<vector<RevEngPoint*> >& other_groups,
			       vector<RevEngPoint*>& single)
//===========================================================================
{
  vector<vector<RevEngPoint*> > groups;
  vector<DirectionCone> cones;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point norm = group_points_[ki]->getTriangNormal();
      size_t kj;
      for (kj=0; kj<cones.size(); ++kj)
	{
	  double ang = cones[kj].unionAngle(norm);
	  if (ang <= lim_cone)
	    {
	      cones[kj].addUnionWith(norm);
	      groups[kj].push_back(group_points_[ki]);
	      break;
	    }
	}
      if (kj == cones.size())
	{
	  DirectionCone curr(norm);
	  cones.push_back(curr);
	  vector<RevEngPoint*> other;
	  other.push_back(group_points_[ki]);
	  groups.push_back(other);
	}
    }
#ifdef DEBUG_PLANAR
  std::ofstream of("planar_groups.g2");
  for (size_t ki=0; ki<groups.size(); ++ki)
    {
      of << "400 1 0 4 100 155 0 255" << std::endl;
      of << groups[ki].size() << std::endl;
      for (size_t kj=0; kj<groups[ki].size(); ++kj)
	of << groups[ki][kj]->getPoint() << std::endl;
    }
#endif
  if (groups.size() == 1)
    return;
  
  vector<RevEngPoint*> remain;
  for (int ka=(int)groups.size()-1; ka>=0; --ka)
    {
      if ((int)groups[ka].size() < min_point_reg)
	{
	  remain.insert(remain.end(), groups[ka].begin(), groups[ka].end());
	  groups.erase(groups.begin()+ka);
	}
    }

  size_t num_groups = groups.size();
  for (size_t ki=0; ki<num_groups; ++ki)
    {
      vector<vector<RevEngPoint*> > connected;
      vector<RevEngPoint*> dummy;
      connectedGroups(groups[ki], connected, false, dummy);
      for (size_t kj=1; kj<connected.size(); ++kj)
	if (connected[kj].size() > connected[0].size())
	  std::swap(connected[kj], connected[0]);

      groups[ki] = connected[0];
      for (size_t kj=1; kj<connected.size(); ++kj)
	{
	  if (connected[kj].size() >= min_point_reg)
	    groups.push_back(connected[kj]);
	  else
	    remain.insert(remain.end(), connected[kj].begin(), connected[kj].end());
	}
    }

  for (size_t ki=1; ki<groups.size(); ++ki)
    if (groups[ki].size() > groups[0].size())
      std::swap(groups[ki], groups[0]);

  std::swap(group_points_, groups[0]);
  updateInfo();

  if (groups.size() > 0)
    other_groups.insert(other_groups.end(), groups.begin()+1, groups.end());

  if (remain.size() > 0)
    {
      vector<vector<RevEngPoint*> > connected2;
      vector<RevEngPoint*> dummy2;
      connectedGroups(remain, connected2, false, dummy2);
      for (size_t kj=0; kj<connected2.size(); ++kj)
	{
	  if (connected2[kj].size() == 1)
	    {
	      connected2[kj][0]->unsetRegion();
	      single.push_back(connected2[kj][0]);
	    }
	  else
	    other_groups.push_back(connected2[kj]);
	}
    }
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::updateRegion(double approx_tol, double anglim,
				vector<RevEngRegion*>& adapted_regions,
				vector<shared_ptr<RevEngRegion> >& outdiv_regions)
//===========================================================================
{
  if (group_points_.size() < 5)
    return;

#ifdef DEBUG_UPDATE
  std::ofstream of("region_grow.g2");
  writeRegionInfo(of);
#endif
  double eps = 1.0e-6;
  shared_ptr<SplineSurface> surf = surfApprox(group_points_, bbox_);
#ifdef DEBUG_UPDATE
  std::ofstream of2("approx_sf.g2");
  surf->writeStandardHeader(of2);
  surf->write(of2);
#endif  
  // Check if the accuracy is sufficient as a basis for growth
  double avdist = 0.0, maxdist = 0.0;
  vector<RevEngPoint*> in, out;
  approximationAccuracy(group_points_, surf, approx_tol, anglim, maxdist, avdist,
			in, out);
  if (in.size() < group_points_.size()/2 || avdist > approx_tol)
    return;

  setBaseSf(surf, maxdist, avdist, (int)in.size(), (int)in.size());
  
  // Grow
  vector<RevEngPoint*> visited;
  int numpt = (int)group_points_.size();
  int numfac = 10;
  // std::set<RevEngPoint*> tmpset0(in.begin(), in.end());
  // if (tmpset0.size() != in.size())
  //   std::cout << "Point number mismatch, init. " << tmpset0.size() << " " << in.size() << std::endl;
  size_t ki = 0;
  for (int ka=0; ka<10; ++ka)
    {
      size_t in_size = in.size();
      for (; ki<in.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = in[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (curr->moved())
		continue;  // Already in an extended region
	      if (curr->isEdge(edge_class_type_))
		continue;  // Do not include points labelled as edge
	      if (curr->region() == this)
		continue;
	      if (curr->visited())
		continue;

	      // if (std::find(in.begin(), in.end(), curr) != in.end())
	      // 	{
	      // 	  std::cout << "Double point" << curr << std::endl;
	      // 	}
	      curr->setVisited();
	      visited.push_back(curr);
	      int numpt2 = 0;
	      if (curr->region())
		numpt2 = curr->region()->numPoints();
	      if (numpt2 > numfac*numpt)
		continue;
	      Vector3D xyz = curr->getPoint();
	      Point pos(xyz[0], xyz[1], xyz[2]);
	      double upar, vpar, dist;
	      Point close;
	      surf->closestPoint(pos, upar, vpar, close, dist, eps);
	      if (dist < approx_tol)
		{
		  Point normal;
		  surf->normal(normal, upar, vpar);
		  double ang = normal.angle(curr->getLocFuncNormal());
		  ang = std::min(ang, M_PI-ang);
		  if (ang < anglim)
		    in.push_back(curr);
		  else
		    out.push_back(curr);
		}
	      else
		out.push_back(curr);
	    }
	}
#ifdef DEBUG_UPDATE
      std::ofstream of3("in_out.g2");
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << in.size() << std::endl;
      for (size_t kj=0; kj<in.size(); ++kj)
	of3 << in[kj]->getPoint() << std::endl;
      of3 << "400 1 0 4 0 255 0 255" << std::endl;
      of3 << out.size() << std::endl;
      for (size_t kj=0; kj<out.size(); ++kj)
	of3 << out[kj]->getPoint() << std::endl;
#endif
      if (in.size() == in_size)
	break;

      vector<RevEngPoint*> in_out(in.begin(), in.end());
      in_out.insert(in_out.end(), out.begin(), out.end());
      shared_ptr<SplineSurface> surf2 = surfApprox(in_out, bbox_);
#ifdef DEBUG_UPDATE
      std::ofstream of4("approx_sf2.g2");
      surf2->writeStandardHeader(of4);
      surf2->write(of4);
#endif
      //int nmb_in2 = 0;
      double avdist2 = 0.0;
      double maxdist2 = 0.0;
      in.clear();
      out.clear();
      double wgt = 1.0/(double)in_out.size();
      for (size_t ki=0; ki<in_out.size(); ++ki)
	{
	  Vector3D xyz = in_out[ki]->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  double upar, vpar, dist;
	  Point close;
	  surf2->closestPoint(pos, upar, vpar, close, dist, eps);
	  avdist2 += wgt*dist;
	  maxdist2 = std::max(maxdist2,dist);
	  if (dist < approx_tol)
	    {
	      Point normal;
	      surf2->normal(normal, upar, vpar);
	      double ang = normal.angle(in_out[ki]->getLocFuncNormal());
	      ang = std::min(ang, M_PI-ang);
	      if (ang < anglim)
		in.push_back(in_out[ki]);
	      else
		out.push_back(in_out[ki]);
	    }
	  else
	    out.push_back(in_out[ki]);
	}

      if (in.size() < out.size() || avdist2 > approx_tol)
	{
	  break;
	}

      setBaseSf(surf2, maxdist2, avdist2, (int)in.size(), (int)in.size());
      surf = surf2;
      // std::set<RevEngPoint*> tmpset(in.begin(), in.end());
      // if (tmpset.size() != in.size())
      // 	std::cout << "Point number mismatch. " << ka << " " << tmpset.size() << " " << in.size() << std::endl;
      int stop_break0 = 1;
    }

  //std::cout << "Ready to move points " << std::endl;
  // Move in points
  set<RevEngRegion*> affected_reg;
  for (size_t ki=0; ki<in.size(); ++ki)
    {
      in[ki]->setMoved();
      
      if (in[ki]->region() == this)
	continue;

      RevEngRegion *other_reg = in[ki]->region();
      if (other_reg)
	{
	  other_reg->removePoint(in[ki]);
	  affected_reg.insert(other_reg);
	}
      in[ki]->setRegion(this);
      in[ki]->addMove();
      group_points_.push_back(in[ki]);
    }

  for (size_t ki=0; ki<visited.size(); ++ki)
    visited[ki]->unsetVisited();
      
  vector<RevEngRegion*> affected_reg2(affected_reg.begin(), affected_reg.end());
  for (size_t ki=0; ki<affected_reg2.size(); ++ki)
    {
      if (affected_reg2[ki]->numPoints() == 0)
	{
	  for (auto it=affected_reg2[ki]->adjacent_regions_.begin();
	       it!=affected_reg2[ki]->adjacent_regions_.end(); ++it)
	    {
	      // if ((*it) != this)
	      //   {
	      // 	addAdjacentRegion(*it);
	      (*it)->removeAdjacentRegion(affected_reg2[ki]);
	      // }
	    }
	  // removeAdjacentRegion(affected_reg2[ki]);
	  adapted_regions.push_back(affected_reg2[ki]);
	}
      else
	{
	  vector<vector<RevEngPoint*> > separate_groups;
	  affected_reg2[ki]->splitRegion(separate_groups);
	  int classtype = affected_reg2[ki]->getClassificationType();
	  for (size_t ki=0; ki<separate_groups.size(); ++ki)
	    {
	      shared_ptr<RevEngRegion> reg(new RevEngRegion(classtype,
							    edge_class_type_,
							    separate_groups[ki]));
	      outdiv_regions.push_back(reg);
	    }
	  affected_reg2[ki]->updateInfo(approx_tol, anglim);
	}
    }

  updateInfo(approx_tol, anglim);
  int stop_break = 1;
  
}

//===========================================================================
bool RevEngRegion::segmentByPlaneAxis(Point mainaxis[3], int min_point_in,
				      int min_pt_reg,
				      double tol, double angtol, int prefer_elementary,
				      vector<RevEngRegion*>& adj_planar,
				      vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      vector<shared_ptr<RevEngRegion> >& added_reg,
				      vector<HedgeSurface*>& prevsfs,
				      vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
#ifdef DEBUG_SEGMENT
  std::ofstream of("adj_planar.g2");
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    adj_planar[ki]->writeRegionPoints(of);
#endif
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem(adj_planar.size());
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      if (!adj_planar[ki]->hasSurface())
	return false;
      shared_ptr<ParamSurface> surf = adj_planar[ki]->getSurface(0)->surface();
      shared_ptr<ElementarySurface> elemsf =
	dynamic_pointer_cast<Plane,ParamSurface>(surf);
      if (!elemsf.get())
	return false;
      adj_elem[ki] = std::make_pair(elemsf, adj_planar[ki]);
    }

  // Get candidate axis directions
  vector<Point> dir;
  size_t kr;
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      Point dir1 = adj_elem[ki].first->direction();
      for (kr=0; kr<dir.size(); ++kr)
	{
	  double ang = dir[kr].angle(dir1);
	  ang = std::min(ang, M_PI-ang);
	  if (ang <= angtol)
	    break;
	}
      if (kr == dir.size())
	dir.push_back(dir1);
      for (size_t kj=ki+1; kj<adj_elem.size(); ++kj)
	{
	  Point dir2 = adj_elem[kj].first->direction();
	  double ang = dir1.angle(dir2);
	  ang = std::min(ang, M_PI-ang);
	  if (ang > angtol)
	    {
	      Point dir3 = dir1.cross(dir2);
	      for (kr=0; kr<dir.size(); ++kr)
		{
		  ang = dir[kr].angle(dir3);
		  ang = std::min(ang, M_PI-ang);
		  if (ang <= angtol)
		    break;
		}
	      if (kr == dir.size())
		dir.push_back(dir3);
	    }
	}
    }

  vector<vector<RevEngPoint*> > group1, group2;
  vector<RevEngPoint*> remaining;
  double axisang = 0.1*M_PI;
  double planeang = 0.05;
  bool divided = sortByAxis(dir, tol, axisang, planeang, group1, group2, remaining);

  vector<RevEngPoint*> all;
  for (size_t ki=0; ki<group1.size(); ++ki)
    all.insert(all.end(), group1[ki].begin(), group1[ki].end());
  for (size_t ki=0; ki<group2.size(); ++ki)
    all.insert(all.end(), group2[ki].begin(), group2[ki].end());
  all.insert(all.end(), remaining.begin(), remaining.end());
  std::set<RevEngPoint*> tmpset1(all.begin(), all.end());
  if (tmpset1.size() != all.size())
    std::cout << "Point number mismatch,  all. " << tmpset1.size() << " " << all.size() << " " << group_points_.size() << std::endl;
  
  vector<RevEngPoint*> remain1;
  bool found1 = integratePlanarPoints(dir, group1, adj_elem, tol, angtol,
				      remain1);
  std::set<RevEngPoint*> tmpset2(remain1.begin(), remain1.end());
  if (tmpset2.size() != remain1.size())
    std::cout << "Point number mismatch,  remain1. " << tmpset2.size() << " " << remain1.size() << std::endl;
  
#ifdef DEBUG_SEGMENT
  std::ofstream ofr1("remain1.g2");
  ofr1 << "400 1 0 4 50 155 50 255" << std::endl;
  ofr1 << remain1.size() << std::endl;
  for (size_t ki=0; ki<remain1.size(); ++ki)
    ofr1 << remain1[ki]->getPoint() << std::endl;
#endif

  if (remain1.size() > 0)
    remaining.insert(remaining.end(), remain1.begin(), remain1.end());
  
  std::set<RevEngPoint*> tmpset3(remaining.begin(), remaining.end());
  if (tmpset3.size() != remaining.size())
    std::cout << "Point number mismatch,  remaining2. " << tmpset3.size() << " " << remaining.size() << std::endl;
  
#ifdef DEBUG_SEGMENT
  std::ofstream of2("adj_planar2.g2");
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    adj_planar[ki]->writeRegionPoints(of2);
#endif

  vector<RevEngPoint*>  remain2;
  bool found2 = defineCylindricalRegs(mainaxis, group2, min_point_in, 
				      min_pt_reg, tol, angtol, 
				      added_reg, hedgesfs, remain2);
  std::set<RevEngPoint*> tmpset4(remain2.begin(), remain2.end());
  if (tmpset4.size() != remain2.size())
    std::cout << "Point number mismatch,  remain2. " << tmpset4.size() << " " << remain2.size() << std::endl;
  
#ifdef DEBUG_SEGMENT
  std::ofstream ofr2("remain2.g2");
  ofr2 << "400 1 0 4 50 50 155 255" << std::endl;
  ofr2 << remain2.size() << std::endl;
  for (size_t ki=0; ki<remain2.size(); ++ki)
    ofr2 << remain2[ki]->getPoint() << std::endl;
#endif
  if (remain2.size() > 0)
    remaining.insert(remaining.end(), remain2.begin(), remain2.end());
  
  std::set<RevEngPoint*> tmpset5(remaining.begin(), remaining.end());
  if (tmpset5.size() != remaining.size())
    std::cout << "Point number mismatch,  remaining3. " << tmpset5.size() << " " << remaining.size() << std::endl;

  std::vector<RevEngPoint*> dummy;
  connectedGroups(remaining, added_groups, false, dummy);
  int max_num = 0;
  int max_ix = -1;
  for (size_t ki=0; ki<added_groups.size(); ++ki)
    if ((int)added_groups[ki].size() > max_num)
      {
	max_num = (int)added_groups[ki].size();
	max_ix = (int)ki;
      }

  if (max_ix >= 0)
    {
      std::swap(group_points_, added_groups[max_ix]);
      std::swap(added_groups[max_ix], added_groups[added_groups.size()-1]);
      added_groups.pop_back();
      updateInfo(tol, angtol);
    }
  
#ifdef DEBUG_SEGMENT
  std::ofstream of3("added_reg_surf.g2");
  for (size_t ki=0; ki<added_reg.size(); ++ki)
    {
      added_reg[ki]->writeRegionPoints(of3);
      if (added_reg[ki]->hasSurface())
	added_reg[ki]->writeSurface(of3);
    }
#endif
  return (added_reg.size() > 0);//(found1 || found2);
}

//===========================================================================
bool RevEngRegion::defineCylindricalRegs(Point mainaxis[3],
					 vector<vector<RevEngPoint*> >& groups,
					 int min_point, int min_pt_reg,
					 double tol, double angtol,
					 vector<shared_ptr<RevEngRegion> >& added_reg,
					 vector<shared_ptr<HedgeSurface> >& hedgesfs,
					 vector<RevEngPoint*>& remaining)
//===========================================================================
{
  bool found = false;
  vector<HedgeSurface*> prevsfs;   // Not feasible
  for (size_t ki=0; ki<groups.size(); ++ki)
    {
      if (groups[ki].size() == 0)
	continue;
      
      std::set<RevEngPoint*> tmpset4(groups[ki].begin(), groups[ki].end());
      if (tmpset4.size() != groups[ki].size())
	std::cout << "Point number mismatch,  axis group  " << ki << " " << tmpset4.size() << " " << groups[ki].size() << std::endl;
      // Split into connected groups
      vector<vector<RevEngPoint*> > connected;
      std::vector<RevEngPoint*> dummy;
      connectedGroups(groups[ki], connected, false, dummy);
      for (size_t kj=0; kj<connected.size(); ++kj)
	{
	  // Identify adjacent regions
	  std::set<RevEngRegion*> adj_reg;
	  for (size_t kr=0; kr<connected[kj].size(); ++kr)
	    {
	      vector<ftSamplePoint*> next = connected[kj][kr]->getNeighbours();
	      for (size_t kh=0; kh<next.size(); ++kh)
		{
		  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kh]);
		  RevEngRegion *adj = curr->region();
		  if (adj && adj != this)
		    adj_reg.insert(adj);
		}
	    }

	  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem;
	  for (auto it=adj_reg.begin(); it!=adj_reg.end(); ++it)
	    {
	      if ((*it)->hasSurface())
		{
		  shared_ptr<ParamSurface> surf = (*it)->getSurface(0)->surface();
		  shared_ptr<ElementarySurface> elemsf =
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
		  if (elemsf.get())
		    adj_elem.push_back(std::make_pair(elemsf, (*it)));
		}
	    }

	  if (adj_elem.size() > 1)
	    {
	      // Try to fit a cylinder from context
	      shared_ptr<RevEngRegion> tmp_reg(new RevEngRegion(classification_type_,
								edge_class_type_,
								connected[kj]));
#ifdef DEBUG_SEGMENT
	      std::ofstream of("cylinder_cand.g2");
	      tmp_reg->writeRegionPoints(of);
#endif
	      vector<vector<RevEngPoint*> > added_groups;
	      bool found1 = tmp_reg->contextCylinder(mainaxis, tol, min_point,
						     min_pt_reg, angtol, 0,
						     adj_elem, hedgesfs, prevsfs,
						     added_groups);

	      for (size_t kr=0; kr<added_groups.size(); ++kr)
		{
		  for (size_t kh=0; kh<added_groups[kr].size(); ++kh)
		    added_groups[kr][kh]->setRegion(this);
		  remaining.insert(remaining.end(), added_groups[kr].begin(),
				   added_groups[kr].end());
		}
	      
	      if (found1)
		{
		  found = true;
		  tmp_reg->setRegionAdjacency();
		  added_reg.push_back(tmp_reg);
		}
	      else
		{
		  vector<RevEngPoint*> curr_points(tmp_reg->pointsBegin(),
						   tmp_reg->pointsEnd());
		  for (size_t kh=0; kh<curr_points.size(); ++kh)
		    curr_points[kh]->setRegion(this);
		  remaining.insert(remaining.end(), curr_points.begin(),
		  		   curr_points.end());
		}
	    }
	  else
	    remaining.insert(remaining.end(), connected[kj].begin(),
	    		     connected[kj].end());
	}
    }
  
  return found;
}

//===========================================================================
bool RevEngRegion::integratePlanarPoints(vector<Point>& dir,
					 vector<vector<RevEngPoint*> >& groups,
					 vector<pair<shared_ptr<ElementarySurface>,RevEngRegion*> >& adj_elem,
					 double tol, double angtol,
					 vector<RevEngPoint*>& remaining)
//===========================================================================
{
  bool found = false;
  for (size_t ki=0; ki<groups.size(); ++ki)
    {
      if (groups[ki].size() == 0)
	continue;

      // Split into connected groups
      vector<vector<RevEngPoint*> > connected;
      std::vector<RevEngPoint*> dummy;
      connectedGroups(groups[ki], connected, false, dummy);
      for (size_t kj=0; kj<connected.size(); ++kj)
	{
	  // Identify compatible adjacent plane
	  std::set<RevEngRegion*> adj_reg;
	  for (size_t kr=0; kr<connected[kj].size(); ++kr)
	    {
	      vector<ftSamplePoint*> next = connected[kj][kr]->getNeighbours();
	      for (size_t kh=0; kh<next.size(); ++kh)
		{
		  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kh]);
		  RevEngRegion *adj = curr->region();
		  if (adj != this)
		    adj_reg.insert(adj);
		}
	    }

	  bool integrated = false;
	  for (auto it=adj_reg.begin(); it!=adj_reg.end(); ++it)
	    {
	      for (size_t kr=0; kr<adj_elem.size(); ++kr)
		{
		  if (adj_elem[kr].second == (*it))
		    {
		      Point curr_dir = dir[ki/2];
		      if (ki%2 == 0)
			curr_dir *= -1;
		      Point norm = adj_elem[kr].first->direction();
		      Point vec = adj_elem[kr].second->getMeanNormal();
		      if (norm*vec < 0.0)
			norm *= -1;
		      double ang = curr_dir.angle(norm);
		      if (ang < angtol)
			{
			  // Check distance
			  Point loc = adj_elem[kr].first->location();
			  //double avdist = 0.0;
			  shared_ptr<ParamSurface> surf = adj_elem[kr].first;
			  double maxd, avd;
			  int nmb_in, nmb2_in;
			  vector<RevEngPoint*> in, out;
			  vector<pair<double,double> > dist_ang;
			  vector<double> parvals;
			  RevEngUtils::distToSurf(connected[kj].begin(), connected[kj].end(),
						  surf, tol, maxd, avd, nmb_in, nmb2_in,
						  in, out, parvals, dist_ang, angtol);
			  if (avd <= tol)
			    {
			      // Integrate group in adjacent plane
			      for (size_t kh=0; kh<connected[kj].size(); ++kh)
				{
				  connected[kj][kh]->setPar(Vector2D(parvals[2*kh],
								     parvals[2*kh+1]));
				  connected[kj][kh]->setSurfaceDist(dist_ang[kh].first,
								    dist_ang[kh].second);
				  adj_elem[kr].second->addPoint(connected[kj][kh]);
				}
			      integrated = true;
			      break;
			    }
			}
		    }
		  if (integrated)
		    break;
		}
	      if (integrated)
		break;
	    }

	  if (integrated)
	    found = true;
	  else
	    remaining.insert(remaining.end(), connected[kj].begin(), connected[kj].end());
	}
    }
  return found;
}

//===========================================================================
void RevEngRegion::growPlaneOrCyl(Point mainaxis[3], int min_pt_reg,
				  double tol, double angtol,
				  vector<RevEngRegion*>& grown_regions,
				  vector<HedgeSurface*>& adj_surfs,
				  vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
  if (associated_sf_.size() == 0)
    return;  // No surface from which to grow
  
  int sfcode;
  ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  if (classtype != Class_Plane && classtype != Class_Cylinder &&
      classtype != Class_Cone)
    return;

  if (surfflag_ != ACCURACY_OK)
    return;  // Grow only surfaces with good accuracy

  vector<RevEngRegion*> adj_reg;
  adj_reg.insert(adj_reg.end(), adjacent_regions_.begin(), adjacent_regions_.end());
  int num_points = (int)group_points_.size();
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      // Get seed points
      vector<RevEngPoint*> adj_pts =
	adj_reg[ki]->extractNextToAdjacent(this);

      // Grow
      growFromNeighbour(mainaxis, min_pt_reg, adj_pts, tol,
			angtol, adj_reg[ki], false);
      if (adj_reg[ki]->numPoints() == 0)
	{
	  for (auto it=adj_reg[ki]->adjacent_regions_.begin();
	       it != adj_reg[ki]->adjacent_regions_.end(); ++it)
	    (*it)->removeAdjacentRegion(adj_reg[ki]);
	  if (adj_reg[ki]->hasSurface())
	    {
	      int num_sf = adj_reg[ki]->numSurface();
	      for (int kb=0; kb<num_sf; ++kb)
		adj_surfs.push_back(adj_reg[ki]->getSurface(kb));
	    }
	  removeAdjacentRegion(adj_reg[ki]);
	  grown_regions.push_back(adj_reg[ki]);
	}
      else
	{
	  // Make sure that the adjacent region is connected
	  vector<vector<RevEngPoint*> > separate_groups;
	  adj_reg[ki]->splitRegion(separate_groups);
	  if (separate_groups.size() > 0)
	    {
	      added_groups.insert(added_groups.end(), separate_groups.begin(),
				  separate_groups.end());
	      if (adj_reg[ki]->hasSurface())
		adj_reg[ki]->checkReplaceSurf(mainaxis, min_pt_reg, tol,
					       angtol);
	    }
	}
    }
  if ((int)group_points_.size() > num_points)
    checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
      
#ifdef DEBUG_GROWNEIGHBOUR
  std::ofstream of("grown_group.g2");
  writeRegionPoints(of);
#endif
  if (false) //(int)group_points_.size() > num_points)
    {
      // Compute distance
      shared_ptr<ParamSurface> surf = getSurface(0)->surface();
      double maxdist, avdist;
      int num_in, num2_in;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      vector<RevEngPoint*> inpt, outpt; 
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf, tol, maxdist, avdist, num_in, num2_in,
			  inpt, outpt, parvals, dist_ang, angtol);

      int sf_flag = defineSfFlag(0, tol, num_in, num2_in,
				 avdist, (classtype == Class_Cylinder));
      for (size_t kh=0; kh<group_points_.size(); ++kh)
	{
	  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	}
      setAccuracy(maxdist, avdist, num_in, num2_in);
      setSurfaceFlag(sf_flag);
      checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
    }
}

//===========================================================================
bool RevEngRegion::segmentByAdjSfContext(Point mainaxis[3], int min_point_in, 
					 int min_pt_reg, double tol, double angtol,
					 vector<RevEngRegion*>& adj_planar,
					 vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
#ifdef DEBUG_SEGMENT
  std::ofstream of("adj_planar.g2");
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    adj_planar[ki]->writeRegionPoints(of);
#endif
  //double lim = 0.1;
  int num_pt = numPoints();
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      if ((double)adj_planar[ki]->numPoints() > min_point_in/*lim*num_pt*/)
	{
	  // Get seed points
	  vector<RevEngPoint*> adj_pts = extractNextToAdjacent(adj_planar[ki]);
	  
	  // Grow adjacent
	  adj_planar[ki]->growFromNeighbour(mainaxis, min_pt_reg,
					    adj_pts, tol, angtol, this);
	}
#ifdef DEBUG_SEGMENT
      std::ofstream of2("updated_planar.g2");
      writeRegionInfo(of2);
      adj_planar[ki]->writeRegionPoints(of2);
#endif
    }
  if (numPoints() > 0 && numPoints() != num_pt)
    {
      splitRegion(added_groups);
      updateInfo(tol, angtol);
      updateRegionAdjacency();
// #ifdef DEBUG_SEGMENT
//       if (!isConnected())
// 	std::cout << "Disconnect, segmentByAdjSfContext" << std::endl;
// #endif
      return (numPoints() > num_pt/10) ? true : false;
      //return true;
    }
  else if (numPoints() == 0)
    {
      removeFromAdjacent();
      clearRegionAdjacency();
    }

  return false;
}

//===========================================================================
bool RevEngRegion::segmentByDirectionContext(int min_point_in, double tol,
					     const Point& dir, double angtol, 
					     vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  vector<vector<RevEngPoint*> > pnt_groups(3);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point norm = group_points_[ki]->getLocFuncNormal();
      double ang = dir.angle(norm);
      ang = std::min(ang, M_PI-ang);
      if (ang < angtol)
	pnt_groups[0].push_back(group_points_[ki]);
      else if (fabs(pihalf-ang) < angtol)
	pnt_groups[1].push_back(group_points_[ki]);
      else
	pnt_groups[2].push_back(group_points_[ki]);
    }

  int num_pnt_groups = (pnt_groups[0].size() > 0) + (pnt_groups[1].size() > 0) +
    (pnt_groups[2].size() > 0);
  if (num_pnt_groups <= 1)
    return false;
#ifdef DEBUG_SEGMENT
  std::ofstream of("point_groups.g2");
  for (int ka=0; ka<3; ++ka)
    {
      if (pnt_groups[ka].size() > 0)
	{
	  of << "400 1 0 0" << std::endl;
	  of << pnt_groups[ka].size() << std::endl;
	  for (size_t ki=0; ki<pnt_groups[ka].size(); ++ki)
	    of << pnt_groups[ka][ki]->getPoint() << std::endl;
	}
    }
#endif
  vector<vector<RevEngPoint*> > sep_groups;
  size_t g_ix = 0;
  for (int ka=0; ka<3; ++ka)
    {
      if (pnt_groups[ka].size() > 1)
	{
	  shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
							edge_class_type_,
							pnt_groups[ka]));
	  vector<vector<RevEngPoint*> > curr_sep_groups;
	  reg->splitRegion(curr_sep_groups);
	  sep_groups.push_back(reg->getPoints());
	  if (curr_sep_groups.size() > 0)
	    sep_groups.insert(sep_groups.end(), curr_sep_groups.begin(),
			      curr_sep_groups.end());
	  for (; g_ix<sep_groups.size(); ++g_ix)
	    for (size_t kj=0; kj<sep_groups[g_ix].size(); ++kj)
	      sep_groups[g_ix][kj]->unsetRegion();
	}
      else if (pnt_groups[ka].size() == 1)
	pnt_groups[ka][0]->unsetRegion();
    }

#ifdef DEBUG_SEGMENT
  std::ofstream of2("sep_groups.g2");
  for (size_t kj=0; kj<sep_groups.size(); ++kj)
    {
      of2 << "400 1 0 0" << std::endl;
      of2 << sep_groups[kj].size() << std::endl;
      for (size_t ki=0; ki<sep_groups[kj].size(); ++ki)
	of2 << sep_groups[kj][ki]->getPoint() << std::endl;
    }
#endif

  int max_ix = -1;
  int num_pnts = 0;
  for (size_t ki=0; ki<sep_groups.size(); ++ki)
    {
      if ((int)sep_groups[ki].size() > num_pnts)
	{
	  num_pnts = (int)sep_groups[ki].size();
	  max_ix = (int)ki;
	}
    }

  for (size_t ki=0; ki<sep_groups.size(); ++ki)
    {
      if ((int)ki == max_ix)
	{
	  group_points_.clear();
	  group_points_ = sep_groups[max_ix];
	  for (size_t kj=0; kj<group_points_.size(); ++kj)
	    group_points_[kj]->setRegion(this);
	  updateInfo(tol, angtol);
	}
      else
	{
	  for (size_t kj=0; kj<sep_groups[ki].size(); ++kj)
	    sep_groups[ki][kj]->unsetRegion();
	  added_groups.push_back(sep_groups[ki]);
	}
    }

  return (added_groups.size() > 0);
}

//===========================================================================
void RevEngRegion::getPCA(double lambda[3], Point& eigen1, Point& eigen2,
			  Point& eigen3)
//===========================================================================
{
  return getPCA(group_points_, lambda, eigen1, eigen2, eigen3);
}

//===========================================================================
void RevEngRegion::getPCA(vector<RevEngPoint*>& points, double lambda[3], 
			  Point& eigen1, Point& eigen2, Point& eigen3)
//===========================================================================
{
  // PCA analysis of group points
  Vector3D xyz = points[0]->getPoint();
  Point pos0(xyz[0], xyz[1], xyz[2]);
  vector<Point> pnts;
  pnts.reserve(points.size());
  //double wgt = 1.0/(double)points.size();
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      pnts.push_back(pos);
    }

  double eigenvec[3][3];
  RevEngUtils::principalAnalysis(pos0, pnts, lambda, eigenvec);
  eigen1 = Point(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2]);
  eigen2 = Point(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
  eigen3 = Point(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
}


//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::surfApprox(vector<RevEngPoint*>& points,
						   const BoundingBox& bbox)
//===========================================================================
{
  // PCA analysis of given points to orient plane for parametrerization
  Vector3D xyz = points[0]->getPoint();
  Point pos0(xyz[0], xyz[1], xyz[2]);
  vector<Point> pnts;
  pnts.reserve(points.size());
  Point vec(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(group_points_.size());
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      pnts.push_back(pos);
      Point vec0 = points[ki]->minCurvatureVec();
      if (vec0*vec < 0.0)
	vec0 *= -1;
      vec += wgt*vec0;
    }
  double lambda[3];
  double eigenvec[3][3];
  RevEngUtils::principalAnalysis(pos0, pnts, lambda, eigenvec);
  Point eigen3(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
  Point ydir = eigen3.cross(vec);
  ydir.normalize_checked();
  Point xdir = ydir.cross(eigen3);
  xdir.normalize_checked();

  // Parameterize points based on planar surface
  pnts.push_back(pos0);
  vector<double> data;
  vector<double> param;
  RevEngUtils::parameterizeWithPlane(pnts, bbox, xdir, ydir,
				     data, param);

  // Surface approximation
  int order = 3;
  double belt = 0.1*bbox.low().dist(bbox.high());
  shared_ptr<SplineSurface> surf = RevEngUtils::surfApprox(data, 3, param, order,
							   order, order, order, belt);
return surf;
}

//===========================================================================
void RevEngRegion::collect(RevEngPoint *pt, RevEngRegion *prev)
//===========================================================================
{
  if (pt->hasRegion() && pt->region() != this)
    return;  // Cannot grow
  group_points_.push_back(pt);
  if (classification_type_ == CLASSIFICATION_UNDEF)
    return; // Cannot grow
  int type = pt->surfaceClassification(classification_type_);
  double meancurv = pt->meanCurvature();
  double gausscurv = pt->GaussCurvature();
  if (type == C1_UNDEF)  // SI_UNDEF == C1_UNDEF
    return; // Cannot grow

  vector<RevEngPoint*> grouped;
  grouped.push_back(pt);
  Vector3D xyz = pt->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  bbox_.addUnionWith(xyz2);
  if (normalcone_.dimension() == 0)
    normalcone_ = DirectionCone(pt->getLocFuncNormal());
  else
    normalcone_.addUnionWith(pt->getLocFuncNormal());
  if (normalcone2_.dimension() == 0)
    normalcone2_ = DirectionCone(pt->getTriangNormal());
  else
    normalcone2_.addUnionWith(pt->getTriangNormal());
  bool planar = planartype();
  for (size_t kj=0; kj<grouped.size(); ++kj)
    {
      vector<ftSamplePoint*> next = grouped[kj]->getNeighbours();
      for (size_t ki=0; ki<next.size(); ++ki)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[ki]);
	  if (!curr)
	    continue;  // Should not happen
	  if (curr->isOutlier())
	    continue;
	  if (curr->hasRegion() && curr->region() != prev)
	    continue;  // Already belonging to a segment
	  if (curr->isEdge(edge_class_type_))
	    continue;  // An edge point
	  int type2 = curr->surfaceClassification(classification_type_);
	  if (type2 != type)
	    continue;   // Different classification
	  if (classification_type_ == CLASSIFICATION_SHAPEINDEX &&
	      (!planar) && meancurv*curr->meanCurvature() < 0.0)
	    continue;  // Not compatible with one primary surface
	  double curv1 = meancurv*gausscurv;
	  double curv2 = curr->meanCurvature()*curr->GaussCurvature();
	  if (classification_type_ == CLASSIFICATION_SHAPEINDEX &&
	      (!planar) && curv1*curv2 < 0.0)
	    continue;  // Not compatible with one primary surface

	  // Add to region
	  curr->setRegion(this);
	  group_points_.push_back(curr);
	  xyz = curr->getPoint();
	  xyz2 = Point(xyz[0], xyz[1], xyz[2]);
	  bbox_.addUnionWith(xyz2);
	  normalcone_.addUnionWith(curr->getLocFuncNormal());
	  normalcone2_.addUnionWith(curr->getTriangNormal());

	  // Continue growing from this point
	  grouped.push_back(curr);
	}
    }
#ifdef DEBUG_COLLECT
  std::ofstream of("collected_group.g2");
  writeRegionPoints(of);
  vector<RevEngPoint*> bd_pts = extractBdPoints();
  of << "400 1 0 4 0 0 0 255" << std::endl;
  of << bd_pts.size() << std::endl;
  for (size_t kj=0; kj<bd_pts.size(); ++kj)
    of << bd_pts[kj]->getPoint() << std::endl;
#endif
  // Principal curvature summary
  updateInfo();
}

//===========================================================================
BoundingBox RevEngRegion::getParameterBox()
//===========================================================================
{
  BoundingBox parbox(2);
  if (associated_sf_.size() == 0)
    return parbox;   // No surface means that the points are not parameterized

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector2D par = group_points_[ki]->getPar();
      Point par2(par[0], par[1]);
      parbox.addUnionWith(par2);
    }
  return parbox;
}


//===========================================================================
RevEngPoint* RevEngRegion::seedPointPlane(int min_next, double rfac, double angtol)
//===========================================================================
{
  double min_in = 0;
  int min_ix = -1;
  int multfac = 10;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double local_len = group_points_[ki]->getMeanEdgLen();
      Point normal = group_points_[ki]->getLocFuncNormal();
      double radius = 2.0*rfac*local_len;
      vector<RevEngPoint*> nearpts;
      group_points_[ki]->fetchClosePoints2(radius, min_next, 5*min_next,
					   nearpts, this);

      // Count deviant points
      int deviant = 0;
      for (size_t kj=0; kj<nearpts.size(); ++kj)
      {
	Point curr_norm = nearpts[kj]->getLocFuncNormal();
	double ang = curr_norm.angle(normal);
	if (nearpts[kj]->region() != this || ang > angtol)
	  ++deviant;
      }

      if (deviant == 0)
	return group_points_[ki];   // Seed point found

      if ((int)nearpts.size() - deviant > min_in)
	{
	  min_in = (int)nearpts.size() - deviant;
	  min_ix = (int)ki;
	  if (min_in < multfac*min_next)
	    break;
	}
    }
  return (min_ix >= 0) ? group_points_[min_ix] : 0;
}

//===========================================================================
void RevEngRegion::growLocalPlane(Point mainaxis[3],
				  double tol, vector<RevEngPoint*>& plane_pts,
				  shared_ptr<Plane>& plane_out)
//===========================================================================
{
  // Get seed point for local grow
  int min_next = std::max(10, (int)group_points_.size()/100);
  double rfac = 3;
  double angtol = 0.1;
  RevEngPoint* seed = seedPointPlane(min_next, rfac, angtol);
  if (!seed)
    return;
 
#ifdef DEBUG_SEGMENT
  std::ofstream of("seed.g2");
  of << "400 1 0 4 200 0 55 255" << std::endl;
  of << "1" << std::endl;
  of << seed->getPoint() << std::endl;
#endif
  // Fetch nearby points belonging to the same region
  double local_len = seed->getMeanEdgLen();
  double radius = 2.0*rfac*local_len;
  vector<RevEngPoint*> nearpts;
  seed->fetchClosePoints2(radius, min_next, 5*min_next, nearpts, this);
  nearpts.insert(nearpts.begin(), seed);
  if ((int)nearpts.size() < min_next/2)
    return;

  // Approximate with plane
  shared_ptr<Plane> plane = computePlane(nearpts, avnorm_, mainaxis);

  double eps = 1.0e-6;
  double maxdist, avdist;
  int num_inside, num2_inside;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(nearpts.begin(), nearpts.end(),
			  plane, tol, maxdist, avdist, num_inside, num2_inside,
			  inpt, outpt, parvals, dist_ang, angtol);

  if (maxdist > tol)
    return;   // For the time being

  
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    nearpts[ki]->setVisited();
  size_t prev_size = nearpts.size();
  while (true)
    {
      for (size_t ki=0; ki<nearpts.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = nearpts[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	      if ((!pt->hasRegion()) || pt->region() != this)
		continue;
	      if (pt->visited())
		continue;
	      Vector3D xyz = pt->getPoint();

	      double upar, vpar, dist;
	      Point close;
	      plane->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
	      pt->setVisited();
	      if (dist < tol)
		nearpts.push_back(pt);
	    }
	}

      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  if (std::find(nearpts.begin(), nearpts.end(), group_points_[ki]) == nearpts.end())
	    group_points_[ki]->unsetVisited();
	}

      if (nearpts.size() == prev_size)
	break;
      
      shared_ptr<Plane> plane2 = computePlane(nearpts, avnorm_, mainaxis);

      double maxdist2, avdist2;
      int num_inside2, num2_inside2;
      vector<pair<double, double> > dist_ang2;
      vector<double> parvals2;
     RevEngUtils::distToSurf(nearpts.begin(), nearpts.end(),
			     plane2, tol, maxdist2, avdist2, num_inside2,
			     num2_inside2, inpt,
			     outpt, parvals2, dist_ang2, angtol);
      if (maxdist2 > tol)
	break;
      
      prev_size = nearpts.size();
      plane = plane2;
      maxdist = maxdist2;
      avdist = avdist2;
      num_inside = num_inside2;
      parvals = parvals2;
      dist_ang = dist_ang2;
    }

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    group_points_[ki]->unsetVisited();
#ifdef DEBUG_EXTRACT  
  std::ofstream of2("plane_pts.g2");
  of2 << "400 1 0 4 100 155 0 255" << std::endl;
  of2 << nearpts.size() << std::endl;
  for (size_t ki=0; ki<nearpts.size(); ++ki)
    of2 << nearpts[ki]->getPoint() << std::endl;

  plane->writeStandardHeader(of2);
  plane->write(of2);
#endif
  plane_pts = nearpts;
  plane_out = plane;

  int stop_break = 1;

}

//===========================================================================
void RevEngRegion::segmentByPlaneGrow(Point mainaxis[3], double tol,
				      double angtol, int min_pt, 
				      vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      vector<HedgeSurface*>& prevsfs,
				      vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // Extract group of planar points and associated plane
  int min_nmb_pts = std::max((int)group_points_.size()/100, 10);
  vector<RevEngPoint*> plane_pts;
  shared_ptr<Plane> plane;
  growLocalPlane(mainaxis, tol, plane_pts, plane);
  if ((int)plane_pts.size() < min_nmb_pts || plane_pts.size() == group_points_.size())
    return;

  // Fetch accuracy information
   double maxd, avd;
   int num_in, num2_in;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(plane_pts.begin(), plane_pts.end(),
			  plane, tol, maxd, avd, num_in, num2_in, inpt,
			  outpt, parvals, dist_ang, angtol);
  
  // Extract remaining points
  for (size_t ki=0; ki<plane_pts.size(); ++ki)
    plane_pts[ki]->setVisited();

  vector<RevEngPoint*> remaining_pts;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (!group_points_[ki]->visited())
	{
	  remaining_pts.push_back(group_points_[ki]);
	  removePoint(group_points_[ki]);
	}
    }

  for (size_t ki=0; ki<plane_pts.size(); ++ki)
    plane_pts[ki]->unsetVisited();

  if (hasSurface())
    {
      // No longer valid
      prevsfs.insert(prevsfs.end(), associated_sf_.begin(), associated_sf_.end());
      clearSurface();
    }
  
  if ((int)plane_pts.size() >= min_pt)
    {
      // Register current plane
      int sf_flag = defineSfFlag(0, tol, num_in, num2_in,
				 avd, false);
      for (size_t ki=0; ki<plane_pts.size(); ++ki)
	{
	  plane_pts[ki]->setPar(Vector2D(parvals[2*ki],parvals[2*ki+1]));
	  plane_pts[ki]->setSurfaceDist(dist_ang[ki].first, dist_ang[ki].second);
	}
      setAccuracy(maxd, avd, num_in, num2_in);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(plane, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }

  // Store as base surface and primary surface
  setBaseSf(plane, maxd, avd, num_in, num2_in);

  updateInfo(tol, angtol);

  // Distribute remaining points into groups and create regions
  shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						 edge_class_type_,
						 remaining_pts));
  vector<vector<RevEngPoint*> > connected;
  reg->splitRegion(connected);
  out_groups.push_back(reg->getPoints());
  if (connected.size() > 0)
    out_groups.insert(out_groups.end(), connected.begin(), connected.end());
}



//===========================================================================
bool RevEngRegion::extractPlane(Point mainaxis[3],
				double tol, int min_pt, int min_pt_reg,
				double angtol, int prefer_elementary,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

#ifdef DEBUG
  if (normalcone_.greaterThanPi())
    {
      std::cout << "Greater than pi" << std::endl;
      std::ofstream ofpi("pi_region.g2");
      writeRegionInfo(ofpi);
    }
#endif
  
  bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;

  // vector<RevEngPoint*> in, out;
  // analysePlaneProperties(normal, angtol, in, out);
  
  // shared_ptr<Plane> surf1(new Plane(pos, normal1));
  // shared_ptr<Plane> surf(new Plane(pos, normal));

  // Point pos2 = pos;
  // Point normal2;
  // vector<pair<vector<RevEngPoint*>::iterator,
  // 	      vector<RevEngPoint*>::iterator> > group;
  // group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  // RevEngUtils::computePlane(group, pos2, normal2);
  // shared_ptr<Plane> surf2(new Plane(pos2, normal2));
  
  // vector<RevEngPoint*> in3, out3;
  // analysePlaneProperties(normal3, angtol, in3, out3);

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);
  
  shared_ptr<Plane> surf3 = computePlane(group_points_, avnorm_, mainaxis); //(new Plane(pos3, normal3));

  // Check accuracy
  // double maxdist, avdist;
  // double maxdist1, avdist1;
  // double maxdist2, avdist2;
  double maxdist3, avdist3;
  int num_inside3, num2_inside3; //, num_inside1, num_inside2, num_inside3;
  vector<RevEngPoint*> inpt, outpt;
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf, tol, maxdist, avdist, num_inside);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf1, tol, maxdist1, avdist1, num_inside1);
  // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  // 			  surf2, tol, maxdist2, avdist2, num_inside2);
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf3, tol, maxdist3, avdist3, num_inside3, num2_inside3,
			  inpt, outpt, parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_plane.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif
  
  if (inpt.size() > group_points_.size()/3 && (int)inpt.size() > min_nmb &&
      (normalcone_.angle() > angtol ||
       normalcone_.centre().angle(surf3->getNormal()) > angtol) &&
      (!(basesf_.get() && basesf_->instanceType() == Class_Cylinder)))
    {
      vector<RevEngPoint*> ang_points;
      double dtol = std::min(0.5*tol, 1.5*avdist3);
      identifyAngPoints(dist_ang, angtol, dtol, ang_points);
      if (ang_points.size() < group_points_.size()/2)
	{
	  extractSpesPoints(ang_points, out_groups, true);
	  
	  if (out_groups.size() > 0)
	    {
	      // Ensure connected region
	      vector<vector<RevEngPoint*> > separate_groups;
	      splitRegion(separate_groups);
	      if (separate_groups.size() > 0)
		{
		  out_groups.insert(out_groups.end(), separate_groups.begin(),
				    separate_groups.end());
		}
	  
	      shared_ptr<Plane> plane_in =
		computePlane(group_points_, avnorm_, mainaxis);
	      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
	      vector<pair<double, double> > dist_ang_in;
	      vector<double> parvals_in;
	      double maxd_in, avd_in;
	      int num2_in, num2_in2;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      plane_in, tol, maxd_in, avd_in, num2_in, num2_in2,
				      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
#ifdef DEBUG_EXTRACT  
  	      std::ofstream ofd2("in_out_plane2.g2");
	      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	      ofd2 << inpt_in.size() << std::endl;
	      for (size_t kr=0; kr<inpt_in.size(); ++kr)
		ofd2 << inpt_in[kr]->getPoint() << std::endl;
	      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	      ofd2 << outpt_in.size() << std::endl;
	      for (size_t kr=0; kr<outpt_in.size(); ++kr)
		ofd2 << outpt_in[kr]->getPoint() << std::endl;
#endif
	      // if (num2_in > num_inside3 || avd_in < avdist3)
	      // 	{
	      std::swap(surf3, plane_in);
	      std::swap(num_inside3, num2_in);
	      std::swap(num2_inside3, num2_in2);
	      std::swap(avdist3, avd_in);
	      std::swap(maxdist3, maxd_in);
	      std::swap(parvals, parvals_in);
	      std::swap(dist_ang, dist_ang_in);
	      // 	  std::cout << "Plane swap" << std::endl;
	      // 	  std::cout << group_points_.size() << " " << num2_in << " ";
	      // std::cout<< num_inside3 << " " << avd_in << " " << avdist3;
	      // 	  std::cout << " " << maxd_in << " " << maxdist3 << std::endl;
	      // }
	    }
	}
    }

  Point low = bbox_.low();
  Point high = bbox_.high();
  //double len = low.dist(high);
  //surf3->setParameterBounds(-len, -len, len, len);
#ifdef DEBUG_EXTRACT  
  std::ofstream plane("curr_plane.g2");
  surf3->writeStandardHeader(plane);
  surf3->write(plane);
#endif
  //int num = (int)group_points_.size();
  int sf_flag = defineSfFlag(0, tol, num_inside3, num2_inside3,
			     avdist3, false);
  if (sf_flag < ACCURACY_POOR)
    {
      found = true;
      for (size_t kh=0; kh<group_points_.size(); ++kh)
	{
	  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	}
      setAccuracy(maxdist3, avdist3, num_inside3, num2_inside3);

#ifdef DEBUG
      std::cout << "Plane. N1: " << num << ", N2: " << num_inside3 << ", max: " << maxdist3 << ", av: " << avdist3 << std::endl;
#endif
      
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf3, this));
      for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	prevsfs.push_back(associated_sf_[kh]);
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }
  if (!basesf_.get() ||
      (num_inside3 >= num_in_base_ && avdist3 < avdist_base_))
    setBaseSf(surf3, maxdist3, avdist3, num_inside3, num2_inside3);

  return found;
}


//===========================================================================
shared_ptr<Plane> RevEngRegion::computePlane(vector<RevEngPoint*>& points,
					     const Point& norm_dir, Point mainaxis[3])
//===========================================================================
{
  Point normal1 = normalcone_.centre();
  Point normal(0.0, 0.0, 0.0);
  Point pos(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point curr = points[ki]->getLocFuncNormal();
      Vector3D xyz = points[ki]->getPoint();
      normal += wgt*curr;
      pos += wgt*Point(xyz[0], xyz[1], xyz[2]);
    }
  
  ImplicitApprox impl;
  impl.approx(points, 1);
  Point pos3, normal3;
  bool found = impl.projectPoint(pos, normal, pos3, normal3);
  if (!found)
    {
      pos3 = pos;
      normal3 = normal;
      normal3.normalize();
    }
  if (normal3*norm_dir < 0.0)
    normal3 *= -1.0;

  // Define x-axis
  int ix = -1;
  double minang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(normal3);
      ang = std::min(ang, M_PI-ang);
      if (ang < minang)
	{
	  minang = ang;
	  ix = ka;
	}
    }

  Point Cy = mainaxis[(ix+1)%3].cross(normal3);
  Point Cx = normal3.cross(Cy);
  
  shared_ptr<Plane> surf(new Plane(pos3, normal3, Cx));
  Point low = bbox_.low();
  Point high = bbox_.high();
  //double len = low.dist(high);
  //surf->setParameterBounds(-len, -len, len, len);
  
  bool plane_project = false;
  if (plane_project)
    {
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > group;
      group.push_back(std::make_pair(points.begin(), points.end()));
      Point vec1, vec2;
      surf->getSpanningVectors(vec1, vec2);
      vector<Point> projected1, projected2;
      double maxdp1, avdp1, maxdp2, avdp2;
      RevEngUtils::projectToPlane(group, vec1, pos3, projected1, maxdp1, avdp1);
      RevEngUtils::projectToPlane(group, vec2, pos3, projected2, maxdp2, avdp2);

#ifdef DEBUG_EXTRACT
      std::ofstream ofp3("plane_project.g2");
      ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
      ofp3 << projected1.size() << std::endl;
      for (size_t kr=0; kr<projected1.size(); ++kr)
	ofp3 << projected1[kr] << std::endl;
      ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
      ofp3 << projected2.size() << std::endl;
      for (size_t kr=0; kr<projected2.size(); ++kr)
	ofp3 << projected2[kr] << std::endl;
#endif
    }
  return surf;
}

//===========================================================================
void RevEngRegion::getDistAndAng(vector<pair<double,double> >& distang)
//===========================================================================
{
  distang.resize(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double dist, ang;
      group_points_[ki]->getSurfaceDist(dist, ang);
      distang[ki] = std::make_pair(dist, ang);
    }
}

//===========================================================================
bool RevEngRegion::possiblePlane(double angtol, double inlim)
//===========================================================================
{
  if (planartype())
    return true;

  double wgt = 1.0/(double)group_points_.size();
  Point avnorm(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getLocFuncNormal();
      avnorm += wgt*curr;
    }

  int nmb_in = 0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getLocFuncNormal();
      double ang = curr.angle(avnorm);
      if (ang <= angtol)
	++nmb_in;
    }

  double frac = (double)nmb_in/(double)group_points_.size();
  return (frac >= inlim);
}

//===========================================================================
vector<RevEngRegion*> RevEngRegion::fetchAdjacentPlanar()
//===========================================================================
{
  vector<RevEngRegion*> adjacent_planar;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  if (sf->instanceType() == Class_Plane)
	    adjacent_planar.push_back(*it);
	}
    }

  return adjacent_planar;
}

//===========================================================================
vector<RevEngRegion*> RevEngRegion::fetchAdjacentCylindrical()
//===========================================================================
{
  vector<RevEngRegion*> adjacent_cylindrical;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  if (sf->instanceType() == Class_Cylinder)
	    adjacent_cylindrical.push_back(*it);
	}
    }

  return adjacent_cylindrical;
}

//===========================================================================
Point RevEngRegion::directionFromAdjacent(double angtol)
//===========================================================================
{
  Point dir;
  vector<std::pair<Point,int> > all_dir;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  if (sf->instanceType() == Class_Plane || sf->instanceType() == Class_Cylinder)
	    {
	      shared_ptr<ElementarySurface> elem =
		dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	      Point curr_dir = elem->direction();
	      size_t kj;
	      for (kj=0; kj<all_dir.size(); ++kj)
		{
		  double ang = curr_dir.angle(all_dir[kj].first);
		  if (M_PI-ang < ang)
		    {
		      curr_dir *= -1.0;
		      ang = M_PI - ang;
		    }
		  if (ang < angtol)
		    {
		      int num = (*it)->numPoints();
		      double fac = 1.0/(double)(all_dir[kj].second+num);
		      curr_dir = fac*(all_dir[kj].second*all_dir[kj].first + num*curr_dir);
		      break;
		    }
		}
	      if (kj == all_dir.size())
		all_dir.push_back(std::make_pair(curr_dir, (*it)->numPoints()));
	    }
	}
    }

  for (size_t ki=0; ki<all_dir.size(); ++ki)
    for (size_t kj=ki+1; kj<all_dir.size(); ++kj)
      if (all_dir[kj].second > all_dir[ki].second)
	std::swap(all_dir[ki], all_dir[kj]);

  if (all_dir.size() > 0)
    dir = all_dir[0].first;

  return dir;
}

//===========================================================================
bool RevEngRegion::potentialBlend(double angtol)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      ClassType type1 = adj_elem[ki].first->instanceType();
      Point dir1 = adj_elem[ki].first->direction();
      if (type1 != Class_Plane && type1 != Class_Cylinder)
	continue;
      for (size_t kj=ki+1; kj<adj_elem.size(); ++kj)
	{
	  ClassType type2 = adj_elem[kj].first->instanceType();
	  if (type2 != Class_Plane && type2 != Class_Cylinder)
	    continue;
	  if (type1 == Class_Cylinder && type2 == Class_Cylinder)
	    continue;
	  Point dir2 = adj_elem[kj].first->direction();
	  double ang = dir1.angle(dir2);
	  if (type1 == Class_Plane && type2 == Class_Plane)
	    ang = std::min(ang, M_PI-ang);
	  else
	    ang = fabs(pihalf - ang);
	  if (ang > angtol)
	    return true;
	}
    }
  return false;
}

//===========================================================================
void RevEngRegion::neighbourBlends(vector<shared_ptr<CurveOnSurface> >& cvs,
				    double width, double tol,
				    vector<RevEngRegion*>& new_blends)
//===========================================================================
{
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasAssociatedBlend())
	continue;
      if ((*it)->hasRevEdges())
	continue;
      int in_blend = 0;
      if ((*it)->isInBlend(cvs, width, tol, in_blend))
	new_blends.push_back(*it);
    }
}

//===========================================================================
bool RevEngRegion::isInBlend(vector<shared_ptr<CurveOnSurface> >& cvs,
			     double width, double tol, int& in_blend)
//===========================================================================
{
  double lim = 0.9;
  in_blend = 0;
  double tmin = cvs[0]->startparam();
  double tmax = cvs[cvs.size()-1]->endparam();
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double tpar, dist=std::numeric_limits<double>::max();;
      Vector3D xyz = group_points_[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      Point close;
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  cvs[kj]->closestPoint(pt, cvs[kj]->startparam(), cvs[kj]->endparam(),
				tpar, close, dist);
	  if (tpar > tmin && tpar < tmax && dist < width+tol)
	    in_blend++;
	}
    }
  if ((double)in_blend >= lim*(double)group_points_.size())
    return true;
  else
    return false;
}

//===========================================================================
void RevEngRegion::getAdjacentElemInfo(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem,
				       vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem_base)
//===========================================================================
{
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      bool found = false;
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> sf = (*it)->getSurface(0)->surface();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	  if (elem.get())
	    {
	      adj_elem.push_back(std::make_pair(elem,*it));
	      found = true;
	    }
	}
      
      if ((*it)->hasBaseSf() && (!found))
	{
	  shared_ptr<ParamSurface> sf = (*it)->getBase();
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(sf);
	  if (elem.get())
	    {
	      adj_elem_base.push_back(std::make_pair(elem, *it));
	      found = true;
	    }
	}
      
    }

#ifdef DEBUG_ADJACENT
  std::ofstream of("adj_sf_pts.g2");
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    adj_elem[ki].second->writeRegionPoints(of);
#endif
}

//===========================================================================
bool RevEngRegion::possibleCylinder(double angtol, double inlim)
//===========================================================================
{
  // if (cylindertype())
  //   return true;

  double wgt = 1.0/(double)group_points_.size();
  Point avvec = wgt*group_points_[0]->minCurvatureVec();
  for (size_t ki=1; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      if (curr*avvec < 0.0)
	curr *= -1;
      avvec += wgt*curr;
    }

  int nmb_in = 0;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      double ang = curr.angle(avvec);
      ang = std::min(ang, M_PI-ang);
      if (ang <= angtol)
	++nmb_in;
    }

  double frac = (double)nmb_in/(double)group_points_.size();
  return (frac >= inlim);
}

//===========================================================================
bool RevEngRegion::possibleCone(double tol, double inlim)
//===========================================================================
{
  return true;  // For the time being
}

//===========================================================================
bool RevEngRegion::possibleTorus(double tol, double inlim)
//===========================================================================
{
  // if (planartype() || cylindertype())
  //   return false;
  // Compute curvature variations
  double k1mean = 0.0, k2mean = 0.0;
  double k1_1 = 0.0, k1_2 = 0.0, k2_1 = 0.0, k2_2 = 0.0;
  double d1 = 0.0, d2 = 0.0;
  double wgt = 1.0/(double)group_points_.size();
  Point vec(0.0, 0.0, 0.0);

  double eps = 1.0e-3;
  vector<Vector3D> cneg, cpos, czero;
  double k1min = std::numeric_limits<double>::max();
  double k1max = std::numeric_limits<double>::lowest();
  double k2min = std::numeric_limits<double>::max();
  double k2max = std::numeric_limits<double>::lowest();
  int k1pos = 0, k1neg = 0, k2pos = 0, k2neg = 0;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmin = group_points_[kr]->minPrincipalCurvature();
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      d1 += kmin;
      d2 += kmax;
      k1_1 += kmin*kmin;
      k2_1 += kmax*kmax;
      k1_2 += fabs(kmin);
      k2_2 += fabs(kmax);
      k1mean += wgt*kmin;
      k2mean += wgt*kmax;
      k1min = std::min(k1min, kmin);
      k1max = std::max(k1max, kmin);
      k2min = std::min(k2min, kmax);
      k2max = std::max(k2max, kmax);
      if (kmin < 0.0)
	k1neg++;
      else
	k1pos++;
      if (kmax < 0.0)
	k2neg++;
      else
	k2pos++;

      Point norm = group_points_[kr]->getLocFuncNormal();
      Point minvec = group_points_[kr]->minCurvatureVec();
      Point temp = norm.cross(minvec);
      temp.normalize_checked();
      vec += wgt*temp;
      if (group_points_[kr]->GaussCurvature()*group_points_[kr]->meanCurvature() < -eps)
	cneg.push_back(group_points_[kr]->getPoint());
      else if (group_points_[kr]->GaussCurvature()*group_points_[kr]->meanCurvature() > eps)
	cpos.push_back(group_points_[kr]->getPoint());
      else
	czero.push_back(group_points_[kr]->getPoint());
    }
  //double vark1 = ((double)group_points_.size()*k1_1)/fabs(d1) - fabs(k1_2);
  //double vark2 = ((double)group_points_.size()*k2_1)/fabs(d2) - fabs(k2_2);

  double klim1 = 0.1*fabs(k1mean);
  double klim2 = 0.1*fabs(k2mean);
  int nmb1=0, nmb2=0;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      double kmin = group_points_[kr]->minPrincipalCurvature();
      double kmax = group_points_[kr]->maxPrincipalCurvature();
      if (fabs(kmin - k1mean) < klim1)
	++nmb1;
     if (fabs(kmax - k2mean) < klim2)
	++nmb2;
    }

  //double frac1 = (double)nmb1/(double)group_points_.size();
  //double frac2 = (double)nmb2/(double)group_points_.size();

#ifdef DEBUG_EXTRACT
  std::ofstream of("posnegcurv.g2");
  of << "400 1 0 4 0 100 155 255" << std::endl;
  of << cneg.size() << std::endl;
  for (size_t kr=0; kr<cneg.size(); ++kr)
    of << cneg[kr] << std::endl;
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << cpos.size() << std::endl;
  for (size_t kr=0; kr<cpos.size(); ++kr)
    of << cpos[kr] << std::endl;
  of << "400 1 0 4 200 55 0 255" << std::endl;
  of << czero.size() << std::endl;
  for (size_t kr=0; kr<czero.size(); ++kr)
    of << czero[kr] << std::endl;
#endif
  // Something with variation in curvature
  return true;  // For the time being
}

//===========================================================================
void RevEngRegion::analysePlaneProperties(Point avnorm, double angtol,
					  vector<RevEngPoint*>& in,
					  vector<RevEngPoint*> out)
//===========================================================================
{
  vector<double> midang(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->getLocFuncNormal();
      double ang = curr.angle(avnorm);
      midang[ki] = ang;
      if (ang <= angtol)
	in.push_back(group_points_[ki]);
      else
	out.push_back(group_points_[ki]);
    }

  std::sort(midang.begin(), midang.end());
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::analyseNormals(double tol, Point& normal, Point& centre,
				  double& radius) //double& beta)
//===========================================================================
{
#ifdef DEBUG0
  std::ofstream of2("curr_normals.g2");
  of2 << "400 1 0 4 100  0 155 255" << std::endl;
  of2 << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getLocFuncNormal();
      of2 << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of2);
  sph.write(of2);
#endif
  
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);

  Point pnt(0.0, 0.0, 0.0);
  Point pos(0.0, 0.0, 0.0);
  vector<Point> vec(group_points_.size());
  double wgt = 1.0/(double)(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      pos += wgt*Point(xyz[0], xyz[1], xyz[2]);
      vec[ki] = group_points_[ki]->getLocFuncNormal();
      pnt += wgt*vec[ki];
    }
  Point tmpnorm;
  double maxang = 0.0;
  size_t ix = 0;
  for (size_t kj=1; kj<vec.size(); ++kj)
    {
      double ang = vec[0].angle(vec[kj]);
      ang = std::min(ang, M_PI-ang);
      if (ang > maxang)
	{
	  maxang = ang;
	  ix = kj;
	}
    }
  tmpnorm = vec[0].cross(vec[ix]);
  tmpnorm.normalize();

  ImplicitApprox impl;
  impl.approxPoints(vec, 1);

  //Point pos; //, normal;
  bool found = impl.projectPoint(pos, tmpnorm, pnt, normal);
  if (!found)
    {
      pnt = pos;
      normal = tmpnorm;
    }
  shared_ptr<Plane> surf(new Plane(pnt, axis));

#ifdef DEBUG0
  Point origo(0.0, 0.0,0.0);
  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  of2 << "1" << std::endl;
  of2 << origo << " " << axis << std::endl;
#endif
  
  double maxdist, avdist;
  int num_inside;
  //int num = (int)vec.size();
  vector<double> distance;
  RevEngUtils::distToSurf(vec, surf, tol, maxdist, avdist, num_inside,
			  distance);

  //double radius;
  RevEngUtils::computeRadius(vec, axis, Cx, Cy, radius);
  centre = pnt;
  // double r2 = std::min(radius, 1.0);
  // double phi = acos(r2);
  // double delta = r2*tan(phi);
  // double alpha = (delta > 1.0e-10) ? atan((1.0-r2)/delta) : M_PI;
  // beta = M_PI - alpha;
  int stop_break = 1;

}

//===========================================================================
bool RevEngRegion::feasiblePlane(double zero_H, double zero_K) const
//===========================================================================
{
  double angtol = 0.2;
  double in_lim = 0.8;
  double in_lim2 = 0.5;
  double fac = 5.0;
  if (basesf_ && basesf_->instanceType() == Class_Plane)
    return true;

  if (normalcone_.angle() < angtol && (!normalcone_.greaterThanPi()))
    return true;

  int nmb_in = 0;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal1 = group_points_[kr]->getLocFuncNormal();
      Point normal2 = group_points_[kr]->getTriangNormal();
      if (avnorm_.angle(normal1) <= angtol ||
	  avnorm2_.angle(normal2) <= angtol)
	nmb_in++;
    }
  double in_frac = (double)nmb_in/(double)group_points_.size();
  if (in_frac > in_lim)
    return true;

  if (in_frac < in_lim2)
    return false;
  
  if (std::max(MAH_, MAK_) > fac*std::max(zero_H, zero_K))
    return false;

  return true;
}

//===========================================================================
bool RevEngRegion::feasibleCylinder(double zero_H, double zero_K) const
//===========================================================================
{
  double fac1 = 5.0;
  double fac2 = 0.5;
  if (basesf_ && basesf_->instanceType() == Class_Cylinder)
    return true;
  if (MAK_ > fac1*zero_H)
    return false;
  if (MAK_ > fac2*MAH_)
    return false;
  // if (fabs(avK_) < fac2*MAK_)
  //   return false;

  return true;
}


//===========================================================================
bool RevEngRegion::extractCylinder(double tol, int min_pt, int min_pt_reg,
				   double angtol, int prefer_elementary,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups,
				   bool& repeat)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
    return false;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  vector<vector<RevEngPoint*> > configs;
  shared_ptr<Cone> cone;
  shared_ptr<Cylinder> cyl = computeCylinder(group_points_, tol);
#ifdef DEBUG_CYL
  std::ofstream ofs("cyl.g2");
  shared_ptr<Cylinder> cyl2(cyl->clone());
  double diag = bbox_.low().dist(bbox_.high());
  cyl2->setParamBoundsV(-0.5*diag,0.5*diag);
  cyl2->writeStandardHeader(ofs);
  cyl2->write(ofs);
#endif
  
  // Check accuracy
  double maxd, avd; 
  int num2, num2_2; 
  vector<RevEngPoint*> inpt, outpt; 
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  cyl, tol, maxd, avd, num2, num2_2, inpt, outpt, 
			  parvals, dist_ang, angtol);

#ifdef DEBUG_CYL
  std::ofstream ofd("in_out_cyl.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif

  int num_in_lin, num_in_cub;
  double avd_lin, avd_cub;
  analyseCylRotate(cyl, tol, avd, num2, avd_lin, num_in_lin,
		   avd_cub, num_in_cub, cone);
  analyseCylProject(cyl, tol, configs);
  if (cone.get())
    {
      double maxd2, avd2;
      int num3, num3_2;
      vector<RevEngPoint*> inpt2, outpt2;
      vector<pair<double, double> > dist_ang2;
      vector<double> parvals2;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      cone, tol, maxd2, avd2, num3, num3_2, inpt2, outpt2,
			      parvals2, dist_ang2, angtol);
#ifdef DEBUG_CYL
      std::ofstream ofsc("cone_cyl.g2");
      shared_ptr<Cone> cone2(cone->clone());
      double diag = bbox_.low().dist(bbox_.high());
      cone2->setParamBoundsV(-0.5*diag,0.5*diag);
      cone2->writeStandardHeader(ofsc);
      cone2->write(ofsc);
      std::ofstream ofdc("in_out_cone_cyl.g2");
      ofdc << "400 1 0 4 155 50 50 255" << std::endl;
      ofdc << inpt2.size() << std::endl;
      for (size_t kr=0; kr<inpt2.size(); ++kr)
	ofdc << inpt2[kr]->getPoint() << std::endl;
      ofdc << "400 1 0 4 50 155 50 255" << std::endl;
      ofdc << outpt.size() << std::endl;
      for (size_t kr=0; kr<outpt2.size(); ++kr)
	ofdc << outpt2[kr]->getPoint() << std::endl;
#endif
      int stop_cone = 1;
    }
  
  if (configs.size() > 1 && (!hasSurface()))
    {
      // Split point group and try again
      int keep_ix = -1;
      int keep_nmb = 0;
      for (size_t ki=0; ki<configs.size(); ++ki)
	if ((int)configs[ki].size() > keep_nmb)
	  {
	    keep_nmb = (int)configs[ki].size();
	    keep_ix = (int)ki;
	  }
      
      for (size_t ki=0; ki<configs.size(); ++ki)
	{
	  if ((int)ki == keep_ix)
	    continue;

	  extractSpesPoints(configs[ki], out_groups);
	}

      // Check that the remaing point cloud is connected
      splitRegion(out_groups);

      if (hasSurface())
	{
	  for (size_t ki=0; ki<associated_sf_.size(); ++ki)
	    prevsfs.push_back(associated_sf_[ki]);
	  clearSurface();
	}
      
      updateInfo(tol, angtol);
      repeat = true;
    }
  else
    {
      int num = (int)group_points_.size();
      if (num2 > min_pt && num2 > num/2)
	{
	  // Check for deviant points at the boundary
	  vector<RevEngPoint*> dist_points;
	  identifyDistPoints(dist_ang, tol, maxd, avd, dist_points);
	  extractSpesPoints(dist_points, out_groups, true);
	  if (out_groups.size() > 0)
	    {
	      // Ensure connected region
	      vector<vector<RevEngPoint*> > separate_groups;
	      splitRegion(separate_groups);
	      if (separate_groups.size() > 0)
		{
		  out_groups.insert(out_groups.end(), separate_groups.begin(),
				    separate_groups.end());
		}
	  
	      shared_ptr<Cylinder> cyl_in = computeCylinder(group_points_, tol);
	      vector<RevEngPoint*> inpt_in, outpt_in; 
	      vector<pair<double, double> > dist_ang_in;
	      vector<double> parvals_in;
	      double maxd_in, avd_in;
	      int num2_in, num2_in2;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      cyl_in, tol, maxd_in, avd_in, num2_in, num2_in2,
				      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);

	      
#ifdef DEBUG_CYL
	      std::ofstream ofd2("in_out_cylinder2.g2");
	      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	      ofd2 << inpt_in.size() << std::endl;
	      for (size_t kr=0; kr<inpt_in.size(); ++kr)
		ofd2 << inpt_in[kr]->getPoint() << std::endl;
	      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	      ofd2 << outpt_in.size() << std::endl;
	      for (size_t kr=0; kr<outpt_in.size(); ++kr)
		ofd2 << outpt_in[kr]->getPoint() << std::endl;
#endif
	      std::swap(cyl, cyl_in);
	      std::swap(num2, num2_in);
	      std::swap(num2_2, num2_in2);
	      std::swap(avd, avd_in);
	      std::swap(maxd, maxd_in);
	      std::swap(parvals, parvals_in);
	      std::swap(dist_ang, dist_ang_in);
	    }
	}
	  
      double maxd_init, avd_init;
      int num_init, num_init2;
      getAccuracy(maxd_init, avd_init, num_init, num_init2);
      int sf_flag = defineSfFlag(0, tol, num2, num2_2,
				 avd, true);
      if (sf_flag < ACCURACY_POOR)
	{
	  bool OK = true;
	  double acc_fac = 1.5;
	  if (associated_sf_.size() > 0)
	    {
	      int sfcode;
	      int sftype = associated_sf_[0]->instanceType(sfcode);
	      double ang = (sftype == Class_Plane) ?
		normalcone_.angle() : M_PI;
	      double ang_lim = 0.1*M_PI;
	      
	      // Check with current approximating surface
	      if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
		OK = false;
	      else if (prefer_elementary == ALWAYS_ELEM ||
		       prefer_elementary == PREFER_ELEM)
		{
		  if (!(ang > ang_lim && (num2 < num_init ||
					  (avd < avd_init &&
					   num2 < acc_fac*num_init))))
		    OK = false;
		}
	      else
		{
		  if (!(num2 < num_init ||
			(avd < avd_init && num2 < acc_fac*num_init)))
		    OK = true;
		}
	      if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
		OK = true;
	    }

	  if (OK)
	    {
	      found = true;
	      for (size_t kh=0; kh<group_points_.size(); ++kh)
		{
		  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
		  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
		}
	      setAccuracy(maxd, avd, num2, num2_2);
      
	      // // Limit cylinder with respect to bounding box
	      // double gap = 1.0e-6;
	      // Point xdir(1.0, 0.0, 0.0);
	      // Point ydir(0.0, 1.0, 0.0);
	      // Point zdir(0.0, 0.0, 1.0);
	      // Point bbdiag = high - low;
	      // CompositeModelFactory factory(gap, gap, 10.0*gap, 0.01, 0.05);
	      // shared_ptr<SurfaceModel> boxmod(factory.createFromBox(low-2.0*bbdiag, xdir, ydir, 
	      // 							    5*bbdiag[0], 5*bbdiag[1],
	      // 							    5*bbdiag[2]));
	      // vector<shared_ptr<ParamSurface> > sfs;
	      // sfs.push_back(cyl);
	      // shared_ptr<SurfaceModel> cylmod(new SurfaceModel(gap, gap, 10.0*gap, 0.01,
	      // 						       0.05, sfs));
	      // vector<shared_ptr<SurfaceModel> > divcyl = cylmod->splitSurfaceModels(boxmod);

#ifdef DEBUG
	      std::cout << "Cylinder. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
#endif
	      shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl, this));
	      for (size_t kh=0; kh<associated_sf_.size(); ++kh)
		prevsfs.push_back(associated_sf_[kh]);
	      setHedge(hedge.get());
	      hedgesfs.push_back(hedge);
	      setSurfaceFlag(sf_flag);
	      // for (int ka=0; ka<divcyl[0]->nmbEntities(); ++ka)
	      // 	{
	      // 	  shared_ptr<ParamSurface> cyl2 = divcyl[0]->getSurface(ka);
	      // 	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
	      // 	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	      // 	    prevsfs.push_back(associated_sf_[kh]);
	      // 	  associated_sf_.push_back(hedge.get());
	      // 	  hedgesfs.push_back(hedge);
	  
	    }
	}
    }
  if (!basesf_.get() ||
      (num2 >= num_in_base_ && avd < avdist_base_))
    setBaseSf(cyl, maxd, avd, num2, num2_2);

  return found;
}

//===========================================================================
shared_ptr<Cylinder>
RevEngRegion::computeCylinder(vector<RevEngPoint*>& points, double tol)
//===========================================================================
{
  // Cylinder orientation by covariance matrix of normal vectors
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);
  
  Point low = bbox_.low();
  Point high = bbox_.high();
  //double len = low.dist(high);
  double rad;
  Point pnt;
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pnt, axis, Cy));
#ifdef DEBUG_CYL
  std::ofstream of("cylinder_compute.g2");
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif
  return cyl;
}

//===========================================================================
void RevEngRegion::analyseCylProject(shared_ptr<Cylinder> cyl, double tol,
				     vector<vector<RevEngPoint*> >& configs)
//===========================================================================
{
#ifdef DEBUG_CYL
  std::ofstream of("projected_pts_cyl.g2");
#endif
  
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  Point pnt = cyl->location();
  Point axis, Cx, Cy;
  cyl->getCoordinateAxes(Cx, Cy, axis);
  double rad = cyl->getRadius();
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pnt, projected, maxdp, avdp);
  shared_ptr<Circle> circ(new Circle(rad, pnt, axis, Cx));
  shared_ptr<SplineCurve> spl;
  Point xpos;
  vector<double> param;
  curveApprox(projected, tol, circ, param, spl, xpos);
  
#ifdef DEBUG_CYL
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    of << projected[kr] << std::endl;
   circ->writeStandardHeader(of);
  circ->write(of);
  spl->writeStandardHeader(of);
  spl->write(of);
#endif
  
  double maxdc, avdc, maxds, avds;
  int num_inc, num_ins;
  RevEngUtils::distToCurve(projected, circ, tol, maxdc, avdc, num_inc);
  RevEngUtils::distToCurve(projected, spl, tol, maxds, avds, num_ins);

  double len = bbox_.low().dist(bbox_.high());
  if (maxds < maxdc && num_ins > num_inc && 
      num_ins > (int)group_points_.size()/4 && rad < 2.0*len)
    {
      // Investigate point configuration
      configSplit(group_points_, param, cyl, spl, maxds, configs);
#ifdef DEBUG_CYL
      std::ofstream ofconf("conf_groups.g2");
      for (size_t ki=0; ki<configs.size(); ++ki)
	{
	  ofconf << "400 1 0 0" << std::endl;
	  ofconf << configs[ki].size() << std::endl;
	  for (size_t kr=0; kr<configs[ki].size(); ++kr)
	    ofconf << configs[ki][kr]->getPoint() << std::endl;
	}
#endif
    }
  
  if (configs.size() <= 1 && (num_ins >= num_inc && avds <= avdc) &&
      (num_ins > (int)group_points_.size()/2 && avds < tol) &&
      (!sweep_.get() || sweep_->type_ != 1 ||
       (sweep_->num_in_ < num_ins && sweep_->avdist_ > avds)))
    {
      // Possible linear sweep
      Point pt1 = pnt - len*axis; 
      Point pt2 = pnt + len*axis; 
      sweep_ = shared_ptr<SweepData>(new SweepData(1, spl, pt1, pt2, maxds, avds, num_ins));
    }
  int stop_break = 1;
 }


//===========================================================================
bool RevEngRegion::defineConeFromCyl(shared_ptr<Cylinder> cyl, double tol,
				     double angtol, int min_pt_reg,
				     double avdist, int num_in, int num2_in,
				     vector<shared_ptr<HedgeSurface> >& hedgesfs,
				     vector<vector<RevEngPoint*> >& out_groups,
				     vector<RevEngPoint*>& single_pts)
//===========================================================================
{
#ifdef DEBUG_CYL
  std::ofstream of("rotated_pts.g2");
#endif
  
  Point pnt = cyl->location();
  Point axis, Cx, Cy;
  cyl->getCoordinateAxes(Cx, Cy, axis);
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  double len = bbox_.low().dist(bbox_.high());
  RevEngUtils::rotateToPlane(group, Cy, axis, pnt, rotated);
  shared_ptr<Line> line(new Line(pnt, axis));
  line->setParameterInterval(-len, len);
  Point pt1 = line->ParamCurve::point(-len);
  Point pt2 = line->ParamCurve::point(len);
  shared_ptr<SplineCurve> line_cv(new SplineCurve(pt1, -len, pt2, len));
  shared_ptr<SplineCurve> crv;
  RevEngUtils::curveApprox(rotated, line_cv, 2, 2, crv);
  double maxdcrv, avdcrv;
  int num_in_crv;
  vector<double> dist;
  RevEngUtils::distToCurve(rotated, crv, tol, maxdcrv, avdcrv, num_in_crv, dist);
#ifdef DEBUG_CYL
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;

  crv->writeStandardHeader(of);
  crv->write(of);
#endif

  if (avdcrv > avdist || num_in_crv < num2_in)
    return false; // No gain in replacing the cylinder with a cone

  // Identify distant points
  double dlim = std::max(2.0*tol, 2.0*avdcrv);
  vector<RevEngPoint*> dist_pts;
  vector<RevEngPoint*> remaining;
  vector<Point> rotated2;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (dist[ki] > dlim)
	dist_pts.push_back(group_points_[ki]);
      else
	{
	  remaining.push_back(group_points_[ki]);
	  rotated2.push_back(rotated[ki]);
	}
    }
#ifdef DEBUG_CYL
  std::ofstream of3("split_group.g2");
  if (remaining.size() > 0)
    {
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << remaining.size() << std::endl;
      for (size_t kr=0; kr<remaining.size(); ++kr)
	of3 << remaining[kr]->getPoint() << std::endl;
    }
  if (dist_pts.size() > 0)
    {
      of3 << "400 1 0 4 0 255 0 255" << std::endl;
      of3 << dist_pts.size() << std::endl;
      for (size_t kr=0; kr<dist_pts.size(); ++kr)
	of3 << dist_pts[kr]->getPoint() << std::endl;
    }
#endif
  if (remaining.size() < 0.5*group_points_.size())
    return false;
  
  // Define cone
  shared_ptr<SplineCurve> crv2;
  RevEngUtils::curveApprox(rotated2, line_cv, 2, 2, crv2);
  double tclose, dclose;
  Point ptclose;
  crv2->closestPoint(pnt, crv2->startparam(), crv2->endparam(), tclose,
		     ptclose, dclose);
  vector<Point> der(2);
  crv2->point(der, tclose, 1);
  double phi = der[1].angle(axis);

  // Check sign of angle
  Point pnt2 = pnt + 0.1*len*axis;
  double tclose2, dclose2;
  Point ptclose2;
  crv2->closestPoint(pnt2, crv2->startparam(), crv2->endparam(), tclose2,
		     ptclose2, dclose2);
  if (dclose2 < dclose)
    phi *= -1;
      
  shared_ptr<Cone> cone(new Cone(dclose, pnt, axis, Cy, phi));
#ifdef DEBUG_CYL
  std::ofstream of2("cone_from_rotate.g2");
  double t1 = crv2->startparam();
  double t2 = crv2->endparam();
  cone->setParamBoundsV(t1-0.1*(t2-t1), t2+0.1*(t2-t1));
  cone->writeStandardHeader(of2);
  cone->write(of2);
#endif

  // Check accuracy
  bool found = false;
  double maxdistco, avdistco;
  int num_insideco, num2_insideco;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inco, outco;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  cone, tol, maxdistco, avdistco,
			  num_insideco, num2_insideco,
			  inco, outco, parvals, dist_ang, angtol);
  int sf_flag = defineSfFlag(remaining.size(), 0, tol, num_insideco, 
			     num2_insideco, avdistco, true);
  if (sf_flag < ACCURACY_POOR)
    {
      bool OK = true;
      double acc_fac = 1.5;
      double frac = (double)remaining.size()/(double)group_points_.size();
      if (associated_sf_.size() > 0)
	{
	  double maxd_init, avd_init;
	  int num_init, num_init2;
	  getAccuracy(maxd_init, avd_init, num_init, num_init2);
	  if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	    OK = false;
	  else
	    {
	      if (!((double)num_insideco < frac*(double)num_init ||
		    (avdistco < avd_init &&
		     (double)num_insideco < acc_fac*frac*(double)num_init)))
		OK = true;
	    }
	  if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
	    OK = true;
	}

      if (OK)
	{
	  found = true;
	  std::swap(group_points_, remaining);
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxdistco, avdistco, num_insideco, num2_insideco);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cone, this));

	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	  
	  // Make sure that the distant points are connected
	  vector<vector<RevEngPoint*> > separate_groups;
	  vector<RevEngPoint*> dummy;
	  connectedGroups(dist_pts, separate_groups, false, dummy);
	  if (separate_groups.size() > 0)
	    {
	      for (size_t kh=0; kh<separate_groups.size(); ++kh)
		{
#ifdef DEBUG_CYL
		  std::ofstream of4("curr_sep_group.g2");
		  of4 << "400 1 0 4 0 0 255 255" << std::endl;
		  of4 << separate_groups[kh].size() << std::endl;
		  for (size_t kr=0; kr<separate_groups[kh].size(); ++kr)
		    of4 << separate_groups[kh][kr]->getPoint() << std::endl;
#endif
		  // Check if the group can be added to the current region given
		  // the updated cone
		  vector<RevEngPoint*> in2, out2;
		  vector<pair<double,double> > dist_ang2;
		  vector<double> parvals2;
		  int nmb_in2, nmb2_in2;
		  double maxd2, avd2;
		  RevEngUtils::distToSurf(separate_groups[kh].begin(),
					  separate_groups[kh].end(),
					  cone, tol, maxd2, avd2, nmb_in2, nmb2_in2,
					  in2, out2, parvals2, dist_ang2, angtol);
		  int sf_flag2 = defineSfFlag(separate_groups[kh].size(), 0, tol,
					      nmb_in2, nmb2_in2, avd2, true);

		  if (sf_flag2 <= sf_flag)
		    {
		      for (size_t kr=0; kr<separate_groups[kh].size(); ++kr)
			{
			  separate_groups[kh][kr]->setPar(Vector2D(parvals2[2*kr],
								   parvals2[2*kr+1]));
			  separate_groups[kh][kr]->setSurfaceDist(dist_ang2[kr].first,
								  dist_ang2[kr].second);
			}
		      (void)addPointsToGroup(separate_groups[kh], tol, angtol);
		    }
		  else if (separate_groups[kh].size() == 1)
		    {
		      separate_groups[kh][0]->unsetRegion();
		      single_pts.push_back(separate_groups[kh][0]);
		    }
		  else
		    {
		      for (size_t kr=0; kr<separate_groups[kh].size(); ++kr)
			separate_groups[kh][kr]->unsetRegion();
		      out_groups.push_back(separate_groups[kh]);
		    }
		}
	    }
	}
    }
  return found;
 }
  
//===========================================================================
void RevEngRegion::analyseCylRotate(shared_ptr<Cylinder> cyl, double tol,
				    double avdist, int num_in,
				    double& avdist_lin, int& num_in_lin,
				    double& avdist_cub, int& num_in_cub,
				    shared_ptr<Cone>& cone)
//===========================================================================
{
#ifdef DEBUG_CYL
  std::ofstream of("rotated_pts_cyl.g2");
#endif
  
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  Point pnt = cyl->location();
  Point axis, Cx, Cy;
  cyl->getCoordinateAxes(Cx, Cy, axis);
  //double rad = cyl->getRadius();
  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cy, axis, pnt, rotated);
  double len = bbox_.low().dist(bbox_.high());
#ifdef DEBUG_CYL
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << pnt-0.5*len*axis << " " << pnt+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << pnt-0.5*len*Cy << " " << pnt+0.5*len*Cy << std::endl;
#endif

  shared_ptr<Line> line(new Line(pnt, axis));
  line->setParameterInterval(-len, len);
  shared_ptr<SplineCurve> line2;
  RevEngUtils::curveApprox(rotated, line, 2, 2, line2);
  vector<Point> der(2);
  line2->point(der, 0.0, 1);

  Point pt1 = line->ParamCurve::point(-len);
  Point pt2 = line->ParamCurve::point(len);
  shared_ptr<SplineCurve> line_cv(new SplineCurve(pt1, -len, pt2, len));
  shared_ptr<SplineCurve> cv1, cv2;
  RevEngUtils::curveApprox(rotated, line_cv, 2, 2, cv1);
  RevEngUtils::curveApprox(rotated, line_cv, 4, 8, cv2);
  double maxdcv1, maxdcv2;
  RevEngUtils::distToCurve(rotated, cv1, tol, maxdcv1, avdist_lin, num_in_lin);
  RevEngUtils::distToCurve(rotated, cv2, tol, maxdcv2, avdist_cub, num_in_cub);
#ifdef DEBUG_CYL
  cv1->writeStandardHeader(of);
  cv1->write(of);
  cv2->writeStandardHeader(of);
  cv2->write(of);
#endif
  
  // Check if a cone is feasible
  double cfac = 1.2;
  if (avdist_lin < cfac*std::min(avdist_cub, avdist) &&
      cfac*num_in_lin > (double)std::max(num_in_cub, num_in))
    {
      double tclose, dclose;
      Point ptclose;
      cv1->closestPoint(pnt, cv1->startparam(), cv1->endparam(), tclose,
			ptclose, dclose);
      vector<Point> der(2);
      cv1->point(der, tclose, 1);
      double phi = der[1].angle(axis);

      // Check sign of angle
      Point pnt2 = pnt + 0.1*len*axis;
      double tclose2, dclose2;
      Point ptclose2;
      cv1->closestPoint(pnt2, cv1->startparam(), cv1->endparam(), tclose2,
			ptclose2, dclose2);
      if (dclose2 < dclose)
	phi *= -1;
      
      cone = shared_ptr<Cone>(new Cone(dclose, pnt, axis, Cy, phi));
#ifdef DEBUG_CYL
      std::ofstream of2("cone_from_rotate.g2");
      double t1 = cv1->startparam();
      double t2 = cv1->endparam();
      cone->setParamBoundsV(t1-0.1*(t2-t1), t2+0.1*(t2-t1));
      cone->writeStandardHeader(of2);
      cone->write(of2);
#endif
    }

  int stop_break = 1;
}



//===========================================================================
bool RevEngRegion::extractLinearSweep(double tol, int min_pt, int min_pt_reg,
				      double angtol, int prefer_elementary,
				      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      std::vector<HedgeSurface*>& prevsfs)
//===========================================================================
{
  if (!sweep_.get())
    return false;

  if (sweep_->type_ != 1)
    return false;
  
  bool found = false;
  shared_ptr<SplineSurface> surf;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);
  
  SweepSurfaceCreator sweep;
  shared_ptr<SplineCurve> along(new SplineCurve(sweep_->location_, sweep_->added_info_));
  Point mid = 0.5*(sweep_->location_ + sweep_->added_info_);
  surf = shared_ptr<SplineSurface>(sweep.linearSweptSurface(*along,
							  *sweep_->profile_, mid));
  if (!surf.get())
    return false;
#ifdef DEBUG_EXTRACT
  std::ofstream of3("sweep.g2");
  along->writeStandardHeader(of3);
  along->write(of3);
  sweep_->profile_->writeStandardHeader(of3);
  sweep_->profile_->write(of3);
  surf->writeStandardHeader(of3);
  surf->write(of3);
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << sweep_->location_ << std::endl;
#endif
  
  // Check accuracy
  double maxd, avd; 
  int num2, num2_2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  surf, tol, maxd, avd, num2, num2_2, inpt, outpt,
			  parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_surf.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif
  //int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  int sf_flag = defineSfFlag(0, tol, num2, num2_2, avd, false);
  if (sf_flag < ACCURACY_POOR)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2, num2_2);
#ifdef DEBUG
	  std::cout << "Linear swept surface. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
#endif
	  
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf, this));
	  hedge->setLinearSweepInfo(sweep_->profile_, sweep_->location_, sweep_->added_info_);
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	}
    }

   return found;
 }

//===========================================================================
shared_ptr<SplineSurface>
RevEngRegion::computeLinearSwept(double tol, shared_ptr<SplineCurve>& profile,
				 Point& pt1, Point& pt2)
//===========================================================================
{
  // Axis by covariance matrix of normal vectors
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::computeAxis(group, axis, Cx, Cy);

  // Point on axis
  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  double rad;
  Point pnt;
  RevEngUtils::computeCylPosRadius(group, low, high,
				   axis, Cx, Cy, pnt, rad);

  // Project the points onto the defined plane 
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pnt, projected, maxdp, avdp);

  // Approximate the projected point cloud with a spline curve
  shared_ptr<Circle> circ(new Circle(rad, pnt, axis, Cx));
  Point xpos;
  vector<double> param;
  curveApprox(projected, tol, circ, param, profile, xpos);

  pt1 = pnt - len*axis;
  pt2 = pnt + len*axis;
  Point mid = 0.5*(pt1 + pt2);
  
  SweepSurfaceCreator sweep;
  shared_ptr<SplineCurve> along(new SplineCurve(pt1, pt2));
  shared_ptr<SplineSurface> swept_surf(sweep.linearSweptSurface(*along, *profile, mid));
				       
  return swept_surf;
}

//===========================================================================
bool RevEngRegion::extractSphere(Point mainaxis[3],
				 double tol, int min_pt, int min_pt_reg,
				 double angtol, int prefer_elementary,
				 vector<shared_ptr<HedgeSurface> >& hedgesfs,
				 vector<HedgeSurface*>& prevsfs,
				 vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  Point adj_axis;
  int num_pt_adj = 0;
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      int type = adj_elem[ki].first->instanceType();
      if (type == Class_Cylinder || type == Class_Cone)
	{
	  if (adj_elem[ki].second->numPoints() > num_pt_adj)
	    {
	      adj_axis = adj_elem[ki].first->direction();
	      num_pt_adj = adj_elem[ki].second->numPoints();
	    }
	}
    }
  
  shared_ptr<Sphere> sphere = computeSphere(mainaxis, adj_axis, group_points_);
  if (!sphere.get())
    return false;
#ifdef DEBUG_EXTRACT
  std::ofstream ofs("sph.g2");
  sphere->writeStandardHeader(ofs);
  sphere->write(ofs);
#endif  
  // Check accuracy
  double maxd, avd; 
  int num2, num2_2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  sphere, tol, maxd, avd, num2, num2_2, inpt, outpt,
			  parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_sphere.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif
  
  if (inpt.size() > group_points_.size()/2 && (int)inpt.size() > min_nmb)
    {
      shared_ptr<Sphere> sphere_in = computeSphere(mainaxis, adj_axis, inpt);
      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
      vector<pair<double, double> > dist_ang_in;
      vector<double> parvals_in;
      double maxd_in, avd_in;
      int num2_in, num2_in2;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      sphere_in, tol, maxd_in, avd_in, num2_in, num2_in2,
			      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
#ifdef DEBUG_EXTRACT
      std::ofstream ofd2("in_out_sphere2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in.size() << std::endl;
      for (size_t kr=0; kr<inpt_in.size(); ++kr)
	ofd2 << inpt_in[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in.size() << std::endl;
      for (size_t kr=0; kr<outpt_in.size(); ++kr)
	ofd2 << outpt_in[kr]->getPoint() << std::endl;
#endif
      
      if (num2_in > num2 || avd_in < avd)
	{
	  std::swap(sphere, sphere_in);
	  std::swap(num2, num2_in);
	  std::swap(num2_2, num2_in2);
	  std::swap(avd, avd_in);
	  std::swap(maxd, maxd_in);
	  std::swap(parvals, parvals_in);
	  std::swap(dist_ang, dist_ang_in);
	  // std::cout << "Sphere swap" << std::endl;
	  // std::cout << group_points_.size() << " " << num2_in;
	  // std::cout << " " << num2 << " " << avd_in << " " << avd;
	  // std::cout << " " << maxd_in << " " << maxd << std::endl;
	}
    }

  //int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  int sf_flag = defineSfFlag(0, tol, num2, num2_2,
			     avd, false);
  if (sf_flag < ACCURACY_POOR)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	    OK = false;
	  else if (surfflag_ < ACCURACY_POOR && sf_flag >= ACCURACY_POOR)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	  if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
	    OK = true;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2, num2_2);
#ifdef DEBUG
	  std::cout << "Sphere. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
#endif
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(sphere, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	}
    }
  if (!basesf_.get() ||
      (num2 >= num_in_base_ && avd < avdist_base_))
    setBaseSf(sphere, maxd, avd, num2, num2_2);

  return found;
}

//===========================================================================
shared_ptr<Sphere> RevEngRegion::computeSphere(Point mainaxis[3],
					       Point adj_axis,
					       vector<RevEngPoint*>& points)
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  Point centre;
  double radius;
  try {
    RevEngUtils::computeSphereProp(group, centre, radius);
  }
  catch (...)
    {
      shared_ptr<Sphere> dummy;
      return dummy;
    }

  Point z_axis;
  int ix = -1;
  double minang = M_PI;
  if (adj_axis.dimension() > 0)
    z_axis = adj_axis;
  else
    {
      Point norm = normalcone_.centre();
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = mainaxis[ka].angle(norm);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < minang)
	    {
	      minang = ang;
	      ix = ka;
	    }
	}
      z_axis = mainaxis[ix];
    }
  
  // Define x-axis
  ix = -1;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(z_axis);
      ang = std::min(ang, M_PI-ang);
      if (ang < minang)
	{
	  minang = ang;
	  ix = ka;
	}
    }

  Point Cy = mainaxis[(ix+1)%3].cross(z_axis);
  Point Cx = z_axis.cross(Cy);
  
  shared_ptr<Sphere> sph(new Sphere(radius, centre, z_axis, Cx));
  //cyl->setParamBoundsV(-len, len);
  return sph;
}

//===========================================================================
bool RevEngRegion::extractCone(double tol, int min_pt, int min_pt_reg,
			       double angtol, int prefer_elementary,
			       vector<shared_ptr<HedgeSurface> >& hedgesfs,
			       vector<HedgeSurface*>& prevsfs,
			       vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
 bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  double minang = 0.05;
  if ((int)group_points_.size() < min_nmb)
    return false;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);
  
  Point apex;
  shared_ptr<Cone> cone = computeCone(group_points_, apex);
#ifdef DEBUG_EXTRACT
  std::ofstream ofs("conesf.g2");
  cone->writeStandardHeader(ofs);
  cone->write(ofs);
#endif
  
  // Check accuracy
  double maxd, avd; 
  int num2, num2_2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  cone, tol, maxd, avd, num2, num2_2, inpt, outpt,
			  parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_cone.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif
  
  if (inpt.size() > group_points_.size()/2 && (int)inpt.size() > min_nmb)
    {
      Point apex2;
      shared_ptr<Cone> cone_in = computeCone(inpt, apex2);
      vector<RevEngPoint*> inpt_in, outpt_in; //, inpt2, outpt2;
      vector<pair<double, double> > dist_ang_in;
      vector<double> parvals_in;
      double maxd_in, avd_in;
      int num2_in, num2_in2;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      cone_in, tol, maxd_in, avd_in, num2_in, num2_in2, 
			      inpt_in, outpt_in, parvals_in, dist_ang_in, angtol);
  
#ifdef DEBUG_EXTRACT
      std::ofstream ofd2("in_out_cone2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in.size() << std::endl;
      for (size_t kr=0; kr<inpt_in.size(); ++kr)
	ofd2 << inpt_in[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in.size() << std::endl;
      for (size_t kr=0; kr<outpt_in.size(); ++kr)
	ofd2 << outpt_in[kr]->getPoint() << std::endl;
#endif
      
      if (num2_in > num2 || avd_in < avd)
	{
	  std::swap(cone, cone_in);
	  std::swap(num2, num2_in);
	  std::swap(num2_2, num2_in2);
	  std::swap(avd, avd_in);
	  std::swap(maxd, maxd_in);
	  std::swap(parvals, parvals_in);
	  std::swap(dist_ang, dist_ang_in);
	  std::swap(apex, apex2);
	  // std::cout << "Cone swap" << std::endl;
	  // std::cout << group_points_.size() << " " << num2_in;
	  // std::cout << " " << num2 << " " << avd_in << " " << avd;
	  // std::cout << " " << maxd_in << " " << maxd << std::endl;
	}
    }

  //int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  int sf_flag = defineSfFlag(0, tol, num2, num2_2,
			     avd, true);
  double cone_ang = cone->getConeAngle();
  if (0.5*M_PI-fabs(cone_ang) > minang && (!bbox_.containsPoint(apex, tol)) &&
      sf_flag < ACCURACY_POOR)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	  if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
	    OK = true;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2, num2_2);
#ifdef DEBUG
	  std::cout << "Cone. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
#endif
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(cone, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	}
    }
  if (!basesf_.get() ||
      (num2 >= num_in_base_ && avd < avdist_base_))
    setBaseSf(cone, maxd, avd, num2, num2_2);

  return found;
}


//===========================================================================
shared_ptr<Cone> RevEngRegion::computeCone(vector<RevEngPoint*>& points, Point& apex)
//===========================================================================
{
  Point low = bbox_.low();
  Point high = bbox_.high();

  Point mid(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)points.size();
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      Vector3D xyz = points[kr]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      mid += wgt*pnt;
    }
  Point axis;
  Point Cx;
  Point Cy;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::coneAxis(group, axis, Cx, Cy);

  //Point apex;
  double phi;
  RevEngUtils::coneApex(group, axis, apex, phi);
  
  vector<Point> rotated;
  RevEngUtils::rotateToPlane(group, Cx, axis, apex, rotated);
#ifdef DEBUG_EXTRACT
  double len = low.dist(high);
  std::ofstream of("rotated_pts_cone.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << apex-0.5*len*axis << " " << apex+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << apex-0.5*len*Cx << " " << apex+0.5*len*Cx << std::endl;
#endif
  
  Point pos0 = apex + ((mid - apex)*axis)*axis;
  double dd = apex.dist(pos0);
  double rad0 = dd*tan(phi);
  shared_ptr<Cone> cone(new Cone(fabs(rad0), pos0, axis, Cx, phi));
  //cone->setParamBoundsV(-0.2*len, 0.2*len);
  
  shared_ptr<Circle> circ(new Circle(rad0, pos0, axis, Cx));
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group, axis, pos0, projected, maxdp, avdp);
#ifdef DEBUG_EXTRACT
  std::ofstream ofp3("projected_pts_cone.g2");
  ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
  ofp3 << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp3 << projected[kr] << std::endl;
  circ->writeStandardHeader(ofp3);
  circ->write(ofp3);
#endif
  
  return cone;
}

//===========================================================================
bool RevEngRegion::adjacentToCylinder(Point mainaxis[3],
				      double tol, int min_pt, int min_pt_reg,
				      double angtol, int prefer_elementary,
				      vector<shared_ptr<HedgeSurface> >& hedgesfs,
				      vector<HedgeSurface*>& prevsfs,
				      vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  // Check for one and only one cylinder axis
  Point pos, axis;
  size_t cyl_ix;
  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      if (adj_elem[ki].first->instanceType() == Class_Cylinder)
	{
	  Point loc = adj_elem[ki].first->location();
	  Point dir = adj_elem[ki].first->direction();
	  if (pos.dimension() == 0)
	    {
	      pos = loc;
	      axis = dir;
	      cyl_ix = ki;
	    }
	  else
	    {
	      // Check compatibility
	      double ang = axis.angle(dir);
	      ang = std::min(ang, M_PI-ang);
	      if (ang > angtol)
		return false;

	      Point axis_pt = pos + ((loc - pos)*axis)*axis;
	      double dist = loc.dist(axis_pt);
	      if (dist > tol)
		return false;
	    }
	}
    }
  if (pos.dimension() == 0)
    return false;  // No adjacent cylinder

  double min_ang = M_PI;
  int ix = -1;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(axis);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  ix = ka;
	}
    }

  if (ix < 0)
    return false; // Should not happen
  
  Point Cx = mainaxis[(ix+1)%3];
  Point Cy = axis.cross(Cx);
  Cy.normalize();
  Cx = Cy.cross(axis);
  Cx.normalize();
  
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, Cy, axis, pos, rotated);
  double len = bbox_.low().dist(bbox_.high());

#ifdef DEBUG_CYL
  std::ofstream of("rotated_pts.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*Cy << " " << pos+0.5*len*Cy << std::endl;
#endif
      
  // Approximate rotated points with a circle
  Point centre;
  double radius;
  RevEngUtils::computeCircPosRadius(rotated, Cx, Cy, axis, centre, radius);
  shared_ptr<Circle> circ(new Circle(radius, centre, Cx, Cy));
#ifdef DEBUG_CYL
  circ->writeStandardHeader(of);
  circ->write(of);
#endif

  // Check accuracy
  bool found = false;
  double maxd, avd;
  int num_in;
  RevEngUtils::distToCurve(rotated, circ, tol, maxd, avd, num_in);
  if (accuracyOK(min_pt, tol, num_in, avd))
    {
      // A surface can be fitted. Check if it is a torus or a sphere
      Point axis_pt = pos + ((centre - pos)*axis)*axis;
      double dist = centre.dist(axis_pt);
      shared_ptr<ElementarySurface> surf;
      if (dist < tol)
	// Sphere
	surf = shared_ptr<ElementarySurface>(new Sphere(radius, centre,
							axis, Cx));
      else
	// Torus
	surf = shared_ptr<ElementarySurface>(new Torus(dist, radius,
						       axis_pt, axis, Cx));

#ifdef DEBUG_CYL
      std::ofstream of2("rot_sf.g2");
      surf->writeStandardHeader(of2);
      surf->write(of2);
      adj_elem[cyl_ix].first->writeStandardHeader(of2);
      adj_elem[cyl_ix].first->write(of2);
#endif
      // Check accuracy of surface
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      double maxd, avd;
      int num_in, num2_in;
      vector<double> parvals;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      surf, tol, maxd, avd, num_in, num2_in,
			      inpt, outpt, parvals, dist_ang,
			      angtol);
#ifdef DEBUG_CYL
      std::ofstream ofd2("in_out_elem.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt.size() << std::endl;
      for (size_t kr=0; kr<inpt.size(); ++kr)
	ofd2 << inpt[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt.size() << std::endl;
      for (size_t kr=0; kr<outpt.size(); ++kr)
	ofd2 << outpt[kr]->getPoint() << std::endl;
#endif

      int sf_flag = defineSfFlag(0, tol, num_in, num2_in,
				 avd, false);
      if (sf_flag < ACCURACY_POOR)
	{
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num_in, num2_in);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(surf, this));
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	  found = true;
	}
      
      if (!basesf_.get() ||
	  (num2_in >= num_in_base_ && avd < avdist_base_))
	setBaseSf(surf, maxd, avd, num_in, num2_in);
    }
  return found;
}



//===========================================================================
bool RevEngRegion::contextCylinder(Point mainaxis[3],
				   double tol, int min_pt, int min_pt_reg,
				   double angtol, int prefer_elementary,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  return contextCylinder(mainaxis, tol, min_pt, min_pt_reg, angtol, 
			 prefer_elementary, adj_elem, hedgesfs, prevsfs,
			 out_groups);
}


//===========================================================================
bool RevEngRegion::contextCylinder(Point mainaxis[3],
				   double tol, int min_pt,  int min_pt_reg,
				   double angtol, int prefer_elementary,
				   vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups,
				   int mode)
//===========================================================================
{
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  double angtol2 = 0.1;
  Point pos, axis, Cx;
  double rad;
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_planar;
  vector<RevEngPoint*> cyl_pts;
  bool foundaxis = analyseCylinderContext(adj_elem, tol, angtol2, mainaxis,
					  mode, adj_planar, pos, axis, Cx,
					  rad, cyl_pts);
  if (!foundaxis)
    return false;

  shared_ptr<Cylinder> cyl(new Cylinder(rad, pos, axis, Cx));

#ifdef DEBUG_CYLCONTEXT
  std::ofstream of("context_cyl.g2");
  cyl->writeStandardHeader(of);
  cyl->write(of);
#endif

  // // Fetch remaining points
  // vector<RevEngPoint*> remaining;
  // getRemainingPoints(cyl_pts, remaining);
  
  // Move points as appropriate
  // First check accuracy with respect to cylinder
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  RevEngUtils::distToSurf(cyl_pts.begin(), cyl_pts.end(),
			  cyl, tol, maxd, avd, num_in, num2_in,
			  inpt, outpt, parvals, dist_ang,
			  angtol);
  int sf_flag_cyl = defineSfFlag((int)cyl_pts.size(), 0, tol, num_in, num2_in,
				 avd, true);
  if (sf_flag_cyl >= ACCURACY_POOR)
    return false;  // Not a good subset

  // Check if parts of this regions represents a likely cylinder
  if (hasSurface())
    {
      double fac = 1.5;
      if (surfflag_ < sf_flag_cyl || avdist_ < fac*avd ||
  	  (avdist_ <= tol && avd > tol))
  	return false;  // Not an improvement
    }

  double cyl_lim = 0.3;
  if (num2_in < (int)(cyl_lim*(double)group_points_.size()))
      return false;

  // Extract points from cylinder
  vector<vector<RevEngPoint*> > planar(adj_planar.size());
  vector<RevEngPoint*> remaining;
  vector<Point> norm(adj_planar.size());
  vector<Point> loc(adj_planar.size());
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      norm[ki] = adj_planar[ki].first->direction();
      Point avnorm = adj_planar[ki].second->getMeanNormal();
      if (norm[ki]*avnorm < 0.0)
	norm[ki] *= -1;
      loc[ki] = adj_planar[ki].first->location();
    }
  
  double cyl_dom[4];   // Cylinder domain
  cyl_dom[0] = cyl_dom[2] = std::numeric_limits<double>::max();
  cyl_dom[1] = cyl_dom[3] = std::numeric_limits<double>::lowest();
  double ang_fac = 1.1;
  double pi4 = 0.25*M_PI;
  vector<double> ptdist(adj_planar.size());
  vector<double> ptang(adj_planar.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point ptnorm = group_points_[ki]->getLocFuncNormal();
      size_t kj;
      for (kj=0; kj<ptdist.size(); ++kj)
	{
	  ptdist[kj] = fabs((ptpos-loc[kj])*norm[kj]);
	  ptang[kj] = norm[kj].angle(ptnorm);
	  ptang[kj] = std::min(ptang[kj], M_PI-ptang[kj]);
	}
      double dist3 = dist_ang[ki].first;
      Point axis_pt = pos + ((ptpos - pos)*axis)*axis;
      double dist4 = ptpos.dist(axis_pt);
      double ang3 = dist_ang[ki].second;
      for (kj=0; kj<ptdist.size(); ++kj)
	{
	  if (ptdist[kj] < dist3 && ptang[kj] < ang_fac*ang3 &&
	      (dist4 >= rad || axis.angle(norm[kj]) > pi4))
	    {
	      planar[kj].push_back(group_points_[ki]);
	      break;
	    }
	}
      if (kj == ptdist.size())
	{
	  remaining.push_back(group_points_[ki]);
	  cyl_dom[0] = std::min(cyl_dom[0], parvals[2*ki]);
	  cyl_dom[1] = std::max(cyl_dom[1], parvals[2*ki]);
	  cyl_dom[2] = std::min(cyl_dom[2], parvals[2*ki+1]);
	  cyl_dom[3] = std::max(cyl_dom[3], parvals[2*ki+1]);
	}
    }
  
  // Planes to cylinder
  double eps = 1.0e-2;
  vector<vector<RevEngPoint*> > plane2cyl(adj_planar.size());
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      for (auto it=adj_planar[ki].second->pointsBegin();
	   it!=adj_planar[ki].second->pointsEnd(); ++it)
	{
	  Vector3D xyz = (*it)->getPoint();
	  Point ptpos(xyz[0], xyz[1], xyz[2]);
	  Point ptnorm = (*it)->getLocFuncNormal();
	  double dist1 = (*it)->getSurfaceDist();
	  double ang1 = norm[ki].angle(ptnorm);
	  ang1 = std::min(ang1, M_PI-ang1);
	  Point vec = ptpos - pos;
	  double dist2 = (vec - (vec*axis)*vec).length();
	  dist2 = fabs(dist2 - rad);
	  //double dist2 = fabs((ptpos - pos)*axis);
	  if (dist2 < dist1)
	    {
	      double upar, vpar, tdist;
	      Point close;
	      cyl->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	      Point cnorm;
	      cyl->normal(cnorm, upar, vpar);
	      double ang2 = cnorm.angle(ptnorm);
	      if (ang2 <= angtol2 && vpar >= cyl_dom[2] && vpar <= cyl_dom[3])
		plane2cyl[ki].push_back(*it);
	    }
	}
    }
  
#ifdef DEBUG_CYLCONTEXT
  std::ofstream of1("cyl2plane.g2");
  for (size_t ki=0; ki<planar.size(); ++ki)
    if (planar[ki].size() > 0)
      {
	of1 << "400 1 0 0" << std::endl;
	of1 << planar[ki].size() << std::endl;
	for (size_t kj=0; kj<planar[ki].size(); ++kj)
	  of1 << planar[ki][kj]->getPoint() << std::endl;
      }

  std::ofstream of2("plane2cyl.g2");
  for (size_t ki=0; ki<plane2cyl.size(); ++ki)
    if (plane2cyl[ki].size() > 0)
      {
	of2 << "400 1 0 0" << std::endl;
	of2 << plane2cyl[ki].size() << std::endl;
	for (size_t kj=0; kj<plane2cyl[ki].size(); ++kj)
	  of2 << plane2cyl[ki][kj]->getPoint() << std::endl;
      }
#endif
  
  // Recompute accuracy
  vector<RevEngPoint*> inpt2, outpt2;
  vector<pair<double, double> > dist_ang2;
  double maxd2, avd2;
  int num_in2, num2_in2;
  vector<double> parvals2;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  cyl, tol, maxd2, avd2, num_in2, num2_in2,
			  inpt2, outpt2, parvals2, dist_ang2,
			  angtol);

  if (num2_in2 > min_pt && num2_in2 > (int)remaining.size()/2)
    {
      // Check for deviant points at the boundary
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						    edge_class_type_, remaining));
      
      vector<RevEngPoint*> dist_points;
      reg->identifyDistPoints(dist_ang2, tol, maxd2, avd2, dist_points);
      if (dist_points.size() > 0)
	{
	  reg->extractSpesPoints(dist_points, out_groups, true);

	  vector<RevEngPoint*> remaining2;
	  remaining2 = reg->getPoints();
	  for (size_t ki=0; ki<remaining2.size(); ++ki)
	    remaining2[ki]->setRegion(this);

	  std::swap(remaining, remaining2);
	}
      else
	{
	  for (size_t ki=0; ki<remaining.size(); ++ki)
	    remaining[ki]->setRegion(this);
	}
    }

  for (size_t ki=0; ki<plane2cyl.size(); ++ki)
    if (plane2cyl[ki].size() > 0)
      {
	for (size_t kj=0; kj<plane2cyl[ki].size(); ++kj)
	  plane2cyl[ki][kj]->setRegion(this);
	remaining.insert(remaining.end(), plane2cyl[ki].begin(), plane2cyl[ki].end()); // Is now double
      }
  
  // Recompute accuracy
  Point x_axis, y_axis, z_axis;
  cyl->getCoordinateAxes(x_axis, y_axis, z_axis);
  Point pos2;
  double rad2;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  points.push_back(make_pair(remaining.begin(), remaining.end()));
  Point low = bbox_.low();
  Point high = bbox_.high();
  RevEngUtils::computeCylPosRadius(points, low, high, z_axis,
				   x_axis, y_axis, pos2, rad2);
  shared_ptr<Cylinder> cyl2(new Cylinder(rad2, pos2, z_axis, x_axis));
#ifdef DEBUG_CYLCONTEXT
  std::ofstream of3("context_cyl2.g2");
  cyl2->writeStandardHeader(of3);
  cyl2->write(of3);
#endif
  
  vector<RevEngPoint*> inpt3, outpt3;
  vector<pair<double, double> > dist_ang3;
  double maxd3, avd3;
  int num_in3, num2_in3;
  vector<double> parvals3;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  cyl, tol, maxd3, avd3, num_in3, num2_in3,
			  inpt3, outpt3, parvals3, dist_ang3,
			  angtol);
  
  vector<RevEngPoint*> inpt4, outpt4;
  vector<pair<double, double> > dist_ang4;
  double maxd4, avd4;
  int num_in4, num2_in4;
  vector<double> parvals4;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  cyl2, tol, maxd4, avd4, num_in4, num2_in4,
			  inpt4, outpt4, parvals4, dist_ang4,
			  angtol);
  if (avd4 < avd3 && num_in4+num2_in4 < num_in3+num2_in3)
    {
      std::swap(maxd4, maxd3);
      std::swap(avd4, avd3);
      std::swap(num_in4, num_in3);
      std::swap(num2_in4, num2_in3);
      std::swap(parvals4, parvals3);
      std::swap(dist_ang4, dist_ang3);
      std::swap(inpt4, inpt3);
      std::swap(outpt4, outpt3);
    }
  
  // Move points
#ifdef DEBUG_CYLCONTEXT
  std::ofstream of4("cyl_pts.g2");
  writeRegionPoints(of4);
#endif
  
  int sf_flag = defineSfFlag((int)remaining.size(), 0, tol, num_in3, num2_in3, avd3, true);
  bool OKsurf = (sf_flag < ACCURACY_POOR);
  bool hasSurf = (associated_sf_.size() > 0);
  if (hasSurf)
    {
      double acc_fac = 1.5;
      double maxd_init, avd_init;
      int num_init, num_init2;
      getAccuracy(maxd_init, avd_init, num_init, num_init2);
      int sfcode;
      int sftype = associated_sf_[0]->instanceType(sfcode);
      double ang = (sftype == Class_Plane) ?
	normalcone_.angle() : M_PI;
      double ang_lim = 0.1*M_PI;
	      
      // Check with current approximating surface
      if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	OKsurf = false;
      else if (num_init > num_in3 && num_init2 > num2_in3 && avd_init < avd3)
	OKsurf = false;
      else if (prefer_elementary == ALWAYS_ELEM ||
	       prefer_elementary == PREFER_ELEM)
	{
	  if (!(ang > ang_lim && (num_in3 < num_init ||
				  (avd3 < avd_init &&
				   num_in3 < acc_fac*num_init))))
	    OKsurf = false;
	}
      else
	{
	  if (!(num_in3 < num_init ||
		(avd3 < avd_init && num_in3 < acc_fac*num_init)))
	    OKsurf = true;
	}
      if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
	OKsurf = true;
    }


  if (OKsurf || hasSurf == false)
    {
      std::swap(group_points_, remaining);
      for (size_t ki=0; ki<plane2cyl.size(); ++ki)
	if (plane2cyl[ki].size() > 0)
	  {
	    for (size_t kj=0; kj<plane2cyl[ki].size(); ++kj)
	      {
		adj_planar[ki].second->removePoint(plane2cyl[ki][kj]);
		plane2cyl[ki][kj]->setRegion(this);
	      }
	    adj_planar[ki].second->updateInfo();
	  }
    }


  if (OKsurf)
    {
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals3[2*ki],parvals3[2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang3[ki].first, dist_ang3[ki].second);
	}
      setAccuracy(maxd3, avd3, num_in3, num2_in3);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }
  
  if (OKsurf || hasSurf == false)
    {
      for (size_t kj=0; kj<planar.size(); ++kj)
	{
	  if (planar[kj].size() > 0)
	    {
	      // Parameterize
	      vector<RevEngPoint*> inpt4, outpt4;
	      vector<pair<double, double> > dist_ang4;
	      double maxd4, avd4;
	      int num_in4, num2_in4;
	      vector<double> parvals4;
	      RevEngUtils::distToSurf(planar[kj].begin(), planar[kj].end(),
				      adj_planar[kj].first, tol, maxd4, avd4,
				      num_in4, num2_in4, inpt4, outpt4, parvals4,
				      dist_ang4, angtol);
	      for (size_t ki=0; ki<planar[kj].size(); ++ki)
		{
		  planar[kj][ki]->setPar(Vector2D(parvals4[2*ki],parvals4[2*ki+1]));
		  planar[kj][ki]->setSurfaceDist(dist_ang4[ki].first, dist_ang4[ki].second);
		  planar[kj][ki]->setRegion(adj_planar[kj].second);
		  adj_planar[kj].second->addPoint(planar[kj][ki]);
		}
	      adj_planar[kj].second->updateInfo();
	  
	      double maxd5, avd5;
	      int num_5, num2_5;
	      adj_planar[kj].second->getAccuracy(maxd5, avd5, num_5, num2_5);
	      int sf_flag2 =
		adj_planar[kj].second->defineSfFlag(0, tol, num_5,
						    num2_5, avd5, false);
	  
	      adj_planar[kj].second->setSurfaceFlag(sf_flag2); 
	    }
	}
  

      for (size_t ki=0; ki<adj_planar.size(); ++ki)
	{
	  int num = adj_planar[ki].second->numPoints();
	  for (int ka=0; ka<num; ++ka)
	    {
	      RevEngRegion *reg = adj_planar[ki].second->getPoint(ka)->region();
	      if (reg != adj_planar[ki].second)
		std::cout << "contextCylinder. Region mismatch " << reg << " " << adj_planar[ki].second << std::endl;
	    }
	}
      
      // Ensure that the point group is connected
      vector<vector<RevEngPoint*> > added_groups;
      splitRegion(added_groups);
      if (added_groups.size() > 0)
	out_groups.insert(out_groups.end(), added_groups.begin(), added_groups.end());
      
      updateInfo(tol, angtol);
    }
  else
    {
      for (size_t ki=0; ki<plane2cyl.size(); ++ki)
	if (plane2cyl[ki].size() > 0)
	  {
	    for (size_t kj=0; kj<plane2cyl[ki].size(); ++kj)
	      {
		plane2cyl[ki][kj]->setRegion(adj_planar[ki].second);
	      }
	  }
      for (size_t ki=0; ki<out_groups.size(); ++ki)
	for (size_t kj=0; kj<out_groups[ki].size(); ++kj)
	  out_groups[ki][kj]->setRegion(this);
      out_groups.clear();
    }
    
  return OKsurf;
  //return false;
}

//===========================================================================
bool RevEngRegion::contextTorus(Point mainaxis[3],
				double tol, int min_pt,  int min_pt_reg,
				double angtol, int prefer_elementary,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  //bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  Point pos, axis, Cx;
  vector<size_t> adj_ix;
  double R1, R2, R1_tor;
  double cyl_dom[4];
  bool outer;
  int plane_ix, cyl_ix;
  bool analyse_rotated = false;
  bool foundaxis = analyseTorusContext(adj_elem, tol, angtol, adj_ix,
				       plane_ix, cyl_ix, pos, axis, Cx,
				       R1, R2, cyl_dom, outer, analyse_rotated);
  if (!foundaxis)
    return false;

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);
  if (analyse_rotated)
    {
      vector<Point> rotated0;
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > group0;
      group0.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
      RevEngUtils::rotateToPlane(group0, Cx, axis, pos, rotated0);
#ifdef DEBUG_TORUSCONTEXT
      std::ofstream of3("rotated_pts_contexttor.g2");
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << rotated0.size() << std::endl;
      for (size_t kr=0; kr<rotated0.size(); ++kr)
	of3 << rotated0[kr] << std::endl;
      of3 << "410 1 0 4 0 0 255 255" << std::endl;
      of3 << "1" << std::endl;
      of3 << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
      of3 << "410 1 0 4 0 255 0 255" << std::endl;
      of3 << "1" << std::endl;
      of3 << pos-0.5*len*Cx << " " << pos+0.5*len*Cx << std::endl;
#endif
    }
  
  double eps = 1.0e-6;
  Point pos2 = pos - R2*axis;
  Point Cy = axis.cross(Cx);
  shared_ptr<Torus> tor0(new Torus(R1, R2, pos2, axis, Cx));
  R1_tor = R1;
  
#ifdef DEBUG_TORUSCONTEXT
  std::ofstream oft("torus_context.g2");
  tor0->writeStandardHeader(oft);
  tor0->write(oft);
#endif

  shared_ptr<ElementarySurface> cyl = adj_elem[cyl_ix].first;
  shared_ptr<ElementarySurface> plane = adj_elem[plane_ix].first;
  double rad = cyl->radius(0.0, 0.0);

  // Extract the points feasible for a torus adjacent to the
  // identified cylinder, and check if they constitute a point cloud
  // of sufficiant size
  double limit = 1.1*R2;
  vector<RevEngPoint*> inside, outside;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point onaxis = pos + ((xyz2-pos)*axis)*axis;
      double dd = xyz2.dist(onaxis);
      if (dd >= rad-limit && dd <= rad+limit)
	inside.push_back(group_points_[kr]);
      else
	outside.push_back(group_points_[kr]);
    }
  
  if ((int)inside.size() < min_pt_reg)
    return false;
  
  //#ifdef DEBUG_TORUSCONTEXT
  // Check accuracy of inside points with respect to the current
  // torus
  vector<RevEngPoint*> inpt0, outpt0;
  vector<pair<double, double> > dist_ang0;
  double maxd0, avd0;
  int num_in0, num2_in0;
  vector<double> parvals0;
  RevEngUtils::distToSurf(inside.begin(), inside.end(),
			  tor0, tol, maxd0, avd0, num_in0, num2_in0,
			  inpt0, outpt0, parvals0, dist_ang0,
			  angtol);
  
  //#endif

  // Identify points belonging to the adjacent cylinder
  // First rotate current point cloud according to the found axis
  vector<RevEngPoint*> cyl_pts, remaining;
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(inside.begin(), inside.end()));
  RevEngUtils::rotateToPlane(group, Cx, axis, pos, rotated);
  double angtol2 = angtol;
  if (adj_elem[cyl_ix].second->getSurfaceFlag() == PROBABLE_HELIX)
    angtol2 = 0.5*M_PI;  // Cannot expect angular accuracy for a helix
  RevEngUtils::extractLinearPoints(inside, rotated,
				   len, pos, axis, rad, axis, false,
				   tol, angtol2, dist_ang0,
				   cyl_pts, true, remaining);
  
#ifdef DEBUG_TORUSCONTEXT
  std::ofstream of3("rotated_pts_contexttor.g2");
  of3 << "400 1 0 4 255 0 0 255" << std::endl;
  of3 << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of3 << rotated[kr] << std::endl;
  of3 << "410 1 0 4 0 0 255 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
  of3 << "410 1 0 4 0 255 0 255" << std::endl;
  of3 << "1" << std::endl;
  of3 << pos-0.5*len*Cx << " " << pos+0.5*len*Cx << std::endl;

  std::ofstream ofl("lin_extract1.g2");
  if (cyl_pts.size() > 0)
    {
      ofl << "400 1 0 4 55 100 100 255" << std::endl;
      ofl << cyl_pts.size() << std::endl;
      for (size_t kr=0; kr<cyl_pts.size(); ++kr)
	ofl << cyl_pts[kr]->getPoint() << std::endl;
    }
  
  if (remaining.size() > 0)
    {
      ofl << "400 1 0 4 100 100 55 255" << std::endl;
      ofl << remaining.size() << std::endl;
      for (size_t kr=0; kr<remaining.size(); ++kr)
	ofl << remaining[kr]->getPoint() << std::endl;
    }
  
  if (outside.size() > 0)
    {
      ofl << "400 1 0 4 10 10 10 255" << std::endl;
      ofl << outside.size() << std::endl;
      for (size_t kr=0; kr<outside.size(); ++kr)
	ofl << outside[kr]->getPoint() << std::endl;
    }
#endif

   // Recompute torus
  double cpdist = R2;
  double up, vp, dd;
  Point cl;
  for (size_t kr=0; kr<cyl_pts.size(); ++kr)
    {
      Vector3D xyz = cyl_pts[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      adj_elem[plane_ix].first->closestPoint(xyz2, up, vp, cl, dd, eps);
      cpdist = std::min(cpdist, dd);
    }
  Point pos3 = pos - cpdist*axis;
  int sgn = (outer) ? -1 : 1;
  R1_tor = rad+sgn*cpdist;
  shared_ptr<Torus> tor(new Torus(R1_tor, cpdist, pos3, axis, Cx));

#ifdef DEBUG_TORUSCONTEXT
  tor->writeStandardHeader(oft);
  tor->write(oft);
#endif

  // Check accuracy of identified cylinder points with respect to the
  // adjacent cylinder
  vector<pair<double, double> > dist_ang_cyl;
  double maxd_cyl, avd_cyl;
  int num_in_cyl, num2_in_cyl;
  vector<double> parvals_cyl;
  if (cyl_pts.size() > 0)
    {
      vector<RevEngPoint*> inpt_cyl, outpt_cyl;
      RevEngUtils::distToSurf(cyl_pts.begin(), cyl_pts.end(),
			      cyl, tol, maxd_cyl, avd_cyl, num_in_cyl, num2_in_cyl,
			      inpt_cyl, outpt_cyl, parvals_cyl, dist_ang_cyl,
			      angtol);
      int sf_flag_cyl = defineSfFlag((int)cyl_pts.size(), 0, tol, num_in_cyl,
				     num2_in_cyl, avd_cyl, true);
      if (sf_flag_cyl >= ACCURACY_POOR)
	{
	  // Something unexpected happened. Return to previous stage
	  R1_tor = R1;
	  tor = tor0;
	  remaining = inside;
	  cyl_pts.clear();
	}
    }

  // Compute torus accuracy
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  double maxd, avd;
  int num_in, num2_in;
  vector<double> parvals;
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  tor, tol, maxd, avd, num_in, num2_in,
			  inpt, outpt, parvals, dist_ang,
			  angtol);

  if (hasSurface())
    {
      double fac = 1.5;
      int sf_flag_tor = defineSfFlag((int)remaining.size(), 0, tol,
				     num_in, num2_in, avd, false);
      if (surfflag_ < sf_flag_tor || avdist_ < fac*avd ||
	  (avdist_ <= tol && avd > tol))
	return false;  // Not an improvement
    }
  
  // Check if parts of this regions represents a likely torus
  double tor_lim = 0.3;
  if (num2_in < (int)(tor_lim*(double)remaining.size()))
    {
      // Check if a candidate segmentation strategy is found
      if (outside.size() > group_points_.size()/10 &&
	  outside.size() > min_pt_reg)
	{
	  shared_ptr<SegmentData> seg_info(new SegmentData(1, pos, axis,
							   rad-limit,
							   rad+limit));
	  seg_info_.push_back(seg_info);
	}
      
      return false;
    }
  
  // Extract planar points and additional cylindrical points from torus
  vector<RevEngPoint*> planar;
  vector<RevEngPoint*> remaining2;
  Point axis2 = cyl->direction();
  double pihalf = 0.5*M_PI;
  double tor_dom[4];   // Torus domain
  tor_dom[0] = tor_dom[2] = std::numeric_limits<double>::max();
  tor_dom[1] = tor_dom[3] = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<remaining.size(); ++ki)
    {
      Vector3D xyz = remaining[ki]->getPoint();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = remaining[ki]->getLocFuncNormal();
      Point norm2 = remaining[ki]->getTriangNormal();
      double dist = pos.dist(ptpos);
      double ang1 = axis.angle(norm);
      double ang2 = axis.angle(norm2);
      double ang3 = axis2.angle(norm);
      double ang4 = axis2.angle(norm2);
      if (((dist < R1_tor && outer) || (dist > R1_tor && (!outer))) &&
	  (ang1 <= angtol2 || ang2 <= angtol2))
	planar.push_back(remaining[ki]);
      else if (fabs(pihalf-ang3) <= dist_ang[ki].second ||
	       fabs(pihalf-ang4) <= dist_ang[ki].second ||
	       dist_ang[ki].first > tol)
	{
	  double upar, vpar, tdist;
	  Point close;
	  cyl->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  bool OKcyl = false;
	  if (tdist <= dist_ang[ki].first && upar >= cyl_dom[0] &&
	      upar <= cyl_dom[1] && vpar >= cyl_dom[2] && vpar <= cyl_dom[3])
	    cyl_pts.push_back(remaining[ki]);
	  else
	    {
	      remaining2.push_back(remaining[ki]);
	      tor_dom[0] = std::min(tor_dom[0], parvals[2*ki]);
	      tor_dom[1] = std::max(tor_dom[1], parvals[2*ki]);
	      tor_dom[2] = std::min(tor_dom[2], parvals[2*ki+1]);
	      tor_dom[3] = std::max(tor_dom[3], parvals[2*ki+1]);
	    }
	}
      else
	{
	  remaining2.push_back(remaining[ki]);
	  tor_dom[0] = std::min(tor_dom[0], parvals[2*ki]);
	  tor_dom[1] = std::max(tor_dom[1], parvals[2*ki]);
	  tor_dom[2] = std::min(tor_dom[2], parvals[2*ki+1]);
	  tor_dom[3] = std::max(tor_dom[3], parvals[2*ki+1]);
	}
    }
  
  // Plane to torus
  vector<RevEngPoint*>  plan2tor;
  for (auto it=adj_elem[plane_ix].second->pointsBegin();
       it!=adj_elem[plane_ix].second->pointsEnd(); ++it)
    {
      Vector3D xyz = (*it)->getPoint();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = (*it)->getLocFuncNormal();
      double dist = pos.dist(ptpos);
      double ang = axis.angle(norm);
      ang = std::min(ang, M_PI-ang);
      if ((dist > R1 && outer) || (dist < R1 && (!outer)))
	{
	  double upar, vpar, tdist;
	  Point close;
	  tor->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  Point tnorm;
	  tor->normal(tnorm, upar, vpar);
	  double ang2 = tnorm.angle(norm);
	  if (tdist < (*it)->getSurfaceDist() && ang2 <= angtol2 &&
	      upar >= tor_dom[0] && upar <= tor_dom[1])
	    plan2tor.push_back(*it);
	}
    }

   // Cylinder to torus
  vector<RevEngPoint*>  cyl2tor;
  for (auto it=adj_elem[cyl_ix].second->pointsBegin();
       it!=adj_elem[cyl_ix].second->pointsEnd(); ++it)
    {
      Vector3D xyz = (*it)->getPoint();
      Vector2D uv = (*it)->getPar();
      Point ptpos(xyz[0], xyz[1], xyz[2]);
      Point norm = (*it)->getLocFuncNormal();
      if (uv[0] < cyl_dom[0] || uv[0] > cyl_dom[1] ||
	  uv[1] < cyl_dom[2] || uv[1] > cyl_dom[3])
	{
	  double upar, vpar, tdist;
	  Point close;
	  tor->closestPoint(ptpos, upar, vpar, close, tdist, eps);
	  Point tnorm;
	  tor->normal(tnorm, upar, vpar);
	  double ang2 = tnorm.angle(norm);
	  if (tdist < (*it)->getSurfaceDist() && ang2 <= angtol2)
	    cyl2tor.push_back(*it);
	}
    }
  
#ifdef DEBUG_TORUSCONTEXT
  if (planar.size() > 0)
    {
      std::ofstream ofp("planar_move.g2");
      ofp << "400 1 0 4 255 0 0 255" << std::endl;
      ofp << planar.size() << std::endl;
      for (size_t kr=0; kr<planar.size(); ++kr)
	ofp << planar[kr]->getPoint() << std::endl;
    }

  if (plan2tor.size() > 0)
    {
      std::ofstream oftp("plan_tor_move.g2");
      oftp << "400 1 0 4 0 255 0 255" << std::endl;
      oftp << plan2tor.size() << std::endl;
      for (size_t kr=0; kr<plan2tor.size(); ++kr)
	oftp << plan2tor[kr]->getPoint() << std::endl;
    }
  
  if (cyl_pts.size() > 0)
    {
      std::ofstream ofc("cyl_move.g2");
      ofc << "400 1 0 4 255 0 0 255" << std::endl;
      ofc << cyl_pts.size() << std::endl;
      for (size_t kr=0; kr<cyl_pts.size(); ++kr)
	ofc << cyl_pts[kr]->getPoint() << std::endl;
    }

  if (cyl2tor.size() > 0)
    {
      std::ofstream oftc("cyl_tor_move.g2");
      oftc << "400 1 0 4 0 255 0 255" << std::endl;
      oftc << cyl2tor.size() << std::endl;
      for (size_t kr=0; kr<cyl2tor.size(); ++kr)
	oftc << cyl2tor[kr]->getPoint() << std::endl;
    }
#endif

  // Recompute accuracy
  vector<RevEngPoint*> inpt2, outpt2;
  vector<pair<double, double> > dist_ang2;
  double maxd2, avd2;
  int num_in2, num2_in2;
  vector<double> parvals2;
  RevEngUtils::distToSurf(remaining2.begin(), remaining2.end(),
			  tor, tol, maxd2, avd2, num_in2, num2_in2,
			  inpt2, outpt2, parvals2, dist_ang2,
			  angtol);
#ifdef DEBUG_TORUSCONTEXT
  std::ofstream ofd2("in_out_tor_context.g2");
  ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
  ofd2 << inpt2.size() << std::endl;
  for (size_t kr=0; kr<inpt2.size(); ++kr)
    ofd2 << inpt2[kr]->getPoint() << std::endl;
  ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
  ofd2 << outpt2.size() << std::endl;
  for (size_t kr=0; kr<outpt2.size(); ++kr)
    ofd2 << outpt2[kr]->getPoint() << std::endl;
#endif
      
  if (num2_in2 > min_pt && num2_in2 > (int)remaining2.size()/2)
    {
      // Check for deviant points at the boundary
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						    edge_class_type_, remaining2));
      
      vector<RevEngPoint*> dist_points;
      reg->identifyDistPoints(dist_ang2, tol, maxd2, avd2, dist_points);
      if (dist_points.size() > 0)
	{
	  reg->extractSpesPoints(dist_points, out_groups, true);

	  vector<RevEngPoint*> remaining3;
	  remaining3 = reg->getPoints();
	  for (size_t ki=0; ki<remaining3.size(); ++ki)
	    remaining3[ki]->setRegion(this);

	  std::swap(remaining2, remaining3);
	}
      else
	{
	  for (size_t ki=0; ki<remaining2.size(); ++ki)
	    remaining2[ki]->setRegion(this);
	}
    }

  if (plan2tor.size() > 0)
    {
      remaining2.insert(remaining2.end(), plan2tor.begin(), plan2tor.end());
    }
  
  if (cyl2tor.size() > 0)
     {
      remaining2.insert(remaining2.end(), cyl2tor.begin(), cyl2tor.end());
    }
   
  // Recompute accuracy
  vector<RevEngPoint*> inpt3, outpt3;
  vector<pair<double, double> > dist_ang3;
  double maxd3, avd3;
  int num_in3, num2_in3;
  vector<double> parvals3;
  RevEngUtils::distToSurf(remaining2.begin(), remaining2.end(),
			  tor, tol, maxd3, avd3, num_in3, num2_in3,
			  inpt3, outpt3, parvals3, dist_ang3,
			  angtol);
#ifdef DEBUG_TORUSCONTEXT
  std::ofstream ofd3("in_out_tor_context3.g2");
  ofd3 << "400 1 0 4 155 50 50 255" << std::endl;
  ofd3 << inpt3.size() << std::endl;
  for (size_t kr=0; kr<inpt3.size(); ++kr)
    ofd3 << inpt3[kr]->getPoint() << std::endl;
  ofd3 << "400 1 0 4 50 155 50 255" << std::endl;
  ofd3 << outpt3.size() << std::endl;
  for (size_t kr=0; kr<outpt3.size(); ++kr)
    ofd3 << outpt3[kr]->getPoint() << std::endl;
#endif
      
  // Check if this surface is better than an eventual previous surface. Move points
  int sf_flag = defineSfFlag((int)remaining2.size(), 0, tol, num_in3, num2_in3, avd3, false);
  bool OKsurf = (sf_flag < ACCURACY_POOR);
  bool hasSurf = (associated_sf_.size() > 0);
  if (hasSurf)
    {
      double acc_fac = 1.5;
      double maxd_init, avd_init;
      int num_init, num_init2;
      getAccuracy(maxd_init, avd_init, num_init, num_init2);
	      
      // Check with current approximating surface
      if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	OKsurf = false;
     else if (prefer_elementary == ALWAYS_ELEM ||
	       prefer_elementary == PREFER_ELEM)
	{
	  if ((double)num_init/(double)group_points_.size() >
	      (double)num_in3/(double)remaining.size() &&
	      (avd_init <= avd3 || num_in3 >= acc_fac*num_init))
	    OKsurf = false;
	}
     else
       {
	 if (!(num_in3 < num_init ||
	       (avd3 < avd_init && num_in3 < acc_fac*num_init)))
	   OKsurf = true;
       }
      if (sf_flag == ACCURACY_OK && surfflag_ > ACCURACY_OK)
	OKsurf = true;
    }

  if (OKsurf || hasSurf == false)
    {
      std::swap(group_points_, remaining2);
      if (plan2tor.size() > 0)
	{
	  for (size_t kj=0; kj<plan2tor.size(); ++kj)
	    {
	      adj_elem[plane_ix].second->removePoint(plan2tor[kj]);
	      plan2tor[kj]->setRegion(this);
	    }
	  adj_elem[plane_ix].second->updateInfo(tol, angtol);
	}
      
      if (cyl2tor.size() > 0)
	{
	  for (size_t kj=0; kj<cyl2tor.size(); ++kj)
	    {
	      adj_elem[cyl_ix].second->removePoint(cyl2tor[kj]);
	      cyl2tor[kj]->setRegion(this);
	    }
	  adj_elem[cyl_ix].second->updateInfo(tol, angtol);
	}
    }
  
  if (OKsurf)
    {
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals3[2*ki],parvals3[2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang3[ki].first, dist_ang3[ki].second);
	}
      setAccuracy(maxd3, avd3, num_in3, num2_in3);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }
  if (OKsurf || hasSurf == false)
    {
      if (planar.size() > 0)
	{
	  // Parameterize
	  vector<RevEngPoint*> inpt4, outpt4;
	  vector<pair<double, double> > dist_ang4;
	  double maxd4, avd4;
	  int num_in4, num2_in4;
	  vector<double> parvals4;
	  RevEngUtils::distToSurf(planar.begin(), planar.end(),
				  plane, tol, maxd4, avd4, num_in4, num2_in4,
				  inpt4, outpt4, parvals4, dist_ang4,
				  angtol);
	  for (size_t ki=0; ki<planar.size(); ++ki)
	    {
	      planar[ki]->setPar(Vector2D(parvals4[2*ki],parvals4[2*ki+1]));
	      planar[ki]->setSurfaceDist(dist_ang4[ki].first, dist_ang4[ki].second);
	      adj_elem[plane_ix].second->addPoint(planar[ki]);
	    }
	  adj_elem[plane_ix].second->updateInfo(tol, angtol);

	  double maxd5, avd5;
	  int num_in5, num2_in5;
	  adj_elem[plane_ix].second->getAccuracy(maxd5, avd5, num_in5, num2_in5);
	  int sf_flag5 = adj_elem[plane_ix].second->defineSfFlag(0, tol,
								 num_in5, num2_in5,
								 avd5, false);
	  adj_elem[plane_ix].second->setSurfaceFlag(sf_flag5);
	}
  
      if (cyl_pts.size() > 0)
	{
	  // Parameterize
	  vector<RevEngPoint*> inpt4, outpt4;
	  vector<pair<double, double> > dist_ang4;
	  double maxd4, avd4;
	  int num_in4, num2_in4;
	  vector<double> parvals4;
	  RevEngUtils::distToSurf(cyl_pts.begin(), cyl_pts.end(),
				  cyl, tol, maxd4, avd4, num_in4, num2_in4,
				  inpt4, outpt4, parvals4, dist_ang4,
				  angtol);
	  for (size_t ki=0; ki<cyl_pts.size(); ++ki)
	    {
	      cyl_pts[ki]->setPar(Vector2D(parvals4[2*ki],parvals4[2*ki+1]));
	      cyl_pts[ki]->setSurfaceDist(dist_ang4[ki].first, dist_ang4[ki].second);
	      adj_elem[cyl_ix].second->addPoint(cyl_pts[ki]);
	    }
	  adj_elem[cyl_ix].second->updateInfo(tol, angtol);

	  double maxd5, avd5;
	  int num_in5, num2_in5;
	  adj_elem[cyl_ix].second->getAccuracy(maxd5, avd5, num_in5, num2_in5);
	  int sf_flag5 = adj_elem[cyl_ix].second->defineSfFlag(0, tol,
							       num_in5, num2_in5,
							       avd5, true);
	  adj_elem[cyl_ix].second->setSurfaceFlag(sf_flag5);
	}

      if (numPoints() > 0)//group_points_.size() > 0)
	{
	  // Ensure that the point group is connected
	  vector<vector<RevEngPoint*> > added_groups;
	  splitRegion(added_groups);
	  if (added_groups.size() > 0)
	    out_groups.insert(out_groups.end(), added_groups.begin(), added_groups.end());
	  
	  updateInfo(tol, angtol);
	}
      else
	{
	  removeFromAdjacent();
	  clearRegionAdjacency();
	}
    }


  if (outside.size() > 0)
    {
      vector<std::vector<RevEngPoint*> > sep_outside;
      std::vector<RevEngPoint*> dummy;
      connectedGroups(outside, sep_outside, false, dummy);
      out_groups.insert(out_groups.end(), sep_outside.begin(),
			sep_outside.end());
    }
  
  return OKsurf;
 
 }

//===========================================================================
bool RevEngRegion::addPointsToGroup(vector<RevEngPoint*>& points,
				    double tol, double angtol, bool compute_accuracy)
//===========================================================================
{
  if (!hasSurface())
    return false;

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();

  if (compute_accuracy)
    {
      // Compute accuracy
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      double maxd, avd;
      int num_in, num2_in;
      vector<double> parvals;
      RevEngUtils::distToSurf(points.begin(), points.end(),
			      surf, tol, maxd, avd, num_in, num2_in,
			      inpt, outpt, parvals, dist_ang,
			      angtol);

      // Update info in points
      for (size_t kh=0; kh<points.size(); ++kh)
	{
	  points[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	  points[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	}
    }
  for (size_t kh=0; kh<points.size(); ++kh)
      points[kh]->setRegion(this);
  group_points_.insert(group_points_.end(), points.begin(), points.end());
  updateInfo(tol, angtol);

  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  int sf_flag = defineSfFlag(0, tol, num_inside_, num_inside2_, avdist_, cyllike);
  setSurfaceFlag(sf_flag);

  return true;
}



//===========================================================================
void RevEngRegion::growInDomain(RevEngRegion *adjacent, double tol,
				double angtol)
//===========================================================================
{
  if (!hasSurface())
    return;   // No growing is possible
  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2;
  if (adjacent->hasSurface())
    surf2 = adjacent->getSurface(0)->surface();

  // Get seed points
  vector<RevEngPoint*> pts = adjacent->extractNextToAdjacent(this);

  vector<RevEngPoint*> move;
  double eps = 1.0e-6;
  for (size_t ki=0; ki<pts.size(); ++ki)
    {
      pts[ki]->setVisited();
      Vector3D xyz = pts[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pos, upar, vpar, close, dist, eps);
      bool in_domain = surf->inDomain(upar, vpar, eps);
      bool at_boundary = surf->onBoundary(upar, vpar, eps);
      if (in_domain && (!at_boundary))
	{
	  bool do_move = true;
	  if (surf2.get())
	    {
	      Vector2D uv = pts[ki]->getPar();
	      if (surf2->inDomain(uv[0], uv[1], eps) &&
		  dist >= pts[ki]->getSurfaceDist())
		do_move = false;
	    }
	  if (do_move)
	    {
	      move.push_back(pts[ki]);
	      vector<ftSamplePoint*> next = pts[ki]->getNeighbours();
	      for (size_t kj=0; kj<next.size(); ++kj)
		{
		  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
		  if (curr->visited())
		    continue;
		  RevEngRegion *adj_reg = curr->region();
		  if (adj_reg != adjacent)
		    continue;
		  curr->setVisited();
		  pts.push_back(curr);
		}
	    }
	}
    }

  for (size_t ki=0; ki<pts.size(); ++ki)
    pts[ki]->unsetVisited();

  if (move.size() > 0)
    {
      adjacent->removePoints(move);
      adjacent->updateInfo(tol, angtol);

      (void)addPointsToGroup(move, tol, angtol);
    }
      

}

//===========================================================================
void RevEngRegion::growFromNeighbour(Point mainaxis[3], int min_pt_reg,
				     vector<RevEngPoint*>& seed, double tol,
				     double angtol, RevEngRegion *neighbour,
				     bool do_update)
//===========================================================================
{
  if (!hasSurface())
    return;   // No growing is possible

  bool check_pts = (neighbour->hasSurface());

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  Point axis;
  if (elem.get())
    axis = elem->direction();
  double axis_ang = 0.0;
  if (surf->instanceType() == Class_Cylinder)
    axis_ang = 0.5*M_PI;   // To be extended as appropriate
  else if (surf->instanceType() == Class_Cone)
    {
      double cone_ang = ((Cone*)(elem.get()))->getConeAngle();
      axis_ang = 0.5*M_PI - cone_ang;
    }

#ifdef DEBUG_GROWNEIGHBOUR
  std::ofstream ofs1("seed1.g2");
  ofs1 << "400 1 0 4 0 200 55 255" << std::endl;
  ofs1 << seed.size() << std::endl;
  for (size_t ki=0; ki<seed.size(); ++ki)
    {
      Vector3D xyz = seed[ki]->getPoint();
      ofs1 << xyz << std::endl;
    }
#endif
  
  // Check initial points
  double eps = 1.0e-6;
  double tol2 = std::max(tol, 0.75*maxdist_);
  tol2 = std::min(tol2, 2.0*tol);
  double tol3 = 0.5*tol;
  //double angtol2 = 0.5*angtol;
  vector<RevEngPoint*> next_pts;
  double upar, vpar, dist, ang, ang2;
  Point close;
  double dfac = 1.5;
  for (size_t ki=0; ki<seed.size(); ++ki)
    {
      seed[ki]->setVisited();
      Vector3D xyz = seed[ki]->getPoint();
      Point norm = seed[ki]->getLocFuncNormal();
      Point norm2 = seed[ki]->getTriangNormal();
      surf->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
      if (!elem.get())
	surf->normal(axis, upar, vpar);
      ang = axis.angle(norm);
      ang = fabs(ang-axis_ang);
      ang2 = axis.angle(norm2);
      ang2 = fabs(ang2-axis_ang);
      ang = std::min(std::min(ang, M_PI-ang), std::min(ang2, M_PI-ang2));
      if (dist <= tol3 || (dist <= tol2 && ang <= angtol))
	{
	  if (check_pts)
	    {
	      double init_dist, init_ang;
	      seed[ki]->getSurfaceDist(init_dist, init_ang);
	      if ((dist < init_dist && ang < init_ang) ||
		  (dfac*dist < init_dist && ang < dfac*init_ang))
		next_pts.push_back(seed[ki]);
	    }
	  else
	    next_pts.push_back(seed[ki]);
	}
    }

#ifdef DEBUG_GROWNEIGHBOUR
  std::ofstream of0("curr_region_from.g2");
  neighbour->writeRegionPoints(of0);
#endif

#ifdef DEBUG_GROWNEIGHBOUR
  std::ofstream ofs2("seed2.g2");
  ofs2 << "400 1 0 4 0 200 55 255" << std::endl;
  ofs2 << next_pts.size() << std::endl;
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      Vector3D xyz = next_pts[ki]->getPoint();
      ofs2 << xyz << std::endl;
    }
#endif
  
  // Grow
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      vector<ftSamplePoint*> next2 = next_pts[ki]->getNeighbours();
      for (size_t kj=0; kj<next2.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next2[kj]);
	  if (curr->visited())
	    continue;
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != neighbour)
	    continue;
	  curr->setVisited();
	  Vector3D xyz = curr->getPoint();
	  Point norm = curr->getLocFuncNormal();
	  Point norm2 = curr->getTriangNormal();
	  surf->closestPoint(Point(xyz[0],xyz[1],xyz[2]), upar, vpar, close, dist, eps);
	  if (!elem.get())
	    surf->normal(axis, upar, vpar);
	  ang = axis.angle(norm);
	  ang = fabs(ang-axis_ang);
	  ang2 = axis.angle(norm2);
	  ang2 = fabs(ang2-axis_ang);
	  ang = std::min(std::min(ang, M_PI-ang), std::min(ang2, M_PI-ang2));
	  if (dist <= tol3 || (dist <= tol2 && ang <= angtol))
	    {
	      if (check_pts)
		{
		  double init_dist, init_ang;
		  curr->getSurfaceDist(init_dist, init_ang);
		  if ((dist < init_dist && ang < init_ang) ||
		      (dfac*dist < init_dist && ang < dfac*init_ang))
		    next_pts.push_back(curr);
		}
	      else
		next_pts.push_back(curr);
	    }
	}
    }
  for (int ka=0; ka<neighbour->numPoints(); ++ka)
    neighbour->getPoint(ka)->unsetVisited();

  int classtype = surf->instanceType();
  if (classtype == Class_Cylinder)
    {
      double dom[4];
      Vector2D uv = group_points_[0]->getPar();
      dom[0] = dom[1] = uv[0];
      dom[2] = dom[3] = uv[1];
      double dist, ang, avang;
      double fac = 1.0/(double)group_points_.size();
      group_points_[0]->getSurfaceDist(dist, ang);
      avang = fac*ang;
      for (size_t ki=1; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  group_points_[ki]->getSurfaceDist(dist, ang);
	  dom[0] = std::min(dom[0], uv[0]);
	  dom[1] = std::max(dom[1], uv[0]);
	  dom[2] = std::min(dom[2], uv[1]);
	  dom[3] = std::max(dom[3], uv[1]);
	  avang += fac*ang;
	}
      
      vector<RevEngPoint*> adjpts;
      vector<double> par_and_dist;
      double avd, ava;
      int nn;
      getAdjInsideDist(surf, dom, tol, neighbour, avd, ava, nn, adjpts,
		       par_and_dist);

#ifdef DEBUG_GROWNEIGHBOUR2
      std::ofstream of2("in_cyl_pts.g2");
      of2 << "400 1 0 4 75 75 75 255" << std::endl;
      of2 << adjpts.size() << std::endl;
      for (size_t kh=0; kh<adjpts.size(); ++kh)
	of2 << adjpts[kh]->getPoint() << std::endl;
#endif
      int stop_cyl = 1;
    }
  
  // Move points
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      //addPoint(next_pts[ki]);
      neighbour->removePoint(next_pts[ki]);
    }

#ifdef DEBUG_GROWNEIGHBOUR2
  std::ofstream of("updated_region_adj.g2");
  writeRegionInfo(of);
#endif
  
  // Update this region
  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  double frac = 0.1;
  int num_init = (int)group_points_.size();
  if (next_pts.size() > 0)
    {
      // Compute distance
      shared_ptr<ParamSurface> surf = getSurface(0)->surface();
      double maxdist, avdist;
      int num_in, num2_in;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      vector<RevEngPoint*> inpt, outpt; 
      RevEngUtils::distToSurf(next_pts.begin(), next_pts.end(),
			  surf, tol, maxdist, avdist, num_in, num2_in,
			  inpt, outpt, parvals, dist_ang, angtol);

      for (size_t kh=0; kh<next_pts.size(); ++kh)
	{
	  next_pts[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	  next_pts[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	  next_pts[kh]->setRegion(this);
	}
      group_points_.insert(group_points_.end(), next_pts.begin(),
			   next_pts.end());
      updateInfo(tol, angtol);
      int sf_flag = defineSfFlag(0, tol, num_inside_,
				 num_inside2_, avdist_, cyllike);
      setSurfaceFlag(sf_flag);
      if (do_update && surf_adaption_ == INITIAL &&
	  (double)next_pts.size() > frac*(double)num_init)
	checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
      for (size_t kj=0; kj<rev_edges_.size(); ++kj)
	rev_edges_[kj]->increaseExtendCount();

      neighbour->updateInfo(tol, angtol);
    }

  int stop_break = 1;
}

//===========================================================================
vector<RevEngRegion*> RevEngRegion::commonAdjacent(RevEngRegion* adj_reg)
//===========================================================================
{
  vector<RevEngRegion*> common;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      bool adjacent = (*it)->isAdjacent(adj_reg);
      if (adjacent)
	common.push_back(*it);
    }
  return common;
}

//===========================================================================
vector<RevEngPoint*> RevEngRegion::extractBdPoints()
//===========================================================================
{
  vector<RevEngPoint*> result;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != this)
	    {
	      result.push_back(group_points_[ki]);
	      break;
	    }
	}
    }
  return result;
}

//===========================================================================
vector<RevEngPoint*>
RevEngRegion::extractBdPoints(vector<RevEngRegion*> regions)
//===========================================================================
{
  vector<RevEngPoint*> result;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != this)
	    {
	      size_t kr;
	      for (kr=0; kr<regions.size(); ++kr)
		if (adj_reg == regions[kr])
		  break;

	      if (kr < regions.size())
		result.push_back(group_points_[ki]);
	      break;
	    }
	}
    }
  return result;
}

//===========================================================================
vector<RevEngPoint*> RevEngRegion::extractBranchPoints()
//===========================================================================
{
  vector<RevEngPoint*> result;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      int num_other = 0;
      RevEngRegion *other = 0;
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != this)
	    {
	      if (num_other == 0)
		{
		  num_other++;
		  other = adj_reg;
		}
	      else if (adj_reg != other)
		{
		  result.push_back(group_points_[ki]);
		  break;
		}
	    }
	}
    }
  return result;
}

//===========================================================================
vector<RevEngPoint*> RevEngRegion::extractNextToAdjacent(RevEngRegion* reg)
//===========================================================================
{
  vector<RevEngPoint*> result;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg == reg)
	    {
	      result.push_back(group_points_[ki]);
	      break;
	    }
	}
    }
  return result;
}

struct DistanceInfo
{
  double mind_, maxd_, avd_;
  double pmin_, pmax_;
  int nmb_;

  DistanceInfo()
  {
    mind_ = maxd_ = avd_ = pmin_ = pmax_ = 0.0;
    nmb_ = 0;
  }
  
  DistanceInfo(double mind, double maxd, double avd,
	       double pmin, double pmax, int nmb)
  {
    mind_ = mind;
    maxd_ = maxd;
    avd_ = avd;
    pmin_ = pmin;
    pmax_ = pmax;
    nmb_ = nmb;
  }

  void setInfo(double mind, double maxd, double avd,
	       double pmin, double pmax, int nmb)
  {
    mind_ = mind;
    maxd_ = maxd;
    avd_ = avd;
    pmin_ = pmin;
    pmax_ = pmax;
    nmb_ = nmb;
  }
};

void RevEngRegion::testBlendGap(vector<shared_ptr<CurveOnSurface> >& cvs,
				double tmin, double tmax, double tdel, double width,
				vector<pair<double,double> >& not_gap)
{
  double cdel0 = tmax - tmin;
  int ncheck = std::max((int)(1.5*cdel0/tdel), 1);
  double cdel = cdel0/(double)ncheck;
  double cpar;
  size_t kr;
  vector<size_t> inside;
  for (kr=0, cpar=tmin+0.5*cdel; kr<ncheck; ++kr, cpar+=cdel)
    {
      size_t kh;
      for (kh=0; kh<cvs.size(); ++kh)
	if (cvs[kh]->startparam() <= cpar && cvs[kh]->endparam() >= cpar)
	  break;
      if (kh >= cvs.size())
	{
	  std::cout << cvs.size() << " " << kh << std::endl;
	}
      Point cpt = cvs[kh]->ParamCurve::point(cpar);
      double dist;
      RevEngPoint* rpt = closestPoint(cpt, dist);
      if (dist <= width)
	inside.push_back(kr);
      else
	{
	  Vector3D xyz = rpt->getPoint();
	  Point rpt2(xyz[0], xyz[1], xyz[2]);
	  double tpar2, dist2;
	  Point close;
	  cvs[kh]->closestPoint(rpt2, std::max(cvs[kh]->startparam(),cpar-0.25*cdel),
				std::min(cvs[kh]->endparam(),cpar+0.25*cdel),
				tpar2, close, dist2);
	  if (dist2 <= width)
	    inside.push_back(kr);
	}
    }
  if ((int)inside.size() > ncheck/2)
    not_gap.push_back(std::make_pair(tmin, tmax));
  else if (inside.size() > 0)
    {
      double t1 = tmin + inside[0]*cdel;
      double t2 = t1 + cdel;
      for (size_t kj=1; kj<inside.size(); ++kj)
	{
	  if (inside[kj] == inside[kj-1]+1)
	    t2 += cdel;
	  else
	    {
	      not_gap.push_back(std::make_pair(t1, t2));
	      t1 = tmin + inside[kj]*cdel;
	      t2 = t1 + cdel;
	    }
	}
      not_gap.push_back(std::make_pair(t1, t2));
    }
}

//===========================================================================
void RevEngRegion::estimateBlendDimensions(vector<shared_ptr<CurveOnSurface> >& cvs,
					   vector<RevEngPoint*>& bd_points,
					   double tol, double distlim,
					   vector<pair<double,double> >& t1_t2,
					   vector<double>& width, int& num_in_lim)
//===========================================================================
{
  // if (cvs.size() != 1)
  //   return;  // Must consider this later. Will it happen?
  
  // Parameterize points with respect to the given curve
  double eps = 1.0e-9;
  vector<double> param(bd_points.size());
  vector<double> distance(bd_points.size());
  
  double tmin = cvs[0]->startparam();
  double tmax = cvs[cvs.size()-1]->endparam();
  num_in_lim = 0;
  for (size_t ki=0; ki<bd_points.size(); ++ki)
    {
      double tpar, dist=std::numeric_limits<double>::max();;
      Vector3D xyz = bd_points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  double tpar0, dist0;
	  Point close;
	  cvs[kj]->closestPoint(pt, cvs[kj]->startparam(), cvs[kj]->endparam(),
				tpar0, close, dist0);
	  if (dist0 < dist)
	    {
	      tpar = tpar0;
	      dist = dist0;
	    }
	}
      param[ki] = tpar;
      distance[ki] = dist;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
      if (dist <= distlim)
	num_in_lim++;
    }

  if (tmax < tmin+eps)
    return;
  
  // Sort
  for (size_t ki=0; ki<bd_points.size(); ++ki)
    for (size_t kj=ki+1; kj<bd_points.size(); ++kj)
      if (param[kj] < param[ki])
	{
	  std::swap(bd_points[ki], bd_points[kj]);
	  std::swap(param[ki], param[kj]);
	  std::swap(distance[ki], distance[kj]);
	}

  // Limit point sequence
  size_t first = 0, last = param.size()-1;
  eps = std::min(1.0e-9, 1.0e-6*(tmax-tmin));
  for (; first<last; ++first)
    if (param[first+1]-param[0] > eps)
      break;
  for (; last>first; --last)
    if (param[param.size()-1]-param[last-1] > eps)
      break;

  double mean_out = 0.0;
  double min_out = std::numeric_limits<double>::max();
  int nmb_out = 0;
  size_t outlim = distance.size()/10;
  for (size_t ki=0; ki<std::max(first,outlim); ++ki)
    {
      min_out = std::min(min_out, distance[ki]);
      mean_out += distance[ki];
      nmb_out++;
    }
  for (size_t ki=std::min(distance.size()-outlim,last+1);
       ki<distance.size(); ++ki)
    {
      min_out = std::min(min_out, distance[ki]);
      mean_out += distance[ki];
      nmb_out++;
    }
  if (nmb_out > 0)
    mean_out /= (double)nmb_out;
  else
    min_out = 10.0*distlim;
  
  // Find start point for search
  double max_dist = 0.0;
  double min_dist = std::numeric_limits<double>::max();
  int min_ix = -1;
  for (size_t ki=first; ki<=last; ++ki)
    {
      if (distance[ki] < min_dist)
	{
	  min_dist = distance[ki];
	  min_ix = (int)ki;
	}
      max_dist = std::max(max_dist,distance[ki]);
    }

  int ndel = std::max(20, (int)(last-first)/100);
  ndel = std::max(1, std::min(ndel, (int)(last-first)/5));
  vector<DistanceInfo> distinfo(ndel);
  double dlim = std::max(2.0*distlim, 0.5*(min_out + mean_out));
  double tdel = (tmax-tmin)/(double)ndel;
  double tpar;
  int ka, kb, ix;
  for (ix=0, ka=first, tpar=tmin+tdel; ix<ndel && ka<=(int)last;
       ++ix, ka=kb, tpar+=tdel)
    {
      int nmb = 0;
      double maxd=0.0, avd=0.0;
      double mind = 2.0*dlim;
      for (kb=ka; kb<=(int)last && param[kb]<=tpar; ++kb)
	{
	  if (distance[kb]>dlim)
	    continue;
	  ++nmb;
	  maxd = std::max(maxd, distance[kb]);
	  mind = std::min(mind, distance[kb]);
	  avd += distance[kb];
	}
      if (nmb > 0)
	avd /= (double)nmb;
      
      distinfo[ix].setInfo(mind, maxd, avd, tpar-tdel, tpar, nmb);
    }
  
  double del = 0.0;
  for (size_t ki=0; ki<distinfo.size(); ++ki)
    del += (distinfo[ki].avd_ - distinfo[ki].mind_);
  del /= (double)(distinfo.size());
  del = std::min(del, mean_out);

  size_t kj, kr;
  int pt_out_lim = 5;
  double eps2 = 0.01*(mean_out - min_dist);
  for (kj=0; kj<distinfo.size(); kj=kr+1)
    {
      double curr_width = 0.0;
      double tmin2 = distinfo[kj].pmin_;
      double tmax2 = distinfo[kj].pmax_;
      int num_pt = 0;
      double min_width = max_dist;
      for (kr=kj; kr<distinfo.size(); ++kr)
	{
	  if (distinfo[kr].nmb_ == 0)
	    break;
	  if (distinfo[kr].mind_ > distlim)
	    break;
	  num_pt += distinfo[kr].nmb_;

	  min_width = std::min(min_width, distinfo[kr].mind_);
	  curr_width = std::max(curr_width, std::min(distinfo[kr].maxd_,
						     distinfo[kr].mind_+del));
	  tmax2 = distinfo[kr].pmax_;
	}
      if (num_pt > 0 && min_width <= distlim)
	{
	  curr_width = std::min(curr_width, mean_out + 0.1*(max_dist-mean_out));
	  curr_width = std::max(curr_width, 0.2*distlim);
	  curr_width = std::max(curr_width, min_width+eps2);
	  t1_t2.push_back(std::make_pair(tmin2, tmax2));
	  width.push_back(curr_width);
	}
    }

  // The boundary points may be too sparse. Check if the gaps between
  // curve pieces are real
  if (t1_t2.size() > 0)
    {
      size_t nmb = t1_t2.size();
      vector<pair<double,double> > t3_t4;
      vector<double> width2;
      if (t1_t2[0].first > tmin+eps)
	{
	  vector<pair<double,double> > not_gap;
	  testBlendGap(cvs, tmin, t1_t2[0].first, tdel, width[0], not_gap);
	  if (not_gap.size() > 0)
	    {
	      t3_t4.insert(t3_t4.end(), not_gap.begin(), not_gap.end());
	      for (size_t kh=0; kh<not_gap.size(); ++kh)
		width2.push_back(width[0]);
	    }
	}
      t3_t4.push_back(t1_t2[0]);
      width2.push_back(width[0]);
      for (size_t kj=1; kj<t1_t2.size(); ++kj)
	{
	  vector<pair<double,double> > not_gap;
	  double wdt = 0.5*(width[kj-1]+width[kj]);
	  testBlendGap(cvs, t1_t2[kj-1].second, t1_t2[kj].first, tdel, wdt, not_gap);
	  if (not_gap.size() > 0)
	    {
	      t3_t4.insert(t3_t4.end(), not_gap.begin(), not_gap.end());
	      for (size_t kh=0; kh<not_gap.size(); ++kh)
		width2.push_back(wdt);
	    }
	  t3_t4.push_back(t1_t2[kj]);
	  width2.push_back(width[kj]);
	}
     if (tmax > t1_t2[nmb-1].second)
	{
	  vector<pair<double,double> > not_gap;
	  testBlendGap(cvs, t1_t2[nmb-1].second, tmax, tdel, width[nmb-1], not_gap);
	  if (not_gap.size() > 0)
	    {
	      t3_t4.insert(t3_t4.end(), not_gap.begin(), not_gap.end());
	      for (size_t kh=0; kh<not_gap.size(); ++kh)
		width2.push_back(width[nmb-1]);
	    }
	}
      if (t3_t4.size() > t1_t2.size())
	{
	  double ltol = 0.5*tdel;
	  for (int ka=(int)t3_t4.size()-1; ka>0; --ka)
	    {
	      if (t3_t4[ka].first - t3_t4[ka-1].second <= ltol)
		{
		  t3_t4[ka-1].second = t3_t4[ka].second;
		  width2[ka-1] = 0.5*(width2[ka-1] + width2[ka]);
		  t3_t4.erase(t3_t4.begin()+ka);
		  width2.erase(width2.begin()+ka);
		}
	    }
	  std::swap(t1_t2, t3_t4);
	  std::swap(width, width2);
	}
    }
  
  int stop_break = 1;
}


//===========================================================================
void RevEngRegion::getNearPoints(shared_ptr<CurveOnSurface>& cv,
				 double& tmin, double& tmax,
				 double width, double angtol, 
				 vector<RevEngPoint*>& nearpoints)
//===========================================================================
{
  // if (cvs.size() != 1)
  //   return;  // Must consider this later. Will it happen?

  double pihalf = 0.5*M_PI;
  double tmin2 = tmax;
  double tmax2 = tmin;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      double tpar, dist=std::numeric_limits<double>::max();
      Point close;
      Vector3D xyz = group_points_[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      cv->closestPoint(pt, cv->startparam(), cv->endparam(),
		       tpar, close, dist);

      vector<Point> der(2);
      cv->point(der, tpar, 1);
      Point vec = pt - close;
      double ang = vec.angle(der[1]);
      Point norm1 = group_points_[ki]->getLocFuncNormal();
      Point norm2 = group_points_[ki]->getTriangNormal();
      double ang1 = norm1.angle(der[1]);
      double ang2 = norm2.angle(der[1]);
      ang1 = std::min(ang1, M_PI-ang1);
      ang2 = std::min(ang2, M_PI-ang2);
      if (tpar > tmin && tpar < tmax && dist <= width &&
	  fabs(pihalf-ang) <= angtol && std::min(ang1,ang2) > angtol)
	{
	  nearpoints.push_back(group_points_[ki]);
	  tmin2 = std::min(tmin2, tpar);
	  tmax2 = std::max(tmax2, tpar);
	}
    }
  tmin = tmin2;
  tmax = tmax2;
}

//===========================================================================
void RevEngRegion::getNearPoints2(vector<RevEngPoint*>& points,
				  shared_ptr<CurveOnSurface>& cv,
				  double width,
				  vector<RevEngPoint*>& nearpoints)
//===========================================================================
{
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      double tpar, dist;
      Point close;
      cv->closestPoint(pt, cv->startparam(), cv->endparam(),
		       tpar, close, dist);

      if (dist <= width)
	nearpoints.push_back(points[ki]);
    }
}

  
  
//===========================================================================
vector<RevEngPoint*>
RevEngRegion::removeOutOfSurf(vector<RevEngPoint*>& points,
			      double tol, double angtol, bool outer,
			      double& min_dist)
//===========================================================================
{
  vector<RevEngPoint*> result;
  double diag = bbox_.low().dist(bbox_.high());
  min_dist = 10.0*diag;
  if (!hasSurface())
    return result;

  int code;
  int classtype = getSurface(0)->instanceType(code);
  if (classtype != Class_Plane && classtype != Class_Cylinder &&
      classtype != Class_Cone)
    return result;  // For the time being

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);

  Point vec = elem->direction();
  Point loc = elem->location();
  double rad = elem->radius(0.0, 0.0);
  double alpha = 0.0;
  if (elem->instanceType() == Class_Cone)
    {
      shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ElementarySurface>(elem);
      alpha = cone->getConeAngle();
    }

  // Check orientation. Apply to plane
  bool plane = (elem->instanceType() == Class_Plane);
  if (plane)
    {
      Point normal = group_points_[0]->getTriangNormal();
      if (vec*normal < 0.0)
	vec *= -1;
      if (!outer)
	vec *= -1;
    }

  double upar, vpar, dist;
  Point close;
  double eps = 1.0e-6;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      if (plane)
	{
	  elem->closestPoint(pt, upar, vpar, close, dist, eps);
	  if (dist < tol || (pt-close)*vec > 0.0)
	    {
	      min_dist = std::min(min_dist, dist);
	      result.push_back(points[ki]);
	    }
	}
      else
	{
	  double vval = (pt-loc)*vec;
	  Point pt2 = pt - vval*vec;
	  dist = loc.dist(pt2);
	  double dist2 = vval*tan(alpha);
	  dist -= dist2;
	  if ((outer && dist < rad) || (!outer && dist > rad))
	    {
	      min_dist = std::min(min_dist, dist);
	      result.push_back(points[ki]);
	    }
	}
    }
  return result;
}

//===========================================================================
bool
RevEngRegion::analyseCylinderContext(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
				     double tol, double angtol, Point mainaxis[3], int mode,
				     vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_planar,
				     Point& pos, Point& axis, Point& Cx, double& rad,
				     vector<RevEngPoint*>& cyl_pts)
//===========================================================================
{
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].first->instanceType() == Class_Plane)
	adj_planar.push_back(adj[ki]);
      if (mode > 1 && adj[ki].first->instanceType() == Class_Cylinder)
	adj_planar.push_back(adj[ki]);
    }
  
  if (adj_planar.size() < 2)
    return false;
  
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_planar.size(); ++kj)
      if (adj_planar[kj].second->numPoints() > adj_planar[ki].second->numPoints())
	std::swap(adj_planar[ki], adj_planar[kj]);

#ifdef DEBUG_CYLCONTEXT
  std::ofstream ofadj("adj_planes.g2");
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      adj_planar[ki].first->writeStandardHeader(ofadj);
      adj_planar[ki].first->write(ofadj);
    }
#endif

  // Find feasible axis directions
  vector<Point> dir;
  double acute_lim = 0.25*M_PI;
  for (size_t ki=0; ki<adj_planar.size(); ++ki)
    {
      Point dir1 = adj_planar[ki].first->direction();
      if (mode > 1)
	dir.push_back(dir1);
      for (size_t kj=ki+1; kj<adj_planar.size(); ++kj)
	{
	  Point dir2 = adj_planar[kj].first->direction();
	  double ang = dir1.angle(dir2);
	  ang = std::min(ang, M_PI-ang);
	  if (ang > angtol && ang < acute_lim)
	    continue;

	  if (mode == 1 && ang <= angtol)
	    {
	      if (dir1*dir2 < 0.0)
		dir2 *= -1;
	      dir.push_back(0.5*(dir1 + dir2));
	    }
	  if (ang > angtol)
	    dir.push_back(dir1.cross(dir2));
	}
    }

  // Clean up in directions
  for (size_t ki=0; ki<dir.size(); ++ki)
    for (size_t kj=ki+1; kj<dir.size(); )
      {
	double ang = dir[ki].angle(dir[kj]);
	ang = std::min(ang, M_PI-ang);
	if (ang < angtol)
	  dir.erase(dir.begin()+kj);
	else
	  ++kj;
      }
  if (dir.size() < 1)
    return false;

  // Find best candidate direction
  double ang_lim = 0.1*M_PI;
  double pihalf = 0.5*M_PI;
  for (size_t ki=0; ki<dir.size(); ++ki)
    dir[ki].normalize();

  // Check feasability of cylinder
  vector<vector<RevEngPoint*> > dir_pts(dir.size());
  for (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      Point norm = group_points_[kj]->getLocFuncNormal();
      Point norm2 = group_points_[kj]->getTriangNormal();
      for (size_t ki=0; ki<dir.size(); ++ki)
	{
	  double ang = dir[ki].angle(norm);
	  double ang2 = dir[ki].angle(norm2);
	  if (std::min(fabs(pihalf-ang2), fabs(pihalf-ang)) < ang_lim)
	    dir_pts[ki].push_back(group_points_[kj]);
	}
    }

  vector<double> MAH(dir.size(), 0.0);
  vector<double> MAK(dir.size(), 0.0);
  for (size_t ki=0; ki<dir.size(); ++ki)
    {
      for (size_t kj=0; kj<dir_pts[ki].size(); ++kj)
	{
	  MAK[ki] += fabs(dir_pts[ki][kj]->GaussCurvature());
	  MAH[ki] += fabs(dir_pts[ki][kj]->meanCurvature());
	}
      MAK[ki] /= (double)dir_pts[ki].size();
      MAH[ki] /= (double)dir_pts[ki].size();
    }

  int min_frac = std::numeric_limits<double>::max();
  int ix = -1;
  double eps = 1.0e-9;
  double num_frac = 0.25;
  double min_num = num_frac*(double)group_points_.size();
  for (size_t ki=0; ki<dir.size(); ++ki)
    {
      if ((double)dir_pts[ki].size() < min_num)
	continue;
      double frac = (MAH[ki]<eps) ? 1.0e8 : MAK[ki]/MAH[ki];
      if (frac < min_frac)
	{
	  min_frac = frac;
	  ix = (int)ki;
	}
    }

  if (ix < 0)
    return false;

#ifdef DEBUG_CYLCONTEXT
  std::ofstream of2("cyl_cand_pts.g2");
  of2 << "400 1 0 4 0 255 0 255" << std::endl;
  of2 << dir_pts[ix].size() << std::endl;
  for (size_t kj=0; kj<dir_pts[ix].size(); ++kj)
    of2 << dir_pts[ix][kj]->getPoint() << std::endl;
#endif

  vector<vector<RevEngPoint*> > conn_groups;
  vector<RevEngPoint*> dummy;
  connectedGroups(dir_pts[ix], conn_groups, false, dummy);
  int min_conn = 0;
  int conn_ix = -1;
  for (size_t ki=0; ki<conn_groups.size(); ++ki)
    if ((int)conn_groups[ki].size() > min_conn)
      {
	min_conn = (int)conn_groups[ki].size();
	conn_ix = (int)ki;
      }
  
  if ((double)conn_groups[conn_ix].size() >= min_num)
    {
      // Candidate found
      // Compute cylinder centre and radius. Project points
      axis = dir[ix];
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > group;
      group.push_back(std::make_pair(conn_groups[conn_ix].begin(),
				     conn_groups[conn_ix].end()));
      Point mid = 0.5*(bbox_.low() + bbox_.high());
      vector<Point> projected;
      double maxdp, avdp;
      RevEngUtils::projectToPlane(group, axis, mid, projected, maxdp, avdp);
      
      int min_ix = -1;
      double min_ang = M_PI;
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = mainaxis[ka].angle(axis);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < min_ang)
	    {
	      min_ix = ka;
	      min_ang = ang;
	    }
	}
      Cx = mainaxis[(min_ix+1)%3];
      Point Cy = axis.cross(Cx);
      Cy.normalize();
      Cx = Cy.cross(axis);
      Cx.normalize();

      RevEngUtils::computeCircPosRadius(projected, axis, Cx, Cy, pos, rad);

      cyl_pts.insert(cyl_pts.end(), conn_groups[conn_ix].begin(),
		     conn_groups[conn_ix].end());
      return true;
    }

  return false;
}

//===========================================================================
bool
RevEngRegion::identifySignificantAxis(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
				      Point& pos, Point& axis, Point& axis2)
//===========================================================================
{
  // Currently applicable only for point groups with an associated plane
  if (!hasSurface())
    return false;
  int code;
  int classtype = getSurface(0)->instanceType(code);
  if (classtype != Class_Plane && classtype != Class_Sphere)
    return false;

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (!elem.get())
    return false;  // Should not happen

  Point norm = elem->direction();
  double ang1 = norm.angle(normalcone_.centre());
  double ang2 = norm.angle(normalcone2_.centre());
  double ang3 = norm.angle(avnorm_);
  double ang4 = norm.angle(avnorm2_);
  double anglim = 0.1*M_PI;
  if (ang1 > anglim && ang2 > anglim && ang3 > anglim && ang4 > anglim)
    return false;   // Surface orientation is incompatible with point cloud

  double min_ang = M_PI;
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      if (adj[ki].first->instanceType() != Class_Cylinder &&
	  adj[ki].first->instanceType() != Class_Cone &&
	  adj[ki].first->instanceType() != Class_Torus)
	continue;
      if (adj[ki].second->getSurfaceFlag() >= ACCURACY_POOR)
	continue;
      Point dir = adj[ki].first->direction();
      double ang = norm.angle(dir);
      ang = std::min(M_PI,ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  axis = dir;
	  axis2 = adj[ki].first->direction2();
	  pos = adj[ki].first->location();
	}
    }

  if (axis.dimension() == 0 || min_ang > anglim)
    return false;
  else
    return true;
}


//===========================================================================
void
RevEngRegion::analyseRotated(Point& pos, Point& axis, Point& axis2)
//===========================================================================
{
  
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, axis2, axis, pos, rotated);

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

#ifdef DEBUG_VALIDATE
  std::ofstream of("rotated_points.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*axis2 << " " << pos+0.5*len*axis2 << std::endl;
#endif

  int stop_break = 1;
}

//===========================================================================
shared_ptr<Torus>
RevEngRegion::analyseRotatedTorus(Point& pos, Point& Cx, Point& axis,
				  double tol, double angtol)
//===========================================================================
{
  // Compute minor circle of torus given rotational axis
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(group_points_.begin(), group_points_.end()));
  RevEngUtils::rotateToPlane(group, Cx, axis, pos, rotated);

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

#ifdef DEBUG_EXTRACT
  std::ofstream of("rotated_points_tor.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << rotated.size() << std::endl;
  for (size_t kr=0; kr<rotated.size(); ++kr)
    of << rotated[kr] << std::endl;
  of << "410 1 0 4 0 0 255 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*axis << " " << pos+0.5*len*axis << std::endl;
  of << "410 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos-0.5*len*Cx << " " << pos+0.5*len*Cx << std::endl;
#endif

  Point cpos;
  double crad;
  Point Cy = axis.cross(Cx);
  RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, axis, cpos, crad);
  shared_ptr<Circle> circ(new Circle(crad, cpos, Cy, Cx));
#ifdef DEBUG_EXTRACT
  circ->writeStandardHeader(of);
  circ->write(of);
#endif
  Point cvec = cpos - pos;
  double R1 = (cvec - (cvec*axis)*axis).length();
  double R2 = (cvec*axis)*axis.length();
  shared_ptr<Torus> tor(new Torus(R1, crad, pos+R2*axis, axis, Cy));
  
#ifdef DEBUG_EXTRACT
  tor->writeStandardHeader(of);
  tor->write(of);
#endif

  // Check accuracy
  double maxd, avd;
  int num, num2;
  vector<RevEngPoint*> in, out;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor, tol, maxd, avd, num, num2, in, out,
			  parvals, dist_ang, angtol);
  int sf_flag = defineSfFlag(0, tol, num, num2, avd, false);
  if (sf_flag == ACCURACY_OK)
    {
      // Finished
      return tor;
    }

  vector<RevEngPoint*> cyl_pts1, cyl_pts2, remaining1, remaining2;
  int sgn = (avH_ < 0.0) ? 1 : -1;
  RevEngUtils::extractLinearPoints(group_points_, rotated,
				   len, pos, axis, R1+sgn*crad, axis, false,
				   tol, angtol, dist_ang,
				   cyl_pts1, true, remaining1);
  RevEngUtils::extractLinearPoints(group_points_, rotated,
				   len, pos, axis, R1+sgn*crad, axis, false,
				   tol, angtol, dist_ang,
				   cyl_pts2, false, remaining2);

#ifdef DEBUG_EXTRACT
  if (cyl_pts1.size() > 0)
    {
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << cyl_pts1.size() << std::endl;
      for (size_t ki=0; ki<cyl_pts1.size(); ++ki)
	of << cyl_pts1[ki]->getPoint() << std::endl;
    }
  if (cyl_pts2.size() > 0)
    {
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << cyl_pts2.size() << std::endl;
      for (size_t ki=0; ki<cyl_pts2.size(); ++ki)
	of << cyl_pts2[ki]->getPoint() << std::endl;
    }
#endif

  vector<RevEngPoint*> distpt, angpt;
  identifyDistPoints(dist_ang, tol, maxd, avd, distpt);
  identifyAngPoints(dist_ang, angtol, tol, angpt);
#ifdef DEBUG_EXTRACT
  std::ofstream of2("dist_ang_out.g2");
  if (distpt.size() > 0)
    {
      of2 << "400 1 0 4 55 200 0 255" << std::endl;
      of2 << distpt.size() << std::endl;
      for (size_t ki=0; ki<distpt.size(); ++ki)
	of2 << distpt[ki]->getPoint() << std::endl;
    }
  if (angpt.size() > 0)
    {
      of2 << "400 1 0 4 55 0 200 255" << std::endl;
      of2 << angpt.size() << std::endl;
      for (size_t ki=0; ki<angpt.size(); ++ki)
	of2 << angpt[ki]->getPoint() << std::endl;
    }
#endif
  int stop_break = 1;
  return tor;
}

//===========================================================================
bool
RevEngRegion::analyseTorusContext(vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
				  double tol, double angtol, vector<size_t>& adjacent_ix,
				  int &plane_ix, int& cyl_ix,
				  Point& pos, Point& axis, Point& Cx, double& R1,
				  double& R2, double dom[4], bool& outer,
				  bool& analyse_rotated)
//===========================================================================
{
  for (size_t ki=0; ki<adj.size(); ++ki)
    for (size_t kj=ki+1; kj<adj.size(); ++kj)
      if (adj[kj].second->numPoints() > adj[ki].second->numPoints())
	std::swap(adj[ki], adj[kj]);

  analyse_rotated = false;
  plane_ix=-1;
  cyl_ix=-1;
  //vector<Point> dir;
  vector<pair<size_t, size_t> > combo;
  vector<pair<int, int> > num_common;
  vector<double> angle;
  for (size_t ki=0; ki<adj.size(); ++ki)
    {
      shared_ptr<ElementarySurface> elem = adj[ki].first;
      if (elem->instanceType() != Class_Plane)
	continue;
      Point plane_axis = elem->direction();
       for (size_t kj=0; kj<adj.size(); ++kj)
	{
	  shared_ptr<ElementarySurface> elem2 = adj[kj].first;
	  if (elem2->instanceType() != Class_Cylinder)
	    continue;
	  Point cyl_axis = elem2->direction();
	  double ang = plane_axis.angle(cyl_axis);
	  ang = std::min(ang, M_PI-ang);
	  if (ang < angtol)
	    {
	      combo.push_back(std::make_pair(ki, kj));
	      vector<RevEngPoint*> adj_pt1 = extractNextToAdjacent(adj[ki].second);
	      vector<RevEngPoint*> adj_pt2 = extractNextToAdjacent(adj[kj].second);
	      num_common.push_back(std::make_pair((int)adj_pt1.size(), (int)adj_pt2.size()));
	      angle.push_back(ang);
	    }
	}
    }
     
  // for (size_t ki=0; ki<adj.size(); ++ki)
  //   {
  //     shared_ptr<ElementarySurface> elem = adj[ki].first;
  //     if (!(elem->instanceType() == Class_Plane || elem->instanceType() == Class_Cylinder))
  // 	continue;

  //     Point curr_dir = elem->direction();
  //     dir.push_back(curr_dir);
  //     adjacent_ix.push_back(ki);
  //     if (elem->instanceType() == Class_Plane)
  // 	{
  // 	  plane = true;
  // 	  plane_ix = (int)ki;
  // 	}
  //     else
  // 	{
  // 	  cyl = true;
  // 	  pos = elem->location();
  // 	  Cx = elem->direction2();
  // 	  radius = elem->radius(0.0, 0.0);
  // 	  cyl_ix = (int)ki;
  // 	}

  //     for (size_t kj=ki+1; kj<adj.size(); ++kj)
  // 	{
  // 	  if (!(adj[kj].first->instanceType() == Class_Plane ||
  // 		adj[kj].first->instanceType() == Class_Cylinder))
  // 	    continue;
  // 	  curr_dir = adj[kj].first->direction();
  // 	  double ang = dir[0].angle(curr_dir);
  // 	  ang = std::min(ang, M_PI-ang);
  // 	  if (ang < angtol)
  // 	    {
  // 	      dir.push_back(curr_dir);
  // 	      adjacent_ix.push_back(kj);
  // 	      if (adj[kj].first->instanceType() == Class_Plane)
  // 		{
  // 		  plane = true;
  // 		  if (plane_ix < 0)
  // 		    plane_ix = (int)kj;
  // 		}
  // 	      else
  // 		{
  // 		  if (!cyl)
  // 		    {
  // 		      pos = adj[kj].first->location();
  // 		      Cx = adj[kj].first->direction2();
  // 		      radius = adj[kj].first->radius(0.0, 0.0);
  // 		      cyl_ix = (int)kj;
  // 		    }
  // 		  cyl = true;
  // 		}
  // 	    }
  // 	}

  //     if (plane && cyl)
  // 	{
  // 	  combo.push_back(plane_ix, cyl_ix);
	  
  // 	break;
  //     else
  // 	{
  // 	  plane = false;
  // 	  cyl = false;
  // 	  dir.clear();
  // 	  adjacent_ix.clear();
  // 	  plane_ix = cyl_ix = -1;
  // 	}
  //   }

  if (combo.size() == 0)
    return false;

  double fac = 1.1;
  int max_num = 0;
  int ix = -1;
  double min_ang = M_PI;
  for (size_t ki=0; ki<combo.size(); ++ki)
    {
      if (num_common[ki].first+num_common[ki].second > max_num /*&& angle[ki] < fac*min_ang*/)
	{
	  max_num = num_common[ki].first+num_common[ki].second;
	  min_ang = angle[ki];
	  ix = (int)ki;
	}
    }

  if (ix < 0)
    return false;
  
  plane_ix = (int)combo[ix].first;
  cyl_ix = (int)combo[ix].second;
  pos = adj[cyl_ix].first->location();
  double radius = adj[cyl_ix].first->radius(0.0, 0.0);

  // if (plane && cyl)
  //   {
  //     // Compute axis
  //     for (size_t ki=1; ki<dir.size(); ++ki)
  // 	if (dir[0]*dir[ki] < 0.0)
  // 	  dir[ki] *= -1;

      // axis = Point(0.0, 0.0, 0.0);
      // double frac = 1.0/(double)dir.size();
      // for (size_t ki=0; ki<dir.size(); ++ki)
      // 	axis += frac*dir[ki];
      outer = true;
      axis = adj[plane_ix].first->direction();
      Point axis2 = adj[cyl_ix].first->direction();

#ifdef DEBUG_TORUSCONTEXT
      std::ofstream of("adj_elem.g2");
      for (size_t ki=0; ki<adjacent_ix.size(); ++ki)
	{
	  adj[adjacent_ix[ki]].first->writeStandardHeader(of);
	  adj[adjacent_ix[ki]].first->write(of);
	  adj[adjacent_ix[ki]].second->writeRegionInfo(of);
	}
#endif
      
      Cx = adj[cyl_ix].first->direction2();
      Point Cy = axis.cross(Cx);
      Cx = Cy.cross(axis);  // Project to plane
      
#ifdef DEBUG_TORUSCONTEXT
      // Check accuracy of alternative plane
      Point loc = adj[cyl_ix].first->location();
      Point mid(0.0, 0.0, 0.0);
      double wgt = 1.0/(double)(adj[plane_ix].second->numPoints());
      for (auto it=adj[plane_ix].second->pointsBegin();
	   it!=adj[plane_ix].second->pointsEnd(); ++it)
	{
	  Vector3D pos0 = (*it)->getPoint();
	  Point xpos(pos0[0], pos0[1], pos0[2]);
	  Point vec = xpos - loc;
	  Point pos2 = loc + (vec*axis2)*axis2;
	  mid += wgt*pos2;
	}
      shared_ptr<Plane> plane(new Plane(mid, axis2));
      vector<RevEngPoint*> inptp, outptp;
      vector<pair<double, double> > dist_angp;
      double maxdp, avdp;
      int num_inp, num2_inp;
      vector<double> parvalsp;
      RevEngUtils::distToSurf(adj[plane_ix].second->pointsBegin(),
			      adj[plane_ix].second->pointsEnd(),
			      plane, tol, maxdp, avdp, num_inp, num2_inp,
			      inptp, outptp, parvalsp, dist_angp, -1.0);

      std::ofstream ofd1("in_out_alt_plane.g2");
      ofd1 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd1 << inptp.size() << std::endl;
      for (size_t kr=0; kr<inptp.size(); ++kr)
	ofd1 << inptp[kr]->getPoint() << std::endl;
      ofd1 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd1 << outptp.size() << std::endl;
      for (size_t kr=0; kr<outptp.size(); ++kr)
	ofd1 << outptp[kr]->getPoint() << std::endl;
      
      // Check accuracy of alternative cylinder
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      points.push_back(std::make_pair(adj[cyl_ix].second->pointsBegin(),
				      adj[cyl_ix].second->pointsEnd()));
      BoundingBox bb = adj[cyl_ix].second->getBbox();
      Point xpos;
      double rad;
      Point low = bb.low();
      Point high = bb.high();
      RevEngUtils::computeCylPosRadius(points, low, high, axis, Cx, Cy,
				       xpos, rad);
      shared_ptr<Cylinder> cyl(new Cylinder(rad, xpos, axis, Cy));
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      double maxd, avd;
      int num_in, num2_in;
      vector<double> parvals;
      RevEngUtils::distToSurf(points[0].first, points[0].second,
			      cyl, tol, maxd, avd, num_in, num2_in,
			      inpt, outpt, parvals, dist_ang, -1.0);

      std::ofstream ofd2("in_out_alt_cyl.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt.size() << std::endl;
      for (size_t kr=0; kr<inpt.size(); ++kr)
	ofd2 << inpt[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt.size() << std::endl;
      for (size_t kr=0; kr<outpt.size(); ++kr)
	ofd2 << outpt[kr]->getPoint() << std::endl;
#endif
      
      // Bound cylinder
      int num_pt = adj[cyl_ix].second->numPoints();
      double dist, ang;
      int ka;
      int num_in_cyl = 0;
      for (ka=0; ka<num_pt; ++ka)
	{
	  Vector2D uv = adj[cyl_ix].second->getPoint(ka)->getPar();
	  adj[cyl_ix].second->getPoint(ka)->getSurfaceDist(dist, ang);
	  if (ang <= angtol)
	    {
	      num_in_cyl++;
	      dom[0] = dom[1] = uv[0];
	      dom[2] = dom[3] = uv[1];
	      break;
	    }
	}
      for (; ka<num_pt; ++ka)
	{
	  Vector2D uv = adj[cyl_ix].second->getPoint(ka)->getPar();
	  adj[cyl_ix].second->getPoint(ka)->getSurfaceDist(dist, ang);
	  if (ang > angtol)
	    continue;
	  num_in_cyl++;
	  dom[0] = std::min(dom[0], uv[0]);
	  dom[1] = std::max(dom[1], uv[0]);
	  dom[2] = std::min(dom[2], uv[1]);
	  dom[3] = std::max(dom[3], uv[1]);
	}
      if (dom[1] - dom[2] > 2*M_PI)
	{
	  dom[0] = 0.0;
	  dom[1] = 2*M_PI;
	}

      if (num_in_cyl < 10)
	return false;
      
#ifdef DEBUG_TORUSCONTEXT
      shared_ptr<ElementarySurface> elem2(adj[cyl_ix].first->clone());
      elem2->setParameterBounds(dom[0], dom[2], dom[1], dom[3]);
      std::ofstream ofcyl("bd_cyl.g2");
      elem2->writeStandardHeader(ofcyl);
      elem2->write(ofcyl);
#endif

      // Compute distance from cylinder to plane
      Point pos1 = adj[cyl_ix].first->point(0.5*(dom[0]+dom[1]), dom[2]);
      Point pos2 = adj[cyl_ix].first->point(0.5*(dom[0]+dom[1]), dom[3]);
      double upar1, upar2, upar3, vpar1, vpar2, vpar3, dist1, dist2, dist3;
      Point close1, close2, close3;
      double eps = 1.0e-6;
      adj[plane_ix].first->closestPoint(pos1, upar1, vpar1, close1, dist1, eps);
      adj[plane_ix].first->closestPoint(pos2, upar2, vpar2, close2, dist2, eps);
      Point pos3 = (dist1 < dist2) ? pos + dom[2]*axis2 : pos + dom[3]*axis2;
      adj[plane_ix].first->closestPoint(pos3, upar3, vpar3, close3, dist3, eps);

      if ((close1-pos1)*axis < 0.0)
	axis *= -1;
      
      if ((close1-pos1)*(close2-pos2) < 0.0)
	{
	  analyse_rotated = true;

	  // Check axis direction
	  int num1 = 0, num2 = 0;
	  for (size_t ki=0; ki<group_points_.size(); ++ki)
	    {
	      Vector3D xyz = group_points_[ki]->getPoint();
	      Point gpt(xyz[0], xyz[1], xyz[2]);
	      if ((gpt-pos)*axis < 0.0)
		num1++;
	      else if ((gpt-pos)*axis > 0.0)
		num2++;
	    }
	  if (num2 < num1)
	    axis *= -1;
	}
      
      R2 = std::min(dist1, dist2);
      R1 = radius - R2;
      pos = close3;

      // Check if the torus in inwards or outwards
      double cdist;
      //RevEngPoint *closest_planar =
      (void)adj[plane_ix].second->closestPoint(pos, cdist);
      if (cdist > 0.99*radius)
	{
	  outer = false;
	  R1 = radius + R2;
	}

#ifdef DEBUG_TORUSCONTEXT
      // Bound plane
      double dom2[4];
      int num_pt2 = adj[plane_ix].second->numPoints();
      Vector2D uv2 = adj[plane_ix].second->getPoint(0)->getPar();
      dom2[0] = dom2[1] = uv2[0];
      dom2[2] = dom2[3] = uv2[1];

      for (int ka=0; ka<num_pt2; ++ka)
	{
	  Vector2D uv2 = adj[plane_ix].second->getPoint(ka)->getPar();
	  dom2[0] = std::min(dom2[0], uv2[0]);
	  dom2[1] = std::max(dom2[1], uv2[0]);
	  dom2[2] = std::min(dom2[2], uv2[1]);
	  dom2[3] = std::max(dom2[3], uv2[1]);
	}
      shared_ptr<ElementarySurface> elem3(adj[plane_ix].first->clone());
      elem3->setParameterBounds(dom2[0], dom2[2], dom2[1], dom2[3]);
      std::ofstream ofpl("bd_plane.g2");
      elem3->writeStandardHeader(ofpl);
      elem3->write(ofpl);
#endif
      
      return true;
  //   }
  
  // return false;  // Requested configuration not found   
}

//===========================================================================
RevEngPoint* RevEngRegion::closestPoint(const Point& pos, double& dist)
//===========================================================================
{
  dist = std::numeric_limits<double>::max();
  int ix = -1;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      double curr_dist = pos.dist(pt);
      if (curr_dist < dist)
	{
	  dist = curr_dist;
	  ix = (int)ki;
	}
    }

  if (ix >= 0)
    return group_points_[ix];
  else
    return 0;
}

//===========================================================================
RevEngPoint* RevEngRegion::closestParPoint(const Point& parpt, double& dist)
//===========================================================================
{
  if (!hasSurface())
    return 0;
  
  dist = std::numeric_limits<double>::max();
  int ix = -1;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector2D uv = group_points_[ki]->getPar();
      Point pt(uv[0], uv[1]);
      double curr_dist = parpt.dist(pt);
      if (curr_dist < dist)
	{
	  dist = curr_dist;
	  ix = (int)ki;
	}
    }

  if (ix >= 0)
    return group_points_[ix];
  else
    return 0;
}

//===========================================================================
bool RevEngRegion::extractTorus(Point mainaxis[3], double tol, int min_pt,
				int min_pt_reg, double angtol,
				int prefer_elementary,
				vector<shared_ptr<HedgeSurface> >& hedgesfs,
				vector<HedgeSurface*>& prevsfs,
				vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  bool found = false;
  int min_nmb = 20;
  if ((int)group_points_.size() < min_nmb)
    return false;
  
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);

  vector<Point> adj_axis;
  axisFromAdjacent(angtol, adj_axis);
  size_t num_adj = adj_axis.size();
  for (int ka=0; ka<3; ++ka)
    {
      size_t kj;
      for (kj=0; kj<num_adj; ++kj)
	{
	  double ang = mainaxis[ka].angle(adj_axis[kj]);
	  ang = std::min(ang, M_PI-ang);
	  if (ang <= angtol)
	    break;
	}
      if (kj == num_adj)
	adj_axis.push_back(mainaxis[ka]);
    }
  shared_ptr<Torus> tor1 = computeTorus(group_points_, adj_axis, tol,
					angtol);
  if (!tor1.get())
    return false;
  if (tor1->getMajorRadius() <= tor1->getMinorRadius())
    return false;
  
  // Check accuracy
  double maxd1, avd1, maxd2=100*tol, avd2=100*tol;
  int num1, num1_2, num2=0, num2_2=0;
  //shared_ptr<ParamSurface> surf = cyl;
  vector<RevEngPoint*> in1, out1,in2, out2;
  vector<pair<double, double> > dist_ang1;
  vector<pair<double, double> > dist_ang2;
  vector<double> parvals1;
  vector<double> parvals2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  tor1, tol, maxd1, avd1, num1, num1_2, in1, out1,
			  parvals1, dist_ang1, angtol);
  int sf_flag1 = defineSfFlag(0, tol, num1, num1_2, avd1, false);
  
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_tor.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << in1.size() << std::endl;
  for (size_t kr=0; kr<in1.size(); ++kr)
    ofd << in1[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << out1.size() << std::endl;
  for (size_t kr=0; kr<out1.size(); ++kr)
    ofd << out1[kr]->getPoint() << std::endl;
  // ofd << "400 1 0 4 155 50 50 255" << std::endl;
  // ofd << in2.size() << std::endl;
  // for (size_t kr=0; kr<in2.size(); ++kr)
  //   ofd << in2[kr]->getPoint() << std::endl;
  // ofd << "400 1 0 4 50 155 50 255" << std::endl;
  // ofd << out2.size() << std::endl;
  // for (size_t kr=0; kr<out2.size(); ++kr)
  //   ofd << out2[kr]->getPoint() << std::endl;
#endif

  shared_ptr<Torus> tor2;
    if (sf_flag1 > ACCURACY_OK)
    {
      Point pos = tor1->location();
      Point axis = tor1->direction();
      Point Cx = tor1->direction2();
      tor2 = analyseRotatedTorus(pos, Cx, axis, tol, angtol);
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      tor2, tol, maxd2, avd2, num2, num2_2, in2, out2,
			      parvals2, dist_ang2, angtol);
    }

  size_t insize = std::max(in1.size(), in2.size());
  if (insize > group_points_.size()/2 && (int)insize > min_nmb)
    {
      shared_ptr<Torus> tor_in2;
      shared_ptr<Torus> tor_in1 =
	computeTorus(in1.size() > in2.size() ? in1 : in2, adj_axis,
		     tol, angtol);
if (tor_in1.get())
{
      vector<RevEngPoint*> inpt_in1, outpt_in1, inpt_in2, outpt_in2; 
      vector<pair<double, double> > dist_ang_in1, dist_ang_in2;
      double maxd_in1, avd_in1, maxd_in2, avd_in2;
      int num1_in, num1_in2, num2_in, num2_in2;
      vector<double> parvals_in1;
      vector<double> parvals_in2;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      tor_in1, tol, maxd_in1, avd_in1, num1_in, num1_in2,
			      inpt_in1, outpt_in1, parvals_in1, dist_ang_in1,
			      angtol);
  
      // RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
      // 			      tor_in2, tol, maxd_in2, avd_in2, num2_in, num2_in2,
      // 			      inpt_in2, outpt_in2, parvals_in2,
      // 			      dist_ang_in2, angtol);
  
#ifdef DEBUG_EXTRACT
      std::ofstream ofd2("in_out_tor2.g2");
      ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      ofd2 << inpt_in1.size() << std::endl;
      for (size_t kr=0; kr<inpt_in1.size(); ++kr)
	ofd2 << inpt_in1[kr]->getPoint() << std::endl;
      ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      ofd2 << outpt_in1.size() << std::endl;
      for (size_t kr=0; kr<outpt_in1.size(); ++kr)
	ofd2 << outpt_in1[kr]->getPoint() << std::endl;
      // ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
      // ofd2 << inpt_in2.size() << std::endl;
      // for (size_t kr=0; kr<inpt_in2.size(); ++kr)
      // 	ofd2 << inpt_in2[kr]->getPoint() << std::endl;
      // ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
      // ofd2 << outpt_in2.size() << std::endl;
      // for (size_t kr=0; kr<outpt_in2.size(); ++kr)
      // 	ofd2 << outpt_in2[kr]->getPoint() << std::endl;
#endif
      
      if (num1_in > num1 || avd_in1 < avd1)
	{
	  std::swap(tor1, tor_in1);
	  std::swap(num1, num1_in);
	  std::swap(num1_2, num1_in2);
	  std::swap(avd1, avd_in1);
	  std::swap(maxd1, maxd_in1);
	  std::swap(parvals1, parvals_in1);
	  std::swap(dist_ang1, dist_ang_in1);
	  // std::cout << "Torus swap 1" << std::endl;
	  // std::cout << group_points_.size() << " " << num1_in << " ";
	  // std::cout << num1 << " " << avd_in1 << " " << avd1;
	  // std::cout << " " << maxd_in1 << " " << maxd1 << std::endl;
	}
      // if (num2_in > num2 || avd_in2 < avd2)
      // 	{
      // 	  std::swap(tor2, tor_in2);
      // 	  std::swap(num2, num2_in);
      // 	  std::swap(num2_2, num2_in2);
      // 	  std::swap(avd2, avd_in2);
      // 	  std::swap(maxd2, maxd_in2);
      // 	  std::swap(parvals2, parvals_in2);
      // 	  std::swap(dist_ang2, dist_ang_in2);
      // 	  // std::cout << "Torus swap 2" << std::endl;
      // 	  // std::cout << group_points_.size() << " " << num2_in << " ";
      // 	  // std::cout << num2 << " " << avd_in2 << " " << avd2;
      // 	  // std::cout << " " << maxd_in2 << " " << maxd2 << std::endl;
      // 	}
    }
    }
  
//int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  int sf_flag2 = defineSfFlag(0, tol, num2, num2_2, avd2, false);
  if (tor1->getMajorRadius() > tor1->getMinorRadius() &&
      sf_flag1 < ACCURACY_POOR && num1 > num2)
    {
      bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag1 > surfflag_)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num1 < acc_fac1*num_init ||
		    avd1 > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num1 < num_init || avd1 > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals1[2*kh],parvals1[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang1[kh].first, dist_ang1[kh].second);
	    }
	  setAccuracy(maxd1, avd1, num1, num1_2);      
      
#ifdef DEBUG
	  std::cout << "Torus 1. N1: " << num << ", N2: " << num1 << ", max: " << maxd1 << ", av: " << avd1 << std::endl;
#endif
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  //shared_ptr<ParamSurface> tor1_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor1, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag1);
	}
	// }
	}
  else if (tor2->getMajorRadius() > tor2->getMinorRadius() &&
	   sf_flag2 < ACCURACY_POOR && num2 >= num1) 
    {
      bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag2 > surfflag_)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd2 > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num2 < num_init || avd2 > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals2[2*kh],parvals2[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang2[kh].first, dist_ang2[kh].second);
	    }
	  setAccuracy(maxd2, avd2, num2, num2_2);
      
#ifdef DEBUG
	  std::cout << "Torus 2. N1: " << num << ", N2: " << num2 << ", max: " << maxd2 << ", av: " << avd2 << std::endl;
#endif
	  // associated_sf_.clear();
	  // for (int ka=0; ka<divmod[0]->nmbEntities(); ++ka)
	  // 	{
	  // 	  shared_ptr<ParamSurface> tor2_2 = divmod[0]->getSurface(ka);
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(tor2, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag2);
	  // }
	}
    }
  if (!basesf_.get() ||
      (num1 >= num_in_base_ && avd1 < avdist_base_))
    setBaseSf(tor1, maxd1, avd1, num1, num1_2);
  if (!basesf_.get() ||
      (num2 >= num_in_base_ && avd2 < avdist_base_))
    setBaseSf(tor2, maxd2, avd2, num2, num2_2);

  return found;
}
 
//===========================================================================
  shared_ptr<Torus> RevEngRegion::computeTorus(vector<RevEngPoint*>& points,
					       vector<Point>& adj_axis,
					       double tol, double angtol)
//===========================================================================
{
#ifdef DEBUG_EXTRACT
  std::ofstream of("torus_compute.g2");
#endif
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)points.size();
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      double kmax = points[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(points.size());
  Point mid(0.0, 0.0, 0.0);
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      double kmax = points[kr]->maxPrincipalCurvature();
      k2mean += wgt*kmax;

      Vector3D xyz = points[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = points[kr]->getLocFuncNormal();
      centr[kr] = xyz2 + rd*norm;
      mid += wgt*centr[kr];
    }
  
#ifdef DEBUG_EXTRACT
  of << "400 1 0 4 155 100 0 255" << std::endl;
  of << points.size() << std::endl;
  for (size_t kr=0; kr<points.size(); ++kr)
    {
      of << centr[kr] << std::endl;
    }
#endif

  Point low = bbox_.low();
  Point high = bbox_.high();
  double len = low.dist(high);

  ImplicitApprox impl;
  impl.approxPoints(centr, 1);

  double val;
  Point grad;
  impl.evaluate(mid, val, grad);
  grad.normalize_checked();
  
  shared_ptr<Torus> dummy;
  Point pos, normal;
  bool found = impl.projectPoint(mid, grad, pos, normal);
  if (!found)
    return dummy;
  double eps1 = 1.0e-8;
  if (normal.length() < eps1)
      return dummy;

  for (size_t kr=0; kr<adj_axis.size(); ++kr)
    {
      double ang = adj_axis[kr].angle(normal);
      ang = std::min(ang, M_PI-ang);
      if (ang < angtol)
	{
	  normal = adj_axis[kr];
	  break;
	}
    }
  
#ifdef DEBUG_EXTRACT
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos << std::endl;
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos << " " << pos+0.5*len*normal << std::endl;
#endif

  // std::ofstream ofimpl("impl_plane.g2");
  // impl->visualize(points, ofimpl);
  
  Vector3D xyz = points[0]->getPoint();
  Point xyz2(xyz[0], xyz[1], xyz[2]);
  Point Cx = centr[0] - xyz2;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  try {
    RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  }
  catch (...)
    {
      shared_ptr<Torus> dummy;
      return dummy;
    }
  pnt -= ((pnt - pos)*normal)*normal;
  shared_ptr<Circle> circ(new Circle(rad, pnt, normal, Cx));
#ifdef DEBUG_EXTRACT
  circ->writeStandardHeader(of);
  circ->write(of);
#endif
  
//   vector<Point> rotated;
//   vector<pair<vector<RevEngPoint*>::iterator,
// 	      vector<RevEngPoint*>::iterator> > group;
//   group.push_back(std::make_pair(points.begin(), points.end()));
//   RevEngUtils::rotateToPlane(group, Cx, normal, pnt, rotated);
// #ifdef DEBUG_EXTRACT
//   std::ofstream of3("rotated_pts_tor.g2");
//   of3 << "400 1 0 4 255 0 0 255" << std::endl;
//   of3 << rotated.size() << std::endl;
//   for (size_t kr=0; kr<rotated.size(); ++kr)
//     of3 << rotated[kr] << std::endl;
//   of3 << "410 1 0 4 0 0 255 255" << std::endl;
//   of3 << "1" << std::endl;
//   of3 << pnt-0.5*len*normal << " " << pnt+0.5*len*normal << std::endl;
//   of3 << "410 1 0 4 0 255 0 255" << std::endl;
//   of3 << "1" << std::endl;
//   of3 << pnt-0.5*len*Cx << " " << pnt+0.5*len*Cx << std::endl;
// #endif 
//   Point cpos;
//   double crad;
//   RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, normal, cpos, crad);
//   shared_ptr<Circle> circ2(new Circle(crad, cpos, Cy, Cx));
// #ifdef DEBUG_EXTRACT
//   circ2->writeStandardHeader(of3);
//   circ2->write(of3);
// #endif
//   shared_ptr<SplineCurve> spl;
//   Point xpos;
//   vector<double> param;
//   try {
//     curveApprox(rotated, tol, circ2, param, spl, xpos);
//   }
//   catch (...)
//     {
//     }
// #ifdef DEBUG_EXTRACT
// if (spl.get())
//     {
//       spl->writeStandardHeader(of3);
//       spl->write(of3);
//     }
//  #endif
  
//   vector<Point> projected;
//   double maxdp, avdp;
//   RevEngUtils::projectToPlane(group, normal, pnt, projected, maxdp, avdp);

// #ifdef DEBUG_EXTRACT
//   std::ofstream ofp3("projected_pts_tor.g2");
//   ofp3 << "400 1 0 4 255 0 0 255" << std::endl;
//   ofp3 << projected.size() << std::endl;
//   for (size_t kr=0; kr<projected.size(); ++kr)
//     ofp3 << projected[kr] << std::endl;
//   circ->writeStandardHeader(ofp3);
//   circ->write(ofp3);
// #endif
  shared_ptr<Torus> tor1(new Torus(rad, fabs(rd), pnt, normal, Cy));

#ifdef DEBUG_EXTRACT
  tor1->writeStandardHeader(of);
  tor1->write(of);
#endif

 //  Point cvec = cpos - pnt;
//   double R1 = (cvec - (cvec*normal)*normal).length();
//   double R2 = (cvec*normal)*normal.length();
//   shared_ptr<Torus> tor2(new Torus(R1, crad, pnt+R2*normal, normal, Cy));
// #ifdef DEBUG_EXTRACT
//   tor2->writeStandardHeader(of);
//   tor2->write(of);
//   //std::cout << "Torus small radius: " << fabs(rd) << ", " << crad << std::endl;
// #endif
  
//   torus2 = tor2;
  return tor1;
}

//===========================================================================
bool RevEngRegion::tryOtherSurf(int prefer_elementary, bool replace)
//===========================================================================
{
  if (associated_sf_.size() > 0 && (!replace))
    return false;
  
  if (prefer_elementary == BEST_ACCURACY ||
      associated_sf_.size() == 0)
    return true;

  if (prefer_elementary != BEST_ACCURACY && surfflag_ == ACCURACY_OK)
    return false;

  if (prefer_elementary != BEST_ACCURACY && surf_adaption_ > INITIAL)
    return false;

  if (surfflag_ == ACCURACY_POOR)
    return true;

  if (prefer_elementary == ALWAYS_ELEM)
    {
      int sfcode;
      ClassType type = associated_sf_[0]->instanceType(sfcode);
      double fac = 5.0;
      if ((type == Class_Plane || type == Class_Sphere) && MAH_ > fac*MAK_)
	return true;
      else
	return false;
    }
  
  // Check if the current accuracy is sufficient
  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  if ((double)num_init > 0.75*num)
    return false;
  else
    return true;
  
}


//===========================================================================
bool RevEngRegion::extractFreeform(double tol, int min_pt, int min_pt_reg,
				   double angtol, int prefer_elementary,
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<HedgeSurface*>& prevsfs,
				   vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // std::ofstream of("curr_region.g2");
  // writeRegionInfo(of);
  
  // std::ofstream of2("curr_normals.g2");
  // writeUnitSphereInfo(of2);

  bool found = false;
  int min_nmb = 20;
  //double eps = 1.0e-6;
  if ((int)group_points_.size() < min_nmb)
    return false;

  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*> > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);
  
  shared_ptr<SplineSurface> spl = computeFreeform(group_points_, tol);
  if (!spl.get())
    return false;
#ifdef DEBUG_EXTRACT
  std::ofstream ofs("spl.g2");
  spl->writeStandardHeader(ofs);
  spl->write(ofs);
#endif
  
  // Check accuracy
  double maxd, avd; 
  int num2, num2_2; 
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
  			  spl, tol, maxd, avd, num2, num2_2, inpt, outpt,
			  parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
  std::ofstream ofd("in_out_spl.g2");
  ofd << "400 1 0 4 155 50 50 255" << std::endl;
  ofd << inpt.size() << std::endl;
  for (size_t kr=0; kr<inpt.size(); ++kr)
    ofd << inpt[kr]->getPoint() << std::endl;
  ofd << "400 1 0 4 50 155 50 255" << std::endl;
  ofd << outpt.size() << std::endl;
  for (size_t kr=0; kr<outpt.size(); ++kr)
    ofd << outpt[kr]->getPoint() << std::endl;
#endif

  int num = (int)group_points_.size();
  double maxd_init, avd_init;
  int num_init, num_init2;
  getAccuracy(maxd_init, avd_init, num_init, num_init2);
  if (associated_sf_.size() > 0 ||  (num2 > min_pt && num2 > num/2))
    {
      extractOutPoints(dist_ang, tol, angtol, out_groups);
      if (out_groups.size() > 0)
	{
	  // Some points has been removed from the group. Redo surface
	  // generation
	  spl = computeFreeform(group_points_, tol);
	  if (!spl.get())
	    return false;
#ifdef DEBUG_EXTRACT
	  std::ofstream ofs2("spl2.g2");
	  spl->writeStandardHeader(ofs2);
	  spl->write(ofs2);
#endif
	  
	  num = (int)group_points_.size();
	  dist_ang.clear();
	  parvals.clear();
	  inpt.clear();
	  outpt.clear();
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  spl, tol, maxd, avd, num2, num2_2, inpt, outpt,
				  parvals, dist_ang, angtol);
#ifdef DEBUG_EXTRACT
	  std::ofstream ofd2("in_out_spl2.g2");
	  ofd2 << "400 1 0 4 155 50 50 255" << std::endl;
	  ofd2 << inpt.size() << std::endl;
	  for (size_t kr=0; kr<inpt.size(); ++kr)
	    ofd2 << inpt[kr]->getPoint() << std::endl;
	  ofd2 << "400 1 0 4 50 155 50 255" << std::endl;
	  ofd2 << outpt.size() << std::endl;
	  for (size_t kr=0; kr<outpt.size(); ++kr)
	    ofd2 << outpt[kr]->getPoint() << std::endl;
#endif
	}
      int stop_break_out = 1;
    }

  int sf_flag = defineSfFlag(0, tol, num2, num2_2, avd, false);
  if (sf_flag < ACCURACY_POOR)
    {
     bool OK = true;
      if (associated_sf_.size() > 0)
	{
	  // Check with current approximating surface
	  double acc_fac1 = 1.25;
	  double acc_fac2 = 0.75;
	  if (surfflag_ == ACCURACY_OK && sf_flag > surfflag_)
	    OK = false;
	  else if (prefer_elementary == ALWAYS_ELEM)
	    OK = false;
	  else if (prefer_elementary == PREFER_ELEM &&
		   ((double)num2 < acc_fac1*num_init ||
		    avd > acc_fac2*avd_init))
	    OK = false;
	  else if (prefer_elementary == BEST_ACCURACY &&
		   (num2 < num_init || avd > avd_init))
	    OK = false;
	}

      if (OK)
	{
	  found = true;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
	    }
	  setAccuracy(maxd, avd, num2, num2_2);
#ifdef DEBUG
	  std::cout << "Spline. N1: " << num << ", N2: " << num2 << ", max: " << maxd << ", av: " << avd << std::endl;
#endif
	  shared_ptr<HedgeSurface> hedge(new HedgeSurface(spl, this));
	  for (size_t kh=0; kh<associated_sf_.size(); ++kh)
	    prevsfs.push_back(associated_sf_[kh]);
	  setHedge(hedge.get());
	  hedgesfs.push_back(hedge);
	  setSurfaceFlag(sf_flag);
	}
    }
  int stop_break0 = 1;
  return found;
}


//===========================================================================
 void RevEngRegion::removeAndUpdatePoints(vector<RevEngPoint*>& points)
//===========================================================================
 {
   for (size_t ki=0; ki<points.size(); ++ki)
     {
       points[ki]->unsetRegion();
       points[ki]->unsetSurfInfo();
       vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
						     group_points_.end(),
						     points[ki]);
       if (it != group_points_.end())
	 {
	   std::swap(*it, group_points_[group_points_.size()-1]);
	   group_points_.pop_back();
	     //group_points_.erase(it);
	 }
      }
 }

//===========================================================================
 void RevEngRegion::extractOutOfEdge(shared_ptr<CurveOnSurface>& cv,
				     vector<shared_ptr<CurveOnSurface> >& intcv,
				     double radius, double tol, double angtol,
				     vector<RevEngPoint*>& out_points)
//===========================================================================
 {
   if (!hasSurface())
     return; // Points are not parameterized
   shared_ptr<ParamCurve> pcv = cv->parameterCurve();
   if (!pcv.get())
     return;  // Not possible to check configuration in parameter domain

   Point dir, norm0, norm;
   if (pcv->isLinear(dir, angtol))
       norm0 = Point(-dir[1],dir[0]);

   // Check orientation of trimming curve
   Point midp(0.5*(domain_[0]+domain_[1]), 0.5*(domain_[2]+domain_[3]));
   double tpar1, dist1;
   Point pclose1;
   pcv->closestPoint(midp, pcv->startparam(), pcv->endparam(), tpar1, pclose1, dist1);
   Point vec = pclose1 - midp;
   if (norm0.dimension() > 0)
     norm = norm0;
   else
     {
       vector<Point> der(2);
       pcv->point(der, tpar1, 1);
       norm = Point(-der[1][1], der[1][0]);
     }
   int sgn = 0;
   double dist2 = std::numeric_limits<double>::max();
   for (size_t ki=0; ki<intcv.size(); ++ki)
     {
       shared_ptr<ParamCurve> pcv2 = intcv[ki]->parameterCurve();
       if (!pcv2.get())
	 continue;
       double tpar2, dist2_2;
       Point pclose2;
       pcv2->closestPoint(midp, pcv2->startparam(), pcv2->endparam(), tpar2, pclose2, dist2_2);
       if (dist2_2 < dist2)
	 {
	   dist2 = dist2_2;
	   if (dist2 < dist1)
	     sgn = (norm*vec < 0.0) ? 1 : -1;
	   else
	     sgn = (norm*vec < 0.0) ? -1 : 1;
	 }
     }
   
   if (sgn == 0)
     sgn = (norm*vec < 0.0) ? -1 : 1;

   // Check with point distribution
   double mdist, cdist;
   RevEngPoint *mclose = closestParPoint(midp, mdist);
   RevEngPoint *cclose = closestParPoint(pclose1, cdist);

   vector<RevEngPoint*> remaining;
   for (size_t ki=0; ki<group_points_.size(); ++ki)
     {
       Vector2D uv = group_points_[ki]->getPar();
       Point par(uv[0], uv[1]);

       // Check if the parameter point lies to the left of the parameter curve
       double tpar, dist;
       Point close;
       pcv->closestPoint(par, pcv->startparam(), pcv->endparam(), tpar,
			 close, dist);

       if (norm0.dimension() > 0)
	 norm = norm0;
       else
	 {
	   vector<Point> der(2);
	   pcv->point(der, tpar, 1);
	   norm = Point(-der[1][1], der[1][0]);
	 }

       vec = sgn*(close - par);
       double ang = norm.angle(vec);
       ang = std::min(ang, M_PI-ang);
       if (norm*vec < 0.0 && ang <= angtol)
	 {
	   Vector3D xyz = group_points_[ki]->getPoint();
	   Point pos(xyz[0], xyz[1], xyz[2]);
	   double dist3 = std::numeric_limits<double>::max();
	   for (size_t kj=0; kj<intcv.size(); ++kj)
	     {
	       shared_ptr<ParamCurve> space = intcv[kj]->spaceCurve();
	       if (!space.get())
		 continue;
	       double tpar3, dist3_2;
	       Point close3;
	       space->closestPoint(pos, space->startparam(), space->endparam(),
				   tpar3, close3, dist3_2);
	       dist3 = std::min(dist3, dist3_2);
	     }

	   if (fabs(dist3-radius) < tol)
	     out_points.push_back(group_points_[ki]);
	   else
	     remaining.push_back(group_points_[ki]);
	 }
       else
	 remaining.push_back(group_points_[ki]);
     }

   if (out_points.size() > 0 && remaining.size() > 0)
     {
       // If all points are outside, none will be extracted
       std::swap(group_points_, remaining);
       updateInfo(tol, angtol);
       
       shared_ptr<ParamSurface> surf = getSurface(0)->surface();
       bool cyllike = (surf->instanceType() == Class_Cylinder ||
		       surf->instanceType() == Class_Cone);
       int sf_flag = defineSfFlag(0, tol, num_inside_, num_inside2_, avdist_, cyllike);
       setSurfaceFlag(sf_flag);
     }
 }

//===========================================================================
 void RevEngRegion::extractOutOfEdge2(vector<shared_ptr<CurveOnSurface> >& intcv,
				      double tol, double angtol,
				     vector<RevEngPoint*>& out_points)
//===========================================================================
 {
   if (!hasSurface())
     return; // Points are not parameterized
   if (intcv.size() == 0)
     return;

   // Check orientation of trimming curve
   double midt = 0.5*(intcv[0]->startparam() + intcv[intcv.size()-1]->endparam());
   BoundingBox pbox(3);
   int ix = -1;
   for (size_t ki=0; ki<intcv.size(); ++ki)
     {
       if (intcv[ki]->startparam() <= midt && intcv[ki]->endparam() >= midt)
	 {
	   ix = (int)ki;
	   break;
	 }
     }
   if (ix < 0)
     return;

   vector<Point> der(2);
   intcv[ix]->point(der, midt, 1);
   Point sfpar = intcv[ix]->faceParameter(midt);
   Point sfnorm;
   shared_ptr<ParamSurface> surf = intcv[ix]->underlyingSurface();
   surf->normal(sfnorm, sfpar[0], sfpar[1]);
   Point vec = der[1].cross(sfnorm);
   vec.normalize();
   double diag = bbox_.low().dist(bbox_.high());
   double dfac = 0.1;
   Point pt1 = der[0] + dfac*diag*vec;
   Point pt2 = der[0] - dfac*diag*vec;
   double mdist1, mdist2;
   RevEngPoint *rpt1, *rpt2;
   rpt1 = closestPoint(pt1, mdist1);
   rpt2 = closestPoint(pt2, mdist2);
   int sgn = (mdist1 <= mdist2) ? -1 : 1;

   vector<RevEngPoint*> remaining;
   for (size_t ki=0; ki<group_points_.size(); ++ki)
     {
       Vector3D xyz = group_points_[ki]->getPoint();
       Point pos(xyz[0], xyz[1], xyz[2]);

       // Check if the point lies to the left of the curve
       double tpar;
       Point close, norm;
       double dist = std::numeric_limits<double>::max();
       int ix2 = -1;
       for (size_t ki=0; ki<intcv.size(); ++ki)
	 {
	   double tpar2, dist2;
	   Point close2;
	   intcv[ki]->closestPoint(pos, intcv[ki]->startparam(), intcv[ki]->endparam(), tpar2, close2, dist2);
	   if (dist2 < dist)
	     {
	       tpar = tpar2;
	       dist = dist2;
	       close = close2;
	       ix2 = (int)ki;
	     }
	 }

       if (ix2 < 0)
	 remaining.push_back(group_points_[ki]);
       else
	 {
	   vector<Point> der1(2);
	   intcv[ix2]->point(der1, tpar, 1);
	   Point sfpar2 = intcv[ix2]->faceParameter(tpar);
	   Point sfnorm2;
	   surf->normal(sfnorm2, sfpar2[0], sfpar2[1]);
	   Point vec2 = der[1].cross(sfnorm);
	   vector<Point> der2(2);
	   intcv[ix]->point(der2, tpar, 1);
	   Point vec3 = sgn*(der2[0] - pos);
	   double ang = der2[1].angle(vec3);
	   ang = fabs(0.5*M_PI - ang);
	   
	   if (vec2*vec3 < 0.0 &&  ang <= angtol)
	     out_points.push_back(group_points_[ki]);
	   else
	     remaining.push_back(group_points_[ki]);
	 }
     }

   if (out_points.size() > 0 && remaining.size() > 0)
     {
       // If all points are outside, none will be extracted
       std::swap(group_points_, remaining);
       updateInfo(tol, angtol);
       
       shared_ptr<ParamSurface> surf = getSurface(0)->surface();
       bool cyllike = (surf->instanceType() == Class_Cylinder ||
		       surf->instanceType() == Class_Cone);
       int sf_flag = defineSfFlag(0, tol, num_inside_, num_inside2_, avdist_, cyllike);
       setSurfaceFlag(sf_flag);
     }
 }

 //===========================================================================
 void RevEngRegion::extractOutPoints(vector<pair<double, double> >& dist_ang,
				     double tol, double angtol,
				     vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
 {
   vector<RevEngPoint*> move;
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].first > tol &&
	   (associated_sf_.size() == 0 ||
	    group_points_[kr]->getSurfaceDist() > tol))
	 move.push_back(group_points_[kr]);
     }
   extractSpesPoints(move, out_groups, true);
   if (out_groups.size() > 0)
     updateInfo(tol, angtol);
 }

//===========================================================================
 void RevEngRegion::identifyAngPoints(vector<pair<double, double> >& dist_ang,
				      double tol, double disttol, 
				      vector<RevEngPoint*>& ang_points)
//===========================================================================
 {
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].second > tol && dist_ang[kr].first > disttol)
	 ang_points.push_back(group_points_[kr]);
     }
 }
//===========================================================================
 void RevEngRegion::identifyAngPoints(vector<pair<double, double> >& dist_ang,
				      double tol, double dtol, double dtol2,
				      vector<RevEngPoint*>& ang_points,
				      vector<RevEngPoint*>& remaining)
//===========================================================================
 {
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if ((dist_ang[kr].second > tol && dist_ang[kr].first > dtol) ||
	   dist_ang[kr].first > dtol2)
	 ang_points.push_back(group_points_[kr]);
       else
	 remaining.push_back(group_points_[kr]);
     }
 }

//===========================================================================
 void RevEngRegion::identifyDistPoints(vector<pair<double, double> >& dist_ang,
				       double tol, double maxd, double avd,
				       vector<RevEngPoint*>& dist_points)
//===========================================================================
 {
   double tol2 = (avd > 0.9*tol) ? 3.0*tol : 2.0*tol;
   for (size_t kr=0; kr<dist_ang.size(); ++kr)
     {
       if (dist_ang[kr].first > tol2)
	 dist_points.push_back(group_points_[kr]);
     }
 }

//===========================================================================
bool RevEngRegion::isConnected()
//===========================================================================
{
  if (group_points_.size() == 0)
    return false;
  vector<RevEngPoint*> sub_group;
  group_points_[0]->fetchConnected(this, (int)group_points_.size(), sub_group);

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    group_points_[ki]->unsetVisited();
#ifdef DEBUG_CHECK
  if (sub_group.size() < group_points_.size())
    {
      vector<vector<RevEngPoint*> > conn_groups;
      std::vector<RevEngPoint*> dummy;
      connectedGroups(group_points_, conn_groups, false, dummy);
      std::ofstream of("disconnect_groups.g2");
      for (size_t ki=0; ki<conn_groups.size(); ++ki)
	{
	  of << "400 1 0 4 0 155 100 255" << std::endl;
	  of << conn_groups[ki].size() << std::endl;
	  for (size_t kj=0; kj<conn_groups[ki].size(); ++kj)
	    of << conn_groups[ki][kj]->getPoint() << std::endl;
	}
    }
#endif

  return (group_points_.size() == sub_group.size());
}

//===========================================================================
 void RevEngRegion::connectedGroups(vector<RevEngPoint*>& move,
				    vector<vector<RevEngPoint*> >& conn_groups,
				    bool outer, vector<RevEngPoint*>& inner)
//===========================================================================
 {
   // Note: derived information in this group is not updated!!
   shared_ptr<RevEngRegion> dummy_reg(new RevEngRegion(edge_class_type_));
      
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->setRegion(dummy_reg.get());  // Preliminary

   vector<vector<RevEngPoint*> > sep_move;
   for (size_t kr=0; kr<move.size(); ++kr)
     {
       if (move[kr]->visited())
	 continue;
       vector<RevEngPoint*> curr_move;
       move[kr]->fetchConnected(dummy_reg.get(), (int)move.size(), curr_move);
       sep_move.push_back(curr_move);
     }
      
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->unsetVisited();
      
#ifdef DEBUG_EXTRACT
   std::ofstream m1("move_group.g2");
   std::ofstream m2("not_move_group.g2");
#endif
   for (size_t kr=0; kr<sep_move.size(); ++kr)
     {
       size_t kh;
       for (kh=0; kh<sep_move[kr].size(); ++kh)
	 {
	   vector<ftSamplePoint*> next = sep_move[kr][kh]->getNeighbours();
	   size_t kh1;
	   for (kh1=0; kh1<next.size(); ++kh1)
	     {
	       RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kh1]);
		  
	       RevEngRegion *reg = pt->region();
	       if (reg != this && reg != dummy_reg.get())
		 break;
	     }
	   if (kh1 < next.size())
	     break;
	 }
       if (kh < sep_move[kr].size() || (!outer))
	 {
	   conn_groups.push_back(sep_move[kr]);
#ifdef DEBUG_EXTRACT
	   m1 << "400 1 0 0" << std::endl;
	   m1 << sep_move[kr].size() << std::endl;
	   for (size_t kh1=0; kh1<sep_move[kr].size(); ++kh1)
	     m1 << sep_move[kr][kh1]->getPoint() << std::endl;
#endif
	 }
       else
	 {
#ifdef DEBUG_EXTRACT
	   m2 << "400 1 0 0" << std::endl;
	   m2 << sep_move[kr].size() << std::endl;
	   for (size_t kh1=0; kh1<sep_move[kr].size(); ++kh1)
	     m2 << sep_move[kr][kh1]->getPoint() << std::endl;
#endif
	   inner.insert(inner.end(), sep_move[kr].begin(), sep_move[kr].end());
	 }
     }
   for (size_t kr=0; kr<move.size(); ++kr)
     move[kr]->setRegion(this);
 }

//===========================================================================
 void RevEngRegion::extractSpesPoints(vector<RevEngPoint*>& move,
				      vector<vector<RevEngPoint*> >& out_groups,
				      bool outer)
//===========================================================================
 {
   // Identify connected groups
   std::vector<RevEngPoint*> dummy;
   connectedGroups(move, out_groups, outer, dummy);
   
   // Extract group of out points from current group
   for (size_t ki=0; ki<out_groups.size(); ++ki)
     for (size_t kj=0; kj<out_groups[ki].size(); ++kj)
       {
	 out_groups[ki][kj]->unsetRegion();
	 out_groups[ki][kj]->addMove();
	 vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
						 group_points_.end(),
						 out_groups[ki][kj]);
	 if (it != group_points_.end())
	   {
	     std::swap(*it, group_points_[group_points_.size()-1]);
	     group_points_.pop_back();
	     //group_points_.erase(it);
	   }
       }

   int stop_break = 1;
 }
 
//===========================================================================
void RevEngRegion::removePoints(vector<RevEngPoint*>& remove)
//===========================================================================
 {
   // Identify other points
   std::vector<RevEngPoint*>::iterator end = group_points_.end();
   std::vector<RevEngPoint*>::iterator last = group_points_.end();
   --last;
   for (size_t ki=0; ki<remove.size(); ++ki)
     {
       auto it = std::find(group_points_.begin(), end, remove[ki]);
       if (it != end)
	 {
	   std::swap((*it), (*last));
	   --last;
	   --end;
	 }
     }
   group_points_.erase(end, group_points_.end());
 }

//===========================================================================
 void RevEngRegion::removeOtherPoints(vector<RevEngPoint*>& keep,
				      vector<HedgeSurface*>& prevsfs,
				      vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
 {
   // Identify other points
   vector<RevEngPoint*> remove;
   for (size_t ki=0; ki<group_points_.size(); ++ki)
     {
       auto it = std::find(keep.begin(), keep.end(), group_points_[ki]);
       if (it == keep.end())
	 remove.push_back(group_points_[ki]);
     }

   int num = (int)group_points_.size();
   if (remove.size() > 0)
     {
       extractSpesPoints(remove, out_groups, false);
       updateInfo();
     }

   if (hasSurface() && (int)group_points_.size() < std::max(20, num/10))
     {
       prevsfs.insert(prevsfs.end(), associated_sf_.begin(), associated_sf_.end());
       clearSurface();
     }

   for (size_t ki=0; ki<out_groups.size(); ++ki)
     for (size_t kj=0; kj<out_groups[ki].size(); ++kj)
       out_groups[ki][kj]->unsetRegion();
   
   updateRegionAdjacency();
   int stop_break = 1;
 }
 
//===========================================================================
shared_ptr<CurveOnSurface>
RevEngRegion::constParSfCv(shared_ptr<ParamSurface> surf, int dir,
			   double par, int bd, double t1, double t2)
//===========================================================================
{
  double eps = std::max(0.0001*(t2-t1), 1.0e-4);
  Point pdir = (dir == 0) ? Point(0.0, 1.0) : Point(1.0, 0.0);
  vector<shared_ptr<ParamCurve> > cvs1 = surf->constParamCurves(par, (dir==1));
  if (!cvs1[0]->isBounded())
    {
      shared_ptr<ElementaryCurve> elemcv =
	dynamic_pointer_cast<ElementaryCurve,ParamCurve>(cvs1[0]);
      elemcv->setParamBounds(t1, t2);
    }

  if (cvs1[0]->startparam() < t1-eps || cvs1[0]->endparam() > t2+eps)
    {
      shared_ptr<ParamCurve> sub(cvs1[0]->subCurve(t1, t2));
      cvs1[0] = sub;
    }
  
  Point ppos = (dir == 0) ? Point(par, 0.0) : Point(0.0, par);
  shared_ptr<ElementaryCurve> pcrv(new Line(ppos, pdir));
  pcrv->setParamBounds(cvs1[0]->startparam(), cvs1[0]->endparam());
#ifdef DEBUG_TRIM
  Point sp1, sp2, ssp1, ssp2;
  Point pp1, pp2;
  cvs1[0]->point(sp1, cvs1[0]->startparam());
  cvs1[0]->point(sp2, cvs1[0]->endparam());
  pcrv->point(pp1, pcrv->startparam());
  pcrv->point(pp2, pcrv->endparam());
  surf->point(ssp1, pp1[0], pp1[1]);
  surf->point(ssp2, pp2[0], pp2[1]);
#endif
  shared_ptr<CurveOnSurface> sfcv(new CurveOnSurface(surf, pcrv, cvs1[0], false,
						     3, dir+1, par, bd, true));
  return sfcv;
 }
//===========================================================================
bool RevEngRegion::trimSurface(double tol)
//===========================================================================
{
  if (!hasSurface())
    return false;

  if (trim_edgs_.size() == 0)
    return false;   // Not ready for trimming
#ifdef DEBUG_TRIM
  std::ofstream of1("par_edgs.g2");
  std::ofstream of1_2("space_edgs.g2");
 for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      if (!sfcv.get())
	continue;
      shared_ptr<ParamCurve> pcurve = sfcv->parameterCurve();
      if (pcurve.get())
	SplineDebugUtils::writeSpaceParamCurve(pcurve, of1, 0.0);
      shared_ptr<ParamCurve> scurve = sfcv->spaceCurve();
      if (scurve.get())
	{
	  scurve->writeStandardHeader(of1_2);
	  scurve->write(of1_2);
	}
    }
#endif
  int status = 0;
  HedgeSurface *hedge = getSurface(0);
  shared_ptr<ParamSurface> surf = hedge->surface();
  RectDomain dom = surf->containingDomain();
  double dom2[4];
  dom2[0] = dom.umin();
  dom2[1] = dom.umax();
  dom2[2] = dom.vmin();
  dom2[3] = dom.vmax();
  bool closed_u = false, closed_v = false;
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (elem.get())
    elem->isClosed(closed_u, closed_v);
  else
    {
      shared_ptr<SplineSurface> splsf =
	dynamic_pointer_cast<SplineSurface,ParamSurface>(surf);
      if (splsf.get())
	{
	  int md1 = GeometryTools::analyzePeriodicityDerivs(*splsf, 0, 0, tol);
	  int md2 = GeometryTools::analyzePeriodicityDerivs(*splsf, 1, 0, tol);
	  closed_u = (md1 >= 0);
	  closed_v = (md2 >= 0);
	}
    }

  //double diag = bbox_.low().dist(bbox_.high());
  if (closed_u)
    {
      // Add seem along v-direction to the edge collection
      double tdel = domain_[3] - domain_[2];
      double t1 = domain_[2] - 0.1*tdel;
      double t2 = domain_[3] + 0.1*tdel;
      shared_ptr<CurveOnSurface> seemcv1 = constParSfCv(surf, 0, dom2[0], 0, t1, t2);
      shared_ptr<CurveOnSurface> seemcv2 = constParSfCv(surf, 0, dom2[1], 1, t1, t2);
      shared_ptr<ftEdge> edg1(new ftEdge(hedge, seemcv1, seemcv1->startparam(),
					 seemcv1->endparam()));
      shared_ptr<ftEdge> edg2(new ftEdge(hedge, seemcv2, seemcv2->startparam(),
					 seemcv2->endparam()));
      int stat = 0;
      edg1->connectTwin(edg2.get(), stat);
      addTrimEdge(edg1);
      addTrimEdge(edg2);
    }
  if (closed_v)
    {
      // Add seem along u-direction to the edge collection
      double tdel = domain_[1] - domain_[0];
      double t1 = domain_[0] - 0.1*tdel;
      double t2 = domain_[1] + 0.1*tdel;
      shared_ptr<CurveOnSurface> seemcv1 = constParSfCv(surf, 1, dom2[2], 2, t1, t2);
      shared_ptr<CurveOnSurface> seemcv2 = constParSfCv(surf, 1, dom2[3], 3, t1, t2);
      shared_ptr<ftEdge> edg1(new ftEdge(hedge, seemcv1, seemcv1->startparam(),
					 seemcv1->endparam()));
      shared_ptr<ftEdge> edg2(new ftEdge(hedge, seemcv2, seemcv2->startparam(),
					 seemcv2->endparam()));
      int stat = 0;
      edg1->connectTwin(edg2.get(), stat);
      addTrimEdge(edg1);
      addTrimEdge(edg2);
    }
  
  bool do_bound = false;  // Necessary to trim only if not all curves are boundary curves
  vector<int> adjusted(trim_edgs_.size(), 0);
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      if (!sfcv.get())
	continue;
      if (sfcv->isConstantCurve())
	{
	  adaptOneEdge(trim_edgs_[ki], dom2);
	  bool same;
	  int bd = sfcv->whichBoundary(tol, same);
#ifdef DEBUG_TRIM
	  bool same_orient = sfcv->sameOrientation();
	  bool same_trace = sfcv->sameTrace(tol);
	  bool same_cv = sfcv->sameCurve(tol);
	  if ((!same_orient) || (!same_trace) || (!same_cv))
	    std::cout << "Surface curve mismatch " << ki << " " << sfcv << " " << same_orient << " " << same_trace << " " << same_cv << std::endl;
	    
#endif
	  if (bd < 0)
	    do_bound = true;
	}
      else
	{
	  do_bound = true;
	    
	  double t1 = sfcv->startparam();
	  double t2 = sfcv->endparam();

	  ftEdgeBase *twin0 = trim_edgs_[ki]->twin();
	  if (twin0)
	    {
	      ftEdge *twin = twin0->geomEdge();
	      shared_ptr<CurveOnSurface> sfcv2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(twin->geomCurve());
	      double t3 = sfcv2->startparam();
	      double t4 = sfcv2->endparam();
	      if (sfcv2->isConstantCurve() && (t3 > t1 || t4 < t2))
		{
		  int mod = 0;
		  if (t3 > t1)
		    mod += 1;
		  if (t4 < t2)
		    mod += 2;
		  double tmin = std::max(t1, t3);
		  double tmax = std::min(t2, t4);
		  shared_ptr<CurveOnSurface> sub(sfcv->subCurve(tmin, tmax));
		  shared_ptr<ftEdge> subedge(new ftEdge(sub, tmin, tmax));
		  ftEdgeBase *twin = trim_edgs_[ki]->twin();
		  trim_edgs_[ki]->disconnectTwin();
		  subedge->connectTwin(twin, status);
		  trim_edgs_[ki] = subedge;
		  adjusted[ki] = mod;
		}
	      int stop_break1 = 1;
	    }
	}
    }
  
  if (do_bound == false && trim_edgs_.size() != 4)
    {
      // Check if the parameter bound in the surface and in the existing trimming edges
      // is constent
      int stop_check = 1;
    }
  
#ifdef DEBUG_TRIM
  std::ofstream of2("par_edgs2.g2");
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      if (!sfcv.get())
	continue;
      shared_ptr<ParamCurve> pcurve = sfcv->parameterCurve();
      if (!pcurve.get())
	continue;
      SplineDebugUtils::writeSpaceParamCurve(pcurve, of2, 0.0);
    }
#endif

  if (do_bound)
    {
      bool OK = arrangeEdgeLoop(tol, adjusted);
      if (!OK)
	return false;
      if (trim_edgs_.size() == 0)
	return false;

      // Make bounded surface.
      // First collect trimming curves
      vector<shared_ptr<CurveOnSurface> > trim_cvs(trim_edgs_.size());
      // NB! Could be missing curves
      bool missing_par_cvs = false;
      for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> tmp_sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
	  trim_cvs[ki] = shared_ptr<CurveOnSurface>(tmp_sfcv->clone());  // To avoid inconsistencies
	  // if two edges share the same curve in geometry space
	  if (!trim_cvs[ki].get())
	    {
	      missing_par_cvs = true;
	      continue;  // Should not happen
	    }
	  if (!trim_cvs[ki]->hasParameterCurve())
	    missing_par_cvs = true;
	}

      if (missing_par_cvs)
	return false;  // Will fail
  
      double eps = 1.0e-9;
      vector<vector<shared_ptr<CurveOnSurface> > > trim_loops;

      // Check for separate and disconnected loops
      shared_ptr<ParamCurve> pcv1 = trim_cvs[0]->parameterCurve();
      Point startpt = pcv1->point(pcv1->startparam());
      int first_ix = 0;
      for (int ka=0; ka<(int)trim_cvs.size()-1; ++ka)
	{
	  int kb = (ka+1)%(int)trim_cvs.size();
	  shared_ptr<ParamCurve> pcv2 = trim_cvs[kb]->parameterCurve();
	  Point end = pcv1->point(pcv1->endparam());
	  Point start = pcv2->point(pcv2->startparam());
	  double next_dd = end.dist(start);
	  double close_dd = end.dist(startpt);
	  if (close_dd <= next_dd+eps && close_dd <= tol)
	    {
	      vector<shared_ptr<CurveOnSurface> > curr_cvs(trim_cvs.begin()+first_ix,
							   (kb!=0) ? trim_cvs.begin()+kb :
							   trim_cvs.end());
	      trim_loops.push_back(curr_cvs);
	      startpt = start;
	      first_ix = (kb!=0) ? kb : (int)trim_cvs.size();
	    }
	  else if (next_dd > tol)
	    {
	      // Open segment. Dismiss
	      startpt = start;
	      first_ix = (kb!=0) ? kb : (int)trim_cvs.size();
	    }
	      
	  // if (end.dist(start) > tol)
	  //   {
	  //     // Check for a closed loop
	  //     if (end.dist(startpt) <= tol)
	  // 	{
	  // 	  vector<shared_ptr<CurveOnSurface> > curr_cvs(trim_cvs.begin()+first_ix,
	  // 						       (kb!=0) ? trim_cvs.begin()+kb :
	  // 						       trim_cvs.end());
	  // 	  trim_loops.push_back(curr_cvs);
	  // 	  startpt = start;
	  // 	  first_ix = (kb!=0) ? kb : (int)trim_cvs.size();
	  // 	}
	  //     else
	  // 	{
	  // 	  MESSAGE("RevEngRegion::trimSurface(): Obsolete code. Why do we end here?");
		  
	  // 	}
	  //   }
	  pcv1 = pcv2;
	}

      if (first_ix < (int)trim_cvs.size())
	{
	  vector<shared_ptr<CurveOnSurface> > curr_cvs(trim_cvs.begin()+first_ix,
						       trim_cvs.end());
	  trim_loops.push_back(curr_cvs);
	}
      
#ifdef DEBUG_TRIM
      std::ofstream of3("par_edgs3.g2");
      std::ofstream of4("space_edgs3.g2");
      for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
	  if (!sfcv.get())
	    continue;
	  shared_ptr<ParamCurve> pcurve = sfcv->parameterCurve();
	  if (!pcurve.get())
	    continue;
	  SplineDebugUtils::writeSpaceParamCurve(pcurve, of3, 0.0);
	}
      for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
	  if (!sfcv.get())
	    continue;
	  shared_ptr<ParamCurve> scurve = sfcv->spaceCurve();
	  scurve->writeStandardHeader(of4);
	  scurve->write(of4);
	}
      std::ofstream of3_2("par_edgs3_2.g2");
      std::ofstream of4_2("space_edgs3_2.g2");
      for (size_t ki=0; ki<trim_loops.size(); ++ki)
	for (size_t kj=0; kj<trim_loops[ki].size(); ++kj)
	  {
	    shared_ptr<CurveOnSurface> sfcv = trim_loops[ki][kj];
	    shared_ptr<ParamCurve> pcurve = sfcv->parameterCurve();
	    if (pcurve.get())
	      SplineDebugUtils::writeSpaceParamCurve(pcurve, of3_2, 0.0);
	    shared_ptr<ParamCurve> scurve = sfcv->spaceCurve();
	    scurve->writeStandardHeader(of4_2);
	    scurve->write(of4_2);
	  }
	
#endif
      // Check orientation
      vector<CurveLoop> loops;
      if (trim_loops.size() == 1)
	{
	  loops.push_back(CurveLoop(trim_loops[0], tol));
	  bool is_CCW = LoopUtils::paramIsCCW(trim_loops[0], tol, tol);
	  if (!is_CCW)
	    {
	      loops[0].turnOrientation();
	      // Must also turn edges
	    }
	}
      else
	{
	  // Sort curves. Expects one outer loop and one or more inner loops
	  vector<Point> pt_on_cv(trim_loops.size());
	  for (size_t ki=0; ki<trim_loops.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> pcv = trim_loops[ki][0]->parameterCurve();
	      pt_on_cv[ki] = pcv->point(0.5*(pcv->startparam()+pcv->endparam()));
	    }

	  size_t ki;
	  for (ki=0; ki<trim_loops.size(); ++ki)
	    {
	      // Check if the point from all other loops lies inside current
	      vector<shared_ptr<ParamCurve> > par_cvs(trim_loops[ki].size());
	      for (size_t kr=0; kr<trim_loops[ki].size(); ++kr)
		{
		  par_cvs[kr] = trim_loops[ki][kr]->parameterCurve();
		}

	      shared_ptr<CurveLoop> cvloop(new CurveLoop(par_cvs, tol, false));
	      CurveBoundedDomain cvdom(cvloop);
	      size_t kj;
	      for (kj=0; kj<pt_on_cv.size(); ++kj)
		{
		  if (kj == ki)
		    continue;
		  Vector2D ppos(pt_on_cv[kj][0], pt_on_cv[kj][1]);
		  bool inside = cvdom.isInDomain(ppos, tol);
		  if (!inside)
		    break;
		}
	      if (kj == pt_on_cv.size())
		break; // Outer boundary found
	    }
	  if (ki < trim_loops.size())
	    std::swap(trim_loops[0], trim_loops[ki]);
	  
	  loops.push_back(CurveLoop(trim_loops[0], tol));
	  bool is_CCW = LoopUtils::paramIsCCW(trim_loops[0], tol, tol);
	  if (!is_CCW)
	    {
	      loops[0].turnOrientation();
	    }
	  for (size_t kr=1; kr<trim_loops.size(); ++kr)
	    {
	      loops.push_back(CurveLoop(trim_loops[kr], tol));
	      bool is_CCW = LoopUtils::paramIsCCW(trim_loops[kr], tol, tol);
	      if (is_CCW)
		{
		  loops[kr].turnOrientation();
		}
	    }
	}

      shared_ptr<BoundedSurface> bdsurf(new BoundedSurface(surf, loops));
#ifdef DEBUG_TRIM
      int valid_state;
      bool valid = bdsurf->isValid(valid_state);
      std::cout << "BoundedSurf is valid? " << valid << " " << valid_state << std::endl;
      std::ofstream of5("bounded_surf.g2");
      bdsurf->writeStandardHeader(of5);
      bdsurf->write(of5);
#endif
      associated_sf_[0]->replaceSurf(bdsurf);

      // // Add existing trimming curves. Later
      // for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
      // 	trim_edgs_[ki]->setFace(associated_sf_[0]);
      // vector<shared_ptr<ftEdgeBase> > tmp_trim(trim_edgs_.begin(), trim_edgs_.end());
      // shared_ptr<Loop> edge_loop(new Loop(associated_sf_[0], tmp_trim, tol));
      // associated_sf_[0]->addOuterBoundaryLoop(edge_loop);
      associated_sf_[0]->clearInitialEdges();
      (void)associated_sf_[0]->createInitialEdges();
#ifdef DEBUG_TRIM
      for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
	  Point startpt1 = sfcv->ParamCurve::point(sfcv->startparam());
	  Point startpt2 = trim_edgs_[ki]->point(trim_edgs_[ki]->tMin());
	  int stop_breakc = 1;
	}
#endif
    }  // End do_bound
  
  return true;
}

struct CloseCvInfo
{
  double par1_, par2_, dist_;
  double par3_, par4_, dist2_;

  CloseCvInfo()
  {
    par1_ = par2_ = par3_ = par4_ = 0.0;
    dist_ = dist2_ = -1.0;
  }
  
  CloseCvInfo(double par1, double par2, double dist)
  {
    par1_ = par1;
    par2_ = par2;
    dist_ = dist;
    par3_ = par4_ = 0.0;
    dist2_ = -1.0;
  }
  
  CloseCvInfo(double par1, double par2, double dist, 
	      double par3, double par4, double dist2)
  {
    par1_ = par1;
    par2_ = par2;
    dist_ = dist;
    par3_ = par3;
    par4_ = par4;
    dist2_ = dist2;
  }
};

CloseCvInfo getCloseInfo(double tol, shared_ptr<ParamCurve>& pcurve1, int adjusted1,
			 shared_ptr<ParamCurve>& pcurve2, int adjusted2)
{
  Point close1, close2;
  int status = 0;
  double tmin1 = pcurve1->startparam();
  double tmax1 = pcurve1->endparam();
  double seed1 = 0.5*(tmin1 + tmax1);
  Point pos1 = pcurve1->point(tmin1);
  Point pos2 = pcurve1->point(tmax1);
  double tmin2 = pcurve2->startparam();
  double tmax2 = pcurve2->endparam();
  double seed2 = 0.5*(tmin2 + tmax2);
  Point pos3 = pcurve2->point(tmin2);
  Point pos4 = pcurve2->point(tmax2);
  double par1[9], par2[9];
  par1[0] = par1[1] = tmin1;
  par1[2] = par1[3] = tmax1;
  par1[5] = tmin1;
  par1[6] = tmax1;
  par2[7] = tmin2;
  par2[8] = tmax2;
  par2[0] = par2[2] = tmin2;
  par2[1] = par2[3] = tmax2;
  par1[4] = par1[7] = par1[8] = 0.5*(tmin1+tmax1);
  par2[4] = par2[5] = par2[6] = 0.5*(tmin2+tmax2);
  double dist[9];
  dist[0] = pos1.dist(pos3);
  dist[1] = pos1.dist(pos4);
  dist[2] = pos2.dist(pos3);
  dist[3] = pos2.dist(pos4);
  dist[4] = dist[5] = dist[6] = dist[7] = dist[8] = std::numeric_limits<double>::max();
  if (adjusted1 < 3 && adjusted2 < 3)
    ClosestPoint::closestPtCurves2D(pcurve1.get(), pcurve2.get(), tol,
				    tmin1, tmax1, tmin2, tmax2, seed1,
				    seed2, 1, false, par1[4], par2[4],
				    dist[4], close1, close2, status);
  if ((adjusted1 == 1 || adjusted1 == 3) && adjusted2 < 3)
    pcurve2->closestPoint(pos1, tmin2, tmax2, par2[5], close2, dist[5]);
  if ((adjusted2 == 1 || adjusted1 == 3) && adjusted2 < 3)
    pcurve2->closestPoint(pos2, tmin2, tmax2, par2[6], close2, dist[6]);
  if ((adjusted2 == 1 || adjusted2 == 3) && adjusted1 < 3)
    pcurve1->closestPoint(pos3, tmin1, tmax1, par1[7], close1, dist[7]);
  if ((adjusted2 == 2 || adjusted2 == 3) && adjusted1 < 3)
    pcurve1->closestPoint(pos4, tmin1, tmax1, par1[8], close1, dist[8]);

  double eps = std::min(0.1*tol, 1.0e-4);
  double eps2 = 1.0e-9;
  for (int ka=4; ka<9; ++ka)
    {
      if (par1[ka]-tmin1 < eps2)
	par1[ka] = tmin1;
      if (tmax1-par1[ka] < eps2)
	par1[ka] = tmax1;
      if (par2[ka]-tmin2 < eps2)
	par2[ka] = tmin2;
      if (tmax2-par2[ka] < eps2)
	par2[ka] = tmax2;
    }
  
  if (adjusted2 == 1 || adjusted2 == 3)
    {
      if (par1[4] - tmin1 > eps2 && par1[4] - tmin1 < eps)
	dist[4] = std::numeric_limits<double>::max();
      if (par1[7] - tmin1 > eps2 && par1[7] - tmin1 < eps)
	dist[7] = std::numeric_limits<double>::max();
      if (par1[8] - tmin1 > eps2 && par1[8] - tmin1 < eps)
	dist[8] = std::numeric_limits<double>::max();
    }
  if (adjusted2 == 2 || adjusted2 == 3)
    {
      if (tmax1 - par1[4] > eps2 && tmax1 - par1[4] < eps)
	dist[4] = std::numeric_limits<double>::max();
      if (tmax1 - par1[7] > eps2 && tmax1 - par1[7] < eps)
	dist[7] = std::numeric_limits<double>::max();
      if (tmax1 - par1[8] > eps2 && tmax1 - par1[8] < eps)
	dist[8] = std::numeric_limits<double>::max();
    }

  if (adjusted1 == 1 || adjusted1 == 3)
    {
      if (par2[4] - tmin2 > eps2 && par2[4] - tmin2 < eps)
	dist[4] = std::numeric_limits<double>::max();
      if (par2[5] - tmin2 > eps2 && par2[5] - tmin2 < eps)
	dist[5] = std::numeric_limits<double>::max();
      if (par2[6] - tmin2 > eps2 && par2[6] - tmin2 < eps)
	dist[6] = std::numeric_limits<double>::max();
    }
  if (adjusted1 == 2 || adjusted1 == 3)
    {
      if (tmax2 - par2[4] > eps2 && tmax2 - par2[4] < eps)
	dist[4] = std::numeric_limits<double>::max();
      if (tmax2 - par2[5] > eps2 && tmax2 - par2[5] < eps)
	dist[5] = std::numeric_limits<double>::max();
      if (tmax2 - par2[6] > eps2 && tmax2 - par2[6] < eps)
	dist[6] = std::numeric_limits<double>::max();
    }


  int ka, kb;
  for (ka=0; ka<9; ++ka)
    for (kb=ka+1; kb<9; ++kb)
      {
	// if ((kb < 4 && dist[kb] < dist[ka]) || (kb == 4 && dist[kb] < dist[ka]-tol))
	if (dist[kb] < dist[ka])
	  {
	    std::swap(dist[ka], dist[kb]);
	    std::swap(par1[ka], par1[kb]);
	    std::swap(par2[ka], par2[kb]);
	  }
      }

  if (fabs(par1[1] - par1[0]) < eps || fabs(par2[1]-par2[0]) < eps)
    {
      par1[1] = par1[2];
      par2[1] = par2[2];
      dist[1] = dist[2];
    }
  CloseCvInfo curr_info(par1[0], par2[0], dist[0], par1[1],
			par2[1], dist[1]);
  return curr_info;
 }


//===========================================================================
bool RevEngRegion::arrangeEdgeLoop(double tol, vector<int>& adjusted)
//===========================================================================
{
  vector<vector<CloseCvInfo> > info(trim_edgs_.size());
  vector<shared_ptr<CurveOnSurface> > sfcv(trim_edgs_.size());
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      info[ki].resize(trim_edgs_.size());
      sfcv[ki] =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
    }

  // Collect distance info
  double tmin1, tmax1; 
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      if (!sfcv[ki].get())
	continue;
      shared_ptr<ParamCurve> pcurve1 = sfcv[ki]->parameterCurve();
      if (!pcurve1.get())
	return false;
      
      tmin1 = pcurve1->startparam();
      tmax1 = pcurve1->endparam();
      Point pos1 = pcurve1->point(tmin1);
      Point pos2 = pcurve1->point(tmax1);
      double dist0 = pos1.dist(pos2);
      info[ki][ki] = CloseCvInfo(tmin1, tmax1, dist0);
      for (size_t kj=ki+1; kj<trim_edgs_.size(); ++kj)
	{
	  if (!sfcv[kj].get())
	    continue;
	  shared_ptr<ParamCurve> pcurve2 = sfcv[kj]->parameterCurve();
	  if (!pcurve2.get())
	    return false;

	  CloseCvInfo curr_info = getCloseInfo(tol, pcurve1, adjusted[ki],
					       pcurve2, adjusted[kj]);
	  info[ki][kj] = info[kj][ki] = curr_info;
 	}
    }

  //#if 0
  // Identify gaps
  double tol2 = 10.0*tol;
  vector<pair<int, int> > adj_gap;
  for (size_t ki=0; ki<info.size(); ++ki)
    {
      double ix1 = -1, ix2 = -1;
      for (size_t kj=0; kj<info[ki].size(); ++kj)
	{
	  // if (kj == ki)
	  //   continue;
	  double dd = info[ki][kj].dist_;
	  if (dd <= tol2)
	    {
	      int num_small = 0;
	      for (size_t kr=0; kr<info[kj].size(); ++kr)
		if (info[kj][kr].dist_ < dd)
		  num_small++;
	      if (num_small < 2)
		{
		  if (ix1 < 0)
		    ix1 = (int)kj;
		  else if (ix2 < 0)
		    ix2 = (int)kj;
		}
	    }
	}
      if (ix1 < 0 || (ix2 < 0 && ix1 != (int)ki))
	adj_gap.push_back(std::make_pair((int)ki, ix1));
    }

  vector<double> limit_par;
  if (adj_gap.size() > 0)
    {
      for (size_t ki=0; ki<adj_gap.size(); ++ki)
	{
	  // Fetch adjacent face
	  int cv_ix = adj_gap[ki].first;
	  ftEdgeBase* twin0 = trim_edgs_[cv_ix]->twin();
	  if (!twin0)
	    continue;
	  ftEdge* twin = twin0->geomEdge();
	  ftFaceBase *twin_face = twin->face();
	  if (!twin_face)
	    continue;
	  HedgeSurface *adj_hedge = dynamic_cast<HedgeSurface*>(twin_face);
	  if (!adj_hedge)
	    continue;
	  RevEngRegion* other_reg = adj_hedge->getRegion(0);
	  vector<RevEngPoint*> edge_pts;
	  if (other_reg == this)
	    {
	      // Fetch points at seem
	      int dir;
	      double parval;
	      sfcv[cv_ix]->isConstantCurve(tol2, dir, parval);
	      bool udir = (dir == 1);   // Should be sufficient
	      vector<RevEngPoint*> seam_pts1, seam_pts2;
	      extractPointsAtSeam(seam_pts1, seam_pts2, udir);
	      int ix = 1 - dir;
	      if (parval-domain_[2*ix] < domain_[2*ix+1]-parval)
		edge_pts = seam_pts1;
	      else
		edge_pts = seam_pts2;
	      int stop_b1 = 1;
	    }
	  else
	    {
	      vector<RevEngRegion*> others;
	      others.push_back(other_reg);
	      edge_pts = extractBdPoints(others);
	      int stop_b2 = 1;
	    }

	  // Identify first and last point along edge
	  RevEngPoint *first_pt=0, *last_pt=0;
	  double t1, t2;
	  RevEngUtils::identifyEndPoints(edge_pts, sfcv[cv_ix], first_pt, t1, last_pt, t2);
#ifdef DEBUG_TRIM
	  if (t2 > t1)
	    {
	      std::ofstream ofpt("edge_ends.g2");
	      if (first_pt)
		{
		  ofpt << "400 1 0 4 100 100 55 255" << std::endl;
		  ofpt << "1" << std::endl;
		  ofpt << first_pt->getPar() << " 0.0" << std::endl;
		}
	      if (last_pt)
		{
		  ofpt << "400 1 0 4 100 100 55 255" << std::endl;
		  ofpt << "1" << std::endl;
		  ofpt << last_pt->getPar() << " 0.0" << std::endl;
		}
	    }
#endif
	  if (t2 < t1)
	    std::swap(t1, t2);
	  if (adj_gap[ki].second < 0)
	    {
	      limit_par.push_back(t1);
	      limit_par.push_back(t2);
	    }
	  else
	    {
	      int ix2 = adj_gap[ki].second;
	      double t3 = (cv_ix < ix2) ? info[cv_ix][ix2].par1_ : info[cv_ix][ix2].par2_;
	      if (fabs(t3 - t1) < fabs(t3 - t2))
		limit_par.push_back(t2);
	      else limit_par.push_back(t1);
	    }
	  int stop_break0 = 1;
	}

      // Sort gaps
      vector<pair<pair<shared_ptr<CurveOnSurface>,double>,
		  pair<shared_ptr<CurveOnSurface>,double> > > gap_cv_bounds;
      vector<bool> use_lower;
      if (adj_gap.size() == 1 && limit_par.size() == 2)
	{
	  int cv_ix = adj_gap[0].first;
	  gap_cv_bounds.push_back(make_pair(make_pair(sfcv[cv_ix],limit_par[0]),
					    make_pair(sfcv[cv_ix],limit_par[1])));
	  use_lower.push_back(true);  // Not used
	}
      else if (adj_gap.size() == 2 && limit_par.size() == 2)
	{
	  int cv_ix1 = adj_gap[0].first;
	  int cv_ix2 = adj_gap[1].first;
	  double t3;
	  if (trim_edgs_[cv_ix1]->twin() == trim_edgs_[cv_ix2].get())
	    {
	      // Seam curve. Unify end parameters
	      int ix2 = adj_gap[0].second;
	      t3 = (cv_ix1 < ix2) ? info[cv_ix1][ix2].par1_ : info[cv_ix1][ix2].par2_;
	      if (fabs(t3-limit_par[0]) < fabs(t3-limit_par[1]))
		limit_par[0] = limit_par[1];
	      else
		limit_par[1] = limit_par[0];
	    }
	  gap_cv_bounds.push_back(make_pair(make_pair(sfcv[cv_ix1],limit_par[0]),
					    make_pair(sfcv[cv_ix2],limit_par[1])));
	  use_lower.push_back((limit_par[0]<t3));
	}
      else
	MESSAGE("RevEngRegion::arrangeEdgeLoop(): Gap configuration not implemented");

      HedgeSurface *hedge = getSurface(0);
      shared_ptr<ParamSurface> surf = hedge->surface();
      size_t num = trim_edgs_.size();
      for (size_t ki=0; ki<gap_cv_bounds.size(); ++ki)
	{
	  shared_ptr<CurveOnSurface> sfcv1 = gap_cv_bounds[ki].first.first;
	  shared_ptr<CurveOnSurface> sfcv2 = gap_cv_bounds[ki].second.first;
	  int dir1, dir2;
	  double parval1, parval2;
	  sfcv1->isConstantCurve(tol, dir1, parval1);
	  sfcv2->isConstantCurve(tol, dir2, parval2);

	  shared_ptr<ParamCurve> pcv1 = sfcv1->parameterCurve();
	  shared_ptr<ParamCurve> pcv2 = sfcv2->parameterCurve();
	  Point ppos1 = pcv1->point(gap_cv_bounds[ki].first.second);
	  Point ppos2 = pcv2->point(gap_cv_bounds[ki].second.second);
	  if (dir1 >= 0 && dir1 == dir2)
	    {
	      // Set surface limit based on parameter domain of points
	      int p_ix = 2 - dir1;
	      // if (std::min(fabs(ppos1[p_ix]-domain_[2*p_ix]), fabs(ppos1[p_ix]-domain_[2*p_ix])) <
	      // 	  std::min(fabs(ppos1[p_ix]-domain_[2*p_ix+1]), fabs(ppos1[p_ix]-domain_[2*p_ix+1])))
	      if (use_lower[ki])
		ppos1[p_ix] = ppos2[p_ix] = domain_[2*p_ix];
	      else
		ppos1[p_ix] = ppos2[p_ix] = domain_[2*p_ix+1];
	    }
	  shared_ptr<SplineCurve> gap_par(new SplineCurve(ppos1, ppos2));
	  shared_ptr<CurveOnSurface> gap_sfcv(new CurveOnSurface(surf, gap_par, true));
	  gap_sfcv->ensureSpaceCrvExistence(tol);
	  shared_ptr<ftEdge> gap_edge(new ftEdge(hedge, gap_sfcv, gap_sfcv->startparam(),
						 gap_sfcv->endparam()));
	  sfcv.push_back(gap_sfcv);
	  trim_edgs_.push_back(gap_edge);
	  adjusted.push_back(true);
	}

      if (sfcv.size() > num)
	{
	  // Extend meeting information
	  info.resize(sfcv.size());
	  for (size_t ki=0; ki<sfcv.size(); ++ki)
	    info[ki].resize(sfcv.size());

	  for (size_t ki=num; ki<trim_edgs_.size(); ++ki)
	    {
	      if (!sfcv[ki].get())
		continue;
	      shared_ptr<ParamCurve> pcurve1 = sfcv[ki]->parameterCurve();
	      if (!pcurve1.get())
		continue;
      
	      tmin1 = pcurve1->startparam();
	      tmax1 = pcurve1->endparam();
	      Point pos1 = pcurve1->point(tmin1);
	      Point pos2 = pcurve1->point(tmax1);
	      double dist0 = pos1.dist(pos2);
	      info[ki][ki] = CloseCvInfo(tmin1, tmax1, dist0);
	      for (size_t kj=0; kj<trim_edgs_.size(); ++kj)
		{
		  if (kj >= num && kj <= ki)
		    continue;
		  if (!sfcv[kj].get())
		    continue;
		  shared_ptr<ParamCurve> pcurve2 = sfcv[kj]->parameterCurve();
		  if (!pcurve2.get())
		    continue;

		  CloseCvInfo curr_info = (ki > kj) ?
		    getCloseInfo(tol, pcurve2, adjusted[kj], pcurve1, adjusted[ki]) :
		    getCloseInfo(tol, pcurve1, adjusted[ki], pcurve2, adjusted[kj]);;
		  info[ki][kj] = info[kj][ki] = curr_info;
		}
	    }

	}
      int stop_break = 1;
    }
  //#endif  
  // Check if any curves must be reduced, and record previous and next curve
  double eps = 1.0e-9;
  int status = 0;
  vector<pair<int, int> > seq(trim_edgs_.size());
  vector<pair<double, double> > param(trim_edgs_.size());
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      int ix1=-1, ix2=-1;
      double td1 = std::numeric_limits<double>::max();
      double td2 = std::numeric_limits<double>::max();
      double td0 = info[ki][ki].dist_;
      double t1, t2;
      t1 = t2 = 2.0*fabs(std::max(sfcv[ki]->startparam(), sfcv[ki]->endparam()));
      double delta = 0.001*fabs(info[ki][ki].par2_ - info[ki][ki].par1_);
      for (size_t kj=0; kj<info[ki].size(); ++kj)
	{
	  if (ki == kj)
	    continue;
	  double dd = info[ki][kj].dist_;
	  double tp1 = (ki < kj) ? info[ki][kj].par1_ : info[ki][kj].par2_;
	  if (ix1 < 0 || dd < td1)
	    {
	      ix1 = (int)kj;
	      td1 = dd;
	      t1 = tp1;
	    }
	}

      
      for (size_t kj=0; kj<info[ki].size(); ++kj)
	{
	  if (ki == kj)
	    continue;
	  double dd = info[ki][kj].dist_;
	  double dd2 = info[ki][kj].dist2_;
	  double tp1 = (ki < kj) ? info[ki][kj].par1_ : info[ki][kj].par2_;
	  double tp2 = (ki < kj) ? info[ki][kj].par3_ : info[ki][kj].par4_;
	  //double td12 = std::max(td1, td2);
	  if (ix2 < 0 || dd < td2)
	    {
	      if (fabs(tp1 - t1) > delta)
		{
		  ix2 = (int)kj;
		  td2 = dd;
		  t2 = tp1;
		}
	      else if ((ix2 < 0 || dd2 < td2) && fabs(tp2 - t1) > delta)
		{
		  ix2 = (int)kj;
		  td2 = dd2;
		  t2 = tp2;
		}
	    }
	}

      
	  // if (ix1 < 0 || dd < td1)
	  //   {
	  //     if (ix1 < 0 || (ix2 >= 0 && td1 >= td2 && (!fabs(tp1 - t2) < eps)))
	  // 	{
	  // 	  ix1 = (int)kj;
	  // 	  td1 = dd;
	  // 	  t1 = tp1;
	  // 	  set1 = true;
	  // 	}
	  //     else if (ix1 < 0 || (ix2 >= 0 && dd2 < td1 &&
	  // 			   td1 >= td2 && (!fabs(tp2 - t2) < eps)))
	  // 	{
	  // 	  ix1 = (int)kj;
	  // 	  td1 = dd2;
	  // 	  t1 = tp2;
	  // 	  set2 = true;
	  // 	}
	  //   }
	  // if (ix2 < 0 || dd < td2)
	  //   {
	  //     if ((!set1) && ix1 >= 0 && td2 > td1 && (!fabs(tp1 - t1) < eps))
	  // 	{
	  // 	  ix2 = (int)kj;
	  // 	  td2 = dd;
	  // 	  t2 = tp1;
	  // 	}
	  //     else if ((!set2) && ix1 >= 0 && dd2 < td2 &&
	  // 	       td2 > td1 && (!fabs(tp2 - t1) < eps))
	  // 	{
	  // 	  ix2 = (int)kj;
	  // 	  td2 = dd2;
	  // 	  t2 = tp2;
	  // 	}
	  //   }
	  // if (ix1 < 0 || ((dd < td12-eps ||
	  // 		   (dd < td12 && fabs(t2-tp1) > eps)) && td2 <= td1))
	  //   {
	  //     ix1 = (int)kj;
	  //     td1 = dd;
	  //     t1 = tp1;
	  //   }
	  // else if (ix2 < 0 || ((dd < td12-eps ||
	  // 			(dd < td12 && fabs(t1-tp1) > eps)) && td1 < td2))
	  //   {
	  //     ix2 = (int)kj;
	  //     td2 = dd;
	  //     t2 = tp1;
	  //   }

	  // td12 = std::max(td1, td2);
	  // if (dd2 < 0.0)
	  //   dd2 = td12;
	  // if (ix1 < 0 || ((dd2 < td12-eps ||
	  // 		  (dd2 < td12 && fabs(t2-tp2) > eps)) && td2 <= td1))
	  //   {
	  //     ix1 = (int)kj;
	  //     td1 = dd2;
	  //     t1 = tp2;
	  //   }
	  // else if (ix2 < 0 || ((dd2 < td12-eps ||
	  // 			(dd2 < td12 && fabs(t1-tp2) > eps)) && td1 < td2))
	  //   {
	  //     ix2 = (int)kj;
	  //     td2 = dd2;
	  //     t2 = (ki < kj) ? info[ki][kj].par3_ : info[ki][kj].par4_;
	  //   }
      //}
      if (/*td1 > tol && td2 > tol &&*/ td0 < std::min(td1, td2))
	{
	  // Closed loop
	  ix1 = ix2 = (int)ki;
	  t1 = sfcv[ki]->startparam();
	  t2 = sfcv[ki]->endparam();
	}
      
      if (t2 < t1)
	{
	  std::swap(t1, t2);
	  std::swap(ix1, ix2);
	}
      seq[ki] = std::make_pair(ix1, ix2);
      param[ki] = std::make_pair(t1, t2);

      if ((t1 > sfcv[ki]->startparam()+eps || t2 < sfcv[ki]->endparam()-eps) && t2 > t1)
	{
	  if ((td1 > tol || td2 > tol) &&
	      std::min(t1-sfcv[ki]->startparam(), sfcv[ki]->endparam()-t2) < tol)
	    {
	      // Extra testing
	      shared_ptr<ParamCurve> pcrv = sfcv[ki]->parameterCurve();
	      double tpar = (t1-sfcv[ki]->startparam() > sfcv[ki]->endparam()-t2) ? t1 : t2;
	      double t3 = 0.75*tpar + 0.25*pcrv->startparam();
	      double t4 = 0.75*tpar + 0.25*pcrv->endparam();
	      Point ppar1 = pcrv->point(t3);
	      Point ppar2  = pcrv->point(t4);
	      Point close1(std::max(domain_[0], std::min(domain_[1], ppar1[0])),
			   std::max(domain_[2], std::min(domain_[3], ppar1[1])));
	      Point close2(std::max(domain_[0], std::min(domain_[1], ppar2[0])),
			   std::max(domain_[2], std::min(domain_[3], ppar2[1])));
	      double d1 = ppar1.dist(close1);
	      double d2 = ppar2.dist(close2);
	      if (tpar == t1 && d1<d2)
		{
		  t2 = sfcv[ki]->startparam();
		  std::swap(t1, t2);
		  std::swap(ix1, ix2);
		  seq[ki] = std::make_pair(ix1, ix2);
		}
	      else if (tpar == t2 && d2<d1)
		{
		  t1 = sfcv[ki]->endparam();
		  std::swap(t1, t2);
		  std::swap(ix1, ix2);
		  seq[ki] = std::make_pair(ix1, ix2);
		}
	      int stop_d = 1;
	    }
	  shared_ptr<CurveOnSurface> sub(sfcv[ki]->subCurve(t1, t2));
	  shared_ptr<ftEdge> subedge(new ftEdge(sub, t1, t2));
	  ftEdgeBase *twin = trim_edgs_[ki]->twin();
	  if (twin)
	    {
	      trim_edgs_[ki]->disconnectTwin();
	      subedge->connectTwin(twin, status);
	    }
	  trim_edgs_[ki] = subedge;
	}
    }

  // Could be necessary to remove extra entrances of closed edges
  
  vector<int> seq_ix;
  vector<bool> turn;
  seq_ix.push_back(0);
  turn.push_back(false);
  vector<size_t> start_ix;
  start_ix.push_back(0);
  while (seq_ix.size() < trim_edgs_.size())
    {
      size_t ix = seq_ix.size() - 1;
      for (size_t ki=ix; ki<seq_ix.size(); ++ki)
	{
	  // Select next curve
	  int ix1 = seq[seq_ix[ki]].first;
	  int ix2 = seq[seq_ix[ki]].second;
	  double t1 = param[seq_ix[ki]].first; 
	  double t2 = param[seq_ix[ki]].second;
	  Point pt1 = sfcv[seq_ix[ki]]->parameterCurve()->point(t1);
	  Point pt2 = sfcv[seq_ix[ki]]->parameterCurve()->point(t2);
	  for (size_t kr=0; kr<seq_ix.size(); ++kr)
	    {
	      if (seq_ix[kr] == ix1)
		ix1 = -1;
	      if (seq_ix[kr] == ix2)
		ix2 = -1;
	    }

	  if (ix1 >= 0 && turn[ki])
	    {
	      seq_ix.push_back(ix1);
	      Point pt3 = sfcv[ix1]->parameterCurve()->point(param[ix1].first);
	      Point pt4 = sfcv[ix1]->parameterCurve()->point(param[ix1].second);
	      if (pt1.dist(pt3) > pt1.dist(pt4))
		turn.push_back(true);
	      else
		turn.push_back(false);
	    }
	  else if (ix2 >= 0 && (!turn[ki]))
	    {
	      seq_ix.push_back(ix2);
	      Point pt3 = sfcv[ix2]->parameterCurve()->point(param[ix2].first);
	      Point pt4 = sfcv[ix2]->parameterCurve()->point(param[ix2].second);
	      if (pt2.dist(pt3) > pt2.dist(pt4))
		turn.push_back(true);
	      else
		turn.push_back(false);
	    }
	}

      for (size_t ki=seq_ix.size()-1; ki<seq_ix.size(); ++ki)
	{
	  // Select previous curve
	  int ix1 = seq[seq_ix[ix]].first;
	  int ix2 = seq[seq_ix[ix]].second;
	  double t1 = param[seq_ix[ix]].first;
	  double t2 = param[seq_ix[ix]].second;
	  Point pt1 = sfcv[seq_ix[ix]]->parameterCurve()->point(t1);
	  Point pt2 = sfcv[seq_ix[ix]]->parameterCurve()->point(t2);
	  for (size_t kr=0; kr<seq_ix.size(); ++kr)
	    {
	      if (seq_ix[kr] == ix1)
		ix1 = -1;
	      if (seq_ix[kr] == ix2)
		ix2 = -1;
	    }

	  if (ix1 >= 0 && (!turn[ix]))
	    {
	      seq_ix.insert(seq_ix.begin()+ix, ix1);
	      Point pt3 = sfcv[ix1]->parameterCurve()->point(param[ix1].first);
	      Point pt4 = sfcv[ix1]->parameterCurve()->point(param[ix1].second);
	      if (pt1.dist(pt4) < pt1.dist(pt3)) 
		turn.insert(turn.begin()+ix, false);
	      else
		turn.insert(turn.begin()+ix, true);
	    }
	  else if (ix2 >= 0 && turn[ix])
	    {
	      seq_ix.push_back(ix2);
	      Point pt3 = sfcv[ix2]->parameterCurve()->point(param[ix2].first);
	      Point pt4 = sfcv[ix2]->parameterCurve()->point(param[ix2].second);
	      if (pt2.dist(pt3) > pt2.dist(pt4))
		turn.insert(turn.begin()+ix, false);
	      else
		turn.insert(turn.begin()+ix, true);
	    }
	}
      // Select unused edge
      size_t kr, kh;
      for (kr=0; kr<trim_edgs_.size(); ++kr)
	{
	  for (kh=0; kh<seq_ix.size(); ++kh)
	    if (seq_ix[kh] == (int)kr)
	      break;

	  if (kh == seq_ix.size())
	    {
	      start_ix.push_back(seq_ix.size());
	      seq_ix.push_back((int)kr);
	      turn.push_back(false);
	      break;
	    }
	}
    }
  start_ix.push_back(seq_ix.size());
  // vector<int> seq_ix(trim_edgs_.size(), -1);
  // vector<bool> turn(trim_edgs_.size(), false);
  // seq_ix[0] = 0;
  // for (size_t ki=0; ki<trim_edgs_.size()-1; ++ki)
  //   {
  //     // Select next curve
  //     int ix1 = seq[seq_ix[ki]].first;
  //     int ix2 = seq[seq_ix[ki]].second;
  //     double t1 = param[seq_ix[ki]].first; //(seq_ix[ki] < ix1) ? info[seq_ix[ki]][ix1].par1_ : info[seq_ix[ki]][ix1].par2_;
  //     double t2 = param[seq_ix[ki]].second; //(seq_ix[ki] < ix2) ? info[seq_ix[ki]][ix2].par1_ : info[seq_ix[ki]][ix2].par2_;
  //     Point pt1 = sfcv[ki]->ParamCurve::point(t1);
  //     Point pt2 = sfcv[ki]->ParamCurve::point(t2);
  //     for (size_t kr=0; kr<ki; ++kr)
  // 	{
  // 	  if (seq_ix[kr] == ix1)
  // 	    ix1 = -1;
  // 	  if (seq_ix[kr] == ix2)
  // 	    ix2 = -1;
  // 	}
  //     if (ix1 != ix2 && ix1 >= 0 && ((turn[ki] && t1<t2) || ((!turn[ki]) && t2<t1)))
  // 	{
  // 	  seq_ix[ki+1] = ix1;
  // 	  //double t3 = (seq_ix[ki] < ix1) ? info[seq_ix[ki]][ix1].par2_ : info[seq_ix[ki]][ix1].par1_;
  // 	  //if (fabs(sfcv[ix1]->endparam()-t3) < fabs(t3-sfcv[ix1]->startparam()))
  // 	  Point pt3 = sfcv[ix1]->ParamCurve::point(sfcv[ix1]->startparam());
  // 	  Point pt4 = sfcv[ix1]->ParamCurve::point(sfcv[ix1]->endparam());
  // 	  if (pt1.dist(pt4)+pt2.dist(pt3) < pt1.dist(pt3)+pt2.dist(pt4))
  // 	    turn[ki+1] = true;
  // 	}
  //     else if (ix1 != ix2 && ix2 >= 0)
  // 	{
  // 	  seq_ix[ki+1] = ix2;
  // 	  //double t3 = (seq_ix[ki] < ix2) ? info[seq_ix[ki]][ix2].par2_ : info[seq_ix[ki]][ix2].par1_;
  // 	  //if (fabs(sfcv[ix2]->endparam()-t2) < fabs(t3-sfcv[ix2]->startparam()))
  // 	  Point pt3 = sfcv[ix2]->ParamCurve::point(sfcv[ix2]->startparam());
  // 	  Point pt4 = sfcv[ix2]->ParamCurve::point(sfcv[ix2]->endparam());
  // 	  if (pt1.dist(pt4)+pt2.dist(pt3) < pt1.dist(pt3)+pt2.dist(pt4))
  // 	    turn[ki+1] = true;
  // 	}
  //     else
  // 	{
  // 	  // Select unused edge
  // 	  size_t kr, kh;
  // 	  for (kr=0; kr<trim_edgs_.size(); ++kr)
  // 	    {
  // 	      for (kh=0; kh<=ki; ++kh)
  // 		if (seq_ix[kh] == (int)kr)
  // 		  break;
  // 	      if (kh > ki)
  // 		{
  // 		  seq_ix[ki+1] = (int)kr;
  // 		  turn[ki+1] = false;
  // 		  break;
  // 		}
  // 	    }
  // 	}
  //   }

  vector<shared_ptr<ftEdge> > tmp_edgs(trim_edgs_.size());
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    tmp_edgs[ki] = trim_edgs_[seq_ix[ki]];
  std::swap(trim_edgs_, tmp_edgs);
  
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    if (turn[ki])
      trim_edgs_[ki]->reverseGeomCurve();

#ifdef DEBUG_TRIM
  std::ofstream of0("trim_edgs_arrange0.g2");
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      if (!sfcv.get())
	continue;
      shared_ptr<ParamCurve> space = sfcv->spaceCurve();
      if (!space.get())
	continue;
      space->writeStandardHeader(of0);
      space->write(of0);
      Point startpt = space->point(space->startparam());
      of0 << "400 1 0 4 255 0 0 255" << std::endl;
      of0 << "1" << std::endl;
      of0 << startpt << std::endl;
    }
#endif
  
  HedgeSurface *hedge = getSurface(0);
  shared_ptr<ParamSurface> surf = hedge->surface();
  double lenfac = 0.6;
  double tol4 = 4.0*tol;
  for (size_t ki=1; ki<start_ix.size(); ++ki)
    {
      size_t kj, kr;
      for (kj=start_ix[ki-1]; kj<start_ix[ki]; ++kj)
	{
	  kr = (start_ix[ki]-start_ix[ki-1] == 1) ? kj :
	    ((kj == start_ix[ki]-1) ? start_ix[ki-1] : kj+1);

	  // Check distance between adjacent curve segments
	  Point pt1 = trim_edgs_[kj]->point(trim_edgs_[kj]->tMax());
	  Point pt2 = trim_edgs_[kr]->point(trim_edgs_[kr]->tMin());
	  double glen = pt1.dist(pt2);
	  if (glen > tol)
	    {
	      // Define missing edge
	      shared_ptr<CurveOnSurface> sfcv1 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[kj]->geomCurve());
	      shared_ptr<CurveOnSurface> sfcv2 =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[kr]->geomCurve());
	      if ((!sfcv1.get()) || (!sfcv2.get()))
		continue;

	      double clen = sfcv1->estimatedCurveLength();
	      if (kj == kr && (!trim_edgs_[kj]->twin()))
		{
		  // Single curve
		  if (glen > lenfac*clen)
		    {
		      // Remove trimming curve
		      trim_edgs_.erase(trim_edgs_.begin()+kj);
		      for (size_t kh=ki+1; kh<start_ix.size(); ++kh)
			start_ix[kh]--;
		      start_ix.erase(start_ix.begin()+ki);
		      --ki;
		      break;
		    }
		  
		}

	      if (clen - glen < tol4)
		continue;   // Don't duplicate an almost straight curve
	      
	      shared_ptr<ParamCurve> pcv1 = sfcv1->parameterCurve();
	      shared_ptr<ParamCurve> pcv2 = sfcv2->parameterCurve();
	      if ((!pcv1.get()) || (!pcv2.get()))
		continue;

	      Point ppos1 = pcv1->point(pcv1->endparam());
	      Point ppos2 = pcv2->point(pcv2->startparam());
	      if (ppos1.dist(ppos2) < eps)
		continue;
	      shared_ptr<SplineCurve> gap_par(new SplineCurve(ppos1, ppos2));
	      shared_ptr<CurveOnSurface> gap_sfcv(new CurveOnSurface(surf,
								     gap_par,
								     true));
	      gap_sfcv->ensureSpaceCrvExistence(tol);
	      shared_ptr<ftEdge> gap_edge(new ftEdge(hedge, gap_sfcv,
						     gap_sfcv->startparam(),
						     gap_sfcv->endparam()));
	      trim_edgs_.insert(trim_edgs_.begin()+kr, gap_edge);

	      for (size_t kh=ki; kh<start_ix.size(); ++kh)
		start_ix[kh]++;
	      ++kj;
	    }
	}
    }
  
#ifdef DEBUG_TRIM
  std::ofstream of1("trim_edgs_arrange.g2");
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      if (!sfcv.get())
	continue;
      shared_ptr<ParamCurve> space = sfcv->spaceCurve();
      if (!space.get())
	continue;
      space->writeStandardHeader(of1);
      space->write(of1);
      Point startpt = space->point(space->startparam());
      of1 << "400 1 0 4 255 0 0 255" << std::endl;
      of1 << "1" << std::endl;
      of1 << startpt << std::endl;
    }
#endif
  
  int stop_break = 1;
  return true;
}

//===========================================================================
void RevEngRegion::extractPointsAtSeam(vector<RevEngPoint*>& seam_pts1,
				       vector<RevEngPoint*>& seam_pts2, bool along_udir)
//===========================================================================
{
  // Investigate only the points close to the domain boundary
  int p_ix = along_udir ? 0 : 1;
  double start = domain_[2*p_ix];
  double end = start + 0.1*(domain_[2*p_ix+1] - domain_[2*p_ix]);
  double end2 = domain_[2*p_ix+1];
  double start2 = end2 - 0.1*(domain_[2*p_ix+1] - domain_[2*p_ix]);
  std::set<RevEngPoint*> tmp_pts2;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector2D uv = group_points_[ki]->getPar();
      if (uv[p_ix] >= start && uv[p_ix] <= end)
	{
	  bool found = false;
	  vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	      Vector2D uv2 = curr->getPar();
	      if (uv2[p_ix] >= start2 && uv2[p_ix] <= end2)
		{
		  tmp_pts2.insert(curr);
		  found = true;
		}
	    }
	  if (found)
	    seam_pts1.push_back(group_points_[ki]);
	}
    }
  if (tmp_pts2.size() > 0)
    {
      vector<RevEngPoint*> tmp2_pts2(tmp_pts2.begin(), tmp_pts2.end());
      seam_pts2 = tmp2_pts2;
    }
}

//===========================================================================
void RevEngRegion::getAdjCandMerge(vector<RevEngRegion*>& adj_surf,
				   vector<RevEngRegion*>& adj_nosurf)
//===========================================================================
{
  vector<RevEngRegion*> adj_reg;
  adj_reg.insert(adj_reg.end(), adjacent_regions_.begin(), adjacent_regions_.end());
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (adj_reg[kj]->numPoints() > adj_reg[ki]->numPoints())
	std::swap(adj_reg[ki], adj_reg[kj]);
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    
    {
      if (adj_reg[ki]->prev_region_ && adj_reg[ki]->prev_region_ == this)
	continue;
      if (commonRevEdge(adj_reg[ki]) || adj_reg[ki]->hasAssociatedBlend() ||
	  adj_reg[ki]->hasBlendEdge())
	continue;

      if (adj_reg[ki]->hasSurface())
	{
	  // Check compatibility (todo)
	  adj_surf.push_back(adj_reg[ki]);
	}
      else
	{
	  // Check compatibility (todo)
	  adj_nosurf.push_back(adj_reg[ki]);
	}
    }
}

//===========================================================================
bool RevEngRegion::commonRevEdge(RevEngRegion *other)
//===========================================================================
{
  vector<RevEngEdge*> other_edgs = other->getAllRevEdges();
  for (size_t ki=0; ki<rev_edges_.size(); ++ki)
    for (size_t kj=0; kj<other_edgs.size(); ++kj)
      if (rev_edges_[ki] == other_edgs[kj])
	return true;

  return false;
}

//===========================================================================
bool RevEngRegion::commonTrimEdge(RevEngRegion *other)
//===========================================================================
{
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      ftEdgeBase *twin0 = trim_edgs_[ki]->twin();
      if (!twin0)
	continue;
      ftEdge *twin = twin0->geomEdge();
      ftFaceBase *face0 = twin->face();
      if (!face0)
	continue;
      HedgeSurface *hedge = dynamic_cast<HedgeSurface*>(face0);
      if (!hedge)
	continue;
      RevEngRegion *reg = hedge->getRegion(0);
      if (reg == other)
	return true;
    }
  return false;
}

//===========================================================================
void RevEngRegion::adaptEdges()
//===========================================================================
{
  if (associated_sf_.size() > 0)
    {
      RectDomain dom = associated_sf_[0]->surface()->containingDomain();
      double dom2[4];
      dom2[0] = dom.umin();
      dom2[1] = dom.umax();
      dom2[2] = dom.vmin();
      dom2[3] = dom.vmax();
      for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	{
	  adaptOneEdge(trim_edgs_[ki], dom2);
	}
    }
}

//===========================================================================
void RevEngRegion::adaptOneEdge(shared_ptr<ftEdge>& edge, double dom[4])
//===========================================================================
{
  int status = 0;
  double eps = 1.0e-9;
  shared_ptr<CurveOnSurface> sfcv =
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(edge->geomCurve());
  if (!sfcv.get())
    return;
  if (sfcv->isConstantCurve())
    {
      int dir, bd;
      double val;
      bool orient;
      sfcv->getConstantCurveInfo(dir, val, bd, orient);
      shared_ptr<ParamCurve> pcurve = sfcv->parameterCurve();
      if (pcurve.get())
	{
	  double t1 = pcurve->startparam();
	  double t2 = pcurve->endparam();
	  Point pt1 = pcurve->point(t1);
	  Point pt2 = pcurve->point(t2);
	  double tpar, dist;
	  Point close;
	  int ix = (dir == 1) ? 1 : 0;
	  double tmax = std::max(pt1[ix], pt2[ix]);
	  double tmin = std::min(pt1[ix], pt2[ix]);
	  if (tmax > dom[2*ix+1]+eps)
	    {
	      Point pos;
	      pos = (ix == 0) ? Point(dom[2*ix+1], val) : Point(val, dom[2*ix+1]);
	      pcurve->closestPoint(pos, t1, t2, tpar, close, dist);
	      t2 = tpar;
	    }
	  if (tmin < dom[2*ix]-eps)
	    {
	      Point pos;
	      pos = (ix == 0) ? Point(dom[2*ix], val) : Point(val, dom[2*ix]);
	      pcurve->closestPoint(pos, t1, t2, tpar, close, dist);
	      t1 = tpar;
	    }

	  if (t1 > pcurve->startparam() || t2 < pcurve->endparam())
	    {
	      shared_ptr<CurveOnSurface> sub(sfcv->subCurve(t1, t2));
	      shared_ptr<ftEdge> subedge(new ftEdge(edge->face(), sub, t1, t2));
	      ftEdgeBase *twin = edge->twin();
	      edge->disconnectTwin();
	      subedge->connectTwin(twin, status);
	      edge = subedge;
	    }
	}
      else
	MESSAGE("adaptCurve, missing parameter curve");
    }
  else
    MESSAGE("adaptCurve for non-constant trimming curves is not implemented");
}

//===========================================================================
void RevEngRegion::getAdjacentBlends(vector<RevEngRegion*>& adj_blends)
//===========================================================================
{
  for (auto it=adjacent_regions_.begin(); it!= adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasBlendEdge())
	adj_blends.push_back(*it);
    }
}

//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::updateFreeform(vector<RevEngPoint*>& points,
						       double tol)
//===========================================================================
{
  shared_ptr<SplineSurface> dummy;
  if (associated_sf_.size() == 0)
    return dummy;

  shared_ptr<ParamSurface> surf0 = associated_sf_[0]->surface();
  shared_ptr<SplineSurface> surf =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(surf0);
  if (!surf.get())
    return dummy;

  vector<double> data, param;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Vector2D uv = points[ki]->getPar();
      data.insert(data.end(), xyz.begin(), xyz.end());
      param.insert(param.end(), uv.begin(), uv.end());
    }

  int dim = surf->dimension();
  ApproxSurf approx(surf, data, param, dim, tol);
  //ApproxSurf approx(surf3, data, param, dim, tol, 0, false, true, 0, true);
  approx.setMBA(true);
  approx.setFixBoundary(false);
  int max_iter = 1;
  double maxd, avd;
  int num_out;
  shared_ptr<SplineSurface> surf2;
  try {
    surf2 = approx.getApproxSurf(maxd, avd, num_out, max_iter);
  }
  catch (...)
  {
    std::cout << "Surface update failed" << std::endl;
  }
  
 return surf2;
}
  

//===========================================================================
shared_ptr<SplineSurface> RevEngRegion::computeFreeform(vector<RevEngPoint*>& points,
							double tol)
//===========================================================================
{
// Parameterize
  vector<double> data, param;
  int inner1=0, inner2=0;
      
  // Compute PCA axes
  double lambda[3];
  Point eigen1, eigen2, eigen3;
  getPCA(lambda, eigen1, eigen2, eigen3);

  //bool usebasesf = false;
  bool done = false;
  bool close1 = false, close2 = false;
  if (associated_sf_.size() > 0)
    {
      done = parameterizeOnSurf(group_points_,
				associated_sf_[0]->surface(), data,
				param, inner1, inner2, close1, close2);
    }
  if (!done && basesf_.get() && avdist_base_ <= 5.0*tol)
    {
      done = parameterizeOnSurf(group_points_, basesf_, data, param, 
				inner1, inner2, close1, close2);
    }
  if (!done)
    {
     // Parameterize on plane
      RevEngUtils::parameterizeWithPlane(points, bbox_, eigen1,
					 eigen2, data, param);
    }
#ifdef DEBUG_EXTRACT
  std::ofstream ofpar("parpoints.g2");
  int nmbpar = (int)param.size()/2;
  ofpar << "400 1 0 0" << std::endl;
  ofpar << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar << param[2*ka] << " " << param[2*ka+1] << "  0.0" << std::endl;
#endif
  
  vector<double> param2;
  double umin, umax, vmin, vmax;
  bool repar = reparameterize(param, param2, umin, umax, vmin, vmax);
  //std::cout << "repar: " << repar << std::endl;
  if (repar)
    {
      vector<double> p_edg1, p_edg2;
      for (size_t ki=0; ki<points.size(); ++ki)
	{
	  vector<ftSamplePoint*> next = points[ki]->getNeighbours();
	  for (size_t kj=0; kj<next.size(); ++kj)
	    {
	      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	      if (pt->region() != this)
		continue;
	      std::vector<RevEngPoint*>::iterator it = std::find(points.begin(),
								 points.end(), pt);
	      if (it != points.end())
		{
		  int ix = (int)(it - points.begin());
		  p_edg1.push_back(param[2*ki]);
		  p_edg1.push_back(param[2*ki+1]);
		  p_edg2.push_back(param2[2*ki]);
		  p_edg2.push_back(param2[2*ki+1]);
		  p_edg1.push_back(param[2*ix]);
		  p_edg1.push_back(param[2*ix+1]);
		  p_edg2.push_back(param2[2*ix]);
		  p_edg2.push_back(param2[2*ix+1]);
		}
	    }
	}

#ifdef DEBUG_EXTRACT
      std::ofstream ofe1("par_edgs1.g2");
      std::ofstream ofe2("par_edgs2.g2");
      ofe1 << "410 1 0 4 255 0 0 255" << std::endl;
      ofe2 << "410 1 0 4 255 0 0 255" << std::endl;
      ofe1 << p_edg1.size()/4 << std::endl;
      ofe2 << p_edg2.size()/4 << std::endl;
      for (size_t ki=0; ki<p_edg1.size(); ki+=4)
	{
	  ofe1 << p_edg1[ki] << " " << p_edg1[ki+1] << " 0.0 " << p_edg1[ki+2];
	  ofe1 << " " << p_edg1[ki+3] << " 0.0" << std::endl;
	  ofe2 << p_edg2[ki] << " " << p_edg2[ki+1] << " 0.0 " << p_edg2[ki+2];
	  ofe2 << " " << p_edg2[ki+3] << " 0.0" << std::endl;
	}
#endif
      
      double plenfac = 5.0;
      if (umax - umin > plenfac*(vmax - vmin))
	{
	  if (inner1 <= inner2)
	    inner1 = inner2 + 1;
	}
      else if (vmax - vmin > plenfac*(umax - umin))
	{
	  if (inner2 <= inner1)
	    inner2 = inner1 + 1;
	}
      std::swap(param, param2);
    }

  // Extend with extra points in corners far from the point cloud
  //size_t nmb_prev_extend = param.size();
  if ((!close1) && (!close2))
    {
      try {
	extendInCorner(data, param, umin, umax, vmin, vmax);
      }
      catch (...)
	{
	  std::cout << "Corner extend failed" << std::endl;
	}
    }
  
  // Approximate
  double maxd, avd;
  int num_out;
  int max_iter = (associated_sf_.size() > 0) ? 3 : 6; //2;
  int dim = bbox_.low().dimension();
  int order = 4;
  shared_ptr<SplineSurface> surf;
  double del = 0.01;
  vector<double> parvals;
  try {
    surf = RevEngUtils::surfApprox(data, dim, param, order, order, 
				   order+inner1,  order+inner2,
				   close1, close2, max_iter, 
				   tol, maxd, avd, num_out, parvals, del);
  }
  catch (...)
    {
#ifdef DEBUG_EXTRACT      
       std::cout << "Surface approximation failed" << std::endl;
#endif
   }

#ifdef DEBUG_EXTRACT
  if (surf.get())
    {
      std::ofstream of("spline_approx.g2");
      surf->writeStandardHeader(of);
      surf->write(of);
    }
  else
    int no_surf = 1;
#endif
  
  return surf;
}
											 
int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//===========================================================================
bool RevEngRegion::reparameterize(vector<double>& param, vector<double>& param2,
				  double& umin, double& umax, double& vmin, double& vmax)
//===========================================================================
{
  // Define raster
  int nmb_div = 15;
  vector<vector<int> > raster;
  vector<double> param_copy(param.begin(), param.end());  // Raster definition
  // changes sequence of parameter points
  defineRaster(param_copy, nmb_div, raster, umin, umax, vmin, vmax);
  int div2 = (int)raster.size();
  if (div2 == 0)
    return false;
  int div1 = (int)raster[0].size();
  double udel = (umax - umin)/(double)(div1);
  double vdel = (vmax - vmin)/(double)(div2);

  // Count number of empty cells
  int nmb_zero = 0;
  for (int kb=0; kb<div2; ++kb)
    for (int ka=0; ka<div1; ++ka)
      if (raster[kb][ka] == 0)
	++nmb_zero;

  if (nmb_zero < (int)(0.5*div1*div2))
    return false;

  // Identify points on edges
  double eps = 1.0e-10;
  vector<double> edge_pts;
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      if (param[kr]-umin < eps || umax-param[kr] < eps)
	{
	  edge_pts.push_back((param[kr]-umin < eps) ? 0.0 : (double)div1);
	  edge_pts.push_back(0.5*(int)(2*(param[kr+1]-vmin)/vdel));
	}
      if (param[kr+1]-vmin < eps || vmax-param[kr+1] < eps)
	{
	  edge_pts.push_back(0.5*(int)(2*(param[kr]-umin)/udel));
	  edge_pts.push_back((param[kr+1]-vmin < eps) ? 0.0 : (double)div2);
	}
    }

  // Find candidates for endpoints of a mid curve through the paramaeter points
  vector<Point> cand_par;
  vector<double> wd;
  int i1, i2, i3, i4, ix, kc;
  for (kc=0, ix=0; kc<2; ++kc, ix=div2-1)
    for (i1=i2=0; i1<div1; i1=i2)
      {
	for (; i1<div1; ++i1)
	  if (raster[ix][i1] > 0)
	    break;
	for (i2=i1; i2<div1; ++i2)
	  if (raster[ix][i2] == 0)
	    break;
	if (i2 > i1)
	  {
	    size_t kr;
	    for (kr=0; kr<edge_pts.size(); kr+=2)
	      {
		if (((kc == 0 && edge_pts[kr+1] == 0) ||
		     (kc == 1 && edge_pts[kr+1] == div2)) &&
		    (edge_pts[kr] >= i1 && edge_pts[kr] <= i2))
		  break;
	      }
	    
	    if (kr < edge_pts.size())
	      {
		Point currpar(2);
		currpar[0] = edge_pts[kr]; //0.5*(i1+i2);
		currpar[1] = (kc == 0) ? 0 : (double)div2;
		cand_par.push_back(currpar);
		wd.push_back(0.5*udel*(i2-i1));
	      }
	  }
      }

  for (kc=0, ix=0; kc<2; ++kc, ix=div1-1)
    for (i1=i2=0; i1<div2; i1=i2)
      {
	for (; i1<div2; ++i1)
	  if (raster[i1][ix] > 0)
	    break;
	for (i2=i1; i2<div2; ++i2)
	  if (raster[i2][ix] == 0)
	    break;
	if (i2 > i1)
	  {
	    size_t kr;
	    for (kr=0; kr<edge_pts.size(); kr+=2)
	      {
		if (((kc == 0 && edge_pts[kr] == 0) ||
		     (kc == 1 && edge_pts[kr] == div1)) &&
		    (edge_pts[kr+1] >= i1 && edge_pts[kr+1] <= i2))
		  break;
	      }
	    
	    if (kr < edge_pts.size())
	      {
		Point currpar(2);
		currpar[0] = (kc == 0) ? 0 : (double)div1;
		currpar[1] = edge_pts[kr+1]; //0.5*(i1+i2);
		cand_par.push_back(currpar);
		wd.push_back(0.5*vdel*(i2-i1));
	      }
	  }
      }

  if (cand_par.size() < 2)
    return false;
  
  // Decide on endpoints
  int ix1=-1, ix2=-1;
  double maxd = std::numeric_limits<double>::lowest();
  for (size_t kr=0; kr<cand_par.size(); ++kr)
    for (size_t kh=kr+1; kh<cand_par.size(); ++kh)
      {
	double dd = cand_par[kr].dist(cand_par[kh]);
	dd -= (wd[kr] + wd[kh]);
	if (dd > maxd)
	  {
	    ix1 = (int)kr;
	    ix2 = (int)kh;
	    maxd = dd;
	  }
      }

  // Compute points along a sceleton in the parameter point cloud
  double curr1[2], curr2[2];
  double wd1, wd2;
  vector<double> ptpar;
  vector<double> width;
  ptpar.push_back(cand_par[ix1][0]);
  ptpar.push_back(cand_par[ix1][1]);
  width.push_back(wd[ix1]);

  Point prev;
  int sgn1=0, sgn2=0;
  bool loop = false;
  for (kc=0; kc<(int)ptpar.size(); kc+=2)
    {
#ifdef DEBUG_REPAR
      std::ofstream ofc("midpol.g2");
      ofc << "410 1 0 4 0 0 0 255" << std::endl;
      ofc << ptpar.size()/2 << std::endl;
      ofc << umin+udel*ptpar[0] << " " << vmin+vdel*ptpar[1] << "  0.0" << std::endl;
      for (size_t kr=2; kr<ptpar.size(); kr+=2)
	{
	  ofc << umin+udel*ptpar[kr] << " " << vmin+vdel*ptpar[kr+1] << "  0.0" << std::endl;
	  ofc << umin+udel*ptpar[kr] << " " << vmin+vdel*ptpar[kr+1] << "  0.0" << std::endl;
	}
      ofc << umin+udel*cand_par[ix2][0] << " " << vmin+vdel*cand_par[ix2][1] << " 0.0" << std::endl;
#endif
      curr1[0] = curr2[0] = ptpar[kc];
      curr1[1] = curr2[1] = ptpar[kc+1];
      if (sgn1 == 0 || curr1[0] == prev[0])
	sgn1 = (cand_par[ix2][0] > curr1[0]) ? 1 : -1;
      if (sgn2 == 0 || curr2[1] == prev[1])
	sgn2 = (cand_par[ix2][1] > curr2[1]) ? 1 : -1;
      double fac1 = (curr1[0] == 0.0) ? 0.5 : 1.0;
      double fac2 = (curr2[1] == 0.0) ? 0.5 : 1.0;
      curr1[0] += fac1*sgn1;
      curr2[1] += fac2*sgn2;
      prev = Point(ptpar[kc], ptpar[kc+1]);
      if (prev.dist(cand_par[ix2]) <= 1.0)
	break;
      curr1[0] = std::min(curr1[0], div1-0.5);
      curr1[1] = std::min(curr1[1], div2-0.5);
      curr2[0] = std::min(curr2[0], div1-0.5);
      curr2[1] = std::min(curr2[1], div2-0.5);
      curr1[0] = std::max(curr1[0], 0.0);
      curr1[1] = std::max(curr1[1], 0.0);
      curr2[0] = std::max(curr2[0], 0.0);
      curr2[1] = std::max(curr2[1], 0.0);

      // Extent with constant first parameter direction
      getParExtent(curr1, 0, raster, i1, i2);
      curr1[1] = 0.5*(double)(i1+i2+1);
      wd1 = 0.5*(i1 - i2);

      // Constant second parameter direction
      getParExtent(curr2, 1, raster, i3, i4);
      curr2[0] = 0.5*(double)(i3+i4+1);
      wd2 = 0.5*(i3 - i4);

      // Select continuation
      Point lastpar(ptpar[ptpar.size()-2], ptpar[ptpar.size()-1]);
      Point c1(curr1[0], curr1[1]);
      Point c2(curr2[0], curr2[1]);
      if (i1-i2 < i3-i4 || c2.dist(lastpar) < eps) 
	{
	  ptpar.push_back(curr1[0]);
	  ptpar.push_back(curr1[1]);
	  width.push_back(udel*wd1);
	  lastpar = c1;
	}
      else
	{
	  ptpar.push_back(curr2[0]);
	  ptpar.push_back(curr2[1]);
	  width.push_back(vdel*wd2);
	  lastpar = c2;
	}
	
      // Check for loops
      for (size_t kh=2; kh<ptpar.size()-2; kh+=2)
	{
	  Point currpar(ptpar[kh], ptpar[kh+1]);
	  if (currpar.dist(lastpar) < eps)
	    {
	      loop = true;
	      break;
	    }
	}
      if (loop)
	break;
      
      int stop_break0 = 1;      
    }
  if (loop)
    return false;
  
  Point last(ptpar[ptpar.size()-2], ptpar[ptpar.size()-1]);
  if (last.dist(cand_par[ix2]) > 0.0)
    {
      ptpar.push_back(cand_par[ix2][0]);
      ptpar.push_back(cand_par[ix2][1]);
      width.push_back(wd[ix2]);
    }

  if (ptpar.size() <= 2)
    return false;
  
  // Translate to real coordinates
  vector<double> ptpar2(ptpar.size());
  for (size_t kr=0; kr<ptpar.size(); kr+=2)
    {
      ptpar2[kr] = umin+udel*ptpar[kr];
      ptpar2[kr+1] = vmin+vdel*ptpar[kr+1];
    }

  // Extend
  double npar1[2], npar2[2];
  for (int ka=0; ka<2; ++ka)
    {
      npar1[ka] = ptpar2[ka] - 0.1*(ptpar2[2+ka]-ptpar2[ka]);
      npar2[ka] = ptpar2[ptpar2.size()-2+ka] + 0.1*(ptpar2[ptpar2.size()-2+ka]-
						   ptpar2[ptpar2.size()-4+ka]);
    }
  ptpar2.insert(ptpar2.begin(), &npar1[0], &npar1[0]+2);
  ptpar2.insert(ptpar2.end(), &npar2[0], &npar2[0]+2);
  width.insert(width.begin(), width[0]);
  width.push_back(width[width.size()-1]);
  

#ifdef DEBUG_REPAR
  std::ofstream ofp("parpt.g2");
  ofp << "400 1 0 4 200 55 0 255" << std::endl;
  ofp << ptpar2.size()/2 << std::endl;
  for (size_t kr=0; kr<ptpar2.size(); kr+=2)
    ofp << ptpar2[kr] << " " << ptpar2[kr+1] << " 0.0" << std::endl;
#endif
  
  // Define "mid" curve
  vector<double> par2(ptpar2.size()/2);
  par2[0] = 0.0;
  for (size_t kr=2; kr<ptpar2.size(); kr+=2)
    {
      double dd = Utils::distance_squared(&ptpar2[kr-2], &ptpar2[kr], &ptpar2[kr]);
      par2[kr/2] = par2[(kr-2)/2] + sqrt(dd);
    }

  int in = 4, ik = 4;
  int dim = 2;
  double tol = 0.1;
  ApproxCurve approx(ptpar2, par2, dim, tol, in, ik);
  
  double maxdist, avdist;
  shared_ptr<SplineCurve> midcv = approx.getApproxCurve(maxdist, avdist, 1);

  ApproxCurve approx2(width, par2, 1, tol, in, ik);
  
  double maxdistw, avdistw;
  shared_ptr<SplineCurve> widthcv = approx2.getApproxCurve(maxdistw, avdistw, 1);

#ifdef DEBUG_REPAR
  std::ofstream ofmid("midcv.g2");
  midcv->writeStandardHeader(ofmid);
  midcv->write(ofmid);
#endif
  
  double avw = 0.0;
  for (size_t kr=0; kr<width.size(); ++kr)
    avw += width[kr];
  avw /= (double)width.size();

  param2.resize(param.size());
  double tmin = midcv->startparam();
  double tmax = midcv->endparam();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      Point currpar(param[kr], param[kr+1]);
      double tpar, dist;
      Point close;
      midcv->closestPoint(currpar, tmin, tmax, tpar, close, dist);

      vector<Point> der(2);
      midcv->point(der, tpar, 1);
      Point vec = close - currpar;
      Point vec2(-vec[1], vec[0]);
      //double tmp = fabs(der[1]*vec);
      // if (tmp > 1.0e-6)
      // 	std::cout << "inner: " << tmp << std::endl;
      if (vec2*der[1] < 0.0)
	dist *= -1.0;
      Point wdt = widthcv->ParamCurve::point(tpar);
      //dist *= (avw/wdt[0]);
      param2[kr] = tpar;
      param2[kr+1] = dist;
      int stop_break0 = 1;
    }
  
#ifdef DEBUG_REPAR
  std::ofstream ofpar2("parpoints_repar.g2");
  int nmbpar = (int)param2.size()/2;
  ofpar2 << "400 1 0 0" << std::endl;
  ofpar2 << nmbpar << std::endl;
  for (int ka=0; ka<nmbpar; ++ka)
    ofpar2 << param2[2*ka] << " " << param2[2*ka+1] << "  0.0" << std::endl;
#endif
  
  vector<vector<int> > raster2;
  double umin2, umax2, vmin2, vmax2;
  vector<double> param2_copy(param2.begin(), param2.end());  // Raster definition
  // changes sequence of parameter points
  defineRaster(param2_copy, nmb_div, raster2, umin2, umax2, vmin2, vmax2);

  // Count number of empty cells
  int nmb_zero2 = 0;
  for (size_t kr=0; kr<raster2.size(); ++kr)
    for (size_t kh=0; kh<raster2[kr].size(); ++kh)
      if (raster2[kr][kh] == 0)
	++nmb_zero2;

  if (nmb_zero2 >= nmb_zero)
    return false;

  umin = umin2;
  umax = umax2;
  vmin = vmin2;
  vmax = vmax2;
  return true;
}

//===========================================================================
void RevEngRegion::getParExtent(double curr[2], int pdir, vector<vector<int> >& raster,
				int& i1, int& i2)
//===========================================================================
{
  int div = (pdir == 0) ? (int)raster.size() : (int)raster[0].size();
  int j1 = (int)curr[pdir];
  int pdir2 = 1 - pdir;
  
  for (i1=(int)curr[pdir2]; i1<div; ++i1)
    {
      int r1 = (pdir == 0) ? raster[i1][j1] : raster[j1][i1];
      if (r1 > 0)
	break;
    }
  if (i1 == div)
    i1=(int)curr[pdir2];
  for (; i1<div; ++i1)
    {
      int r1 = (pdir == 0) ? raster[i1][j1] : raster[j1][i1];
      if (r1 == 0)
	break;
    }
  for (i2=(int)curr[pdir2]; i2>=0; --i2)
    {
      int r1 = (pdir == 0) ? raster[i2][j1] : raster[j1][i2];
      if (r1 > 0)
	break;
    }
  if (i2 < 0)
    i2=(int)curr[pdir2];
  for (; i2>=0; --i2)
    {
      int r1 = (pdir == 0) ? raster[i2][j1] : raster[j1][i2];
      if (r1 == 0)
      break;
    }
  }


//===========================================================================
void RevEngRegion::defineRaster(vector<double>& param, int nmb_div,
				vector<vector<int> >& raster, double& umin,
				double& umax, double& vmin, double& vmax)
//===========================================================================
{
  umin = std::numeric_limits<double>::max();
  umax = std::numeric_limits<double>::lowest();
  vmin = std::numeric_limits<double>::max();
  vmax = std::numeric_limits<double>::lowest();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      umin = std::min(umin, param[kr]);
      umax = std::max(umax, param[kr]);
      vmin = std::min(vmin, param[kr+1]);
      vmax = std::max(vmax, param[kr+1]);
    }
  
  int nm = nmb_div*nmb_div;
  double dom = (umax-umin)*(vmax-vmin);
  double c1 = std::pow((double)nm/dom, 1.0/2.0);
  int div1, div2;
  div1 = (int)(c1*(umax-umin));
  ++div1;
  div2 = (int)(c1*(vmax-vmin));
  ++div2;
  double udel = (umax - umin)/(double)(div1);
  double vdel = (vmax - vmin)/(double)(div2);

#ifdef DEBUG_REPAR
  std::ofstream of("division_lines.g2");
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << div1+div2+2 << std::endl;
  for (int ki=0; ki<=div1; ++ki)
    {
      Point p1(umin+ki*udel, vmin, 0.0);
      Point p2(umin+ki*udel, vmax, 0.0);
      of << p1 << " " << p2 << std::endl;
    }
  for (int ki=0; ki<=div2; ++ki)
    {
      Point p1(umin, vmin+ki*vdel, 0.0);
      Point p2(umax, vmin+ki*vdel, 0.0);
      of << p1 << " " << p2 << std::endl;
    }
#endif
  
  raster.resize(div2);
  for (size_t kr=0; kr<raster.size(); ++kr)
    {
      raster[kr].resize(div1, 0);
    }

  int nmb_pts = (int)param.size()/2;
  qsort(&param[0], nmb_pts, 2*sizeof(double), compare_v_par);
  int pp0, pp1;
  int ka, kb;
  double upar, vpar;
  for (vpar=vmin+vdel, pp0=0, kb=0; kb<div2; ++kb, vpar+=vdel)
    {
      for (pp1=pp0; pp1<2*nmb_pts && param[pp1+1]<=vpar; pp1+=2);
      qsort(&param[pp0], (pp1-pp0)/2, 2*sizeof(double), compare_u_par);

      int pp2, pp3;
      for (upar=umin+udel, pp2=pp0, ka=0; ka<div1; ++ka, upar+=udel)
	{
	  for (pp3=pp2; pp3<pp1 && param[pp3]<=upar; pp3+=2);
	  raster[kb][ka] = (pp3-pp2)/2;
	  pp2 = pp3;
	}
      pp0 = pp1;
    }
}

//===========================================================================
void RevEngRegion::extendInCorner(vector<double>& data, vector<double>& param,
				  double umin, double umax, double vmin, double vmax)
//===========================================================================
{
  Point corner[4];
  corner[0] = Point(umin, vmin);
  corner[1] = Point(umax, vmin);
  corner[2] = Point(umin, vmax);
  corner[3] = Point(umax, vmax);
  double cdist[4];
  cdist[0] = cdist[1] = cdist[2] = cdist[3] = std::numeric_limits<double>::max();
  for (size_t kr=0; kr<param.size(); kr+=2)
    {
      Point currpar(param[kr], param[kr+1]);
      for (int ka=0; ka<4; ++ka)
	{
	  double dd = corner[ka].dist(currpar);
	  cdist[ka] = std::min(cdist[ka], dd);
	}
    }

  int div = 15;
  double minc = std::min((umax-umin)/(double)div, (vmax-vmin)/(double)div);
  vector<double> corner_data;
  vector<double> corner_par;
  int min_nmb = 20;
  for (int ka=0; ka<4; ++ka)
    {
      if (cdist[ka] > minc)
      {
	// Estimate corner point
	// Collect points in the vicinity of the corner
	vector<double> data2;
	vector<double> param2;
	double umin2 = std::numeric_limits<double>::max();
	double umax2 = std::numeric_limits<double>::lowest();
	double vmin2 = std::numeric_limits<double>::max();
	double vmax2 = std::numeric_limits<double>::lowest();
	double lim = 1.1*cdist[ka];  //2.0
	for (size_t kr=0; kr<param.size()/2; ++kr)
	  {
	    Point currpar(param[2*kr], param[2*kr+1]);
	    if (currpar.dist(corner[ka]) < lim)
	      {
		data2.insert(data2.end(), data.begin()+3*kr, data.begin()+3*(kr+1));
		param2.insert(param2.end(), param.begin()+2*kr, param.begin()+2*(kr+1));
		umin2 = std::min(umin2, param[2*kr]);
		umax2 = std::max(umax2, param[2*kr]);
		vmin2 = std::min(vmin2, param[2*kr+1]);
		vmax2 = std::max(vmax2, param[2*kr+1]);
	      }
	  }

	while ((int)data2.size()/3 < min_nmb)
	  {
	    data2.clear();
	    param2.clear();
	    lim += 0.5*cdist[ka];
	    for (size_t kr=0; kr<param.size()/2; ++kr)
	      {
		Point currpar(param[2*kr], param[2*kr+1]);
		if (currpar.dist(corner[ka]) < lim)
		  {
		    data2.insert(data2.end(), data.begin()+3*kr, data.begin()+3*(kr+1));
		    param2.insert(param2.end(), param.begin()+2*kr, param.begin()+2*(kr+1));
		    umin2 = std::min(umin2, param[2*kr]);
		    umax2 = std::max(umax2, param[2*kr]);
		    vmin2 = std::min(vmin2, param[2*kr+1]);
		    vmax2 = std::max(vmax2, param[2*kr+1]);
		  }
	      }
	  }

	// Approximate sub cloud with a planar surface
	umin2 = std::min(umin2, corner[ka][0]);
	umax2 = std::max(umax2, corner[ka][0]);
	vmin2 = std::min(vmin2, corner[ka][1]);
	vmax2 = std::max(vmax2, corner[ka][1]);
	umin2 -= 0.1*(umax2 - umin2);
	umax2 += 0.1*(umax2 - umin2);
	vmin2 -= 0.1*(vmax2 - vmin2);
	vmax2 += 0.1*(vmax2 - vmin2);

	shared_ptr<SplineSurface> plane =
	  RevEngUtils::surfApprox(data2, 3, param2, 2, 2, 2, 2, umin, umax,
				  vmin, vmax);
	Point pos = plane->ParamSurface::point(corner[ka][0], corner[ka][1]);
	corner_data.insert(corner_data.end(), pos.begin(), pos.end());
	corner_par.insert(corner_par.end(), corner[ka].begin(), corner[ka].end());
      }
    }

  if (corner_data.size() > 0)
    {
      data.insert(data.end(), corner_data.begin(), corner_data.end());
      param.insert(param.end(), corner_par.begin(), corner_par.end());

#ifdef DEBUG_REPAR
      std::ofstream of("added_corners.g2");
      for (size_t ki=0; ki<corner_data.size()/3; ++ki)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  for (int ka=0; ka<3; ++ka)
	    of << corner_data[3*ki+ka] << " ";
	  of << std::endl;
	}
#endif
    }
}


//===========================================================================
bool RevEngRegion::parameterizeOnSurf(shared_ptr<ParamSurface> surf, 
				      vector<double>& data, vector<double>& param,
				      int& inner1, int& inner2, bool& close1, bool& close2)
//===========================================================================
{
  return parameterizeOnSurf(group_points_, surf, data, param, inner1, inner2,
			    close1, close2);
}

//===========================================================================
bool RevEngRegion::parameterizeOnSurf(vector<RevEngPoint*>& points,
				      shared_ptr<ParamSurface> surf, 
				      vector<double>& data, vector<double>& param,
				      int& inner1, int& inner2, bool& close1, bool& close2)
//===========================================================================
{
  ClassType classtype = surf->instanceType();
  if (classtype == Class_Cone || classtype == Class_Torus)
    {
      // Check angle. An almost plane cone is not appropriate for
      // parametrization
      shared_ptr<ParamSurface> sf = surf;
      shared_ptr<BoundedSurface> bdsf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
      if (bdsf.get())
	sf = bdsf->underlyingSurface();
      if (classtype == Class_Cone)
	{
	  shared_ptr<Cone> cone =
	    dynamic_pointer_cast<Cone,ParamSurface>(sf);
	  double angle = cone->getConeAngle();
	  if (angle > 0.3*M_PI)
	    return false;
	}
      else  if (classtype == Class_Torus)
	{
	  shared_ptr<Torus> torus =
	    dynamic_pointer_cast<Torus,ParamSurface>(sf);
	  double rad1 = torus->getMajorRadius();
	  double rad2 = torus->getMinorRadius();
	  if (rad2 > 0.9*rad1)
	    return false;
	}
    }
  
  bool OK = RevEngUtils::parameterizeOnPrimary(points, surf, data, param,
					       inner1, inner2, close1, close2);
  return OK;
}


//===========================================================================
void RevEngRegion::analyseCylinderProperties(Point avvec, double angtol,
					  vector<RevEngPoint*>& in,
					  vector<RevEngPoint*> out)
//===========================================================================
{
  vector<double> midang(group_points_.size());
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point curr = group_points_[ki]->minCurvatureVec();
      double ang = curr.angle(avvec);
      ang = std::min(ang, M_PI-ang);
      midang[ki] = ang;
      if (ang <= angtol)
	in.push_back(group_points_[ki]);
      else
	out.push_back(group_points_[ki]);
    }

  std::sort(midang.begin(), midang.end());
  int stop_break = 1;
}


//===========================================================================
bool RevEngRegion::sortByAxis(vector<Point>& axis, double tol,
			      double axisang, double planeang,
			      vector<vector<RevEngPoint*> >& groups1,
			      vector<vector<RevEngPoint*> >& groups2,
			      vector<RevEngPoint*>& remaining)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  groups1.resize(2*axis.size());
  groups2.resize(axis.size());
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal = group_points_[kr]->getLocFuncNormal();
      Point normal2 = group_points_[kr]->getTriangNormal();
      int min_ix1 = -1, min_ix2 = -1;
      double min_ang1 = pihalf, min_ang2 = pihalf;
      for (int ka=0; ka<(int)axis.size(); ++ka)
	{
	  double ang = axis[ka].angle(normal);
	  double ang2 = axis[ka].angle(normal2);
	  ang = std::min(fabs(pihalf-ang), fabs(pihalf-ang2));
	  if (ang < min_ang2)
	    {
	      min_ang2 = ang;
	      min_ix2 = ka;
	    }
	}
      for (int ka=0; ka<(int)axis.size(); ++ka)
	{
	  double ang = axis[ka].angle(normal);
	  double ang2 = axis[ka].angle(normal2);
	  if (std::min(ang,ang2) < min_ang1)
	    {
	      min_ang1 = std::min(ang,ang2);
	      min_ix1 = 2*ka+1;
	    }
	  else if (std::min(M_PI-ang,M_PI-ang2) < min_ang1)
	    {
	      min_ang1 = std::min(M_PI-ang,M_PI-ang2);
	      min_ix1 = 2*ka;
	    }
	}
      if (min_ang1 > planeang && min_ang2 < min_ang1 &&
	  min_ang2 < axisang && min_ix2 >= 0 )
	groups2[min_ix2].push_back(group_points_[kr]);
      else if (min_ang1 < axisang && min_ix1 >= 0)
	groups1[min_ix1].push_back(group_points_[kr]);
      else
	remaining.push_back(group_points_[kr]);
    }
      
#ifdef DEBUG_SEGMENT
  std::ofstream of1("axis_groups1.g2");
  for (size_t ki=0; ki<groups1.size(); ++ki)
    {
      if (groups1[ki].size() > 0)
	{
	  of1 << "400 1 0 4 255 0 0 255" << std::endl;
	  of1 <<  groups1[ki].size() << std::endl;
	  for (int ka=0; ka<(int)groups1[ki].size(); ++ka)
	    of1 << groups1[ki][ka]->getPoint() << std::endl;
	}
    }
  std::ofstream of2("axis_groups2.g2");
  for (size_t ki=0; ki<groups2.size(); ++ki)
    {
      if (groups2[ki].size() > 0)
	{
	  of2 << "400 1 0 4 255 0 0 255" << std::endl;
	  of2 <<  groups2[ki].size() << std::endl;
	  for (int ka=0; ka<(int)groups2[ki].size(); ++ka)
	    of2 << groups2[ki][ka]->getPoint() << std::endl;
	}
    }
  std::ofstream of3("remaining.g2");
  if (remaining.size() > 0)
    {
      of3 << "400 1 0 4 0 255 0 255" << std::endl;
      of3 <<  remaining.size() << std::endl;
      for (int ka=0; ka<(int)remaining.size(); ++ka)
	of3 << remaining[ka]->getPoint() << std::endl;
    }
#endif
  int num_groups = 0;
  for (size_t ki=0; ki<groups1.size(); ++ki)
    if (groups1.size() > 0)
      num_groups++;
  for (size_t ki=0; ki<groups2.size(); ++ki)
    if (groups2.size() > 0)
      num_groups++;
  if (remaining.size() > 0)
    num_groups++;
    
  return (num_groups > 1);
}

//===========================================================================
void RevEngRegion::removeLowAccuracyPoints(int min_pt_reg, double tol, double angtol,
					    vector<vector<RevEngPoint*> >& added_groups)
//===========================================================================
{
  shared_ptr<ParamSurface> surf= getSurface(0)->surface();
  int code;
  ClassType classtype = getSurface(0)->instanceType(code);
  bool cyllike = (classtype == Class_Cylinder || classtype == Class_Cone);
  double angtol2 = (surfflag_ == PROBABLE_HELIX) ? 0.5*M_PI : 0.25*M_PI;
  
  vector<vector<RevEngPoint*> > out_groups;
  vector<RevEngPoint*> remain;
  vector<pair<double,double> > distang;
  getDistAndAng(distang);
  identifyOutPoints(distang, tol, angtol, angtol2, out_groups,
		    remain);
#ifdef DEBUG_AXIS
  std::ofstream of1("dist_separated.g2");
  of1 << "400 1 0 4 255 0 0 255" << std::endl;
  of1 << remain.size() << std::endl;
  for (size_t kr=0; kr<remain.size(); ++kr)
    of1 << remain[kr]->getPoint() << std::endl;
  for (size_t kj=0; kj<out_groups.size(); ++kj)
    {
      of1 << "400 1 0 4 0 255 0 255" << std::endl;
      of1 << out_groups[kj].size() << std::endl;
      for (size_t kr=0; kr<out_groups[kj].size(); ++kr)
	of1 << out_groups[kj][kr]->getPoint() << std::endl;
    }
#endif

  double maxdist, avdist;
  int num_in, num2_in;
  //int surf_flag2 = NOT_SET;
  if (remain.size() < group_points_.size())
    {
      maxdist = avdist = 0.0;
      num_in = num2_in = 0;
      double frac = 1.0/(double)remain.size();
      for (size_t kj=0; kj<remain.size(); ++kj)
	{
	  double dist, ang;
	  remain[kj]->getSurfaceDist(dist, ang);
	  maxdist = std::max(maxdist, dist);
	  avdist += frac*dist;
	  if (dist <= tol)
	    {
	      num2_in++;
	      if (ang <= angtol)
		num_in++;
	    }
	}
      int surf_flag = defineSfFlag((int)remain.size(), 0, tol,
				   num_in, num2_in, avdist, cyllike);
      if (surf_flag <= surfflag_ && avdist < avdist_)
	{
	  // Extract points
	  std::swap(group_points_, remain);
	  added_groups.insert(added_groups.end(), out_groups.begin(),
			      out_groups.end());

	  updateInfo(tol, angtol);
	  setSurfaceFlag(surf_flag);
	  updateRegionAdjacency();
	  
	  int stop_break = 1;
	}
	  
    }
  
}

//===========================================================================
void RevEngRegion::setPlaneParam(int min_pt_reg, Point mainaxis[3],
				 double tol, double angtol)
//===========================================================================
{
  if (!hasSurface())
    return;
  shared_ptr<ParamSurface> surf_init = getSurface(0)->surface();
  if (surf_init->instanceType() != Class_Plane)
    return;

  shared_ptr<Plane> plane = dynamic_pointer_cast<Plane, ParamSurface>(surf_init);
  int ix = -1;
  double min_ang = M_PI;
  double pihalf = 0.5*M_PI;
  Point dir = plane->direction();
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(dir);
      ang = fabs(pihalf-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  ix = ka;
	}
    }

  Point dir2 = dir.cross(mainaxis[ix]);
  dir2.normalize();
  std::vector<double> rot_mat = GeometryTools::getRotationMatrix(dir2, min_ang);
  Point rotated(3);
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      rotated[ka] += rot_mat[ka*3+kb]*dir[kb];

  bool updated = updateSurfaceWithAxis(min_pt_reg, rotated, mainaxis, -1, tol,
				       angtol, plane->location());
  
  shared_ptr<Plane> plane2(new Plane(plane->location(), plane->direction(), mainaxis[ix]));

  // Reparameterize. The accuracy should stay the same
  double maxd, avd;
  int num_in, num2_in;
  vector<pair<double, double> >  dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  plane, tol, maxd, avd, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);

  int sf_flag = defineSfFlag(0, tol, num_in, num2_in, avd, false);
  for (size_t kh=0; kh<group_points_.size(); ++kh)
    {
      group_points_[kh]->setPar(Vector2D(parvals[2*kh], parvals[2*kh+1]));
      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
    }
  HedgeSurface *hedge = getSurface(0);
  hedge->replaceSurf(plane2);
  setSurfaceFlag(sf_flag);
  updateInfo(tol, angtol);
  surf_adaption_ =  AXIS_ADAPTED;
 
}

//===========================================================================
bool RevEngRegion::updateSurfaceWithAxis(int min_pt_reg, Point adj_axis, 
					 Point mainaxis[3], int ix, double tol, 
					 double angtol, Point pos)
//===========================================================================
{
  bool updated = false;
  if (!hasSurface())
    return false;
  vector<Point> curr_axis;
  if (adj_axis.dimension() != 0)
    curr_axis.push_back(adj_axis);
  if (ix >=0 && ix < 3)
    curr_axis.push_back(mainaxis[ix]);
  if (curr_axis.size() == 0)
    return false;

  shared_ptr<ParamSurface> surf_init = getSurface(0)->surface();
  int code;
  ClassType classtype = getSurface(0)->instanceType(code);
  bool cyllike = (classtype == Class_Cylinder || classtype == Class_Cone);
  //double angtol2 = (surfflag_ == PROBABLE_HELIX) ? 0.5*M_PI : 0.25*M_PI;

  // Ensure updated information
  double maxdist1 = 0.0, avdist1 = 0.0;
  int num_in1 = 0, num2_in1 = 0;
  double frac = 1.0/(double)group_points_.size();
  for (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double dist, ang;
      group_points_[kj]->getSurfaceDist(dist, ang);
      maxdist1 = std::max(maxdist1, dist);
      avdist1 += frac*dist;
      if (dist <= tol)
	{
	  num2_in1++;
	  if (ang <= angtol)
	    num_in1++;
	}
    }
  int sf_flag1 = defineSfFlag(0, tol, num_in1, num2_in1,
			      avdist1, cyllike);
  
  vector<double> maxdist(curr_axis.size());
  vector<double> avdist(curr_axis.size());
  vector<int> num_in(curr_axis.size());
  vector<int> num2_in(curr_axis.size());
  vector<vector<pair<double, double> > > dist_ang(curr_axis.size());
  vector<vector<double> > parvals(curr_axis.size());
  vector<int> sf_flag(curr_axis.size());
  vector<shared_ptr<ParamSurface> > surf(curr_axis.size());
   for (size_t ki=0; ki<curr_axis.size(); ++ki)
    {
      // Update surface
      surf[ki] = surfaceWithAxis(group_points_, curr_axis[ki], pos,
				 mainaxis);

      // Check accuracy
      vector<RevEngPoint*> inpt, outpt;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      surf[ki], tol, maxdist[ki], avdist[ki], num_in[ki],
			      num2_in[ki], inpt, outpt, parvals[ki],
			      dist_ang[ki], angtol);

      sf_flag[ki] = defineSfFlag(0, tol, num_in[ki], num2_in[ki],
				 avdist[ki], cyllike);

    }

   // Prefer axis adapted
   int numpt = numPoints();
   double avfrac = 0.9;
   double fracfrac = 0.9;
   double tol2 = 0.5*tol;
   int alt_ix = (curr_axis.size() > 1) ? 1 : 0;
   if (alt_ix > 0)
     {
       double num_frac1 = (double)(num_in[0]+num2_in[0])/(double)numpt;
       double num_frac2 = (double)(num_in[1]+num2_in[1])/(double)numpt;
       if ((avdist[0] <= avfrac*avdist[1] || avdist[0] <= tol2) &&
	   fracfrac*num_frac1 >= num_frac2)
	 alt_ix = 0;
     }

   //double num_frac = (double)(num_in1+num2_in1)/(double)numpt;
   //double num_frac_alt = (double)(num_in[alt_ix]+num2_in[alt_ix])/(double)numpt;
   // if ((avfrac*avdist[alt_ix] <= avdist1 || avdist[alt_ix] <= tol2) &&
   //     num_frac_alt >= fracfrac*num_frac)
   if (sf_flag[alt_ix] < ACCURACY_POOR || sf_flag[alt_ix] <= sf_flag1)
     {
       // Replace
       setBaseSf(surf_init, maxdist1, avdist1, num_in1, num2_in1);
       for (size_t kh=0; kh<group_points_.size(); ++kh)
	 {
	   group_points_[kh]->setPar(Vector2D(parvals[alt_ix][2*kh],
					      parvals[alt_ix][2*kh+1]));
	   group_points_[kh]->setSurfaceDist(dist_ang[alt_ix][kh].first,
					     dist_ang[alt_ix][kh].second);
	 }
       HedgeSurface *hedge = getSurface(0);
       hedge->replaceSurf(surf[alt_ix]);
       if (!surf[alt_ix]->isBounded())
	 {
	   double diag = bbox_.low().dist(bbox_.high());
	   hedge->limitSurf(2*diag);
	 }
       setSurfaceFlag(sf_flag[alt_ix]);
       updateInfo(tol, angtol);
       for (size_t kj=0; kj<rev_edges_.size(); ++kj)
	 {
	   rev_edges_[kj]->replaceSurf(this, surf[alt_ix], tol);
	   rev_edges_[kj]->increaseSurfChangeCount();
	 }
       surf_adaption_ = (adj_axis.dimension() != 0 && alt_ix == 0) ?
	 ADJACENT_ADAPTED : AXIS_ADAPTED;
       updated = true;
     }
   else
     {
       // Store alternative as base surface
       setBaseSf(surf[alt_ix], maxdist[alt_ix], avdist[alt_ix],
		 num_in[alt_ix], num2_in[alt_ix]);

       setAccuracy(maxdist1, avdist1, num_in1, num2_in1);
     }
   
   int stop_break = 1;

   return updated;
}

//===========================================================================
void RevEngRegion::updateSurfaceAndInfo(shared_ptr<ParamSurface> surf,
					double tol, double angtol,
					vector<double>& parvals,
					vector<pair<double,double> >& dist_ang,
					vector<RevEngEdge*>& nopar_edgs)
//===========================================================================
{
  if (parvals.size() != 2*group_points_.size() ||
      dist_ang.size() != group_points_.size())
    THROW("RevEngRegion::updateSurfaceAndInfo: Inconsistent input");

    shared_ptr<ParamSurface> surf_init = getSurface(0)->surface();

    // Replace
    setBaseSf(surf_init, maxdist_, avdist_, num_inside_, num_inside2_);
    for (size_t kh=0; kh<group_points_.size(); ++kh)
      {
	group_points_[kh]->setPar(Vector2D(parvals[2*kh],
					   parvals[2*kh+1]));
	group_points_[kh]->setSurfaceDist(dist_ang[kh].first,
					  dist_ang[kh].second);
      }
    HedgeSurface *hedge = getSurface(0);
    hedge->replaceSurf(surf);
    if (!surf->isBounded())
      {
	double diag = bbox_.low().dist(bbox_.high());
	hedge->limitSurf(2*diag);
      }
    updateInfo(tol, angtol);
    bool cyl_like = ((surf->instanceType() == Class_Cylinder ||
		      surf->instanceType() == Class_Cone));
    
    int sf_flag = defineSfFlag((int)group_points_.size(), 0, tol, num_inside_,
			       num_inside2_, avdist_, cyl_like);
    setSurfaceFlag(sf_flag);
    for (size_t kj=0; kj<rev_edges_.size(); ++kj)
      {
	rev_edges_[kj]->replaceSurf(this, surf, tol);
	rev_edges_[kj]->increaseSurfChangeCount();
	int missing = rev_edges_[kj]->missingParCrv();
	if (missing > 0)
	  nopar_edgs.push_back(rev_edges_[kj]);
      }
}

//===========================================================================
void RevEngRegion::rangeAlongAxis(const Point& pos, const Point& axis, 
				  double& tmin, double& tmax)
//===========================================================================
{
  tmin = std::numeric_limits<double>::max();
  tmax = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point point(xyz[0], xyz[1], xyz[2]);
      double tpar = (point - pos)*axis;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
    }
}

//===========================================================================
void RevEngRegion::identifyOutPoints(vector<pair<double,double> >& distang,
				     double tol, double angtol, double angtol2,
				     vector<vector<RevEngPoint*> >& out_groups,
				     vector<RevEngPoint*>& remaining)
//===========================================================================
{
  // Identify points with poor accuracy
  double fac = 2.0;
  vector<RevEngPoint*> points_poor;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    if ((distang[ki].second > angtol && distang[ki].first > tol) ||
	distang[ki].first > fac*tol || distang[ki].second > angtol2)
      points_poor.push_back(group_points_[ki]);
    else
      remaining.push_back(group_points_[ki]);

  // Release points in the inner of the group
  if (points_poor.size() > 0)
    {
      vector<RevEngPoint*> inner;
      connectedGroups(points_poor, out_groups, true, inner);
      if (inner.size() > 0)
	remaining.insert(remaining.end(), inner.begin(), inner.end());
    }
}

//===========================================================================
shared_ptr<ParamSurface> RevEngRegion::surfaceWithAxis(vector<RevEngPoint*>& points,
						       Point axis, Point pos,
						       Point mainaxis[3])
//===========================================================================
{
  shared_ptr<ParamSurface> surf;
  if (!hasSurface())
    return surf;

  shared_ptr<ParamSurface> init_surf = getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(init_surf);
  if (!elem.get())
    return surf;

  Point low = bbox_.low();
  Point high = bbox_.high();
  Point loc = elem->location();
  if (hasBlendEdge())
    {
      // Maintain radius while updating surface
      surf = init_surf;  // For the time being
    }
  else if (elem->instanceType() == Class_Plane)
    {
      Point init_norm = elem->direction();
      surf = RevEngUtils::planeWithAxis(points, axis, loc, mainaxis);
      if (axis*init_norm < 0.0)
	surf->swapParameterDirection();
    }
  else if (elem->instanceType() == Class_Cylinder)
    surf = RevEngUtils::cylinderWithAxis(points, axis, low, high,
					 mainaxis);
  else if (elem->instanceType() == Class_Torus)
    surf = RevEngUtils::torusWithAxis(points, axis,
				      pos.dimension() == 0 ? loc : pos,
				      mainaxis);
  else if (elem->instanceType() == Class_Sphere)
    surf = RevEngUtils::sphereWithAxis(points, axis, mainaxis);
  else if (elem->instanceType() == Class_Cone)
    surf = RevEngUtils::coneWithAxis(points, axis, low, high,
				     mainaxis);

  return surf;
}


//===========================================================================
void RevEngRegion::axisFromAdjacent(double angtol, vector<Point>& axis)
//===========================================================================
{
  vector<pair<shared_ptr<ElementarySurface>, RevEngRegion*>  > adj_elem, adj_elem_base;
  getAdjacentElemInfo(adj_elem, adj_elem_base);

  for (size_t ki=0; ki<adj_elem.size(); ++ki)
    {
      if (adj_elem[ki].first->instanceType() == Class_Plane ||
	  adj_elem[ki].first->instanceType() == Class_Cylinder ||
	  adj_elem[ki].first->instanceType() == Class_Torus)
	{
	  Point dir = adj_elem[ki].first->direction();
	  size_t kj;
	  for (kj=0; kj<axis.size(); ++kj)
	    {
	      double ang = dir.angle(axis[kj]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang < angtol)
		break;
	    }
	  if (kj == axis.size())
	    axis.push_back(dir);
	}
    }
}

//===========================================================================
bool RevEngRegion::divideWithSegInfo(int seg_ix, int min_pt_reg, 
				     vector<vector<RevEngPoint*> >& sep_groups,
				     vector<RevEngPoint*>& single_pts)
//===========================================================================
{
  if (seg_ix < 0 || seg_ix >= (int)seg_info_.size())
    return false;

  if (seg_info_[seg_ix]->type_ != 1)
    return false;  // No strategy

  // Extract points outside a belt around the specified axis
  Point axis = seg_info_[seg_ix]->axis_;
  Point pos = seg_info_[seg_ix]->loc_;
  double lim1 = seg_info_[seg_ix]->min_dist_;
  double lim2 = seg_info_[seg_ix]->max_dist_;
  vector<RevEngPoint*> inside, outside;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point onaxis = pos + ((xyz2-pos)*axis)*axis;
      double dd = xyz2.dist(onaxis);
      if (dd >= lim1 && dd <= lim2)
	inside.push_back(group_points_[kr]);
      else
	outside.push_back(group_points_[kr]);
    }

  // Check for a reasonably balanced segmentation
  int min_num = std::min(min_pt_reg, (int)group_points_.size()/10);
  if ((int)inside.size() < min_num || (int)outside.size() < min_num)
    return false;

  // Ensure connected point groups
  vector<std::vector<RevEngPoint*> > sep1, sep2;
  std::vector<RevEngPoint*> dummy;
  connectedGroups(inside, sep1, false, dummy);
  connectedGroups(outside, sep2, false, dummy);
  sep1.insert(sep1.end(), sep2.begin(), sep2.end());

  // Identify largest group
  int max_num = 0;
  int ix = -1;
  for (size_t ki=0; ki<sep1.size(); ++ki)
    if ((int)sep1[ki].size() > max_num)
      {
	max_num = (int)sep1[ki].size();
	ix = (int)ki;
      }

  if (ix >= 0)
    {
      group_points_.clear();
      group_points_ = sep1[ix];
      
      // Update bounding box and principal curvature summary
      updateInfo();
      
      for (size_t kj=0; kj<sep1.size(); ++kj)
	{
	  if ((int)kj == ix)
	    continue;
	  if (sep1[kj].size() == 1)
	    {
	      sep1[kj][0]->unsetRegion();
	      single_pts.push_back(sep1[kj][0]);
	    }
	  else
	    sep_groups.push_back(sep1[kj]);
	}
    }
  return true;
}


//===========================================================================
bool RevEngRegion::extractCylByAxis(Point mainaxis[3], int min_point,
				    int min_pt_reg, double tol, 
				    double angtol, int prefer_elementary,
				    vector<shared_ptr<HedgeSurface> >& hedgesfs,
				    vector<shared_ptr<RevEngRegion> >& added_reg,
				    vector<vector<RevEngPoint*> >& out_groups,
				    vector<RevEngPoint*>& single_pts)
//===========================================================================
{
  // Get axis from adjacent
  vector<Point> axis;
  axisFromAdjacent(angtol, axis);
  
  double pihalf = 0.5*M_PI;
  double axisang = 0.1*M_PI;
  vector<vector<RevEngPoint*> > axis_groups(3);
  vector<RevEngPoint*> remaining;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal = group_points_[kr]->getLocFuncNormal();
      int min_ix = -1;
      double min_ang = pihalf;
      for (int ka=0; ka<3; ++ka)
	{
	  double ang = mainaxis[ka].angle(normal);
	  ang = fabs(pihalf-ang);
	  if (ang < min_ang)
	    {
	      min_ang = ang;
	      min_ix = ka;
	    }
	}
      if (min_ang < axisang && min_ix >= 0)
	axis_groups[min_ix].push_back(group_points_[kr]);
      else
	remaining.push_back(group_points_[kr]);
    }
      
#ifdef DEBUG_SEGMENT
  std::ofstream of("axis_groups.g2");
  for (size_t ki=0; ki<axis_groups.size(); ++ki)
    {
      if (axis_groups[ki].size() > 0)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of <<  axis_groups[ki].size() << std::endl;
	  for (int ka=0; ka<(int)axis_groups[ki].size(); ++ka)
	    of << axis_groups[ki][ka]->getPoint() << std::endl;
	}
    }
  if (remaining.size() > 0)
    {
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of <<  remaining.size() << std::endl;
      for (int ka=0; ka<(int)remaining.size(); ++ka)
	of << remaining[ka]->getPoint() << std::endl;
    }
#endif

  size_t remain_size = remaining.size();
  for (int ka=0; ka<(int)axis_groups.size(); ++ka)
    if (axis_groups[ka].size() < remain_size)
      {
	remaining.insert(remaining.end(), axis_groups[ka].begin(), axis_groups[ka].end());
	axis_groups[ka].clear();
      }
  
  double remain_fac = 0.9;
  if ((double)remaining.size()/(double)group_points_.size() > remain_fac)
    return false;  // Not a feasible segmentation
  for (int ka=0; ka<3; ++ka)
    if (axis_groups[ka].size() == group_points_.size())
      return false;

  // Split into connected groups
  for (size_t ki=0; ki<axis_groups.size(); ++ki)
    {
      if ((int)axis_groups[ki].size() < min_point)
	{
	  if (axis_groups[ki].size() > 0)
	    remaining.insert(remaining.end(), axis_groups[ki].begin(),
			     axis_groups[ki].end());
	  continue;
	}
      
      vector<vector<RevEngPoint*> > grouped;
      std::vector<RevEngPoint*> dummy;
      connectedGroups(axis_groups[ki], grouped, false, dummy);
      
      for (size_t kj=0; kj<grouped.size(); ++kj)
	{
	  if ((int)grouped[kj].size() < min_point)
	    {
	      remaining.insert(remaining.end(), grouped[kj].begin(),
			       grouped[kj].end());
	      continue;
	    }

	  // Try to fit a cylinder
	  shared_ptr<RevEngRegion> tmp_reg(new RevEngRegion(classification_type_,
							    edge_class_type_,
							    grouped[kj]));
	  tmp_reg->setRegionAdjacency();
	  vector<HedgeSurface*> prevsfs;
	  bool repeat = false;
	  bool found = tmp_reg->extractCylinder(tol, min_point, min_pt_reg,
						angtol,
						prefer_elementary, hedgesfs,
						prevsfs, out_groups, repeat);
	  added_reg.push_back(tmp_reg);
	}
    }
	  
  vector<vector<RevEngPoint*> > grouped2;
  std::vector<RevEngPoint*> dummy;
  connectedGroups(remaining, grouped2, false, dummy);
      
  int max_num = 0;
  int max_ix = -1;
  for (size_t ki=0; ki<grouped2.size(); ++ki)
    if ((int)grouped2[ki].size() > max_num)
      {
	max_num = (int)grouped2[ki].size();
	max_ix = (int)ki;
      }

  if (max_ix >= 0)
    {
      std::swap(group_points_, grouped2[max_ix]);
      std::swap(grouped2[max_ix], grouped2[grouped2.size()-1]);
      grouped2.pop_back();
      updateInfo(tol, angtol);
    }
  if (grouped2.size() > 0)
    out_groups.insert(out_groups.end(), grouped2.begin(), grouped2.end());
 
  return true;
}

//===========================================================================
void RevEngRegion::computeFracNorm(double angtol, Point mainaxis[3],
				   int nmb_axis[3], double& in_frac1,
				   double& in_frac2)
//===========================================================================
{
  double pihalf = 0.5*M_PI;
  double axisang = 0.1*M_PI;
  nmb_axis[0] = nmb_axis[1] = nmb_axis[2] = 0;
  if (normalcone_.angle() <= angtol)
    in_frac1 = 1.0;
  else
    {
      int nmb_in = 0;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	{
	  Point normal = group_points_[kr]->getLocFuncNormal();
	  if (avnorm_.angle(normal) <= angtol)
	    nmb_in++;
	  for (int ka=0; ka<3; ++ka)
	    {
	      double ang = mainaxis[ka].angle(normal);
	      if (fabs(pihalf-ang) < axisang)
		nmb_axis[ka]++;
	    }
	      
	}
      in_frac1 = (double)nmb_in/(double)group_points_.size();
    }
  frac_norm_in_ = in_frac1;

  if (normalcone2_.angle() <= angtol)
    in_frac2 = 1.0;
  else
    {
      int nmb_in = 0;
      for (size_t kr=0; kr<group_points_.size(); ++kr)
	{
	  Point normal = group_points_[kr]->getTriangNormal();
	  if (avnorm2_.angle(normal) <= angtol)
	    nmb_in++;
	  for (int ka=0; ka<3; ++ka)
	    {
	      double ang = mainaxis[ka].angle(normal);
	      if (fabs(pihalf-ang) < axisang)
		nmb_axis[ka]++;
	    }
	      
	}
      in_frac2 = (double)nmb_in/(double)group_points_.size();
    }
  frac_norm_in2_ = in_frac2;
}

//===========================================================================
void
RevEngRegion::initPlaneCyl(int min_point, int min_pt_reg, double tol, 
			   double angtol, Point mainaxis[3],
			   double zero_H, double zero_K,
			   vector<shared_ptr<HedgeSurface> >& hedgesfs,
			   vector<vector<RevEngPoint*> >& out_groups,
			   vector<RevEngPoint*>& single_pts, bool& repeat)
//===========================================================================
{
#ifdef DEBUG_SEGMENT
  std::ofstream of1("init_sf_region.g2");
  writeRegionPoints(of1);
#endif

  repeat = false;
  
  // Count fraction of normals closer to the centre than the tolerance
  double angtol2 = 2.0*angtol;
  //double dfrac = 0.75;
  double in_frac, in_frac2;
  int nmb_axis[3];
  computeFracNorm(angtol2, mainaxis, nmb_axis, in_frac, in_frac2);

  int min_pt_reg2 = 4*min_pt_reg/5; // Allow some slack

  // Initial approximation with plane and cylinder
  vector<shared_ptr<ElementarySurface> > elemsf(6);
  vector<double> maxd(6), avd(6);
  vector<int> num_in(6), num2_in(6);
  vector<vector<double> > param(6);
  vector<vector<pair<double,double> > > d_a(6);
  vector<int> surfflag(6, NOT_SET);
  
  double in_lim = 0.9;
  shared_ptr<Plane> plane = computePlane(group_points_, avnorm_, mainaxis);
  vector<RevEngPoint*> inpt1, outpt1, inpt2, outpt2;
  int num_ang_in1=0, num_ang_in2=0;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(), plane, 
			  tol, maxd[4], avd[4], num_in[4], num2_in[4],
			  inpt1, outpt1, param[4], d_a[4], angtol);
  for (size_t kr=0; kr<d_a[4].size(); ++kr)
    {
      if (d_a[4][kr].second < angtol)
	num_ang_in1++;
    }
  elemsf[4] = plane;
  surfflag[4] = defineSfFlag(0, tol, num_in[4], num2_in[4], avd[4],
			      false);

  shared_ptr<Cylinder> cyl = computeCylinder(group_points_, tol);
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(), cyl, 
			  tol, maxd[5], avd[5], num_in[5], num2_in[5],
			  inpt2, outpt2, param[5], d_a[5], angtol);

  for (size_t kr=0; kr<d_a[5].size(); ++kr)
    {
      if (d_a[5][kr].second < angtol)
	num_ang_in2++;
    }
  elemsf[5] = cyl;
  surfflag[5] = defineSfFlag(0, tol, num_in[5], num2_in[5], avd[5],
			      true);

  // Check with reduced point set
  vector<vector<RevEngPoint*> > remaining(4);
  vector<vector<vector<RevEngPoint*> > > extracted(6);  // Two last is empty
    
  int min_ang_in = (int)group_points_.size()/20;
  int min_ang_in2 = (int)group_points_.size()/5;
  double dfac = 3.0;
  int nfac = 2;
  if (accuracyOK(min_point, tol, num_in[4], avd[4]) == false &&
      avd[4] < avd[5] && avd[4] < dfac*tol && num2_in[4] > num2_in[5] && 
      num_ang_in1 > min_ang_in && num_ang_in1 > min_point)
    {
      // Search for a (reasonable) connected planar component in the point cloud
      Point vec = plane->direction();
      elemsf[0] = planarComponent(vec, min_point, min_pt_reg2,
				  tol, angtol, mainaxis,
				  remaining[0], extracted[0]);
    }

    if (accuracyOK(min_point, tol, num_in[5], avd[5]) == false &&
      normalcone_.greaterThanPi() && avd[5] < avd[4] && avd[5] < dfac*tol &&
      num2_in[5] >= nfac*num_in[5])
    {
      // Potential helical surface. Cleanup in points and register occurance
      elemsf[1] = helicalComponent(cyl, tol, angtol, min_point, min_pt_reg2,
				   avd[5], num2_in[5], d_a[5],
				   remaining[1], extracted[1],
				   maxd[1], avd[1], num_in[1], num2_in[1],
				   param[1], d_a[1]);
    }

  vector<vector<RevEngPoint*> > out_gr2;
  if ((in_frac >= in_lim && (!(MAH_ > 10.0*MAK_ && MAH_ > 0.1))) ||
      accuracyOK(min_point, tol, num_in[4], avd[4]) ||
      num_ang_in1 >= min_ang_in2)
    {
      vector<RevEngPoint*> remain1;
      vector<vector<RevEngPoint*> > out_gr1;
      identifyOutPoints(d_a[4], tol, angtol, angtol2, out_gr1, remain1);

#ifdef DEBUG_SEGMENT
      if (out_gr1.size() > 0)
	{
	  std::ofstream ofa("ang_points1.g2");
	  ofa << "400 1 0 4 255 0 0 255" << std::endl;
	  ofa << remain1.size() << std::endl;
	  for (size_t kr=0; kr<remain1.size(); ++kr)
	    ofa << remain1[kr]->getPoint() << std::endl;
	  for (size_t kh=0; kh<out_gr1.size(); ++kh)
	    {
	      ofa << "400 1 0 4 50 50 155 255" << std::endl;
	      ofa << out_gr1[kh].size() << std::endl;
	      for (size_t kr=0; kr<out_gr1[kh].size(); ++kr)
		ofa << out_gr1[kh][kr]->getPoint() << std::endl;
	    }
	}
#endif

      if ((int)remain1.size() > min_pt_reg2)
	{
	  vector<vector<RevEngPoint*> > connected;
	  vector<RevEngPoint*> inner;
	  connectedGroups(remain1, connected, false, inner);
	  for (size_t kr=1; kr<connected.size(); ++kr)
	    if (connected[kr].size() > connected[0].size())
	      std::swap(connected[0], connected[kr]);
	  std::swap(remain1, connected[0]);
	  out_gr1.insert(out_gr1.end(), connected.begin()+1, connected.end());

	  elemsf[2] = computePlane(remain1, avnorm_, mainaxis);
	  remaining[2] = remain1;
	  extracted[2] = out_gr1;
	}
    }
  
  if (MAH_ > 3*MAK_ || accuracyOK(min_point, tol, num_in[5], avd[5]) ||
      num_ang_in2 >= min_ang_in2)
    {
      vector<RevEngPoint*> remain2;
      identifyOutPoints(d_a[5], tol, angtol, angtol2, out_gr2, remain2);

#ifdef DEBUG_SEGMENT
      if (out_gr2.size() > 0)
	{
	  std::ofstream ofa("ang_points2.g2");
	  ofa << "400 1 0 4 255 0 0 255" << std::endl;
	  ofa << remain2.size() << std::endl;
	  for (size_t kr=0; kr<remain2.size(); ++kr)
	    ofa << remain2[kr]->getPoint() << std::endl;
	  for (size_t kh=0; kh<out_gr2.size(); ++kh)
	    {
	      ofa << "400 1 0 4 50 50 155 255" << std::endl;
	      ofa << out_gr2[kh].size() << std::endl;
	      for (size_t kr=0; kr<out_gr2[kh].size(); ++kr)
		ofa << out_gr2[kh][kr]->getPoint() << std::endl;
	    }
	}
#endif

       double rfac = 0.5;
       //size_t remain2_size = remain2.size();
       if ((double)remain2.size() < rfac*(double)group_points_.size())
       {
	 vector<RevEngPoint*> ang_pts;
	 vector<RevEngPoint*> remain_ang;
	 double dtol = 2.0*avd[5];
	 identifyAngPoints(d_a[5], angtol, tol, dtol, ang_pts, remain_ang);
      
#ifdef DEBUG_SEGMENT
	 std::ofstream ofa("ang_pts_cyl.g2");
	 ofa << "400 1 0 4 155 50 50 255" << std::endl;
	 ofa << ang_pts.size() << std::endl;
	 for (size_t kr=0; kr<ang_pts.size(); ++kr)
	   ofa << ang_pts[kr]->getPoint() << std::endl;
#endif
	 if (remain_ang.size() > remain2.size())
	   {
	     std::swap(remain_ang, remain2);
	     vector<vector<RevEngPoint*> > connected2;
	     vector<RevEngPoint*> dummy_inner;
	     connectedGroups(ang_pts, connected2, false, dummy_inner);
	     std::swap(out_gr2, connected2);
	   }
       }
	 
       if ((int)remain2.size() > min_pt_reg2)
	 {
	   shared_ptr<Cylinder> cyl2 = computeCylinder(remain2, tol);
 
	   double len = bbox_.low().dist(bbox_.high());
	   Point axis = cyl2->direction();
	   Point Cx = cyl2->direction2();
	   Point loc = cyl2->location();
	   shared_ptr<Cone> cone =
	     RevEngUtils::coneWithAxis(remain2, axis, Cx, loc, len);

	   // Check accuracy with respect to reduced point cloud
	   // For cylinder
	   double maxdist3, avdist3;
	   int num_inside3, num2_inside3=0;
	   vector<RevEngPoint*> inpt3, outpt3;
	   vector<double> parvals3;
	   vector<pair<double, double> > dist_ang3;
	   RevEngUtils::distToSurf(remain2.begin(), remain2.end(), cyl2, 
				   tol, maxdist3, avdist3, num_inside3,
				   num2_inside3, inpt3, outpt3, parvals3,
				   dist_ang3, angtol);
	   int sf_flag3 = defineSfFlag((int)remain2.size(), 0, tol, num_inside3,
				       num2_inside3, avdist3, true);

	   // For cone
	   int sf_flag4 = NOT_SET;
	   double maxdist4, avdist4;
	   int num_inside4, num2_inside4=0;
	   vector<RevEngPoint*> inpt4, outpt4;
	   vector<double> parvals4;
	   vector<pair<double, double> > dist_ang4;
	   double cone_fac = 0.1*angtol;
	   if (sf_flag3 >= ACCURACY_POOR ||
	       fabs(cone->getConeAngle()) > cone_fac*angtol)
	     {
	       vector<RevEngPoint*> inpt4, outpt4;
	       RevEngUtils::distToSurf(remain2.begin(), remain2.end(), cone, 
				       tol, maxdist4, avdist4, num_inside4,
				       num2_inside4, inpt4, outpt4, parvals4,
				       dist_ang4, angtol);
	       sf_flag4 = defineSfFlag((int)remain2.size(), 0, tol, num_inside4,
				       num2_inside4, avdist4, true);
	     }

	   if (std::min(sf_flag3, sf_flag4) < ACCURACY_POOR)
	     {
	       remaining[3] = remain2;
	       extracted[3] = out_gr2;

	       if (sf_flag3 < sf_flag4 ||
		   (sf_flag3 == sf_flag4 &&
		    (avdist3 < avdist4 || num_inside3 > num_inside4)))
		 {
		   elemsf[3] = cyl2;
		   maxd[3] = maxdist3;
		   avd[3] = avdist3;
		   num_in[3] = num_inside3;
		   num2_in[3] = num2_inside3;
		   param[3] = parvals3;
		   d_a[3] = dist_ang3;
		   surfflag[3] = sf_flag3;
		 }
	       else
		 {
		   elemsf[3] = cone;
		   maxd[3] = maxdist4;
		   avd[3] = avdist4;
		   num_in[3] = num_inside4;
		   num2_in[3] = num2_inside4;
		   param[3] = parvals4;
		   d_a[3] = dist_ang4;
		   surfflag[3] = sf_flag4;
		 }
	     }
	 }
    }

  // Compute accuracy for remaining candidates
  for (size_t ki=0; ki<elemsf.size(); ++ki)
    {
      if (!elemsf[ki].get())
	continue;  // Not a candidate
      if (surfflag[ki] < NOT_SET)
	continue;  // Already computed

      if (param[ki].size() == 0)
	{
	   vector<RevEngPoint*> inpt3, outpt3;
	   RevEngUtils::distToSurf(remaining[ki].begin(), remaining[ki].end(),
				   elemsf[ki], tol, maxd[ki], avd[ki],
				   num_in[ki], num2_in[ki], inpt3, outpt3,
				   param[ki], d_a[ki], angtol);
	}
      bool cyllike = (elemsf[ki]->instanceType() == Class_Cylinder ||
		      elemsf[ki]->instanceType() == Class_Cone);
      surfflag[ki] = defineSfFlag((int)remaining[ki].size(), 0, tol,
				  num_in[ki], num2_in[ki], avd[ki], cyllike);
    }

  // Select surface
  int perm[6] = {4, 5, 2, 3, 0, 1};
  size_t num_pts_in[6];
  for (size_t kr=0; kr<4; ++kr)
    num_pts_in[kr] = remaining[kr].size();
  num_pts_in[4] = num_pts_in[5] = group_points_.size();
  double size_fac = 0.9;
  double avd_fac = 0.95;
  
  int ix = -1;
  for (int ka=0; ka<6; ++ka)
    {
      if (surfflag[perm[ka]] >= ACCURACY_POOR)
	continue;
      if (ix < 0)
	ix = ka;
      else
	{
	  if (surfflag[perm[ka]] == surfflag[perm[ix]] ||
	      (surfflag[perm[ka]] == ANGULAR_DEVIATION &&
	       surfflag[perm[ix]] == PROBABLE_HELIX) ||
	      (surfflag[perm[ix]] == ANGULAR_DEVIATION &&
	       surfflag[perm[ka]] == PROBABLE_HELIX))
	    {
	      //int curr_num = std::max(num_pts_in[perm[ix]], num_pts_in[perm[ka]]);
	      int curr_num = ((int)num_pts_in[perm[ix]] +
			      (int)num_pts_in[perm[ka]])/2;
	      double frac1 = (double)(num_in[perm[ix]] + num2_in[perm[ix]])/
		(double)curr_num;
	      double frac2 = (double)(num_in[perm[ka]] + num2_in[perm[ka]])/
		(double)curr_num;

	      bool prefer_first = false;
	      if (perm[ix] % 2 == 0 && perm[ka] % 2 == 1)
		{
		  // Replace plane by rotational surface?
		  double r = elemsf[perm[ka]]->radius(0.0, 0.0);
		  double d = bbox_.low().dist(bbox_.high());
		  double x2 = 4*r*r - d*d;
		  double x = (x2 > 0.0) ? r - 0.5*sqrt(x2) : r;
		  double rfac = 0.05;
		  if (x < rfac*r)
		    prefer_first = true;
		}
	      

	      if (prefer_first)
		;   // Keep current
	      else if (avd[perm[ka]] <= tol && avd_fac*avd[perm[ix]] > tol)
		ix = ka;
	      else if (size_fac*(double)num_pts_in[perm[ka]] >
		  (double)num_pts_in[perm[ix]])
		ix = ka;
	      else if (avd[perm[ka]] < avd_fac*avd[perm[ix]] && frac2 > frac1)
		ix = ka;
		  
	    }
	  else if (surfflag[perm[ka]] < surfflag[perm[ix]])
	    ix = ka;
	}
    }

  if (ix >= 0)
    {
      // Surface found
      int ix2 = perm[ix];

      // Pre check for connected group
      vector<vector<RevEngPoint*> > grouped;
      std::vector<RevEngPoint*> dummy;
      connectedGroups((ix2<4) ? remaining[ix2] : group_points_, grouped,
		      false, dummy);
      int ix3 = 0;
      for (size_t kr=1; kr<grouped.size(); ++kr)
	if (grouped[kr].size() > grouped[ix3].size())
	  ix3 = (int)kr;

      if (grouped[ix3].size() < min_pt_reg2 || (!elemsf[ix2].get()))
	ix = -1;
    }
     
  if (ix >= 0)
    {
      // Surface found
      int ix2 = perm[ix];

      if (ix2 < 4)
	std::swap(group_points_, remaining[ix2]);

      for (size_t kh=0; kh<group_points_.size(); ++kh)
	{
	  group_points_[kh]->setPar(Vector2D(param[ix2][2*kh],
					     param[ix2][2*kh+1]));
	  group_points_[kh]->setSurfaceDist(d_a[ix2][kh].first,
					    d_a[ix2][kh].second);
	}
      setAccuracy(maxd[ix2], avd[ix2], num_in[ix2], num2_in[ix2]);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(elemsf[ix2], this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(surfflag[ix2]);
      
      // Ensure connected region
      size_t num_extract = extracted[ix2].size();
      for (size_t kr=0; kr<extracted[ix2].size(); ++kr)
	for (size_t kh=0; kh<extracted[ix2][kr].size(); ++kh)
	  extracted[ix2][kr][kh]->unsetRegion();
	
      vector<vector<RevEngPoint*> > separate_groups;
      splitRegion(separate_groups);
      if (separate_groups.size() > 0)
	{
	  extracted[ix2].insert(extracted[ix2].end(), separate_groups.begin(),
				separate_groups.end());
      
	  for (size_t kr=num_extract; kr<extracted[ix2].size(); ++kr)
	    for (size_t kh=0; kh<extracted[ix2][kr].size(); ++kh)
	      extracted[ix2][kr][kh]->unsetRegion();
      
	  if (hasSurface())
	    checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
	}
	
	for (size_t kr=0; kr<extracted[ix2].size(); )
	  {
	    if (extracted[ix2][kr].size() == 1)
	      {
		extracted[ix2][kr][0]->unsetRegion();
		single_pts.push_back(extracted[ix2][kr][0]);
		extracted[ix2].erase(extracted[ix2].begin()+kr);
	      }
	    else
	      ++kr;
	  }
	out_groups = extracted[ix2];
  
#ifdef DEBUG_SEGMENT
	std::ofstream ofb("res_group.g2");
	writeRegionPoints(ofb);
	writeSurface(ofb);
#endif
	int stop_break = 1;
    }
}




//===========================================================================
void RevEngRegion::splitRegion(vector<vector<RevEngPoint*> >& separate_groups)
//===========================================================================
{
  vector<vector<RevEngPoint*> > connected;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (group_points_[ki]->visited())
	continue;
      vector<RevEngPoint*> curr_group;
      group_points_[ki]->fetchConnected(this, (int)group_points_.size(), curr_group);
      connected.push_back(curr_group);
    }

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      group_points_[ki]->unsetVisited();
    }

#ifdef DEBUG_SEGMENT
  std::ofstream of("curr_region_split2.g2");
  for (size_t kj=0; kj<connected.size(); ++kj)
    {
      of << "400 1 0 0" << std::endl;
      of << connected[kj].size() << std::endl;
      for (size_t kr=0; kr<connected[kj].size(); ++kr)
	of << connected[kj][kr]->getPoint() << std::endl;
    }
#endif
  
  if (connected.size() <= 1)
    return;
  
  int max_num = 0;
  int ix = -1;
  for (size_t ki=0; ki<connected.size(); ++ki)
    if ((int)connected[ki].size() > max_num)
      {
	max_num = (int)connected[ki].size();
	ix = (int)ki;
      }

  if (ix >= 0)
    {
      group_points_.clear();
      group_points_ = connected[ix];
      
      // Update bounding box and principal curvature summary
      updateInfo();
      
      for (size_t kj=0; kj<connected.size(); ++kj)
	{
	  if ((int)kj == ix)
	    continue;
	  for (size_t kr=0; kr<connected[kj].size(); ++kr)
	    connected[kj][kr]->unsetRegion();
	  separate_groups.push_back(connected[kj]);
	}
    }
}

//===========================================================================
shared_ptr<Plane>
RevEngRegion::planarComponent(Point vec, int min_point, int min_pt_reg,
			      double tol, double angtol, Point mainaxis[3],
			      vector<RevEngPoint*>& remaining,
			      vector<vector<RevEngPoint*> >& extracted)
//===========================================================================
{
  shared_ptr<Plane> dummy_plane;
  vector<vector<RevEngPoint*> > vec_pts(2);
  int min_size = std::max(10, std::min(min_point, (int)group_points_.size()/20));
  vector<RevEngPoint*> pts_out;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point normal = group_points_[ki]->getLocFuncNormal();
      Point normal2 = group_points_[ki]->getTriangNormal();
      double ang = vec.angle(normal);
      double ang2 = vec.angle(normal2);
      if (ang < angtol || ang2 < angtol)
	vec_pts[0].push_back(group_points_[ki]);
      else if (M_PI-ang < angtol || M_PI-ang < angtol)
	vec_pts[1].push_back(group_points_[ki]);
      else
	pts_out.push_back(group_points_[ki]);
    }

  double eps = 1.0e-6;
  double *seed = 0;
  vector<vector<RevEngPoint*> > conn_groups;
  for (int ka=0; ka<2; ++ka)
    {
      if ((int)vec_pts[ka].size() < min_size)
	continue;
      
      // Check if the group can be represented by a plane
      shared_ptr<Plane> plane = computePlane(vec_pts[ka], vec, mainaxis);
      
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      double maxdist, avdist;
      int num_inside, num2_inside;
      RevEngUtils::distToSurf(vec_pts[ka].begin(), vec_pts[ka].end(), plane, 
			      tol, maxdist, avdist, num_inside, num2_inside,
			      inpt, outpt, parvals, dist_ang, angtol);
      int sf_flag0 = defineSfFlag((int)vec_pts[ka].size(), 0, tol, num2_inside,
				 num_inside, avdist, false);
      if (sf_flag0 < ACCURACY_POOR)
	{
	  // Move points from out group if feasible
	  for (int kb=(int)pts_out.size()-1; kb>=0; --kb)
	    {
	      Vector3D xyz = pts_out[kb]->getPoint();
	      Point pos(xyz[0], xyz[1], xyz[2]);
	      double upar, vpar, dist;
	      Point close;
	      plane->closestPoint(pos, upar, vpar, close, dist, eps, 0, seed);
	      if (dist <= tol)
		{
		  vec_pts[ka].push_back(pts_out[kb]);
		  pts_out.erase(pts_out.begin()+kb);
		}
	    }
	      
	  std::vector<RevEngPoint*> dummy;
	  connectedGroups(vec_pts[ka], conn_groups, false, dummy);
	}
    }
  
  // Identify largest group
  if (conn_groups.size() == 0)
    return dummy_plane;  // No large planar component is found
  
  int max_group = 0;
  int ix = -1;
  for (size_t ki=0; ki<conn_groups.size(); ++ki)
    {
      if ((int)conn_groups[ki].size() > max_group)
	{
	  max_group = (int)conn_groups[ki].size();
	  ix = (int)ki;
	}
    }

  if (max_group < min_size)
    return dummy_plane;   // No large planar component is found

   // Update plane
  shared_ptr<Plane> plane2 = computePlane(conn_groups[ix], vec, mainaxis);
	  
#ifdef DEBUG_PLANAR
  std::ofstream of("planar_component.g2");
  writeRegionPoints(of);
  shared_ptr<Plane> plane3(plane2->clone());
  double len = bbox_.low().dist(bbox_.high());
  plane3->setParameterBounds(-len, -len, len, len);
  plane3->writeStandardHeader(of);
  plane3->write(of);
#endif

  // Check angular behaviour of points
  vector<RevEngPoint*> inpt2, outpt2;
  vector<pair<double, double> > dist_ang2;
  vector<double> parvals2;
  double maxdist2, avdist2;
  int num_inside2, num2_inside2;
  RevEngUtils::distToSurf(conn_groups[ix].begin(), conn_groups[ix].end(),
			  plane2, tol, maxdist2, avdist2, num_inside2, num2_inside2,
			  inpt2, outpt2, parvals2, dist_ang2, angtol);
  int sf_flag = defineSfFlag((int)conn_groups[ix].size(), 0, tol, num_inside2,
			     num2_inside2, avdist2, false);
  if (sf_flag >= ACCURACY_POOR)
    return dummy_plane;
  
  vector<RevEngPoint*> ang_out;
  if (sf_flag > ACCURACY_OK)
    {
      // Check if any points should be removed
      size_t num_in = conn_groups[ix].size();
      double dtol = std::min(0.5*tol, 1.5*avdist2);
      for (size_t kr=0; kr<conn_groups[ix].size(); ++kr)
	if (dist_ang2[kr].second > 2.0*angtol && dist_ang2[kr].first > dtol)
	  {
	    std::swap(conn_groups[ix][kr], conn_groups[ix][num_in-1]);
	    --num_in;
	  }
      
      if (num_in < conn_groups[ix].size())
	{
	  ang_out.insert(ang_out.end(), conn_groups[ix].begin()+num_in,
			 conn_groups[ix].end());
	  conn_groups[ix].erase(conn_groups[ix].begin()+num_in,
				conn_groups[ix].end());
#ifdef DEBUG_PLANAR
	  std::ofstream ofa("ang_points_pc.g2");
	  ofa << "400 1 0 4 50 50 155 255" << std::endl;
	  ofa << ang_out.size() << std::endl;
	  for (size_t kr=0; kr<ang_out.size(); ++kr)
	    ofa << ang_out[kr]->getPoint() << std::endl;
#endif
	}
    }

  if ((int)conn_groups[ix].size() < max_group)
    return dummy_plane;

  // A result is reached. 
  vector<vector<RevEngPoint*> > connected;
  vector<RevEngPoint*> dummy_inner;
  connectedGroups(conn_groups[ix], connected, false, dummy_inner);
  if (connected.size() > 1)
    {
      for (size_t kr=0; kr<connected.size(); ++kr)
	for (size_t kh=kr+1; kh<connected.size(); ++kh)
	  if (connected[kh].size() > connected[kr].size())
	    std::swap(connected[kh], connected[kr]);
      conn_groups[ix] = connected[0];
      conn_groups.insert(conn_groups.end(), connected.begin()+1,
			 connected.end());
    }
  
  remaining = conn_groups[ix];
  for (size_t kr=0; kr<conn_groups.size(); ++kr)
    {
      if ((int)kr == ix)
	continue;
      extracted.push_back(conn_groups[kr]);
    }
  if (ang_out.size() > 0)
    extracted.push_back(ang_out);

  vector<vector<RevEngPoint*> > connected2;
  connectedGroups(pts_out, connected2, false, dummy_inner);
  extracted.insert(extracted.end(), connected2.begin(), connected2.end());
  
  return plane2;
 }


//===========================================================================
bool RevEngRegion::planarComponent(Point vec, int min_point, int min_pt_reg,
				   double tol, double angtol, Point mainaxis[3],
				   vector<shared_ptr<HedgeSurface> >& hedgesfs,
				   vector<vector<RevEngPoint*> >& out_groups,
				   vector<RevEngPoint*>& single_pts,
				   bool create_surface)
//===========================================================================
{
  vector<vector<RevEngPoint*> > vec_pts(2);
  vector<RevEngPoint*> remaining;
  int min_size = std::max(10, std::min(min_point, (int)group_points_.size()/20));
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Point normal = group_points_[ki]->getLocFuncNormal();
      Point normal2 = group_points_[ki]->getTriangNormal();
      double ang = vec.angle(normal);
      double ang2 = vec.angle(normal2);
      if (ang < angtol || ang2 < angtol)
	vec_pts[0].push_back(group_points_[ki]);
      else if (M_PI-ang < angtol || M_PI-ang < angtol)
	vec_pts[1].push_back(group_points_[ki]);
      else
	remaining.push_back(group_points_[ki]);
    }

  double eps = 1.0e-6;
  double *seed = 0;
  vector<vector<RevEngPoint*> > conn_groups;
  for (int ka=0; ka<2; ++ka)
    {
      if ((int)vec_pts[ka].size() < min_size)
	continue;

	  // Check if the group can be represented by a plane
	  shared_ptr<Plane> plane = computePlane(vec_pts[ka], vec, mainaxis);
	  
	  vector<RevEngPoint*> inpt, outpt;
	  vector<pair<double, double> > dist_ang;
	  vector<double> parvals;
	  double maxdist, avdist;
	  int num_inside, num2_inside;
	  RevEngUtils::distToSurf(vec_pts[ka].begin(), vec_pts[ka].end(),
				  plane, tol, maxdist, avdist, num_inside, num2_inside,
				  inpt, outpt, parvals, dist_ang, -1);
	  if (num_inside > min_point && num_inside > (int)vec_pts[ka].size()/2 &&
	      avdist <= tol)
	    {
	      // Move points from remaining group if feasible
	      for (size_t kr=0; kr<remaining.size(); ++kr)
		{
		  Vector3D xyz = remaining[kr]->getPoint();
		  Point pos(xyz[0], xyz[1], xyz[2]);
		  double upar, vpar, dist;
		  Point close;
		  plane->closestPoint(pos, upar, vpar, close, dist, eps, 0, seed);
		  if (dist < tol)
		    vec_pts[ka].push_back(remaining[kr]);
		}
	      
	      std::vector<RevEngPoint*> dummy;
	      connectedGroups(vec_pts[ka], conn_groups, false, dummy);

	      int stop_break = 1;
	    }
    }
  
  // Identify largest group
  if (conn_groups.size() == 0)
    return false;  // No large planar component is found

  int max_group = 0;
  int ix = -1;
  for (size_t ki=0; ki<conn_groups.size(); ++ki)
    {
      if ((int)conn_groups[ki].size() > max_group)
	{
	  max_group = (int)conn_groups[ki].size();
	  ix = (int)ki;
	}
    }

  if (max_group < min_size)
    return false;   // No large planar component is found

  // Extract identified points
  if ((int)group_points_.size() > max_group)
    {
      for (size_t ki=0; ki<conn_groups[ix].size(); ++ki)
	removePoint(conn_groups[ix][ki]);
      std::swap(group_points_, conn_groups[ix]);
    }
  else
    conn_groups[ix].clear();

  // Update plane
  shared_ptr<Plane> plane2 = computePlane(group_points_, vec, mainaxis);
	  
#ifdef DEBUG_PLANAR
  std::ofstream of("planar_component.g2");
  writeRegionPoints(of);
  shared_ptr<Plane> plane3(plane2->clone());
  double len = bbox_.low().dist(bbox_.high());
  plane3->setParameterBounds(-len, -len, len, len);
  plane3->write(of);
#endif
  
  vector<RevEngPoint*> inpt2, outpt2;
  vector<pair<double, double> > dist_ang2;
  vector<double> parvals2;
  double maxdist2, avdist2;
  int num_inside2, num2_inside2;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  plane2, tol, maxdist2, avdist2, num_inside2, num2_inside2,
			  inpt2, outpt2, parvals2, dist_ang2, angtol);
  int sf_flag = defineSfFlag(0, tol, num_inside2, num2_inside2,
			     avdist2, false);
  if (sf_flag > ACCURACY_OK && sf_flag < NOT_SET)
    {
      // Check if any points should be removed
      vector<RevEngPoint*> ang_points;
      double dtol = std::min(0.5*tol, 1.5*avdist2);
      identifyAngPoints(dist_ang2, 2.0*angtol, dtol, ang_points);
#ifdef DEBUG_PLANAR
      if (ang_points.size() > 0)
	{
	  std::ofstream ofa("ang_points_pc.g2");
	  ofa << "400 1 0 4 50 50 155 255" << std::endl;
	  ofa << ang_points.size() << std::endl;
	  for (size_t kr=0; kr<ang_points.size(); ++kr)
	    ofa << ang_points[kr]->getPoint() << std::endl;
	}
#endif

      double del_frac = 0.01;
      if ((double)ang_points.size() > del_frac*(double)group_points_.size())
	{
	  vector<vector<RevEngPoint*> > separate_groups;
	  extractSpesPoints(ang_points, separate_groups, true);
	  if (separate_groups.size() > 0)
	    {
	      updateInfo(tol, angtol);
	      inpt2.clear();
	      outpt2.clear();
	      parvals2.clear();
	      dist_ang2.clear();
	      plane2 = computePlane(group_points_, avnorm_, mainaxis);
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      plane2, tol, maxdist2, avdist2,
				      num_inside2, num2_inside2,
				      inpt2, outpt2, parvals2, dist_ang2, angtol);
	    }
	  for (size_t kr=0; kr<separate_groups.size(); ++kr)
	    {
	      if (separate_groups[kr].size() == 1)
		{
		  separate_groups[kr][0]->unsetRegion();
		  single_pts.push_back(separate_groups[kr][0]);
		}
	      else
		out_groups.push_back(separate_groups[kr]);
	    }
	  sf_flag = defineSfFlag(0, tol, num_inside2, num2_inside2,
				 avdist2, false);
	}
    }

  if (create_surface)
    {
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals2[2*ki],parvals2[2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang2[ki].first, dist_ang2[ki].second);
	}
      setAccuracy(maxdist2, avdist2, num_inside2, num2_inside2);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(plane2, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }
  updateInfo(tol, angtol);
    
   
  // Split remaining points into connected components
  if (conn_groups[ix].size() > 0)
    {
      shared_ptr<RevEngRegion> reg(new RevEngRegion(classification_type_,
						    edge_class_type_,
						    conn_groups[ix]));
      vector<vector<RevEngPoint*> > connected;
      reg->splitRegion(connected);
      out_groups.push_back(reg->getPoints());
      for (size_t ki=0; ki<connected.size(); ++ki)
	{
	  if (connected[ki].size() == 1)
	    {
	      connected[ki][0]->unsetRegion();
	      single_pts.push_back(connected[ki][0]);
	    }
	  else
	    out_groups.push_back(connected[ki]);
	}
    }


  return true;
}

//===========================================================================
shared_ptr<ElementarySurface>
RevEngRegion::helicalComponent(shared_ptr<Cylinder> cyl,  double tol,
			       double angtol, int min_point,
			       int min_pt_reg, double avdist, int num_in2,
			       vector<pair<double,double> >& dist_ang,
			       vector<RevEngPoint*>& remaining,
			       vector<vector<RevEngPoint*> >& extracted,
			       double& maxd_out, double& avd_out,
			       int& num_in_out, int& num2_in_out,
			       vector<double>& parvals_out,
			       vector<pair<double,double> >& distang_out)
//===========================================================================
{
  shared_ptr<ElementarySurface> dummy_surf;
  // Compute rotated info
  int num_in_lin1, num_in_cub1;
  double avdist_lin1, avdist_cub1;
  shared_ptr<Cone> cone;
  analyseCylRotate(cyl, tol, avdist, num_in2, avdist_lin1, num_in_lin1,
		   avdist_cub1, num_in_cub1, cone);

  // Restrict to linear cases
  double lfac = 1.2;
  if (avdist_lin1 > lfac*avdist_cub1)
    return dummy_surf;
  
  // Identify distant points
  double eps = 1.0e-6;
  double dlim = 2.0*avdist;
  vector<RevEngPoint*> dist_pts;
  vector<RevEngPoint*> in_pts;
  BoundingBox bb(3);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      bool within = true;
      Vector3D xyz = group_points_[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      if (cone.get())
	{
	  double upar, vpar, dist;
	  Point close;
	  cone->closestPoint(pos, upar, vpar, close, dist, eps);
	  if (dist > dlim)
	    within = false;	  
	}
      else
	{
	  if (dist_ang[ki].first > dlim)
	    within = false;
	}

      if (within)
	{
	  in_pts.push_back(group_points_[ki]);
	  bb.addUnionWith(pos);
	}
      else
	dist_pts.push_back(group_points_[ki]);
    }

  if ((int)in_pts.size() < min_point)
    return dummy_surf;

#ifdef DEBUG_CYL
  std::ofstream of("helical_split.g2");
  of << "400 1 0 4 0 0 255 255" << std::endl;
  of << in_pts.size() << std::endl;
  for (size_t ki=0; ki<in_pts.size(); ++ki)
    of << in_pts[ki]->getPoint() << std::endl;
  
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << dist_pts.size() << std::endl;
  for (size_t ki=0; ki<dist_pts.size(); ++ki)
    of << dist_pts[ki]->getPoint() << std::endl;
#endif
  
  // Check accuracy with respect to reduced point cloud
  // For cylinder
  double maxdist2, avdist2;
  int num_inside2, num2_inside2=0;
  vector<RevEngPoint*> inpt2, outpt2;
  vector<double> parvals2;
  vector<pair<double, double> > dist_ang2;
  shared_ptr<Cylinder> cyl2 = computeCylinder(in_pts, tol);
  RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), cyl2, 
			  tol, maxdist2, avdist2, num_inside2, num2_inside2,
			  inpt2, outpt2, parvals2, dist_ang2, angtol);
  int sf_flag2 = defineSfFlag((int)in_pts.size(), 0, tol, num_inside2,
			      num2_inside2, avdist2, true);

  // For cone
  double len = bb.low().dist(bb.high());
  Point axis = cyl2->direction();
  Point Cx = cyl2->direction2();
  Point loc = cyl2->location();
  shared_ptr<Cone> cone2 =
    RevEngUtils::coneWithAxis(in_pts, axis, Cx, loc, len);
  
  double maxdist3, avdist3;
  int num_inside3, num2_inside3=0;
  vector<RevEngPoint*> inpt3, outpt3;
  vector<double> parvals3;
  vector<pair<double, double> > dist_ang3;
  RevEngUtils::distToSurf(in_pts.begin(), in_pts.end(), cone2,
			  tol, maxdist3, avdist3, num_inside3, num2_inside3,
			  inpt3, outpt3, parvals3, dist_ang3, angtol);
  int sf_flag3 = defineSfFlag((int)in_pts.size(), 0, tol, num_inside3,
			      num2_inside3, avdist3, true);

  if (std::min(sf_flag2, sf_flag3) >= ACCURACY_POOR)
    return dummy_surf;

  remaining = in_pts;
  vector<vector<RevEngPoint*> > connected;
  vector<RevEngPoint*> dummy_inner;
  connectedGroups(dist_pts, connected, false, dummy_inner);
  extracted = connected;

  if (sf_flag2 < sf_flag3 ||
      (sf_flag2 == sf_flag3 &&
       (avdist2 < avdist3 || num_inside2 > num_inside3)))
    {
      maxd_out = maxdist2;
      avd_out = avdist2;
      num_in_out = num_inside2;
      num2_in_out = num2_inside2;
      parvals_out = parvals2;
      distang_out = dist_ang2;
      return cyl2;
    }
  else
    {
      maxd_out = maxdist3;
      avd_out = avdist3;
      num_in_out = num_inside3;
      num2_in_out = num2_inside3;
      parvals_out = parvals3;
      distang_out = dist_ang3;
      return cone2;
    }
}


//===========================================================================
bool RevEngRegion::defineHelicalInfo(shared_ptr<Cylinder> cyl,  double tol,
				     double angtol, int min_point,
				     int min_pt_reg, double avdist,
				     int num_in1, int num_in2,
				     vector<pair<double,double> >& dist_ang,
				     vector<shared_ptr<HedgeSurface> >& hedgesfs,
				     vector<vector<RevEngPoint*> >& out_groups,
				     vector<RevEngPoint*>& single_pts)
//===========================================================================
{
  // Compute rotated info
  int num_in_lin1, num_in_cub1;
  double avdist_lin1, avdist_cub1;
  shared_ptr<Cone> cone;
  analyseCylRotate(cyl, tol, avdist, num_in2, avdist_lin1, num_in_lin1,
		   avdist_cub1, num_in_cub1, cone);

  // Restrict to linear cases
  double lfac = 1.2;
  if (avdist_lin1 > lfac*avdist_cub1)
    return false;

  // Identify distant points
  double dlim = 2.0*avdist;
  vector<RevEngPoint*> dist_pts;
  vector<RevEngPoint*> remaining;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      if (dist_ang[ki].first > dlim)
	dist_pts.push_back(group_points_[ki]);
      else
	remaining.push_back(group_points_[ki]);
    }

  // Check accuracy with respect to reduced point cloud
  double maxdist2, avdist2;
  int num_inside2, num2_inside2=0;
  vector<RevEngPoint*> inpt2, outpt2;
  vector<double> parvals2;
  vector<pair<double, double> > dist_ang2;
  shared_ptr<Cylinder> cyl2 = computeCylinder(remaining, tol);
  RevEngUtils::distToSurf(remaining.begin(), remaining.end(),
			  cyl2, tol, maxdist2, avdist2, num_inside2, num2_inside2,
			  inpt2, outpt2, parvals2, dist_ang2, angtol);

  // Check level of change
  double ang = cyl->direction().angle(cyl2->direction());
  ang = std::min(ang, M_PI-ang);
  double distfac = 0.5;
  double angmin = 0.1*M_PI;
  if (ang > angmin && avdist2 > tol &&
      (double)num2_inside2/(double)num_in2 > distfac)
    return false;  // Not likely that the computed cylinder is a good representation
  
  int nfac = 2;
  if (num2_inside2 < nfac*num_inside2 && (num_inside2 < (int)remaining.size()/2 ||
					   avdist2 > tol))
    return false;   // Probably not a helical surface and not sufficient accuracy for
  // a cylinder

#ifdef DEBUG_CYL
  std::ofstream of("helical_split.g2");
  of << "400 1 0 4 0 0 255 255" << std::endl;
  of << remaining.size() << std::endl;
  for (size_t ki=0; ki<remaining.size(); ++ki)
    of << remaining[ki]->getPoint() << std::endl;
  
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << dist_pts.size() << std::endl;
  for (size_t ki=0; ki<dist_pts.size(); ++ki)
    of << dist_pts[ki]->getPoint() << std::endl;
#endif
  // Set info to remaining
  for (size_t ki=0; ki<remaining.size(); ++ki)
    {
      remaining[ki]->setPar(Vector2D(parvals2[2*ki],parvals2[2*ki+1]));
      remaining[ki]->setSurfaceDist(dist_ang2[ki].first, dist_ang2[ki].second);
    }
  std::swap(group_points_, remaining);
  updateInfo(tol, angtol);

  // Make sure that the remaining points are connected
  vector<vector<RevEngPoint*> > separate_groups;
  splitRegion(separate_groups);
  if (separate_groups.size() > 0)
    {
      // Recompute
      cyl2 = computeCylinder(group_points_, tol);
      parvals2.clear();
      inpt2.clear();
      outpt2.clear();
      
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      cyl2, tol, maxdist2, avdist2,
			      num_inside2, num2_inside2,
			      inpt2, outpt2, parvals2, dist_ang2, angtol);
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals2[2*ki],parvals2[2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang2[ki].first, dist_ang2[ki].second);
	}
      updateInfo(tol, angtol);
    }
  
  for (size_t ki=0; ki<separate_groups.size(); ++ki)
    {
      if (separate_groups[ki].size() == 1)
	{
	  separate_groups[ki][0]->unsetRegion();
	  single_pts.push_back(separate_groups[ki][0]);
	}
      else
	out_groups.push_back(separate_groups[ki]);
    }

  // Surface flag
  int sf_flag = defineSfFlag(0, tol, num_inside2, num2_inside2,
			     avdist2, true);

  if (sf_flag < ACCURACY_POOR)
    {
      setAccuracy(maxdist2, avdist2, num_inside2, num2_inside2);
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(cyl2, this));
      setHedge(hedge.get());
      hedgesfs.push_back(hedge);
      setSurfaceFlag(sf_flag);
    }
  else if (num2_inside2 >= nfac*num_inside2)
    {
      shared_ptr<SplineCurve> dummy;
      sweep_ = shared_ptr<SweepData>(new SweepData(3, dummy, cyl2->location(),
						   cyl2->direction(), maxdist2,
						   avdist2, num_inside2,
						   cyl2->getRadius()));
    }

  setBaseSf(cyl2, maxdist2, avdist2, num_inside2, num2_inside2);
  
  // Split remaining points into connected components
  vector<vector<RevEngPoint*> > connected;
  std::vector<RevEngPoint*> dummy;
  connectedGroups(dist_pts, connected, false, dummy);
  
  for (size_t ki=0; ki<connected.size(); ++ki)
    {
      if (connected[ki].size() == 1)
	{
	  connected[ki][0]->unsetRegion();
	  single_pts.push_back(connected[ki][0]);
	}
      else
	out_groups.push_back(connected[ki]);
    }

  return true;
}

//===========================================================================
void RevEngRegion::setAssociatedSurface(shared_ptr<ParamSurface>& surf,
					double tol, double angtol,
					int min_pt_reg,
					shared_ptr<HedgeSurface>& hedge)
//===========================================================================
{
  // Compute accuracy
  double maxdist, avdist;
  int num_in, num2_in;
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  surf, tol, maxdist, avdist, num_in, num2_in,
			  inpt, outpt, parvals, dist_ang, angtol);

  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  int sf_flag = defineSfFlag(0, tol, num_in, num2_in, avdist, cyllike);
  
  for (size_t kh=0; kh<group_points_.size(); ++kh)
    {
      group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
      group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
    }
  setAccuracy(maxdist, avdist, num_in, num2_in);
  hedge = shared_ptr<HedgeSurface>(new HedgeSurface(surf, this));
  setHedge(hedge.get());
  setSurfaceFlag(sf_flag);
}

//===========================================================================
void RevEngRegion::extractOutPoints(int dir, double tol, double angtol,
				    double angtol2,
				    vector<vector<RevEngPoint*> >& out_groups)
//===========================================================================
{
  // Compute real domain
  double dom[4];
  int ka;
  int num_pt = (int)group_points_.size();
  double dist, ang;
  for (ka=0; ka<num_pt; ++ka)
    {
      Vector2D uv = group_points_[ka]->getPar();
      group_points_[ka]->getSurfaceDist(dist, ang);
      if (ang <= angtol2)
	{
	  dom[0] = dom[1] = uv[0];
	  dom[2] = dom[3] = uv[1];
	  break;
	}
    }
  for (; ka<num_pt; ++ka)
    {
      Vector2D uv = group_points_[ka]->getPar();
      group_points_[ka]->getSurfaceDist(dist, ang);
      if (ang > angtol2)
	continue;
      dom[0] = std::min(dom[0], uv[0]);
      dom[1] = std::max(dom[1], uv[0]);
      dom[2] = std::min(dom[2], uv[1]);
      dom[3] = std::max(dom[3], uv[1]);
    }

  vector<RevEngPoint*> out1, out2, remain;
  for (ka=0; ka<num_pt; ++ka)
    {
      Vector2D uv = group_points_[ka]->getPar();
      if (uv[dir] < dom[2*dir])
	out1.push_back(group_points_[ka]);
      else if (uv[dir] > dom[2*dir+1])
	out2.push_back(group_points_[ka]);
      else
	remain.push_back(group_points_[ka]);
    }

  if (remain.size() < group_points_.size())
    {
      std::swap(remain, group_points_);
      updateInfo(tol, angtol);

      if (out1.size() > 0)
	{
	  vector<vector<RevEngPoint*> > con;
	  vector<RevEngPoint*> dummy;
	  connectedGroups(out1, con, false, dummy);
	  for (size_t ki=0; ki<con.size(); ++ki)
	    {
	      for (size_t kj=0; kj<con[ki].size(); ++kj)
		con[ki][kj]->unsetRegion();
	    }
	  out_groups.insert(out_groups.end(), con.begin(), con.end());
	}
      
      if (out2.size() > 0)
	{
	  vector<vector<RevEngPoint*> > con;
	  vector<RevEngPoint*> dummy;
	  connectedGroups(out2, con, false, dummy);
	  for (size_t ki=0; ki<con.size(); ++ki)
	    {
	      for (size_t kj=0; kj<con[ki].size(); ++kj)
		con[ki][kj]->unsetRegion();
	    }
	  out_groups.insert(out_groups.end(), con.begin(), con.end());
	}
    }
  }

//===========================================================================
int RevEngRegion::defineSfFlag(int min_point, double tol, int num_in, 
			       int num_in2, double avd, bool type_cyl)
//===========================================================================
{
  int sf_flag = NOT_SET;
  double fac = 2.0;
  bool OK = accuracyOK(0, tol, num_in, avd);
  int nfac = 2;
  if (OK)
    sf_flag = ACCURACY_OK;
  else
    {
      OK = accuracyOK(min_point, tol, num_in2, avd);
      if (OK)
	sf_flag = ANGULAR_DEVIATION;
      else
	{
	  OK = accuracyOK(min_point, fac*tol, num_in2, avd);
	  if (OK)
	    sf_flag = ACCURACY_POOR;
	}
    }
  
  if ((sf_flag == ANGULAR_DEVIATION ||
       (num_in2 >= nfac*num_in && sf_flag < ACCURACY_POOR)) && type_cyl)
    sf_flag = PROBABLE_HELIX;

  return sf_flag;
}
 
//===========================================================================
int RevEngRegion::defineSfFlag(int num_points, int min_point, double tol, int num_in, 
			       int num_in2, double avd, bool type_cyl)
//===========================================================================
{
  int sf_flag = NOT_SET;
  double fac = 2.0;
  bool OK = accuracyOK(num_points, min_point, tol, num_in, avd);
  int nfac = 2;
  if (OK)
    sf_flag = ACCURACY_OK;
  else
    {
      OK = accuracyOK(num_points, min_point, tol, num_in2, avd);
      if (OK)
	sf_flag = ANGULAR_DEVIATION;
      else
	{
	  OK = accuracyOK(num_points, min_point, fac*tol, num_in2, avd);
	  if (OK)
	    sf_flag = ACCURACY_POOR;
	}
    }
  
  if ((sf_flag == ANGULAR_DEVIATION ||
       (num_in2 >= nfac*num_in && sf_flag < ACCURACY_POOR)) && type_cyl)
    sf_flag = PROBABLE_HELIX;

  return sf_flag;
}
 
//===========================================================================
void RevEngRegion::updateInfo(double tol, double angtol)
//===========================================================================
{
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  MAH_ = MAK_ = avH_ = avK_ = 0.0;
  bbox_ = BoundingBox(3);
  double fac = 1.0/(double)group_points_.size();
  for  (size_t kj=0; kj<group_points_.size(); ++kj)
    {
      double k1 = group_points_[kj]->minPrincipalCurvature();
      double k2 = group_points_[kj]->maxPrincipalCurvature();
      double H = group_points_[kj]->meanCurvature();
      double K = group_points_[kj]->GaussCurvature();
      mink1_ = std::min(mink1_, fabs(k1));
      maxk1_ = std::max(maxk1_, fabs(k1));
      mink2_ = std::min(mink2_, fabs(k2));
      maxk2_ = std::max(maxk2_, fabs(k2));
      avH_ += fac*H;
      avK_ += fac*K;
      MAH_ += fac*fabs(H);
      MAK_ += fac*fabs(K);
      Vector3D point = group_points_[kj]->getPoint();
      Point point2(point[0], point[1], point[2]);
      bbox_.addUnionWith(point2);
    }

  if (group_points_.size() > 0)
    {
      normalcone_ = DirectionCone(group_points_[0]->getLocFuncNormal());
      normalcone2_ = DirectionCone(group_points_[0]->getTriangNormal());
      avnorm_ = Point(0.0, 0.0, 0.0);
      avnorm2_ = Point(0.0, 0.0, 0.0);
      for  (size_t kj=0; kj<group_points_.size(); ++kj)
	{
	  Point norm = group_points_[kj]->getLocFuncNormal();
	  normalcone_.addUnionWith(norm);
	  avnorm_ += fac*norm;
	  Point norm2 = group_points_[kj]->getTriangNormal();
	  normalcone2_.addUnionWith(norm2);
	  avnorm2_ += fac*norm2;
	}
      // (void)avnorm_.normalize_checked();
      // (void)avnorm2_.normalize_checked();

      double anglim = 0.2;
      if (normalcone_.angle() <= anglim)
	frac_norm_in_ = 1.0;
      else
	{
	  int nmb_in = 0;
	  for (size_t kr=0; kr<group_points_.size(); ++kr)
	    {
	      Point normal = group_points_[kr]->getLocFuncNormal();
	      if (avnorm_.angle(normal) <= anglim)
		nmb_in++;
	    }
	  frac_norm_in_ = (double)nmb_in/(double)group_points_.size();
	}
  
      if (normalcone2_.angle() <= anglim)
	frac_norm_in2_ = 1.0;
      else
	{
	  int nmb_in = 0;
	  for (size_t kr=0; kr<group_points_.size(); ++kr)
	    {
	      Point normal = group_points_[kr]->getTriangNormal();
	      if (avnorm2_.angle(normal) <= anglim)
		nmb_in++;
	    }
	  frac_norm_in2_ = (double)nmb_in/(double)group_points_.size();
	}
    }
  
  if (hasSurface() && group_points_.size() > 0)
    {
      if (tol > 0.0)
	{
	  maxdist_ = avdist_ = 0.0;
	  num_inside_ = num_inside2_ = 0;
	}
      Vector2D par = group_points_[0]->getPar();
      domain_[0] = domain_[1] = par[0];
      domain_[2] = domain_[3] = par[1];
      for  (size_t kj=0; kj<group_points_.size(); ++kj)
	{
	  Vector2D par = group_points_[kj]->getPar();
	  domain_[0] = std::min(domain_[0], par[0]);
	  domain_[1] = std::max(domain_[1], par[0]);
	  domain_[2] = std::min(domain_[2], par[1]);
	  domain_[3] = std::max(domain_[3], par[1]);
	  if (tol > 0.0)
	    {
	      double dist, ang;
	      group_points_[kj]->getSurfaceDist(dist, ang);
	      maxdist_ = std::max(maxdist_, dist);
	      avdist_ += fac*dist;
	      if (dist <= tol)
		{
		  num_inside2_++;
		  if (ang <= angtol)
		    num_inside_++;
		}
	    }
 	}
      if (tol > 0.0)
	{
	  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
	  bool cyllike = (surf->instanceType() == Class_Cylinder ||
			  surf->instanceType() == Class_Cone);
	  surfflag_ = defineSfFlag(0, tol, num_inside_, num_inside2_, avdist_,
				   cyllike);
	}
    }
}

//===========================================================================
void RevEngRegion::addPoint(RevEngPoint* point)
//===========================================================================
{
#ifdef DEBUG_CHECK
  if (std::find(group_points_.begin(), group_points_.end(), point) != group_points_.end())
    std::cout << "addPoint: point exists already. " << point << std::endl;
#endif
  int nmb = (int)group_points_.size();
  group_points_.push_back(point);
  point->setRegion(this);
  Vector3D point2 = point->getPoint();
  Point point3(point2[0], point2[1], point2[2]);
  bbox_.addUnionWith(point3);
  Point norm = point->getLocFuncNormal();
  normalcone_.addUnionWith(norm);
  Point norm2 = point->getTriangNormal();
  normalcone2_.addUnionWith(norm2);
  double k1 = point->minPrincipalCurvature();
  double k2 = point->maxPrincipalCurvature();
  mink1_ = std::min(mink1_, fabs(k1));
  maxk1_ = std::max(maxk1_, fabs(k1));
  mink2_ = std::min(mink2_, fabs(k2));
  maxk2_ = std::max(maxk2_, fabs(k2));
  avnorm_ = (nmb*avnorm_ + norm)/(double)(nmb+1);
  double H = point->meanCurvature();
  double K = point->GaussCurvature();
  avH_ = (nmb*avH_ + H)/(double)(nmb+1);
  avK_ = (nmb*avK_ + K)/(double)(nmb+1);
  MAH_ = (nmb*MAH_ + fabs(H))/(double)(nmb+1);
  MAK_ = (nmb*MAK_ + fabs(K))/(double)(nmb+1);

  if (hasSurface())
    {
      Vector2D par = point->getPar();
      domain_[0] = std::min(domain_[0], par[0]);
      domain_[1] = std::max(domain_[1], par[0]);
      domain_[2] = std::min(domain_[2], par[1]);
      domain_[3] = std::max(domain_[3], par[1]);
    }
 
}

//===========================================================================
void RevEngRegion::computeDomain()
//===========================================================================
{
  if (hasSurface() && group_points_.size() > 0)
    {
      Vector2D par = group_points_[0]->getPar();
      domain_[0] = domain_[1] = par[0];
      domain_[2] = domain_[3] = par[1];
      for  (size_t kj=1; kj<group_points_.size(); ++kj)
	{
	  Vector2D par = group_points_[kj]->getPar();
	  domain_[0] = std::min(domain_[0], par[0]);
	  domain_[1] = std::max(domain_[1], par[0]);
	  domain_[2] = std::min(domain_[2], par[1]);
	  domain_[3] = std::max(domain_[3], par[1]);
 	}
    }
}

//===========================================================================
void RevEngRegion::removePoint(RevEngPoint* point)
//===========================================================================
{
  auto found = std::find(group_points_.begin(), group_points_.end(), point);
  if (found != group_points_.end())
    {
      std::swap(*found, group_points_[group_points_.size()-1]);
      group_points_.pop_back();
    }
  // for (size_t ki=0; ki<group_points_.size(); ++ki)
  //   if (group_points_[ki] == point)
  //     {
  // 	group_points_.erase(group_points_.begin()+ki);
  // 	break;
  //     }
}

//===========================================================================
void RevEngRegion::updateWithPointsInOut(vector<RevEngPoint*>& points_out,
					 vector<RevEngPoint*>& points_in,
					 double tol, double angtol)
//===========================================================================
{
  for (size_t ki=0; ki<points_out.size(); ++ki)
    {
      removePoint(points_out[ki]);
      points_out[ki]->unsetRegion();
    }

  shared_ptr<ParamSurface> surf;
  bool cyllike = false;
  if (hasSurface())
    {
      surf = getSurface(0)->surface();
      cyllike = (surf->instanceType() == Class_Cylinder ||
		 surf->instanceType() == Class_Cone);
      double maxd, avd;
      int num_inside, num2_inside;
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      RevEngUtils::distToSurf(points_in.begin(), points_in.end(),
			      surf, tol, maxd, avd, num_inside, num2_inside,
			      inpt, outpt, parvals, dist_ang, angtol);
      for (size_t ki=0; ki<points_in.size(); ++ki)
	{
	  points_in[ki]->setPar(Vector2D(parvals[2*ki],parvals[2*ki+1]));
	  points_in[ki]->setSurfaceDist(dist_ang[ki].first, dist_ang[ki].second);
	}
    }
  for (size_t ki=0; ki<points_in.size(); ++ki)
    addPoint(points_in[ki]);

  updateInfo(tol, angtol);
  int sf_flag = NOT_SET;
  if (hasSurface())
    {
      double maxd2, avd2;
      int num_in2, num2_in2;
      getAccuracy(maxd2, avd2, num_in2, num2_in2);
      sf_flag = defineSfFlag(0, tol, num_in2, num2_in2, avd2, cyllike);
      setSurfaceFlag(sf_flag);
    }
}

//===========================================================================
void RevEngRegion::sortBlendPoints(vector<RevEngPoint*>& points,
				   vector<shared_ptr<CurveOnSurface> >& cvs,
				   double distance, bool in_blend,
				   vector<RevEngPoint*>& blend_points)
//===========================================================================
{
  int num_points = (int)points.size();
  for (int ki=num_points-1; ki>=0; --ki)
    {
      double tpar, dist;
      Point close;
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      //int ix = -1;
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  cvs[kj]->closestPoint(pt, cvs[kj]->startparam(), cvs[kj]->endparam(),
				tpar, close, dist);
	  if ((in_blend && dist < distance) ||
	      (in_blend == false && dist >= distance))
	    {
	      blend_points.push_back(points[ki]);
	      std::swap(points[ki],points[num_points-1]);
	      --num_points;
	    }
	}
    }

  if (num_points < (int)points.size())
    points.erase(points.begin()+num_points, points.end());
}

//===========================================================================
void RevEngRegion::sortBlendPoints(vector<RevEngPoint*>& points,
				   vector<shared_ptr<CurveOnSurface> >& cvs,
				   double distance, RevEngRegion* other,
				   vector<RevEngPoint*>& blend_points1,
				   vector<RevEngPoint*>& blend_points2)
//===========================================================================
{
  if ((!hasSurface()) || (!other->hasSurface()))
    return;  // Not as expected

  double eps = 1.0e-6;
  int num_points = (int)points.size();
  shared_ptr<ParamSurface> surf1 = getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = other->getSurface(0)->surface();
  double upar1, upar2, vpar1, vpar2, dist1, dist2;
  Point close1, close2;
  for (int ki=num_points-1; ki>=0; --ki)
    {
      double tpar, dist;
      Point close;
      Vector3D xyz = points[ki]->getPoint();
      Point pt(xyz[0], xyz[1], xyz[2]);
      int ix = -1;
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  cvs[kj]->closestPoint(pt, cvs[kj]->startparam(), cvs[kj]->endparam(),
				tpar, close, dist);
	  if (dist >= distance)
	    {
	      // Check which adjacent region the point belongs to
	      surf1->closestPoint(pt, upar1, vpar1, close1, dist1, eps);
	      surf2->closestPoint(pt, upar2, vpar2, close2, dist2, eps);
	      if (dist1 <= dist2)
		blend_points1.push_back(points[ki]);
	      else
		blend_points2.push_back(points[ki]);
	      std::swap(points[ki],points[num_points-1]);
	      --num_points;
	    }
	}
    }

  if (num_points < (int)points.size())
    points.erase(points.begin()+num_points, points.end());
}


//===========================================================================
void RevEngRegion::setRegionAdjacency()
//===========================================================================
{
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg == this)
	    continue;
	  if (adj_reg)
	    {
	      adj_reg->addAdjacentRegion(this);
	      addAdjacentRegion(adj_reg);
	    }
	}
    }
}

//===========================================================================
void RevEngRegion::updateRegionAdjacency()
//===========================================================================
{
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    (*it)->removeAdjacentRegion(this);
  adjacent_regions_.clear();

  setRegionAdjacency();
}

struct integrate_info
{
  RevEngRegion *adjacent;
  int nmb_pt_adj;
  double maxd, avd, maxd_adj, avd_adj;
  int nmb_in, nmb_in_adj;
  double min_ang, max_ang, min_dist, max_dist;
  bool outlier;

  integrate_info()
  {
    adjacent = 0;
    nmb_pt_adj = nmb_in = nmb_in_adj = 0;
    maxd = avd = maxd_adj = avd_adj = 0.0;
    min_ang = max_ang = min_dist = max_dist = -1.0;
    outlier = false;
  }

  void setAdjacent(RevEngRegion* adj)
  {
    adjacent = adj;
  }

  void setOutlier()
  {
    outlier = true;
  }

  void setNmbAdj(int nmb)
  {
    nmb_pt_adj = nmb;
  }

  void setAccuracy(double maxd1, double avd1, int nmb_in1, double maxd2,
		   double avd2, int nmb_in2)
  {
    maxd = maxd1;
    avd = avd1;
    nmb_in = nmb_in1;
    maxd_adj = maxd2;
    avd_adj = avd2;
    nmb_in_adj = nmb_in2;
  }

  void setMinMax(double min_a, double max_a, double min_d, double max_d)
  {
    min_ang = min_a;
    max_ang = max_a;
    min_dist = min_d;
    max_dist = max_d;
  }
};
  
//===========================================================================
bool RevEngRegion::integrateInAdjacent(double mean_edge_len, int min_next,
				       int max_next, double tol, double angtol,
				       int max_nmb_outlier,
				       RevEngRegion* taboo)
//===========================================================================
{
  // if (adjacent_regions_.size() != 1)
  //   return false;   // To be removed

  if ((int)group_points_.size() > max_next/2)
    return false;

#ifdef DEBUG_INTEGRATE
  std::ofstream of("curr_integrate.g2");
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kh=0; kh<group_points_.size(); ++kh)
    of << group_points_[kh]->getPoint() << std::endl;
#endif 
  size_t kj=0;
  double local_len = group_points_[0]->getMeanEdgLen(10.0*mean_edge_len);
  double radius = 3.0*(local_len + mean_edge_len);
  radius = std::min(radius, 20.0*mean_edge_len);
  size_t adjsize = adjacent_regions_.size();
  vector<integrate_info> info(adjsize);
  
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it, ++kj)
    info[kj].setAdjacent(*it);
  
  if ((int)group_points_.size() <= max_nmb_outlier)
    {
      // Simplified test
      vector<RevEngRegion*> adj_reg;
      vector<pair<double,double> > adj_info;
      double lentol = 2.0*mean_edge_len;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  vector<RevEngRegion*> pt_adj_reg;
	  vector<RevEngPoint*> pt_adj_pt;
	  group_points_[ki]->getAdjInfo(mean_edge_len, pt_adj_reg, pt_adj_pt);
	  Point monge1 = group_points_[ki]->getLocFuncNormal();
	  for (size_t kj=0; kj<pt_adj_pt.size(); ++kj)
	    {
	      if (pt_adj_reg[kj] == this)
		continue;
	      double len = group_points_[ki]->pntDist(pt_adj_pt[kj]);
	      if (len > lentol)
		continue;
	      Point monge2 = pt_adj_pt[kj]->getLocFuncNormal();
	      if (monge1*monge2 < 0.0 || monge1.angle(monge2) > angtol)
		continue;
	      adj_reg.push_back(pt_adj_reg[kj]);
	      adj_info.push_back(std::make_pair(len, monge1.angle(monge2)));
	    }
	}

      for (size_t kj=0; kj<info.size(); ++kj)
	{
	  RevEngRegion *reg2 = info[kj].adjacent;
	  for (size_t ki=0; ki<adj_reg.size(); ++ki)
	    {
	      if (adj_reg[ki] == reg2)
		{
		  if (info[kj].max_dist < 0.0)
		    {
		      info[kj].min_dist = info[kj].max_dist = adj_info[ki].first;
		      info[kj].min_ang = info[kj].max_ang = adj_info[ki].second;
		    }
		  else
		    {
		      info[kj].min_dist = std::min(info[kj].min_dist,adj_info[ki].first);
		      info[kj].max_dist = std::max(info[kj].max_dist,adj_info[ki].first);
		      info[kj].min_ang = std::min(info[kj].min_ang,adj_info[ki].second);
		      info[kj].max_ang = std::max(info[kj].max_ang,adj_info[ki].second);
		    }
		}
	    }
	}
    }

  kj = 0;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it, ++kj)
    {
      double maxd1, maxd2, avd1, avd2;
      maxd1 = maxd2 = avd1 = avd2 = std::numeric_limits<double>::max();
      int nmb_in1 = 0, nmb_in1_2 = 0, nmb_in2 = 0;
      bool local_approx = (info[kj].max_dist < 0.0);
      bool outlier = false;
      int nmb_pt_adj;
      bool computed = computeIntegrateInfo(group_points_, *it, tol, angtol, radius,
					   local_approx, min_next, max_next, max_nmb_outlier, 
					   outlier, nmb_pt_adj, maxd2, avd2, nmb_in2, maxd1, 
					   avd1, nmb_in1, nmb_in1_2);
      info[kj].setNmbAdj(nmb_pt_adj);
      if (outlier)
	info[kj].setOutlier();
      if (!computed)
	continue;
      // shared_ptr<ParamSurface> surf;
      // if ((*it)->hasSurface())
      // 	{
      // 	  surf = (*it)->getSurface(0)->surface();
      // 	  (*it)->getAccuracy(maxd1, avd1, nmb_in1);
      // 	}
      // else if ((*it)->hasBaseSf())
      // 	(*it)->getBase(surf, maxd1, avd1, nmb_in1);

      // if (surf.get())
      // 	{
      // 	  // Fetch accuracy of current surface
      // 	  info[kj].setNmbAdj((*it)->numPoints());
	  
      // 	  // Check accuracy of new points
      // 	  vector<RevEngPoint*> in2, out2;
      // 	  vector<pair<double,double> > dist_ang2;
      // 	  vector<double> parvals2;
      // 	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
      // 				  surf, tol, maxd2, avd2, nmb_in2, in2, out2,
      // 				  parvals2, dist_ang2, angtol);
      // 	}
      // else if (info[kj].max_dist < 0.0)
      // 	{
      // 	  // Fetch nearby points
      // 	  vector<RevEngPoint*> nearpts;
      // 	  if ((*it)->numPoints() < min_next)
      // 	    nearpts.insert(nearpts.end(), (*it)->pointsBegin(), (*it)->pointsEnd());
      // 	  else
      // 	    {
      // 	      for (size_t ki=0; ki<group_points_.size(); ++ki)
      // 		{
      // 		  nearpts.clear();
      // 		  group_points_[ki]->fetchClosePoints2(radius, min_next,
      // 						       max_next-(int)group_points_.size(),
      // 						       nearpts, *it);
      // 		  if (nearpts.size() > max_nmb_outlier)
      // 		    break;
      // 		}
      // 	    }
      // 	  size_t nearnmb = nearpts.size();
      // 	  info[kj].setNmbAdj((int)nearnmb);
      // 	  if (((int)(nearnmb+group_points_.size()) <= max_nmb_outlier) ||
      // 	      ((int)nearnmb <= max_nmb_outlier && (*it)->numPoints() > min_next))
      // 	    {
      // 	      if ((int)group_points_.size() <= max_nmb_outlier)
      // 		info[kj].setOutlier();
      // 	      continue;
      // 	    }

      // 	  nearpts.insert(nearpts.end(), group_points_.begin(), group_points_.end());
      // 	  BoundingBox bbox = bbox_;
      // 	  for (size_t ki=0; ki<nearnmb; ++ki)
      // 	    {
      // 	      Vector3D xyz = nearpts[ki]->getPoint();
      // 	      Point xyz2(xyz[0], xyz[1], xyz[2]);
      // 	      bbox.addUnionWith(xyz2);
      // 	    }
      // 	  surf = surfApprox(nearpts, bbox);

      // 	  // Check accuracy
      // 	  vector<RevEngPoint*> in1, in2, out1, out2;
      // 	  vector<pair<double,double> > dist_ang1, dist_ang2;
      // 	  vector<double> parvals1, parvals2;
      // 	  RevEngUtils::distToSurf(nearpts.begin(), nearpts.begin()+nearnmb, surf,
      // 				  tol, maxd1, avd1, nmb_in1, in1, out1,
      // 				  parvals1, dist_ang1, angtol);
      // 	  RevEngUtils::distToSurf(nearpts.begin()+nearnmb, nearpts.end(), surf,
      // 				  tol, maxd2, avd2, nmb_in2, in2, out2,
      // 				  parvals2, dist_ang2, angtol);
      // 	}
#ifdef DEBUG_INTEGRATE
      int num = (*it)->numPoints();
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of << (*it)->getPoint(ka)->getPoint() << std::endl;

      // if (surf.get())
      // 	{
      // 	  surf->writeStandardHeader(of);
      // 	  surf->write(of);
      // 	}
#endif
      info[kj].setAccuracy(maxd2, avd2, nmb_in2, maxd1, avd1, nmb_in1);
    }

  // Select adjacent
  bool outlier = true;
  int ix = -1;
  double score = 0.0;
  int num = (int)group_points_.size();
  for (size_t kj=0; kj<info.size(); ++kj)
    {
      if (!info[kj].outlier)
	outlier = false;
      if (taboo && info[kj].adjacent == taboo /*&& info.size()>1*/)
	continue;
      if (prev_region_ && info[kj].adjacent == prev_region_ && info.size()>1)
	continue;
      if (info[kj].avd > tol || (info[kj].nmb_in < num/2 && info[kj].maxd > tol))
	continue;
      if (info[kj].adjacent->hasAssociatedBlend())
	continue;  // Not to be extended
      double frac1 = (double)info[kj].nmb_in/(double)num;
      double frac2 = (double)info[kj].nmb_in_adj/(double)info[kj].nmb_pt_adj;
      if (frac1 < 0.9*frac2 && info[kj].maxd > tol)
	continue;
      double avH = avH_*info[kj].adjacent->avH_;
      double avK = avK_*info[kj].adjacent->avK_;
      double curr_score = (tol/std::max(1.0e-6, info[kj].avd)) + /* frac2/frac1 +*/
	(info[kj].adjacent->hasSurface()) + (avH > 0.0) + (avK > 0.0);
      if (curr_score > score)
	{
	  ix = (int)kj;
	  score = curr_score;
	}
    }

  if (ix < 0 && (!outlier))
    {
      double div = std::numeric_limits<double>::max();
      for (size_t kj=0; kj<info.size(); ++kj)
	{
	  if (info[kj].max_dist < 0)
	    continue;
	  if (taboo && info[kj].adjacent == taboo /*&& info.size()>1*/)
	    continue;
	  if (prev_region_ && info[kj].adjacent == prev_region_ && info.size()>1)
	    continue;
	  double curr_div = info[kj].min_dist + info[kj].min_ang;
	  if (curr_div < div)
	    {
	      div = curr_div;
	      ix = (int)kj;
	    }
	}
    }
  
  if (outlier)
    {
     for (size_t ki=0; ki<group_points_.size(); ++ki)
       {
	group_points_[ki]->unsetRegion();
	group_points_[ki]->setOutlier();
	group_points_[ki]->addMove();
       }
     for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
       (*it)->removeAdjacentRegion(this);
     //updateInfo();
     return true;
    }
  else if (ix >= 0)
    {
      if (info[ix].adjacent->hasSurface())
	{
	  // Set parameter value, distance and angle difference
	  shared_ptr<ParamSurface> surf = info[ix].adjacent->getSurface(0)->surface();
	  double maxd2, avd2;
	  int nmb_in2, nmb2_in2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double,double> > dist_ang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  surf, tol, maxd2, avd2, nmb_in2, nmb2_in2, in2, out2,
				  parvals2, dist_ang2, angtol);
	  bool cyllike = (surf->instanceType() == Class_Cylinder ||
			  surf->instanceType() == Class_Cone);
	  int sf_flag = defineSfFlag(0, tol, nmb_in2, nmb2_in2, avd2, cyllike);
	  for (size_t ki=0; ki<group_points_.size(); ++ki)
	    {
	      group_points_[ki]->setPar(Vector2D(parvals2[2*ki],parvals2[2*ki+1]));
	      group_points_[ki]->setSurfaceDist(dist_ang2[ki].first, dist_ang2[ki].second);
	    }
	  setSurfaceFlag(sf_flag);
	}
     for (size_t ki=0; ki<group_points_.size(); ++ki)
       group_points_[ki]->addMove();
      vector<pair<double, double> > dummy;
      info[ix].adjacent->addRegion(this, dummy);
      for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
	(*it)->removeAdjacentRegion(this);
      info[ix].adjacent->updateInfo(tol, angtol);
      return true;
      }
  return false;
}

//===========================================================================
bool RevEngRegion::computeIntegrateInfo(vector<RevEngPoint*>& points, RevEngRegion *adj_reg,
					double tol, double angtol, double radius,
					bool local_approx, int min_next, int max_next,
					int max_nmb_outlier, bool& outlier, int& nmb_pt_adj,
					double& maxdist, double& avdist, int& nmb_in,
					double& maxdist_adj, double& avdist_adj,
					int& nmb_in_adj, int& nmb_in_adj2)
//===========================================================================
{
  outlier = false;
  nmb_pt_adj = 0;
  
  shared_ptr<ParamSurface> surf;
  if (adj_reg->hasSurface())
    {
      surf = adj_reg->getSurface(0)->surface();
      adj_reg->getAccuracy(maxdist_adj, avdist_adj, nmb_in_adj, nmb_in_adj2);
    }
  else if (adj_reg->hasBaseSf())
    {
      adj_reg->getBase(surf, maxdist_adj, avdist_adj, nmb_in_adj, nmb_in_adj2);
      nmb_in_adj2 = nmb_in_adj;
    }

  if (surf.get())
    {
      // Fetch accuracy of current surface
      nmb_pt_adj = adj_reg->numPoints();
	  
      // Check accuracy of new points
      vector<RevEngPoint*> in, out;
      vector<pair<double,double> > dist_ang;
      vector<double> parvals;
      int nmb2_in;
      RevEngUtils::distToSurf(points.begin(), points.end(),
			      surf, tol, maxdist, avdist, nmb_in, nmb2_in, in, out,
			      parvals, dist_ang, angtol);
    }
  else if (local_approx)
    {
      // Fetch nearby points
      vector<RevEngPoint*> nearpts;
      if (adj_reg->numPoints() < min_next)
	nearpts.insert(nearpts.end(), adj_reg->pointsBegin(), adj_reg->pointsEnd());
      else
	{
	  for (size_t ki=0; ki<points.size(); ++ki)
	    {
	      nearpts.clear();
	      points[ki]->fetchClosePoints2(radius, min_next,
					    max_next-(int)points.size(),
					    nearpts, adj_reg);
	      if ((int)nearpts.size() > max_nmb_outlier)
		break;
	    }
	}
      size_t nearnmb = nearpts.size();
      nmb_pt_adj = (int)nearnmb;
      if (((int)(nearnmb+points.size()) <= max_nmb_outlier) ||
	  ((int)nearnmb <= max_nmb_outlier && adj_reg->numPoints() > min_next))
	{
	  if ((int)points.size() <= max_nmb_outlier)
	    outlier = true;
	  return false;
	}

      nearpts.insert(nearpts.end(), points.begin(), points.end());
      BoundingBox bbox = bbox_;
      for (size_t ki=0; ki<nearnmb; ++ki)
	{
	  Vector3D xyz = nearpts[ki]->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  bbox.addUnionWith(xyz2);
	}
      surf = surfApprox(nearpts, bbox);

      // Check accuracy
      vector<RevEngPoint*> in1, in2, out1, out2;
      vector<pair<double,double> > dist_ang1, dist_ang2;
      vector<double> parvals1, parvals2;
      int nmb2_in, nmb2_in_adj;
      RevEngUtils::distToSurf(nearpts.begin(), nearpts.begin()+nearnmb, surf,
			      tol, maxdist_adj, avdist_adj, nmb_in_adj,
			      nmb2_in_adj, in1, out1,
			      parvals1, dist_ang1, angtol);
      RevEngUtils::distToSurf(nearpts.begin()+nearnmb, nearpts.end(), surf,
			      tol, maxdist, avdist, nmb_in, nmb2_in, in2, out2,
			      parvals2, dist_ang2, angtol);
    }

  return true;
}

//===========================================================================
bool RevEngRegion::adjustWithCylinder(Point mainaxis[3],
				      double tol, double angtol, int min_pt_reg,
				      vector<vector<RevEngPoint*> >& out_groups,
				      vector<RevEngRegion*>& grown_regions,
				      vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (!hasSurface())
    return false;

  HedgeSurface *hedge = getSurface(0);
  int code;
  ClassType type = hedge->instanceType(code);
  if (type != Class_Cylinder &&
      (!(type == Class_SplineSurface && code == LINEARSWEPT_SURF)))
    return false;
  shared_ptr<ParamSurface> surf = hedge->surface();

#ifdef DEBUG_ADJUST
  std::ofstream of("cylinder_adjust.g2");
  writeRegionInfo(of);
#endif
  
  // Make bounding parameter domain
  int ixd = (type == Class_Cylinder) ? 0 : 2;
  //int ixd2 = (ixd == 0) ? 2 : 0;
  int ixp = ixd/2;
  double dom[4];
  Vector2D uv = group_points_[0]->getPar();
  dom[0] = dom[1] = uv[0];
  dom[2] = dom[3] = uv[1];
  double dist, ang, avang;
  double fac = 1.0/(double)group_points_.size();
  group_points_[0]->getSurfaceDist(dist, ang);
  avang = fac*ang;
  for (size_t ki=1; ki<group_points_.size(); ++ki)
    {
      Vector2D uv = group_points_[ki]->getPar();
      group_points_[ki]->getSurfaceDist(dist, ang);
      dom[0] = std::min(dom[0], uv[0]);
      dom[1] = std::max(dom[1], uv[0]);
      dom[2] = std::min(dom[2], uv[1]);
      dom[3] = std::max(dom[3], uv[1]);
      avang += fac*ang;
    }

#ifdef DEBUG_ADJUST
  std::ofstream ofp("projected_pts.g2");
#endif
  
  Point axis, pnt;
  shared_ptr<ParamCurve> crv;
  if (type == Class_Cylinder)
    {
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
      axis = cyl->getAxis();
      pnt = cyl->getLocation();
      shared_ptr<Circle> circ(new Circle(cyl->radius(0,0), pnt, axis, cyl->direction2()));
      crv = circ;
    }
  else
    {
      vector<Point> der(3);
      double upar = 0.5*(dom[0]+dom[1]);
      double vpar = 0.5*(dom[2]+dom[3]);
      surf->point(der, upar, vpar, 1);
      axis = der[1];
      axis.normalize();
      pnt = der[0];
      vector<shared_ptr<ParamCurve> > cvs = surf->constParamCurves(vpar, true);
      crv = cvs[0];
    }
  
  vector<Point> projected;
  double maxdp, avdp;
  RevEngUtils::projectToPlane(group_points_, axis, pnt, projected, maxdp, avdp);
#ifdef DEBUG_ADJUST
  ofp << "400 1 0 4 255 0 0 255" << std::endl;
  ofp << projected.size() << std::endl;
  for (size_t kr=0; kr<projected.size(); ++kr)
    ofp << projected[kr] << std::endl;
  crv->writeStandardHeader(ofp);
  crv->write(ofp);
#endif
  
  // Reduce domain from the start
  vector<RevEngPoint*> outer;
  double del = 0.02;
  double del2 = del*(dom[ixd+1] - dom[ixd]);
  double par;
  int knmb = 10;
  int ka;
  double dfac = 2.0;
  //double afac = 2.0;
  double pfac = 0.5;
  int part = (int)(del*(double)group_points_.size());
  for (ka=0, par=dom[ixd]+del2; ka<knmb; ++ka, par+=del2)
    {
      vector<RevEngPoint*> curr_pts;
      double avd = 0.0, ava = 0.0;
      int nn = 0;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  if (uv[ixp] >= par-del2 && uv[ixp]<par)
	    {
	      curr_pts.push_back(group_points_[ki]);
	      group_points_[ki]->getSurfaceDist(dist, ang);
	      avd += dist;
	      ava += ang;
	      ++nn;
	    }
	}
      avd /= (double)nn;
      ava /= (double)nn;

      // Check with adjacent regions
#ifdef DEBUG_ADJUST
      std::ofstream ofo("outer_cand.g2");
      ofo << "400 1 0 4 0 255 0 255" << std::endl;
      ofo << curr_pts.size() << std::endl;
     for (size_t ki=0; ki<curr_pts.size(); ++ki)
       ofo << curr_pts[ki]->getPoint() << std::endl;
#endif
     
      vector<RevEngRegion*> next_reg;
      vector<int> nmb_next;
      for (size_t ki=0; ki<curr_pts.size(); ++ki)
	{
	  vector<RevEngRegion*> adjr;
	  curr_pts[ki]->adjacentRegions(adjr);
	  for (size_t kj=0; kj<adjr.size(); ++kj)
	    {
	      size_t kr=0;
	      for (kr=0; kr<next_reg.size(); ++kr)
		if (next_reg[kr] == adjr[kj])
		  break;
	      if (kr == next_reg.size())
		{
		  next_reg.push_back(adjr[kj]);
		  nmb_next.push_back(1);
		}
	      else
		nmb_next[kr]++;
	    }
	}

      
      if (avd > dfac*avdist_ && (ava > dfac*avang || nn < pfac*part))
	{
	  outer.insert(outer.end(), curr_pts.begin(), curr_pts.end());
	  dom[ixd] = par;
	}
      else
	break;
      int stop_break = 1;
    }

  // Reduce domain from the end
  for (ka=0, par=dom[ixd+1]-del2; ka<knmb; ++ka, par-=del2)
    {
      vector<RevEngPoint*> curr_pts;
       double avd = 0.0, ava = 0.0;
      int nn = 0;
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  Vector2D uv = group_points_[ki]->getPar();
	  if (uv[ixp] > par && uv[ixp]<=par+del2)
	    {
	      curr_pts.push_back(group_points_[ki]);
	      group_points_[ki]->getSurfaceDist(dist, ang);
	      avd += dist;
	      ava += ang;
	      ++nn;
	    }
	}
      avd /= (double)nn;
      ava /= (double)nn;
      // Check with adjacent regions
#ifdef DEBUG_ADJUST
      std::ofstream ofo("outer_cand.g2");
      ofo << "400 1 0 4 0 255 0 255" << std::endl;
      ofo << curr_pts.size() << std::endl;
     for (size_t ki=0; ki<curr_pts.size(); ++ki)
       ofo << curr_pts[ki]->getPoint() << std::endl;
#endif
     
      vector<RevEngRegion*> next_reg;
      vector<int> nmb_next;
      for (size_t ki=0; ki<curr_pts.size(); ++ki)
	{
	  vector<RevEngRegion*> adjr;
	  curr_pts[ki]->adjacentRegions(adjr);
	  for (size_t kj=0; kj<adjr.size(); ++kj)
	    {
	      size_t kr=0;
	      for (kr=0; kr<next_reg.size(); ++kr)
		if (next_reg[kr] == adjr[kj])
		  break;
	      if (kr == next_reg.size())
		{
		  next_reg.push_back(adjr[kj]);
		  nmb_next.push_back(1);
		}
	      else
		nmb_next[kr]++;
	    }
	}

      
      if (avd > dfac*avdist_ && (ava > dfac*avang || nn < pfac*part))
	{
	  outer.insert(outer.end(), curr_pts.begin(), curr_pts.end());
	  dom[ixd+1] = par;
	}
      else
	break;
      int stop_break = 1;
    }

  // Integrate points from adjacent regions
  vector<vector<RevEngPoint*> > adjpts;
  vector<RevEngRegion*> adjreg;
  vector<vector<double> > par_and_dist;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
#ifdef DEBUG_ADJUST
      std::ofstream of2("adj_group.g2");
      int num = (*it)->numPoints();
      of2 << "400 1 0 4 255 0 0 255" << std::endl;
      of2 << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of2 << (*it)->getPoint(ka)->getPoint() << std::endl;
#endif
      
      vector<Point> projected2;
      double maxdp2, avdp2;
      vector<RevEngPoint*> curr_pts = (*it)->getPoints();
      RevEngUtils::projectToPlane(curr_pts, axis, pnt, projected2, maxdp2, avdp2);
#ifdef DEBUG_ADJUST
      std::ofstream ofp2("projected_pts_adj.g2");
      ofp2 << "400 1 0 4 100 155 0 255" << std::endl;
      ofp2 << projected2.size() << std::endl;
      for (size_t kr=0; kr<projected2.size(); ++kr)
	ofp2 << projected2[kr] << std::endl;
#endif
      
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> surf2((*it)->getSurface(0)->surface()->clone());
	  double upar2, vpar2, dist;
	  Point close;
	  surf2->closestPoint(pnt, upar2, vpar2, close, dist, tol);
	  if (!surf2->isBounded())
	    {
	      BoundingBox bb = getBbox();
	      double len = bb.low().dist(bb.high());
	      shared_ptr<Cylinder> elem1 =
		dynamic_pointer_cast<Cylinder,ParamSurface>(surf2);
	      shared_ptr<Plane> elem2 =
		dynamic_pointer_cast<Plane,ParamSurface>(surf2);
	      if (elem1.get())
		elem1->setParamBoundsV(-len, len);
	      else if (elem2.get())
		elem2->setParameterBounds(-len, -len, len, len);
	    }

#ifdef DEBUG_ADJUST
	  if (surf2->isBounded())
	    {
	      vector<shared_ptr<ParamCurve> > cvs2_1 = surf2->constParamCurves(upar2, false);
	      vector<shared_ptr<ParamCurve> > cvs2_2 = surf2->constParamCurves(vpar2, true);
	      cvs2_1[0]->writeStandardHeader(ofp2);
	      cvs2_1[0]->write(ofp2);
	      cvs2_2[0]->writeStandardHeader(ofp2);
	      cvs2_2[0]->write(ofp2);
	    }
#endif
	}

      vector<RevEngPoint*> curr_adjpts;
       vector<double> curr_par_and_dist;
       double avd, ava;
       int nn;
       getAdjInsideDist(surf, dom, tol, *it, avd, ava, nn, curr_adjpts, curr_par_and_dist);

      if (avd < dfac*avdist_ /*&& (ava < dfac*avang || nn == num)*/)
	{
	  adjpts.push_back(curr_adjpts);
	  adjreg.push_back(*it);
	  par_and_dist.push_back(curr_par_and_dist);
	}
     }

#ifdef DEBUG_ADJUST
  std::ofstream of2("move2pts.g2");
  for (size_t ki=0; ki<adjpts.size(); ++ki)
    {
      of2 << "400 1 0 4 75 75 75 255" << std::endl;
      of2 << adjpts[ki].size() << std::endl;
      for (size_t kh=0; kh<adjpts[ki].size(); ++kh)
	of2 << adjpts[ki][kh]->getPoint() << std::endl;
    }
#endif
  
  // Adjust point groups
  if (outer.size() > 0)
    {
      extractSpesPoints(outer, out_groups);
      splitRegion(out_groups);
      }

  for (size_t ki=0; ki<adjpts.size(); ++ki)
    {
      for (size_t kh=0; kh<adjpts[ki].size(); ++kh)
	{
	  adjpts[ki][kh]->setMoved();
	  adjreg[ki]->removePoint(adjpts[ki][kh]);
	  adjpts[ki][kh]->setRegion(this);
	  adjpts[ki][kh]->setPar(Vector2D(par_and_dist[ki][4*kh],
					  par_and_dist[ki][4*kh+1]));
	  adjpts[ki][kh]->setSurfaceDist(par_and_dist[ki][4*kh+2], par_and_dist[ki][4*kh+3]);
	  addPoint(adjpts[ki][kh]);
	}

      if (adjreg[ki]->numPoints() == 0)
	{
	  for (auto it=adjreg[ki]->adjacent_regions_.begin();
	       it!=adjreg[ki]->adjacent_regions_.end(); ++it)
	    {
	      if (*it != this)
		{
		  (*it)->addAdjacentRegion(this);
		  (*it)->removeAdjacentRegion(adjreg[ki]);
		}
	    }
	  grown_regions.push_back(adjreg[ki]);
	  int num_sf = adjreg[ki]->numSurface();
	  for (int kb=0; kb<num_sf; ++kb)
	    adj_surfs.push_back(adjreg[ki]->getSurface(kb));
	      
	  removeAdjacentRegion(adjreg[ki]);
	}
      else
	{
	  if (adjreg[ki]->hasSurface())
	    {
	      if (adjreg[ki]->numPoints() >= min_pt_reg)
		adjreg[ki]->checkReplaceSurf(mainaxis, min_pt_reg,
					     tol, angtol);
	      else
		{
		  int num_sf = adjreg[ki]->numSurface();
		  for (int kb=0; kb<num_sf; ++kb)
		    adj_surfs.push_back(adjreg[ki]->getSurface(kb));
		  adjreg[ki]->clearSurface();
		}
	    }
	  else
	    adjreg[ki]->updateInfo(tol, angtol);
	}
    }
  
  //updateInfo();
  checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);

   
  return (outer.size() > 0 || adjpts.size() > 0);
}

//===========================================================================
void RevEngRegion::getAdjInsideDist(shared_ptr<ParamSurface> surf, double dom[4],
				    double tol, RevEngRegion* reg,
				    double& avd, double& ava, int& nn,
				    vector<RevEngPoint*>& adjpts,
				    vector<double>& par_and_dist)
//===========================================================================
{
  double maxdd, avdd;
  int num_inside, num2_inside;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> in, out;
  RevEngUtils::distToSurf(reg->pointsBegin(),
			  reg->pointsEnd(), surf, tol, maxdd, avdd, 
			  num_inside, num2_inside, in, out, parvals, dist_ang);
  avd = 0.0;
  ava = 0.0;
  nn = 0;
  for (int kh=0; kh<(int)dist_ang.size(); ++kh)
    {
      if (parvals[2*kh] >= dom[0] && parvals[2*kh] <= dom[1] &&
	  parvals[2*kh+1] >= dom[2] && parvals[2*kh+1] <= dom[3])
	{
	  adjpts.push_back(reg->getPoint(kh));
	  par_and_dist.push_back(parvals[2*kh]);
	  par_and_dist.push_back(parvals[2*kh+1]);
	  par_and_dist.push_back(dist_ang[kh].first);
	  par_and_dist.push_back(dist_ang[kh].second);
	  avd += dist_ang[kh].first;
	  ava += dist_ang[kh].second;
	  ++nn;
	}
    }
  if (nn > 0)
    {
      avd /= (double)nn;
      ava /= (double)nn;
    }
}

//===========================================================================
void RevEngRegion::addRegion(RevEngRegion* reg,
			     vector<pair<double, double> >& dist_ang,
			     double maxd, double avd, int num_inside,
			     int num_inside2)
//===========================================================================
{
  int num = reg->numPoints();
#ifdef DEBUG_CHECK
  for (int ki=0; ki<num; ++ki)
    {
      auto it = std::find(group_points_.begin(), group_points_.end(), reg->getPoint(ki));
      if (it != group_points_.end())
	std::cout << "addRegion: point exists already. " << it-group_points_.begin() <<" " << reg << " ki= " << ki << " point: " << reg->getPoint(ki) << std::endl;
    }
#endif
  
  bbox_.addUnionWith(reg->boundingBox());
  normalcone_.addUnionWith(reg->getNormalCone());
  normalcone2_.addUnionWith(reg->getNormalConeTriang());
  if (num_inside >= 0)
    {
      maxdist_ = std::max(maxdist_, maxd);
      double div = (double)((int)group_points_.size() + num);
      avdist_ = ((double)(group_points_.size())*avdist_ + num*avd)/div;
      num_inside_ += num_inside;
      num_inside2_ += num_inside2;
    }

  double mink1, maxk1, mink2, maxk2;
  reg->getPrincipalCurvatureInfo(mink1, maxk1, mink2, maxk2);
  mink1_ = std::min(mink1_, mink1);
  maxk1_ = std::max(maxk1_, maxk1);
  mink2_ = std::min(mink2_, mink2);
  maxk2_ = std::max(maxk2_, maxk2);

  int num_all = numPoints() + reg->numPoints();
  double fac1 = (double)numPoints()/(double)num_all;
  double fac2 = (double)reg->numPoints()/(double)num_all;
  double avH, avK, MAH, MAK;
  reg->getAvCurvatureInfo(avH, avK, MAH, MAK);
  Point avnorm = reg->getMeanNormal();
  Point avnorm2 = reg->getMeanNormalTriang();
  avH_ = fac1*avH_ + fac2*avH;
  avK_ = fac1*avK_ + fac2*avK;
  MAH_ = fac1*MAH_ + fac2*MAH;
  MAK_ = fac1*MAK_ + fac2*MAK;
  avnorm_ = fac1*avnorm_ + fac2*avnorm;
  avnorm2_ = fac1*avnorm2_ + fac2*avnorm2;

  size_t kr=0;
  for (auto it=reg->pointsBegin(); it != reg->pointsEnd(); ++it, ++kr)
    {
      if (dist_ang.size() > 0)
	(*it)->setSurfaceDist(dist_ang[kr].first, dist_ang[kr].second);
      (*it)->setRegion(this);
    }
  group_points_.insert(group_points_.end(), reg->pointsBegin(),
		       reg->pointsEnd());
  removeAdjacentRegion(reg);

  if (hasSurface())
    computeDomain();
  
}


//===========================================================================
bool RevEngRegion::includeAdjacent(RevEngRegion* adj, Point mainaxis[3], 
				   double tol, double angtol,
				   vector<RevEngRegion*>& grown_regions,
				   vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  if (!hasSurface())
    return false;

  // Check accuracy
  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  double maxdist, avdist;
  int num_in, num2_in;
  vector<RevEngPoint*> inpt, outpt;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  RevEngUtils::distToSurf(adj->pointsBegin(), adj->pointsEnd(), surf,
			  tol, maxdist, avdist, num_in, num2_in, inpt, outpt,
			  parvals, dist_ang, angtol);
  
  bool type_cyl = (surf->instanceType() == Class_Cylinder ||
		   surf->instanceType() == Class_Cone);
  int num_pt_adj = adj->numPoints();
  int sf_flag = defineSfFlag(num_pt_adj, 0, tol, num_in, num2_in,
			     avdist, type_cyl);
  int adj_flag = adj->getSurfaceFlag();
  if (sf_flag < NOT_SET &&
      (sf_flag <= adj_flag ||
       (sf_flag == PROBABLE_HELIX && adj_flag == ANGULAR_DEVIATION) ||
       (num_in > num_pt_adj/2 && avdist <= tol)))
    {
      // Include
      vector<RevEngRegion*> added_adjacent;
      includeAdjacentRegion(adj, maxdist, avdist, num_in, num2_in,
			    parvals, dist_ang, added_adjacent);
      grown_regions.push_back(adj);
      int num_sf = adj->numSurface();
      for (int kb=0; kb<num_sf; ++kb)
	adj_surfs.push_back(adj->getSurface(kb));
      removeAdjacentRegion(adj);
      setSurfaceFlag(sf_flag);
      
      checkReplaceSurf(mainaxis, 0, tol, angtol);
      return true;
    }
  return false;
}

//===========================================================================
void RevEngRegion::growWithSurf(Point mainaxis[3], int min_pt_reg, double tol,
				double angtol, vector<RevEngRegion*>& grown_regions,
				vector<HedgeSurface*>& adj_surfs,
				vector<RevEngEdge*>& adj_edgs, bool use_base)
//===========================================================================
{
  //double eps = 1.0e-6;
  if (associated_sf_.size() == 0)
    return;  // No surface to grow

  if (hasAssociatedBlend())
    return;

  size_t num_group_points = group_points_.size();
  
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  int classtype = surf->instanceType();
  bool cyllike = (classtype == Class_Cylinder || classtype == Class_Cone);

  vector<RevEngRegion*> adj_reg;
  adj_reg.insert(adj_reg.end(), adjacent_regions_.begin(), adjacent_regions_.end());
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (adj_reg[kj]->numPoints() > adj_reg[ki]->numPoints())
	std::swap(adj_reg[ki], adj_reg[kj]);

  vector<grow_cand> cands;
  bool changed = false;
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      if ((adj_reg[ki]->prev_region_ && adj_reg[ki]->prev_region_ == this) ||
	  adj_reg[ki]->hasAssociatedBlend())
	{
	  continue;
	}

      int num = adj_reg[ki]->numPoints();
      
#ifdef DEBUG_GROW
      int write_extend = 1;
      std::ofstream of("curr_extend.g2");
       if (write_extend)
	{
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << group_points_.size() << std::endl;
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    of << group_points_[kh]->getPoint() << std::endl;
      
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << num << std::endl;
	  for (int ka=0; ka<num; ++ka)
	    of << adj_reg[ki]->getPoint(ka)->getPoint() << std::endl;
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
#endif

       double maxd, avd;
       int num_in, num2_in;
       vector<RevEngPoint*> in, out;
       vector<pair<double, double> > dist_ang;
       vector<double> parvals;
       RevEngUtils::distToSurf(adj_reg[ki]->pointsBegin(),
			       adj_reg[ki]->pointsEnd(), surf, tol,
			       maxd, avd, num_in, num2_in,
			       in, out, parvals, dist_ang, angtol);
       int num_ang_in = 0;
       double avang = 0.0;
       double nfac = 1.0/(double)num;
       for (size_t kh=0; kh<dist_ang.size(); ++kh)
	 {
	   avang += nfac*dist_ang[kh].second;
	   if (dist_ang[kh].second <= angtol)
	     num_ang_in++;
	 }
       
       int adj_sf_flag = adj_reg[ki]->defineSfFlag(0, tol, num_in, num2_in, 
						   avd, cyllike);

       double tol_fac = 2.0;
       double tol_fac2 = 5.0;
       double ang_lim = 0.9;
       double num2_lim = 0.25;
       if (adj_sf_flag == ACCURACY_OK ||
	   (adj_sf_flag <= surfflag_ &&
	    avd <= tol_fac*tol && (double)num_ang_in >= ang_lim*(double)num &&
	    (double)num2_in >= num2_lim*(double)num2_in))
	 {
	   // Accepted. Should possibly have a test on number of points and number
	   // of neighbours
	   changed = true;
	   
	   // Include adjacent region in present
	   vector<RevEngRegion*> added_adjacent;
	   includeAdjacentRegion(adj_reg[ki], maxd, avd, num_in, 
				 num2_in, parvals, dist_ang,
				 added_adjacent);
	   grown_regions.push_back(adj_reg[ki]);

	   int num_sf = adj_reg[ki]->numSurface();
	   for (int kb=0; kb<num_sf; ++kb)
	     adj_surfs.push_back(adj_reg[ki]->getSurface(kb));
	   if (adj_reg[ki]->numRevEdges() > 0)
	     {
	       vector<RevEngEdge*> rev_edgs = adj_reg[ki]->getAllRevEdges();
	       for (size_t kr=0; kr<rev_edgs.size(); ++kr)
		 {
		   RevEngRegion *adj1, *adj2;
		   rev_edgs[kr]->getAdjacent(adj1, adj2);
		   RevEngRegion *other = (adj1 == adj_reg[ki]) ? adj2 : adj1;
		   other->removeRevEngEdge(rev_edgs[kr]);
		 }
	       adj_edgs.insert(adj_edgs.end(), rev_edgs.begin(), rev_edgs.end());
	     }
	 }
       else if (avd <= tol_fac2*tol && (double)num_ang_in >= ang_lim*(double)num)
	 {
	   grow_cand curr_cand(adj_reg[ki], maxd, avd, avang, num_in, num2_in,
			       num_ang_in);
	   cands.push_back(curr_cand);
	 }
    }
       
  // Check if the surface should be updated
  if (changed)
    checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
	   
  // Integrate candidate neighbours if feasible
  if (cands.size() > 0)
    {
      integrateGrowCand(cands, mainaxis, tol, angtol, grown_regions,
			adj_surfs);
      int stop_break = 1;
    }

  double fac = 1.5;
  if ((int)group_points_.size() > (int)(fac*(double)num_group_points) ||
      adj_surfs.size() > 0)
    {
      for (size_t kj=0; kj<rev_edges_.size(); ++kj)
	rev_edges_[kj]->increaseExtendCount();
    }
}


//===========================================================================
int RevEngRegion::getGrowAccuracy(RevEngRegion *other, double tol,
				  double angtol, double& maxdist,
				  double& avdist, int& num_inside, 
				  int& num2_inside, double& maxdist2,
				  double& avdist2, int& num_in2,
				  int& num2_in2,vector<double>& parvals,
				  vector<pair<double,double> >& distang)
//===========================================================================
{
  if (!hasSurface())
    return NOT_SET;

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  
  vector<RevEngPoint*> in, out;
  RevEngUtils::distToSurf(other->pointsBegin(),
			  other->pointsEnd(), surf, tol,
			  maxdist2, avdist2, num_in2, num2_in2,
			  in, out, parvals, distang, angtol);

  maxdist = std::max(maxdist_, maxdist2);
  int all_pts = (int)group_points_.size() + other->numPoints();
  double fac1 = (double)(group_points_.size())/(double)all_pts;
  double fac2 = (double)(other->numPoints())/(double)all_pts;
  avdist = fac1*avdist_ + fac2*avdist2;
  num_inside = num_inside_ + num_in2;
  num2_inside = num_inside2_ + num2_in2;

  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  int sf_flag = other->defineSfFlag(0, tol, num_in2, num2_in2, avdist2,
				    cyllike);
  return sf_flag;
 }

//===========================================================================
void RevEngRegion::getDomainBoundaries(double tol, double angtol,
				       vector<pair<int, double> >& bd_par1,
				       vector<pair<int, double> >& bd_par2)
//===========================================================================
{
  if (!hasSurface())
    return;

  // Check domain
  double angtol2 = 1.1*angtol;
  double dom[4];
  int ka;
  int num_pt = (int)group_points_.size();
  double dist, ang;
  for (ka=0; ka<num_pt; ++ka)
    {
      Vector2D uv = group_points_[ka]->getPar();
      group_points_[ka]->getSurfaceDist(dist, ang);
      if (ang <= angtol2)
	{
	  dom[0] = dom[1] = uv[0];
	  dom[2] = dom[3] = uv[1];
	  break;
	}
    }
  for (; ka<num_pt; ++ka)
    {
      Vector2D uv = group_points_[ka]->getPar();
      group_points_[ka]->getSurfaceDist(dist, ang);
      if (ang > angtol2)
	continue;
      dom[0] = std::min(dom[0], uv[0]);
      dom[1] = std::max(dom[1], uv[0]);
      dom[2] = std::min(dom[2], uv[1]);
      dom[3] = std::max(dom[3], uv[1]);
    }

  
  // boundary index: 0=umin, 1=umax, 2=vmin, 3=vmax
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(trim_edgs_[ki]->geomCurve());
      int dir;
      double parval;
      if (sfcv->isConstantCurve(tol, dir, parval))
	{
	  bool same;
	  int bd = sfcv->whichBoundary(tol, same);
	  if (bd >= 0)
	    bd_par1.push_back(std::make_pair(bd, parval));
	}
    }

  for (size_t ki=0; ki<bd_par1.size(); ki++)
    for (size_t kj=1; kj<bd_par1.size(); ++kj)
      if (bd_par1[kj].first < bd_par1[ki].first)
	std::swap(bd_par1[ki], bd_par1[kj]);

  int ix=0;
  size_t kj;
  for (size_t ki=0; ki<bd_par1.size(); ki=kj)
    {
      for (kj=ki+1; kj<bd_par1.size(); ++kj)
	if (bd_par1[ki].first != bd_par1[kj].first)
	  break;
      if (kj > ki+1)
	{
	  // Check consistence
	  double tmin, tmax;
	  tmin = tmax = bd_par1[ki].second;
	  for (size_t kr=ki+1; kr<kj; ++kr)
	    {
	      tmin = std::min(tmin, bd_par1[kr].second);
	      tmax = std::min(tmax, bd_par1[kr].second);
	    }
	  if (tmax - tmin > tol)
	    {
	      bd_par1.erase(bd_par1.begin()+ki, bd_par1.begin()+kj);
	      kj = ki;
	    }
	  else
	    {
	      bd_par1[ki] = std::make_pair(bd_par1[ki].first, 0.5*(tmin+tmax));
	      bd_par1.erase(bd_par1.begin()+ki+1, bd_par1.begin()+kj);
	      kj = ki+1;
	    }
	}
      if (bd_par1[ki].first != ix)
	bd_par2.push_back(std::make_pair(ix, dom[ix]));
      ++ix;
    }

  for (; ix<4; ++ix)
    bd_par2.push_back(std::make_pair(ix, dom[ix]));
     
}

//===========================================================================
void RevEngRegion::growBlendSurf(vector<RevEngRegion*>& next_blend, double tol,
				 double angtol, vector<RevEngRegion*>& grown_regions,
				 vector<vector<RevEngPoint*> >& added_regions)
//===========================================================================
{
  if (!hasSurface())
    return;

  if (!hasBlendEdge())
    return;
  RevEngRegion *adj1, *adj2;
  blend_edge_->getAdjacent(adj1, adj2);
  
  // Set absolute and vague domain boundaries
  vector<pair<int, double> > bd_par1;
  vector<pair<int, double> > bd_par2;
  getDomainBoundaries(tol, angtol, bd_par1, bd_par2);
  if (bd_par1.size() == 0)
    return;  // Not as expected

  vector<RevEngRegion*> adj_reg;
  adj_reg.insert(adj_reg.end(), adjacent_regions_.begin(), adjacent_regions_.end());
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    {
      if (adj_reg[ki]->hasSurface())
	continue;
#ifdef DEBUG_BLEND
      std::ofstream of("curr_blend_adj.g2");
      adj_reg[ki]->writeRegionPoints(of);
#endif
      // Compute distance to current surface
      shared_ptr<ParamSurface> surf = getSurface(0)->surface();
      double maxd, avd;
      int num_in, num2_in;
      vector<RevEngPoint*> in, out;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      RevEngUtils::distToSurf(adj_reg[ki]->pointsBegin(),
			      adj_reg[ki]->pointsEnd(), surf, tol,
			      maxd, avd, num_in, num2_in,
			      in, out, parvals, dist_ang, angtol);

      vector<RevEngRegion*> other_adj;
      adj_reg[ki]->getAdjacentRegions(other_adj);
      if (other_adj.size() == 1 && other_adj[0] == this)
	{
	  // Include adjacent region in current
	  vector<RevEngRegion*> added_adjacent;
	  includeAdjacentRegion(adj_reg[ki], maxd, avd, num_in, 
				num2_in, parvals, dist_ang,
				added_adjacent);
	  adj_reg[ki]->removeFromAdjacent();
	  adj_reg[ki]->clearRegionAdjacency();
	  grown_regions.push_back(adj_reg[ki]);

	}
      else
	{
	  // Distribute points according to identified boundaries
	  vector<vector<int> > out1_ix(bd_par1.size());
	  vector<int> in1_ix;
	  double eps = 1.0e-9;
	  int num = adj_reg[ki]->numPoints();
	  vector<RevEngPoint*> adj_pts = adj_reg[ki]->getPoints();
	  for (int ka=0; ka<num; ++ka)
	    {
	      size_t kj=0;
	      for (kj=0; kj<bd_par1.size(); ++kj)
		{
		  int ix = (bd_par1[kj].first <= 1) ? 0 : 1;
		  int sgn = (bd_par1[kj].first%2 == 0) ? -1 : 1;
		  if ((sgn < 0 && parvals[2*ka+ix] < bd_par1[kj].second+eps) ||
		      (sgn > 0 && parvals[2*ka+ix] > bd_par1[kj].second-eps))
		    {
		      out1_ix[kj].push_back(ka);
		      break;
		    }
		}
	      if (kj == bd_par1.size())
		in1_ix.push_back(ka);
	    }

	  if (other_adj.size() == 2)
	    {
	      RevEngRegion* other_reg = (other_adj[0] == this) ?
		other_adj[1] : other_adj[0];
	      if (other_reg == adj1 || other_reg == adj2)
		{
		  if ((int)in1_ix.size() == num)
		    {
		      vector<RevEngRegion*> added_adjacent;
		      includeAdjacentRegion(adj_reg[ki], maxd, avd, num_in, 
					    num2_in, parvals, dist_ang,
					    added_adjacent);
		    }
		  else
		    {
		      for (size_t kj=0; kj<in1_ix.size(); ++kj)
			{
			  RevEngPoint *curr = adj_pts[in1_ix[kj]];
			  int ix = in1_ix[kj];
			  curr->setPar(Vector2D(parvals[2*ix],
						parvals[2*ix+1]));
			  curr->setSurfaceDist(dist_ang[ix].first, dist_ang[ix].second);
			  adj_reg[ki]->removePoint(curr);
			  addPoint(curr);
			}
		      shared_ptr<ParamSurface> surf2 = other_reg->getSurface(0)->surface();
		      double maxd2, avd2;
		      int num_in2, num2_in2;
		      vector<RevEngPoint*> in2, out2;
		      vector<pair<double, double> > dist_ang2;
		      vector<double> parvals2;
		      RevEngUtils::distToSurf(adj_reg[ki]->pointsBegin(),
					      adj_reg[ki]->pointsEnd(), surf2,
					      tol, maxd2, avd2, num_in2,
					      num2_in2, in2, out2, parvals2,
					      dist_ang2, angtol);
		      
		      vector<RevEngRegion*> added_adjacent;
		      other_reg->includeAdjacentRegion(adj_reg[ki], maxd2, avd2,
						       num_in2, num2_in2,
						       parvals2, dist_ang2,
						       added_adjacent);
		      
		    }
		  adj_reg[ki]->removeFromAdjacent();
		  adj_reg[ki]->clearRegionAdjacency();
		  grown_regions.push_back(adj_reg[ki]);
		  continue;
		}
	    }
	  
	  if (in1_ix.size() == 0)
	    continue;

 	  vector<vector<int> > out2_ix;
	  vector<int> in2_ix;
	  if (bd_par2.size() > 0)
	    {
	      out2_ix.resize(bd_par2.size());
	      for (size_t kr=0; kr<in1_ix.size(); ++kr)
		{
		  int ka = in1_ix[kr];
		  size_t kj=0;
		  for (kj=0; kj<bd_par2.size(); ++kj)
		    {
		      int ix = (bd_par2[kj].first <= 1) ? 0 : 1;
		      int sgn = (bd_par2[kj].first%2 == 0) ? -1 : 1;
		      if ((sgn < 0 && parvals[2*ka+ix] < bd_par2[kj].second+eps) ||
			  (sgn > 0 && parvals[2*ka+ix] > bd_par2[kj].second-eps))
			{
			  out2_ix[kj].push_back(ka);
			  break;
			}
		    }
		  if (kj == bd_par2.size())
		    in2_ix.push_back(ka);
		}
	    }
	  else
	    in2_ix = in1_ix;

	  // Identify grow candidates from adjacent
	  vector<RevEngRegion*> adj_grow;
	  vector<RevEngRegion*> remain_adj;
	  if (in2_ix.size() < in1_ix.size())
	    adj_grow.push_back(this);
	  for (size_t kj=0; kj<other_adj.size(); ++kj)
	    {
	      if (other_adj[kj] == adj1 || other_adj[kj] == adj2)
		adj_grow.push_back(other_adj[kj]);
	      else
		{
		  size_t kh;
		  for (kh=0; kh<next_blend.size(); ++kh)
		    if (next_blend[kh] == other_adj[kj])
		      break;
		  if (kh < next_blend.size())
		    adj_grow.push_back(other_adj[kj]);
		  else if (other_adj[kj] != this)
		    remain_adj.push_back(other_adj[kj]);
		}
	    }
	  int glast = (int)adj_grow.size()-1;
	  for (size_t kj=0; kj<adj_grow.size(); ++kj)
	    if (adj_grow[kj] == adj1)
	      {
		if ((int)kj < glast)
		  std::swap(adj_grow[kj], adj_grow[glast]);
		glast--;
	      }
	  for (size_t kj=0; kj<adj_grow.size(); ++kj)
	    if ((int)kj < glast && adj_grow[kj] == adj2)
	      std::swap(adj_grow[kj], adj_grow[glast]);

	  vector<vector<RevEngPoint*> > in_pts(adj_grow.size());
	  for (size_t kj=0; kj<adj_grow.size(); ++kj)
	    {
	      for (size_t kh=0; kh<out2_ix.size(); ++kh)
		if (out2_ix[kh].size() > 0)
		  adj_grow[kj]->blendGrowFromAdjacent(adj_reg[ki], out2_ix[kh],
						      tol, angtol, in_pts[kj]);
	      if (adj_grow[kj] != this)
		{
		  for (size_t kh=0; kh<out1_ix.size(); ++kh)
		    if (out1_ix[kh].size() > 0)
		      adj_grow[kj]->blendGrowFromAdjacent(adj_reg[ki], out1_ix[kh],
							  tol, angtol, in_pts[kj]);
		}
	    }

	  for (int ka=0; ka<num; ++ka)
	    {
	      adj_reg[ki]->getPoint(ka)->unsetVisited();
	    }

	  for (size_t kh=0; kh<in2_ix.size(); ++kh)
	    {
	      int ix = in2_ix[kh];
	      adj_pts[ix]->setPar(Vector2D(parvals[2*ix], parvals[2*ix+1]));
	      adj_pts[ix]->setSurfaceDist(dist_ang[ix].first, dist_ang[ix].second);
	      adj_reg[ki]->removePoint(adj_pts[ix]);
	      addPoint(adj_pts[ix]);
	    }
	      
	  for (size_t kj=0; kj<adj_grow.size(); ++kj)
	    {
	      if (in_pts[kj].size() > 0)
		{
		  shared_ptr<ParamSurface> surf2 = adj_grow[kj]->getSurface(0)->surface();
		  double maxd2, avd2;
		  int num_in2, num2_in2;
		  vector<RevEngPoint*> in2, out2;
		  vector<pair<double, double> > dist_ang2;
		  vector<double> parvals2;
		  RevEngUtils::distToSurf(in_pts[kj].begin(), in_pts[kj].end(),
					  surf2, tol, maxd2, avd2, num_in2,
					  num2_in2, in2, out2, parvals2,
					  dist_ang2, angtol);
		  for (size_t kh=0; kh<in_pts[kj].size(); ++kh)
		    {
		      in_pts[kj][kh]->setPar(Vector2D(parvals[2*kh],
		  				  parvals[2*kh+1]));
		      in_pts[kj][kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
		      adj_reg[ki]->removePoint(in_pts[kj][kh]);
		      adj_grow[kj]->addPoint(in_pts[kj][kh]);
		    }
		  
		}
	      adj_grow[kj]->computeDomain();
	    }

	  if (adj_reg[ki]->numPoints() == 0)
	    {
	      adj_reg[ki]->removeFromAdjacent();
	      adj_reg[ki]->clearRegionAdjacency();
	      grown_regions.push_back(adj_reg[ki]);
	    }
	  else
	    {
	      vector<vector<RevEngPoint*> > sep_groups;
	      adj_reg[ki]->splitRegion(sep_groups);
	      if (sep_groups.size() > 0)
	  	{
	  	  added_regions.insert(added_regions.end(), sep_groups.begin(),
	  			       sep_groups.end());
	  	  adj_reg[ki]->updateRegionAdjacency();
	  	}
	    }

	  if (remain_adj.size() > 0)
	    {
	      for (size_t kj=0; kj<remain_adj.size(); ++kj)
		{
		  size_t kh;
		  for (kh=0; kh<adj_reg.size(); ++kh)
		    if (adj_reg[kh] == remain_adj[kj])
		      break;
		  if (kh == adj_reg.size())
		    adj_reg.push_back(remain_adj[kj]);
		}
#ifdef DEBUG_BLEND
	      std::ofstream of2("updated_blend_grow.g2");
	      adj_reg[ki]->writeRegionPoints(of2);
	      for (size_t kr=0; kr<adj_grow.size(); ++kr)
		adj_grow[kr]->writeRegionPoints(of2);
#endif
	    }
	  int stop_break = 1;
	}
    }
#ifdef DEBUG_BLEND
  std::ofstream of3("updated_blend_grow2.g2");
  writeRegionPoints(of3);
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    adj_reg[ki]->writeRegionPoints(of3);
#endif
  int stop_break2 = 1;
}

  
//===========================================================================
void RevEngRegion::blendGrowFromAdjacent(RevEngRegion* adjacent,
					 vector<int>& pt_ix, double tol,
					 double angtol,
					 vector<RevEngPoint*>& grow_pt)
//===========================================================================
{
  if (!hasSurface())
    return;

  vector<pair<int, double> > bd_par1;
  vector<pair<int, double> > bd_par2;
  if (hasBlendEdge())
    getDomainBoundaries(tol, angtol, bd_par1, bd_par2);

  double eps = 1.0e-6;
  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  vector<RevEngPoint*> adj_pts = adjacent->getPoints();
  double fac = 1.1;
  double tol2 = std::min(fac*tol, avdist_); //fac*std::max(tol, avdist_);
  vector<RevEngPoint*> next_pts;
  double par[2];
  double dist;
  Point close;
  Point norm, norm2, norm3;
  vector<int> remove_ix;
  for (int ka=(int)pt_ix.size()-1; ka>=0; --ka)
    {
      RevEngPoint* curr = adj_pts[pt_ix[ka]];
      curr->setMarkIx(ka);
      if (curr->isNeighbour(this))
	{
	  curr->setVisited();
	  Vector3D xyz = curr->getPoint();
	  Point pnt(xyz[0], xyz[1], xyz[2]);
	  surf->closestPoint(pnt, par[0], par[1], close, dist, eps);
	  
	  size_t kh=0;
	  for (kh=0; kh<bd_par1.size(); ++kh)
	    {
	      int ix = (bd_par1[kh].first <= 1) ? 0 : 1;
	      int sgn = (bd_par1[kh].first%2 == 0) ? -1 : 1;
	      if ((sgn < 0 && par[ix] < bd_par1[kh].second+eps) ||
		  (sgn > 0 && par[ix] > bd_par1[kh].second-eps))
		break;
	    }

	  size_t kr=0;
	  for (kr=0; kr<bd_par2.size(); ++kr)
	    {
	      int ix = (bd_par2[kr].first <= 1) ? 0 : 1;
	      int sgn = (bd_par2[kr].first%2 == 0) ? -1 : 1;
	      if ((sgn < 0 && par[ix] < bd_par2[kr].second+eps) ||
		  (sgn > 0 && par[ix] > bd_par2[kr].second-eps))
		break;
	    }

	      
	  if (kh == bd_par1.size() &&
	      (dist <= tol2 || (bd_par2.size() > 0 && kr == bd_par2.size())))
	    {
	      surf->normal(norm, par[0], par[1]);
	      norm2 = curr->getLocFuncNormal();
	      norm3 = curr->getTriangNormal();
	      double ang = norm.angle(norm2);
	      double ang2 = norm.angle(norm3);
	      ang = std::min(std::min(ang,M_PI-ang), std::min(ang2,M_PI-ang2));
	      if (dist <= tol || ang < fac*angtol ||
		  (bd_par2.size() > 0 && kr == bd_par2.size() && ang <= angtol))
		{
		  next_pts.push_back(curr);
		  grow_pt.push_back(curr);
		  remove_ix.push_back(ka);
		}
	    }

	}
    }

  if (next_pts.size() == 0)
    {
      for (int ka=(int)pt_ix.size()-1; ka>=0; --ka)
	{
	  RevEngPoint* curr = adj_pts[pt_ix[ka]];
	  curr->unsetMarkIx();
	}
      return;
    }
  
#ifdef DEBUG_GROW
  std::ofstream of3("seed_blendgrow.g2");
  if (next_pts.size() > 0)
    {
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << next_pts.size() << std::endl;
      for (size_t kr=0; kr<next_pts.size(); ++kr)
	of3 << next_pts[kr]->getPoint() << std::endl;
    }
  std::ofstream of4("cand_blendgrow.g2");
  if (pt_ix.size() > 0)
    {
      of4 << "400 1 0 4  0 255 0 255" << std::endl;
      of4 << pt_ix.size() << std::endl;
      for (size_t kr=0; kr<pt_ix.size(); ++kr)
	of4 << adj_pts[pt_ix[kr]]->getPoint() << std::endl;
    }
#endif
	  
  // Grow
  for (size_t ki=0; ki<next_pts.size(); ++ki)
    {
      vector<ftSamplePoint*> next2 = next_pts[ki]->getNeighbours();
      for (size_t kj=0; kj<next2.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next2[kj]);
	  if (curr->visited())
	    continue;
	  RevEngRegion *adj_reg = curr->region();
	  if (adj_reg != adjacent)
	    continue;
	  if (curr->getMarkIx() < 0)
	    continue;
	  curr->setVisited();
	  Vector3D xyz = curr->getPoint();
	  surf->closestPoint(Point(xyz[0],xyz[1],xyz[2]), par[0], par[1], close, dist, eps);
	  size_t kh=0;
	  for (kh=0; kh<bd_par1.size(); ++kh)
	    {
	      int ix = (bd_par1[kh].first <= 1) ? 0 : 1;
	      int sgn = (bd_par1[kh].first%2 == 0) ? -1 : 1;
	      if ((sgn < 0 && par[ix] < bd_par1[kh].second+eps) ||
		  (sgn > 0 && par[ix] > bd_par1[kh].second-eps))
		break;
	    }

	  size_t kr=0;
	  for (kr=0; kr<bd_par2.size(); ++kr)
	    {
	      int ix = (bd_par2[kr].first <= 1) ? 0 : 1;
	      int sgn = (bd_par2[kr].first%2 == 0) ? -1 : 1;
	      if ((sgn < 0 && par[ix] < bd_par2[kr].second+eps) ||
		  (sgn > 0 && par[ix] > bd_par2[kr].second-eps))
		break;
	    }

	  if (kh == bd_par1.size() &&
	      (dist <= tol2 || (bd_par2.size() > 0 && kr == bd_par2.size())))
	    {
	      surf->normal(norm, par[0], par[1]);
	      norm2 = curr->getLocFuncNormal();
	      norm3 = curr->getTriangNormal();
	      double ang = norm.angle(norm2);
	      double ang2 = norm.angle(norm3);
	      ang = std::min(std::min(ang,M_PI-ang), std::min(ang2,M_PI-ang2));
	      if (dist <= tol || ang < fac*angtol ||
		  (bd_par2.size() > 0 && kr == bd_par2.size() && ang <= angtol))
		{
		  grow_pt.push_back(curr);
		  next_pts.push_back(curr);
		  remove_ix.push_back(curr->getMarkIx());
		}
	    }
	}
    }

  for (int ka=(int)pt_ix.size()-1; ka>=0; --ka)
    {
      RevEngPoint* curr = adj_pts[pt_ix[ka]];
      curr->unsetMarkIx();
    }

  if (remove_ix.size() > 1)
    std::sort(remove_ix.begin(), remove_ix.end());
  for (int ka=(int)remove_ix.size()-1; ka>=0; --ka)
    pt_ix.erase(pt_ix.begin()+remove_ix[ka]);
  // for (size_t kj=0; kj<remove_ix.size(); ++kj)
  //   {
  //     auto it = std::find(pt_ix.begin(), pt_ix.end(), remove_ix[kj]);
  //     if (it != pt_ix.end())
  // 	{
  // 	  std::swap(*it, pt_ix[pt_ix.size()-1]);
  // 	  pt_ix.pop_back();
  // 	}
  //   }
}


//===========================================================================
void RevEngRegion::integrateGrowCand(vector<grow_cand>& cand,
				     Point mainaxis[3], double tol,
				     double angtol, vector<RevEngRegion*>& grown_regions,
				     vector<HedgeSurface*>& adj_surfs)
//===========================================================================
{
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  ClassType classtype = surf->instanceType();
  bool cyllike;
  double tol3 = 3.0*tol;
  double anglim = 0.75;

#ifdef DEBUG_GROW
  std::ofstream of1("source_reg.g2");
  writeRegionPoints(of1);
#endif
  // Collect all points
  // Better to keep them separate, but don't change all that code now
  vector<RevEngPoint*> all_points;
  all_points.insert(all_points.end(), group_points_.begin(), group_points_.end());
  while (cand.size() > 0)
    {
      vector<int> all_num(cand.size()+1);
      int num;
      all_num[0] = num = (int)group_points_.size();
      for (size_t ki=0; ki<cand.size(); ++ki)
	{
	  all_points.insert(all_points.end(), cand[ki].cand_->pointsBegin(),
			    cand[ki].cand_->pointsEnd());
	  all_num[ki+1] = cand[ki].cand_->numPoints();
	  num += all_num[ki+1];
	}

      // Compute updated surface
      shared_ptr<ParamSurface> merged1, merged2;
      computeSurface(all_points, mainaxis, tol, angtol, classtype,
		     merged1, merged2, cyllike);

      // Check accuracy
      vector<double> maxd(cand.size()+1), avd(cand.size()+1), avang(cand.size()+1, 0.0);
      vector<int> num_in(cand.size()+1), num2_in(cand.size()+1), ang_in(cand.size()+1, 0);
      vector<vector<double> > parvals(cand.size()+1);
      vector<vector<pair<double,double> > > dist_ang(cand.size()+1);
      vector<RevEngPoint*> inpt, outpt;

      if (!merged1.get())
	{
	  cand.clear();
	  break;
	}
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(), merged1,
			      tol, maxd[0], avd[0], num_in[0], num2_in[0], inpt,
			      outpt, parvals[0], dist_ang[0], angtol);
      if (merged2.get())
	{
	  double maxd2, avd2;
	  int num_in2, num2_in2;
	  vector<double> parvals2;
	  vector<pair<double,double> > dist_ang2;
	  vector<RevEngPoint*> inpt2, outpt2;
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(), merged2,
				  tol, maxd2, avd2, num_in2, num2_in2, inpt,
				  outpt2, parvals2, dist_ang2, angtol);
	  if (avd2 < avd[0] && num_in2+num2_in2 > num_in[0]+num2_in[0])
	    {
	      std::swap(merged1, merged2);
	      std::swap(maxd[0], maxd2);
	      std::swap(avd[0], avd2);
	      std::swap(num_in[0], num_in2);
	      std::swap(num2_in[0], num2_in2);
	      std::swap(parvals[0], parvals2);
	      std::swap(dist_ang[0], dist_ang2);
	    }
	}
      double frac = 1.0/(double)group_points_.size();
      for (size_t kh=0; kh<dist_ang[0].size(); ++kh)
	  {
	    avang[0] += frac*dist_ang[0][kh].second;
	    if (dist_ang[0][kh].second <= angtol)
	      ang_in[0]++;
	  }

      double maxdist = maxd[0];
      double avdist = (double)all_num[0]*avd[0]/(double)num;
      int num_inside = num_in[0];
      int num2_inside = num2_in[0];
      double meanang = (double)all_num[0]*avang[0]/(double)num;
      int inside_ang = ang_in[0];
      
      for (size_t ki=0; ki<cand.size(); ++ki)
	{
	  inpt.clear();
	  outpt.clear();
	  RevEngUtils::distToSurf(cand[ki].cand_->pointsBegin(),
				  cand[ki].cand_->pointsEnd(), merged1,
				  tol, maxd[ki+1], avd[ki+1], num_in[ki+1], num2_in[ki+1], 
				  inpt, outpt, parvals[ki+1], dist_ang[ki+1], angtol);
	  frac = 1.0/(double)cand[ki].cand_->numPoints();
	  for (size_t kh=0; kh<dist_ang[ki+1].size(); ++kh)
	    {
	      avang[ki+1] += frac*dist_ang[ki+1][kh].second;
	      if (dist_ang[ki+1][kh].second <= angtol)
		ang_in[ki+1]++;
	    }

	  maxdist = std::max(maxdist, maxd[ki+1]);
	  avdist += (double)all_num[ki+1]*avd[ki+1]/(double)num;
	  num_inside += num_in[ki+1];
	  num2_inside += num2_in[ki+1];
	  meanang += (double)all_num[ki+1]*avang[ki+1]/(double)num;
	  inside_ang += ang_in[ki+1];
	}
      int sf_flag = defineSfFlag(num, 0, tol, num_inside, num2_inside, avdist,
				 cyllike);
#ifdef DEBUG_GROW
      std::ofstream of2("adj_source_reg.g2");
      for (size_t ki=0; ki<cand.size(); ++ki)
	{
	  cand[ki].cand_->writeRegionPoints(of2);
	}
#endif
      
      vector<size_t> out_cand;
      for (size_t ki=0; ki<cand.size(); ++ki)
	{
	  if (avd[ki+1] > tol3 ||
	      (avd[ki+1] > tol && (double)ang_in[ki+1] < anglim*(double)all_num[ki+1]))
	    out_cand.push_back(ki);
	}

      if (sf_flag >= ACCURACY_POOR && out_cand.size() == 0)
	{
	  double max_av = avd[1];
	  size_t ix = 0;
	  for (size_t ki=1; ki<cand.size(); ++ki)
	    if (avd[ki+1] > max_av)
	      {
		max_av = avd[ki+1];
		ix = ki;
	      }

	  out_cand.push_back(ix);
	}

      if (out_cand.size() == 0)
	{
	  // Integrate
	  for (size_t kh=0; kh<group_points_.size(); ++kh)
	    {
	      group_points_[kh]->setPar(Vector2D(parvals[0][2*kh],parvals[0][2*kh+1]));
	      group_points_[kh]->setSurfaceDist(dist_ang[0][kh].first,
						dist_ang[0][kh].second);
	    }
	  for (size_t ki=0; ki<cand.size(); ++ki)
	    {
	      vector<RevEngRegion*> added_adjacent;
	      includeAdjacentRegion(cand[ki].cand_, maxd[ki+1], avd[ki+1],
				    num_in[ki+1], num2_in[ki+1],
				    parvals[ki+1], dist_ang[ki+1], added_adjacent);
	      grown_regions.push_back(cand[ki].cand_);
	      int num_sf = cand[ki].cand_->numSurface();
	      for (int kb=0; kb<num_sf; ++kb)
		adj_surfs.push_back(cand[ki].cand_->getSurface(kb));
	    }
	  associated_sf_[0]->replaceSurf(merged1);
	  if (!merged1->isBounded())
	    {
	      double diag = bbox_.low().dist(bbox_.high());
	      associated_sf_[0]->limitSurf(2*diag);
	    }
	  updateInfo(tol, angtol);
	  setSurfaceFlag(sf_flag);
	  for (size_t kj=0; kj<rev_edges_.size(); ++kj)
	    rev_edges_[kj]->replaceSurf(this, merged1, tol);

	  break;
	}
      else
	{
	  for (int ka=(int)out_cand.size()-1; ka>=0; --ka)
	    cand.erase(cand.begin()+ka);
	}
      int stop_break0 = 1;
    }
  
  int stop_break = 1;
}


//===========================================================================
bool RevEngRegion::mergePlanarReg(double zero_H, double zero_K, double tol,
				  Point mainaxis[3],
				  vector<RevEngRegion*>& grown_regions)
//===========================================================================
{
  if (!feasiblePlane(zero_H, zero_K))
    return false;
#ifdef DEBUG_MERGE
  std::ofstream of1("source_planar_reg.g2");
  writeRegionPoints(of1);
#endif
  
  double anglim = 0.1*M_PI; //0.05;
  vector<RevEngRegion*> merge_cand;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if (!(*it)->feasiblePlane(zero_H, zero_K))
	continue;
      if ((*it)->hasSurface())
	continue;
      if ((*it)->hasAssociatedBlend())
	continue;

      Point norm = (*it)->getMeanNormal();
      double ang = avnorm_.angle(norm);
      if (ang <= anglim)
	merge_cand.push_back(*it);
    }

  if (merge_cand.size() == 0)
    return false;
  
#ifdef DEBUG_MERGE
  std::ofstream of2("adj_planar_reg.g2");
  for (size_t ki=0; ki<merge_cand.size(); ++ki)
    merge_cand[ki]->writeRegionPoints(of2);
#endif

  // Check accuracy
  vector<RevEngPoint*> all_pts;
  all_pts.insert(all_pts.end(), group_points_.begin(), group_points_.end());
  for (size_t ki=0; ki<merge_cand.size(); ++ki)
    all_pts.insert(all_pts.end(), merge_cand[ki]->pointsBegin(),
		   merge_cand[ki]->pointsEnd());
  shared_ptr<Plane> plane = computePlane(all_pts, avnorm_, mainaxis);

  double angtol = -1.0;
  int min_pt = 2;  // Not crucial here
  double maxdist, avdist;
  int num_in, num2_in;
  vector<pair<double, double> > dist_ang;
  vector<double> parvals;
  vector<RevEngPoint*> inpt, outpt;
  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			  plane, tol, maxdist, avdist, num_in, num2_in, inpt,
			  outpt, parvals, dist_ang, angtol);
  if (!accuracyOK(min_pt, tol, num_in, avdist))
    return false;

  for (int kj=(int)merge_cand.size()-1; kj>=0; --kj)
    {
      double maxdist2, avdist2;
      int num_in2, num2_in2;
      vector<pair<double, double> > dist_ang2;
      vector<double> parvals2;
      vector<RevEngPoint*> inpt2, outpt2;
      RevEngUtils::distToSurf(merge_cand[kj]->pointsBegin(),
			      merge_cand[kj]->pointsEnd(),
			      plane, tol, maxdist2, avdist2, num_in2, num2_in2, inpt2,
			      outpt2, parvals2, dist_ang2, angtol);
      if (!merge_cand[kj]->accuracyOK(0, tol, num_in2, avdist2))
	merge_cand.erase(merge_cand.begin()+kj);	  
    }
  if (merge_cand.size() == 0)
    return false;

  // Include next layer of adjacent regions
  vector<RevEngRegion*> adj_reg;
  size_t cand_size = merge_cand.size();
  Point dir = plane->direction();
  for (size_t ki=0; ki<cand_size; ++ki)
    {
      for (auto it=merge_cand[ki]->adjacent_regions_.begin();
  	     it!=merge_cand[ki]->adjacent_regions_.end(); ++it)
  	{
  	  if ((*it) == this)
  	    continue;
	  if (!(*it)->feasiblePlane(zero_H, zero_K))
	    continue;
	  if ((*it)->hasSurface())
	    continue;
	  if ((*it)->hasAssociatedBlend())
	    continue;
      
  	  size_t kr;
  	  for (kr=0; kr<merge_cand.size(); ++kr)
  	    if ((*it) == merge_cand[kr])
  	      break;
	  if (kr < merge_cand.size())
	    continue;

  	  Point norm = (*it)->getMeanNormal();
  	  double ang = dir.angle(norm);
  	  if (ang > anglim)
  	    continue;
  	  double maxdist3, avdist3;
  	  int num_in3, num2_in3;
  	  vector<pair<double, double> > dist_ang3;
  	  vector<double> parvals3;
  	  vector<RevEngPoint*> inpt3, outpt3;
  	  RevEngUtils::distToSurf((*it)->pointsBegin(),
  				  (*it)->pointsEnd(),
  				  plane, tol, maxdist3, avdist3, num_in3, num2_in3, inpt3,
  				  outpt3, parvals3, dist_ang3, angtol);
  	  if ((*it)->accuracyOK(0, tol, num_in3, avdist3))
  	    merge_cand.push_back(*it);
  	}
    }

#ifdef DEBUG_MERGE
  std::ofstream of3("adj_planar2_reg.g2");
  for (size_t ki=0; ki<merge_cand.size(); ++ki)
    merge_cand[ki]->writeRegionPoints(of3);
#endif

  // Integrate identified regions
  for (size_t ki=0; ki<merge_cand.size(); ++ki)
    {
      grown_regions.push_back(merge_cand[ki]);
      for (auto it=merge_cand[ki]->adjacent_regions_.begin();
	       it!=merge_cand[ki]->adjacent_regions_.end(); ++it)
	{
	  if (*it != this)
	    {
	      addAdjacentRegion(*it);
	      (*it)->addAdjacentRegion(this);
	      (*it)->removeAdjacentRegion(merge_cand[ki]);
	    }
	}
      // for (size_t kj=ki+1; kj<merge_cand.size(); ++kj)
      // 	if (merge_cand[kj]->isAdjacent(merge_cand[ki]))
      // 	  merge_cand[kj]->removeAdjacentRegion(merge_cand[ki]);
      std::vector<std::pair<double, double> > dummy;
      addRegion(merge_cand[ki], dummy);
      removeAdjacentRegion(merge_cand[ki]);
    }

  return true;
}

//===========================================================================
void RevEngRegion::mergeAdjacentSimilar(double tol, double angtol,
					vector<RevEngRegion*>& grown_regions,
					vector<HedgeSurface*>& adj_surfs,
					vector<RevEngEdge*>& adj_edgs)
//===========================================================================
{
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  if (hasAssociatedBlend())
    return;

#ifdef DEBUG_MERGE
  std::ofstream of1("source_group.g2");
  writeRegionPoints(of1);
#endif
  int sfcode;
  ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  bool cyllike = (classtype == Class_Cylinder || classtype == Class_Cone);
  //HedgeSurface *hedge = associated_sf_[0];
  vector<RevEngRegion*> adj_reg;
  vector<double> score;
  double frac2 = 0.75;
  double frac3 = 2.0;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if ((*it)->hasSurface() && (!(*it)->hasAssociatedBlend()))
	{
	  double curr_score;
	  bool compatible = associated_sf_[0]->isCompatible((*it)->getSurface(0),
							    angtol, tol,
							    classtype, curr_score);
	  if (compatible)
	    {
	      adj_reg.push_back(*it);
	      score.push_back(curr_score);
	    }
	}
    }

  if (adj_reg.size() == 0)
    return; // Nothing with which to merge

#ifdef DEBUG_MERGE
  std::ofstream of2("adj_group.g2");
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    adj_reg[ki]->writeRegionPoints(of2);
#endif
  for (size_t ki=0; ki<adj_reg.size(); ++ki)
    for (size_t kj=ki+1; kj<adj_reg.size(); ++kj)
      if (score[kj] > score[ki])
	{
	  std::swap(score[ki], score[kj]);
	  std::swap(adj_reg[ki], adj_reg[kj]);
	}

  shared_ptr<ParamSurface> surf;
  int ka;
  double maxdist=0.0, avdist=0.0;
  int num_in = 0;
  vector<double> maxd(adj_reg.size()+1, 0.0);
  vector<double> avd(adj_reg.size()+1, 0.0);
  vector<int> ninside(adj_reg.size()+1, 0);
  vector<int> ninside2(adj_reg.size()+1, 0);
  vector<vector<double> > parvals(adj_reg.size()+1);
  vector<vector<pair<double,double> > > dist_ang(adj_reg.size()+1);
  for (ka=(int)adj_reg.size(); ka>=1; --ka)
    {
      // Create surface from combined point set
      vector<pair<vector<RevEngPoint*>::iterator,
		  vector<RevEngPoint*>::iterator> > points;
      BoundingBox bbox(3);
      vector<int> nmbpts;

      nmbpts.push_back(numPoints());
      points.push_back(std::make_pair(pointsBegin(), pointsEnd()));
      bbox = boundingBox();
      
      for (int kb=0; kb<ka; ++kb)
	{
	  nmbpts.push_back(adj_reg[kb]->numPoints());
	  points.push_back(std::make_pair(adj_reg[kb]->pointsBegin(),
					  adj_reg[kb]->pointsEnd()));
	  bbox.addUnionWith(adj_reg[kb]->boundingBox());
	}
      
      int num_all = 0;
      for (int kb=0; kb<=ka; ++kb)
	num_all += nmbpts[kb];
      double frac = 1.0/(double)num_all;
      if (classtype == Class_Plane)
	{
	  surf = RevEngUtils::doMergePlanes(points, bbox, nmbpts);
	}
      else if (classtype == Class_Cylinder)
	{
	  surf = RevEngUtils::doMergeCylinders(points, bbox, nmbpts);
	}
      else if (classtype == Class_Sphere)
	{
	  Point normal = frac*numPoints()*avnorm_;
	  
	  for (int kb=1; kb<=ka; ++kb)
	    normal += frac*adj_reg[kb-1]->numPoints()*adj_reg[kb-1]->getMeanNormal();
	  normal.normalize_checked();
	  surf = RevEngUtils::doMergeSpheres(points, bbox, nmbpts, normal);
	}
      else if (classtype == Class_Torus)
	{
	  surf = RevEngUtils::doMergeTorus(points, bbox, nmbpts);
	}
      if (!surf.get())
	continue;

      // Check accuracy
#ifdef DEBUG_GROW
       std::ofstream of("in_out_adj.g2");
#endif
       maxdist = avdist = 0.0;
       num_in = 0;
       int num2_in = 0;
       vector<int> sfflag(ka+1, NOT_SET);
       bool flagOK = true;
      for (int kb=0; kb<=ka; ++kb)
	{
	  maxd[kb] = avd[kb] = 0.0;
	  ninside[kb] = ninside2[kb] = 0;
	  parvals[kb].clear();
	  dist_ang[kb].clear();
	  vector<RevEngPoint*> in, out;
	  RevEngUtils::distToSurf(points[kb].first, points[kb].second, surf, tol,
				  maxd[kb], avd[kb], ninside[kb], ninside2[kb],
				  in, out, parvals[kb],
				  dist_ang[kb], angtol);
	  maxdist = std::max(maxdist, maxd[kb]);
	  avdist += frac*nmbpts[kb]*avd[kb];
	  num_in += ninside[kb];
	  num2_in += ninside2[kb];
	  RevEngRegion *curr = (kb==0) ? this : adj_reg[kb-1];
	  sfflag[kb] = curr->defineSfFlag(0, tol, ninside[kb], ninside2[kb],
					  avd[kb], cyllike);
	  if (sfflag[kb] == NOT_SET)
	    flagOK = false;

#ifdef DEBUG_GROW
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << in.size() << std::endl;
	  for (size_t kj=0; kj<in.size(); ++kj)
	    of << in[kj]->getPoint() << std::endl;
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << out.size() << std::endl;
	  for (size_t kj=0; kj<out.size(); ++kj)
	    of << out[kj]->getPoint() << std::endl;
#endif
	}

      int num = numPoints();
      double init_max = frac*num*maxdist_;
      double init_av = frac*num*avdist_;
      int init_in = num_inside_;
      int init_in2 = num_inside2_;
      for (int kb=0; kb<ka; ++kb)
	{
	  int num2 = adj_reg[kb]->numPoints();
	  double max2, av2;
	  int num_in2, num2_in2;
	  adj_reg[kb]->getAccuracy(max2, av2, num_in2, num2_in2);
	  init_max += frac*num2*max2;
	  init_av += frac*num2*av2;
	  init_in += num_in2;
	  init_in2 += num2_in2;
	}
      if (flagOK && num2_in > num_all/2 && avdist < tol &&
	  (double)num2_in > frac2*(double)init_in2 /*&&
						     (avdist < frac3*init_av || avdist < frac2*tol)*/)
	break;

      // Swap adjacent regions to skip the least accurate region
      if (ka > 1 && avd[1] > avd[ka])  // The test should be made more accurate
	{
	  std::swap(adj_reg[0], adj_reg[ka-1]);
	  std::swap(score[0], score[ka-1]);
	}
    }

  if (!surf.get())
    return;
  if (ka >= 1)
    {
      setAccuracy(maxd[0], avd[0], ninside[0], ninside2[0]);
      for (size_t ki=0; ki<group_points_.size(); ++ki)
	{
	  group_points_[ki]->setPar(Vector2D(parvals[0][2*ki],parvals[0][2*ki+1]));
	  group_points_[ki]->setSurfaceDist(dist_ang[0][ki].first, dist_ang[0][ki].second);
	}
      
      for (int kb=0; kb<ka; ++kb)
	{
	  vector<RevEngRegion*> added_adjacent;
	  includeAdjacentRegion(adj_reg[kb], maxd[kb+1], avd[kb+1], ninside[kb+1],
				ninside2[kb+1], parvals[kb+1],
			dist_ang[kb+1], added_adjacent);
	  grown_regions.push_back(adj_reg[kb]);
	  int num_sf = adj_reg[kb]->numSurface();
	  for (int kc=0; kc<num_sf; ++kc)
	    adj_surfs.push_back(adj_reg[kb]->getSurface(kc));
	  vector<RevEngEdge*> rev_edgs = adj_reg[kb]->getAllRevEdges();
	  for (size_t kr=0; kr<rev_edgs.size(); ++kr)
	    {
	      RevEngRegion *adj1, *adj2;
	      rev_edgs[kr]->getAdjacent(adj1, adj2);
	      RevEngRegion *other = (adj1 == adj_reg[kb]) ? adj2 : adj1;
	      other->removeRevEngEdge(rev_edgs[kr]);
	    }
	  adj_edgs.insert(adj_edgs.end(), rev_edgs.begin(), rev_edgs.end());
	  removeAdjacentRegion(adj_reg[kb]);
	}
      updateInfo(tol, angtol);
      double maxdist2, avdist2;
      int num_in2, num2_in2;
      getAccuracy(maxdist2, avdist2, num_in2, num2_in2);
      int surfflag = defineSfFlag(0, tol, num_in2, num2_in2, avdist2, cyllike);
      setSurfaceFlag(surfflag);
      associated_sf_[0]->replaceSurf(surf);
      computeDomain();
    }
  
}

//===========================================================================
void
RevEngRegion::includeAdjacentRegion(RevEngRegion* reg, double maxd, double avd,
				    int num_inside, int num_inside2,
				    vector<double>& parvals,
				    vector<pair<double, double> >& dist_ang,
				    vector<RevEngRegion*>& added_adjacent)
//===========================================================================
{
  // First update parameter values
  size_t kr=0;
  for (auto it1=reg->pointsBegin(); it1!=reg->pointsEnd(); ++it1, kr+=2)
    {
      (*it1)->addMove();
      (*it1)->setPar(Vector2D(parvals[kr],parvals[kr+1]));
    }
  addRegion(reg, dist_ang, maxd, avd, num_inside, num_inside2);


  // Update adjacent regions
  for (auto it2=reg->adjacent_regions_.begin(); it2!=reg->adjacent_regions_.end(); ++it2)
    {
      if (*it2 != this)
	{
	  size_t nmb_adj = adjacent_regions_.size();
	  addAdjacentRegion(*it2);
	  if (adjacent_regions_.size() > nmb_adj)
	    {
	      added_adjacent.push_back(*it2);
	      (*it2)->addAdjacentRegion(this);
	    }
	  (*it2)->removeAdjacentRegion(reg);
	}
    }
}

//===========================================================================
void RevEngRegion::checkEdgeAssociation(double tol, int min_point_reg,
					vector<HedgeSurface*>& removed_sfs)
//===========================================================================
{
  vector<RevEngRegion*> adj_regs(adjacent_regions_.begin(), adjacent_regions_.end());
  std::set<RevEngEdge*> rev_edgs;
  for (size_t ki=0; ki<adj_regs.size(); ++ki)
    if (adj_regs[ki]->hasRevEdges())
      {
	vector<RevEngEdge*> curr_edgs = adj_regs[ki]->getAllRevEdges();
	rev_edgs.insert(curr_edgs.begin(), curr_edgs.end());
      }
  vector<RevEngEdge*> edges(rev_edgs.begin(), rev_edgs.end());
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      vector<RevEngRegion*> blend_regs;
      edges[ki]->getAllBlendRegs(blend_regs);
      size_t kj = 0;
      for (kj=0; kj<blend_regs.size(); ++kj)
	if (blend_regs[kj]->isAdjacent(this))
	  {
	    //std::cout << "Potential edge" << std::endl;
	    vector<shared_ptr<CurveOnSurface> > cvs;
	    edges[ki]->getCurve(cvs);
	    double width = edges[ki]->getDistance();
	    int num_in_blend = 0;
	    if (isInBlend(cvs, 2.0*width, tol, num_in_blend))
	      {
		edges[ki]->addBlendRegion(this);
		setAssociatedBlend(edges[ki]);
	      }
	    else if (num_in_blend > 0 && (int)group_points_.size() < min_point_reg)
	      {
		// Dubious surface
		if (associated_sf_.size() > 0)
		  removed_sfs.insert(removed_sfs.end(), associated_sf_.begin(),
				     associated_sf_.end());
		clearSurface();
	      }
	    break;
	  }
      if (kj < blend_regs.size())
	break;
    }
  int stop_break = 1;
}

//===========================================================================

void endPoints(vector<RevEngPoint*>& seq_pts, RevEngPoint*& end1,
		    RevEngPoint*& end2)
{
  size_t ix1 = 0, ix2 = 1;
  ix2 = std::min(ix2, seq_pts.size()-1);
  double maxlen = seq_pts[ix1]->pntDist(seq_pts[ix2]);
  for (size_t ki=0; ki<seq_pts.size(); ++ki)
    for (size_t kj=ki+1; kj<seq_pts.size(); ++kj)
      {
	double len = seq_pts[ki]->pntDist(seq_pts[kj]);
	if (len > maxlen)
	  {
	    maxlen = len;
	    ix1 = ki;
	    ix2 = kj;
	  }
      }
  end1 = seq_pts[ix1];
  end2 = seq_pts[ix2];
}

void getBranches(vector<vector<RevEngPoint*> >& bd_pts,
		 vector<vector<pair<RevEngPoint*,int> > >& branches)
{
  RevEngPoint *dummy = 0;
  branches.resize(bd_pts.size());
  for (size_t ki=0; ki<bd_pts.size(); ++ki)
    {
      for (size_t kr=ki+1; kr<bd_pts.size(); ++kr)
	{
	  for (size_t kj=0; kj<bd_pts[ki].size(); ++kj)
	    {
	      for (size_t kh=0; kh<bd_pts[kr].size(); ++kh)
		{
		  if (bd_pts[ki][kj]->isNeighbour(bd_pts[kr][kh]))
		    {
		      branches[ki].push_back(std::make_pair(bd_pts[ki][kj],(int)kr));
		      branches[kr].push_back(std::make_pair(bd_pts[kr][kh],(int)ki));
		    }
		}
	    }
	}
      if (branches[ki].size() == 0)
	branches[ki].push_back(std::make_pair(dummy,-1));
    }
}

int getNextSeq(vector<int>& prev, int curr,
	       vector<vector<pair<RevEngPoint*,int> > >& branches,
	       vector<vector<RevEngPoint*> >& bd_pts,
	       RevEngPoint*& pnt)
{
  vector<pair<RevEngPoint*,int> > cand;
  for (size_t ki=0; ki<branches[curr].size(); ++ki)
    {
      size_t kj;
      for (kj=0; kj<prev.size(); ++kj)
	if (branches[curr][ki].second == prev[kj])
	  break;
      if (kj < prev.size())
	continue;
      for (kj=0; kj<cand.size(); ++kj)
	if (cand[kj].first == branches[curr][ki].first &&
	    cand[kj].second == branches[curr][ki].second)
	  break;
      if (kj == cand.size())
	cand.push_back(branches[curr][ki]);
    }

  // Sort candidates
  for (size_t ki=0; ki<cand.size(); ++ki)
    for (size_t kj=ki+1; kj<cand.size(); ++kj)
      if (cand[kj].second < cand[ki].second)
	std::swap(cand[kj], cand[ki]);

  // Select one candidate for each adjacent curve
  vector<double> dist;
  for (size_t ki=0; ki<cand.size(); )
    {
      size_t kj;
      for (kj=ki+1;
	   kj<cand.size() && cand[kj].second == cand[ki].second; ++kj);

      int kr = cand[ki].second;
      double mindist = std::numeric_limits<double>::max();
      int ix = -1;
      if (kr >= 0)
	{
	  for (size_t kh=ki; kh<kj; ++kh)
	    {
	      for (size_t kv=0; kv<branches[kr].size(); ++kv)
		{
		  if (branches[kr][kv].second != curr)
		    continue;
		  double dd = cand[kh].first->pntDist(branches[kr][kv].first);
		  if (dd < mindist)
		    {
		      mindist = dd;
		      ix = (int)kh;
		    }
		}
	    }
	}
      for (int ka=(int)kj-1; ka>=(int)ki; --ka)
	{
	  if (ka != ix)
	    cand.erase(cand.begin()+ka);
	}
      if (ix >= 0)
	{
	  dist.push_back(mindist);
	  ++ki;
	}
    }

  if (cand.size() == 0)
    return -1;
  
  // Select adjacent branch
  int ix = 0;
  double mindist = dist[0];
  double eps = 1.0e-9;
  for (size_t ki=1; ki<cand.size(); ++ki)
    {
      if (fabs(dist[ki] - mindist) < eps)
	{
	  int num1 = (int)bd_pts[ix].size();
	  int num2 = (int)bd_pts[ki].size();
	  if (num2 > num1)
	    {
	      ix = (int)ki;
	      mindist= dist[ki];
	    }
	}
      else if (dist[ki] < mindist)
	{
	  ix = (int)ki;
	  mindist= dist[ki];
	}
    }
  pnt = cand[ix].first;
  return cand[ix].second;
}

void splitAtSeam(shared_ptr<ParamSurface>& surf, shared_ptr<ParamCurve>& cv,
		 double tol, double angtol, double maxlim,
		 vector<shared_ptr<ParamCurve> >& subcvs)
{
  vector<Point> pos;
  vector<Point> norm;
  vector<Point> axis;
  vector<int> type;  // 1=cylinder, 2=cone, 3=torus big, 4=torus small, 5=sphere
  shared_ptr<ElementarySurface> elemsf =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  double radius = elemsf->radius(0,0);
  RectDomain dom = surf->containingDomain();
  if (surf->instanceType() == Class_Cylinder ||
      surf->instanceType() == Class_Cone ||
      surf->instanceType() == Class_Sphere)
    {
      vector<Point> seam_pt(3);
      surf->point(seam_pt, dom.umin(), 0.5*(dom.vmin()+dom.vmax()), 1);
      pos.push_back(seam_pt[0]);
      seam_pt[1].normalize();
      norm.push_back(seam_pt[1]);
      axis.push_back(elemsf->direction());
      type.push_back((surf->instanceType() == Class_Cylinder) ? 1 :
		     ((surf->instanceType() == Class_Cone) ? 2 : 3));
    }
  else if (surf->instanceType() == Class_Torus)
    {
      shared_ptr<Torus> tor = dynamic_pointer_cast<Torus,ParamSurface>(surf);
      Point centre = tor->location();
      Point dir = tor->direction();
      vector<Point> seam_pt(3);
      surf->point(seam_pt, dom.umin(), dom.vmin(), 1);
      seam_pt[1].normalize();

      pos.push_back(centre);
      norm.push_back(dir);
      axis.push_back(seam_pt[1]);
      type.push_back(3);

      pos.push_back(seam_pt[0]);
      norm.push_back(seam_pt[1]);
      axis.push_back(dir);
      type.push_back(4);
    }

  double eps = std::min(tol, 1.0e-5);
  vector<double> splitpar;
  for (size_t ki=0; ki<pos.size(); ++ki)
    {
      vector<double> intpar;
      vector<pair<double,double> > intcvs;
      intersectCurvePlane(cv.get(), pos[ki], norm[ki], eps,
			  intpar, intcvs);
      for (size_t kr=0; kr<intcvs.size(); ++kr)
	{
	  intpar.push_back(intcvs[kr].first);
	  intpar.push_back(intcvs[kr].second);
	}
      
      for (size_t kr=0; kr<intpar.size(); ++kr)
	{
	  // Check if the intersection is at the seam
	  Point pt = cv->point(intpar[kr]);
	  Point dir = axis[ki];
	  double rad = radius;
	  if (type[ki] == 2)
	    {
	      shared_ptr<Cone> cone = dynamic_pointer_cast<Cone, ParamSurface>(surf);
	      double upar, vpar, dist;
	      Point close;
	      cone->closestPoint(pt, upar, vpar, close, dist, eps);
	      rad = cone->radius(vpar, 0.0);
	      shared_ptr<Line> line = cone->getLine(upar);
	      dir = line->getDirection();
	    }
	  double dd = pt.dist(pos[ki]);
	  double dd2 = fabs((pt - pos[ki])*dir);
	  double dd3 = (pt - pos[ki])*norm[ki];
	  if (dd-dd2 < rad && dd3 < tol)
	    splitpar.push_back(intpar[kr]);  // Will need to consider for
	  // surfaces different from cylinder
	}
    }
  if (splitpar.size() > 1)
    {
      std::sort(splitpar.begin(), splitpar.end());
      for (int ka=(int)splitpar.size()-1; ka>0; --ka)
	if (splitpar[ka] - splitpar[ka-1] < eps)
	  {
	    splitpar[ka-1] = 0.5*(splitpar[ka-1] + splitpar[ka]);
	    splitpar.erase(splitpar.begin()+ka);
	  }
    }
  
  shared_ptr<ParamCurve> tmpcv = cv;
  for (size_t ki=0; ki<splitpar.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > split = tmpcv->split(splitpar[ki]);
      subcvs.push_back(split[0]);
      tmpcv = split[1];
    }
  subcvs.push_back(tmpcv);
  int stop_all = 1;
      
  // shared_ptr<ElementarySurface> elemsf =
  //   dynamic_pointer_cast<ElementarySurface, ParamSurface>(surf);
  // if (elemsf.get())
  //   {
  //     bool close_u, close_v;
  //     elemsf->isClosed(close_u, close_v);
  //     vector<shared_ptr<ParamCurve> > seam_cvs;
  //     vector<int> dir;
  //     if (close_u)
  // 	{
  // 	  vector<shared_ptr<ParamCurve> > tmp_cvs =
  // 	    elemsf->constParamCurves(dom.umin(), false);
  // 	  seam_cvs.push_back(tmp_cvs[0]);
  // 	  dir.push_back(1);
  // 	}
  //     if (close_v)
  // 	{
  // 	  vector<shared_ptr<ParamCurve> > tmp_cvs =
  // 	    elemsf->constParamCurves(dom.vmin(), true);
  // 	  seam_cvs.push_back(tmp_cvs[0]);
  // 	  dir.push_back(2);
  // 	}

  //     double eps = std::max(1.0e-4, 0.001*(cv->endparam()-cv->startparam()));
  //     subcvs.push_back(cv);
  //     for (size_t kj=0; kj<subcvs.size(); )
  // 	{
  // 	  bool atseam = false;
  // 	  for (size_t ki=0; ki<seam_cvs.size(); ++ki)
  // 	    {
  // 	      double par1, par2, dist;
  // 	      Point close1, close2;
  // 	      ClosestPoint::closestPtCurves(subcvs[kj].get(), seam_cvs[ki].get(),
  // 					    par1, par2, dist, close1, close2);
  // 	      if (dist > maxlim)
  // 		continue;
  // 	      if (par1 < subcvs[kj]->startparam()+eps ||
  // 		  par1 > subcvs[kj]->endparam()-eps)
  // 		continue;
	      
  // 	      // Check angle with surface normal
  // 	      double sf_u = (dir[ki] == 1) ? dom.umin() : par2;
  // 	      double sf_v = (dir[ki] == 1) ? par2 : dom.vmin();
  // 	      Point norm;
  // 	      surf->normal(norm, sf_u, sf_v);
  // 	      Point vec = close2 - close1;
  // 	      double ang = norm.angle(vec);
  // 	      ang = std::min(ang, M_PI-ang);
  // 	      if (ang < angtol)
  // 		{
  // 		  vector<shared_ptr<ParamCurve> > split =
  // 		    subcvs[kj]->split(par1);
  // 		  subcvs[kj] = split[0];
  // 		  subcvs.insert(subcvs.end(), split.begin()+1, split.end());
  // 		  atseam = true;
  // 		  break;
  // 		}
  // 	    }
  // 	  if (!atseam)
  // 	    ++kj;
  // 	}
  //   }
  // else
  //   subcvs.push_back(cv);
}

      
void RevEngRegion::extendBoundaries(double mean_edge_len, int min_point_reg,
				    double tol, double angtol, Point mainaxis[3])
//===========================================================================
{

#ifdef DEBUG_ADJUST
  std::ofstream of1("curr_regions_adjust.g2");
  std::ofstream of11("par_cvs_in.g2");
  writeRegionPoints(of1);
  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = trim_edgs_[ki]->geomCurve();
      shared_ptr<CurveOnSurface> sfcv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(of1);
	  sfcv->spaceCurve()->write(of1);
	  if (sfcv->hasParameterCurve())
	    SplineDebugUtils::writeSpaceParamCurve(sfcv->parameterCurve(),of11);
	}
    }
#endif
  
  int min_nmb_adj = min_point_reg/10;
  size_t min_bd = min_point_reg/50;
  double min_len = 20.0*mean_edge_len;
#ifdef DEBUG_ADJUST
  std::ofstream of5("all_adj_adjust.g2");
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if (commonRevEdge(*it))
	continue;

      if (commonTrimEdge(*it))
	continue;

      if ((*it)->hasBlendEdge())
	continue;   // Trimming curves should exist already
      // if ((*it)->numPoints() < min_nmb)
      // 	continue;
      (*it)->writeRegionPoints(of5);
    }
#endif

  vector<vector<RevEngPoint*> > bd_pts1; //, bd_pts2;
  vector<RevEngRegion*> adj_bd;
  vector<pair<RevEngPoint*,RevEngPoint*> > end_pts;
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      if (commonRevEdge(*it))
	continue;

      if (commonTrimEdge(*it))
	continue;

      if ((*it)->hasBlendEdge())
	continue;   // Trimming curves should exist already
      
      // The boundary towards adjacent regions with a surface should preferably be
      // handled by other tools, but are currently included
      // if ((*it)->hasSurface())
      // 	continue;
      // if ((*it)->numPoints() < min_nmb)
      // 	continue;

// #ifdef DEBUG_ADJUST
//       std::ofstream of2("adj_regions_adjust.g2");
//       (*it)->writeRegionPoints(of2);
// #endif
      
      // Get boundary points
      // Dismiss boundary points that are too close to existing trimming curves
      vector<RevEngPoint*> adj_pts1 = extractNextToAdjacent(*it);
      int nmb_lim = (int)adj_pts1.size()/10;
      int nmb_in = 0;
      for (size_t kj=0; kj<adj_pts1.size(); ++kj)
	{
	  Vector3D xyz = adj_pts1[kj]->getPoint();
	  Point pnt(xyz[0], xyz[1], xyz[2]);
	  for (size_t ki=0; ki<trim_edgs_.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = trim_edgs_[ki]->geomCurve();
	      double tpar, dist;
	      Point close;
	      cv->closestPoint(pnt, cv->startparam(), cv->endparam(), tpar, close, dist);
	      if (dist <= tol)
		{
		  nmb_in++;
		  break;
		}
	    }
	}
      
      //vector<RevEngPoint*> adj_pts2 = (*it)->extractNextToAdjacent(this);
      if (nmb_in < nmb_lim)
	{
	  RevEngPoint *end1, *end2;
	  RevEngPoint *dummy = 0;
	  vector<RevEngPoint*> tmp_pts;
	  for (size_t kj=0; kj<adj_pts1.size(); ++kj)
	    if (adj_pts1[kj]->getSurfaceDist() <= tol)
	      tmp_pts.push_back(adj_pts1[kj]);
	  bd_pts1.push_back(tmp_pts); //adj_pts1);
	  adj_bd.push_back(*it);
	  if (adj_pts1.size() > 0)
	    {
	      endPoints(adj_pts1, end1, end2);
	      end_pts.push_back(std::make_pair(end1, end2));
	    }
	  else
	    end_pts.push_back(std::make_pair(dummy, dummy));
	}
      //bd_pts2.push_back(adj_pts2);
    }
      // if (adj_pts1.size() < min_bd || adj_pts2.size() < min_bd)
      // 	continue;

  
  vector<vector<pair<RevEngPoint*,int> > > branches;
  getBranches(bd_pts1, branches);
  
#ifdef DEBUG_ADJUST
  std::ofstream of2("bd1.g2");
  std::ofstream of4("end1.g2");
  for (size_t kj=0; kj<bd_pts1.size(); ++kj)
    if (bd_pts1[kj].size() > 0)
      {
	of4 << "400 1 0 4 255 0 0 255" << std::endl;
	of4 << "2" << std::endl;
	of4 << end_pts[kj].first->getPoint() << std::endl;
	of4 << end_pts[kj].second->getPoint() << std::endl;
	of2 << "400 1 0 4 0 255 0 255" << std::endl;
	of2 << bd_pts1[kj].size() << std::endl;
	for (size_t ki=0; ki<bd_pts1[kj].size(); ++ki)
	  of2 << bd_pts1[kj][ki]->getPoint() << std::endl;
      }
  
  // std::ofstream of3("bd2.g2");
  // for (size_t kj=0; kj<bd_pts2.size(); ++kj)
  //   if (bd_pts2[kj].size() > 0)
  //     {
  // 	of3 << "400 1 0 4 0 255 0 255" << std::endl;
  // 	of3 << bd_pts2[kj].size() << std::endl;
  // 	for (size_t ki=0; ki<bd_pts2[kj].size(); ++ki)
  // 	  of3 << bd_pts2[kj][ki]->getPoint() << std::endl;
  //     }

  std::ofstream of6("bd_branches.g2");
  for (size_t ki=0; ki<branches.size(); ++ki)
    if (branches[ki].size() > 0 && branches[ki][0].first)
      {
	of6 << "400 1 0 4 255 0 0 255" << std::endl;
	of6 << branches[ki].size() << std::endl;
	for (size_t kj=0; kj<branches[ki].size(); ++kj)
	  of6 << branches[ki][kj].first->getPoint() << std::endl;
      }
#endif

  // Sort sequences and remember branch points
  vector<vector<int> > ixs;
  vector<vector<RevEngPoint*> > joints;
  vector<int> num_bd;
  for (size_t ki=0; ki<branches.size(); ++ki)
    {
      size_t kr,kh;
      for (kr=0; kr<ixs.size(); ++kr)
	{
	  for (kh=0; kh<ixs[kr].size(); ++kh)
	    {
	      if (ixs[kr][kh] == (int)ki)
		break;
	    }
	  if (kh < ixs[kr].size())
	    break;
	}
      if (kr < ixs.size())
	continue;

      if (branches[ki].size() == 0)
	continue;
      
      vector<int> ixs0;
      vector<RevEngPoint*> joints0;
      ixs0.push_back((int)ki);
      int curr=(int)ki, next=-1;
      RevEngPoint *pt;
      int num0 = 0;
      next = getNextSeq(ixs0, curr, branches, bd_pts1, pt);
      while (next >= 0)
	{
	  ixs0.push_back(next);
	  joints0.push_back(pt);
	  num0 += (int)branches[curr].size();
	  curr = next;
	  next = getNextSeq(ixs0, curr, branches, bd_pts1, pt);
	}
      if (ixs0.size() > 1)
	{
	  curr = ixs0[0];
	  next = getNextSeq(ixs0, curr, branches, bd_pts1, pt);
	  while (next >= 0)
	    {
	      ixs0.insert(ixs0.begin(), next);
	      joints0.insert(joints0.begin(), pt);
	      num0 += (int)branches[curr].size();
	      curr = next;
	      next = getNextSeq(ixs0, curr, branches, bd_pts1, pt);
	    }
	}
      ixs.push_back(ixs0);
      joints.push_back(joints0);
      num_bd.push_back(num0);
    }

  for (size_t ki=0; ki<ixs.size(); ++ki)
    for (size_t kj=ki+1; kj<ixs.size(); ++kj)
      if (num_bd[kj] > num_bd[ki])
	{
	  std::swap(ixs[ki], ixs[kj]);
	  std::swap(joints[ki], joints[kj]);
	  std::swap(num_bd[ki], num_bd[kj]);
	}
  
#ifdef DEBUG_ADJUST
  std::ofstream of7("joined_seqs.g2");
  for (size_t ki=0; ki<ixs.size(); ++ki)
    {
      int num=0;
      for (size_t kj=0; kj<ixs[ki].size(); ++kj)
	num += (int)bd_pts1[ixs[ki][kj]].size();

      of7 << "400 1 0 0" << std::endl;
      of7 << num << std::endl;
      for (size_t kj=0; kj<ixs[ki].size(); ++kj)
	for (size_t kr=0; kr<bd_pts1[ixs[ki][kj]].size(); ++kr)
	  of7 << bd_pts1[ixs[ki][kj]][kr]->getPoint() << std::endl;
    }
#endif

  // Check if the points can be approximated by a circle or a line
  // First collect seqences
  vector<bool> done(ixs.size(), false);
  vector<vector<size_t> > all_bd_ixs;
  vector<int> all_num;
  vector<vector<Point> > all_points;
  vector<int> all_adj_num;
  vector<Point> all_norm;
  for (size_t ki=0; ki<ixs.size(); ++ki)
    {
      if (done[ki])
	continue;

      // Collect points
      vector<size_t> bd_ixs;
      size_t kj=ki;
      int num = 0;
      while (kj < ixs.size())
	{
	  done[kj] = true;
	  for (size_t kr=0; kr<ixs[kj].size(); ++kr)
	    {
	      size_t kh;
	      for (kh=0; kh<bd_ixs.size(); ++kh)
		if (bd_ixs[kh] == ixs[kj][kr])
		  break;
	      if (kh == bd_ixs.size())
		{
		  bd_ixs.push_back(ixs[kj][kr]);
		  num += (int)bd_pts1[ixs[kj][kr]].size();
		}
	    }
	  
	  for (++kj; kj<ixs.size(); ++kj)
	    {
	      size_t kr;
	      for (kr=0; kr<ixs[kj].size(); ++kr)
		{
		  size_t kh;
		  for (kh=0; kh<bd_ixs.size(); ++kh)
		    if (bd_ixs[kh] == ixs[kj][kr])
		      break;
		  if (kh < bd_ixs.size())
		    break;
		}
	      if (kr < ixs[kj].size())
		break;
	    }
	}

      vector<Point> points;
      points.reserve(num);
      Point norm(0.0, 0.0, 0.0);
      double fac = 1.0/(double)num;
      int num_adj = 0;
      for (size_t kr=0; kr<bd_ixs.size(); ++kr)
	{
	  for (size_t kh=0; kh<bd_pts1[bd_ixs[kr]].size(); ++kh)
	    {
	      Vector3D xyz = bd_pts1[bd_ixs[kr]][kh]->getPoint();
	      points.push_back(Point(xyz[0], xyz[1], xyz[2]));
	      Point norm0 = bd_pts1[bd_ixs[kr]][kh]->getTriangNormal();
	      norm += fac*norm0;
	    }
	  num_adj += adj_bd[bd_ixs[kr]]->numPoints();
	}
      all_bd_ixs.push_back(bd_ixs);
      all_num.push_back(num);
      all_points.push_back(points);
      all_norm.push_back(norm);
      all_adj_num.push_back(num_adj);
    }

  int num_lim = 10;
  vector<shared_ptr<ElementaryCurve> > elemcv(2*all_points.size());
  for (size_t ki=0; ki<all_points.size(); ++ki)
  {
    if (all_points[ki].size() < num_lim)
      continue;
    
    // Compute circle
    Point pos0, axis, Cx, Cy;
    RevEngUtils::computePlane(all_points[ki], all_norm[ki], mainaxis, 
			      pos0, axis, Cx, Cy);

    Point centre;
    double radius;
    try {
      RevEngUtils::computeCircPosRadius(all_points[ki], axis, Cx, Cy,
					centre, radius);
    }
    catch (...)
      {
	continue;
      }

    shared_ptr<Circle> circ(new Circle(radius, centre, axis, Cx));
    elemcv[2*ki] = circ;

    // Compute line
    Point mid, dir;
    RevEngUtils::computeLine(all_points[ki], mid, dir);
    shared_ptr<Line> line(new Line(mid, dir));
    elemcv[2*ki+1] = line;
    
    int stop_circ = 1;
  }
#ifdef DEBUG_ADJUST
  std::ofstream of8("circle_line.g2");
  for (size_t ki=0; ki<elemcv.size(); ++ki)
    {
      if (elemcv[ki].get())
	{
	  elemcv[ki]->writeStandardHeader(of8);
	  elemcv[ki]->write(of8);
	}
    }
#endif
  
  // Check distance
  vector<double> maxdist(all_points.size(), std::numeric_limits<double>::max());
  vector<double> avdist(all_points.size(), std::numeric_limits<double>::max());
  vector<int> num_inside(all_points.size(), 0);
  vector<vector<double> > param(all_points.size());
  vector<vector<double> > distance(all_points.size());
  vector<int> match(all_points.size(), -1);

  double fac1 = 0.9;
  for (size_t ki=0; ki<all_points.size(); ++ki)
  {
    if (all_points[ki].size() < num_lim)
      continue;
    for (size_t kj=0; kj<elemcv.size(); ++kj)
      {
	if (!elemcv[kj].get())
	  continue;
	double maxd, avd;
	int in;
	vector<double> parvals, dist;
	RevEngUtils::distToCurve(all_points[ki], elemcv[kj], tol, maxd,
				 avd, in, parvals, dist);
	if ((avd < avdist[ki] &&
	     ((double)in > fac1*(double)num_inside[ki] || num_inside[ki] == 0)) ||
	    (in > num_inside[ki] && fac1*avd < avdist[ki]) )
	  {
	    maxdist[ki] = maxd;
	    avdist[ki] = avd;
	    num_inside[ki] = in;
	    param[ki] = parvals;
	    distance[ki] = dist;
	    match[ki] = (int)kj;
	  }
      }
    int stop_dist = 1;
  }

#ifdef DEBUG_ADJUST
  std::ofstream of9("spline.g2");
#endif

  shared_ptr<ParamSurface> surf = getSurface(0)->surface();
  double smoothwgt = 0.5;
  int in=6, ik=4;
  int maxiter = 2;
#ifdef DEBUG_ADJUST
  std::ofstream of10("proj_cvs.g2");
  surf->writeStandardHeader(of10);
  surf->write(of10);
#endif
  for (size_t kj=0; kj<elemcv.size(); ++kj)
    {
      if (!elemcv[kj].get())
	continue;
      
      // Collect points associated to this curve
      vector<int> pt_ix;
      for (size_t ki=0; ki<match.size(); ++ki)
	if (match[ki] == (int)kj)
	  pt_ix.push_back((int)ki);

      if (pt_ix.size() == 0)
	continue;

      vector<Point> pts;
      vector<double> pts2;
      vector<double> par;
      int curr_adj_num = 0;
      for (size_t ki=0; ki<pt_ix.size(); ++ki)
	{
	  pts.insert(pts.end(), all_points[pt_ix[ki]].begin(),
		     all_points[pt_ix[ki]].end());
	  par.insert(par.end(), param[pt_ix[ki]].begin(), param[pt_ix[ki]].end());
	  curr_adj_num += all_adj_num[pt_ix[ki]];
	}

      // Check significance of current boundary points
      if ((int)pts.size() < min_bd || curr_adj_num < min_nmb_adj)
	continue;
      
      for (size_t ki=0; ki<par.size(); ++ki)
	for (size_t kr=ki+1; kr<par.size(); ++kr)
	  if (par[kr] < par[ki])
	    {
	      std::swap(par[kr], par[ki]);
	      std::swap(pts[kr], pts[ki]);
	    }

      bool close = false;
      if (kj%2 == 0)
	{
	  // Circular base curve. Check for internal seam
	  double closefac = 0.1*M_PI;
	  double avdd = 2*M_PI/(double)(par.size()-1);
	  double ddlim = 20.0*avdd;
	  double maxdd = 0;
	  size_t max_ix = 0;
	  for (size_t ki=1; ki<par.size(); ++ki)
	    if (par[ki] - par[ki-1] > maxdd)
	      {
		maxdd = par[ki] - par[ki-1];
		max_ix = ki;
	      }
	    if (max_ix > 0 && maxdd > ddlim)
	      {
		size_t nn = par.size()-max_ix;
		par.insert(par.end(), par.begin(), par.begin()+max_ix);
		pts.insert(pts.end(), pts.begin(), pts.begin()+max_ix);
		par.erase(par.begin(), par.begin()+max_ix);
		pts.erase(pts.begin(), pts.begin()+max_ix);
		for (size_t kr=nn; kr<par.size(); ++kr)
		  par[kr] += 2*M_PI;
		
	      }

	  if (par[par.size()-1] - par[0] > 2*M_PI-closefac)
	    close = true;
 	}


      double tmin = par[0];
      double tmax = par[par.size()-1];
      for (size_t ki=1; ki<par.size(); ++ki)
	{
	  double dd = pts[ki-1].dist(pts[ki]);
	  par[ki] = par[ki-1] + dd;
	}

      if (close)
	{
	  size_t num = par.size();
	  pts.push_back(pts[0]);  // Repeat point at seam
	  par.push_back(par[num-1] + pts[num-1].dist(pts[num]));
	}


      BoundingBox ptbb(3);
      for (size_t ki=0; ki<pts.size(); ++ki)
	{
	  ptbb.addUnionWith(pts[ki]);
	  pts2.insert(pts2.end(), pts[ki].begin(), pts[ki].end());
	}
      
      
      ApproxCurve approx(pts2, par, 3, tol, in, ik);
      approx.setSmooth(smoothwgt);

      double maxspl, avspl;
      shared_ptr<SplineCurve> cv = approx.getApproxCurve(maxspl, avspl, maxiter);
      BoundingBox cvbb = cv->boundingBox();
#ifdef DEBUG_ADJUST
      cv->writeStandardHeader(of9);
      cv->write(of9);
#endif
      shared_ptr<ParamCurve> selcv = cv;
      double maxlim = maxspl;
      double avelem = avdist[pt_ix[0]];
      int inelem = num_inside[pt_ix[0]];
      double maxelem = maxdist[pt_ix[0]];
      if (pt_ix.size() > 1)
	{
	  size_t num_pt = par.size();
	  avelem = 0;
	  inelem = 0;
	  for (size_t kr=0; kr<pt_ix.size(); ++kr)
	    {
	      avelem += (double)all_points[pt_ix[kr]].size()*avdist[pt_ix[kr]]/
		(double)num_pt;
	      inelem += num_inside[pt_ix[kr]];
	      maxelem = std::max(maxelem, maxdist[pt_ix[kr]]);
	    }
	}

      double bbfac = 5.0;
      if (avelem < avspl || (avelem <= tol && 2*inelem > (int)par.size()) ||
	  cvbb.low().dist(cvbb.high()) > bbfac*ptbb.low().dist(ptbb.high()))
	{
	  if (!close)
	    {
	      if (tmin < -2*M_PI)
		{
		  tmin += 2*M_PI;
		  tmax += 2*M_PI;
		}
	      if (tmax > 2*M_PI)
		{
		  tmin -= 2*M_PI;
		  tmax -= 2*M_PI;
		}
	    }
	    elemcv[kj]->setParamBounds(tmin, tmax);
	  selcv = elemcv[kj];
	  maxlim = maxelem;
	}

      // Check if the curve crosses a seam
      vector<shared_ptr<ParamCurve> > subcvs;
      splitAtSeam(surf, selcv, tol, angtol, maxlim, subcvs);
      
      vector<shared_ptr<CurveOnSurface> > trim_cv(subcvs.size());
      for (size_t kr=0; kr<subcvs.size(); ++kr)
	{
	  shared_ptr<SplineCurve> proj_cv, par_cv;
	  CurveCreators::projectCurve(subcvs[kr], surf, tol, proj_cv, 
				      par_cv);
	  trim_cv[kr] = shared_ptr<CurveOnSurface>(new CurveOnSurface(surf, proj_cv,
								      false));
	  trim_cv[kr]->ensureParCrvExistence(tol);
#ifdef DEBUG_ADJUST
	  proj_cv->writeStandardHeader(of10);
	  proj_cv->write(of10);
	  shared_ptr<ParamCurve> tmp_par = trim_cv[kr]->parameterCurve();
	  if (tmp_par.get())
	    SplineDebugUtils::writeSpaceParamCurve(tmp_par, of11);
#endif
	}

      for (size_t ki=0; ki<trim_cv.size(); ++ki)
	{
	  double len = trim_cv[ki]->estimatedCurveLength();
	  if (len < min_len)
	    continue;
	  shared_ptr<ftEdge> edg(new ftEdge(associated_sf_[0], trim_cv[ki],
					    trim_cv[ki]->startparam(),
					    trim_cv[ki]->endparam()));
	  addTrimEdge(edg);
	}
      int stop_spl = 1;
    }

      
  int stop_break2 = 1;
}

//===========================================================================
vector<RevEngPoint*>  RevEngRegion::extractBdOutPoints(shared_ptr<SplineCurve>& crv,
						       vector<RevEngPoint*>& seq_pts,
						       double tol)
//===========================================================================
{
  vector<RevEngPoint*> left, right, upper, lower;
  vector<RevEngPoint*> points;
  points.insert(points.end(), seq_pts.begin(), seq_pts.end());
  size_t seq_size = seq_pts.size();

  for (size_t ki=0; ki<points.size(); ++ki)
    points[ki]->setVisited();

  double maxdist = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      if (ki == seq_size)
	maxdist *= 2.0;
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double tpar, dist;
      Point close;
      crv->ParamCurve::closestPoint(pos, tpar, close, dist);
      if (dist > maxdist)
	{
	  if (ki < seq_size)
	    maxdist = dist;
	  else
	    continue;
	}
      
      if (tpar <= crv->startparam())
	lower.push_back(points[ki]);
      else if (tpar >= crv->endparam())
	upper.push_back(points[ki]);
      else
	{
	  vector<Point> der(2);
	  crv->point(der, tpar, 1);
	  Point norm = points[ki]->getLocFuncNormal();
	  Point vec = pos - close;
	  Point vec2 = vec.cross(der[1]);
	  if (vec2*norm >= 0)
	    left.push_back(points[ki]);
	  else
	    right.push_back(points[ki]);
	}
      
      vector<ftSamplePoint*> next = points[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *curr = dynamic_cast<RevEngPoint*>(next[kj]);
	  RevEngRegion *adj_reg = curr->region();
	  if (curr->visited())
	    continue;
	  if (adj_reg != this)
	    continue;
	  curr->setVisited();
	  points.push_back(curr);
	}
    }

  for (size_t ki=0; ki<points.size(); ++ki)
    points[ki]->unsetVisited();
  
  return (left.size() < right.size()) ? left : right;
}

//===========================================================================
vector<RevEngPoint*>  RevEngRegion::sortPtsSeq(double mean_edge_len,
					       vector<RevEngPoint*>& seq_pts,
					       vector<RevEngPoint*>& sub_pts)
//===========================================================================
{
  vector<RevEngPoint*> sub_pts2;
  if (sub_pts.size() < 2)
    {
      // Check for maximum distance between input points
      vector<RevEngPoint*> dummy;
      if (seq_pts.size() < 2)
	return dummy;

      double fac = 10.0;
      size_t ix1 = 0, ix2 = 1;
      ix2 = std::min(ix2, seq_pts.size()-1);
      double maxlen = seq_pts[ix1]->pntDist(seq_pts[ix2]);
      for (size_t ki=0; ki<seq_pts.size(); ++ki)
	for (size_t kj=ki+1; kj<seq_pts.size(); ++kj)
	  {
	    double len = seq_pts[ki]->pntDist(seq_pts[kj]);
	    if (len > maxlen)
	      {
		maxlen = len;
		ix1 = ki;
		ix2 = kj;
	      }
	  }

      if (sub_pts.size() == 1)
	{
	  sub_pts2.push_back(sub_pts[0]);
	  double len1 = sub_pts[0]->pntDist(seq_pts[ix1]);
	  double len2 = sub_pts[0]->pntDist(seq_pts[ix2]);
	  if (len1 >= len2)
	    sub_pts2.push_back(seq_pts[ix1]);
	  else
	    sub_pts2.push_back(seq_pts[ix2]);
	}
      else
	{
	  sub_pts2.push_back(seq_pts[ix1]);
	  if (maxlen > fac*mean_edge_len)
	    sub_pts2.push_back(seq_pts[ix2]);
	}
    }
  else
    sub_pts2.insert(sub_pts2.end(), sub_pts.begin(), sub_pts.end());
  vector<RevEngPoint*> seq_pts2(seq_pts.begin(), seq_pts.end());
  vector<vector<RevEngPoint*> > all_seq;
  
  while (sub_pts2.size() > 0)
    {
      size_t ix1 = 0, ix2 = 1;
      ix2 = std::min(ix2, sub_pts2.size()-1);
      double maxlen = sub_pts2[ix1]->pntDist(sub_pts2[ix2]);
      for (size_t ki=0; ki<sub_pts2.size(); ++ki)
	for (size_t kj=ki+1; kj<sub_pts2.size(); ++kj)
	  {
	    double len = sub_pts2[ki]->pntDist(sub_pts2[kj]);
	    if (len > maxlen)
	      {
		maxlen = len;
		ix1 = ki;
		ix2 = kj;
	      }
	  }

      vector<RevEngPoint*> sortseq;
      sortseq.push_back(sub_pts2[ix1]);
      bool more_pts = true;
      RevEngPoint *curr = sub_pts2[ix1];
      vector<RevEngPoint*> prev_pts;
      while (more_pts)
	{
	  vector<RevEngPoint*> next_pts;
	  vector<ftSamplePoint*> next = curr->getNeighbours();
	  for (size_t ki=0; ki<next.size(); ++ki)
	    {
	      size_t kj;
	      for (kj=0; kj<next_pts.size(); ++kj)
		if (next_pts[kj] == next[ki])
		  break;
	      if (kj == next_pts.size())
		{
		  size_t kh;
		  for (kh=0; kh<sortseq.size(); ++kh)
		    if (sortseq[kh] == next[ki])
		      break;
		  if (kh == sortseq.size())
		    {
		      size_t kr;
		      for (kr=0; kr<seq_pts2.size(); ++kr)
			if (seq_pts2[kr] == next[ki])
			  break;
		      if (kr < seq_pts2.size())
			next_pts.push_back(dynamic_cast<RevEngPoint*>(next[ki]));
		    }
		}
	    }
	  for (size_t ki=0; ki<next_pts.size(); ++ki)
	    for (size_t kj=ki+1; kj<next_pts.size(); )
	      {
		if (next_pts[ki] == next_pts[kj])
		  next_pts.erase(next_pts.begin()+kj);
		else
		  ++kj;
	      }

	  if (next_pts.size() == 0)
	    next_pts = prev_pts;
	  
	  if (next_pts.size() == 0)
	    break;
	  if (next_pts.size() == 1)
	    {
	      prev_pts.clear();
	      sortseq.push_back(next_pts[0]);
	    }
	  else
	    {
	      double minlen = curr->pntDist(next_pts[0]);
	      size_t ix3 = 0;
	      for (size_t kr=1; kr<next_pts.size(); ++kr)
		{
		  double len = curr->pntDist(next_pts[kr]);
		  if (len < minlen)
		    {
		      minlen = len;
		      ix3 = kr;
		    }
		}
	      prev_pts.clear();
	      sortseq.push_back(next_pts[ix3]);
	      for (size_t kr=0; kr<next_pts.size(); ++kr)
		if (kr != ix3)
		  prev_pts.push_back(next_pts[kr]);
	    }
	  curr = sortseq[sortseq.size()-1];
	}

#ifdef DEBUG_SEGMENT
      std::ofstream of("sorted_seq.g2");
      for (size_t ki=0; ki<sortseq.size(); ++ki)
	{
	  of << "400 1 0 4 0 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << sortseq[ki]->getPoint() << std::endl;
	}
#endif
      int stop_break = 1;

      for (size_t ki=0; ki<sortseq.size(); ++ki)
	{
	  for (size_t kj=0; kj<sub_pts2.size(); ++kj)
	    {
	      if (sortseq[ki] == sub_pts2[kj])
		{
		  sub_pts2.erase(sub_pts2.begin()+kj);
		  break;
		}
	      for (size_t kj=0; kj<seq_pts2.size(); ++kj)
		{
		  if (sortseq[ki] == seq_pts2[kj])
		    {
		      seq_pts2.erase(seq_pts2.begin()+kj);
		      break;
		    }
		}
	    }
	}
      all_seq.push_back(sortseq);
    }

  // Join sub sequences
  vector<RevEngPoint*> sorted;
  if (all_seq.size() > 0)
    {
      int all_size = (int)all_seq.size();
      while (all_size > 1)
	{
	  double t1min = std::numeric_limits<double>::max();
	  double t2min = std::numeric_limits<double>::max();
	  int x1 = -1, x2 = -1;
	  bool turn1 = false, turn2 = false;
	  for (int kj=1; kj<all_size; ++kj)
	    {
	      double l1 = all_seq[0][0]->pntDist(all_seq[kj][0]);
	      double l2 = all_seq[0][0]->pntDist(all_seq[kj][all_seq[kj].size()-1]);
	      double l3 = all_seq[0][all_seq[0].size()-1]->pntDist(all_seq[kj][0]);
	      double l4 = all_seq[0][all_seq[0].size()-1]->pntDist(all_seq[kj][all_seq[kj].size()-1]);
	      if (std::min(l1, l2) < t1min)
		{
		  t1min = std::min(l1, l2);
		  x1 = (int)kj;
		  turn1 = (l1 < l2);
		}
	    
	      if (std::min(l3, l4) < t2min)
		{
		  t2min = std::min(l3, l4);
		  x2 = (int)kj;
		  turn2 = (l4 < l3);
		}
	    }
	  if (t1min < t2min)
	    {
	      if (turn1)
		{
		  size_t last = all_seq[x1].size()-1;
		  for (size_t kr=0; kr<all_seq[x1].size()/2; ++kr)
		    std::swap(all_seq[x1][kr], all_seq[x1][last-kr]);
		}
	      all_seq[0].insert(all_seq[0].begin(), all_seq[x1].begin(),
				all_seq[x1].end());
	      if (x1 < all_size-1)
		std::swap(all_seq[x1],all_seq[all_size-1]);
	    }
	  else
	    {
	      if (turn2)
		{
		  size_t last = all_seq[x2].size()-1;
		  for (size_t kr=0; kr<all_seq[x2].size()/2; ++kr)
		    std::swap(all_seq[x2][kr], all_seq[x2][last-kr]);
		}
	      all_seq[0].insert(all_seq[0].end(), all_seq[x2].begin(),
				all_seq[x2].end());
	      if (x2 < all_size-1)
		std::swap(all_seq[x2],all_seq[all_size-1]);
	    }
	  all_size--;
	}
    }

  bool include_missing = false;
  if (include_missing)
    {
      // Include missing input points
      for (size_t ki=0; ki<seq_pts2.size(); ++ki)
	{
	  double minlen = seq_pts2[ki]->pntDist(all_seq[0][0]);
	  size_t min_ix = 0;
	  for (size_t kj=1; kj<all_seq[0].size(); ++kj)
	    {
	      double len = seq_pts2[ki]->pntDist(all_seq[0][kj]);
	      if (len < minlen)
		{
		  min_ix = kj;
		  minlen = len;
		}
	    }

	  double len1 = (min_ix == 0) ? 0.0 : seq_pts2[ki]->pntDist(all_seq[0][min_ix-1]);
	  double len2 = (min_ix == all_seq[0].size()-1) ? 0.0 :
	    seq_pts2[ki]->pntDist(all_seq[0][min_ix+1]);
	  size_t ix = min_ix + (len2 > len1);
	  all_seq[0].insert(all_seq[0].begin()+ix, seq_pts2[ki]);
    }
    }
  return all_seq[0];
}

//===========================================================================
void RevEngRegion::adjustWithSurf(Point mainaxis[3], int min_pt_reg,
				  double tol, double angtol)
//===========================================================================
{
  //double eps = 1.0e-6;
  if (associated_sf_.size() == 0)
    return;  // No surface with which to check growt

  //int sfcode;
  //ClassType classtype = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<BoundedSurface> bdsurf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (bdsurf.get())
    surf = bdsurf->underlyingSurface();

#ifdef DEBUG_ADJUST  
  std::ofstream res1("residuals_source.txt");
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    res1 << group_points_[ki]->getPoint() << " " << group_points_[ki]->getSurfaceDist() << std::endl;
#endif
  
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      // Check criteria for growt
      //int num = (*it)->numPoints();
      // if (num > group_points_.size())
      // 	continue;
      
      std::ofstream res2("residuals_adj.txt");
      for (size_t kh=0; kh<(*it)->group_points_.size(); ++kh)
	res2 << (*it)->group_points_[kh]->getPoint() << " " << (*it)->group_points_[kh]->getSurfaceDist() << std::endl;
  
#ifdef DEBUG_SEGMENT
      int num = (*it)->numPoints();
       std::ofstream of("curr_extend.g2");
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << group_points_.size() << std::endl;
      for (size_t kh=0; kh<group_points_.size(); ++kh)
      	of << group_points_[kh]->getPoint() << std::endl;
      
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << num << std::endl;
      for (int ka=0; ka<num; ++ka)
      	of << (*it)->getPoint(ka)->getPoint() << std::endl;
      
      surf->writeStandardHeader(of);
      surf->write(of);
#endif
      
      vector<RevEngPoint*> in, out;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      double maxd, avd;
      int num_inside, num2_inside;
      RevEngUtils::distToSurf((*it)->pointsBegin(),
			      (*it)->pointsEnd(), surf, tol,
			      maxd, avd, num_inside, num2_inside, in, out,
			      parvals, dist_ang);
      vector<RevEngPoint*> better1, worse1;
#ifdef DEBUG_SEGMENT
      std::ofstream res3("residuals_adj2.txt");
      for (size_t kh=0; kh<(*it)->group_points_.size(); ++kh)
	res3 << (*it)->group_points_[kh]->getPoint() << " " << dist_ang[kh].first << std::endl;
#endif
      
      for (size_t kh=0; kh<dist_ang.size(); ++kh)
	{
	  double ptdist = tol, ptang = angtol;
	  if ((*it)->hasSurface())
	    (*it)->group_points_[kh]->getSurfaceDist(ptdist, ptang);
	  if (dist_ang[kh].first <= ptdist && dist_ang[kh].second < angtol) //ptang)
	  //if (dist_ang[kh].first <= tol)
	    better1.push_back((*it)->group_points_[kh]);
	  else
	    worse1.push_back((*it)->group_points_[kh]);
	}
#ifdef DEBUG_SEGMENT
      std::ofstream ofc1("better_worse1.g2");
      ofc1 << "400 1 0 4 155 50 50 255" << std::endl;
      ofc1 << better1.size() << std::endl;
      for (size_t kr=0; kr<better1.size(); ++kr)
	ofc1 << better1[kr]->getPoint() << std::endl;
      ofc1 << "400 1 0 4 50 155 50 255" << std::endl;
      ofc1 << worse1.size() << std::endl;
      for (size_t kr=0; kr<worse1.size(); ++kr)
	ofc1 << worse1[kr]->getPoint() << std::endl;
      
      std::ofstream ofd("in_out_extend.g2");
      ofd << "400 1 0 4 155 50 50 255" << std::endl;
      ofd << in.size() << std::endl;
      for (size_t kr=0; kr<in.size(); ++kr)
	ofd << in[kr]->getPoint() << std::endl;
      ofd << "400 1 0 4 50 155 50 255" << std::endl;
      ofd << out.size() << std::endl;
      for (size_t kr=0; kr<out.size(); ++kr)
	ofd << out[kr]->getPoint() << std::endl;
#endif
      
      vector<RevEngPoint*> better2, worse2;
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> surf2 = (*it)->getSurface(0)->surface();
	  double maxd2, avd2;
	  int num_inside2, num2_inside2;
	  vector<RevEngPoint*> in2, out2;
	  vector<pair<double, double> > dist_ang2;
	  vector<double> parvals2;
	  RevEngUtils::distToSurf(pointsBegin(),
				  pointsEnd(), surf2, tol,
				  maxd2, avd2, num_inside2, num2_inside2, in2, out2,
				  parvals2, dist_ang2);
	  for (size_t kh=0; kh<dist_ang2.size(); ++kh)
	    {
	      double ptdist, ptang;
	      group_points_[kh]->getSurfaceDist(ptdist, ptang);
	      if (dist_ang2[kh].first <= ptdist && dist_ang2[kh].second < angtol) //ptang)
		{
		  better2.push_back(group_points_[kh]);
		}
	      else
		worse2.push_back(group_points_[kh]);
	}
#ifdef DEBUG_SEGMENt
	  std::ofstream ofc2("better_worse2.g2");
	  ofc2 << "400 1 0 4 155 50 50 255" << std::endl;
	  ofc2 << better2.size() << std::endl;
	  for (size_t kr=0; kr<better2.size(); ++kr)
	    ofc2 << better2[kr]->getPoint() << std::endl;
	  ofc2 << "400 1 0 4 50 155 50 255" << std::endl;
	  ofc2 << worse2.size() << std::endl;
	  for (size_t kr=0; kr<worse2.size(); ++kr)
	    ofc2 << worse2[kr]->getPoint() << std::endl;
#endif
	}

      // Move points
      bool move2 = false;
      size_t b1size = 2*better1.size();
      while (better1.size() < b1size)
	{
	  b1size = better1.size();
	  for (size_t ki=0; ki<better1.size(); )
	    {
	      if (better1[ki]->isNeighbour(this))
		{
		  move2 = true;
		  (*it)->removePoint(better1[ki]);
		  better1[ki]->setRegion(this);
		  better1[ki]->addMove();
		  group_points_.push_back(better1[ki]);
		  better1.erase(better1.begin()+ki);
		}
	      else
		++ki;
	    }
	}

      if (false)
	{
      size_t b2size = 2*better2.size();
      while (better2.size() < b2size)
	{
	  b2size = better2.size();
	  for (size_t ki=0; ki<better2.size(); )
	    {
	      if (better2[ki]->isNeighbour(*it))
		{
		  move2 = true;
		  removePoint(better2[ki]);
		  better2[ki]->setRegion(*it);
		  better2[ki]->addMove();
		  (*it)->addPoint(better2[ki]);
		  better2.erase(better2.begin()+ki);
		}
	      else
		++ki;
	    }
	}
	}
     
      if ((*it)->hasSurface() && (*it)->numPoints() > 5 && move2)
	{
	  // Check if the surface should be updated
	  (*it)->checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
	}
    }
  

  // Check if the surface should be updated
  checkReplaceSurf(mainaxis, min_pt_reg, tol, angtol);
#ifdef DEBUG_SEGMENt
  std::ofstream res4("residuals_res.txt");
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    res4 << group_points_[ki]->getPoint() << " " << group_points_[ki]->getSurfaceDist() << std::endl;
#endif
  int stop_break = 1;
  
}
  
//===========================================================================
bool RevEngRegion::getCurveRestriction(vector<shared_ptr<CurveOnSurface> >& cvs,
				       double tol, double anglim,
				       vector<pair<double,double> >& endpars)
//===========================================================================
{
  double int_tol = 1.0e-6;
  Point corner[4];
  corner[0] = Point(domain_[0], domain_[2]);
  corner[1] = Point(domain_[1], domain_[2]);
  corner[2] = Point(domain_[0], domain_[3]);
  corner[3] = Point(domain_[1], domain_[3]);
  
  // Ensure parameter representation of curves
  endpars.resize(cvs.size());
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      double start = std::numeric_limits<double>::max();
      double end = std::numeric_limits<double>::lowest();
      if (!(cvs[ki]->hasParameterCurve()))
	{
	  shared_ptr<ParamSurface> surf = cvs[ki]->underlyingSurface();
	  if (!surf->isBounded())
	    {
	      double diag = bbox_.low().dist(bbox_.high());
	      associated_sf_[0]->limitSurf(2*diag);
	    }
	  bool OK = cvs[ki]->ensureParCrvExistence(tol);
	  if (!OK)
	    {
	      // Curve larger than surface. Get initial restriction
	      shared_ptr<ParamSurface> surf = cvs[ki]->underlyingSurface();
	      CurveLoop cvloop = surf->outerBoundaryLoop();
	      vector<shared_ptr<ParamCurve> > lcvs = cvloop.getCurves();
	      shared_ptr<ParamCurve> scurve = cvs[ki]->spaceCurve();
	      vector<double> parvals;
	      for (size_t kr=0; kr<lcvs.size(); ++kr)
		{
		  vector<pair<double,double> > intpts;
		  intersectParamCurves(scurve.get(), lcvs[kr].get(), int_tol,
				       intpts);
		  for (size_t kj=0; kj<intpts.size(); ++kj)
		    parvals.push_back(intpts[kj].first);
		}
	      std::sort(parvals.begin(), parvals.end());
	      double t1 = cvs[ki]->startparam();
	      double t2 = cvs[ki]->endparam();
	      if (parvals.size() > 1)
		{
		  t1 = parvals[0];
		  t2 = parvals[parvals.size()-1];
		}
	      else if (parvals.size() == 1)
		{
		  Point pt1 = scurve->point(t1);
		  Point pt2 = scurve->point(t2);
		  double u1, u2, v1, v2, d1, d2;
		  Point cl1, cl2;
		  surf->closestPoint(pt1, u1, v1, cl1, d1, int_tol);
		  surf->closestPoint(pt2, u2, v2, cl2, d2, int_tol);
		  if (d1 <= d2)
		    t2 = parvals[0];
		  else
		    t1 = parvals[0];
		}
	      shared_ptr<CurveOnSurface> sub_cv(cvs[ki]->subCurve(t1,t2));
	      cvs[ki] = sub_cv;
	      OK = cvs[ki]->ensureParCrvExistence(tol);
	      if (!OK)
		continue;
	    }
	}
      
      start = std::min(start, cvs[ki]->startparam());
      end = std::max(end, cvs[ki]->endparam());
      endpars[ki] = std::make_pair(start,end);
    }

  // Special treatment of closed curve
  double seam_dist = std::numeric_limits<double>::max();
  if (cvs.size() == 1)
    {
      Point startpt = cvs[0]->ParamCurve::point(cvs[0]->startparam());
      Point endpt = cvs[0]->ParamCurve::point(cvs[0]->endparam());
      if (startpt.dist(endpt) <= tol)
	{
	  closestPoint(startpt, seam_dist);
	}
    }

  // Restrict curve with respect to parameter domain
  double fac = 0.1;
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      vector<double> cvparam(4);
      shared_ptr<ParamCurve> pcurve = cvs[ki]->parameterCurve();
      if (!pcurve.get())
	continue;
      BoundingBox bbp = pcurve->boundingBox();
      double u1 = bbp.low()[0];
      double u2 = bbp.high()[0];
      double v1 = bbp.low()[1];
      double v2 = bbp.high()[1];
      double udel = fac*(domain_[1] - domain_[0]);
      double vdel = fac*(domain_[3] - domain_[2]);
      bool inside = true;
      for (int ka=0; ka<4; ++ka)
	{
	  double dist;
	  Point close;
	  pcurve->closestPoint(corner[ka], pcurve->startparam(),
			       pcurve->endparam(), cvparam[ka], close, dist);
	  if (dist < seam_dist && (u1 < domain_[0]-udel || u2 > domain_[1]+udel ||
				   v1 < domain_[2]-vdel || v2 > domain_[3]+vdel))
	    inside = false;
	}
      std::sort(cvparam.begin(), cvparam.end());
      if (inside)
	{
	  endpars[ki].first = pcurve->startparam();
	  endpars[ki].second = pcurve->endparam();
	}
      else
	{
	  endpars[ki].first = cvparam[0];
	  endpars[ki].second = cvparam[3];
	}
    }


  return true;
}

//===========================================================================
void RevEngRegion::checkReplaceSurf(Point mainaxis[3], int min_pt_reg,
				    double tol, double angtol, bool always)
//===========================================================================
{
  int sfcode;
  ClassType classtype[2];
  classtype[0] = Class_Unknown;
  classtype[1] = associated_sf_[0]->instanceType(sfcode);
  shared_ptr<ParamSurface> primary;
  if (basesf_.get() && avdist_base_ < tol &&
      num_in_base_ > (int)group_points_.size()/2)
    {
      primary = basesf_;
      classtype[0] = primary->instanceType();
    }

  shared_ptr<SplineCurve> profile;
  Point pt1, pt2;
  bool cyllike = false;
  for (int ka=0; ka<2; ++ka)
    {
      if (classtype[ka] == Class_Unknown)
	continue;
  
      shared_ptr<ParamSurface> updated, updated2;
      if (classtype[ka] == Class_SplineSurface && sfcode == 1)
	updated = computeLinearSwept(tol, profile, pt1, pt2);
      else
	computeSurface(group_points_, mainaxis, tol, angtol, classtype[ka], 
		       updated, updated2, cyllike);

      shared_ptr<ParamSurface> replacesurf;
      if (updated.get())
	{
	  double maxd, avd;
	  int num_inside, num2_inside;
	  vector<RevEngPoint*> inpt, outpt;
	  vector<pair<double, double> > dist_ang;
	  vector<double> parvals;
	  RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				  updated, tol, maxd, avd, num_inside, num2_inside,
				  inpt, outpt, parvals, dist_ang, angtol);
	  int sf_flag1 = defineSfFlag(0, tol, num_inside, num2_inside, avd,
				      cyllike);
	  if (updated2.get())
	    {
	      double maxd2, avd2;
	      int num_inside2, num2_inside2;
	      vector<RevEngPoint*> inpt2, outpt2;
	      vector<pair<double, double> > dist_ang2;
	      vector<double> parvals2;
	      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
				      updated2, tol, maxd2, avd2, num_inside2,
				      num2_inside2, inpt2, outpt2, parvals2,
				      dist_ang2, angtol);
	      int sf_flag2 = defineSfFlag(0, tol, num_inside2, num2_inside2, avd2,
					  cyllike);
	      if ((sf_flag2 < sf_flag1 ||
		   (num_inside2 > num_inside ||
		    (num_inside2 == num_inside && avd2 < avd))) &&
		  (sf_flag2 < surfflag_ ||
		    (((num_inside2 > num_inside_ ||
		       (num_inside2 == num_inside_ && avd2 < avdist_))
		      && avd2 < tol)) || always))
		{
		  replacesurf = updated2;
		  for (size_t kh=0; kh<group_points_.size(); ++kh)
		    {
		      group_points_[kh]->setPar(Vector2D(parvals2[2*kh],parvals2[2*kh+1]));
		      group_points_[kh]->setSurfaceDist(dist_ang2[kh].first, dist_ang2[kh].second);
		    }
		  setAccuracy(maxd2, avd2, num_inside2, num2_inside2);
		  setSurfaceFlag(sf_flag2);
		}
	      if (ka == 0 && num_inside2 >= num_in_base_ && avd2 < avdist_base_)
		{
		  setBaseSf(updated2, maxd2, avd2, num_inside2, num2_inside2);
		}
	      int stop_break2 = 1;
	    }
	  if ((!replacesurf.get()) &&
	      (sf_flag1 < surfflag_ ||
	       (num_inside > num_inside_ ||
		(num_inside == num_inside_ && avd < avdist_)) || always))
	    {
	      replacesurf = updated;
	      for (size_t kh=0; kh<group_points_.size(); ++kh)
		{
		  group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
		  group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
		}
	      setAccuracy(maxd, avd, num_inside, num2_inside);
	      setSurfaceFlag(sf_flag1);
	    }
	  if (ka == 0 && num_inside >= num_in_base_ && avd < avdist_base_)
	    {
	      setBaseSf(updated, maxd, avd, num_inside, num2_inside);
	    }
	  int stop_break = 1;
	}

      if (replacesurf.get())
	{
	  associated_sf_[0]->replaceSurf(replacesurf);
	  if (!replacesurf->isBounded())
	    {
	      double diag = bbox_.low().dist(bbox_.high());
	      associated_sf_[0]->limitSurf(2*diag);
	    }
	  
	  surf_adaption_ = INITIAL;
	  if (sfcode == 1 && profile.get())
	    {
	      // A linear swept surface
	      associated_sf_[0]->setLinearSweepInfo(profile, pt1, pt2);
	    }
	  computeDomain();
	  for (size_t kj=0; kj<rev_edges_.size(); ++kj)
	    rev_edges_[kj]->replaceSurf(this, replacesurf, tol);
	}
    }
}


//===========================================================================
int RevEngRegion::checkSurfaceAccuracy(vector<shared_ptr<ElementarySurface> >& sfs,
				       double tol, double angtol, double& maxd,
				       double& avd, int& num_in,
				       int& num2_in, int& sf_flag)
//===========================================================================
{
  vector<double> maxdist(sfs.size());
  vector<double> avdist(sfs.size());
  vector<int> num_inside(sfs.size());
  vector<int> num2_inside(sfs.size());
  vector<int> surf_flag(sfs.size());

  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      bool cyllike = (sfs[ki]->instanceType() == Class_Cylinder ||
		      sfs[ki]->instanceType() == Class_Cone);
      vector<RevEngPoint*> inpt, outpt;
      vector<pair<double, double> > dist_ang;
      vector<double> parvals;
      RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			      sfs[ki], tol, maxdist[ki], avdist[ki],
			      num_inside[ki], num2_inside[ki],
			      inpt, outpt, parvals, dist_ang, angtol);
      surf_flag[ki] = defineSfFlag(0, tol, num_inside[ki], num2_inside[ki], 
				   avdist[ki], cyllike);
    }

  int ix = 0;
  double fac = 1.2;
  for (size_t ki=1; ki<sfs.size(); ++ki)
    {
      if (surf_flag[ki] < surf_flag[ix] ||
	  (surf_flag[ki] == surf_flag[ix] &&
	   ((num2_inside[ki] >= num2_inside[ix] && avdist[ki] <= avdist[ix]) ||
	    (double)num2_inside[ki] >= fac*(double)num2_inside[ix] ||
	    fac*avdist[ki] <= avdist[ix])))
	ix = (int)ki;
    }
  maxd = maxdist[ix];
  avd = avdist[ix];
  num_in = num_inside[ix];
  num2_in = num2_inside[ix];
  sf_flag = surf_flag[ix];
  
  return (sf_flag < NOT_SET) ? ix : -1;
}

//===========================================================================
void RevEngRegion::parameterizePoints(double tol, double angtol)
//===========================================================================
{
  if (!hasSurface())
    return;

    shared_ptr<ParamSurface> surf = getSurface(0)->surface();

    // Compute accuracy
    vector<RevEngPoint*> inpt, outpt;
    vector<pair<double, double> > dist_ang;
    double maxd, avd;
    int num_in, num2_in;
    vector<double> parvals;
    RevEngUtils::distToSurf(group_points_.begin(), group_points_.end(),
			    surf, tol, maxd, avd, num_in, num2_in,
			    inpt, outpt, parvals, dist_ang,
			    angtol);
    
    // Update info in points
    for (size_t kh=0; kh<group_points_.size(); ++kh)
      {
	group_points_[kh]->setPar(Vector2D(parvals[2*kh],parvals[2*kh+1]));
	group_points_[kh]->setSurfaceDist(dist_ang[kh].first, dist_ang[kh].second);
      }
  
    updateInfo(tol, angtol);

    bool cyllike = (surf->instanceType() == Class_Cylinder ||
		    surf->instanceType() == Class_Cone);
    int sf_flag = defineSfFlag(0, tol, num_inside_, num_inside2_,
			       avdist_, cyllike);
    setSurfaceFlag(sf_flag);
}

//===========================================================================
bool RevEngRegion::isCompatible(ClassType classtype, int sfcode)
//===========================================================================
{
  return true;
  double anglim = 0.1;
  if (classtype == Class_Plane)
    {
      if (planartype() || (normalcone_.angle() <= anglim &&
			   (!normalcone_.greaterThanPi())))
	return true;
      else
	return false;
    }
  else
  if (classtype == Class_Cylinder)
    {
      return true;  // Must mike proper criteria
    }
    return false;  // To be extended
}


//===========================================================================
void  RevEngRegion::computeSurface(vector<RevEngPoint*>& points,
				   Point mainaxis[3], double tol,
				   double angtol, ClassType classtype,
				   shared_ptr<ParamSurface>& updated,
				   shared_ptr<ParamSurface>& updated2,
				   bool& cyllike)
//===========================================================================
{
  if (classtype == Class_Plane)
    {
      updated = computePlane(points, avnorm_, mainaxis);
    }
  else if (classtype == Class_Cylinder)
    {
      cyllike = true;
      updated = computeCylinder(points, tol);
    }
  else if (classtype == Class_Sphere)
    {
      Point axis;
      if (hasSurface())
	{
	  shared_ptr<ElementarySurface> elem =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(associated_sf_[0]->surface());
	  if (elem.get())
	    axis = elem->direction();
	}
      updated = computeSphere(mainaxis, axis, points);
    }
  else if (classtype == Class_Cone)
    {
      cyllike = true;
      Point apex;
      updated = computeCone(points, apex);
    }
  else if (classtype == Class_Torus)
    {
      //shared_ptr<Torus> torus2;
      vector<Point> dummy_axis;
      updated = computeTorus(points, dummy_axis, tol, angtol);
      //updated2 = torus2;
    }
  else if (classtype == Class_SplineSurface)
    {
      updated = updateFreeform(points, tol);
      if (!updated.get())
	updated = computeFreeform(points, tol);
    }
#ifdef DEBUG_UPDATE
  if (updated.get())
    {
      std::ofstream ofn("updated.g2");
      updated->writeStandardHeader(ofn);
      updated->write(ofn);
      if (updated2.get())
	{
	  updated2->writeStandardHeader(ofn);
	  updated2->write(ofn);
	}
    }
#endif
}


//===========================================================================
void  RevEngRegion::splitCylinderRad(const Point& pos, const Point& axis,
				     const Point& Cx, const Point& Cy,
				     int nmb_split, vector<Point>& centr,
				     vector<double>& rad)
//===========================================================================
{
  centr.resize(nmb_split);
  rad.resize(nmb_split);
  shared_ptr<Line> line(new Line(pos, axis));
  vector<double> par(group_points_.size());
  double tmin = std::numeric_limits<double>::max();
  double tmax = std::numeric_limits<double>::lowest();
  double diag = bbox_.low().dist(bbox_.high());
  double tmin0 = -diag;
  double tmax0 = diag;
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D xyz = group_points_[ki]->getPoint();
      Point curr(xyz[0], xyz[1], xyz[2]);
      double tpar, dist;
      Point close;
      line->closestPoint(curr, tmin0, tmax0, tpar, close, dist);
      par[ki] = tpar;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
    }

  double tdel = (tmax - tmin)/(double)nmb_split;
  vector<Point> mid(nmb_split);
    double tpar = tmin + 0.5*tdel;
  for (int ka=0; ka<nmb_split; ++ka, tpar+=tdel)
    mid[ka] = line->ParamCurve::point(tpar);

  vector<vector<Point> > proj_pts(nmb_split);
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      Vector3D pnt = group_points_[ki]->getPoint();
      size_t kj;
      for (kj=0, tpar=tmin; kj<mid.size(); ++kj, tpar+=tdel)
	if (par[ki] >= tpar && par[ki] < tpar+tdel)
	  break;
      kj = std::min(kj, mid.size()-1);
      Point curr(pnt[0], pnt[1], pnt[2]);
      Point curr2 = curr - mid[kj];
      curr2 -= ((curr2*axis)*axis);
      curr2 += mid[kj];
      proj_pts[kj].push_back(curr2);
    }

#ifdef DEBUG_SEGMENT
  std::ofstream of("split_circ.g2");
  for (size_t kj=0; kj<proj_pts.size(); ++kj)
    {
      RevEngUtils::computeCircPosRadius(proj_pts[kj], axis, Cx, Cy, centr[kj], rad[kj]);
      shared_ptr<Circle> circ(new Circle(rad[kj], centr[kj], axis, Cx));
      circ->writeStandardHeader(of);
      circ->write(of);
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << proj_pts[kj].size() << std::endl;
      for (size_t kr=0; kr<proj_pts[kj].size(); ++kr)
	of << proj_pts[kj][kr] << std::endl;
    }
#endif

  
  int stop_break = 1;
 }


int compare_t(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

//===========================================================================
void  RevEngRegion::curveApprox(vector<Point>& points, double tol,
				shared_ptr<Circle> circle, vector<double>& parval,
				shared_ptr<SplineCurve>& curve, Point& xpos)
//===========================================================================
{
  double eps = 0.001;
  vector<double> pts;
  vector<double> param;
  double tmin = circle->startparam();
  double tmax = circle->endparam();
  double tdel = tmax - tmin;
  double tmin2 = tmax;
  double tmax2 = tmin;
  vector<double> tmppts;
  tmppts.reserve(4*points.size());
  parval.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      circle->closestPoint(points[ki], tmin, tmax, tpar, close, dist);
      // if (ki > 0 && ((tpar-tmin < tmin2-tpar && tmax-tmax2 < tmin2-tmin) ||
      // 		     (tmax-tpar < tpar-tmax2 && tmin2-tmin < tmax-tmax2)))
      if (ki > 0 && (tpar<tmin2 || tpar>tmax2) &&
	  std::min(fabs(tmin2-tpar+tdel),fabs(tpar+tdel-tmax2)) <
	  std::min(fabs(tpar-tmax2),fabs(tpar-tmin2)))
	//(fabs(tmin2-tpar+tdel) < fabs(tpar-tmax2) || fabs(tpar+tdel-tmax2) < fabs(tpar-tmax2)))
	{
	  if (tpar-tmin < tmax-tpar)
	    tpar += tdel;
	  else
	    tpar -= tdel;
	}

      parval[ki] = tpar;
      tmppts.push_back(tpar);
      tmppts.insert(tmppts.end(), points[ki].begin(), points[ki].end());
      // pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      // param.push_back(tpar);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  qsort(&tmppts[0], points.size(), 4*sizeof(double), compare_t);
  pts.resize(3*points.size());
  param.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      param[ki] = tmppts[4*ki];
      for (size_t ka=0; ka<3; ++ka)
	pts[3*ki+ka] = tmppts[4*ki+ka+1];
    }

  if (tmax2 - tmin2 < 2*M_PI-eps && (tmin2 < -eps || tmax2 > 2*M_PI+eps))
    {
      double tpar = (tmax2 > 2*M_PI+eps) ? 0.5*(tmax2 + tmin2 - 2*M_PI)
		      : 0.5*(tmin2 + tmax2 + 2*M_PI);
      xpos = circle->ParamCurve::point(tpar);
    }
    
  int inner = (int)(2.0*(tmax2 - tmin2)/M_PI);
  int ik = 4;
  int in = ik + inner;
  // double tdel = (tmax2 - tmin2)/(double)(in - ik + 1);
  // double et[12];
  // for (int ka=0; ka<ik; ++ka)
  //   {
  //     et[ka] = tmin2;
  //     et[in+ka] = tmax2;
  //   }
  // for (int ka=ik; ka<in; ++ka)
  //   et[ka] = tmin2 + (ka-ik+1)*tdel;

  double smoothwgt = 1.0e-9; //0.001;
  ApproxCurve approx(pts, param, 3, tol, in, ik);
  approx.setSmooth(smoothwgt);
  int maxiter = 4; //3;
  double maxdist, avdist;
  curve = approx.getApproxCurve(maxdist, avdist, maxiter);
  // vector<double> ecoef(3*in, 0.0);
  // shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, et, &ecoef[0], 3));

  // SmoothCurve smooth(3);
  // vector<int> cfn(in, 0);
  // vector<double> wgts(param.size(), 1.0);
  // smooth.attach(cv, &cfn[0]);

  // double wgt1 = 0.0, wgt2 = 0.1, wgt3 = 0.1;
  // double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  // smooth.setOptim(wgt1, wgt2, wgt3);
  // smooth.setLeastSquares(pts, param, wgts, approxwgt);

  // shared_ptr<SplineCurve> curve0;
  // smooth.equationSolve(curve0);
#ifdef DEBUG
  std::ofstream of("points_and_tangents.txt");
  of << points.size() << std::endl;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      of << pts[3*ki] << " " << pts[3*ki+1] << " ";
      vector<Point> der(2);
      curve->point(der, param[ki], 1);
      of << der[1][0] << " " << der[1][1] << std::endl;
    }
#endif
  int stop_break = 1;
}

//===========================================================================
void RevEngRegion::configSplit(vector<RevEngPoint*>& points,
			       vector<double>& param,
			       shared_ptr<Cylinder> cyl,
			       shared_ptr<SplineCurve> spl, double tol,
			       vector<vector<RevEngPoint*> >& configs)
//===========================================================================
{
  // Check cylinder axis
  Point axis = cyl->direction();
  double angtol = 0.2;
  double inlim = 0.8;
  int nmb_in = 0;
  double mpi2 = 0.5*M_PI;
  vector<Vector3D> low, high;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Point normal = group_points_[kr]->getLocFuncNormal();
      double ang = axis.angle(normal);
      if (ang < mpi2)
	low.push_back(group_points_[kr]->getPoint());
      else
	high.push_back(group_points_[kr]->getPoint());
      if (fabs(ang-mpi2) <= angtol)
	nmb_in++;
    }
  double in_frac = (double)nmb_in/(double)group_points_.size();
  if (in_frac < inlim)
    return;

#ifdef DEBUG_SEGMENT
  std::ofstream of1("low_axis.g2");
  of1 << "400 1 0 4 100 155 0 255" << std::endl;
  of1 << low.size() << std::endl;
  for (size_t kr=0; kr<low.size(); ++kr)
    of1 << low[kr] << std::endl;
  
  std::ofstream of2("high_axis.g2");
  of2 << "400 1 0 4 0 155 100 255" << std::endl;
  of2 << high.size() << std::endl;
  for (size_t kr=0; kr<high.size(); ++kr)
    of2 << high[kr] << std::endl;
#endif
  
  // Compute all intersections between the cylinder and the spline curve
  double eps = std::min(1.0e-6, 0.1*tol);
  vector<double> intpar;
  vector<pair<double,double> > int_cvs;
  intersectCurveCylinder(spl.get(), cyl->getLocation(), cyl->getAxis(),
			 cyl->getRadius(), eps, intpar, int_cvs);
  for (size_t ki=0; ki<int_cvs.size(); ++ki)
    {
      intpar.push_back(int_cvs[ki].first);
      intpar.push_back(int_cvs[ki].second);
    }
  if (intpar.size() == 0)
    return;

  // Define parameter intervals of different configurations
  std::sort(intpar.begin(), intpar.end());
  vector<double> delpar;

  // Check startpoint, endpoint and points in the middle of intervals
  double upar, vpar, dist;
  Point close;
  Point pos = spl->ParamCurve::point(spl->startparam());
  cyl->closestPoint(pos, upar, vpar, close, dist, eps);
  delpar.push_back(spl->startparam()-tol);
  if (dist > tol)
    {
      delpar.push_back(intpar[0]);
    }
  for (size_t ki=1; ki<intpar.size(); ++ki)
    {
      double tpar = 0.5*(intpar[ki-1] + intpar[ki]);
      pos = spl->ParamCurve::point(tpar);
      cyl->closestPoint(pos, upar, vpar, close, dist, eps);
      if (dist > tol)
	delpar.push_back(intpar[ki]);
    }
  pos = spl->ParamCurve::point(spl->endparam());
  cyl->closestPoint(pos, upar, vpar, close, dist, eps);
  if (dist > tol)
    {
      if (delpar.size() > 0 &&
	  intpar[intpar.size()-1] > delpar[delpar.size()-1])
	delpar.push_back(intpar[intpar.size()-1]);
    }
  delpar.push_back(spl->endparam());

  if (delpar.size() == 0)
    return;
  
  // Divide point set according to configuration
  configs.resize(delpar.size()-1);
  for (size_t kj=0; kj<points.size(); ++kj)
    {
      for (size_t ki=1; ki<delpar.size(); ++ki)
	if (param[kj] > delpar[ki-1] && param[kj] <= delpar[ki])
	  {
	    configs[ki-1].push_back(points[kj]);
	    break;
	  }
    }
  int stop_break = 1;
}

//===========================================================================
bool  RevEngRegion::hasEdgeBetween(RevEngRegion* adj)
//===========================================================================
{
  if (adj->numPoints() < (int)group_points_.size())
    return adj->hasEdgeBetween(this);

  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      vector<ftSamplePoint*> next = group_points_[ki]->getNeighbours();
      for (size_t kj=0; kj<next.size(); ++kj)
	{
	  RevEngPoint *pt = dynamic_cast<RevEngPoint*>(next[kj]);
	  if (!pt->isEdge(edge_class_type_))
	    continue;

	  vector<ftSamplePoint*> next2 = pt->getNeighbours();
	  for (size_t kr=0; kr<next2.size(); ++kr)
	    {
	      RevEngPoint *pt2 = dynamic_cast<RevEngPoint*>(next2[kr]);
	      if (pt2->region() == adj)
		return true;
	    }
	}
    }
  return false;
}

//===========================================================================
void RevEngRegion::store(std::ostream& os) const
//===========================================================================
{
  os << Id_ << std::endl;
  os << group_points_.size() << std::endl;
for (size_t ki=0; ki<group_points_.size(); ++ki)
  os << group_points_[ki]->getIndex() << " ";
 os << std::endl;
 os << classification_type_ << " " << surfflag_ << " " << surf_adaption_;
 os << " " << frac_norm_in_ << " " << frac_norm_in2_ << std::endl;
 os << maxdist_ << " " << avdist_ << " " << num_inside_;
 os << " " << num_inside2_ << " " << std::endl;
 int base = basesf_.get() ? 1 : 0;
 os << base << std::endl;
 if (base)
   {
     basesf_->writeStandardHeader(os);
     basesf_->write(os);
     os << maxdist_base_ << " " << avdist_base_ << " " << num_in_base_;
     os << " " << num_in_base2_ << std::endl;
   }

 os << associated_sf_.size() << std::endl;
 for (size_t ki=0; ki<associated_sf_.size(); ++ki)
     os << associated_sf_[ki]->getId() << " ";

 int sweep = sweep_.get() ? 1 : 0;
 os << sweep << std::endl;
 if (sweep_)
   {
     os << sweep_->type_ << std::endl;
     if (sweep_->profile_.get())
       {
	 os << "1" << std::endl;
	 sweep_->profile_->writeStandardHeader(os);
	 sweep_->profile_->write(os);
       }
     else
       os << "0" << std::endl;
     os << sweep_->location_ << " " << sweep_->added_info_ << " ";
     os << sweep_->radius_ << " " << sweep_->angle_ << std::endl;
     os << sweep_->maxdist_ << " " << sweep_->avdist_ << " " << sweep_->num_in_<< std::endl;
   }
}

//===========================================================================
void RevEngRegion::read(std::istream& is,
			shared_ptr<ftPointSet>& tri_sf,
			vector<int>& associated_sf_id)
//===========================================================================
{
  GoTools::init();
  is >> Id_;
  int num_points;
  is >> num_points;
  group_points_.resize(num_points);
  int ix;
  for (int ki=0; ki<num_points; ++ki)
    {
      is >> ix;
      RevEngPoint* pt = dynamic_cast<RevEngPoint*>((*tri_sf)[ix]);
      pt->setRegion(this);
      group_points_[ki] = pt;
      }
  is >> classification_type_ >> surfflag_ >> surf_adaption_;
  is >> frac_norm_in_ >> frac_norm_in2_;
  is >> maxdist_ >> avdist_ >> num_inside_ >> num_inside2_;
  int base;
  is >> base;

  if (base)
    {
      ObjectHeader header;
      header.read(is);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(is);
      basesf_ = dynamic_pointer_cast<ParamSurface,GeomObject>(obj);
      is >> maxdist_base_ >> avdist_base_ >> num_in_base_ >> num_in_base2_;
    }

  int num_sf;
  is >> num_sf;
  for (int ki=0; ki<num_sf; ++ki)
    {
      int sf_id;
      is >> sf_id;
      associated_sf_id.push_back(sf_id);
    }
  if (num_sf > 0)
    computeDomain();

  int sweep;
  is >> sweep;
  if (sweep)
    {
      int stype, snum_in, sprof;
      double rad, ang, maxd, avd;
      Point loc(3), added(3);
      shared_ptr<SplineCurve> profile;
      is >> stype;
      is >> sprof;
      if (sprof)
	{
	  ObjectHeader header;
	  header.read(is);
	  profile = shared_ptr<SplineCurve>(new SplineCurve());
	  profile->read(is);
	}
      is >> loc >> added >> rad >> ang >> maxd >> avd >> snum_in;
      sweep_ = shared_ptr<SweepData>(new SweepData(stype, profile, loc, added, maxd,
						   avd, snum_in, rad, ang));
    }
	
  // Bounding box and principal curvature summary
  maxk2_ = std::numeric_limits<double>::lowest();
  mink2_ = std::numeric_limits<double>::max();
  maxk1_ = std::numeric_limits<double>::lowest();
  mink1_ = std::numeric_limits<double>::max();
  MAH_ = MAK_ = avH_ = avK_ = 0.0;
  bbox_ = BoundingBox(3);
  if (group_points_.size() > 0)
    {
      double fac = 1.0/(double)group_points_.size();
      for  (size_t kj=0; kj<group_points_.size(); ++kj)
	{
	  double k1 = group_points_[kj]->minPrincipalCurvature();
	  double k2 = group_points_[kj]->maxPrincipalCurvature();
	  double H = group_points_[kj]->meanCurvature();
	  double K = group_points_[kj]->GaussCurvature();
	  mink1_ = std::min(mink1_, fabs(k1));
	  maxk1_ = std::max(maxk1_, fabs(k1));
	  mink2_ = std::min(mink2_, fabs(k2));
	  maxk2_ = std::max(maxk2_, fabs(k2));
	  avH_ += fac*H;
	  avK_ += fac*K;
	  MAH_ += fac*fabs(H);
	  MAK_ += fac*fabs(K);
	  Vector3D point = group_points_[kj]->getPoint();
	  Point point2(point[0], point[1], point[2]);
	  bbox_.addUnionWith(point2);
	}
  
      normalcone_ = DirectionCone(group_points_[0]->getLocFuncNormal());
      normalcone2_ = DirectionCone(group_points_[0]->getTriangNormal());
      avnorm_ = Point(0.0, 0.0, 0.0);
      avnorm2_ = Point(0.0, 0.0, 0.0);
      for  (size_t kj=1; kj<group_points_.size(); ++kj)
	{
	  Point norm = group_points_[kj]->getLocFuncNormal();
	  normalcone_.addUnionWith(norm);
	  avnorm_ += fac*norm;
	  Point norm2 = group_points_[kj]->getTriangNormal();
	  normalcone2_.addUnionWith(norm2);
	  avnorm2_ += fac*norm2;
	}
    }
}

//===========================================================================
void RevEngRegion::getRemainingPoints(vector<RevEngPoint*>& curr_pts,
				      vector<RevEngPoint*>& remaining)
//===========================================================================
{
  remaining.clear();
  for (size_t ki=0; ki<group_points_.size(); ++ki)
    {
      auto it = std::find(curr_pts.begin(), curr_pts.end(), group_points_[ki]);
      if (it == curr_pts.end())
	remaining.push_back(group_points_[ki]);
    }
}

//===========================================================================

void RevEngRegion::writeSubTriangulation(std::ostream& of)
{
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      vector<ftSamplePoint*> next = group_points_[kr]->getNeighbours();
      size_t nmb_next = next.size();
      for (int ka=(int)(next.size()-1); ka>=0; --ka)
	{
	  RevEngPoint *pr = dynamic_cast<RevEngPoint*>(next[ka]);
	  if (pr->region() != this)
	    next.erase(next.begin()+ka);
	}

      int bd = (next.size() < nmb_next) ? 1 : 0;
      of << kr << " " << group_points_[kr]->getPoint() << " " << bd << std::endl;
      of << next.size() << " ";
      for (size_t kh=0; kh<next.size(); ++kh)
	{
	  vector<RevEngPoint*>::iterator it = std::find(group_points_.begin(),
							group_points_.end(), next[kh]);
	  if (it == group_points_.end())
	    {
#ifdef DEBUG
	      std::cout << "writeSubTriangulation: missing connection" << std::endl;
#endif
	    }
	  else
	    {
	      size_t ix = it - group_points_.begin();
	      of << ix << " ";
	    }
	}
      of << std::endl;
    }
}

void RevEngRegion::writeSurface(std::ostream& of)
{
  if (associated_sf_.size() == 0)
    return;
  shared_ptr<ParamSurface> surf = associated_sf_[0]->surface();
  shared_ptr<ElementarySurface> elemsf =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (elemsf.get())
    {
      // double umin = std::numeric_limits<double>::max();
      // double umax = std::numeric_limits<double>::lowest();
      // double vmin = std::numeric_limits<double>::max();
      // double vmax = std::numeric_limits<double>::lowest();
      // for (size_t ki=0; ki<group_points_.size(); ++ki)
      // 	{
      // 	  Vector2D par = group_points_[ki]->getPar();
      // 	  umin = std::min(umin, par[0]);
      // 	  umax = std::max(umax, par[0]);
      // 	  vmin = std::min(vmin, par[1]);
      // 	  vmax = std::max(vmax, par[1]);
      // 	}
      double umin = domain_[0];
      double umax = domain_[1];
      double vmin = domain_[2];
      double vmax = domain_[3];
      shared_ptr<ElementarySurface> elemsf2(elemsf->clone());
      if (elemsf2->instanceType() != Class_Plane && umax-umin > 2*M_PI)
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      if (elemsf2->instanceType() != Class_Plane && elemsf2->instanceType() != Class_Sphere &&
	  (umin < -2.0*M_PI || umax > 2.0*M_PI))
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      if (elemsf2->instanceType() == Class_Sphere && (umin < 0.0 || umax > 2.0*M_PI))
	{
	  umin = 0;
	  umax = 2*M_PI;
	}
      
     if (elemsf2->instanceType() == Class_Sphere &&
	 (vmax-vmin > M_PI || vmin < -0.5*M_PI || vmax > 0.5*M_PI))
       {
	 vmin = -0.5*M_PI;
	 vmax = 0.5*M_PI;
       }
      if (elemsf2->instanceType() == Class_Torus &&
	  (vmax - vmin > 2*M_PI || vmin < -2.0*M_PI || vmax > 2.0*M_PI))
       {
	 vmin = 0;
	 vmax = 2*M_PI;
       }

      if (elemsf2->isBounded())
	{
	  RectDomain dom = elemsf2->getParameterBounds();
	  double udel = umax - umin;
	  double vdel = vmax - vmin;
	  if (dom.umax() - dom.umin() < 2.0*udel)
	    {
	      umin = dom.umin();
	      umax = dom.umax();
	    }
	  if (dom.vmax() - dom.vmin() < 2.0*vdel)
	    {
	      vmin = dom.vmin();
	      vmax = dom.vmax();
	    }
	}
      if (umax > umin && vmax > vmin)
	elemsf2->setParameterBounds(umin, vmin, umax, vmax);
      elemsf2->writeStandardHeader(of);
      elemsf2->write(of);
    }
  else
    {
      surf->writeStandardHeader(of);
      surf->write(of);
    }
}

void RevEngRegion::writeRegionInfo(std::ostream& of)
{
  double len = bbox_.low().dist(bbox_.high());
  double ll = 0.05*len;
  of << "400 1 0 4 100 0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
  of << "410 1 0 4 200 55 0 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point norm = group_points_[kr]->getLocFuncNormal();
      of << xyz2 << " " << xyz2 + ll*norm << std::endl;
    }
  
  of << "410 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      Vector3D xyz = group_points_[kr]->getPoint();
      Point xyz2(xyz[0], xyz[1], xyz[2]);
      Point vec = group_points_[kr]->getTriangNormal(); //minCurvatureVec();
      of << xyz2 << " " << xyz2 + ll*vec << std::endl;
    }

  // double lambda[2];
  // Point eigen1, eigen2, eigen3;
  // getPCA(lambda, eigen1, eigen2, eigen3);
  // std::ofstream of2("points_two_tangents.txt");
  // of2 << group_points_.size() << std::endl;
  // for (size_t kr=0; kr<group_points_.size(); ++kr)
  //   {
  //     Vector3D xyz = group_points_[kr]->getPoint();
  //     Point norm = group_points_[kr]->getLocFuncNormal();
  //     Point vec1 = eigen1 - (eigen1*norm)*norm;
  //     vec1.normalize();
  //     Point vec2 = norm.cross(vec1);
  //     vec2.normalize();
  //     of2 << xyz << " " << vec1 << " " << vec2 << std::endl;
  //   }
}

void RevEngRegion::writeRegionPoints(std::ostream& of)
{
  of << "400 1 0 4 100 0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    of << group_points_[kr]->getPoint() << std::endl;
}  

void RevEngRegion::writeAdjacentPoints(std::ostream& of)
{
  for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
    {
      (*it)->writeRegionPoints(of);
      if ((*it)->hasSurface())
	{
	  shared_ptr<ParamSurface> tmp((*it)->getSurface(0)->surface()->clone());
	  if (!tmp->isBounded())
	    {
	      shared_ptr<ElementarySurface> elem =
		dynamic_pointer_cast<ElementarySurface,ParamSurface>(tmp);
	      double dom[4];
	      (*it)->getDomain(dom);
	      try {
	      if (elem.get())
		elem->setParameterBounds(dom[0], dom[2], dom[1], dom[3]);
	      }
	      catch (...)
		{
		}
	    }
	  tmp->writeStandardHeader(of);
	  tmp->write(of);
	}
    }
}

void RevEngRegion::writeUnitSphereInfo(std::ostream& of)
{
  of << "400 1 0 4 100  0 155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point norm = pt->getLocFuncNormal();
      of << norm << std::endl;
    }
  Sphere sph(1.0, Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
	     Point(1.0, 0.0, 0.0));
  sph.writeStandardHeader(of);
  sph.write(of);

  of << "400 1 0 4 0 100  155 255" << std::endl;
  of << group_points_.size() << std::endl;
  for (size_t kr=0; kr<group_points_.size(); ++kr)
    {
      RevEngPoint *pt = dynamic_cast<RevEngPoint*>(group_points_[kr]);
      Point vec = pt->minCurvatureVec();
      of << vec << std::endl;
    }
  of << std::endl;
}
