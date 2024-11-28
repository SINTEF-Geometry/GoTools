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

#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/creators/TrimCrvUtils.h"
#include "GoTools/creators/TrimUtils.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

//#define DEBUG_TRIM

using namespace Go;
using std::vector;

// //===========================================================================
// HedgeSurface::HedgeSurface()
//   : ftSurface()
// //===========================================================================
// {
// }

//===========================================================================
HedgeSurface::HedgeSurface()
  : ftSurface()
//===========================================================================
{
}

//===========================================================================
HedgeSurface::HedgeSurface(shared_ptr<ParamSurface> sf, RevEngRegion *region)
  : ftSurface(sf, -1)
//===========================================================================
{
  regions_.push_back(region);
  bbox_ = region->boundingBox();
  surf_code_ = SURF_TYPE_UNDEF;
}

//===========================================================================
HedgeSurface::HedgeSurface(shared_ptr<ParamSurface> sf,
			   vector<RevEngRegion*>& region)
  : ftSurface(sf, -1), regions_(region)
//===========================================================================
{
  if (region.size() > 0)
    {
      bbox_ = region[0]->boundingBox();
      for (size_t ki=1; ki<region.size(); ++ki)
	bbox_.addUnionWith(region[ki]->boundingBox());
    }
  surf_code_ = SURF_TYPE_UNDEF;
}


//===========================================================================
HedgeSurface::~HedgeSurface()
//===========================================================================
{
}

//===========================================================================
int HedgeSurface::numPoints()
//===========================================================================
{
  int num = 0;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    num += regions_[ki]->numPoints();
  return num;
 }

//===========================================================================
ClassType HedgeSurface::instanceType(int& code)
//===========================================================================
{
  code = surf_code_;
  ClassType type = surface()->instanceType();
  if (type == Class_BoundedSurface)
    {
      shared_ptr<ParamSurface> surf = surface();
      shared_ptr<BoundedSurface> bdsurf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
      type = bdsurf->underlyingSurface()->instanceType();
    }
  return type;
}


//===========================================================================
bool HedgeSurface::updateSurfaceWithAxis(Point axis[3], int ix, double tol,
					 double angtol)
//===========================================================================
{
  bool updated = false;
  int code = -1;
  ClassType type = instanceType(code);
  if (type == Class_Plane)
    updated = updatePlaneWithAxis(axis, ix, tol, angtol);
  else if (type == Class_Cylinder)
    updated = updateCylinderWithAxis(axis, ix, tol, angtol);

  return updated;
}

//===========================================================================
bool HedgeSurface::updatePlaneWithAxis(Point axis[3], int ix, double tol,
				       double angtol)
//===========================================================================
{
  shared_ptr<Plane> init_plane = dynamic_pointer_cast<Plane,ParamSurface>(surf_);
  if (!init_plane.get())
    return false;
  Point loc = init_plane->location();
  Point mid(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)numPoints();
  for (size_t ki=0; ki<regions_.size(); ++ki)
    for (auto it=regions_[ki]->pointsBegin(); it!=regions_[ki]->pointsEnd(); ++it)
      {
	Vector3D pos0 = (*it)->getPoint();
	Point pos(pos0[0], pos0[1], pos0[2]);
	Point vec = pos - loc;
	Point pos2 = loc + (vec*axis[ix])*axis[ix];
	mid += wgt*pos2;
      }

  shared_ptr<Plane> plane(new Plane(mid, axis[ix], axis[(ix+1)%3]));

  // Check accuracy
  bool updated = checkAccuracyAndUpdate(plane, tol, angtol);
  return updated;
}

//===========================================================================
bool HedgeSurface::updateCylinderWithAxis(Point axis[3], int ix, double tol,
					  double angtol)
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points;
  for (size_t ki=0; ki<regions_.size(); ++ki)
    points.push_back(std::make_pair(regions_[ki]->pointsBegin(),
				    regions_[ki]->pointsEnd()));
  
  Point pos;
  double rad;
  Point low = bbox_.low();
  Point high = bbox_.high();
  RevEngUtils::computeCylPosRadius(points, low, high, axis[ix], axis[(ix+1)%3],
				   axis[(ix+2)%3], pos, rad);

  // Select x-axis
  shared_ptr<ElementarySurface> elem =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf_);
  int ix2 = (ix+2)%3;
  if (elem.get())
    {
      Point dir = elem->direction2();
      double ang1 = dir.angle(axis[(ix+1)%3]);
      ang1 = std::min(ang1, M_PI-ang1);
      double ang2 = dir.angle(axis[(ix+2)%3]);
      ang2 = std::min(ang1, M_PI-ang2);
      if (ang1 < ang2)
	ix2 = (ix+1)%3;
    }
  shared_ptr<Cylinder> cyl(new Cylinder(rad, pos, axis[ix], axis[ix2]));
  
  // Check accuracy
  bool updated = checkAccuracyAndUpdate(cyl, tol, angtol);
  return updated;
}

//===========================================================================
bool HedgeSurface::checkAccuracyAndUpdate(shared_ptr<ParamSurface> surf,
					  double tol, double angtol)
//===========================================================================
{
  bool updated = false;
  vector<vector<pair<double, double> > > dist_ang(regions_.size());
  vector<double> maxd(regions_.size()), avd(regions_.size());
  vector<int> num_in(regions_.size());
  vector<int> num2_in(regions_.size());
  vector<vector<double> > parvals(regions_.size());
  int all_in = 0;
  double avd_all = 0.0;
  double fac = 1.0/(double)numPoints();
  vector<int> sfflag(regions_.size(), NOT_SET);
  bool sfOK = true;
  bool cyllike = (surf->instanceType() == Class_Cylinder ||
		  surf->instanceType() == Class_Cone);
  for (size_t ki=0; ki<regions_.size(); ++ki)
    {
      vector<RevEngPoint*> inpt, outpt;
      RevEngUtils::distToSurf(regions_[ki]->pointsBegin(), regions_[ki]->pointsEnd(),
			      surf, tol, maxd[ki], avd[ki], num_in[ki], num2_in[ki],
			      inpt, outpt, parvals[ki], dist_ang[ki], angtol);
      all_in += num_in[ki];
      avd_all += fac*regions_[ki]->numPoints()*avd[ki];
      sfflag[ki] = regions_[ki]->defineSfFlag(0, tol, num_in[ki], num2_in[ki],
					      avd[ki], cyllike);
      if (sfflag[ki] == NOT_SET)
	sfOK = false;
    }

  if (sfOK) //all_in > numPoints()/2 && avd_all <= tol)
    {
      updated = true;
      for (size_t ki=0; ki<regions_.size(); ++ki)
	{
	  size_t kj=0;
	  for (auto it=regions_[ki]->pointsBegin(); it!=regions_[ki]->pointsEnd();
	       ++it, ++kj)
	    {
	      (*it)->setPar(Vector2D(parvals[ki][2*kj],parvals[ki][2*kj+1]));
	      (*it)->setSurfaceDist(dist_ang[ki][kj].first, dist_ang[ki][kj].second);
	    }
	  regions_[ki]->setAccuracy(maxd[ki], avd[ki], num_in[ki], num2_in[ki]);
	  regions_[ki]->computeDomain();
	  regions_[ki]->setSurfaceFlag(sfflag[ki]);
	}

      replaceSurf(surf);
    }
  return updated;
}

//===========================================================================
bool HedgeSurface::isCompatible(HedgeSurface* other, double angtol, double approx_tol, ClassType& type, double& score)
//===========================================================================
{
  score = std::numeric_limits<double>::max();
  int code1 = -1, code2 = -1;
  ClassType type1 = instanceType(code1);
  ClassType type2 = other->instanceType(code2);

  shared_ptr<ParamSurface> surf1 = surface();
  shared_ptr<ElementarySurface> psurf1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  if (!psurf1.get())
    {
      shared_ptr<BoundedSurface> bdsf1 =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
      if (bdsf1.get())
	{
	  surf1 = bdsf1->underlyingSurface();
	  psurf1 = dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
	}
    }
  
  shared_ptr<ParamSurface> surf2 = other->surface();
  shared_ptr<ElementarySurface> psurf2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
  if (!psurf2.get())
    {
      shared_ptr<BoundedSurface> bdsf2 =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
      if (bdsf2.get())
	{
	  surf2 = bdsf2->underlyingSurface();
	  psurf2 = dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);
	}
    }

  if (!psurf1.get())
    {
      RevEngRegion *preg = 0;
      int numpt = 0;
      int nreg = numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  getRegion(ka);
	  if (reg->hasBaseSf())
	    {
	      double maxdp, avdp;
	      int num_inp, num2_inp;
	      int curr_numpt = reg->numPoints();
	      reg->getBaseDist(maxdp, avdp, num_inp, num2_inp);
	      if (avdp < approx_tol && num_inp > curr_numpt/2 && curr_numpt > numpt)
		{
		  preg = reg;
		  numpt = curr_numpt;
		}
	    }
	}
      if (preg)
	{
	  psurf1 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(preg->getBase());
	  type1 = preg->getBase()->instanceType();
	}
    }
  
  if (!psurf2.get())
    {
      RevEngRegion *preg = 0;
      int numpt = 0;
      int nreg = other->numRegions();
      for (int ka=0; ka<nreg; ++ka)
	{
	  RevEngRegion *reg =  other->getRegion(ka);
	  if (reg->hasBaseSf())
	    {
	      double maxdp, avdp;
	      int num_inp, num2_inp;
	      int curr_numpt = reg->numPoints();
	      reg->getBaseDist(maxdp, avdp, num_inp, num2_inp);
	      if (avdp < approx_tol && num_inp > curr_numpt/2 && curr_numpt > numpt)
		{
		  preg = reg;
		  numpt = curr_numpt;
		}
	    }
	}
      if (preg)
	{
	  psurf2 =
	    dynamic_pointer_cast<ElementarySurface,ParamSurface>(preg->getBase());
	  type2 = preg->getBase()->instanceType();
	}
    }
  
  if (!psurf1.get() || !psurf2.get())
    return false;
  if (type1 != type2 || code1 != code2)
    return false;  // Not the same surface type
  type = type1;

  Point loc1 = psurf1->location();
  Point loc2 = psurf2->location();
  Point vec1 = psurf1->direction();
  Point vec2 = psurf2->direction();
  double rad1 = psurf1->radius(0.0, 0.0);  // Parameters must be set properly for cone
  double rad2 = psurf2->radius(0.0, 0.0);
  double smallrad1 = psurf1->radius2(0.0, 0.0);
  double smallrad2 = psurf2->radius2(0.0, 0.0);
  vec1.normalize_checked();
  vec2.normalize_checked();

  int sgn = (type1 == Class_Plane) ? -1 : 1;
  double ang = vec1.angle(vec2);
  ang = std::min(ang, M_PI-ang);

  double dlim = (rad1 < 0.0) ? approx_tol : std::max(0.05*rad1, approx_tol);
  double anglim = 10.0*angtol;
  double eps = 1.0e-8;
  if (ang > anglim)
    return false;
  // if (fabs(rad2-rad1) > dlim && fabs(smallrad2-smallrad1) < eps)
  //   return false;
  if (rad1 >= 0.0 && rad2 >= 0.0 &&
      fabs(rad2-rad1) > std::min(rad1, rad2) && fabs(smallrad2-smallrad1) < eps)
    return false;
  else if (smallrad1 > 0.0 &&
	   (rad1 < rad2-smallrad2 || rad1 > rad2+smallrad2 ||
	    rad2 < rad1-smallrad1 || rad2 > rad1+smallrad1))
    return false;
    
  double pdist1 = 0.0, pdist2 = 0.0;
  if (type1 == Class_Plane)
    {
      Point loc2_0 = loc2 - ((loc2-loc1)*vec1)*vec1;
      Point loc1_0 = loc1 - ((loc2-loc1)*vec2)*vec2;
      pdist1 = loc2.dist(loc2_0);
      pdist2 = loc1.dist(loc1_0);
      pdist1 = pdist2 = std::min(pdist1, pdist2);
      if (pdist1 > 2.0*dlim || pdist2 > 2.0*dlim)    
	return false;
    }
  else if (type1 == Class_Cylinder)
    {
      Point loc2_0 = loc1 + ((loc2-loc1)*vec1)*vec1;
      pdist1 = loc2.dist(loc2_0);
      if (pdist1 + fabs(rad2-rad1) > std::min(rad1, rad2))
      // if (pdist1 > dlim)
	return false;
    }
  else if (type1 == Class_Torus)
    {
      pdist1 = loc1.dist(loc2);
      if (pdist1 > 2.0*dlim || pdist2 > 2.0*dlim)    
	return false;
    }

  score = 10.0*ang + fabs(rad2-rad1) + fabs(smallrad2-smallrad1) +
    pdist1 + pdist2;
  return true;
}

//===========================================================================
bool HedgeSurface::hasBaseSf()
//===========================================================================
{
  for (size_t ki=0; ki<regions_.size(); ++ki)
    if (regions_[ki]->hasBaseSf())
      return true;
  return false;
}

//===========================================================================
void HedgeSurface::ensureSurfaceBounded()
//===========================================================================
{
  shared_ptr<ParamSurface> surf = surface();
  shared_ptr<ElementarySurface> elemsf =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (!elemsf.get())
    return;  // Always bounded
  if (!elemsf->isBounded())
    {
      double diag = bbox_.low().dist(bbox_.high());
      if (elemsf->instanceType() == Class_Plane)
	elemsf->setParameterBounds(-0.6*diag, -0.6*diag, 0.6*diag, 0.6*diag);
      else if (elemsf->instanceType() == Class_Cylinder)
	{
	  shared_ptr<Cylinder> cyl =
	    dynamic_pointer_cast<Cylinder,ElementarySurface>(elemsf);
	  cyl->setParamBoundsV(-0.6*diag, 0.6*diag);
	}
      else if (elemsf->instanceType() == Class_Cone)
	{
	  shared_ptr<Cone> cone =
	    dynamic_pointer_cast<Cone,ElementarySurface>(elemsf);
	  cone->setParamBoundsV(-0.6*diag, 0.6*diag);
	}
    }
}

//===========================================================================
bool HedgeSurface::isTangential(HedgeSurface* surf)
//===========================================================================
{
  return false;  // To be implemented properly
}
 
//===========================================================================
void HedgeSurface::doTrim(vector<shared_ptr<CurveOnSurface> >& int_cvs,
			  shared_ptr<BoundedSurface>& bdsf,
			  double tol,
			  vector<shared_ptr<HedgeSurface> >& added_sfs)
//===========================================================================
{
  vector<shared_ptr<BoundedSurface> > trim_sfs;
  if (int_cvs.size() > 0)
    {
      trim_sfs = BoundedUtils::splitWithTrimSegments(bdsf, int_cvs, tol);
    }

  std::ofstream of1("curr_trim.g2");
  for (size_t ki=0; ki<trim_sfs.size(); ++ki)
    {
      trim_sfs[ki]->writeStandardHeader(of1);
      trim_sfs[ki]->write(of1);
    }

  if (trim_sfs.size() <= 1)
    {
      // Do something to complement the intersection results
      int stop_break1 = 1;
    }

  // Assume, for the time being, that no region is split between sub surfaces
  vector<shared_ptr<BoundedSurface> > valid_trim;
  vector<vector<RevEngRegion*> > regs;
  for (size_t ki=0; ki<trim_sfs.size(); ++ki)
    {
      vector<RevEngRegion*> curr_regs;
      for (size_t kj=0; kj<regions_.size(); ++kj)
	{
	  int numpt = regions_[kj]->numPoints();
	  int inside = 2;
	  for (int ka=0; ka<numpt; ++ka)
	    {
	      RevEngPoint *pt = regions_[kj]->getPoint(ka);
	      Vector2D parpt = pt->getPar();
	      inside = trim_sfs[ki]->inDomain2(parpt[0], parpt[1]);
	      if (inside != 2)
		break;  // Not on a boundary
	    }
	  if (inside == 1)
	    curr_regs.push_back(regions_[kj]);
	}
      if (curr_regs.size() > 0)
	{
	  valid_trim.push_back(trim_sfs[ki]);
	  regs.push_back(curr_regs);
	}
    }

  for (size_t ki=1; ki<valid_trim.size(); ++ki)
    {
      shared_ptr<HedgeSurface> hedge(new HedgeSurface(valid_trim[ki], regs[ki]));
      for (size_t kj=0; kj<regs[ki].size(); ++kj)
	{
	  (void)removeRegion(regs[ki][kj]);
	  regs[ki][kj]->setHedge(hedge.get());
	}
      added_sfs.push_back(hedge);
    }
  if (valid_trim.size() > 0)
    replaceSurf(valid_trim[0]);
  
  
  int stop_break = 1;
}

//===========================================================================
bool HedgeSurface::removeRegion(RevEngRegion* reg)
//===========================================================================
{
  auto it = std::find(regions_.begin(), regions_.end(), reg);
  if (it != regions_.end())
    {
      regions_.erase(it);
      if (regions_.size() > 0)
	{
	  bbox_ = regions_[0]->boundingBox();
	  for (size_t ki=1; ki<regions_.size(); ++ki)
	    bbox_.addUnionWith(regions_[ki]->boundingBox());
	}
    }
  return false;
}

//===========================================================================
void HedgeSurface::addRegion(RevEngRegion* reg)
//===========================================================================
{
    regions_.push_back(reg);
    if (bbox_.dimension() == 0)
      bbox_ = reg->boundingBox();
    else
      bbox_.addUnionWith(reg->boundingBox());
}

//===========================================================================
void HedgeSurface::limitSurf(double diag)
//===========================================================================
{
  shared_ptr<ParamSurface> surf = surface();
  shared_ptr<ElementarySurface> elemsf =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf);
  if (elemsf.get() && !elemsf->isBounded())
    {
      if (diag < 0.0)
	diag = bbox_.low().dist(bbox_.high());
      double domain[4];
      regions_[0]->getDomain(domain);
      double midu = 0.5*(domain[0]+domain[1]);
      double midv = 0.5*(domain[2]+domain[3]);
      double fac = 1.0;
      shared_ptr<Plane> plane = dynamic_pointer_cast<Plane,ParamSurface>(surf);
      shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(surf);
      shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ParamSurface>(surf);
      if (plane.get())
	plane->setParameterBounds(midu-fac*diag, midv-fac*diag,
				  midu+fac*diag, midv+fac*diag);
      else if (cyl.get())
	cyl->setParamBoundsV(midv-fac*diag, midv+fac*diag);
      else if (cone.get())
	cone->setParamBoundsV(midv-fac*diag, midv+fac*diag);
    }
}

//===========================================================================
bool HedgeSurface::trimWithPoints(double aeps)
//===========================================================================
{
  // Extract data points
  int del = 5;
  int num_pts = numPoints();
  if (num_pts < 10)
    return false;
  vector<double> data(del*num_pts);
  int num_reg = numRegions();
  vector<double> extent(2*del);   // Limits for points in all coordinates
  for (int ki=0, kj=0; ki<num_reg; ++ki)
    {
      RevEngRegion *reg = getRegion(ki);
      int npts = reg->numPoints();
      for (int kr=0; kr<npts; ++kr, ++kj)
	{
	  RevEngPoint *pt = reg->getPoint(kr);
	  Vector2D uv = pt->getPar();
	  Vector3D xyz = pt->getPoint();
	  data[del*kj] = uv[0];
	  data[del*kj+1] = uv[1];
	  data[del*kj+2] = xyz[0];
	  data[del*kj+3] = xyz[1];
	  data[del*kj+4] = xyz[2];
	}
    }
  for (int kj=0; kj<del; ++kj)
    extent[2*kj] = extent[2*kj+1] = data[kj];
  for (int ki=1; ki<num_pts; ++ki)
    for (int kj=0; kj<del; ++kj)
      {
	extent[2*kj] = std::min(extent[2*kj],data[ki*del+kj]);
	extent[2*kj+1] = std::max(extent[2*kj+1],data[ki*del+kj]);
      }
  

  // Rough trimming curve
  double tol = 1.0e-5;
  int max_rec = 1;
  int nmb_div = std::min(num_pts/70, 8); //15;

  if (extent[1]-extent[0] < tol || extent[3]-extent[2] < tol)
    return false;
    // Compute trimming seqence
  bool only_outer = true;
  vector<vector<double> > seqs;
  //TrimUtils trimutil(points2, 1, domain);
  TrimUtils trimutil(&data[0], num_pts, del-2, &extent[0]);
  trimutil.computeTrimSeqs(max_rec, nmb_div, seqs, only_outer);
  
  // Compute trimming loop
  // First extract parts of the trimming sequences following iso trim curves
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);

  int nmb_loops = only_outer ? std::min(1,(int)seqs.size()) : (int)seqs.size();
  vector<vector<vector<double> > >  all_seqs(nmb_loops);
  if (nmb_loops == 0)
    return false;
  for (int kh=0; kh<nmb_loops; ++kh)
    all_seqs[kh].push_back(seqs[kh]);

  int nmb_match = 4;
  double eps = std::max(udel, vdel);
  vector<vector<shared_ptr<CurveOnSurface> > > loop(nmb_loops);
  for (int kh=0; kh<nmb_loops; ++kh)
    {
      for (int kj=0; kj<4; ++kj)
      	{
      	  int ix = (kj < 2) ? 0 : 1;
      	  for (int ki=0; ki<(int)all_seqs[kh].size();)
      	    {
      	      vector<vector<double> > split_seqs1 = 
      		TrimCrvUtils::extractConstParSeqs(all_seqs[kh][ki], ix, 
      						  extent[kj], nmb_match, 
      						  tol, eps);
      	      if (split_seqs1.size() > 1 || split_seqs1[0].size() == 4)
      		{
      		  all_seqs[kh].erase(all_seqs[kh].begin()+ki);
      		  all_seqs[kh].insert(all_seqs[kh].begin()+ki, 
      				      split_seqs1.begin(), split_seqs1.end());
      		  ki += (int)split_seqs1.size();
      		}
      	      else
      		++ki;
      	    }
      	}

      // Split the remaining sequences from the outer loop in kinks
      double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
      vector<vector<double> > split_seqs;
      for (size_t ki=0; ki<all_seqs[kh].size(); ++ki)
	{
	  vector<vector<double> > curr_seqs = 
	    TrimCrvUtils::splitCurvePointsInKinks(all_seqs[kh][ki], kink_tol);
	  split_seqs.insert(split_seqs.end(), curr_seqs.begin(), curr_seqs.end());
	}

      // Ensure a closed trimming loop
      TrimCrvUtils::makeConnectedLoop(split_seqs, tol);
  
      // Create trimming curves
      const int par_dim = 2;
      const int max_iter = 5;
      vector<shared_ptr<SplineCurve> > par_cvs;
      for (size_t ki = 0; ki < split_seqs.size(); ++ki)
	{
	  shared_ptr<SplineCurve> spline_cv_appr_2d
	    (TrimCrvUtils::approximateTrimPts(split_seqs[ki], par_dim, eps, 
					      max_iter));
	  par_cvs.push_back(spline_cv_appr_2d);
	}
  
      // The curve should be CCW.
      // Assume one outer and the rest inner (this is not necessarily true)
      const double int_tol = 1e-06;
      vector<shared_ptr<ParamCurve> > par_cvs2(par_cvs.begin(), par_cvs.end());
      bool loop_is_ccw = LoopUtils::loopIsCCW(par_cvs2, eps, int_tol);
      if ((kh==0 && !loop_is_ccw) || (kh>0 && loop_is_ccw))
	{
	  //MESSAGE("We should change direction of the loop cv!");
	  for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	    {
	      par_cvs[ki]->reverseParameterDirection();
	    }
	  reverse(par_cvs.begin(), par_cvs.end());
	}
      loop[kh].resize(par_cvs.size());
      for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	{
	  loop[kh][ki] =
	    shared_ptr<CurveOnSurface>(new CurveOnSurface(surface(), par_cvs[ki], true));
	}
    }
#ifdef DEBUG_TRIM
  std::ofstream fileout("trimcrvs.g2");
  for (size_t kh=0; kh<loop.size(); ++kh)
    {
      for (size_t kr=0; kr<loop[kh].size(); ++kr)
	{
	  loop[kh][kr]->parameterCurve()->writeStandardHeader(fileout);
	  loop[kh][kr]->parameterCurve()->write(fileout);
	}
    }
  if (del != 3)
    {
      fileout << "400 1 0 0" << std::endl;
      fileout << num_pts << std::endl;
      for (int ka=0; ka<num_pts; ++ka)
	{
	  fileout << data[ka*del] << " " << data[ka*del+1] << " 0.0" << std::endl;
	}
    }
#endif
  shared_ptr<BoundedSurface> bdsf(new BoundedSurface(surface(), loop, aeps));
  replaceSurf(bdsf);

  return true;
}

//===========================================================================
void HedgeSurface::store(std::ostream& os) const
//===========================================================================
{
  os << id_ << " " << surf_code_ << std::endl;
  surf_->writeStandardHeader(os);
  surf_->write(os);
  int profile = (profile_.get()) ? 1 : 0;
  os << profile << std::endl;
  if (profile_.get())
    {
      profile_->writeStandardHeader(os);
      profile_->write(os);
      os << sweep1_ << " " << sweep2_ << std::endl;
    }
}

//===========================================================================
void HedgeSurface::read(std::istream& is)
//===========================================================================
{
  GoTools::init();
  is >> id_ >> surf_code_;
  ObjectHeader header;
  header.read(is);
  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
  obj->read(is);
  surf_ =  dynamic_pointer_cast<ParamSurface,GeomObject>(obj);
  int profile;
  is >> profile;
  if (profile)
    {
      ObjectHeader header2;
      header2.read(is);
      shared_ptr<GeomObject> obj2(Factory::createObject(header2.classType()));
      obj2->read(is);
      profile_ =  dynamic_pointer_cast<SplineCurve,GeomObject>(obj);
      sweep1_ = Point(3);
      sweep1_.read(is);
      sweep2_ = Point(3);
      sweep2_.read(is);
    }
}
