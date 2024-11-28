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

#include "GoTools/compositemodel/RevEngEdge.h"
#include "GoTools/compositemodel/HedgeSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/intersections/Identity.h"
#include <fstream>

//#define DEBUG
//#define DEBUG_BLEND

using namespace Go;
using std::vector;
using std::istream;
using std::ostream;

//===========================================================================
RevEngEdge::RevEngEdge()
//===========================================================================
  : adjacent1_(0), adjacent2_(0), defined_blend_(0),
    blend_type_(BLEND_NOT_SET), outer1_(false), outer2_(false),
    distance_(0.0), radius_(0.0), surfchangecount_(0), extendcount_(0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(RevEngRegion* reg1, RevEngRegion* reg2)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0),
    blend_type_(BLEND_NOT_SET), outer1_(false), outer2_(false),
    distance_(0.0), radius_(0.0), surfchangecount_(0), extendcount_(0)
{
}

//===========================================================================
RevEngEdge::RevEngEdge(int type, RevEngRegion* reg1, 
		       vector<shared_ptr<CurveOnSurface> > cvs1,
		       bool out1, RevEngRegion* reg2,
		       vector<shared_ptr<CurveOnSurface> > cvs2,
		       bool out2, double distance, double radius)
//===========================================================================
  : adjacent1_(reg1), adjacent2_(reg2), defined_blend_(0), blend_type_(type),
    outer1_(out1), outer2_(out2), distance_(distance), radius_(radius),
    surfchangecount_(0), extendcount_(0)
{
  cvs1_.insert(cvs1_.end(), cvs1.begin(), cvs1.end());
  cvs2_.insert(cvs2_.end(), cvs2.begin(), cvs2.end());
}

//===========================================================================
RevEngEdge::~RevEngEdge()
//===========================================================================
{
  if (adjacent1_)
    adjacent1_->removeRevEngEdge(this);
  if (adjacent2_)
    adjacent2_->removeRevEngEdge(this);
  for (size_t ki=0; ki<blend_regs_.size(); ++ki)
    blend_regs_[ki]->removeAssociatedBlend();
  if (defined_blend_)
    defined_blend_->removeBlendEdge();
}


//===========================================================================
void RevEngEdge::setReg1(RevEngRegion *reg)
//===========================================================================
{
  adjacent1_ = reg;
  shared_ptr<ParamSurface> surf;
  if (reg->hasSurface())
    surf = reg->getSurface(0)->surface();
  else if (reg->hasBaseSf())
    surf = reg->getBase();

  if (surf.get())
    {
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	cvs1_[ki]->setUnderlyingSurface(surf);
    }
}

//===========================================================================
void RevEngEdge::setReg2(RevEngRegion *reg)
//===========================================================================
{
  adjacent2_ = reg;
  shared_ptr<ParamSurface> surf;
  if (reg->hasSurface())
    surf = reg->getSurface(0)->surface();
  else if (reg->hasBaseSf())
    surf = reg->getBase();

  if (surf.get())
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	cvs2_[ki]->setUnderlyingSurface(surf);
    }
}

//===========================================================================
void RevEngEdge::fixMismatchCurves(double tol)
//===========================================================================
{
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    cvs1_[ki]->fixMismatchCurves(tol);
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    cvs2_[ki]->fixMismatchCurves(tol);
}

//===========================================================================
void RevEngEdge::replaceSurf(RevEngRegion* reg,
			     shared_ptr<ParamSurface>& new_surf, double tol)
//===========================================================================
{
  if (reg == adjacent1_)
    {
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	{
	  cvs1_[ki]->setUnderlyingSurface(new_surf);
	  cvs1_[ki]->unsetParameterCurve();
	  cvs1_[ki]->ensureParCrvExistence(tol);
#ifdef DEBUG
	  if (!cvs1_[ki]->hasParameterCurve())
	    std::cout << "RevEngEdge::replaceSurf: No parameter curve" << std::endl;
#endif
	}
    }
  else if (reg == adjacent2_)
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	{
	  cvs2_[ki]->setUnderlyingSurface(new_surf);
	  cvs2_[ki]->unsetParameterCurve();
	  cvs2_[ki]->ensureParCrvExistence(tol);
#ifdef DEBUG
	  if (!cvs2_[ki]->hasParameterCurve())
	    std::cout << "RevEngEdge::replaceSurf: No parameter curve" << std::endl;
#endif
	}
     }
      
}

//===========================================================================
void RevEngEdge::store(ostream& os)
//===========================================================================
{
  os << Id_ << std::endl;
  int id1 = (adjacent1_) ? adjacent1_->getId() : -1;
  int id2 = (adjacent2_) ? adjacent2_->getId() : -1;
  int id3 = (defined_blend_ == 0) ? -1 : defined_blend_->getId();
  os << id1 << " " << id2 << " " << id3 << std::endl;
  os << cvs1_.size() << " " << cvs2_.size() << std::endl;
  
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      cvs1_[ki]->write(os);
      // bool space = cvs1_[ki]->hasSpaceCurve();
      // bool param = cvs1_[ki]->hasParameterCurve();
      // os << space << " " << param << std::endl;
      // if (space)
      // 	{
      // 	  shared_ptr<ParamCurve> spacecv = cvs1_[ki]->spaceCurve();
      // 	  shared_ptr<ParamCurve> parcv = cvs1_[ki]->parameterCurve();
      // 	  spacecv->writeStandardHeader(os);
      // 	  spacecv->write(os);
      // 	  parcv->writeStandardHeader(os);
      // 	  parcv->write(os);
      // 	}
    }
      
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    {
      cvs2_[ki]->write(os);
      // bool space = cvs2_[ki]->hasSpaceCurve();
      // bool param = cvs2_[ki]->hasParameterCurve();
      // os << space << " " << param << std::endl;
      // if (space)
      // 	{
      // 	  shared_ptr<ParamCurve> spacecv = cvs2_[ki]->spaceCurve();
      // 	  shared_ptr<ParamCurve> parcv = cvs2_[ki]->parameterCurve();
      // 	  spacecv->writeStandardHeader(os);
      // 	  spacecv->write(os);
      // 	  parcv->writeStandardHeader(os);
      // 	  parcv->write(os);
      // 	}
    }

  os << blend_regs_.size() << std::endl;
  for (size_t ki=0; ki < blend_regs_.size(); ++ki)
    os << blend_regs_[ki]->getId() << " ";
  os << std::endl;

  os << blend_type_ << " " << distance_ << " " << radius_ << " ";
  os << outer1_ << " " << outer2_ << " ";
  os << surfchangecount_ << " " << extendcount_ << std::endl;
}


//===========================================================================
void RevEngEdge::read(istream& is, int& reg_id1, int& reg_id2, int& reg_id3,
		      vector<int>& blend_id)
//===========================================================================
{
  is >> Id_;
  is >> reg_id1 >> reg_id2 >> reg_id3;
  int num_cv1, num_cv2;
  is >> num_cv1 >> num_cv2;
  if (num_cv1 > 0)
    cvs1_.resize(num_cv1);
  if (num_cv2 > 0)
    cvs2_.resize(num_cv2);
  for (int ka=0; ka<num_cv1; ++ka)
    {
      cvs1_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface());
      cvs1_[ka]->read(is);
      // int space, param;
      // is >> space >> param;
      // shared_ptr<ParamCurve> spacecv, parcv;
      // shared_ptr<ParamSurface> dummy;
      // if (space)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}

      // if (param)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}
      
      // cvs1_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
      // 								spacecv, false,
      // 								-1, -1, 0.0,
      // 								-1, true));
    }

  for (int ka=0; ka<num_cv2; ++ka)
    {
      cvs2_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface());
      cvs2_[ka]->read(is);
      // int space, param;
      // is >> space >> param;
      // shared_ptr<ParamCurve> spacecv, parcv;
      // shared_ptr<ParamSurface> dummy;
      // if (space)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  spacecv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}

      // if (param)
      // 	{
      // 	  ObjectHeader header;
      // 	  header.read(is);
      // 	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      // 	  obj->read(is);
      // 	  parcv = dynamic_pointer_cast<ParamCurve,GeomObject>(obj);
      // 	}
      
      // cvs2_[ka] = shared_ptr<CurveOnSurface>(new CurveOnSurface(dummy, parcv,
      // 								spacecv, false,
      // 								-1, -1, 0.0,
      // 								-1, true));
    }

  int num_blend_reg;
  is >> num_blend_reg;
  if (num_blend_reg > 0)
    blend_id.resize(num_blend_reg);
  for (int ka=0; ka<num_blend_reg; ++ka)
    is >> blend_id[ka];
  
  is >> blend_type_ >> distance_ >> radius_ >> outer1_ >> outer2_;
  is >> surfchangecount_ >> extendcount_;
}


//===========================================================================
vector<shared_ptr<ParamCurve> > RevEngEdge::getSpaceCurves()
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > curves;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    curves.push_back(cvs1_[ki]->spaceCurve());

  return curves;
}

//===========================================================================
void RevEngEdge::getCrvEndPoints(Point& pos1, Point& pos2)
//===========================================================================
{
  if (cvs1_.size() > 0)
    {
      cvs1_[0]->point(pos1, cvs1_[0]->startparam());
      cvs1_[cvs1_.size()-1]->point(pos2, cvs1_[cvs1_.size()-1]->endparam());
    }
}

//===========================================================================
void RevEngEdge::closestPoint(const Point& pos, double& par, Point& close,
			      double& dist)
//===========================================================================
{
  int ix;
  closestPoint(pos, par, close, dist, ix);
}

//===========================================================================
void RevEngEdge::closestPoint(const Point& pos, double& par, Point& close,
			      double& dist, int& ix)
//===========================================================================
{
  dist = std::numeric_limits<double>::max();
  ix = -1;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      double par1, dist1;
      Point close1;
      cvs1_[ki]->closestPoint(pos, cvs1_[ki]->startparam(),
			      cvs1_[ki]->endparam(), par1, close1, dist1);
      if (dist1 < dist)
	{
	  par = par1;
	  close = close1;
	  dist = dist1;
	  ix = (int)ki;
	}
    }
}

//===========================================================================
int RevEngEdge::closedSfAtEnd(double tol, double& par, Point& pos, bool at_start)
//===========================================================================
{
  if (cvs1_.size() == 0)
    return 0;
  
  shared_ptr<ParamSurface> surf1 = cvs1_[0]->underlyingSurface();
  shared_ptr<ParamSurface> surf2 = cvs2_[0]->underlyingSurface();

  // The parameter interval and the space curve of the two curve sequences
  // correspond
  int ix = (at_start) ? 0 : (int)cvs1_.size()-1;
  par = (at_start) ? cvs1_[ix]->startparam() : cvs1_[ix]->endparam();
  pos = cvs1_[ix]->ParamCurve::point(par);

  double eps = 1.0e-9;
  double u1, u2, v1, v2, d1, d2;
  Point close1, close2;
  surf1->closestBoundaryPoint(pos, u1, v1, close1, d1, eps);
  surf2->closestBoundaryPoint(pos, u2, v2, close2, d2, eps);
  if (d1 > tol && d2 > tol)
    return 0;

  if (d1 <= tol)
    {
      RectDomain dom = surf1->containingDomain();
      if (fabs(u1-dom.umin()) <= eps)
	{
	  Point pos2 = surf1->point(dom.umax(), v1);
	  if (pos.dist(pos2) <= tol)
	    return 1;
	}
      else if (fabs(dom.umax()-u1) <= eps)
	{
	  Point pos2 = surf1->point(dom.umin(), v1);
	  if (pos.dist(pos2) <= tol)
	    return 1;
	}
      if (fabs(v1-dom.vmin()) <= eps)
	{
	  Point pos2 = surf1->point(u1, dom.vmax());
	  if (pos.dist(pos2) <= tol)
	    return 1;
	}
      else if (fabs(dom.vmax()-v1) <= eps)
	{
	  Point pos2 = surf1->point(u1, dom.vmin());
	  if (pos.dist(pos2) <= tol)
	    return 1;
	}
    }

  if (d2 <= tol)
    {
      RectDomain dom = surf2->containingDomain();
      if (fabs(u2-dom.umin()) <= eps)
	{
	  Point pos2 = surf2->point(dom.umax(), v2);
	  if (pos.dist(pos2) <= tol)
	    return 2;
	}
      else if (fabs(dom.umax()-u2) <= eps)
	{
	  Point pos2 = surf2->point(dom.umin(), v2);
	  if (pos.dist(pos2))
	    return 2;
	}
      if (fabs(v2-dom.vmin()) <= eps)
	{
	  Point pos2 = surf2->point(u2, dom.vmax());
	  if (pos.dist(pos2) <= tol)
	    return 2;
	}
      else if (fabs(dom.vmax()-v2) <= eps)
	{
	  Point pos2 = surf2->point(u2, dom.vmin());
	  if (pos.dist(pos2) <= tol)
	    return 2;
	}
    }

  return 0;
}

//===========================================================================
bool RevEngEdge::isClosed(double tol)
//===========================================================================
{
  if (cvs1_.size() == 0)
    return false;
  double par1 = cvs1_[0]->startparam();
  double par2 = cvs1_[cvs1_.size()-1]->endparam();
  Point pos1 = cvs1_[0]->ParamCurve::point(par1);
  Point pos2 = cvs1_[cvs1_.size()-1]->ParamCurve::point(par2);
  double dist = pos1.dist(pos2);
  return (dist <= tol);
}

//===========================================================================
bool RevEngEdge::isAdjacent(RevEngEdge* other, double tol, double& par1, double& par2)
//===========================================================================
{
  if (cvs1_.size() == 0)
    return false;
  double tpar1 = cvs1_[0]->startparam();
  double tpar2 = cvs1_[cvs1_.size()-1]->endparam();
  Point pos1 = cvs1_[0]->ParamCurve::point(tpar1);
  Point pos2 = cvs1_[cvs1_.size()-1]->ParamCurve::point(tpar2);
  double tpar3 = other->startparam();
  double tpar4 = other->endparam();
  Point pos3 = other->point(tpar3);
  Point pos4 = other->point(tpar4);

  double dd1 = pos1.dist(pos3);
  double dd2 = pos1.dist(pos4);
  double dd3 = pos2.dist(pos3);
  double dd4 = pos2.dist(pos4);
  bool adjacent = false;
  if (dd1 <= tol && dd1 <= std::min(dd2, std::min(dd3,dd4)))
    {
      par1 = tpar1;
      par2 = tpar3;
      adjacent = true;
    }
  else if (dd2 <= tol && dd2 <= std::min(dd3, dd4))
    {
      par1 = tpar1;
      par2 = tpar4;
      adjacent = true;
    }
    else if (dd3 <= tol && dd3 <= dd4)
    {
      par1 = tpar2;
      par2 = tpar3;
      adjacent = true;
    }
    else if (dd4 <= tol)
     {
      par1 = tpar2;
      par2 = tpar4;
      adjacent = true;
    }
    return adjacent;
}

//===========================================================================
Point RevEngEdge::point(double par)
//===========================================================================
{
  Point dummy;
  if (cvs1_.size() == 0)
    return dummy;
  
  size_t ix = 0;
  for (ix; ix<cvs1_.size() && cvs1_[ix]->endparam() <= par; ++ix);
  if (ix == cvs1_.size())
    --ix;

  return cvs1_[ix]->ParamCurve::point(par);
}

//===========================================================================
bool RevEngEdge::append(RevEngEdge* other, double tol)
//===========================================================================
{
  double tpar1, tpar2;
  bool adjacent = isAdjacent(other, tol, tpar1, tpar2);
  if (!adjacent)
    return false;   // Cannot append
  
  RevEngRegion *adj3, *adj4;
  other->getAdjacent(adj3, adj4);
  bool first;
  if (adjacent1_ == adj3 && adjacent2_ == adj4)
    first = true;
  else if (adjacent1_ == adj4 && adjacent2_ == adj3)
    first = false;
  else
    return false;  // Not an appendable configuration

  size_t ncv1 = cvs1_.size();
  if (ncv1 == 0)
    return false;
  
  size_t ix1 = 0;
  for (; ix1<ncv1 && cvs1_[ix1]->endparam() <= tpar1; ++ix1);
  if (ix1 == ncv1)
    --ix1;
  Point ppos1 = cvs1_[ix1]->faceParameter(tpar1);
  Point ppos2 = cvs2_[ix1]->faceParameter(tpar1);
  
  vector<shared_ptr<CurveOnSurface> > cvs3, cvs4;
  other->getCurve(cvs3, true);
  other->getCurve(cvs4, false);
  size_t ncv3 = cvs3.size();
  size_t ix2 = 0;
  for (; ix2<ncv3 && cvs3[ix2]->endparam() <= tpar2; ++ix2);
  if (ix2 == ncv3)
    --ix2;
  Point ppos3 = cvs3[ix2]->faceParameter(tpar2);
  Point ppos4 = cvs4[ix2]->faceParameter(tpar2);

  if ((!adjacent1_->hasSurface()) || (!adjacent2_->hasSurface()))
    return false;  // Unexpected
  
  shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  Point par_eps1 = SurfaceTools::getParEpsilon(*surf1, tol);
  Point par_eps2 = SurfaceTools::getParEpsilon(*surf2, tol);
  double epspar1 = 0.5*(par_eps1[0], par_eps1[1]);
  double epspar2 = 0.5*(par_eps2[0], par_eps2[1]);
  double dd1 = ppos1.dist(first ? ppos3 : ppos4);
  double dd2 = ppos2.dist(first ? ppos4 : ppos3);

  if (dd1 > epspar1 || dd2 > epspar2)
    return false;

  // Try to append
  // Copy curves to keep originals in case of failure
  vector<shared_ptr<CurveOnSurface> > cvs1_2(ncv1);
  for (size_t ki=0; ki<ncv1; ++ki)
    cvs1_2[ki] = shared_ptr<CurveOnSurface>(cvs1_[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs2_2(ncv1);
  for (size_t ki=0; ki<ncv1; ++ki)
    cvs2_2[ki] = shared_ptr<CurveOnSurface>(cvs2_[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs3_2(ncv3);
  for (size_t ki=0; ki<ncv3; ++ki)
    cvs3_2[ki] = shared_ptr<CurveOnSurface>(cvs3[ki]->clone());
  vector<shared_ptr<CurveOnSurface> > cvs4_2(ncv3);
  for (size_t ki=0; ki<ncv3; ++ki)
    cvs4_2[ki] = shared_ptr<CurveOnSurface>(cvs4[ki]->clone());

  if (fabs(tpar1-cvs1_2[0]->startparam()) < fabs(cvs1_2[ncv1-1]->endparam()-tpar1))
    {
      for (size_t ki=0; ki<ncv1; ++ki)
	{
	  cvs1_2[ki]->reverseParameterDirection();
	  cvs2_2[ki]->reverseParameterDirection();
	}
      for (size_t ki=0; ki<ncv1/2; ++ki)
	{
	  std::swap(cvs1_2[ki], cvs1_2[ncv1-1-ki]);
	  std::swap(cvs2_2[ki], cvs2_2[ncv1-1-ki]);
	}
    }

  if (fabs(tpar2-cvs3_2[0]->startparam()) > fabs(cvs3_2[ncv3-1]->endparam()-tpar2))
    {
      for (size_t ki=0; ki<ncv3; ++ki)
	{
	  cvs3_2[ki]->reverseParameterDirection();
	  cvs4_2[ki]->reverseParameterDirection();
	}
      for (size_t ki=0; ki<ncv3/2; ++ki)
	{
	  std::swap(cvs3_2[ki], cvs3_2[ncv3-1-ki]);
	  std::swap(cvs4_2[ki], cvs4_2[ncv3-1-ki]);
	}
    }

  // Do append
  double dist1, dist2;
  if (first)
    {
      cvs1_2[ncv1-1]->appendCurve(cvs3_2[0].get(), 1, dist1, false, epspar1);
      cvs2_2[ncv1-1]->appendCurve(cvs4_2[0].get(), 1, dist2, false, epspar2);
    }
  else
    {
      cvs1_2[ncv1-1]->appendCurve(cvs4_2[0].get(), 1, dist1, false, epspar1);
      cvs2_2[ncv1-1]->appendCurve(cvs2_2[0].get(), 1, dist2, false, epspar2);
    }

  if (dist1 > tol || dist2 > tol)
    return false;

  // Check that the parameter curves of the joined curves exists
  if ((!cvs1_2[ncv1-1]->hasParameterCurve()) ||
      (!cvs2_2[ncv1-1]->hasParameterCurve()))
    return false;

  // Replace curves
  cvs1_.clear();
  cvs1_.insert(cvs1_.end(), cvs1_2.begin(), cvs1_2.end());
  if (first && cvs3_2.size() > 1)
    cvs1_.insert(cvs1_.end(), cvs3_2.begin()+1, cvs3_2.end());
  else if ((!first) && cvs4_2.size() > 1)
    cvs1_.insert(cvs1_.end(), cvs4_2.begin()+1, cvs4_2.end());

  cvs2_.clear();
  cvs2_.insert(cvs2_.end(), cvs2_2.begin(), cvs2_2.end());
  if (first && cvs4_2.size() > 1)
    cvs2_.insert(cvs2_.end(), cvs4_2.begin()+1, cvs4_2.end());
  else if ((!first) && cvs4_2.size() > 1)
    cvs2_.insert(cvs2_.end(), cvs3_2.begin()+1, cvs3_2.end());

  for (size_t kj=0; kj<other->blend_regs_.size(); ++kj)
    other->blend_regs_[kj]->setAssociatedBlend(this);
  blend_regs_.insert(blend_regs_.end(), other->blend_regs_.begin(),
		     other->blend_regs_.end());
  other->clearBlendRegions();

  distance_ = 0.5*(distance_ + other->distance_);
  radius_ = 0.5*(radius_ + other->radius_);

  return true;
}

//===========================================================================
int RevEngEdge::missingParCrv()
//===========================================================================
{
  int missing = 0;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    if (!cvs1_[ki]->hasParameterCurve())
      {
	missing += 1;
	break;
      }
  
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    if (!cvs2_[ki]->hasParameterCurve())
      {
	missing += 2;
	break;
      }
  return missing;
}

//===========================================================================
void RevEngEdge::splitAtSeam(double tol,
			     vector<shared_ptr<RevEngEdge> >& added_edgs,
			     vector<shared_ptr<RevEngRegion> >& added_regs,
			     vector<shared_ptr<HedgeSurface> >& added_sfs)
//===========================================================================
{
  if ((!adjacent1_->hasSurface()) || (!adjacent2_->hasSurface()))
    return;  // Something is wrong
#ifdef DEBUG
  std::ofstream of("edge_split.g2");
  adjacent1_->writeRegionPoints(of);
  adjacent2_->writeRegionPoints(of);
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp = cvs1_[ki]->spaceCurve();
      tmp->writeStandardHeader(of);
      tmp->write(of);
    }
#endif
  shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  shared_ptr<ElementarySurface> elem1 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf1);
  shared_ptr<ElementarySurface> elem2 =
    dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf2);

  bool close_u1=false, close_u2=false, close_v1=false, close_v2=false;
  if (elem1.get())
    elem1->isClosed(close_u1, close_v1);
  if (elem2.get())
    elem2->isClosed(close_u2, close_v2);

  RectDomain dom1 = surf1->containingDomain();
  RectDomain dom2 = surf2->containingDomain();
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      if (!cvs1_[ki]->hasParameterCurve())
	{
	  // Try to split curve
	  if (close_u1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.umin(), false);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 1;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	  
	  if (close_v1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.vmin(), true);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 1;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	}
    }

  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    {
      if (!cvs2_[ki]->hasParameterCurve())
	{
	  // Try to split curve
	  if (close_u2)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf2->constParamCurves(dom2.umin(), false);

	      shared_ptr<ParamCurve> space = cvs2_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(), 
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 2;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	  
	  if (close_v1)
	    {
	      vector<shared_ptr<ParamCurve> > seam =
		surf1->constParamCurves(dom1.vmin(), true);

	      shared_ptr<ParamCurve> space = cvs1_[ki]->spaceCurve();
	      double par1, par2, dist;
	      Point ptc1, ptc2;
	      ClosestPoint::closestPtCurves(space.get(), seam[0].get(),
					    par1, par2, dist, ptc1, ptc2);
	      if (dist < tol)
		{
#ifdef DEBUG
		  std::ofstream of2("int_with_seam.g2");
		  space->writeStandardHeader(of2);
		  space->write(of2);
		  seam[0]->writeStandardHeader(of2);
		  seam[0]->write(of2);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc1 << std::endl;
		  of2 << "400 1 0 4 0 255 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << ptc2 << std::endl;
#endif
		  int nr = 2;
		  shared_ptr<RevEngEdge> new_edge =
		    doSplit(ki, nr, par1, tol, added_regs, added_sfs);
		  if (new_edge.get())
		    {
		      added_edgs.push_back(new_edge);
		      new_edge->splitAtSeam(tol, added_edgs, added_regs, added_sfs);
		    }
		}
	    }
	}
    }

  
}


//===========================================================================
shared_ptr<RevEngEdge>
RevEngEdge::doSplit(size_t ix, int side, double par, double tol,
		    vector<shared_ptr<RevEngRegion> >& added_regs,
		    vector<shared_ptr<HedgeSurface> >& added_sfs)
//===========================================================================
{
  double eps = 1.0e-9;
  shared_ptr<RevEngEdge> new_edg;
  if (ix >= cvs1_.size() || cvs1_.size() == 0)
    return new_edg;

  if (par <= cvs1_[ix]->startparam()+eps || par >= cvs1_[ix]->endparam()-eps)
    return new_edg;

  // Split curves
  vector<shared_ptr<ParamCurve> > sub1 = cvs1_[ix]->split(par);
  vector<shared_ptr<ParamCurve> > sub2 = cvs2_[ix]->split(par);
  if (sub1.size() != 2 || sub2.size() != 2)
    return new_edg;
  
  // Split blend regions
  vector<RevEngRegion*> move_reg;
  for (size_t ki=0; ki<blend_regs_.size(); )
    {
      vector<RevEngPoint*> points = blend_regs_[ki]->getPoints();
      vector<RevEngPoint*> keep, move;
      
      for (size_t kj=0; kj<points.size(); ++kj)
	{
	  Vector3D xyz = points[kj]->getPoint();
	  Point pnt(xyz[0], xyz[1], xyz[2]);
	  double par2, dist;
	  Point close;
	  int ix2;
	  closestPoint(pnt, par2, close, dist, ix2);
	  if (ix2 < ix || par2 <= par)
	    keep.push_back(points[kj]);
	  else
	    move.push_back(points[kj]);
	}

      if (move.size() == 0)
	++ki; // Do nothing
      else if (keep.size() == 0)
	{
	  // Move regions to new edge
	  move_reg.push_back(blend_regs_[ki]);
	  blend_regs_.erase(blend_regs_.begin()+ki);
	}
      else
	{
	  // Split region
	  blend_regs_[ki]->removePoints(move);  // This is not the most effective
	  // method, but the simplest to implement
	  blend_regs_[ki]->updateInfo();

	  shared_ptr<RevEngRegion> new_reg(new RevEngRegion(blend_regs_[ki]->getClassificationType(),
							    blend_regs_[ki]->getEdgeClassificationType(),
							    move));
	  added_regs.push_back(new_reg);
	  move_reg.push_back(new_reg.get());
	  ++ki;
	}
    }

  // Distribute curves
  vector<shared_ptr<CurveOnSurface> > cvs1_2, cvs2_2;
  cvs1_2.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub1[1]));
  cvs1_2[0]->ensureParCrvExistence(tol);
  for (size_t ki=ix+1; ki<cvs1_.size(); ++ki)
    cvs1_2.push_back(cvs1_[ki]);
  cvs2_2.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub2[1]));
  cvs2_2[0]->ensureParCrvExistence(tol);
  for (size_t ki=ix+1; ki<cvs2_.size(); ++ki)
    cvs2_2.push_back(cvs2_[ki]);

  cvs1_[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub1[0]);
  cvs1_[ix]->ensureParCrvExistence(tol);
  cvs2_[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub2[0]);
  cvs2_[ix]->ensureParCrvExistence(tol);
  if (ix < cvs1_.size()-1)
    {
      cvs1_.erase(cvs1_.begin()+ix+1, cvs1_.end());
      cvs2_.erase(cvs1_.begin()+ix+1, cvs1_.end());
    }

  if (defined_blend_)
    {
      // Split associated surface
      int stop_blend = 1;
    }

  new_edg = shared_ptr<RevEngEdge>(new RevEngEdge(blend_type_, adjacent1_,
						  cvs1_2, outer1_, adjacent2_,
						  cvs2_2, outer2_, radius_,
						  distance_));
  adjacent1_->addRevEdge(new_edg.get());
  adjacent2_->addRevEdge(new_edg.get());
  if (move_reg.size() > 0)
    {
      for (size_t kr=0; kr<move_reg.size(); ++kr)
	move_reg[kr]->setAssociatedBlend(new_edg.get());
      new_edg->addBlendRegions(move_reg);
    }
  
  return new_edg;
}

//===========================================================================
bool
RevEngEdge::updateCurve(double int_tol, double tol, double len)
//===========================================================================
{
  double eps = 1.0e-6;
  if (surfchangecount_ == 0)
    return false;

  if (cvs1_.size() == 0)
    return false;
  
  shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  if (!surf1->isBounded())
    adjacent1_->getSurface(0)->limitSurf(len);
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  if (!surf2->isBounded())
    adjacent2_->getSurface(0)->limitSurf(len);
  shared_ptr<BoundedSurface> bd1, bd2;
  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
  BoundedUtils::getSurfaceIntersections(surf1, surf2, int_tol,
					int_cvs1, bd1,
					int_cvs2, bd2);
#ifdef DEBUG_BLEND
  std::ofstream of_int("intcurves_edge.g2");
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = cvs1_[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
  for (size_t ki=0; ki<int_cvs1.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = int_cvs1[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
#endif
  //  #if 0
  if (isClosed(tol))
    replaceCurves(int_cvs1, int_cvs2);  // The seam might have moved
      
  else 
    {
      if (cvs1_[0]->isClosed())
	{
	  cvs1_[0] = int_cvs1[0];
	  cvs2_[0] = int_cvs2[0];
	}
      else
	{
	  size_t st = cvs1_.size() - 1;
	  Point pos1, pos2, pos3, pos4;
	  double tdel1 = cvs1_[0]->endparam() - cvs1_[0]->startparam();
	  double tdel2 = cvs1_[st]->endparam() - cvs1_[st]->startparam();
	  cvs1_[0]->point(pos1, cvs1_[0]->startparam());
	  cvs1_[0]->point(pos3, cvs1_[0]->startparam() + 0.1*tdel1);
	  cvs1_[st]->point(pos2, cvs1_[st]->endparam());
	  cvs1_[st]->point(pos4, cvs1_[st]->endparam() - 0.1*tdel2);
	  double tp1, tp2;
	  double td1 = std::numeric_limits<double>::max();
	  double td2 = std::numeric_limits<double>::max();
	  double td3 = std::numeric_limits<double>::max();
	  double td4 = std::numeric_limits<double>::max();
	  int ix1 = -1, ix2 = -1;
	  for (size_t kr=0; kr<int_cvs1.size(); ++kr)
	    {
	      double tp1_2, tp2_2, tp1_3, tp2_3, td1_2, td2_2, td1_3, td2_3;
	      Point cl1_2, cl2_2;
	      double tmin = int_cvs1[kr]->startparam();
	      double tmax = int_cvs1[kr]->endparam();
	      int_cvs1[kr]->closestPoint(pos1, tmin, tmax, tp1_2, cl1_2, td1_2);
	      int_cvs1[kr]->closestPoint(pos2, tmin, tmax, tp2_2, cl2_2, td2_2);
	      int_cvs1[kr]->closestPoint(pos3, tmin, tmax, tp1_3, cl1_2, td1_3);
	      int_cvs1[kr]->closestPoint(pos4, tmin, tmax, tp2_3, cl2_2, td2_3);
	      if (td1_2 < td1-eps || (td1_2 <= td1+eps && td1_3 < td3))
		{
		  ix1 = (int)kr;
		  tp1 = tp1_2;
		  td1 = td1_2;
		  td3 = td1_3;
		}
	      if (td2_2 < td2-eps || (td2_2 <= td2+eps && td2_3 < td4))
		{
		  ix2 = (int)kr;
		  tp2 = tp2_2;
		  td2 = td2_2;
		  td4 = td2_3;
		}
	    }
	  
	  if (ix2 - ix1 + 1 == (int)cvs1_.size() && tp1 < tp2)
	    {
	      for (size_t kr=0; kr<int_cvs1.size(); )
		{
		  double tp3 = std::max(tp1, int_cvs1[kr]->startparam());
		  double tp4 = std::min(tp2, int_cvs1[kr]->endparam());
		  if (ix1 <= ix2 && ((int)kr < ix1 || (int)kr > ix2))
		    {
		      int_cvs1.erase(int_cvs1.begin()+kr);
		      int_cvs2.erase(int_cvs2.begin()+kr);
		      if ((int)kr < ix1)
			ix1--;
		      if ((int)kr < ix2)
			ix2--;
		    }
		  else if (tp4 > tp3 && tp3 < int_cvs1[kr]->endparam() &&
			   tp4 > int_cvs1[kr]->startparam() && 
			   (tp3 > int_cvs1[kr]->startparam() || tp4 < int_cvs1[kr]->endparam()))
		    {
		      shared_ptr<CurveOnSurface> sub1(int_cvs1[kr]->subCurve(tp3, tp4));
		      cvs1_[kr] = sub1;
		      shared_ptr<CurveOnSurface> sub2(int_cvs2[kr]->subCurve(tp3, tp4));
		      cvs2_[kr] = sub2;
#ifdef DEBUG_BLEND
		      shared_ptr<ParamCurve> tmp_cv = sub1->spaceCurve();
		      tmp_cv->writeStandardHeader(of_int);
		      tmp_cv->write(of_int);
#endif
		      ++kr;
		    }
		  else
		    {
		      cvs1_[kr] = int_cvs1[kr];
		      cvs2_[kr] = int_cvs2[kr];
		      ++kr;
		    }
		}
	    }
	}
    }

  resetSurfChangeCount();
  return true;

}

//===========================================================================
void
RevEngEdge::updateParCurve(RevEngRegion *adj, double int_tol)
//===========================================================================
{
  if (adj == adjacent1_)
    {
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	cvs1_[ki]->ensureParCrvExistence(int_tol);
    }
  else if (adj == adjacent2_)
    {
      for (size_t ki=0; ki<cvs2_.size(); ++ki)
	cvs2_[ki]->ensureParCrvExistence(int_tol);
    }
}

//===========================================================================
bool
RevEngEdge::extendCurve(double int_tol, double tol, double anglim, 
			double len, double lenlim, double blendlim,
			vector<shared_ptr<RevEngRegion> >& added_regions,
			vector<vector<RevEngPoint*> >& extract_groups,
			vector<HedgeSurface*>& out_sfs)
//===========================================================================
{
  double eps = 1.0e-6;
  if (extendcount_ == 0)
    return false;

  if (cvs1_.size() == 0)
    return false;
  
  if (isClosed(tol))
    {
      extendcount_ = 0;
      return false;  // Nothing to do here
    }

    shared_ptr<ParamSurface> surf1 = adjacent1_->getSurface(0)->surface();
  if (!surf1->isBounded())
    adjacent1_->getSurface(0)->limitSurf(len);
  shared_ptr<ParamSurface> surf2 = adjacent2_->getSurface(0)->surface();
  if (!surf2->isBounded())
    adjacent2_->getSurface(0)->limitSurf(len);
  shared_ptr<BoundedSurface> bd1, bd2;
  vector<shared_ptr<CurveOnSurface> > int_cvs1, int_cvs2;
  BoundedUtils::getSurfaceIntersections(surf1, surf2, int_tol,
					int_cvs1, bd1,
					int_cvs2, bd2);
#ifdef DEBUG_BLEND
  std::ofstream of_int("intcurves_edge.g2");
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = cvs1_[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
  for (size_t ki=0; ki<int_cvs1.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = int_cvs1[ki]->spaceCurve();
      tmp_cv->writeStandardHeader(of_int);
      tmp_cv->write(of_int);
    }
#endif

  if (int_cvs1.size() == 0)
    return false;

  if (int_cvs1.size() > cvs1_.size())
    {
      // Remove possible extra curves
      vector<Point> mid(cvs1_.size());
      for (size_t ki=0; ki<cvs1_.size(); ++ki)
	cvs1_[ki]->point(mid[ki], 0.5*(cvs1_[ki]->startparam()+cvs1_[ki]->endparam()));
      size_t kr, kh;
      for (kr=0; kr<int_cvs1.size(); )
	{
	  for (kh=0; kh<mid.size(); ++kh)
	    {
	      double tpar, dist;
	      Point close;
	      int_cvs1[kr]->closestPoint(mid[kh], int_cvs1[kr]->startparam(),
					 int_cvs1[kr]->endparam(), tpar, close, dist);
	      if (dist <= tol)
		break;
	    }
	  if (kh == mid.size())
	    {
	      int_cvs1.erase(int_cvs1.begin()+kr);
	      int_cvs2.erase(int_cvs2.begin()+kr);
	    }
	  else
	    ++kr;
	}
    }

  // Limit intersection curves to relevant intervals
  vector<pair<double,double> > t1_t2, t3_t4;
  bool OK1 =
    adjacent1_->getCurveRestriction(int_cvs1, tol, anglim, t1_t2);
  bool OK2 =
    adjacent2_->getCurveRestriction(int_cvs2, tol, anglim, t3_t4);

  for (int ka=(int)int_cvs1.size()-1; ka>=0; --ka)
    {
      double t1 = std::max(t1_t2[ka].first, t3_t4[ka].first);
      double t2 = std::min(t1_t2[ka].second, t3_t4[ka].second);
      if (t2 > t1 && (t1 > int_cvs1[ka]->startparam() ||
		      t2 < int_cvs1[ka]->endparam()))
	{
	  double pmin = std::max(t1, int_cvs1[ka]->startparam());
	  double pmax = std::min(t2, int_cvs1[ka]->endparam());
	  shared_ptr<CurveOnSurface> sub1(int_cvs1[ka]->subCurve(pmin,pmax));
	  int_cvs1[ka] = sub1;
	  shared_ptr<CurveOnSurface> sub2(int_cvs2[ka]->subCurve(pmin,pmax));
	  int_cvs2[ka] = sub2;
	}

      if (t2 <= t1 || int_cvs1[ka]->estimatedCurveLength() < lenlim)
	{
	  int_cvs1.erase(int_cvs1.begin()+ka);
	  int_cvs2.erase(int_cvs2.begin()+ka);
	}
    }
  if (int_cvs1.size() == 0)
    return false;
  
  vector<RevEngRegion*> common_reg2 =
    adjacent1_->commonAdjacent(adjacent2_);
  for (size_t kj=0; kj<common_reg2.size(); )
    {
      // Check if the common region is registered already or is unfeasable for a blend
      size_t kr=0;
      for (kr=0; kr<blend_regs_.size(); ++kr)
	if (common_reg2[kj] == blend_regs_[kr])
	  break;
      
      if (kr < blend_regs_.size() || common_reg2[kj]->hasAssociatedBlend() ||
	  common_reg2[kj]->hasBlendEdge())
  	common_reg2.erase(common_reg2.begin()+kj);
      else
  	++kj;
    }

  // Check extension
  vector<RevEngPoint*> bd_pts1 = adjacent1_->extractBdPoints(); 
  vector<RevEngPoint*> bd_pts2 = adjacent2_->extractBdPoints();
  if (bd_pts1.size() == 0 || bd_pts2.size() == 0)
    return false;

  int num_in_lim1=0, num_in_lim2=0;
  vector<pair<double, double> > t5_t6, t7_t8;
  vector<double> wwd1, wwd2;
  adjacent1_->estimateBlendDimensions(int_cvs1, bd_pts1, tol, blendlim,
				      t5_t6, wwd1, num_in_lim1);

  adjacent2_->estimateBlendDimensions(int_cvs2, bd_pts2, tol, blendlim,
				      t7_t8, wwd2, num_in_lim2);
  if (num_in_lim1 == 0 || num_in_lim2 == 0)
    return false;
  if (t5_t6.size() == 0 || t7_t8.size() == 0)
    return false;

  Point midp = point(0.5*(startparam() + endparam()));

  int ix1 = -1, ix2 = -1;
  double td1 = std::numeric_limits<double>::max();
  double td2 = std::numeric_limits<double>::max();
  for (size_t kj=0; kj<t5_t6.size(); ++kj)
    {
      double tmid = 0.5*(t5_t6[kj].first + t5_t6[kj].second);
      size_t kr;
      for (kr=0; kr<int_cvs1.size(); ++kr)
	if (int_cvs1[kr]->startparam()-eps <= tmid && int_cvs1[kr]->endparam()+eps >= tmid)
	  break;
      if (kr == int_cvs1.size())
	return false;  // Should not happen

      double tpar, tdist;
      Point close;
      int_cvs1[kr]->closestPoint(midp, t5_t6[kj].first, t5_t6[kj].second,
				 tpar, close, tdist);
      if (tdist < td1)
	{
	  ix1 = (int)kj;
	  td1 = tdist;
	}
    }
  
  for (size_t kj=0; kj<t7_t8.size(); ++kj)
    {
      double tmid = 0.5*(t7_t8[kj].first + t7_t8[kj].second);
      size_t kr;
      for (kr=0; kr<int_cvs2.size(); ++kr)
	if (int_cvs2[kr]->startparam()-eps <= tmid && int_cvs2[kr]->endparam()+eps >= tmid)
	  break;
      if (kr == int_cvs1.size())
	return false;  // Should not happen

      double tpar, tdist;
      Point close;
      int_cvs2[kr]->closestPoint(midp, t7_t8[kj].first, t7_t8[kj].second,
				 tpar, close, tdist);
      if (tdist < td2)
	{
	  ix2 = (int)kj;
	  td2 = tdist;
	}
    }

  double t1 = std::max(t5_t6[ix1].first, t7_t8[ix2].first);
  double t2 = std::min(t5_t6[ix1].second, t7_t8[ix2].second);
  if (t2 <= t1+eps)
    return false;  // Nothing left
  
  for (int ka=(int)int_cvs1.size()-1; ka>=0; --ka)
    {
      if (t2 > t1 && (t1 > int_cvs1[ka]->startparam() ||
		      t2 < int_cvs1[ka]->endparam()))
	{
	  double pmin = std::max(t1, int_cvs1[ka]->startparam());
	  double pmax = std::min(t2, int_cvs1[ka]->endparam());
	  shared_ptr<CurveOnSurface> sub1(int_cvs1[ka]->subCurve(pmin,pmax));
	  int_cvs1[ka] = sub1;
	  shared_ptr<CurveOnSurface> sub2(int_cvs2[ka]->subCurve(pmin,pmax));
	  int_cvs2[ka] = sub2;
	}

      if (t2 <= t1 || int_cvs1[ka]->estimatedCurveLength() < lenlim)
	{
	  int_cvs1.erase(int_cvs1.begin()+ka);
	  int_cvs2.erase(int_cvs2.begin()+ka);
	}
    }
  if (int_cvs1.size() == 0)
    return false;

  // Adjust added common regions
  vector<RevEngRegion*> adj_regs;
  if (common_reg2.size() > 0)
    {
      int state = (surf1->instanceType() == Class_Plane ||
		   surf2->instanceType() == Class_Plane) ? 1 : 2;
      double angtol = 5.0*anglim;
      double tol10 = 10.0*tol;
      vector<RevEngPoint*> near_pts;
      for (size_t kr=0; kr<common_reg2.size(); ++kr)
	{
	  vector<RevEngPoint*> adj_near1, adj_near2;
	  double dummy_min = 0.0;
	  for (size_t kh=0; kh<int_cvs1.size(); ++kh)
	    {
	      double tmin3 = int_cvs1[kh]->startparam();
	      double tmax3 = int_cvs1[kh]->endparam();
	      common_reg2[kr]->getNearPoints(int_cvs1[kh], tmin3, tmax3, distance_,
					     angtol, adj_near1);
	      if (state == 1)
		{
		  adj_near2 =
		    adjacent1_->removeOutOfSurf(adj_near1, tol10,
						angtol, outer1_, dummy_min);
		  near_pts = 
		    adjacent2_->removeOutOfSurf(adj_near2, tol10,
						angtol, outer2_, dummy_min);
		}
	      else
		near_pts = adj_near1;
	    }

	  if ((int)near_pts.size() == common_reg2[kr]->numPoints())
	    adj_regs.push_back(common_reg2[kr]);
	  else if ((int)near_pts.size() < common_reg2[kr]->numPoints())
	    {
	      int num_init = common_reg2[kr]->numPoints();
	      vector<vector<RevEngPoint*> > out_groups;
	      vector<HedgeSurface*> out_sfs;
	      vector<vector<RevEngPoint*> > near_groups;
	      common_reg2[kr]->extractSpesPoints(near_pts, near_groups);
	      common_reg2[kr]->updateInfo();
	      common_reg2[kr]->splitRegion(extract_groups);
	      if (common_reg2[kr]->hasSurface() &&
		  common_reg2[kr]->numPoints() < num_init/2)
		{
		  int num_sf = common_reg2[kr]->numSurface();
		  for (int ka=0; ka<num_sf; ++ka)
		    out_sfs.push_back(common_reg2[kr]->getSurface(ka));
		  common_reg2[kr]->clearSurface();
		}

	      // Make new region
	      shared_ptr<RevEngRegion> curr_adj(new RevEngRegion(common_reg2[kr]->getClassificationType(),
								 common_reg2[kr]->getEdgeClassificationType(),
								 near_pts));
	      added_regions.push_back(curr_adj);
	      adj_regs.push_back(curr_adj.get());
	    }
	      
	}
    }

  // Update current edge
  replaceCurves(int_cvs1, int_cvs2);
  if (adj_regs.size() > 0)
    addBlendRegions(adj_regs);
  for (size_t kr=0; kr<adj_regs.size(); ++kr)
    adj_regs[kr]->setAssociatedBlend(this);
  extendcount_ = 0;
  
  return true;
}
  

//===========================================================================
void
RevEngEdge::setTrimCurves(double tol, double angtol,
			  vector<RevEngRegion*>& out_regs,
			  vector<HedgeSurface*>& out_sfs)
//===========================================================================
{
  if (blend_type_ != NOT_BLEND)
    return;

  // Ensure parameter curve existence
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    if (!cvs1_[ki]->hasParameterCurve())
      cvs1_[ki]->ensureParCrvExistence(tol);
  
  for (size_t ki=0; ki<cvs2_.size(); ++ki)
    if (!cvs2_[ki]->hasParameterCurve())
      cvs2_[ki]->ensureParCrvExistence(tol);

  // Distribute points along the curve as appropriate
  // Start with the inbetween points
#ifdef DEBUG_BLEND
  std::ofstream of1("init_points.g2");
  adjacent1_->writeRegionPoints(of1);
  adjacent2_->writeRegionPoints(of1);
  for (size_t kj=0; kj<blend_regs_.size(); ++kj)
    blend_regs_[kj]->writeRegionPoints(of1);
  for (size_t kj=0; kj<cvs1_.size(); ++kj)
    {
      cvs1_[kj]->spaceCurve()->writeStandardHeader(of1);
      cvs1_[kj]->spaceCurve()->write(of1);
    }
#endif
  for (int ka=(int)blend_regs_.size()-1; ka>=0; --ka)
    {
      vector<RevEngPoint*> points = blend_regs_[ka]->getPoints();
      vector<RevEngPoint*> move1, move2;
      adjacent1_->sortBlendPoints(points, cvs1_, 0.0, adjacent2_,
				  move1, move2);
      blend_regs_[ka]->removePoints(move1);
      adjacent1_->addPointsToGroup(move1, tol, angtol);
      adjacent2_->addPointsToGroup(move2, tol, angtol);
      blend_regs_[ka]->removePoints(move2);
      if (blend_regs_[ka]->hasSurface())
	{
	  int num_sf = blend_regs_[ka]->numSurface();
	  for (int kb=0; kb<num_sf; ++kb)
	    out_sfs.push_back(blend_regs_[ka]->getSurface(kb));
	  blend_regs_[ka]->clearSurface();
	}
      
      if (blend_regs_[ka]->numPoints() == 0)
	{
	  blend_regs_[ka]->removeFromAdjacent();
	  blend_regs_[ka]->clearRegionAdjacency();
	  out_regs.push_back(blend_regs_[ka]);
	  blend_regs_.erase(blend_regs_.begin()+ka);
	}
      int stop_break0 = 1;
    }

  vector<RevEngPoint*> move_adj1, move_adj2;
  bool do_move = false;
  if (do_move)
    {
      adjacent1_->extractOutOfEdge2(cvs1_, tol, angtol, move_adj1);
      adjacent2_->extractOutOfEdge2(cvs2_, tol, angtol, move_adj2);
    }
#ifdef DEBUG_BLEND
  std::ofstream of2("move_points.g2");
  adjacent1_->writeRegionPoints(of2);
  if (move_adj1.size() > 0)
    {
      of2 << "400 1 0 4 255 0 0 255" << std::endl;
      of2 << move_adj1.size() << std::endl;
      for (size_t kj=0; kj<move_adj1.size(); ++kj)
	of2 << move_adj1[kj]->getPoint() << std::endl;
    }
  adjacent2_->writeRegionPoints(of2);
  if (move_adj2.size() > 0)
    {
      of2 << "400 1 0 4 0 255 0 255" << std::endl;
      of2 << move_adj2.size() << std::endl;
      for (size_t kj=0; kj<move_adj2.size(); ++kj)
	of2 << move_adj2[kj]->getPoint() << std::endl;
    }
  for (size_t kj=0; kj<cvs1_.size(); ++kj)
    {
      cvs1_[kj]->spaceCurve()->writeStandardHeader(of2);
      cvs1_[kj]->spaceCurve()->write(of2);
    }
#endif
  if (move_adj1.size() > 0)
    adjacent2_->addPointsToGroup(move_adj1, tol, angtol);
  if (move_adj2.size() > 0)
    adjacent1_->addPointsToGroup(move_adj2, tol, angtol);
  
  // Define trim curves
  int stat = 0;
  for (size_t ki=0; ki<cvs1_.size(); ++ki)
    {
      shared_ptr<ftEdge> edg1(new ftEdge(adjacent1_->getSurface(0),
					 cvs1_[ki], cvs1_[ki]->startparam(),
					 cvs1_[ki]->endparam()));
      shared_ptr<ftEdge> edg2(new ftEdge(adjacent2_->getSurface(0),
					 cvs2_[ki], cvs2_[ki]->startparam(),
					 cvs2_[ki]->endparam()));
      edg2->setReversed(true);
      edg1->connectTwin(edg2.get(), stat);
      adjacent1_->addTrimEdge(edg1);
      adjacent2_->addTrimEdge(edg2);
    }

  int stop_break = 1;
}


//===========================================================================
bool
RevEngEdge::contains(RevEngEdge *other, double tol)
//===========================================================================
{
  if (!((adjacent1_ == other->adjacent1_ && adjacent2_ == other->adjacent2_) ||
	(adjacent1_ == other->adjacent2_ && adjacent2_ == other->adjacent1_)))
    return false;  // Not the same adjacent surfaces (groups)

  // Assumes only one curve associated to the edge
  Identity ident;
  int stat = ident.identicalCvs(cvs1_[0], other->cvs1_[0], tol);
  if (stat == 1 || stat == 3)
    return true;
  else
    return false;
}

//===========================================================================
bool
RevEngEdge::integrate(RevEngEdge *other)
//===========================================================================
{
  if (!((adjacent1_ == other->adjacent1_ && adjacent2_ == other->adjacent2_) ||
	(adjacent1_ == other->adjacent2_ && adjacent2_ == other->adjacent1_)))
    return false;  // Not the same adjacent surfaces (groups)

  if (defined_blend_ && defined_blend_ != other->defined_blend_)
    return false;

  for (size_t ki=0; ki<other->blend_regs_.size(); ++ki)
    other->blend_regs_[ki]->setAssociatedBlend(this);

  other->clearBlendRegions();
  
  return true;
}
