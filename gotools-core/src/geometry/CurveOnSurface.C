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

#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/ElementaryCurve.h"
#include "GoTools/geometry/BoundedCurve.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/creators/TrimCurve.h"
#include "GoTools/creators/HermiteAppS.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include <fstream>
#include <cassert>

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::endl;
using std::streamsize;
using std::numeric_limits;
using CoonsPatchGen::blendcoef;

//===========================================================================
CurveOnSurface::CurveOnSurface()
  : prefer_parameter_(true), ccm_(0), constdir_(0),
    constval_(0.0), at_bd_(-1), same_orientation_(true), 
    fix_performed_(0)
//===========================================================================
{
}


//===========================================================================
CurveOnSurface::CurveOnSurface(shared_ptr<ParamSurface> surf,
			       shared_ptr<ParamCurve> curve,
			       bool preferparameter)
  : surface_(surf), ccm_(0), constdir_(0), constval_(0.0), 
    at_bd_(-1), same_orientation_(true), fix_performed_(0)
//===========================================================================
{
  ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");

  ALWAYS_ERROR_IF(curve.get() == 0,"Missing curve.");
  ALWAYS_ERROR_IF((preferparameter && curve->dimension() != 2) ||
		  (!preferparameter && (curve->dimension() != surf->dimension())),
		  "Conflict in dimension of geometric objects");


  prefer_parameter_ = preferparameter;
  if (preferparameter)
  {
    pcurve_ = curve;
    spacecurve_ = shared_ptr<ParamCurve>();
  }
  else
  {
    pcurve_ = shared_ptr<ParamCurve>();
    spacecurve_ = curve;
  }
}

//===========================================================================
CurveOnSurface::CurveOnSurface(shared_ptr<ParamSurface> surf,
			       shared_ptr<ParamCurve> curve,
			       int constdir, double constpar, int boundary)
  : surface_(surf), ccm_(3), constdir_(constdir),
    constval_(constpar), at_bd_(boundary), same_orientation_(true),
    fix_performed_(0)
//===========================================================================
{
  ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");

  ALWAYS_ERROR_IF(curve.get() == 0,"Missing curve.");
  ALWAYS_ERROR_IF(curve->dimension() != surf->dimension(),
		  "Conflict in dimension of geometric objects");


  prefer_parameter_ = false;
  pcurve_ = shared_ptr<ParamCurve>();
  spacecurve_ = curve;
  double t1 = startparam();
  double t2 = endparam();
  Point pt1 = faceParameter(t1);
  Point pt2 = faceParameter(t2);

  // Check for orientation
  Point pnt1 = curve->point(t1);
  Point pnt2 = curve->point(t2);
  Point pnt3 = surf->point(pt1[0], pt1[1]);
  Point pnt4 = surf->point(pt2[0], pt2[1]);
  if (pnt1.dist(pnt3) + pnt2.dist(pnt4) > pnt1.dist(pnt4) + pnt2.dist(pnt3))
    {
      same_orientation_ = false;
      std::swap(pt1, pt2);
    }

  pcurve_ = shared_ptr<ParamCurve>(new SplineCurve(pt1, t1, pt2, t2));
}
//===========================================================================
CurveOnSurface::CurveOnSurface(shared_ptr<ParamSurface> surf,
			       int constdir, double constpar, 
			       double par1, double par2, int boundary)
  : surface_(surf), ccm_(3), constdir_(constdir),
    constval_(constpar), at_bd_(boundary), fix_performed_(0)
//===========================================================================
{
  ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");


  pcurve_ = shared_ptr<ParamCurve>();
  spacecurve_ = shared_ptr<ParamCurve>();

  vector<shared_ptr<ParamCurve> > cvs = surf->constParamCurves(constpar, 
							       constdir == 2);
  size_t ki;
  double tol = 1.0e-10;
  double t1 = std::min(par1, par2);
  double t2 = std::max(par1, par2);
  for (ki=0; ki<cvs.size(); ++ki)
    if (cvs[ki]->startparam() < t1+tol &&
	cvs[ki]->endparam() > t2-tol)
      break;
  if (ki < cvs.size())
    {
      spacecurve_ = shared_ptr<ParamCurve>(cvs[ki]->subCurve(t1, t2));
      if (par2 < par1)
	spacecurve_->reverseParameterDirection();
    }

  Point pt1, pt2;
  if (constdir == 1)
    {
      pt1 = Point(constpar, par1);
      pt2 = Point(constpar, par2);
    }
  else
    {
      pt1 = Point(par1, constpar);
      pt2 = Point(par2, constpar);
    }

  // Check for orientation
  same_orientation_ = (par2 > par1);

  pcurve_ = shared_ptr<ParamCurve>(new SplineCurve(pt1, par1, pt2, par2));
  prefer_parameter_ = (spacecurve_.get()) ? true : false;
}

//===========================================================================
CurveOnSurface::CurveOnSurface(shared_ptr<ParamSurface> surf,
				   shared_ptr<ParamCurve> pcurve,
                                   shared_ptr<ParamCurve> spacecurve,
			       bool preferparameter, int ccm)
    : surface_(surf), pcurve_(pcurve), spacecurve_(spacecurve),
      prefer_parameter_(preferparameter), ccm_(ccm), 
      constdir_(0), constval_(0.0), at_bd_(-1), same_orientation_(true),
      fix_performed_(0)
//===========================================================================
{
  ALWAYS_ERROR_IF(surf.get() == 0, "Missing surface.");
  if (preferparameter) {
      ALWAYS_ERROR_IF(pcurve.get() == 0, "Missing curve.");
  } else {
      ALWAYS_ERROR_IF(spacecurve.get() == 0, "Missing curve.");
  }
  ALWAYS_ERROR_IF((pcurve.get() != 0) && (pcurve->dimension() != 2),
		  "The parametric curve must have dimension 2");

  ALWAYS_ERROR_IF((spacecurve.get() != 0)
	      && (spacecurve->dimension() != surf->dimension()),
		  "Conflict in dimension of geometric objects");


  // @afr: Deactivated this check before committing to CVS.
  // It should be activated, but perhaps not until other
  // code is able to deal with the consequences. For example,
  // code in the IGES reader may have to turn curves or
  // reorder them in order to make proper loops.
  // @afr 2004-10-05: Put it back in. It solves the HSDF04.igs case...
#if 0 //1
  if ((spacecurve.get() != 0) && (pcurve.get() != 0)) {
      Point p1_1 = spacecurve->point(spacecurve->startparam());
      Point temp = pcurve->point(pcurve->startparam());
      Point p1_2 = surf->point(temp[0],temp[1]);
      double dist1 = p1_1.dist(p1_2);
      Point p2_1 = spacecurve->point(spacecurve->endparam());
      temp = pcurve->point(pcurve->endparam());
      Point p2_2 = surf->point(temp[0],temp[1]);
      double dist2 = p2_1.dist(p2_2);
      double releps = 1e-3;
      double abseps = 1e-6;
      double eps = max(releps * p1_1.dist(p2_1), abseps);
      if (dist1 > eps || dist2 > eps) {
	  double ddist1 = p1_1.dist(p2_2);
	  double ddist2 = p2_1.dist(p1_2);
	  if (ddist1 > eps || ddist2 > eps) {
	      cerr << "dist1 = " << dist1 << "  dist2 = " << dist2
		   << "  ddist1 = " << ddist1 << "  ddist2 = " << ddist2
		   << "     eps = " << eps << endl;
	      cerr << p1_1 << p2_1 << p1_2 << p2_2 << endl;
	      THROW("Inconsistent parametric and spatial curves.");
		    
	  }
	  cerr << "Fixed inconsistency between space and parameter curve!"
	       << endl;
	  // Which to turn? Depends on the preference.
	  if (prefer_parameter_) {
	      spacecurve->reverseParameterDirection();
	  } else {
	      pcurve->reverseParameterDirection();
	  }
      }
  }
#endif
}

//===========================================================================
CurveOnSurface::CurveOnSurface(const CurveOnSurface& surface_curve)
  : surface_(surface_curve.surface_), 
    prefer_parameter_(surface_curve.prefer_parameter_),
    ccm_(surface_curve.ccm_), constdir_(surface_curve.constdir_),
    constval_(surface_curve.constval_), at_bd_(surface_curve.at_bd_),
    same_orientation_(surface_curve.same_orientation_), fix_performed_(0) 
//===========================================================================
{
  // Clones the curves, not the surface
  if (surface_curve.pcurve_.get() != 0)
    {
// #ifdef _MSC_VER
//       ParamCurve *pcv
// 	= dynamic_cast<ParamCurve*>(surface_curve.pcurve_->clone());
// #else
      ParamCurve *pcv = surface_curve.pcurve_->clone();
// #endif
//        std::cout << "Par crv: " << pcv;
      pcurve_  = shared_ptr<ParamCurve>(pcv);
    }
  if (surface_curve.spacecurve_.get() != 0)
      {
// #ifdef _MSC_VER
// 	  ParamCurve *scv
// 	      = dynamic_cast<ParamCurve*>(surface_curve.spacecurve_->clone());
// #else
	  ParamCurve *scv = surface_curve.spacecurve_->clone();
// #endif
	  //    std::cout << ". Geom crv: " << scv << std::endl;
	  spacecurve_ = shared_ptr<ParamCurve>(scv);
      }
}

//===========================================================================
CurveOnSurface::CurveOnSurface(shared_ptr<ParamSurface> surf,
			       shared_ptr<ParamCurve> pcurve,
			       shared_ptr<ParamCurve> spacecurve,
			       bool preferparameter, int ccm,
			       int constdir, double constpar, int boundary,
			       bool same_orientation)
    : surface_(surf), pcurve_(pcurve), spacecurve_(spacecurve),
      prefer_parameter_(preferparameter), ccm_(ccm), 
      constdir_(constdir), constval_(constpar), at_bd_(boundary), 
      same_orientation_(same_orientation), fix_performed_(0)
//===========================================================================
{
}

//===========================================================================
CurveOnSurface& CurveOnSurface::operator=(const CurveOnSurface& other)
//===========================================================================
{
  if (&other != this) {
    // Clones the curves, not the surface
    surface_ = other.surface_;
    prefer_parameter_ = other.prefer_parameter_;
    constdir_ = other.constdir_;
    constval_ = other.constval_;
    at_bd_ = other.at_bd_;
    same_orientation_ = other.same_orientation_;
    if (other.pcurve_.get() != 0)
      {
// #ifdef _MSC_VER
// 	ParamCurve *pcv
// 	  = dynamic_cast<ParamCurve*>(other.pcurve_->clone());
// #else
	ParamCurve *pcv = other.pcurve_->clone();
// #endif
	//	std::cout << "Par crv: " << pcv;
	pcurve_  = shared_ptr<ParamCurve>(pcv);
      }
    else
      pcurve_ = shared_ptr<ParamCurve>();
    if (other.spacecurve_.get() != 0)
	{
// #ifdef _MSC_VER
//     ParamCurve *scv
//       = dynamic_cast<ParamCurve*>(other.spacecurve_->clone());
// #else
	  ParamCurve *scv = other.spacecurve_->clone();
// #endif
    //    std::cout << ". Geom crv: " << scv << std::endl;
	  spacecurve_ = shared_ptr<ParamCurve>(scv);
	}
    else
      spacecurve_ = shared_ptr<ParamCurve>();

  }
  return *this;
}

//===========================================================================
CurveOnSurface::~CurveOnSurface()
//===========================================================================
{
}


//===========================================================================
void CurveOnSurface::read(std::istream& is)
//===========================================================================
{
//    streamsize prev = is.precision(15);

    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    // Do not care about surface...
    ALWAYS_ERROR_IF(pcurve_.get() != NULL,
		    "Parameter curve already exists!");

    ALWAYS_ERROR_IF(spacecurve_.get() != NULL,
		    "Space curve already exists!");

    bool prefer_parameter;
    int  prefer_parameter_int;
    int  pcurve_type;
    int  spacecurve_type;
    shared_ptr<ParamCurve> pcurve;
    shared_ptr<ParamCurve> spacecurve;

    is >> prefer_parameter_int;
    if (prefer_parameter_int == 0)
	prefer_parameter = false;
    else if (prefer_parameter_int == 1)
	prefer_parameter = true;
    else 
	THROW("Unknown input for preferred CurveOnSurface parameter");

    is >> pcurve_type;
    is >> spacecurve_type;

    if (pcurve_type == 0) {
	// Do nothing - continue
    }
    else if (pcurve_type == 1) {
	// Interpreting pcurve_type = 1 as SplineCurve for backwards
	// compatibility, and give warning. We should really use the
	// ClassType of the relevant ParamCurve.
// 	MESSAGE("Warning: Read ClassType '1'.\n"
// 		"Interpreting this as a SplineCurve for backward "
// 		"compatibility.\n"
// 		"Continuing...");
	pcurve = shared_ptr<ParamCurve>(new SplineCurve());
	pcurve->read(is);
    }
    else {
	ClassType type = ClassType(pcurve_type); // Needs this conversion
	shared_ptr<GeomObject> goobject(Factory::createObject(type));
	pcurve = dynamic_pointer_cast<ParamCurve, GeomObject>(goobject);
	ALWAYS_ERROR_IF(pcurve.get() == 0,
			"Can not read this instance type");
	pcurve->read(is);
    }
    if (spacecurve_type == 0) {
	// Do nothing - continue
    }
    else if (spacecurve_type == 1) {
	// Interpreting spacecurve_type = 1 as SplineCurve for backwards
	// compatibility, and give warning. We should really use the
	// ClassType of the relevant ParamCurve.
// 	MESSAGE("Warning: Read ClassType '1'.\n"
// 		"Interpreting this as a SplineCurve for backward "
// 		"compatibility.\n"
// 		"Continuing...");
	spacecurve = shared_ptr<ParamCurve>(new SplineCurve());
	spacecurve->read(is);
    }
    else {
	ClassType type = ClassType(spacecurve_type); // Needs this conversion
	shared_ptr<GeomObject> goobject(Factory::createObject(type));
	spacecurve = dynamic_pointer_cast<ParamCurve, GeomObject>(goobject);
	ALWAYS_ERROR_IF(spacecurve.get() == 0,
			"Can not read this instance type");
	spacecurve->read(is);
    }

    prefer_parameter_ = prefer_parameter;
    pcurve_ = pcurve;
    spacecurve_ = spacecurve;

#ifndef NDEBUG
    if ((pcurve_.get() != NULL) && (spacecurve_.get() != NULL))
    {
	double partol =1e-12;
	double diff_start = fabs(pcurve_->startparam() - spacecurve_->startparam());
	double diff_end = fabs(pcurve_->endparam() - spacecurve_->endparam());
	if ((diff_start > partol) || (diff_end > partol))
	{
	    double max_diff = std::max(diff_start, diff_end);
	    MESSAGE("CurveOnSurface::read(): End parameters not tol-equal (" << max_diff << " > 1e-12)!");
	}
    }
#endif

    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

//    is.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
void CurveOnSurface::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    // Currently do not write the surface...
    if (!prefer_parameter_)
	os << "0";
    else
	os << "1";
    os << ' ';

    if (pcurve_.get() == NULL)
	os << "0";
    else
	os << pcurve_->instanceType();
    os << ' ';

    if (spacecurve_.get() == NULL)
	os << "0";
    else
	os << spacecurve_->instanceType();
    os << endl;


    if (pcurve_.get() != NULL)
	pcurve_->write(os);

    if (spacecurve_.get() != NULL)
	spacecurve_->write(os);

    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
BoundingBox CurveOnSurface::boundingBox() const
//===========================================================================
{
    if (spacecurve_.get() != 0)
      return spacecurve_->boundingBox();
    else
      return surface_->boundingBox();
}

//===========================================================================
DirectionCone CurveOnSurface::directionCone() const
//===========================================================================
{
    if (spacecurve_.get() != 0)
      return spacecurve_->directionCone();
    else
      {
	// @@@ VSK 0904. Returns an empty cone for the time being
	DirectionCone cone;
	return cone;
      }
}

//===========================================================================
RectDomain CurveOnSurface::containingDomain() const
//===========================================================================
{
    if (pcurve_.get() != 0) {
	BoundingBox bb = pcurve_->boundingBox();
	return RectDomain(Vector2D(bb.low()[0], bb.low()[1]),
			  Vector2D(bb.high()[0], bb.high()[1]));
    }
    else
	return surface_->containingDomain();
}


//===========================================================================
int CurveOnSurface::dimension() const
//===========================================================================
{
    return surface_->dimension();
}


//===========================================================================
ClassType CurveOnSurface::instanceType() const
//===========================================================================
{
    return classType();
}




//===========================================================================
double CurveOnSurface::startparam() const
//===========================================================================
{
    if (prefer_parameter_ && (pcurve_ != NULL))
    return pcurve_->startparam();
  else
    return spacecurve_->startparam();
}


//===========================================================================
double CurveOnSurface::endparam() const
//===========================================================================
{
  if (prefer_parameter_ && (pcurve_ != NULL))
    return pcurve_->endparam();
  else
    return spacecurve_->endparam();
}

//===========================================================================
SplineCurve* CurveOnSurface::geometryCurve()
//===========================================================================
{
  if (spacecurve_.get() != 0)
    return spacecurve_->geometryCurve();
  else {
    double tol = 1.0e-4; // Arbitrary. Do rather call 
    // ensureSpaceCrve from the application
    bool done = ensureSpaceCrvExistence(tol);
    if (done)
      return spacecurve_->geometryCurve();
    else
      {
	MESSAGE("Could not compute space curve");
	return NULL;
      }
  }
}


//===========================================================================
void CurveOnSurface::reverseParameterDirection(bool switchparam)
//===========================================================================
{
  if (pcurve_.get() != 0)
    pcurve_->reverseParameterDirection(switchparam);
  if (spacecurve_.get()!= 0)
    spacecurve_->reverseParameterDirection();
  if (constdir_)
    {
      same_orientation_ = (same_orientation_) ? false : true;

      if (switchparam)
	{
	  constdir_ = 3 - constdir_;
	  if (at_bd_ >= 0)
	    {
	      int sgn = (at_bd_ <= 1) ? 1 : -1;
	      at_bd_ += sgn*2;
	    }
	      
	}
    }
}


//===========================================================================
void CurveOnSurface::setParameterInterval(double t1, double t2)
//===========================================================================
{
  if (pcurve_.get() != 0)
    pcurve_->setParameterInterval(t1, t2);
  if (spacecurve_.get()!= 0)
    spacecurve_->setParameterInterval(t1, t2);
}

//===========================================================================
bool CurveOnSurface::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    if (prefer_parameter_) {
	// @@ Simple approach: Evaluate start, mid and endpoints
	Point start, mid, end;
	point(start, startparam());
	point(mid, 0.5*startparam() + 0.5*endparam());
	point(end, endparam());
	double len = start.dist(mid) + mid.dist(end);
	if (len > degenerate_epsilon)
	    return false;
	else
	    return true;
    } else {
	return spacecurve_->isDegenerate(degenerate_epsilon);
    }
}

//===========================================================================
void CurveOnSurface::point(Point& pt, double tpar) const
//===========================================================================
{
//    std::cout << "Curve on surface: " << this;
//    std::cout << ", Space curve: " << spacecurve_.get();
//    std::cout << ",count: " << *(spacecurve_.pn) << std::endl;
    if (prefer_parameter_) {
	Point param_pt(2);
	pcurve_->point(param_pt, tpar);
	surface_->point(pt, param_pt[0], param_pt[1]);
    } else {
	spacecurve_->point(pt, tpar);
    }
}




//===========================================================================
double  CurveOnSurface::nextSegmentVal(double par, bool forward, double tol) const
//===========================================================================
{
    if (prefer_parameter_)
    {
	if (fabs(startparam() - pcurve_->startparam()) < tol &&
	    fabs(endparam() - pcurve_->endparam()) < tol)
	    return pcurve_->nextSegmentVal(par, forward, tol);
	else
	    return (forward) ? endparam() : startparam();
    }
    else
	return spacecurve_->nextSegmentVal(par, forward, tol);
}

//===========================================================================
void CurveOnSurface::point(std::vector<Point>& pts,
			     double tpar,
			     int derivs,
			     bool from_right) const
//===========================================================================
{
  if (prefer_parameter_)
  {
    // @@@
    ALWAYS_ERROR_IF(derivs > 2, "Only two derivatives supported.");
    std::vector<Point> param_pts(derivs+1, Point(2));
    std::vector<Point> surf_pts((derivs+1)*(derivs+2)/2,
				  Point(surface_->dimension()));
    pcurve_->point(param_pts, tpar, derivs, from_right);
    surface_->point(surf_pts, param_pts[0][0], param_pts[0][1], derivs);
    pts.clear();
    pts.push_back(surf_pts[0]);
    if (derivs > 0) {
	pts.push_back(surf_pts[1]*param_pts[1][0] + surf_pts[2]*param_pts[1][1]);
    }
    if (derivs > 1) {
	// add second derivative
	Point& p_dt   = param_pts[1];
	Point& p_dtdt = param_pts[2];
	Point& s_du   = surf_pts[1];
	Point& s_dv   = surf_pts[2];
	Point& s_dudu = surf_pts[3];
	Point& s_dudv = surf_pts[4];
	Point& s_dvdv = surf_pts[5];
	Point double_der = 
	    s_du * p_dtdt[0] +
	    s_dv * p_dtdt[1] +
	    s_dudu * p_dt[0] * p_dt[0] +
	    s_dudv * p_dt[0] * p_dt[1] + 
	    s_dvdv * p_dt[1] * p_dt[1];
	pts.push_back(double_der);
    }
  }
  else
    spacecurve_->point(pts, tpar, derivs, from_right);
}


namespace {
    template <typename PairType>
    struct comparepair_second
    {
	bool operator() (const PairType& p1, const PairType& p2)
	{
	    return p1.second < p2.second;
	}
    };
}

//===========================================================================
void CurveOnSurface::closestPoint(const Point&   pt,
				  double           tmin,
				  double           tmax,
				  double&          clo_t,
				  Point&         clo_pt,
				  double&          clo_dist,
				  double const     *seed) const
//===========================================================================
{
    int ki;
    //     // Find some initial guess values for the parameter, and
    //     // pass those to the generic solver.
    //     // We do that by evaluating the distance at the interval endpoints
    //     // and the interval midpoint (simple, may fail).

    double guess_param;
    if (seed != 0) {
	guess_param = *seed;
	if (guess_param < tmin || guess_param > tmax) {
	    MESSAGE("Suggested parameter for closest point "
		       "must lie inside domain!");
	    guess_param = (tmax < guess_param ? tmax : guess_param);
	    guess_param = (tmin > guess_param ? tmin : guess_param);
	}
    } else {
	// Find some initial guess values for the parameter, and
	// pass those to the generic solver.
        // We evaluate in nbm_ctl_pts params.
	std::vector<std::pair<double, double> > par_and_dist;
	// More values may be push_backed here in order to make a
	// better starting guess
	int nmb_sample_pts = 3;
	if ((pcurve_.get() != 0) && (pcurve_->instanceType() == Class_SplineCurve)) {
	    shared_ptr<SplineCurve> tempsc(dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve_));
	    nmb_sample_pts = max(nmb_sample_pts, tempsc->numCoefs());
	}
	if ((spacecurve_.get() != 0) && (spacecurve_->instanceType() == Class_SplineCurve)) {
	    nmb_sample_pts = max(nmb_sample_pts,
				 (dynamic_pointer_cast<SplineCurve, ParamCurve>(spacecurve_))->numCoefs());
	}
	par_and_dist.reserve(nmb_sample_pts);
	double tstep = (tmax - tmin)/(nmb_sample_pts - 1);
	for (ki = 0; ki < nmb_sample_pts; ++ki) {
	    double tpar = tmin + ki*tstep;
	    point(clo_pt, tpar);
	    double dist = pt.dist(clo_pt);
	    par_and_dist.push_back(std::make_pair(tpar, dist));
	}
	std::sort(par_and_dist.begin(), par_and_dist.end(),
		  comparepair_second< std::pair<double, double> >());
	guess_param = par_and_dist[0].first;
    }

    ParamCurve::closestPointGeneric(pt, tmin, tmax, guess_param, 
				    clo_t, clo_pt, clo_dist);

}


//===========================================================================
double CurveOnSurface::length(double tol)
//===========================================================================
{
    shared_ptr<SplineCurve> geom_cv(geometryCurve());
    assert(geom_cv.get() != NULL);
    return geom_cv->length(tol);
}


//===========================================================================
void CurveOnSurface::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    // We're assuming C1 as default.
    int cont = 1;
    double dist_dummy;
    appendCurve(cv, cont, dist_dummy, reparam);
}


//===========================================================================
void CurveOnSurface::appendCurve(ParamCurve* other_curve,
				   int continuity, double& dist, bool reparam)
//===========================================================================
{
    CurveOnSurface* other_cv = dynamic_cast<CurveOnSurface*>(other_curve);
    ALWAYS_ERROR_IF(other_cv == 0,
		"Argument cv was of wrong type!");
    ALWAYS_ERROR_IF(surface_.get() != (other_cv->underlyingSurface()).get(),
		    "Trying to append a curveOnSurface attached to another surface.");


    shared_ptr<ParamCurve> other_pcurve = other_cv->parameterCurve();
    shared_ptr<ParamCurve> other_spacecurve = other_cv->spaceCurve();

    // Either must both param curves exist, or both space curves.
    if ((pcurve_.get() == 0) || (other_pcurve.get() == 0))
	pcurve_.reset();
    if ((spacecurve_.get() == 0) || (other_spacecurve.get() == 0))
	spacecurve_.reset();

    ALWAYS_ERROR_IF((pcurve_.get() == 0) && (spacecurve_.get() == 0),
		    "Mismatch between curve represetations of the two objects.");


    double tol = 1.0e-4;
    if (prefer_parameter_)
      {
	try {
	  pcurve_->appendCurve(other_pcurve.get(), continuity, dist, reparam);
	}
	catch (...)
	  {
	    shared_ptr<SplineCurve> tmp1 = 
	      shared_ptr<SplineCurve>(pcurve_->geometryCurve());
	    shared_ptr<SplineCurve> tmp2 = 
	      shared_ptr<SplineCurve>(other_pcurve->geometryCurve());
	    tmp1->appendCurve(tmp2.get(), continuity, dist, reparam);
	    pcurve_ = tmp1;
	  }
	spacecurve_.reset();
	ensureSpaceCrvExistence(tol);
      }
    else 
      {
	double pardist;
	try {
	  spacecurve_->appendCurve(other_spacecurve.get(), continuity, dist, reparam);
	}
	catch (...)
	  {
	    shared_ptr<SplineCurve> tmp1 = 
	      shared_ptr<SplineCurve>(spacecurve_->geometryCurve());
	    shared_ptr<SplineCurve> tmp2 = 
	      shared_ptr<SplineCurve>(other_spacecurve->geometryCurve());
	    tmp1->appendCurve(tmp2.get(), continuity, dist, reparam);
	    spacecurve_ = tmp1;
	  }

#ifdef DEBUG
	std::ofstream of("par_crvs.g2");
	pcurve_->writeStandardHeader(of);
	pcurve_->write(of);
	other_pcurve->writeStandardHeader(of);
	other_pcurve->write(of);
#endif

	if (continuity < 1 && (!reparam))
	  {
	    try {
	      pcurve_->appendCurve(other_pcurve.get(), continuity, pardist, reparam);
	    }
	    catch (...)
	      {
		shared_ptr<SplineCurve> tmp1 = 
		  shared_ptr<SplineCurve>(pcurve_->geometryCurve());
		shared_ptr<SplineCurve> tmp2 = 
		  shared_ptr<SplineCurve>(other_pcurve->geometryCurve());
		tmp1->appendCurve(tmp2.get(), continuity, dist, reparam);
		pcurve_ = tmp1;
	      }
	  }
	else
	  {
	    Point par1 = pcurve_->point(pcurve_->startparam());
	    Point par2 = pcurve_->point(pcurve_->endparam());
	    Point par3 = other_pcurve->point(other_pcurve->startparam());
	    Point par4 = other_pcurve->point(other_pcurve->endparam());

	    // // Adjust parameter and tolerance
	    // Point pos1 = spacecurve_->point(spacecurve_->startparam());
	    // Point pos2 = spacecurve_->point(spacecurve_->endparam());
	    // double u1, u2, v1, v2, d1, d2;
	    // Point close1, close2;
	    // surface_->closestPoint(pos1, u1, v1, close1, d1, tol, NULL, par1.begin());
	    // surface_->closestPoint(pos2, u2, v2, close2, d2, tol, NULL, par4.begin());
	    // par1 = Point(u1,v1);
	    // par4 = Point(u2,v2);
	    // tol = std::max(tol, std::max(d1,d2));

	    pcurve_.reset();
	    if (false /*reparam*/)
	      ensureParCrvExistence(tol);
	    else
	      makeParameterCurve(tol, par1, par4);
	  }
#ifdef DEBUG
	pcurve_->writeStandardHeader(of);
	pcurve_->write(of);
#endif
      }
    // We do not alter value of prefer_parameter_.

    // @@@ Should perform curve-curve-intersection test to see if created curve
    // is inside legal boundaries.
#ifdef GEOMETRY_DEBUG
    MESSAGE("No test performed to check whether created curve is inside domain.");
#endif // GEOMETRY_DEBUG

}



//===========================================================================
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
// ParamCurve* CurveOnSurface::subCurve(double from_par, 
// 					 double to_par, double fuzzy) const
// #else
CurveOnSurface* CurveOnSurface::subCurve(double from_par, 
					     double to_par, double fuzzy) const
// #endif
//===========================================================================
{
    shared_ptr<ParamCurve> subpcurve;
    shared_ptr<ParamCurve> subspacecurve;

    if (prefer_parameter_) {
	if (pcurve_.get() != 0) {
	    subpcurve = shared_ptr<ParamCurve>(pcurve_->subCurve(from_par, to_par,
								   fuzzy));
	    if (spacecurve_.get() != 0) {
		// Check if parametric and spatial curve are reasonably
		// co-parametrized.
		Point from_par_pt = pcurve_->point(from_par);
		Point to_par_pt = pcurve_->point(to_par);
		Point from_pt = surface_->point(from_par_pt[0], from_par_pt[1]);
		Point to_pt = surface_->point(to_par_pt[0], to_par_pt[1]);
		Point guess_from_sp_pt = spacecurve_->point(from_par);
		Point guess_to_sp_pt = spacecurve_->point(to_par);
		double d1 = guess_from_sp_pt.dist(from_pt);
		double d2 = guess_to_sp_pt.dist(to_pt);
		double dmin = (d1 <= d2 ? d1 : d2);

		// if (dmin > fuzzy) {
		//     MESSAGE("Not well co-parametrized: dmin = " << dmin);
		if (dmin > 1.0e-5) { // Too much noice
		    MESSAGE("Not well co-parametrized: dmin = " << dmin);
		}

		double clo_from, clo_to, clo_dist1, clo_dist2;
		Point clo_pt;
		double seed[1];
		*seed = spacecurve_->startparam() +
		    (spacecurve_->endparam() - spacecurve_->startparam()) *
		    (from_par - pcurve_->startparam()) /
		    (pcurve_->endparam() - pcurve_->startparam());
		spacecurve_->closestPoint(from_pt, spacecurve_->startparam(),
					  spacecurve_->endparam(), clo_from,
					  clo_pt, clo_dist1, &seed[0]);
		*seed = spacecurve_->startparam() +
		    (spacecurve_->endparam() - spacecurve_->startparam()) *
		    (to_par - pcurve_->startparam()) /
		    (pcurve_->endparam() - pcurve_->startparam());
		spacecurve_->closestPoint(to_pt, spacecurve_->startparam(),
					  spacecurve_->endparam(), clo_to,
					  clo_pt, clo_dist2, &seed[0]);
		try {
		    // subspacecurve =
		    // 	shared_ptr<ParamCurve>(spacecurve_->subCurve(clo_from, clo_to,
		    // 						   fuzzy));
		    subspacecurve =
		    	shared_ptr<ParamCurve>(spacecurve_->subCurve(from_par, 
		    						     to_par,
		    						     fuzzy));
		} catch (...) {
		    MESSAGE("Failed extracting spatial part.");
		}
		//subspacecurve->setParameterInterval(from_par, to_par);

		shared_ptr<SplineCurve> tmp_space =
		  dynamic_pointer_cast<SplineCurve,ParamCurve>(subspacecurve);
		if (tmp_space.get())
		  {
		    // The space curve is a spline curve. This allows us
		    // to ensure that the endpoints of the curve lies at
		    // the surface and corresponds to the endpoints of
		    // the parameter curve
		    tmp_space->replaceEndPoint(from_pt, true);
		    tmp_space->replaceEndPoint(to_pt, false);
		  }

	    }
	} else
	    THROW("Missing parametercurve.");
    } else {
	if (spacecurve_.get() != 0) {
	    subspacecurve =
		shared_ptr<ParamCurve>(spacecurve_->subCurve(from_par, to_par,
							       fuzzy));
	    if (pcurve_.get() != 0) {
		Point from_space_pt = spacecurve_->point(from_par);
		Point to_space_pt = spacecurve_->point(to_par);
		double clo_u_from, clo_v_from, clo_u_to, clo_v_to;
		double clo_dist1, clo_dist2;
		Point clo_pt, clo_pt_from, clo_pt_to;
		Point from_seed = pcurve_->point(from_par);
		Point to_seed = pcurve_->point(to_par);
		surface_->closestPoint(from_space_pt, clo_u_from, clo_v_from,
				       clo_pt_from, clo_dist1, fuzzy,
				       NULL, from_seed.begin());
		surface_->closestPoint(to_space_pt, clo_u_to, clo_v_to,
				       clo_pt_to, clo_dist2, fuzzy,
				       NULL, to_seed.begin());
		Point from_par_pt(clo_u_from, clo_v_from);
		Point to_par_pt(clo_u_to, clo_v_to);

		double clo_from, clo_to;
		pcurve_->closestPoint(from_par_pt, pcurve_->startparam(),
				      pcurve_->endparam(), clo_from,
				      clo_pt, clo_dist1,&from_par);
		pcurve_->closestPoint(to_par_pt, pcurve_->startparam(),
				      pcurve_->endparam(), clo_to,
				      clo_pt, clo_dist2, &to_par);
		try {
		    // subpcurve =
		    // 	shared_ptr<ParamCurve>(pcurve_->subCurve(clo_from, clo_to,
		    // 						   fuzzy));
		    subpcurve =
		    	shared_ptr<ParamCurve>(pcurve_->subCurve(from_par, 
		    						 to_par,
		    						   fuzzy));
		} catch (...) {
		    MESSAGE("Failed extracting parametric part.");
		}
		//subpcurve->setParameterInterval(from_par, to_par);

		shared_ptr<SplineCurve> tmp_space =
		  dynamic_pointer_cast<SplineCurve,ParamCurve>(subspacecurve);
		shared_ptr<SplineCurve> tmp_par =
		  dynamic_pointer_cast<SplineCurve,ParamCurve>(subpcurve);
		if (tmp_par.get() && tmp_space.get())
		  {
		    // Both curves are spline curves. This allows us
		    // to ensure that the endpoints of the curve lies at
		    // the surface and corresponds to the endpoints of
		    // the parameter curve
		    tmp_space->replaceEndPoint(clo_pt_from, true);
		    tmp_space->replaceEndPoint(clo_pt_to, false);
		    tmp_par->replaceEndPoint(from_par_pt, true);
		    tmp_par->replaceEndPoint(to_par_pt, false);
		  }
	    }
	} else
	    THROW("Missing spacecurve.");
    }

    CurveOnSurface *sub_cv = new CurveOnSurface(surface_, subpcurve, 
						subspacecurve,
						prefer_parameter_);
    sub_cv->ccm_ = ccm_;
    sub_cv->constdir_ = constdir_;
    sub_cv->constval_ = constval_;
    sub_cv->at_bd_ = at_bd_;
    sub_cv->same_orientation_ = same_orientation_;

    return sub_cv;
}

//===========================================================================
vector<shared_ptr<ParamCurve> >  CurveOnSurface::split(double param,
						       double fuzzy) const
//===========================================================================
{
  double ang_tol = 0.05;  // A rather arbitrary angular tolerance
  vector<shared_ptr<ParamCurve> > pcvs(2);
  vector<shared_ptr<ParamCurve> > spacecvs(2);
  int nmb_replace_param = 0;

  if (prefer_parameter_ && pcurve_.get() != 0) 
    {
      pcvs = pcurve_->split(param, fuzzy);

      if (spacecurve_.get() != 0)
	{
	  spacecvs = spacecurve_->split(param, fuzzy);
	  Point par_pt = pcurve_->point(param);
	  Point sf_pt = surface_->point(par_pt[0], par_pt[1]);

	  // Make sure that the space curve corrsponds with the
	  // parameter curve in the split point.
	  // Should also the tangent be modified to ensure C1 continuity?
	  shared_ptr<SplineCurve> tmp_space =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(spacecvs[0]);
	  if (tmp_space.get())
	    {
	      // The space curve is a spline curve. This allows us
	      // to ensure that the endpoints of the curve lies at
	      // the surface and corresponds to the endpoints of
	      // the parameter curve
	      tmp_space->replaceEndPoint(sf_pt, false);
	    }

	  tmp_space =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(spacecvs[1]);
	  if (tmp_space.get())
	    {
	      // The space curve is a spline curve. This allows us
	      // to ensure that the endpoints of the curve lies at
	      // the surface and corresponds to the endpoints of
	      // the parameter curve
	      tmp_space->replaceEndPoint(sf_pt, true);
	    }
	}
    }
  else if (spacecurve_.get() != 0)
    {
      spacecvs = spacecurve_->split(param, fuzzy);

      if (pcurve_.get() != 0) 
	{
	  pcvs = pcurve_->split(param, fuzzy);
	  
	  Point space = spacecurve_->point(param);
	  Point seed = pcurve_->point(param);

	  double par_u, par_v, sf_dist;
	  Point sf_pt;
	  surface_->closestPoint(space, par_u, par_v,
				 sf_pt, sf_dist, fuzzy,
				 NULL, seed.begin());

	  // Make sure that the prameter curve corrsponds with the
	  // space curve in the split point.
	  // Should also the tangent be modified to ensure C1 continuity?
	  shared_ptr<SplineCurve> tmp_par =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(pcvs[0]);
	  if (tmp_par.get())
	    {
	      // The parameter curve is a spline curve. This allows us
	      // to ensure correspondance in the split point
	      vector<Point> p1(2), p2(2);
	      tmp_par->point(p1, tmp_par->endparam(), 1);
	      tmp_par->replaceEndPoint(Point(par_u, par_v), false);

	      // Check consistence of tangent
	      tmp_par->point(p2, tmp_par->endparam(), 1);
	      if (p1[1].angle(p2[1]) > ang_tol)
		{
		  // Reset curve
		  tmp_par->replaceEndPoint(p1[0], false);
		  nmb_replace_param++;
		}
	    }

	  tmp_par =
	    dynamic_pointer_cast<SplineCurve,ParamCurve>(pcvs[1]);
	  if (tmp_par.get())
	    {
	      // The parameter curve is a spline curve. This allows us
	      // to ensure correspondance in the split point
	      vector<Point> p1(2), p2(2);
	      tmp_par->point(p1, tmp_par->startparam(), 1);
	      tmp_par->replaceEndPoint(Point(par_u, par_v), true);

	      // Check consistence of tangent
	      tmp_par->point(p2, tmp_par->startparam(), 1);
	      if (p1[1].angle(p2[1]) > ang_tol)
		{
		  // Reset curve
		  tmp_par->replaceEndPoint(p1[0], true);
		  nmb_replace_param++;
		}
	    }
	}
    }
  else
    THROW("Missing parameter- and space curve");

  shared_ptr<ParamCurve> dummy;
  vector<shared_ptr<ParamCurve> > sub_cvs(2);
  for (int ki=0; ki<2; ++ki)
    {
      if (nmb_replace_param == 1)
	pcvs[ki] = dummy;  // Something wrong
      sub_cvs[ki] = 
	shared_ptr<ParamCurve>(new CurveOnSurface(surface_,pcvs[ki],
						  spacecvs[ki], 
						  prefer_parameter_,
						  ccm_, constdir_,
						  constval_, at_bd_,
						  same_orientation_));
    }
      

  return sub_cvs;
}

//===========================================================================
bool CurveOnSurface::ensureParCrvExistence(double epsgeo,
					   const RectDomain* domain_of_interest)
//===========================================================================
{
  const Point sf_epspar = SurfaceTools::getParEpsilon(*surface_, epsgeo);
  const double epspar = std::min(sf_epspar[0], sf_epspar[1]);

  if (!pcurve_)
    {
      // Check first for elementary curves and surfaces
      shared_ptr<ElementarySurface> elem_sf =
	dynamic_pointer_cast<ElementarySurface, ParamSurface>(surface_);
      shared_ptr<ElementaryCurve> elem_cv =
	dynamic_pointer_cast<ElementaryCurve, ParamCurve>(spacecurve_);
      if (elem_sf.get() && (!elem_cv.get()))
	{
	  shared_ptr<BoundedCurve> bd_cv =
	    dynamic_pointer_cast<BoundedCurve, ParamCurve>(spacecurve_);
	  if (bd_cv.get())
	    {
	      shared_ptr<ElementaryCurve> elem_cv2 =
		dynamic_pointer_cast<ElementaryCurve, ParamCurve>(bd_cv->underlyingCurve());
	      if (elem_cv2.get())
		elem_cv = shared_ptr<ElementaryCurve>(elem_cv2->clone());
	      elem_cv->setParamBounds(bd_cv->startparam(), bd_cv->endparam());
	    }
	}

      if (elem_sf.get() && elem_cv.get())
	{
	  // The function returns a curve only if the configuration is simple
	  pcurve_ = elem_sf->getElementaryParamCurve(elem_cv.get(), epspar);
	}
    }
	     
  // If the space curve and surface are not elementary geometry or the parameter curve
  // is not a simple elementary curve, use a more general approach
  if (!pcurve_)
    {
      Point startpt = faceParameter(startparam(), domain_of_interest);
      Point endpt = faceParameter(endparam(), domain_of_interest);

#ifndef NDEBUG
      {
	  Point spacecv_startpt = spacecurve_->point(spacecurve_->startparam());
	  Point spacecv_endpt = spacecurve_->point(spacecurve_->endparam());
	  Point lifted_startpt = surface_->point(startpt[0], startpt[1]);
	  Point lifted_endpt = surface_->point(endpt[0], endpt[1]);
	  double dist_start = spacecv_startpt.dist(lifted_startpt);
	  double dist_end = spacecv_endpt.dist(lifted_endpt);
	  if (dist_start > epspar || dist_end > epspar)
	  {
	      MESSAGE("Mismatch between space_cv and proj end pts, dist_start: "
		      << dist_start << ", dist_end: " << dist_end);

	      // We write to file the iso-curves in these parameters.
	      std::ofstream debug("tmp/iso_cvs.g2");
	      spacecurve_->writeStandardHeader(debug);
	      spacecurve_->write(debug);
	      vector<double> pts;
	      pts.insert(pts.end(), spacecv_startpt.begin(), spacecv_startpt.end());
	      pts.insert(pts.end(), lifted_startpt.begin(), lifted_startpt.end());
	      pts.insert(pts.end(), spacecv_endpt.begin(), spacecv_endpt.end());
	      pts.insert(pts.end(), lifted_endpt.begin(), lifted_endpt.end());
	      PointCloud3D pt_cl(pts.begin(), 4);
	      pt_cl.writeStandardHeader(debug);
	      pt_cl.write(debug);
	      vector<double> iso_par(4);
	      iso_par[0] = startpt[0];
	      iso_par[1] = startpt[1];
	      iso_par[2] = endpt[0];
	      iso_par[3] = endpt[1];
	      for (size_t ki = 0; ki < iso_par.size(); ++ki)
	      {
		  bool pardir_is_u = (ki%2 == 1); // Direction of moving parameter.
		  vector<shared_ptr<ParamCurve> > const_cvs = surface_->constParamCurves(iso_par[ki], pardir_is_u);
		  for (size_t kj = 0; kj < const_cvs.size(); ++kj)
		  {
		      const_cvs[kj]->writeStandardHeader(debug);
		      const_cvs[kj]->write(debug);
		  }
	      }
	      double debug_val = 0.0;
	  }
      }
#endif

      vector<Point> start;
      vector<Point> end;
      start.push_back(startpt);
      end.push_back(endpt);

      // First handle unbounded surfaces by replacing infinite bounds with "large" numbers
      const double large_number = 1.0e8;
      shared_ptr<ElementarySurface> es
          = dynamic_pointer_cast<ElementarySurface, ParamSurface>(surface_);
      bool isUnbounded = (es.get() != NULL && !es->isBounded());
      RectDomain dom = surface_->containingDomain();
      if  (isUnbounded) {
          double umin = dom.umin();
          double umax = dom.umax();
          double vmin = dom.vmin();
          double vmax = dom.vmax();
          double inf = numeric_limits<double>::infinity();
          if (umax == inf)
              umax = large_number;
          if (umin == -inf)
              umin = -large_number;
          if (vmax == inf)
              vmax = large_number;
          if (vmin == -inf)
              vmin = -large_number;
          Array<double, 2> ll(umin, vmin);
          Array<double, 2> ur(umax, vmax);
          dom = RectDomain(ll, ur);
          es->setParameterBounds(umin, vmin, umax, vmax);
      }
      
      // Check for closed surfaces. First check if the endpoint lies
      // at a boundary
      Point pos = spacecurve_->ParamCurve::point(startparam());
      Point close;
      double upar, vpar, dist;
      bool notfound = false;
      // We check that the surface is bounded before computing closest
      // points on the boundary
      isUnbounded = (es.get() != NULL && !es->isBounded());
      if  (isUnbounded) {
          MESSAGE("Shouldn't get here!");
          notfound = true;
      }
      else {
          try {
              surface_->closestBoundaryPoint(pos, upar, vpar, close, dist, epspar, 
                  &dom, startpt.begin());
          }
          catch (...)
          {
              notfound = true;
          }
      }
      if (notfound == false && pos.dist(close) < epspar)
	{
	  // The point lies at a boundary. Check the opposite boundary
	  if (startpt[0] - dom.umin() < epspar)
	    {
	      Point pos2 = surface_->point(dom.umax(), startpt[1]);
	      if (pos.dist(pos2) < epspar)
		start.push_back(Point(dom.umax(), startpt[1]));
	    }
	  else if (dom.umax() - startpt[0] < epspar)
	    {
	      Point pos2 = surface_->point(dom.umin(), startpt[1]);
	      if (pos.dist(pos2) < epspar)
		start.push_back(Point(dom.umin(), startpt[1]));
	    }
	  if (startpt[1] - dom.vmin() < epspar)
	    {
	      Point pos2 = surface_->point(startpt[0], dom.vmax());
	      if (pos.dist(pos2) < epspar)
		start.push_back(Point(startpt[0], dom.vmax()));
	    }
	  else if (dom.vmax() - startpt[1] < epspar)
	    {
	      Point pos2 = surface_->point(startpt[0], dom.vmin());
	      if (pos.dist(pos2) < epspar)
		start.push_back(Point(startpt[0], dom.vmin()));
	    }
	}
	  
      notfound = false;
      pos = spacecurve_->ParamCurve::point(endparam());
      if  (isUnbounded) {
          MESSAGE("Shouldn't get here!");
          notfound = true;
      }
      else {
          try {
              surface_->closestBoundaryPoint(pos, upar, vpar, close, dist, epspar, 
                  &dom, endpt.begin());
          }
          catch (...)
          {
              notfound = true;
          }
      }
      if (notfound == false && pos.dist(close) < epspar)
	{
	  // The point lies at a boundary. Check the opposite boundary
	  if (endpt[0] - dom.umin() < epspar)
	    {
	      Point pos2 = surface_->point(dom.umax(), endpt[1]);
	      if (pos.dist(pos2) < epspar)
		end.push_back(Point(dom.umax(), endpt[1]));
	    }
	  else if (dom.umax() - endpt[0] < epspar)
	    {
	      Point pos2 = surface_->point(dom.umin(), endpt[1]);
	      if (pos.dist(pos2) < epspar)
		end.push_back(Point(dom.umin(), endpt[1]));
	    }
	  if (endpt[1] - dom.vmin() < epspar)
	    {
	      Point pos2 = surface_->point(endpt[0], dom.vmax());
	      if (pos.dist(pos2) < epspar)
		end.push_back(Point(endpt[0], dom.vmax()));
	    }
	  else if (dom.vmax() - endpt[1] < epspar)
	    {
	      Point pos2 = surface_->point(endpt[0], dom.vmin());
	      if (pos.dist(pos2) < epspar)
		end.push_back(Point(endpt[0], dom.vmin()));
	    }
	}


      // For multiple start or end parameters we should look at the
      // corresponding tangent to see if we may reduce the number of
      // candidates. This will not handle cases with a tangent
      // along the seem.
      if (start.size() > 1)
      {
	  pickParamPoint(start, startparam(), epspar);
      }
      if (end.size() > 1)
      {
	  pickParamPoint(end, endparam(), epspar);
      }

#ifndef NDEBUG
      {
	  if (start.size() == 0 || end.size() == 0)
	  {
	      MESSAGE("Oops, missing end point(s).");
	  }
      }
#endif

     for (size_t ki=0; ki<start.size(); ++ki)
	{
	  for (size_t kj=0; kj<end.size(); ++kj)
	    {
	      shared_ptr<Point> pt1(new Point(start[ki].begin(), 
					      start[ki].end()));
	      shared_ptr<Point> pt2(new Point(end[kj].begin(), end[kj].end()));
	      shared_ptr<SplineCurve> pcv;
	      try {
		pcv = shared_ptr<SplineCurve>(CurveCreators::projectSpaceCurve(spacecurve_,
									       surface_,
									       pt1, pt2, epsgeo));//epspar));
	      }
	      catch(...)
		{
		}
	      if (pcv.get())
		{
#ifndef NDEBUG
		    {
			double partol =1e-12;
			if ((fabs(pcv->startparam() - spacecurve_->startparam()) > partol) ||
			    (fabs(pcv->endparam() - spacecurve_->endparam()) > partol))
			{
			    MESSAGE("End parameters not tol-equal (1e-12)!");
			}
		    }
#endif
		  pcurve_ = pcv;
		  break;
		}
	    }
	  if (pcurve_)
	    break;
	}
    }
  return pcurve_.get() != 0;
}

//===========================================================================
bool CurveOnSurface::makeParameterCurve(double tol, const Point& par1, 
					 const Point& par2)
//===========================================================================
{
  // No test performed on par1 & par2 prior to projecting (i.e. trying to).
  bool remade = false;
  Point parval1 = faceParameter(startparam());
  Point parval2 = faceParameter(endparam());
  shared_ptr<Point> pt1(new Point(par1.begin(), par1.end()));
  shared_ptr<Point> pt2(new Point(par2.begin(), par2.end()));
  shared_ptr<SplineCurve> pcv;
  try {
    pcv = shared_ptr<SplineCurve>(CurveCreators::projectSpaceCurve(spacecurve_,
								   surface_,
								   pt1, pt2, tol));
  }
  catch(...)
    {
    }
  if (pcv.get())
    {
#ifndef NDEBUG
	{
	    double partol =1e-12;
	    if ((fabs(pcv->startparam() - spacecurve_->startparam()) > partol) ||
		(fabs(pcv->endparam() - spacecurve_->endparam()) > partol))
	    {
		MESSAGE("End parameters not tol-equal (1e-12)!");
	    }
	}
#endif
      pcurve_ = pcv;
      if (constdir_ > 0)
	{
	  if (constdir_ == 1)
	    constval_ += (par1[0] - parval1[0]);
	  else
	    constval_ += (par1[1] - parval1[1]);

	  same_orientation_ = (par2[constdir_-1] > par1[constdir_-1]);
	  
	  RectDomain dom = surface_->containingDomain();
	  double ptol = 1.0e-6;
	  if (constdir_ == 1 && fabs(constval_ - dom.umin()) < ptol)
	    at_bd_ = 0;
	  else if (constdir_ == 1 && fabs(dom.umax() - constval_) < ptol)
	    at_bd_ = 1;
	  else if (constdir_ == 2 && fabs(constval_ - dom.vmin()) < ptol)
	    at_bd_ = 2;
	  else if (constdir_ == 2 && fabs(dom.vmax() - constval_) < ptol)
	    at_bd_ = 3;
	}

      remade = true;
    }
  return remade;
}

//===========================================================================
bool CurveOnSurface::translateParameterCurve(const Point& dir)
//===========================================================================
{
  shared_ptr<SplineCurve> pcrv =
    dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve_);
  shared_ptr<ElementaryCurve> pcrv2 = 
    dynamic_pointer_cast<ElementaryCurve, ParamCurve>(pcurve_);

  if (pcrv.get())
    pcrv->translateCurve(dir);
  else if (pcrv2.get())
    pcrv2->translateCurve(dir);
  else
    return false;

  // Update constant parameter values
  double tp = 0.5*(pcurve_->startparam() + pcurve_->endparam());
  Point par = pcurve_->point(tp);
  RectDomain dom = surface_->containingDomain();
  if (constdir_ == 1)
    constval_ += dir[0];
  if (constdir_ == 2)
    constval_ += dir[1];
  if (at_bd_ >= 0)
    {
      // Check if the boundary information is still valid
      double ptol = 1.0e-4;
      if (at_bd_ == 0 && fabs(par[0] - dom.umin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 1 && fabs(par[0] - dom.umax()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 2 && fabs(par[1] - dom.vmin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 3 && fabs(par[1] - dom.vmax()) > ptol)
	at_bd_ = -1;
    }

  return true;
}

//===========================================================================
bool CurveOnSurface::translateSwapParameterCurve(const Point& dir, double sgn,
						 int pdir)
//===========================================================================
{
  shared_ptr<SplineCurve> pcrv =
    dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve_);
  if (!pcrv.get())
    return false;
  
  pcrv->translateSwapCurve(dir, sgn, pdir);

  // Update constant parameter values
  double tp = 0.5*(pcrv->startparam() + pcrv->endparam());
  Point par = pcurve_->point(tp);
  RectDomain dom = surface_->containingDomain();
  if (constdir_ == 1 && (pdir == 1 || pdir == 3))
    constval_ = (sgn*constval_+dir[0]);
  if (constdir_ == 2 && (pdir == 2 || pdir == 3))
    constval_ = (sgn*constval_+dir[1]);
  if (at_bd_ >= 0)
    {
      // Check if the boundary information is still valid
      double ptol = 1.0e-4;
      if (at_bd_ == 0 && fabs(par[0] - dom.umin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 1 && fabs(par[0] - dom.umax()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 2 && fabs(par[1] - dom.vmin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 3 && fabs(par[1] - dom.vmax()) > ptol)
	at_bd_ = -1;
    }

  return true;
}

//===========================================================================
bool CurveOnSurface::setDomainParCrv(double umin, double umax, 
				     double vmin, double vmax,
				     double uminprev, double umaxprev,
				     double vminprev, double vmaxprev)
//===========================================================================
{
  shared_ptr<SplineCurve> pcrv =
    dynamic_pointer_cast<SplineCurve, ParamCurve>(pcurve_);
  if (!pcrv.get())
    return false;  // Not a spline curve, cannot change coefficients
  
  double old_diff_u = umaxprev - uminprev;
  double old_diff_v = vmaxprev - vminprev;
  double new_diff_u = umax - umin;
  double new_diff_v = vmax - vmin;
  // double umin_diff = umin - uminprev;
  // double vmin_diff = vmin - vminprev;

  if (pcrv->rational())
    {
      vector<double>::iterator iter = pcrv->rcoefs_begin();
      while (iter != pcrv->rcoefs_end()) {
	double w1 = iter[2];
	double u = iter[0]/w1;
	u = (u - uminprev)*new_diff_u/old_diff_u + umin;
	iter[0] = u*w1;
	double v = iter[1]/w1;
	v = (v - vminprev)*new_diff_v/old_diff_v + vmin;
	iter[1] = v*w1;
	iter+=3;
      }
      pcrv->updateCoefsFromRcoefs();
    }
  else
    {
      vector<double>::iterator iter = pcrv->coefs_begin();
      while (iter != pcrv->coefs_end()) {
	iter[0] = (iter[0]-uminprev)*new_diff_u/old_diff_u + umin;
	iter[1] = (iter[1]-vminprev)*new_diff_v/old_diff_v + vmin;
	iter+=2;
      }
    }

   if (constdir_ == 1)
    {
      constval_ = (constval_-uminprev)*new_diff_u/old_diff_u + umin;
    }
  else if (constdir_ == 2)
    {
      constval_ = (constval_-vminprev)*new_diff_v/old_diff_v + vmin;
    }
  if (at_bd_ >= 0)
    {
      // Check if the boundary information is still valid
      double ptol = 1.0e-4;
      double tp = 0.5*(pcrv->startparam() + pcrv->endparam());
      Point par = pcurve_->point(tp);
      RectDomain dom = surface_->containingDomain();
      if (at_bd_ == 0 && fabs(par[0] - dom.umin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 1 && fabs(par[0] - dom.umax()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 2 && fabs(par[1] - dom.vmin()) > ptol)
	at_bd_ = -1;
      else if (at_bd_ == 3 && fabs(par[1] - dom.vmax()) > ptol)
	at_bd_ = -1;
    }
  return true;
}

//===========================================================================
bool CurveOnSurface:: ensureSpaceCrvExistence(double tol)
//===========================================================================
{
  if (!spacecurve_)
    {
      shared_ptr<SplineCurve> spacecv;
      try {
	spacecv = 
	  shared_ptr<SplineCurve>(CurveCreators::liftParameterCurve(pcurve_,
								    surface_,
								    tol));
      }
      catch(...)
	{
	}
      if (spacecv.get())
	spacecurve_ = spacecv;
    }
  return spacecurve_.get() != 0;
}

///===========================================================================
bool CurveOnSurface::updateIsoCurves(int constdir, double constpar, 
				     int boundary)
//===========================================================================
{
  // Create new space curve
  vector<shared_ptr<ParamCurve> > isocurves = 
    surface_->constParamCurves(constpar, (constdir == 2));

  if (isocurves.size() != 1)
    return false;   // More than one curve. Don't know how to update
  
  // Check orientation
  Point s1 = spacecurve_->point(spacecurve_->startparam());
  Point e1 = spacecurve_->point(spacecurve_->endparam());
  Point s2 = isocurves[0]->point(isocurves[0]->startparam());
  Point e2 = isocurves[0]->point(isocurves[0]->endparam());
  if ((e1-s1)*(e2-s2) < 0.0)
    same_orientation_ = false;

  prefer_parameter_ = false;
  pcurve_ = shared_ptr<ParamCurve>();
  double t1 = startparam();
  double t2 = endparam();
  Point pt1 = faceParameter(t1);
  Point pt2 = faceParameter(t2);

  spacecurve_ = isocurves[0];
  int dir = 2 - constdir;
  double tmin = std::min(pt1[dir], pt2[dir]);
  double tmax = std::max(pt1[dir], pt2[dir]);
  spacecurve_ = shared_ptr<ParamCurve>(spacecurve_->subCurve(tmin, tmax));
  if (!same_orientation_)
    spacecurve_->reverseParameterDirection();

  spacecurve_->setParameterInterval(t1, t2);

  pcurve_ = shared_ptr<ParamCurve>(new SplineCurve(pt1, t1, pt2, t2));
  
  constdir_ = constdir;
  constval_ = constpar;
  at_bd_ = boundary;
  same_orientation_ = true;
  ccm_ = 3;

  
  return true;  // Update performed
}


///===========================================================================
bool CurveOnSurface::updateIsoCurves()
//===========================================================================
{
  // Only performed if the curve is set as a constant parameter curve
  if (constdir_ < 1 || constdir_ > 2)
    return false;   

  // Create new space curve
  vector<shared_ptr<ParamCurve> > isocurves = 
    surface_->constParamCurves(constval_, (constdir_ == 2));

  if (isocurves.size() != 1)
    return false;   // More than one curve. Don't know how to update
  
  spacecurve_ = isocurves[0];
  if (!same_orientation_)
    spacecurve_->reverseParameterDirection();
  
  return true;  // Update performed
}


//===========================================================================
bool CurveOnSurface::updateCurves(double epsge)
//===========================================================================
{
  if (constdir_ == 1 || constdir_ == 2)
    return updateIsoCurves();

  // Define evaluator based curve
  shared_ptr<TrimCurve> trim_crv = 
    shared_ptr<TrimCurve>(new TrimCurve(this));

  // Define approximator
  vector<int> dims(2);
  dims[0] = dimension();
  dims[1] = 2;  // Curve in the parameter domain
  shared_ptr<HermiteAppS> approximator =
    shared_ptr<HermiteAppS>(new HermiteAppS(trim_crv.get(), epsge, epsge, dims));
  try {
    // Approximate
    approximator->refineApproximation();
  }
  catch (...)
    {
      return false;
    }

  // Fetch curves
  vector<shared_ptr<SplineCurve> > crvs = approximator->getCurves();

  // Update curves
  spacecurve_ = crvs[0];
  pcurve_ = crvs[1];

  return true;
  
}

//===========================================================================
bool CurveOnSurface::updateCurves(Point vx1, Point vx2, double epsge)
//===========================================================================
{
  if (constdir_ == 1 || constdir_ == 2)
    return updateIsoCurves();

  // Define evaluator based curve
//   shared_ptr<TrimCurve> trim_crv = 
//     shared_ptr<TrimCurve>(new TrimCurve(vx1, vx2, this));
  shared_ptr<TrimCurve> trim_crv = 
    shared_ptr<TrimCurve>(new TrimCurve(this));

  // Define approximator
  vector<int> dims(2);
  dims[0] = dimension();
  dims[1] = 2;  // Curve in the parameter domain
  shared_ptr<HermiteAppS> approximator =
    shared_ptr<HermiteAppS>(new HermiteAppS(trim_crv.get(), epsge, epsge, dims));
  try {
    // Approximate
    approximator->refineApproximation();
  }
  catch (...)
    {
      return false;
    }

  // Fetch curves
  vector<shared_ptr<SplineCurve> > crvs = approximator->getCurves();

  // Update curves
  spacecurve_ = crvs[0];
  pcurve_ = crvs[1];

  return true;
  
}

//===========================================================================
Point CurveOnSurface::faceParameter(double crv_par, 
				    const RectDomain* domain_of_interest) const
//===========================================================================
{
  Point param(2);
  bool same = same_orientation_;
  if (constdir_ > 0)
    {
      if (pcurve_)
	{
	  param = pcurve_->ParamCurve::point(crv_par);
	  crv_par = param[2-constdir_];
	  same = true;
	}

      if (constdir_ == 1)
	{
	  param[0] = constval_;
	  param[1] = (same) ? crv_par : endparam() - (crv_par - startparam());
// 	  param[1] = crv_par; // VSK, 1004. More consistent, can it create problems?
	}
      else if (constdir_ == 2)
	{
	  param[0] = (same) ? crv_par : endparam() - (crv_par - startparam());
	  //  param[0] = crv_par; // VSK, 1004. More consistent, can it create problems?
	  param[1] = constval_;
	}
    }
  else if (pcurve_)
    {
      param = pcurve_->ParamCurve::point(crv_par);
      if (!prefer_parameter_)
	{
	  // Iterate to get to a better position
	  Point pos = spacecurve_->ParamCurve::point(crv_par);
	  
	  double clo_u, clo_v, clo_dist;
	  double eps = 1.0e-6;
	  Point clo_pt;
	  Point seed = param;
	  surface_->closestPoint(pos, clo_u, clo_v, clo_pt, clo_dist, eps,
				 domain_of_interest, seed.begin());
	  param = Point(clo_u, clo_v);
	}
    }	  
  else
    {
      // No parameter curve exist. Perform closest point computation
      Point pos = spacecurve_->ParamCurve::point(crv_par);

      double clo_u, clo_v, clo_dist;
      double eps = 1.0e-6;
      Point clo_pt;
      surface_->closestPoint(pos, clo_u, clo_v, clo_pt, clo_dist, eps,
			     domain_of_interest);
      param = Point(clo_u, clo_v);
    }
  return param;
}

//===========================================================================
int CurveOnSurface::whichBoundary(double tol, bool& same_orientation) const
//===========================================================================
{
  if (constdir_ == 1 || constdir_ == 2)
    {
      same_orientation = same_orientation_;
      return at_bd_;
    }
  else 
    {
      // Check if the curve is iso parametric
      if (!pcurve_.get())
	return -1;   // Can't check

      BoundingBox box2D = pcurve_->boundingBox();
      Point low = box2D.low();
      Point high = box2D.high();
      if (high[0] - low[0] > tol && high[1]-low[1] > tol)
	{
	  // 2D curve different from a line
	  return -1;
	}

      RectDomain domain = surface_->containingDomain();
      Array<double,2> p1(low[0], low[1]);
      Array<double,2> p2(high[0], high[1]);
      int bd = domain.whichBoundary(p1, p2, tol);

      // Check orientation
      Point pt1 = pcurve_->point(pcurve_->startparam());
      Point pt2 = pcurve_->point(pcurve_->endparam());
      int idx = (bd <= 1) ? 1 : 0;
      same_orientation = (pt2[idx] > pt1[idx]);
      
      return bd;
    }
}

//===========================================================================
bool CurveOnSurface::isConstantCurve(double tol, int& pardir, double& parval) const
//===========================================================================
{
  if (isConstantCurve())
    {
      pardir = constdir_;
      parval = constval_;
      return true;
    }
  else if (!pcurve_.get())
    return false;
  else
    {
      // Make box around curve
      BoundingBox box2d = pcurve_->boundingBox();
      Point high = box2d.high();
      Point low = box2d.low();
      if (high[0] - low[0] > tol && high[1]-low[1] > tol)
	return false;
      else if (high[0] - low[0] <= tol)
	{
	  pardir = 1;
	  parval = 0.5*(high[0] + low[0]);
	  return true;
	}
      else
	{
	  pardir = 2;
	  parval = 0.5*(high[1] + low[1]);
	  return true;
	}
    }
}

//===========================================================================
bool CurveOnSurface::sameParameterDomain() const
//===========================================================================
{
    double tol = 1.0e-12;

    if (pcurve_.get() && spacecurve_.get())
    {
	// Both versions exist
	
	// Check start parameter
	if (fabs(pcurve_->startparam() - spacecurve_->startparam()) > tol)
	    return false;

	// Check end parameter
	if (fabs(pcurve_->endparam() - spacecurve_->endparam()) > tol)
	    return false;
    }
	
    return true;
}

//===========================================================================
bool CurveOnSurface::sameOrientation() const
//===========================================================================
{
    if (!pcurve_.get() || !spacecurve_.get())
	return true;  // Only one curve, no consistency problems

    // Evaluate endpoints in both versions of the curve
    Point geomstart, geomend, parstart, parend;
    Point par1, par2;

    // Space curve
    spacecurve_->point(geomstart, spacecurve_->startparam());
    spacecurve_->point(geomend, spacecurve_->endparam());

    // If space end pts are identical then we need to compare pts
    // along the curve to check if orientation coincides. If cv is deg
    // then the orientation is ok by definition.
    double space_dist = geomstart.dist(geomend);
    double epsgeo = 1e-06;
    if (space_dist < epsgeo)
	return true;

    // Parameter curve
    pcurve_->point(par1, pcurve_->startparam());
    pcurve_->point(par2, pcurve_->endparam());
    surface_->point(parstart, par1[0], par1[1]);
    surface_->point(parend, par2[0], par2[1]);

    if (geomstart.dist(parstart) + geomend.dist(parend) <
	geomstart.dist(parend) + geomend.dist(parstart))
	return true;
    else
	return false;
}

//===========================================================================
bool CurveOnSurface::sameTrace(double tol, int nmb_samples) const
//===========================================================================
{
    if (!pcurve_.get() || !spacecurve_.get())
	return true;  // Only one curve, no consistency problems

    double tp1 =  pcurve_->startparam();
    double tp2 =  pcurve_->endparam();
    double tg1 = spacecurve_->startparam();
    double tg2 = spacecurve_->endparam();
    if (!sameOrientation())
      std::swap(tg1, tg2);

//     int nmb_sample = 5;
    double tdelp = (tp2-tp1)/(double)(nmb_samples-1);
    double tdelg = (tg2-tg1)/(double)(nmb_samples-1);

    int ki;
    double tp, tg;
    for (ki=0, tp=tp1, tg=tg1; ki<nmb_samples; ++ki, tp+=tdelp, tg+=tdelg)
      {
	Point pntp, pntg;
	Point par;
	double tg_close, dist;
	pcurve_->point(par, tp);
	surface_->point(pntp, par[0], par[1]);
	spacecurve_->closestPoint(pntp, spacecurve_->startparam(),
				  spacecurve_->endparam(), tg_close,
				  pntg, dist, &tg);
	if (dist > tol)
	  return false;
      }

    return true;
}

//===========================================================================
bool CurveOnSurface::sameCurve(double tol, int nmb_samples) const
//===========================================================================
{
    if (!pcurve_.get() || !spacecurve_.get())
	return true;  // Only one curve, no consistency problems


    double tp1 =  pcurve_->startparam();
    double tp2 =  pcurve_->endparam();
    double tg1 = spacecurve_->startparam();
    double tg2 = spacecurve_->endparam();
    if (!sameOrientation())
      std::swap(tg1, tg2);

//     int nmb_sample = 5;
    double tdelp = (tp2-tp1)/(double)(nmb_samples-1);
    double tdelg = (tg2-tg1)/(double)(nmb_samples-1);

    int ki;
    double tp, tg;
    for (ki=0, tp=tp1, tg=tg1; ki<nmb_samples; ++ki, tp+=tdelp, tg+=tdelg)
      {
	Point pntp, pntg;
	Point par;
	pcurve_->point(par, tp);
	surface_->point(pntp, par[0], par[1]);
	spacecurve_->point(pntg, tg);
	double dist = pntp.dist(pntg);
	if (dist > tol)
	  return false;
      }

    return true;
}


//===========================================================================
void CurveOnSurface::makeCurvesConsistent(bool prefer_parameter)
//===========================================================================
{
    if (!pcurve_.get() || !spacecurve_.get())
      return;// true;  // Only one curve, no consistency problems

    if (!sameOrientation())
      {
	if (prefer_parameter)
	  spacecurve_->reverseParameterDirection();
	else
	  pcurve_->reverseParameterDirection();
	fix_performed_ += 1;
      }

    if (!sameParameterDomain())
      {
	if (prefer_parameter)
	  spacecurve_->setParameterInterval(pcurve_->startparam(),
					    pcurve_->endparam());
	else
	   pcurve_->setParameterInterval(spacecurve_->startparam(),
					 spacecurve_->endparam());
	fix_performed_ += 2;
      }
    prefer_parameter_ = prefer_parameter;

//     return sameCurve(tol);  // Check if fix succeeded
}


//===========================================================================
int CurveOnSurface::curveCreationMethod() const
//===========================================================================
{
    return ccm_;
}


//===========================================================================
double CurveOnSurface::maxTraceDiff(int nmb_sample) const
//===========================================================================
{
    if (!pcurve_.get() || !spacecurve_.get())
	return true;  // Only one curve, no consistency problems

    double max_dist = -1.0;
    double tp1 =  pcurve_->startparam();
    double tp2 =  pcurve_->endparam();
    double tg1 = spacecurve_->startparam();
    double tg2 = spacecurve_->endparam();
    if (!sameOrientation())
      std::swap(tg1, tg2);

//     int nmb_sample = 5;
    double tdelp = (tp2-tp1)/(double)(nmb_sample-1);
    double tdelg = (tg2-tg1)/(double)(nmb_sample-1);

    int ki;
    double tp, tg;
    for (ki=0, tp=tp1, tg=tg1; ki<nmb_sample; ++ki, tp+=tdelp, tg+=tdelg)
      {
	Point pntp, pntg;
	Point par;
	double tg_close, dist;
	pcurve_->point(par, tp);
	surface_->point(pntp, par[0], par[1]);
	spacecurve_->closestPoint(pntp, spacecurve_->startparam(),
				  spacecurve_->endparam(), tg_close,
				  pntg, dist, &tg);
	if (dist > max_dist)
	  max_dist = dist;
      }

    return max_dist;
}

//===========================================================================
bool CurveOnSurface::isAxisRotational(Point& centre, Point& axis, Point& vec,
				   double& angle)
//===========================================================================
{
  if (spacecurve_.get())
    return spacecurve_->isAxisRotational(centre, axis, vec, angle);
  else
    return false;
}

//===========================================================================
bool CurveOnSurface::isLinear(Point& dir, double tol)
//===========================================================================
{
  if (spacecurve_.get())
    return spacecurve_->isLinear(dir, tol);
  else
    return false;
}

//===========================================================================
bool CurveOnSurface::isInPlane(const Point& loc, const Point& axis,
			       double eps, Point& normal) const
//===========================================================================
{
  if (spacecurve_.get())
    return spacecurve_->isInPlane(loc, axis, eps, normal);
  else
    return false;
}

//===========================================================================
bool CurveOnSurface::isInPlane(const Point& norm,
			       double eps, Point& pos) const
//===========================================================================
{
  if (spacecurve_.get())
    return spacecurve_->isInPlane(norm, eps, pos);
  else
    return false;
}


//===========================================================================
void CurveOnSurface::pickParamPoint(vector<Point>& par_candidates,
				    double tpar, double epsgeo)
//===========================================================================
{
    // @@sbr Possibly make this function public if it works out well.

    // We require the space curve to be defined.
    if (!spacecurve_)
    {
	return;
    }

    Point sf_pt = ParamCurve::point(tpar);

    // We also expect the par_candidates to refer to the same sf pt.
    for (size_t ki = 0; ki < par_candidates.size(); ++ki)
    {
	Point pt = surface_->point(par_candidates[ki][0], par_candidates[ki][1]);
	Point sf_pt2 = surface_->point(par_candidates[1][0], par_candidates[1][1]);
	if (sf_pt.dist(pt) > epsgeo)
	{
	    par_candidates.erase(par_candidates.begin() + ki);
	    --ki;
	}
    }

    if (par_candidates.size() == 0)
    {
	return;
    }

    bool closed_dir_u, closed_dir_v;
    Go::SurfaceTools::checkSurfaceClosed(*surface_, closed_dir_u, closed_dir_v, epsgeo);
    bool is_closed = (closed_dir_u || closed_dir_v);

    // Assuming that our domain is the full domain. For our purposes I
    // suppose it is.
    RectDomain rect_dom = surface_->containingDomain();
    // We consider an angle of more than tang_tol to be not tangential.
    const double tang_tol = 1e-02;
    Point u_dir(1.0, 0.0);
    Point v_dir(0.0, 1.0);
    for (size_t ki = 0; ki < par_candidates.size(); ++ki)
    {
	double upar = par_candidates[ki][0];
	double vpar = par_candidates[ki][1];

	const double knot_diff_tol = 1e-08;
	// If we are at the end of the curve, our test differs slightly.
	bool at_cv_end = (fabs(tpar - endparam()) < knot_diff_tol);
	bool at_cv_start = (fabs(tpar - startparam()) < knot_diff_tol);

	bool at_u_start = (fabs(upar - rect_dom.umin()) < knot_diff_tol);
	bool at_u_end = (fabs(upar - rect_dom.umax()) < knot_diff_tol);
	bool at_v_start = (fabs(vpar - rect_dom.vmin()) < knot_diff_tol);
	bool at_v_end = (fabs(vpar - rect_dom.vmax()) < knot_diff_tol);

	// By at_u_bd we mean that the seem corresponds to a u-parameter.
	bool at_u_bd = (at_u_start || at_u_end);
	bool at_v_bd = (at_v_start || at_v_end);

	if (!at_u_bd && !at_v_bd){
	    continue;
	}

	// We then project the space curve tangent onto the surface.
	Point par_tangent = projectSpaceCurveTangent(par_candidates[ki], tpar);

	vector<Point> sf_pt = surface_->point(par_candidates[ki][0], par_candidates[ki][1], 1);
	vector<Point> cv_pt = spacecurve_->point(tpar, 1);
	double ang_u_space = cv_pt[1].angle(sf_pt[1]);
	double ang_v_space = cv_pt[1].angle(sf_pt[2]);

	// The range of angle2() is [0, 2*M_PI).
	double ang_u = u_dir.angle2(par_tangent);
	// But we want [-M_PI, M_PI).
	if (ang_u >= M_PI)
	{
	    ang_u -= 2*M_PI;
	}
	// To fix problems with uneven scaling of domain directions we check angle for space tangents.
	if (ang_u_space < tang_tol)
	{
	    ang_u = 0.0;
	}
	double ang_v = v_dir.angle2(par_tangent);
	if (ang_v >= M_PI)
	{
	    ang_v -= 2*M_PI;
	}
	if (ang_v_space < tang_tol)
	{
	    ang_v = 0.0;
	}

	// // We do not yet handle corner points of toruses (and its like).
	// if (at_u_bd && at_v_bd){
	//     MESSAGE("Not yet handled!");
	//     continue;
	// }
	// else if (at_u_bd)
	if (at_u_bd && closed_dir_u)
	{
	    if ((fabs(ang_v) < tang_tol) || ((fabs(ang_v + M_PI) < tang_tol)) || ((fabs(ang_v - M_PI) < tang_tol)))
	    {
		; // Do nothing.
	    }
	    else
	    {
		if (at_cv_end) // The end of the space curve.
		{
		    if (((ang_v < 0.0) && at_u_start) ||
			((ang_v > 0.0) && at_u_end))
		    {
			par_candidates.erase(par_candidates.begin() + ki);
			--ki;
			continue;
		    }
		}
		else if (at_cv_start)
		{
		    if (((ang_v < 0.0) && at_u_end) ||
			((ang_v > 0.0) && at_u_start))
		    {
			par_candidates.erase(par_candidates.begin() + ki);
			--ki;
			continue;
		    }
		}
		else
		{
		    MESSAGE("This routine does not handle curves crossing the seem!");
		    continue;
		}
	    }
	}
	// else // at_v_bd
	if (at_v_bd && closed_dir_v)
	{
	    if ((fabs(ang_u) < tang_tol) || ((fabs(ang_u + M_PI) < tang_tol)) || ((fabs(ang_u - M_PI) < tang_tol)))
	    {
		; // Do nothing.
	    }
	    else
	    {
		if (at_cv_end) // The end of the space curve.
		{
		    if (((ang_u > 0.0) && at_v_start) ||
			((ang_u < 0.0) && at_v_end))
		    {
			par_candidates.erase(par_candidates.begin() + ki);
			--ki;
			continue;
		    }
		}
		else if (at_cv_start)
		{
		    if (((ang_u > 0.0) && at_v_end) ||
			((ang_u < 0.0) && at_v_start))
		    {
			par_candidates.erase(par_candidates.begin() + ki);
			--ki;
			continue;
		    }
		}
		else
		{
		    MESSAGE("This routine does not handle curves crossing the seem!");
		    continue;
		}
	    }
	}
    }
}


//===========================================================================
Point CurveOnSurface::projectSpaceCurveTangent(const Point& par_pt, double tpar)
//===========================================================================
{
    vector<Point> space_cv_pt = spacecurve_->point(tpar, 1);
    const int dim = spacecurve_->dimension();
    // We then project the tangent into the parameter domain.
    vector<Point> sf_pt = surface_->ParamSurface::point(par_pt[0], par_pt[1], 1);
    // We describe the curve tangent as a linear combination of the partial derivs.
    double coef1, coef2;
    blendcoef(&sf_pt[1][0], &sf_pt[2][0],
	      &space_cv_pt[1][0], dim, 1, &coef1, &coef2);

    Point pt_dir(2);
    pt_dir[0] = coef1;
    pt_dir[1] = coef2;

    return pt_dir;
}

