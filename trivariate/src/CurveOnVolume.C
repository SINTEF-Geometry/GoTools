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

#include "GoTools/trivariate/CurveOnVolume.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/Factory.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;
using std::max;
using std::min;
using std::streamsize;
using std::endl;

//===========================================================================
CurveOnVolume::CurveOnVolume()
  : prefer_parameter_(true)
//===========================================================================
{
}

//===========================================================================
CurveOnVolume::CurveOnVolume(shared_ptr<ParamVolume> vol,
			     shared_ptr<ParamCurve> curve,
			     bool preferparameter)
  : volume_(vol), prefer_parameter_(preferparameter)
//===========================================================================
{
  if (prefer_parameter_)
    pcurve_ = curve;
  else
    spacecurve_ = curve;
}

//===========================================================================
CurveOnVolume::CurveOnVolume(shared_ptr<ParamVolume> vol,
			     shared_ptr<ParamCurve> pcurve,
			     shared_ptr<ParamCurve> spacecurve,
			     bool preferparameter)
  : volume_(vol), pcurve_(pcurve), spacecurve_(spacecurve),
    prefer_parameter_(preferparameter)
//===========================================================================
{
}

//===========================================================================
CurveOnVolume::CurveOnVolume(const CurveOnVolume& volume_curve)
//===========================================================================
{
  // Clones the curves, not the volume
  volume_ = volume_curve.volume_;
  prefer_parameter_ = volume_curve.prefer_parameter_;
  if (volume_curve.pcurve_.get() != 0)
    {
      ParamCurve *pcv = volume_curve.pcurve_->clone();
      pcurve_  = shared_ptr<ParamCurve>(pcv);
    }
  else
    pcurve_ = shared_ptr<ParamCurve>();
  if (volume_curve.spacecurve_.get() != 0)
    {
      ParamCurve *scv = volume_curve.spacecurve_->clone();
      spacecurve_ = shared_ptr<ParamCurve>(scv);
    }
  else
    spacecurve_ = shared_ptr<ParamCurve>();

}

//===========================================================================
CurveOnVolume& CurveOnVolume::operator=(const CurveOnVolume& other)
//===========================================================================
{
  if (&other != this) {
    // Clones the curves, not the volume
    volume_ = other.volume_;
    prefer_parameter_ = other.prefer_parameter_;
    if (other.pcurve_.get() != 0)
      {
	ParamCurve *pcv = other.pcurve_->clone();
	pcurve_  = shared_ptr<ParamCurve>(pcv);
      }
    else
      pcurve_ = shared_ptr<ParamCurve>();
    if (other.spacecurve_.get() != 0)
	{
	  ParamCurve *scv = other.spacecurve_->clone();
	  spacecurve_ = shared_ptr<ParamCurve>(scv);
	}
    else
      spacecurve_ = shared_ptr<ParamCurve>();

  }
  return *this;
}

//===========================================================================
CurveOnVolume::~CurveOnVolume()
//===========================================================================
{}

//===========================================================================
void CurveOnVolume::read(std::istream& is)
//===========================================================================
{
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
}

//===========================================================================
void CurveOnVolume::write(std::ostream& os) const
//===========================================================================
{
  streamsize prev = os.precision(15);

  // We do not write the volume. It must be reconnected while reading
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
BoundingBox CurveOnVolume::boundingBox() const
//===========================================================================
{
    if (spacecurve_.get() != 0)
      return spacecurve_->boundingBox();
    else
      return volume_->boundingBox();
}

//===========================================================================
DirectionCone CurveOnVolume::directionCone() const
//===========================================================================
{
    if (spacecurve_.get() != 0)
      return spacecurve_->directionCone();
    else
      {
	// @@@ VSK 1611. Returns an empty cone for the time being
	DirectionCone cone;
	return cone;
      }
}

//===========================================================================
int CurveOnVolume::dimension() const
//===========================================================================
{
  if (volume_.get())
    return volume_->dimension();
  else if (spacecurve_.get())
    return spacecurve_->dimension();
  else
    return 3;  // Should not occur
}

//===========================================================================
ClassType CurveOnVolume::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
void CurveOnVolume::point(Point& pt, double tpar) const
//===========================================================================
{
    if (prefer_parameter_) {
	Point param_pt(3);
	pcurve_->point(param_pt, tpar);
	volume_->point(pt, param_pt[0], param_pt[1], param_pt[2]);
    } else {
	spacecurve_->point(pt, tpar);
    }
}

//===========================================================================
void CurveOnVolume::point(std::vector<Point>& pts,
			  double tpar,
			  int derivs,
			  bool from_right) const
//===========================================================================
{
  if (prefer_parameter_)
  {
    // @@@
    ALWAYS_ERROR_IF(derivs > 2, "Only two derivatives are supported.");
    std::vector<Point> param_pts(derivs+1, Point(3));
    std::vector<Point> vol_pts((derivs+1)*(derivs+2)*(derivs+3)/6, 
			       Point(volume_->dimension()));
    pcurve_->point(param_pts, tpar, derivs, from_right);
    volume_->point(vol_pts, param_pts[0][0], param_pts[0][1], param_pts[0][2], 
		   derivs);
    pts.resize(derivs+1);
    pts[0] = vol_pts[0];
    if (derivs >= 1) 
      {
	pts[1] = param_pts[1][0]*vol_pts[1] + param_pts[1][1]*vol_pts[2] +
	  param_pts[1][2]*vol_pts[3];
      }
    if (derivs == 2)
      {
	pts[2] = param_pts[2][0]*vol_pts[1] + param_pts[2][1]*vol_pts[2] +
	  param_pts[2][2]*vol_pts[3] + 
	  param_pts[1][0]*param_pts[1][0]*vol_pts[4] +
	  param_pts[1][0]*param_pts[1][1]*vol_pts[5] +
	  param_pts[1][0]*param_pts[1][2]*vol_pts[6] +
	  param_pts[1][1]*param_pts[1][1]*vol_pts[7] +
	  param_pts[1][1]*param_pts[1][2]*vol_pts[8] +
	  param_pts[1][2]*param_pts[1][2]*vol_pts[8];
      }
  }
  else
    spacecurve_->point(pts, tpar, derivs, from_right);
}


//===========================================================================
double CurveOnVolume::startparam() const
//===========================================================================
{
    if (prefer_parameter_ && (pcurve_.get() != NULL))
    return pcurve_->startparam();
  else
    return spacecurve_->startparam();
}


//===========================================================================
double CurveOnVolume::endparam() const
//===========================================================================
{
  if (prefer_parameter_ && (pcurve_.get() != NULL))
    return pcurve_->endparam();
  else
    return spacecurve_->endparam();
}

//===========================================================================
void CurveOnVolume::reverseParameterDirection(bool switchparam)
//===========================================================================
{
  if (pcurve_.get() != 0)
    pcurve_->reverseParameterDirection();
  if (spacecurve_.get()!= 0)
    spacecurve_->reverseParameterDirection();
}

//===========================================================================
void CurveOnVolume::setParameterInterval(double t1, double t2)
//===========================================================================
{
  if (pcurve_.get() != 0)
    pcurve_->setParameterInterval(t1, t2);
  if (spacecurve_.get()!= 0)
    spacecurve_->setParameterInterval(t1, t2);
}

//===========================================================================
bool CurveOnVolume::isDegenerate(double degenerate_epsilon)
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
SplineCurve* CurveOnVolume::geometryCurve()
//===========================================================================
{
  if (spacecurve_.get() == 0)
    {
      double tol = 1.0e-4; // Arbitrary
      spacecurve_ = VolumeTools::liftVolParamCurve(pcurve_, volume_, tol);
    }
  return spacecurve_->geometryCurve();
}


//===========================================================================
CurveOnVolume* CurveOnVolume::subCurve(double from_par, 
				       double to_par, double fuzzy) const
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
	Point from_pt, to_pt;
	volume_->point(from_pt, from_par_pt[0], from_par_pt[1], from_par_pt[2]);
	volume_->point(to_pt, to_par_pt[0], to_par_pt[1], to_par_pt[2]);
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
	    shared_ptr<ParamCurve>(spacecurve_->subCurve(from_par, to_par,
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
	double clo_u_from, clo_v_from, clo_w_from, clo_u_to, clo_v_to, clo_w_to;
	double clo_dist1, clo_dist2;
	Point clo_pt, clo_pt_from, clo_pt_to;
	Point from_seed = pcurve_->point(from_par);
	Point to_seed = pcurve_->point(to_par);
	volume_->closestPoint(from_space_pt, clo_u_from, clo_v_from,
			      clo_w_from, clo_pt_from, clo_dist1, fuzzy,
			      from_seed.begin());
	volume_->closestPoint(to_space_pt, clo_u_to, clo_v_to,
			      clo_w_to, clo_pt_to, clo_dist2, fuzzy,
			      to_seed.begin());
	Point from_par_pt(clo_u_from, clo_v_from, clo_w_from);
	Point to_par_pt(clo_u_to, clo_v_to, clo_w_to);

	double clo_from, clo_to;
	pcurve_->closestPoint(from_par_pt, pcurve_->startparam(),
			      pcurve_->endparam(), clo_from,
			      clo_pt, clo_dist1, &from_par);
	pcurve_->closestPoint(to_par_pt, pcurve_->startparam(),
			      pcurve_->endparam(), clo_to,
			      clo_pt, clo_dist2, &to_par);
	try {
	  // subpcurve =
	  // 	shared_ptr<ParamCurve>(pcurve_->subCurve(clo_from, clo_to,
	  // 						   fuzzy));
	  subpcurve =
	    shared_ptr<ParamCurve>(pcurve_->subCurve(from_par, to_par,
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
  
  CurveOnVolume *sub_cv = new CurveOnVolume(volume_, subpcurve, 
					    subspacecurve,
					    prefer_parameter_);

  return sub_cv;
}

//===========================================================================
vector<shared_ptr<ParamCurve> >  CurveOnVolume::split(double param,
						      double fuzzy) const
//===========================================================================
{
  // Simple and inefficient solution
  vector<shared_ptr<ParamCurve> > sub_cvs(2);
  sub_cvs[0] = shared_ptr<ParamCurve>(subCurve(startparam(), param));
  sub_cvs[1] = shared_ptr<ParamCurve>(subCurve(param, endparam()));
  return sub_cvs;
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
void CurveOnVolume::closestPoint(const Point&   pt,
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
double CurveOnVolume::length(double tol)
//===========================================================================
{
    shared_ptr<SplineCurve> geom_cv(geometryCurve());
    assert(geom_cv.get() != NULL);
    return geom_cv->length(tol);
}

//===========================================================================
void CurveOnVolume::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
  // Must be implemented
}


//===========================================================================
void CurveOnVolume::appendCurve(ParamCurve* other_curve,
				int continuity, double& dist, bool reparam)
//===========================================================================
{
  // Must be implemented
}

//===========================================================================
RectDomain CurveOnVolume::containingDomain() const
//===========================================================================
{
  // Does not make sense
  RectDomain dummy;
  return dummy;
}

//===========================================================================
double  CurveOnVolume::nextSegmentVal(double par, bool forward, double tol) const
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
Point 
CurveOnVolume::volumeParameter(double crv_par,
			       const RectDomain* domain_of_interest) const
//===========================================================================
{
  Point param;
  if (pcurve_)
    param = pcurve_->point(crv_par);
  else
    {
      Point pos = spacecurve_->point(crv_par);
      double clo_u, clo_v, clo_w, dist;
      Point clo_pt;
      double eps = 1.0e-6;
      volume_->closestPoint(pos, clo_u, clo_v, clo_w, clo_pt, dist, eps);
      param = Point(clo_u, clo_v, clo_w);
    }
  return param;
}
