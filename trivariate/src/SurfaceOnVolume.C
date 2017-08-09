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

#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/trivariate/CurveOnVolume.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/Factory.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;
using std::streamsize;
using std::endl;

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume()
  : prefer_parameter_(false), constdir_(0), constval_(0.0),
    at_bd_(-1), orientation_(-1), swap_(false)
//===========================================================================
{
}

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume(shared_ptr<ParamVolume> vol,
				 shared_ptr<ParamSurface> parsurf,
				 shared_ptr<ParamSurface> spacesurf,
				 bool preferparameter)
  : volume_(vol), psurf_(parsurf), spacesurf_(spacesurf), 
    prefer_parameter_(preferparameter), constdir_(0), constval_(0.0),
    at_bd_(-1), orientation_(-1), swap_(false)
//===========================================================================
{
}

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume(shared_ptr<ParamVolume> vol,
				 shared_ptr<ParamSurface> spacesurf,
				 int constdir, double constpar, int boundary,
				 bool swapped, int orientation)
  : volume_(vol), spacesurf_(spacesurf), prefer_parameter_(false), 
    constdir_(constdir), constval_(constpar), at_bd_(boundary), 
    orientation_(orientation), swap_(swapped)
//===========================================================================
{
  // Make parameter surface
  // @@@ This construction may not be sufficient when trimmed volumes are
  // introduced
  if (constdir > 0)
    {
      Array<double,6> domain = volume_->parameterSpan();
      RectDomain pardom = spacesurf_->containingDomain();
      vector<double> knots1(4);
      vector<double> knots2(4);
      knots1[0] = knots1[1] = (swap_) ? pardom.vmin() : pardom.umin();
      knots1[2] = knots1[3] = (swap_) ? pardom.vmax() : pardom.umax();
      knots2[0] = knots2[1] = (swap_) ? pardom.umin() : pardom.vmin();
      knots2[2] = knots2[3] = (swap_) ? pardom.umax() : pardom.vmax();
      vector<double> coefs(12);
      domain[2*(constdir_-1)] = domain[2*(constdir_-1)+1] = constval_;
      coefs[0] = domain[0];
      coefs[1] = domain[2];
      coefs[2] = domain[4];
      coefs[3] = domain[1];
      coefs[4] = (constdir_ <= 2) ? domain[3] : domain[2];
      coefs[5] = domain[4];
      coefs[6] = domain[0];
      coefs[7] = (constdir_ <= 2) ? domain[2] : domain[3];
      coefs[8] = domain[5];
      coefs[9] = domain[1];
      coefs[10] = domain[3];
      coefs[11] = domain[5];

      psurf_ = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0],
							  &knots2[0], &coefs[0], 3));
      if (swap_)
	{
	  psurf_->swapParameterDirection();
	}
    }
}

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume(shared_ptr<ParamVolume> vol,
				 shared_ptr<ParamSurface> spacesurf,
				 shared_ptr<ParamSurface> parsurf,
				 bool prefer_parameter,
				 int constdir, double constpar, int boundary,
				 bool swapped)
  : volume_(vol), psurf_(parsurf), spacesurf_(spacesurf),
    prefer_parameter_(prefer_parameter), 
    constdir_(constdir), constval_(constpar), at_bd_(boundary), 
    orientation_(0), swap_(swapped)
//===========================================================================
{
}

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume(shared_ptr<ParamVolume> vol,
				 int constdir, double constpar, int boundary)
  : volume_(vol), prefer_parameter_(true), constdir_(constdir), 
    constval_(constpar), at_bd_(boundary), orientation_(0), swap_(false)
//===========================================================================
{
  // Fetch constant parameter surface
  spacesurf_ = 
    shared_ptr<ParamSurface>(volume_->constParamSurface(constval_, constdir_-1));

  // Make parameter surface
  // @@@ This construction may not be sufficient when trimmed volumes are
  // introduced
  Array<double,6> domain = volume_->parameterSpan();
  RectDomain pardom = spacesurf_->containingDomain();
  vector<double> knots1(4);
  vector<double> knots2(4);
  knots1[0] = knots1[1] = pardom.umin();
  knots1[2] = knots1[3] = pardom.umax();
  knots2[0] = knots2[1] = pardom.vmin();
  knots2[2] = knots2[3] = pardom.vmax();
  vector<double> coefs(12);
  domain[2*constdir_] = domain[2*constdir_+1] = constval_;
  coefs[0] = domain[0];
  coefs[1] = domain[2];
  coefs[2] = domain[4];
  coefs[3] = domain[1];
  coefs[4] = (constdir_ <= 2) ? domain[3] : domain[2];
  coefs[5] = domain[4];
  coefs[6] = domain[0];
  coefs[7] = (constdir_ <= 2) ? domain[2] : domain[3];
  coefs[8] = domain[5];
  coefs[9] = domain[1];
  coefs[10] = domain[3];
  coefs[11] = domain[5];

  psurf_ = shared_ptr<ParamSurface>(new SplineSurface(2, 2, 2, 2, &knots1[0],
						      &knots2[0], &coefs[0], 3));
}

//===========================================================================
SurfaceOnVolume::SurfaceOnVolume(const SurfaceOnVolume& other)
  : volume_(other.volume_), prefer_parameter_(other.prefer_parameter_), 
    constdir_(other.constdir_), constval_(other.constval_), 
    at_bd_(other.at_bd_), orientation_(other.orientation_),
    swap_(other.swap_)
//===========================================================================
{
  // Clones the surfaces
  if (other.psurf_.get())
    psurf_ = shared_ptr<ParamSurface>(other.psurf_->clone());

  if (other.spacesurf_.get())
    spacesurf_ = shared_ptr<ParamSurface>(other.spacesurf_->clone());
}

//===========================================================================
SurfaceOnVolume::~SurfaceOnVolume()
//===========================================================================
{}

//===========================================================================
void SurfaceOnVolume::read(std::istream& is)
//===========================================================================
{
  bool is_good = is.good();
  if (!is_good) {
    THROW("Invalid geometry file!");
  }
  // Do not care about surface...
  ALWAYS_ERROR_IF(psurf_.get() != NULL,
		  "Parameter surf already exists!");

  ALWAYS_ERROR_IF(spacesurf_.get() != NULL,
		  "Space surf already exists!");

  bool prefer_parameter;
  int  prefer_parameter_int;
  int  psurf_type;
  int  spacesurf_type;
  shared_ptr<ParamSurface> psurf;
  shared_ptr<ParamSurface> spacesurf;

  is >> prefer_parameter_int;
  if (prefer_parameter_int == 0)
    prefer_parameter = false;
  else if (prefer_parameter_int == 1)
    prefer_parameter = true;
  else 
    THROW("Unknown input for preferred SurfaceOnVolume parameter");

  is >> psurf_type;
  is >> spacesurf_type;

  if (psurf_type == 0) {
    // Do nothing - continue
  }
  else {
    ClassType type = ClassType(psurf_type); // Needs this conversion
    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    psurf = dynamic_pointer_cast<ParamSurface, GeomObject>(goobject);
    ALWAYS_ERROR_IF(psurf.get() == 0,
		    "Can not read this instance type");
    psurf->read(is);
  }

  if (spacesurf_type == 0) {
    // Do nothing - continue
  }
  else {
    ClassType type = ClassType(spacesurf_type); // Needs this conversion
    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    spacesurf = dynamic_pointer_cast<ParamSurface, GeomObject>(goobject);
    ALWAYS_ERROR_IF(spacesurf.get() == 0,
		    "Can not read this instance type");
    spacesurf->read(is);
  }

  prefer_parameter_ = prefer_parameter;
  psurf_ = psurf;
  spacesurf_ = spacesurf;
  is >> constdir_;
  is >> constval_;
  is >> at_bd_;
  is >> orientation_;
  is >> swap_;
}

//===========================================================================
void SurfaceOnVolume::write(std::ostream& os) const
//===========================================================================
{
  streamsize prev = os.precision(15);

  // We do not write the volume. It must be reconnected while reading
    if (!prefer_parameter_)
	os << "0";
    else
	os << "1";
    os << ' ';

    if (psurf_.get() == NULL)
	os << "0";
    else
	os << psurf_->instanceType();
    os << ' ';

    if (spacesurf_.get() == NULL)
	os << "0";
    else
	os << spacesurf_->instanceType();
    os << endl;

    if (psurf_.get() != NULL)
	psurf_->write(os);

    if (spacesurf_.get() != NULL)
	spacesurf_->write(os);

    os << constdir_ << " " << constval_ << " " << at_bd_ << " ";
    os << orientation_ << " " << swap_ << std::endl;
    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
BoundingBox SurfaceOnVolume::boundingBox() const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->boundingBox();
  else
    return volume_->boundingBox();  // Way too big, not a good solution
}

//===========================================================================
int SurfaceOnVolume::dimension() const
//===========================================================================
{
  if (volume_.get())
    return volume_->dimension();
  else if (spacesurf_.get())
    return spacesurf_->dimension();
  else
    return 3;  // Should not occur
}

//===========================================================================
ClassType SurfaceOnVolume::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
const Domain& SurfaceOnVolume::parameterDomain() const
//===========================================================================
{
  if (prefer_parameter_ && (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)))
    return psurf_->parameterDomain();
  else
    return spacesurf_->parameterDomain();
}

//===========================================================================
RectDomain SurfaceOnVolume::containingDomain() const
//===========================================================================
{
  if (prefer_parameter_ && (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)))
    return psurf_->containingDomain();
  else
    return spacesurf_->containingDomain();
}

//===========================================================================
bool SurfaceOnVolume::inDomain(double u, double v, double eps) const
//===========================================================================
{
  if (prefer_parameter_&& (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)) )
    return psurf_->inDomain(u,v,eps);
  else
    return spacesurf_->inDomain(u,v,eps);
}

//===========================================================================
int SurfaceOnVolume::inDomain2(double u, double v, double eps) const
//===========================================================================
{
  if (prefer_parameter_&& (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)) )
    return psurf_->inDomain2(u,v,eps);
  else
    return spacesurf_->inDomain2(u,v,eps);
}

//===========================================================================
bool SurfaceOnVolume::onBoundary(double u, double v, double eps) const
//===========================================================================
{
  if (prefer_parameter_ && (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)))
    return psurf_->onBoundary(u,v,eps);
  else
    return spacesurf_->onBoundary(u,v,eps);
}

//===========================================================================
Point SurfaceOnVolume::closestInDomain(double u, double v) const
//===========================================================================
{
  if (prefer_parameter_ && (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)))
    return psurf_->closestInDomain(u,v);
  else
    return spacesurf_->closestInDomain(u,v);
}

//===========================================================================
CurveLoop SurfaceOnVolume::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
  // NB! Looses back pointer information
  if (spacesurf_.get())
    return spacesurf_->outerBoundaryLoop(degenerate_epsilon);
  else if (constdir_ > 0)
    {
      // Pick constant parameter surface
      shared_ptr<ParamSurface> csf = 
	shared_ptr<ParamSurface>(volume_->constParamSurface(constval_, constdir_-1));
      return csf->outerBoundaryLoop(degenerate_epsilon);
    }
  else
    {
      // Approximate boundary curves. First fetch loop in parameter domain
      CurveLoop ploop = psurf_->outerBoundaryLoop(degenerate_epsilon);

      // Local approximation tolerance
      double eps = std::max(degenerate_epsilon, 1.0e-4);

      vector<shared_ptr<ParamCurve> > loop_cvs(ploop.size());
      for (int ki=0; ki<ploop.size(); ++ki)
	{
	  loop_cvs[ki] = VolumeTools::liftVolParamCurve(ploop[ki], volume_, eps);
	}

      CurveLoop bd_loop(loop_cvs, degenerate_epsilon);
      return bd_loop;
    }
}

//===========================================================================
vector<CurveLoop> SurfaceOnVolume::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
  // NB! Looses back pointer information
  if (spacesurf_.get())
    return spacesurf_->allBoundaryLoops(degenerate_epsilon);
  else if (constdir_ > 0)
    {
      // Pick constant parameter surface
      shared_ptr<ParamSurface> csf = 
	shared_ptr<ParamSurface>(volume_->constParamSurface(constval_, constdir_-1));
      return csf->allBoundaryLoops(degenerate_epsilon);
    }
  else
    {
      // Approximate boundary curves. First fetch loop in parameter domain
      vector<CurveLoop> ploops = psurf_->allBoundaryLoops(degenerate_epsilon);

      // Local approximation tolerance
      double eps = std::max(degenerate_epsilon, 1.0e-4);

      vector<CurveLoop> bd_loops;
      for (size_t kj=0; kj<ploops.size(); ++kj)
	{
	  vector<shared_ptr<ParamCurve> > loop_cvs(ploops[kj].size());
	  for (int ki=0; ki<ploops[kj].size(); ++ki)
	    {
	      loop_cvs[ki] = VolumeTools::liftVolParamCurve(ploops[kj][ki], volume_, eps);
	    }
	  
	  CurveLoop bd_loop(loop_cvs, degenerate_epsilon);
	  bd_loops.push_back(bd_loop);
	}
      return bd_loops;
    }
}

//===========================================================================
DirectionCone SurfaceOnVolume::normalCone() const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->normalCone();
  else
    THROW("normalCone only supported when the space surface is given");
}

//===========================================================================
DirectionCone SurfaceOnVolume::tangentCone(bool pardir_is_u) const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->tangentCone(pardir_is_u);
  else
    THROW("tangentCone only supported when the space surface is given");
}

//===========================================================================
CompositeBox SurfaceOnVolume::compositeBox() const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->compositeBox();
  else
    THROW("compositeBox only supported when the space surface is given");
}

//===========================================================================
void SurfaceOnVolume::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
  if (spacesurf_.get())
    {
      // Evaluate the space surface
      spacesurf_->point(pt, upar, vpar);
    }
  else if (constdir_ > 0)
    {
      // Evaluate the volume
      double tmp[3];
      if (swap_)
	std::swap(upar, vpar);   // Safe due to transfer by value
      tmp[constdir_-1] = constval_;
      tmp[(constdir_==1) ? 1 : 0] = upar;
      tmp[(constdir_<=2) ? 2 : 1] = vpar;
      volume_->point(pt, tmp[0], tmp[1], tmp[2]);
    }
  else 
    {
      // Evaluate the parameter surface and then the volume
      Point tmp;
      psurf_->point(tmp, upar, vpar);
      volume_->point(pt, tmp[0], tmp[1], tmp[2]);
    }
}

//===========================================================================
void SurfaceOnVolume:: point(std::vector<Point>& pts, 
			     double upar, double vpar,
			     int derivs,
			     bool u_from_right,
			     bool v_from_right,
			     double resolution) const
//===========================================================================
{
  if (spacesurf_.get())
      // Evaluate the space surface
      spacesurf_->point(pts, upar, vpar, derivs, u_from_right, v_from_right,
			resolution);
  else if (constdir_ > 0)
    {
      // Evaluate the volume
      double tmp[3];
      if (swap_)
	std::swap(upar, vpar);   // Safe due to transfer by value
      tmp[constdir_-1] = constval_;
      tmp[(constdir_==1) ? 1 : 0] = upar;
      tmp[(constdir_<=2) ? 2 : 1] = vpar;
      bool from_right[3];
      from_right[constdir_-1] = true;
      from_right[(constdir_==1) ? 1 : 0] = u_from_right;
      from_right[(constdir_<=2) ? 2 : 1] = v_from_right;
      vector<Point> tmp_pts((derivs+1)*(derivs+2)*(derivs+3)/6);
      volume_->point(tmp_pts, tmp[0], tmp[1], tmp[2], derivs, from_right[0], 
		     from_right[1], from_right[2], resolution);

      int ki, kj, idx1, idx2;
      int del1 = 1;
      int del2 = 1;
      for (ki=0, idx1=0, idx2=0; ki<=derivs; ki++, idx1+=del1)
	{
	  for (kj=0, del2=1; kj<=ki; kj++)
	    {
	      pts[idx2++] = tmp_pts[idx1];
	      if (constdir_ == 2 || (constdir_ == 3 && kj>0))
		del2++;
	      if (kj < ki)
		idx1 += del2;
	    }
	  if (constdir_ == 1)
	    del1 += (ki+1);
	  else if (constdir_ == 3 && ki>0)
	    del1++;
	}
	
    }
  else 
    {
      // Evaluate the parameter surface and then the volume
      vector<Point> tmp((derivs+1)*(derivs+2)/2);
      psurf_->point(tmp, upar, vpar, derivs);
      bool from_right[3];
      from_right[constdir_-1] = true;
      from_right[(constdir_==1) ? 1 : 0] = u_from_right;
      from_right[(constdir_<=2) ? 2 : 1] = v_from_right;
      vector<Point> Vd((derivs+1)*(derivs+2)*(derivs+3)/6);
      volume_->point(Vd, tmp[0][0], tmp[0][1], tmp[0][2], 
		     derivs, from_right[0], from_right[1], 
		     from_right[2], resolution);

      pts[0] = Vd[0];
      if (derivs >= 1)
	{
	  pts[1] = Vd[1]*tmp[1][0] + Vd[2]*tmp[1][1] + 
	    Vd[3]*tmp[1][2];
	  pts[2] = Vd[1]*tmp[2][0] + Vd[2]*tmp[2][1] + 
	    Vd[3]*tmp[2][2];
	  if (derivs >= 2)
	    {
	      pts[3] = Vd[4]*tmp[1][0]*tmp[1][0] +
		2.0*Vd[5]*tmp[1][0]*tmp[1][1] +
		2.0*Vd[6]*tmp[1][0]*tmp[1][2] +
		Vd[7]*tmp[1][1]*tmp[1][1] +
		2.0*Vd[8]*tmp[1][1]*tmp[1][2] +
		Vd[9]*tmp[1][2]*tmp[1][2] +
		Vd[1]*tmp[4][0] + Vd[4]*tmp[2][1] + 
		Vd[3]*tmp[4][2];

	      pts[4] = Vd[4]*tmp[1][0]*tmp[2][0] +
		Vd[5]*(tmp[1][0]*tmp[2][1] + tmp[2][0]*tmp[1][1]) +
		Vd[6]*(tmp[1][0]*tmp[2][2] + tmp[2][0]*tmp[1][2]) +
		Vd[7]*tmp[1][1]*tmp[2][1] +
		Vd[8]*(tmp[1][1]*tmp[2][2] + tmp[2][1]*tmp[1][2]) +
		Vd[9]*tmp[1][2]*tmp[2][2] +
		Vd[1]*tmp[4][0] + Vd[4]*tmp[4][1] +
		Vd[3]*tmp[4][2];

	      pts[5] = Vd[4]*tmp[2][0]*tmp[2][0] +
		2.0*Vd[5]*tmp[2][0]*tmp[2][1] +
		2.0*Vd[6]*tmp[2][0]*tmp[2][2] +
		Vd[7]*tmp[2][1]*tmp[2][1] +
		2.0*Vd[8]*tmp[2][1]*tmp[2][2] +
		Vd[9]*tmp[2][2]*tmp[2][2] +
		Vd[1]*tmp[5][0] + Vd[2]*tmp[5][1] + 
		Vd[3]*tmp[5][2];
	    }
	}
    }
        
}


//===========================================================================
void SurfaceOnVolume::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
  vector<Point> res(3);
  point(res, upar, vpar, 1);
  n = res[1].cross(res[2]);
}

//===========================================================================
vector<shared_ptr<ParamCurve> > SurfaceOnVolume::constParamCurves(double parameter,
								  bool pardir_is_u) const
//===========================================================================
{
  // Local approximation tolerance
  double eps = 1.0e-4;
  vector<shared_ptr<CurveOnVolume> > vol_cvs;
  if (spacesurf_.get())
    {
      vector<shared_ptr<ParamCurve> > cvs1 = 
	spacesurf_->constParamCurves(parameter, pardir_is_u);
      for (size_t ki=0; ki<cvs1.size(); ++ki)
	{
	  // shared_ptr<ParamCurve> cv2 = 
	  //   VolumeTools::projectVolParamCurve(cvs1[ki], volume_, eps);
	  // shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, cv2,
	  // 						  cvs1[ki],
	  // 						  prefer_parameter_));
	  shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, 
							  cvs1[ki],
							  false));
	  vol_cvs.push_back(crv);
	}
    }
  else if (constdir_ > 0)
    {
      // Pick constant parameter surface
      shared_ptr<ParamSurface> csf = 
	shared_ptr<ParamSurface>(volume_->constParamSurface(constval_, constdir_-1));
      vector<shared_ptr<ParamCurve> > cvs1 = 
	csf->constParamCurves(parameter, pardir_is_u);
      for (size_t ki=0; ki<cvs1.size(); ++ki)
	{
	  // shared_ptr<ParamCurve> cv2 = 
	  //   VolumeTools::projectVolParamCurve(cvs1[ki], volume_, eps);
	  // shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, cv2,
	  // 						  cvs1[ki],
	  // 						  prefer_parameter_));
	  shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, 
							  cvs1[ki],
							  false));
	  vol_cvs.push_back(crv);
	}
    }
  else
    {
      // Pick constant parameter curves in parameter surface
      vector<shared_ptr<ParamCurve> > pcvs = 
	psurf_->constParamCurves(parameter, pardir_is_u);

      // Make approximated geometry curves
      for (size_t ki=0; ki<pcvs.size(); ++ki)
	{
	  // shared_ptr<ParamCurve> cv2 = 
	  //   VolumeTools::liftVolParamCurve(pcvs[ki], volume_, eps);
	  // shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, pcvs[ki],
	  // 						  cv2,
	  // 						  prefer_parameter_));
	  shared_ptr<CurveOnVolume> crv(new CurveOnVolume(volume_, pcvs[ki],
							  true));
	  vol_cvs.push_back(crv);
	}
    }
  vector<shared_ptr<ParamCurve> > return_cvs(vol_cvs.begin(), vol_cvs.end());
  return return_cvs;
}
    

//===========================================================================
vector<shared_ptr<ParamSurface> > SurfaceOnVolume::subSurfaces(double from_upar, 
							       double from_vpar,
							       double to_upar, 
							       double to_vpar,
							       double fuzzy) const
//===========================================================================
{
  vector<shared_ptr<ParamSurface> > sub_sfs;
  vector<shared_ptr<ParamSurface> > p_sfs;
  vector<shared_ptr<ParamSurface> > geom_sfs;

  if (psurf_.get())
    p_sfs = psurf_->subSurfaces(from_upar, from_vpar, to_upar, to_vpar, fuzzy);
  if (spacesurf_.get())
    geom_sfs = spacesurf_->subSurfaces(from_upar, from_vpar, to_upar, to_vpar, 
				       fuzzy);

  size_t nmb = geom_sfs.size();
  for (size_t ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> surf;
      if (p_sfs.size() == nmb)
	 surf = shared_ptr<ParamSurface>(new SurfaceOnVolume(volume_, geom_sfs[ki],
							     p_sfs[ki], 
							     prefer_parameter_,
							     constdir_, constval_,
							     at_bd_, swap_));
      else
	 surf = shared_ptr<ParamSurface>(new SurfaceOnVolume(volume_, geom_sfs[ki],
							     constdir_, constval_,
							     at_bd_, swap_));
      sub_sfs.push_back(surf);
    }
  return sub_sfs;
}



//===========================================================================
double SurfaceOnVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
  if (prefer_parameter_ && psurf_.get())
    return psurf_->nextSegmentVal(dir, par, forward, tol);
  else if (spacesurf_.get())
    return spacesurf_->nextSegmentVal(dir, par, forward, tol);
  else
    {
      RectDomain dom = containingDomain();
      if (forward)
	return (dir == 0) ? dom.umax() : dom.vmax();
      else
	return (dir == 0) ? dom.umin() : dom.vmin();
    }
}

//===========================================================================
 void SurfaceOnVolume::closestPoint(const Point& pt,
				    double&        clo_u,
				    double&        clo_v, 
				    Point&       clo_pt,
				    double&        clo_dist,
				    double         epsilon,
				    const RectDomain* domain_of_interest,
				    double   *seed) const
//===========================================================================
{
  ParamSurface::closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
			     domain_of_interest, seed);
}

//===========================================================================
 void SurfaceOnVolume::closestBoundaryPoint(const Point& pt,
					    double&        clo_u,
					    double&        clo_v, 
					    Point&       clo_pt,
					    double&        clo_dist,
					    double         epsilon,
					    const RectDomain* domain_of_interest,
					    double   *seed) const
//===========================================================================
{
  // Must be implemented properly, for the time being ...
  if (spacesurf_.get())
    spacesurf_->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
				     domain_of_interest, seed);
  else
    {
      // Fetch boundary loops
      vector<CurveLoop> loops = allBoundaryLoops(epsilon);
      
      // Find closest point
      int min_loop = -1;
      int min_cv = -1;
      clo_dist = 1.0e10;  // Initialize with a large number
      double clo_par;
      for (size_t ki=0; ki<loops.size(); ++ki)
	for (int kj=0; kj<loops[ki].size(); ++kj)
	  {
	    double dist, clo_t;
	    Point curr_pt;
	    shared_ptr<ParamCurve> crv = loops[ki][kj];
	    crv->closestPoint(pt, crv->startparam(), crv->endparam(),
			      clo_t, clo_pt, dist);
	    if (dist < clo_dist)
	      {
		clo_dist = dist;
		min_loop = (int)ki;
		min_cv = kj;
		clo_par = clo_t;
		clo_pt = curr_pt;
	      }
	  }

      // Fetch information about closest point
      shared_ptr<ParamCurve> crv = loops[min_loop][min_cv];
      shared_ptr<CurveOnSurface> sf_cv =
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
      if (sf_cv.get() && sf_cv->parameterCurve().get())
	{
	  Point par = sf_cv->parameterCurve()->point(clo_par);
	  clo_u = par[0];
	  clo_v = par[1];
	}
      else
	{
	  Point clo_pt2 = clo_pt;
	  closestPoint(clo_pt2, clo_u, clo_v, clo_pt, clo_dist, epsilon,
		       domain_of_interest, seed);
	}
    }  
}

//===========================================================================
 void SurfaceOnVolume::getBoundaryInfo(Point& pt1, Point& pt2,
				       double epsilon, SplineCurve*& cv,
				       SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
  // Must be implemented
}

//===========================================================================
 void SurfaceOnVolume::turnOrientation()
//===========================================================================
{
  // Must be implemented
}

//===========================================================================
 void SurfaceOnVolume::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
  // Must be implemented
}

//===========================================================================
 void SurfaceOnVolume::swapParameterDirection()
//===========================================================================
{
  if (spacesurf_.get())
    spacesurf_->swapParameterDirection();

  if (psurf_.get())
    psurf_->swapParameterDirection();

  swap_ = !swap_;
}

//===========================================================================
 double SurfaceOnVolume::area(double tol) const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->area(tol);
  else
    THROW("Area not supported in this case");
}

//===========================================================================
 bool SurfaceOnVolume::isDegenerate(bool& b, bool& r,
				    bool& t, bool& l, double tolerance) const
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->isDegenerate(b, r, t, l, tolerance);
  else
    return ParamSurface::isDegenerate(b, r, t, l, tolerance);
}

//===========================================================================
void SurfaceOnVolume::getDegenerateCorners(vector<Point>& deg_corners, 
					   double tol) const
//===========================================================================
{
  if (spacesurf_.get())
    spacesurf_->getDegenerateCorners(deg_corners, tol);
  else
    {
      vector<Point> derivs(3);
      double ang;
      vector<pair<Point,Point> > corners;
      psurf_->getCornerPoints(corners);
      for (size_t ki=0; ki<corners.size(); ++ki)
	{
	  Point parval = corners[ki].second;
	  point(derivs, parval[0], parval[1], 1);
	  ang = derivs[1].angle(derivs[2]);
	  if (fabs(ang) < tol || fabs(M_PI-ang) < tol)
	    deg_corners.push_back(parval);
	}
    }
	  
}

//===========================================================================
void SurfaceOnVolume::getCornerPoints(vector<pair<Point,Point> >& corners) const
//===========================================================================
{
  if (prefer_parameter_ && (!(spacesurf_.get() && 
			      spacesurf_->instanceType() == Class_BoundedSurface)))
    {
      vector<pair<Point,Point> > tmp_corners;
      psurf_->getCornerPoints(tmp_corners);
      corners.reserve(tmp_corners.size());
      for (size_t ki=0; ki<tmp_corners.size(); ++ki)
	{
	  Point par  = tmp_corners[ki].first;
	  Point pos;
	  volume_->point(pos, par[0], par[1], par[2]);
	  corners.push_back(std::make_pair(pos, tmp_corners[ki].second));
	}
    }
  else
    spacesurf_->getCornerPoints(corners);
}

//===========================================================================
 bool SurfaceOnVolume::isIsoTrimmed(double tol) const
//===========================================================================
{
  if (prefer_parameter_)
    return psurf_->isIsoTrimmed(tol);
  else
    return spacesurf_->isIsoTrimmed(tol);
}

//===========================================================================
int SurfaceOnVolume::whichBoundary(double tol, int& orientation, bool& swap) const
//===========================================================================
{
  if (constdir_ == 1 || constdir_ == 2 || constdir_ == 3)
    {
      orientation = orientation_;
      swap = swap_;
      return at_bd_;
    }
  else
    return -1;  // For the time being
}

//===========================================================================
Point SurfaceOnVolume:: volumeParameter(double u_par, double v_par) const
//===========================================================================
{
  Point param(3);
  RectDomain dom = containingDomain();
  if ((constdir_ == 1 || constdir_ == 2 || constdir_ == 3) &&
      orientation_ >= 0)
    {
      int idx1, idx2, idx3;
      idx1 = constdir_ - 1;
      idx2 = (idx1 == 0) ? 1 : 0;
      idx3 = (constdir_ == 3) ? 1 : 2;
      param[idx1] = constval_;
      if (swap_)
	{
	  param[idx2] = (orientation_ == 2 || orientation_ == 3) ?
	    dom.vmax() - (v_par - dom.vmin()) : v_par;
	  param[idx3] = (orientation_ == 1 || orientation_ == 3) ?
	    dom.umax() - (u_par - dom.umin()) : u_par;
	}
      else
	{
	  param[idx3] = (orientation_ == 2 || orientation_ == 3) ?
	    dom.vmax() - (v_par - dom.vmin()) : v_par;
	  param[idx2] = (orientation_ == 1 || orientation_ == 3) ?
	    dom.umax() - (u_par - dom.umin()) : u_par;
	}
    }
  else if (psurf_)
    param = psurf_->ParamSurface::point(u_par, v_par);
  else
    {
      // No parameter surface exist. Perform closest point computation
      Point pos = spacesurf_->ParamSurface::point(u_par, v_par);
      double clo_u, clo_v, clo_w, dist;
      Point clo_pt;
      double eps = 1.0e-6;
      volume_->closestPoint(pos, clo_u, clo_v, clo_w, clo_pt, dist, eps);
      param = Point(clo_u, clo_v, clo_w);
    }
  return param;
}

//===========================================================================
void SurfaceOnVolume::unsetParamSurf()
//===========================================================================
{
  shared_ptr<ParamSurface> dummy;
  psurf_ = dummy;
  prefer_parameter_ = false;

  constdir_ = 0;
  constval_ = 0.0;
  at_bd_ = -1;
  orientation_ = -1;
  swap_ = false;
}

//===========================================================================
bool SurfaceOnVolume::isLinear(Point& dir1, Point& dir2, double tol)
//===========================================================================
{
  if (spacesurf_.get())
    return spacesurf_->isLinear(dir1, dir2, tol);
  else
    {
      dir1.resize(0);
      dir2.resize(0);
      return false;
    }
}

