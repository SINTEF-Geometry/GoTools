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

#include "GoTools/geometry/BoundedCurve.h"
#include "GoTools/geometry/Factory.h"
#include <vector>


using std::vector;
using std::endl;
using std::istream;
using std::ostream;
using std::streamsize;

namespace Go {


//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   bool prefer_parameter,
			   double startpar, double endpar,
			   Point start_pt, Point end_pt)
    : curve_(curve), prefer_parameter_(prefer_parameter),
      startparam_(startpar), endparam_(endpar),
      start_pt_(start_pt), end_pt_(end_pt)
//===========================================================================
{
}


// Constructor. Input is point and direction
//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   Point start_pt,
			   Point end_pt)
    : curve_(curve), prefer_parameter_(false),
      start_pt_(start_pt), end_pt_(end_pt)
//===========================================================================
{
    startparam_ = curve_->startparam();
    endparam_ = curve_->endparam();
}


//===========================================================================
BoundedCurve::BoundedCurve(shared_ptr<ParamCurve> curve,
			   double startpar, double endpar)
  : curve_(curve), prefer_parameter_(true),
    startparam_(startpar), endparam_(endpar)
//===========================================================================
{
}

/// virtual destructor - ensures safe inheritance
//===========================================================================
BoundedCurve::~BoundedCurve()
//===========================================================================
{
}


/// Read object from stream
//===========================================================================
void BoundedCurve::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int curve_type;
    is >> curve_type;

    is >> prefer_parameter_;
    is >> startparam_;
    is >> endparam_;

    int dim;
    is >> dim;
    start_pt_.resize(dim);
    end_pt_.resize(dim);
    int i;
    for (i=0; i<dim; ++i)
	is >> start_pt_[i];    
    for (i=0; i<dim; ++i)
	is >> end_pt_[i];

    ClassType type = ClassType(curve_type); // Needs this conversion
    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    shared_ptr<ParamCurve> tmp_crv =
	dynamic_pointer_cast<ParamCurve, GeomObject>(goobject);
    ALWAYS_ERROR_IF(tmp_crv.get() == 0,
		    "Can not read this instance type");
    tmp_crv->read(is);
    curve_ = tmp_crv;
    
}


/// Write object to stream
//===========================================================================
void BoundedCurve::write(std::ostream& os) const
//===========================================================================
{
    int i;
    streamsize prev = os.precision(15);
    os << curve_->instanceType();
    os << endl;

    if (prefer_parameter_)
	os << 1 << " ";
    else 
	os << 0 << " ";
    os << startparam_ << " " << endparam_ << endl;

    int dim = start_pt_.dimension();
    os << dim << "  ";
    for (i=0; i<dim; ++i)
	os << start_pt_[i] << "  ";
    for (i=0; i<dim; ++i)
	os << end_pt_[i] << "  ";
    os << endl;

    curve_->write(os);
    os.precision(prev);
}


//===========================================================================
BoundingBox BoundedCurve::boundingBox() const
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->boundingBox();
}
    

//===========================================================================
int BoundedCurve::dimension() const
//===========================================================================
{
    return curve_->dimension();
}


//===========================================================================
ClassType BoundedCurve::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType BoundedCurve::classType()
//===========================================================================
{
    return Class_BoundedCurve;
}


//===========================================================================
BoundedCurve* BoundedCurve::clone() const
//===========================================================================
{
    BoundedCurve* bd_cv = new BoundedCurve(curve_, prefer_parameter_,
					   startparam_, endparam_,
					   start_pt_, end_pt_);

    return bd_cv;
}


//===========================================================================
void BoundedCurve::point(Point& pt, double tpar) const
//===========================================================================
{
    curve_->point(pt, tpar);
}
 

//===========================================================================
void BoundedCurve::point(std::vector<Point>& pts, 
			 double tpar,
			 int derivs,
			 bool from_right) const
//===========================================================================
{
    curve_->point(pts, tpar, derivs, from_right);
}

 
//===========================================================================
double BoundedCurve::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double BoundedCurve::endparam() const
//===========================================================================
{
    return endparam_;
}

 
//===========================================================================
void BoundedCurve::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    curve_->reverseParameterDirection(switchparam);
    std::swap(start_pt_, end_pt_);
}
    

//===========================================================================
void BoundedCurve::setParameterInterval(double t1, double t2)
//===========================================================================
{
    setParamBounds(t1, t2);
}


//===========================================================================
SplineCurve* BoundedCurve::geometryCurve()
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));

    return sub_cv->geometryCurve();
}


//===========================================================================
bool BoundedCurve::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->isDegenerate(degenerate_epsilon);
}

     
//===========================================================================
BoundedCurve* BoundedCurve::subCurve(double from_par, double to_par,
				     double fuzzy) const
//===========================================================================
{
    BoundedCurve* bd_cv = clone();
    bd_cv->setParamBounds(from_par, to_par);

    return bd_cv;
}


//===========================================================================
DirectionCone BoundedCurve::directionCone() const
//===========================================================================
{
    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(startparam_, endparam_));
    return sub_cv->directionCone();
}
 

//===========================================================================
void BoundedCurve::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void BoundedCurve::appendCurve(ParamCurve* cv, int continuity, double& dist,
			       bool reparam, double tol)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void BoundedCurve::closestPoint(const Point& pt,
				double tmin,
				double tmax,
				double& clo_t,
				Point& clo_pt,
				double& clo_dist,
				double const *seed) const
//===========================================================================
{
    // Check and fix the parameter bounds
    if (tmin < startparam_) {
	tmin = startparam_;
	MESSAGE("tmin too small. Using startparam_.");
    }
    if (tmax > endparam_) {
	tmax = endparam_;
	MESSAGE("tmax too large. Using endparam_.");
    }

    shared_ptr<ParamCurve> sub_cv(curve_->subCurve(tmin, tmax));
    sub_cv->ParamCurve::closestPoint(pt, clo_t, clo_pt, clo_dist);
}


//===========================================================================
double BoundedCurve::length(double tol)
//===========================================================================
{
    return ParamCurve::length(tol, startparam_, endparam_);
}


//===========================================================================
bool BoundedCurve::isAxisRotational(Point& centre, Point& axis, Point& vec,
				    double& angle)
//===========================================================================
{
  double radius;
  return isAxisRotational(centre, axis, vec, angle, radius);
}

//===========================================================================
bool BoundedCurve::isAxisRotational(Point& centre, Point& axis, Point& vec,
				    double& angle, double& radius)
//===========================================================================
{
  bool rotational = curve_->isAxisRotational(centre, axis, vec, angle, radius);
  if (!rotational)
    return false;

  // The rotational angle may be reduced compared to the underlying curve.
  // Recompute
  Point pt1 = curve_->point(startparam_);
  Point pt2 = curve_->point(endparam_);
  vec = pt1 - centre;
  vec.normalize();
  Point vec2 = pt2 - centre;
  angle = vec.angle(vec2);
  
  return true;
}

//===========================================================================
  bool BoundedCurve::isInPlane(const Point& loc, const Point& axis,
			       double eps, Point& normal) const
//===========================================================================
{
  return curve_->isInPlane(loc, axis, eps, normal);
}

//===========================================================================
  bool BoundedCurve::isLinear(Point& dir, double tol)
//===========================================================================
{
  return curve_->isLinear(dir, tol);
}

//===========================================================================
void BoundedCurve::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
	THROW("First parameter must be strictly less than second.");

    startparam_ = startpar;
    endparam_ = endpar;
}


} // namespace Go
