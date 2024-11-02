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

#include "GoTools/geometry/Hyperbola.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>
#include <limits>


using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::streamsize;
using std::swap;


namespace Go
{


//===========================================================================
Hyperbola::Hyperbola(Point location, Point direction, Point normal,
		     double r1, double r2,
                     bool isReversed)
    : location_(location), vec1_(direction), normal_(normal),
      r1_(r1), r2_(r2)
//===========================================================================
{
    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }

    if (dimension() == 3)
	normal_.normalize();
    setSpanningVectors();

    double inf = numeric_limits<double>::infinity();
    setParamBounds(-inf, inf);
    setParameterInterval(-inf, inf);

    if (isReversed)
        reverseParameterDirection();
}


//===========================================================================
Hyperbola::~Hyperbola()
//===========================================================================
{
}


//===========================================================================
void Hyperbola::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    location_.resize(dim);
    normal_.resize(dim);
    vec1_.resize(dim);
    is >> r1_
       >> r2_
       >> location_
       >> normal_
       >> vec1_;

    if(dim == 3)
        normal_.normalize();
    setSpanningVectors();

    int isBounded; 
    is >> isBounded;
    bool has_param_int = (isBounded >= 10);
    isBounded = isBounded % 10;
    if (isBounded == 0) {
        // Unbounded - don't read parameters
        double inf = numeric_limits<double>::infinity();
        setParamBounds(-inf, inf);
    }
    else if (isBounded == 1) {
      is >> parbound1_ >> parbound2_;
    }
    else {
        THROW("Bounded flag must be 0 or 1");
    }
    startparam_ = parbound1_;
    endparam_ = parbound2_;

    if (has_param_int)
        is >> startparam_ >> endparam_;

    // "Reset" reversion
    isReversed_ = false;

    // Swapped flag
    int isReversed; // 0 or 1
    is >> isReversed;
    if (isReversed == 0) {
        // Do nothing
    }
    else if (isReversed == 1) {
        reverseParameterDirection();
    }
    else {
        THROW("Swapped flag must be 0 or 1");
    }
}


//===========================================================================
void Hyperbola::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    int dim = dimension();
    os << dim << endl
       << r1_ << endl
       << r2_ << endl
       << location_ << endl
       << normal_ << endl
       << vec1_ << endl;

    if (!isBounded()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl
	   << parbound1_ << parbound2_ << endl
           << startparam_ << endparam_ << endl;
	  
    }

    if (!isReversed()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
    }

    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
BoundingBox Hyperbola::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    // A rather inefficient hack...
    Hyperbola* hyperbola = const_cast<Hyperbola*>(this);
    SplineCurve* tmp = hyperbola->geometryCurve();
    if (tmp == NULL)
    {
	MESSAGE("Failed creating spline represention");
	return box;
    }
    box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
int Hyperbola::dimension() const
//===========================================================================
{
    return location_.dimension();
}
    

//===========================================================================
ClassType Hyperbola::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType Hyperbola::classType()
//===========================================================================
{
    return Class_Hyperbola;
}


//===========================================================================
Hyperbola* Hyperbola::clone() const
//===========================================================================
{
    Hyperbola* hyp = new Hyperbola(location_, vec1_, normal_, r1_, r2_,
        isReversed_);
    hyp->setParamBounds(parbound1_, parbound2_);
    hyp->setParameterInterval(startparam_, endparam_);
    return hyp;
}


//===========================================================================
void Hyperbola::point(Point& pt, double tpar) const
//===========================================================================
{
  getReversedParameter(tpar);
  tpar = parbound1_ + 
    (tpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
  pt = location_ + r1_*cosh(tpar)*vec1_ + r2_*sinh(tpar)*vec2_;
}


//===========================================================================
void Hyperbola::point(std::vector<Point>& pts, 
		      double tpar,
		      int derivs,
		      bool from_right) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0, 
		   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1);
    int ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz < totpts, 
		   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
	}
	pts[i].setValue(0.0);
    }

    point(pts[0], tpar);
    if (derivs == 0)
        return;

    // Since the hyperbola is parametrized as:
    // c(t) = location_ + r1_*cosh(t)*dir1_ + r2_*sinh(t)*dir2_,
    // the derivatives follow easily.
    getReversedParameter(tpar);
    double fac = (parbound2_-parbound1_)/(endparam_-startparam_);
    tpar = parbound1_ + fac*(tpar - startparam_);

    double sinh_t = sinh(tpar);
    double cosh_t = cosh(tpar);
    for (int ki = 1; ki < derivs + 1; ++ki) {
	pts[ki] = (ki%2 == 1) ? r1_*sinh_t*vec1_ + r2_*cosh_t*vec2_ :
	    r1_*cosh_t*vec1_ + r2_*sinh_t*vec2_;
        // Take reversion into account
        if (isReversed()) {
            double sgnrev = (ki % 2 == 1) ? -1.0 : 1.0;
            pts[ki] *= (sgnrev*fac);
	    fac *= (parbound2_-parbound1_)/(endparam_-startparam_);
        }
    }
}


//===========================================================================
double Hyperbola::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double Hyperbola::endparam() const
//===========================================================================
{
    return endparam_;
}


//===========================================================================
void Hyperbola::swapParameters2D()
//===========================================================================
{
    if (dimension() == 2) {
        swap(location_[0], location_[1]);
        swap(vec1_[0], vec1_[1]);
        swap(vec2_[0], vec2_[1]);
    }
}


//===========================================================================
void Hyperbola::setParameterInterval(double t1, double t2)
//===========================================================================
{
  startparam_ = t1;
  endparam_ = t2;
}


//===========================================================================
SplineCurve* Hyperbola::geometryCurve()
//===========================================================================
{
    return createSplineCurve();
}


//===========================================================================
SplineCurve* Hyperbola::createSplineCurve() const
//===========================================================================
{
    MESSAGE("createSplineCurve() not yet implemented.");
    return NULL;
}


//===========================================================================
bool Hyperbola::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    // We consider a Hyperbola as degenerate if either radii is smaller
    // than the epsilon.

    return ((r1_*vec1_.length() < degenerate_epsilon) ||
	    (r2_*vec2_.length() < degenerate_epsilon));
}



//===========================================================================
Hyperbola* Hyperbola::subCurve(double from_par, double to_par,
			  double fuzzy) const
//===========================================================================
{
  Hyperbola* hyperbola = clone();
  if (isReversed())
    {
      double start = endparam_ - (to_par - startparam_);
      double end = startparam_ + (endparam_ - from_par);
      if (start > end)
	{
	  std::swap(start, end);
	}
      double bound1 = parbound1_ + 
	(start-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      double bound2 = parbound1_ + 
	(end-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      hyperbola->setParamBounds(bound1, bound2);
      if (from_par > to_par)
	std::swap(from_par, to_par);
      hyperbola->setParameterInterval(from_par, to_par);
    }
  else
    {
      double bound1 = parbound1_ + 
	(from_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      double bound2 = parbound1_ + 
	(to_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
      hyperbola->setParamBounds(bound1, bound2);
      hyperbola->setParameterInterval(from_par, to_par);
    }
  return hyperbola;
}


//===========================================================================
DirectionCone Hyperbola::directionCone() const
//===========================================================================
{
    double tmin = startparam();
    double tmax = endparam();
    vector<Point> pts;
    point(pts, 0.5*(tmin+tmax), 1);
    // We must calculate the angle between the mid point and the end
    // points. As the curvature is monotone this gives the boundaries
    // for the tangents.
    Point start_pt, end_pt;
    point(start_pt, startparam_);
    point(end_pt, endparam_);
    Point dir1 = start_pt - location_;
    Point dir2 = end_pt - location_;
    Point dir3 = pts[0] - location_;
    double ang1 = dir1.angle(dir3);
    double ang2 = dir2.angle(dir3);
    return DirectionCone(pts[1], std::max(fabs(ang1), fabs(ang2)));
}
 

//===========================================================================
void Hyperbola::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    MESSAGE("appendCurve() not implemented!");
}


//===========================================================================
void Hyperbola::appendCurve(ParamCurve* cv,
			    int continuity, double& dist, bool reparam,
			    double tol)
//===========================================================================
{
    MESSAGE("appendCurve() not implemented!");
}


//===========================================================================
void Hyperbola::closestPoint(const Point& pt,
			   double tmin,
			   double tmax,
			   double& clo_t,
			   Point& clo_pt,
			   double& clo_dist,
			   double const *seed) const
//===========================================================================
{
    double guess_param = 0.5*(tmin + tmax);
    ParamCurve::closestPointGeneric(pt, tmin, tmax,
				    guess_param, clo_t, clo_pt, clo_dist);
}


//===========================================================================
double Hyperbola::length(double tol)
//===========================================================================
{
    if (!isBounded())
        return numeric_limits<double>::infinity();

    int num_spans = 4;

    double result = 0.0;
    double tstep = (endparam_ - startparam_)/(double)num_spans;
    for (int ki = 0; ki < num_spans; ++ki)
    {
	double start = ki*tstep;
	double end = (ki + 1)*tstep;
	result += ParamCurve::length(tol, start, end);
    }

    return result;
}


//===========================================================================
void Hyperbola::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
	THROW("First parameter must be strictly less than second.");

    double start =  
      parbound1_ + (startpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
    double end =  
      parbound1_ + (endpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
    parbound1_ = startpar;
    parbound2_ = endpar;
     startparam_ = start;
    endparam_ = end;
}


//===========================================================================
void Hyperbola::translateCurve(const Point& dir)
//===========================================================================
{
    location_ += dir;
}

//===========================================================================
bool Hyperbola::isBounded() const
//===========================================================================
{
    double inf = numeric_limits<double>::infinity();
    return parbound1_ > -inf && parbound2_ < inf;
}


//===========================================================================
void Hyperbola::setSpanningVectors()
//===========================================================================
{
    // In 3D, the spanning vectors vec1_, vec2_, and the vector
    // normal_ defines a right-handed coordinate system. Similar to an
    // axis2_placement_3d entity in STEP.

    int dim = location_.dimension();
    if (dim == 2) {
	vec2_.resize(2);
	vec2_[0] = -vec1_[1];
	vec2_[1] = vec1_[0];
    }
    else if (dim ==3) {
	Point tmp = vec1_ - (vec1_ * normal_) * normal_;
	if (tmp.length() == 0.0) 
	    THROW("X-axis parallel to normal.");
	vec1_ = tmp;
	vec2_ = normal_.cross(vec1_);
    }
    else {
	THROW("Dimension must be 2 or 3");
    }
    vec1_.normalize();
    vec2_.normalize();
}

//===========================================================================
bool Hyperbola::isInPlane(const Point& norm,
		       double eps, Point& pos) const

//===========================================================================
{
  double ang = norm.angle(normal_);
  pos = location_;

  return (ang <= eps || fabs(M_PI-ang) <= eps);
}


}
