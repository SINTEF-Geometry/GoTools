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

#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include <limits>


using std::vector;
using std::endl;
using std::numeric_limits;
using std::streamsize;
using std::swap;


namespace Go {


// Constructor. Input is point and direction
//===========================================================================
Line::Line(Point point, Point direction, bool isReversed)
    : location_(point), dir_(direction),
      parbound1_(-numeric_limits<double>::infinity()),
      parbound2_(numeric_limits<double>::infinity()),
      startparam_(-numeric_limits<double>::infinity()),
      endparam_(numeric_limits<double>::infinity())
//===========================================================================
{
    // Note: dir_ is not normalized.

    if (isReversed)
        reverseParameterDirection();
}


// Constructor to bounded line. Input is point, direction
// and length
//===========================================================================
Line::Line(Point point, Point direction, double length,
    bool isReversed)
  : location_(point), dir_(direction), parbound1_(0.0),
    parbound2_(length), startparam_(0.0), endparam_(length)
//===========================================================================
{
    double len = dir_.length();
    if (len > 0.0)
    {
        dir_.normalize();
    }

    if (isReversed)
        reverseParameterDirection();
}

// Constructor to bounded line. Input is two points with parameters.
//===========================================================================
Line::Line(Point point1, Point point2, double par1, double par2,
    bool isReversed)
  : parbound1_(par1), parbound2_(par2), startparam_(par1), endparam_(par2)
//===========================================================================
{
    // Throws if the two Points are equal.
  dir_ = point2 - point1;
  double len = dir_.length();
  // if (len == 0.0)
  //     THROW("Constructing a Line from two equal Points!");
  if (len > 0)
  {
      dir_.normalize();
      dir_ *= (len/(par2-par1));
  }
  location_ = point1 - par1*dir_;

    if (isReversed)
        reverseParameterDirection();
}

  // Copy constructor
//===========================================================================
Line& Line::operator=(const Line& other)
//===========================================================================
{
  if (&other == this)
    return *this;
  else
    {
      location_ = other.location_;
      dir_ = other.dir_;
      parbound1_ = other.parbound1_;
      parbound2_ = other.parbound2_;
      startparam_ = other.startparam_;
      endparam_ = other.endparam_;
      isReversed_ = other.isReversed_;
      return *this;
    }
}


// Destructor
//===========================================================================
Line::~Line()
//===========================================================================
{
}

//===========================================================================
void Line::read(std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    location_.resize(dim);
    dir_.resize(dim);
    is >> location_
       >> dir_;

    // Note: dir_ is not normalized.

    int isBounded; 
    is >> isBounded;
    bool has_param_int = (isBounded >= 10);
    isBounded = isBounded % 10;
    if (isBounded == 0) {
        // Unbounded
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
void Line::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << location_.dimension() << endl
       << location_ << endl
       << dir_ << endl;
    
    if (!isBounded()) {
        os << "0" << endl;
    }
    else {
        os << "11" << endl;
	os << parbound1_ << " " << parbound2_ << endl;
        os << startparam() << " " << endparam() << endl;
	
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
BoundingBox Line::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    vector<Point> points(2);
    point(points[0], startparam_);
    point(points[1], endparam_);
    box.setFromPoints(points);
    return box;
}


//===========================================================================
int Line::dimension() const
//===========================================================================
{
    return location_.dimension();
}


//===========================================================================
ClassType Line::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType Line::classType()
//===========================================================================
{
    return Class_Line;
}


//===========================================================================
Line* Line::clone() const
//===========================================================================
{
    Line* line = new Line(location_, dir_, isReversed_);
    line->setParamBounds(parbound1_, parbound2_);
    line->setParameterInterval(startparam_, endparam_);
    return line;
}


//===========================================================================
void Line::point(Point& pt, double tpar) const
//===========================================================================
{
    getReversedParameter(tpar);
    if (isBounded())
      tpar = parbound1_ + 
	(tpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_);
    pt = location_ + tpar * dir_;
}



//===========================================================================
void Line::point(vector<Point>& pts, 
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

    // The derivative is just the direction vector
    double fac = (isBounded()) ?
      (parbound2_-parbound1_)/(endparam_-startparam_) : 1.0;
    pts[1] = fac*dir_;

    if (isReversed())
        pts[1] *= -1.0;

    // Second order and higher derivatives vanish. They are already
    // set to zero, so we return.
    return;

}


//===========================================================================
double Line::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double Line::endparam() const
//===========================================================================
{
    return endparam_;
}


//===========================================================================
void Line::swapParameters2D()
//===========================================================================
{
    if (dimension() == 2) {
        swap(location_[0], location_[1]);
        swap(dir_[0], dir_[1]);
    }
}


//===========================================================================
  void Line::setParameterInterval(double t1, double t2)
//===========================================================================
  {
    startparam_ = t1;
    endparam_ = t2;
  }


//===========================================================================
SplineCurve* Line::geometryCurve()
//===========================================================================
{
    return createSplineCurve();
}


//===========================================================================
SplineCurve* Line::createSplineCurve() const
//===========================================================================
{
    double t0 = startparam();
    double t1 = endparam();
    double inf = numeric_limits<double>::infinity();
    double max = 1.0e8; // "Large" number...
    if (t0 == -inf)
        t0 = -max;
    if (t1 == inf)
        t1 = max;

    Point p0, p1;
    point(p0, t0);
    point(p1, t1);
    return new SplineCurve(p0, t0, p1, t1);
}


//===========================================================================
bool Line::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    // We consider a Line as degenerate if the length of the direction
    // vector dir_ is smaller than the epsilon.

    return (endparam_ - startparam_ < degenerate_epsilon) ||
        (dir_.length() < degenerate_epsilon);

}


//===========================================================================
Line* Line::subCurve(double from_par, double to_par,
                     double fuzzy) const
//===========================================================================
{
  Line *line = clone();
  bool bounded = isBounded();
  if (isReversed())
    {
      double start = endparam_ - (to_par - startparam_);
      double end = startparam_ + (endparam_ - from_par);

      if (start > end)
      	{
      	  std::swap(start, end);
      	}
      double bound1 = bounded ? parbound1_ + 
      	(start-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) :
      	start;
      double bound2 = bounded ? parbound1_ + 
      	(end-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) :
      	end;
      line->setParamBounds(bound1, bound2);
      if (from_par > to_par)
      	std::swap(from_par, to_par);
      line->setParameterInterval(from_par, to_par);
    }
  else
    {
      double bound1 = bounded ? parbound1_ + 
	(from_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) :
	from_par;
      double bound2 = bounded ? parbound1_ + 
	(to_par-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) :
	to_par;
      line->setParamBounds(bound1, bound2);
      line->setParameterInterval(from_par, to_par);
    }
  return line;
}


//===========================================================================
DirectionCone Line::directionCone() const
//===========================================================================
{
    Point dir = isReversed() ? -dir_ : dir_;
    return DirectionCone(dir);
}


//===========================================================================
void Line::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
  double dist;
  appendCurve(cv, 0, dist, reparam);
}


//===========================================================================
void Line::appendCurve(ParamCurve* cv,
                       int continuity, double& dist, bool reparam, double tol)
//===========================================================================
{
  // Check input
  if (cv->instanceType() != Class_Line)
    THROW("Inconsistency in curve types for appendCurve"); 

  double eps = 1.0e-5; // Not really used
  double angtol = 0.01;
  vector<Point> der1(2), der2(2);
  point(der1, endparam(), 1);
  cv->point(der2, cv->startparam(), 1);
  double ang = der1[1].angle(der2[1]);
  if (ang > angtol)
    THROW("AppendCurve: Line descriptions not compatible");

  dist = der1[0].dist(der2[0]);
  if (isReversed())
    setParamBounds(parbound1_-cv->length(eps), parbound2_);
  else
    setParamBounds(parbound1_, parbound2_+cv->length(eps));
}


//===========================================================================
void Line::closestPoint(const Point& pt,
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

    if (isBounded())
      {
	double fac = (parbound2_-parbound1_)/(endparam_-startparam_);
	tmin = parbound1_ + fac*(tmin - startparam_);
	tmax = parbound1_ + fac*(tmax - startparam_);
      }

    Point vec = pt - location_;
    Point dirnormal = dir_;
    double dirlen = dir_.length();
    if (dirlen > 0.0)
    {
        dirnormal.normalize();
        clo_t = vec * dirnormal / dirlen;
    }
    else
    {
        clo_t = (seed != nullptr) ? *seed : tmin;
    }
    if (clo_t < tmin)
        clo_t = tmin;
    if (clo_t > tmax)
        clo_t = tmax;
    clo_pt = location_ + clo_t * dir_;
    clo_dist = (clo_pt - pt).length();
    if (isBounded())
      {
	clo_t = startparam_ + 
	  (clo_t-parbound1_)*(endparam_-startparam_)/(parbound2_-parbound1_);
      }
    getReversedParameter(clo_t);
}


//===========================================================================
double Line::length(double tol)
//===========================================================================
{
    if (!isBounded())
        return numeric_limits<double>::infinity();

    double len = parbound2_ - parbound1_;
    return len * dir_.length();
}


//===========================================================================
void Line::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
        THROW("First parameter must be strictly less than second.");

    double start = isBounded() ? 
      parbound1_ + (startpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) : startpar;
    double end =  isBounded() ?
      parbound1_ + (endpar-startparam_)*(parbound2_-parbound1_)/(endparam_-startparam_) : endpar;
    parbound1_ = startpar;
    parbound2_ = endpar;
    startparam_ = start;
    endparam_ = end;
}

//===========================================================================
void Line::translateCurve(const Point& dir)
//===========================================================================
{
  location_ += dir;
}

//===========================================================================
bool Line::isBounded() const
//===========================================================================
{
    return parbound1_ > -numeric_limits<double>::infinity() &&
        parbound2_ < numeric_limits<double>::infinity();
}

//===========================================================================
  bool Line::isLinear(Point& dir, double tol)
//===========================================================================
{
  dir = dir_;
  return true;
}

//===========================================================================
  bool Line::isInPlane(const Point& loc, const Point& axis,
			   double eps, Point& normal) const
//===========================================================================
{
  normal = dir_.cross(axis);
  Point vec = location_ - loc;
  if (vec.length() < eps)
      return true;
  
  if (normal.length() < eps)
    {
      normal = vec.cross(axis);
      return true;
    }

  if (normal.length() > 0.0)
  {
      normal.normalize();
  }
  double dist = (location_ + dir_)*normal;
  return (dist < eps);
}

//===========================================================================
bool Line::isInPlane(const Point& norm,
		       double eps, Point& pos) const

//===========================================================================
{
  double ang = norm.angle(dir_);
  pos = location_;

  return (fabs(0.5*M_PI-ang) <= eps);
}
} // namespace Go
