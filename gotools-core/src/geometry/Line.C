//==========================================================================
//                                                                          
// File: Line.C                                                              
//                                                                          
// Created: Mon Sep  8 14:50:58 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Line.C,v 1.9 2009-02-19 15:52:54 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================



#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>
#include <limits>


using std::vector;
using std::endl;
using std::numeric_limits;
using std::streamsize;


namespace Go {


// Constructor. Input is point and direction
//===========================================================================
Line::Line(Point point, Point direction)
    : location_(point), dir_(direction),
    startparam_(-numeric_limits<double>::infinity()),
    endparam_(numeric_limits<double>::infinity())
//===========================================================================
{
    // Note: dir_ is not normalized.
}


// Constructor to bounded line. Input is point, direction
// and length
//===========================================================================
Line::Line(Point point, Point direction, double length)
  : location_(point), dir_(direction), startparam_(0.0),
    endparam_(length)
//===========================================================================
{
  dir_.normalize();
}

// Constructor to bounded line. Input is point, direction
// and length
//===========================================================================
Line::Line(Point point1, Point point2, double par1, double par2)
  : startparam_(par1), endparam_(par2)
//===========================================================================
{
  dir_ = point2 - point1;
  double len = dir_.length();
  dir_.normalize();
  dir_ *= (len/(par2-par1));
  location_ = point1 - par1*dir_;
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
    if (isBounded == 0) {
        // Unbounded
        startparam_ = -numeric_limits<double>::infinity();
        endparam_ = numeric_limits<double>::infinity();
    }
    else if (isBounded == 1) {
        is >> startparam_
            >> endparam_;
    }
    else {
        THROW("Bounded flag must be 0 or 1");
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
        os << "1" << endl;
        os << startparam() << " " << endparam() << endl;
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
    Line* line = new Line(location_, dir_);
    line->setParamBounds(startparam_, endparam_);
    return line;
}


//===========================================================================
void Line::point(Point& pt, double tpar) const
//===========================================================================
{
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
    pts[1] = dir_;

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
void Line::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    // This function can be implemented in two different ways. One
    // where we simply flip the direction vector, dir -> -dir_, and
    // one where we preserve the parameter interval but switch the
    // roles of the endpoints. The former is more natural to Lines due
    // to their natural parametrization. However, we will use the
    // latter since it is more likely to be the intention of the
    // caller. If direction flip is intended, simply construct a new
    // Line with a negative direction vector.

    dir_ = -dir_;
    if (isBounded()) {
        double x = endparam_ + startparam_;
        location_ -= x * dir_;
    }

    if (dimension() == 2 && switchparam) {
        double tmp = location_[0];
        location_[0] = location_[1];
        location_[1] = tmp;
        tmp = dir_[0];
        dir_[0] = dir_[1];
        dir_[1] = tmp;
    }
}


//===========================================================================
  void Line::setParameterInterval(double t1, double t2)
//===========================================================================
  {
    if (isBounded())
      {
	Point p1 = this->ParamCurve::point(t1);
	Point p2 = this->ParamCurve::point(t2);
	double len = p1.dist(p2);
	dir_.normalize();
	dir_ *= (len/(endparam_-startparam_));
	location_ -= (t1-startparam_)*dir_;
	startparam_ = t1;
	endparam_ = t2;
      }
    else
      {
	// VSK. This really dosn't make sense
	setParamBounds(t1, t2);
      }
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
    double max = 1.0e8; // "Large" number...
    if (t0 == -numeric_limits<double>::infinity())
        t0 = -max;
    if (t1 == numeric_limits<double>::infinity())
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
    Line* line = clone();
    line->setParamBounds(from_par, to_par);
    return line;
}


//===========================================================================
DirectionCone Line::directionCone() const
//===========================================================================
{
    return DirectionCone(dir_);
}


//===========================================================================
void Line::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void Line::appendCurve(ParamCurve* cv,
                       int continuity, double& dist, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
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

    Point vec = pt - location_;
    Point dirnormal = dir_;
    dirnormal.normalize();
    double dirlen = dir_.length();
    clo_t = vec * dirnormal / dirlen;
    if (clo_t < tmin)
        clo_t = tmin;
    if (clo_t > tmax)
        clo_t = tmax;
    clo_pt = location_ + clo_t * dir_;
    clo_dist = (clo_pt - pt).length();
}


//===========================================================================
double Line::length(double tol)
//===========================================================================
{
    if (!isBounded())
        return numeric_limits<double>::infinity();

    double len = endparam_ - startparam_;
    return len * dir_.length();
}


//===========================================================================
void Line::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
        THROW("First parameter must be strictly less than second.");

    startparam_ = startpar;
    endparam_ = endpar;
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
    return startparam_ > -numeric_limits<double>::infinity() &&
        endparam_ < numeric_limits<double>::infinity();
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

  normal.normalize();
  double dist = (location_ + dir_)*normal;
  return (dist < eps);
}

//===========================================================================


} // namespace Go
