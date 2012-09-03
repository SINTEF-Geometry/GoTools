//===========================================================================
//                                                                           
// File: Hyperbola.C                                                         
//                                                                           
// Created: Thu Sep 17 12:49:17 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/Hyperbola.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>
#include <limits>


using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;


namespace Go
{


//===========================================================================
Hyperbola::Hyperbola(Point location, Point direction, Point normal,
		     double r1, double r2)
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
    THROW("read(): Not yet implemented!");
}


//===========================================================================
void Hyperbola::write(std::ostream& os) const
//===========================================================================
{
    THROW("write(): Not yet implemented!");
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
    return new Hyperbola(location_, vec1_, normal_, r1_, r2_);
}


//===========================================================================
void Hyperbola::point(Point& pt, double tpar) const
//===========================================================================
{
    ASSERT((tpar >= startparam_) && ( tpar < endparam_));

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
    double sinh_t = sinh(tpar);
    double cosh_t = cosh(tpar);
    for (int ki = 1; ki < derivs + 1; ++ki) {
	pts[ki] = (ki%2 == 1) ? r1_*sinh_t*vec1_ + r2_*cosh_t*vec2_ :
	    r1_*cosh_t*vec1_ + r2_*sinh_t*vec2_;
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
void Hyperbola::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    if (switchparam) {
	if (dimension() == 2) {
	    Point tmp = vec1_;
	    vec1_ = vec2_;
	    vec2_ = tmp;
	}
	return;
    }

    // Flip
    normal_ = -normal_;
    vec2_ = -vec2_;

    // Rotate to keep parametrization consistent
    double alpha = startparam_ + endparam_;
    if (alpha >= 2.0 * M_PI)
	alpha -= 2.0 * M_PI;
    if (alpha <= -2.0 * M_PI)
	alpha += 2.0 * M_PI;
    if (alpha != 0.0) {
	GeometryTools::rotatePoint(normal_, -alpha, vec1_);
	GeometryTools::rotatePoint(normal_, -alpha, vec2_);
    }
}


//===========================================================================
void Hyperbola::setParameterInterval(double t1, double t2)
//===========================================================================
{
    setParamBounds(t1, t2);
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
    MESSAGE("Not yet implemented.");
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
    if (from_par >= to_par)
	THROW("First parameter must be strictly less than second.");

    Hyperbola* hyperbola = clone();
    hyperbola->setParamBounds(from_par, to_par);
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
    MESSAGE("Not implemented!");
}


//===========================================================================
void Hyperbola::appendCurve(ParamCurve* cv,
			  int continuity, double& dist, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
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
    if (startpar < -2.0 * M_PI || endpar > 2.0 * M_PI)
	THROW("Parameters must be in [-2pi, 2pi].");
    if (endpar - startpar > 2.0 * M_PI)
	THROW("(endpar - startpar) must not exceed 2pi.");

    startparam_ = startpar;
    endparam_ = endpar;
}


//===========================================================================
void Hyperbola::translateCurve(const Point& dir)
//===========================================================================
{
  location_ += dir;
}

//===========================================================================
bool Hyperbola::isBounded()
//===========================================================================
{
    return startparam_ > -numeric_limits<double>::infinity() &&
	endparam_ < numeric_limits<double>::infinity();
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
