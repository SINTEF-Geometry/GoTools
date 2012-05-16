//===========================================================================
//                                                                           
// File: Parabola.C                                                          
//                                                                           
// Created: Thu Sep 17 12:49:22 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/Parabola.h"
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
Parabola::Parabola(Point location, Point direction,
		   Point normal, double focal_dist)
    : location_(location), vec1_(direction), normal_(normal),
      f_(focal_dist)
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
Parabola::~Parabola()
//===========================================================================
{
}


//===========================================================================
void Parabola::read(std::istream& is)
//===========================================================================
{
    THROW("read(): Not yet implemented!");
}


//===========================================================================
void Parabola::write(std::ostream& os) const
//===========================================================================
{
    THROW("write(): Not yet implemented!");
}


//===========================================================================
BoundingBox Parabola::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    // A rather inefficient hack...
    Parabola* parabola = const_cast<Parabola*>(this);
    SplineCurve* tmp = parabola->geometryCurve();
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
int Parabola::dimension() const
//===========================================================================
{
    return location_.dimension();
}
    

//===========================================================================
ClassType Parabola::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType Parabola::classType()
//===========================================================================
{
    return Class_Parabola;
}


//===========================================================================
Parabola* Parabola::clone() const
//===========================================================================
{
    return new Parabola(location_, vec1_, normal_, f_);
}


//===========================================================================
void Parabola::point(Point& pt, double tpar) const
//===========================================================================
{
    ASSERT((tpar >= startparam_) && ( tpar < endparam_));

    pt = location_ + f_*(tpar*tpar*vec1_ + 2*tpar*vec2_);
}


//===========================================================================
void Parabola::point(std::vector<Point>& pts, 
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

    // Since the parabola is parametrized as:
    // c(t) = location_ + focal-dist_*(t^2*vec1_ + 2*t*vec2_),
    // the derivatives follow easily.
    for (int ki = 1; ki < derivs + 1; ++ki) {
	if (ki == 1)
	    pts[ki] = 2.0*f_*(tpar*vec1_ + vec2_);
	else if (ki == 2)
	    pts[ki] = 2.0*f_*vec1_;
	else
	{
	    pts[ki] = Point(dim);
	    pts[ki].setValue(0.0);
	}
    }
}


//===========================================================================
double Parabola::startparam() const
//===========================================================================
{
    return startparam_;
}


//===========================================================================
double Parabola::endparam() const
//===========================================================================
{
    return endparam_;
}


//===========================================================================
void Parabola::reverseParameterDirection(bool switchparam)
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
void Parabola::setParameterInterval(double t1, double t2)
//===========================================================================
{
    setParamBounds(t1, t2);
}


//===========================================================================
SplineCurve* Parabola::geometryCurve()
//===========================================================================
{
    return createSplineCurve();
}

//===========================================================================
SplineCurve* Parabola::createSplineCurve() const
//===========================================================================
{
    MESSAGE("Not yet implemented.");
    return NULL;
}


//===========================================================================
bool Parabola::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    // We consider a Parabola as degenerate if either radii is smaller
    // than the epsilon. Meaning that the curve is either flat or
    // thin.
    return ((f_*vec1_.length()*vec1_.length() < degenerate_epsilon) ||
	    (2.0*f_*vec2_.length() < degenerate_epsilon));
}



//===========================================================================
Parabola* Parabola::subCurve(double from_par, double to_par,
			  double fuzzy) const
//===========================================================================
{
    if (from_par >= to_par)
	THROW("First parameter must be strictly less than second.");

    Parabola* parabola = clone();
    parabola->setParamBounds(from_par, to_par);
    return parabola;
}

 
//===========================================================================
DirectionCone Parabola::directionCone() const
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
void Parabola::appendCurve(ParamCurve* cv, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void Parabola::appendCurve(ParamCurve* cv,
			  int continuity, double& dist, bool reparam)
//===========================================================================
{
    MESSAGE("Not implemented!");
}


//===========================================================================
void Parabola::closestPoint(const Point& pt,
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
double Parabola::length(double tol)
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
void Parabola::setParamBounds(double startpar, double endpar)
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
bool Parabola::isBounded()
//===========================================================================
{
    return startparam_ > -numeric_limits<double>::infinity() &&
	endparam_ < numeric_limits<double>::infinity();
}


//===========================================================================
void Parabola::translateCurve(const Point& dir)
//===========================================================================
{
  location_ += dir;
}

//===========================================================================
void Parabola::setSpanningVectors()
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


}
