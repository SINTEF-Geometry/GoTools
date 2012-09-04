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
using std::streamsize;
using std::swap;


namespace Go
{


//===========================================================================
Parabola::Parabola(Point location, Point direction,
		   Point normal, double focal_dist,
                   bool isReversed_)
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

    double inf = numeric_limits<double>::infinity();
    setParamBounds(-inf, inf);

    if (isReversed())
        reverseParameterDirection();
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
    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    location_.resize(dim);
    normal_.resize(dim);
    vec1_.resize(dim);
    is >> f_
       >> location_
       >> normal_
       >> vec1_;

    if(dim == 3)
        normal_.normalize();
    setSpanningVectors();

    int isBounded; 
    is >> isBounded;
    if (isBounded == 0) {
        // Unbounded - don't read parameters
        double inf = numeric_limits<double>::infinity();
        setParamBounds(-inf, inf);
    }
    else if (isBounded == 1) {
        is >> startparam_ >> endparam_;
    }
    else {
        THROW("Bounded flag must be 0 or 1");
    }

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
void Parabola::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    int dim = dimension();
    os << dim << endl
       << f_ << endl
       << location_ << endl
       << normal_ << endl
       << vec1_ << endl;

    if (!isBounded()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl
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
    Parabola* par = new Parabola(location_, vec1_, normal_, f_,
        isReversed_);
    par->setParamBounds(startparam_, endparam_);
    return par;
}


//===========================================================================
void Parabola::point(Point& pt, double tpar) const
//===========================================================================
{
    getReversedParameter(tpar);
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
    getReversedParameter(tpar);
    if (derivs <= 1) {
        pts[1] = 2.0*f_*(tpar*vec1_ + vec2_);
        if (isReversed()) {
            pts[1] *= -1.0;
        }
    }
    if (derivs <= 2)
        pts[2] = 2.0*f_*vec1_;
    return;
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
void Parabola::swapParameters2D()
//===========================================================================
{
    if (dimension() == 2) {
        swap(location_[0], location_[1]);
        swap(vec1_[0], vec1_[1]);
        swap(vec2_[0], vec2_[1]);
    }
}


//===========================================================================
void Parabola::setParameterInterval(double t1, double t2)
//===========================================================================
{
    MESSAGE("setParameterInterval() doesn't make sense.");
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
    MESSAGE("createSplineCurve() not yet implemented.");
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
    MESSAGE("appendCurve() not implemented!");
}


//===========================================================================
void Parabola::appendCurve(ParamCurve* cv,
			  int continuity, double& dist, bool reparam)
//===========================================================================
{
    MESSAGE("appendCurve() not implemented!");
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
void Parabola::setParamBounds(double startpar, double endpar)
//===========================================================================
{
    if (startpar >= endpar)
	THROW("First parameter must be strictly less than second.");

    startparam_ = startpar;
    endparam_ = endpar;
}


//===========================================================================
bool Parabola::isBounded() const
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

//===========================================================================
bool Parabola::isInPlane(const Point& norm,
		       double eps, Point& pos) const

//===========================================================================
{
  double ang = norm.angle(normal_);
  pos = location_;

  return (ang <= eps || fabs(M_PI-ang) <= eps);
}


}
