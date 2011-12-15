//==========================================================================
//                                                                          
// File: Cone.C                                                              
//                                                                          
// Created: Tue Nov 18 16:15:45 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Cone.C,v 1.12 2009-05-13 07:30:52 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>
#include <limits>


using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::shared_ptr;


namespace Go
{


// Constructor.
//===========================================================================
Cone::Cone(double radius,
	   Point location, Point z_axis, Point x_axis,
	   double cone_angle)
    : radius_(radius),
      location_(location), z_axis_(z_axis), x_axis_(x_axis),
      cone_angle_(cone_angle)
//===========================================================================
{
    // The use of the z- and x-axis in this constructor - and in the
    // definition of a cone - is based on the use of the
    // 'axis2_placement_3d' entity in STEP.

    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }
    setCoordinateAxes();

    Array<double, 2> ll(0.0, -numeric_limits<double>::infinity());
    Array<double, 2> ur(2.0 * M_PI, numeric_limits<double>::infinity());
    domain_ = RectDomain(ll, ur);
}


// Destructor
//===========================================================================
Cone::~Cone()
//===========================================================================
{
}

//===========================================================================
void Cone::read (std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    location_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> radius_
       >> location_
       >> z_axis_
       >> x_axis_
       >> cone_angle_;

    setCoordinateAxes();
}


//===========================================================================
void Cone::write(std::ostream& os) const
//===========================================================================
{
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl
       << cone_angle_ << endl;   
}

//===========================================================================
int Cone::dimension() const
//===========================================================================
{
    return location_.dimension();
}

//===========================================================================
ClassType Cone::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
BoundingBox Cone::boundingBox() const
//===========================================================================
{
    // A rather unefficient hack...
    Cone* cone = const_cast<Cone*>(this);
    SplineSurface* tmp = cone->geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
const RectDomain& Cone::parameterDomain() const
//===========================================================================
{
    return domain_;
}


//===========================================================================
CurveLoop Cone::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty loop.");
    CurveLoop loop;
    return loop;
}


//===========================================================================
std::vector<CurveLoop> 
Cone::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Cone::normalCone() const
//===========================================================================
{
    double umin = domain_.umin();
    double umax = domain_.umax();
    Point dir;
    normal(dir, 0.5*(umin+umax), 0.0);

    return DirectionCone(dir, 0.5*(umax-umin));
}


//===========================================================================
DirectionCone Cone::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    if (pardir_is_u) {
	DirectionCone normals = normalCone();
	Point dir = normals.centre();
	Point tandir = z_axis_.cross(dir);
	double angle = normals.angle();
	return DirectionCone(tandir, angle);
    }
    else {
	return DirectionCone(z_axis_);
    }
}


//===========================================================================
void Cone::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    pt = location_
	+ (radius_ + vpar * tan(cone_angle_)) * (cos(upar) * x_axis_ 
						 + sin(upar) * y_axis_)
	+ vpar * z_axis_;
}


//===========================================================================
void Cone::point(std::vector<Point>& pts, 
		 double upar, double vpar,
		 int derivs,
		 bool u_from_right,
		 bool v_from_right,
		 double resolution) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0,
		   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)/2;
    int ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz< totpts,
		   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
        }
	pts[i].setValue(0.0);
    }

    // Zero'th derivative
    point(pts[0], upar, vpar);
    if (derivs == 0)
        return;

    // First derivatives
    pts[1] = (radius_ + vpar * tan(cone_angle_)) 
	* (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[2] = tan(cone_angle_) * (cos(upar) * x_axis_ 
				 + sin(upar) * y_axis_)  
	+ z_axis_;
    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void Cone::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    double tana = tan(cone_angle_);
    double tana2 = tana * tana;
    n = (cos(upar) * x_axis_ + sin(upar) * y_axis_ - tana * z_axis_)
	/ sqrt(1.0 + tana2);

    double k = radius_ + vpar * tana;
    if (k > 0.0) {
	// Leave n as it is
	return;
    }
    else if (k < 0.0) {
	// Change sign of n
	n *= -1.0;
	return;
    }
    // If we get here, then k == 0.0, and the normal is not defined
    MESSAGE("Normal is not defined.");
    n = Point(0.0, 0.0, 0.0);

}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Cone::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
Cone* Cone::subSurface(double from_upar, double from_vpar,
			       double to_upar, double to_vpar,
			       double fuzzy) const
//===========================================================================
{
    Cone* cone = clone();
    cone->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return cone;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
Cone::subSurfaces(double from_upar, double from_vpar,
		  double to_upar, double to_vpar,
		  double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Cone> cone(subSurface(from_upar, from_vpar,
				     to_upar, to_vpar));
    res.push_back(cone);
    return res;
}


//===========================================================================
double 
Cone::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    MESSAGE("Does not make sense. Return arbitrarily zero.");
    return 0.0;
}


//===========================================================================
void Cone::closestPoint(const Point& pt,
			double&        clo_u,
			double&        clo_v, 
			Point&         clo_pt,
			double&        clo_dist,
			double         epsilon,
			const RectDomain* domain_of_interest,
			double   *seed) const
//===========================================================================
{
    // Set the domain
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    if (domain_of_interest != NULL) {
      umin = std::max(umin, domain_of_interest->umin());
      umax = std::min(umax, domain_of_interest->umax());
      vmin = std::max(vmin, domain_of_interest->vmin());
      vmax = std::min(vmax, domain_of_interest->vmax());
    }

    // Identify the two values of the v-parameter where an unbounded
    // cone is orthogonal to the cone-to-point vector.
    double rad = radius_;
    if (radius_ < epsilon)
	rad = 1.0;
    Point loc = location_;
    Circle circle(rad, loc, z_axis_, x_axis_);
    circle.closestPoint(pt, 0.0, 2.0*M_PI, clo_u, clo_pt, clo_dist);
    shared_ptr<Line> line = getLine(clo_u);
    double vvalmin, vvalmax, tmp;
    line->closestPoint(pt, vmin, vmax, vvalmin, clo_pt, clo_dist);
    line = getLine(clo_u - M_PI);
    line->closestPoint(pt, vmin, vmax, vvalmax, clo_pt, clo_dist);
    if (vvalmin > vvalmax) {
	tmp = vvalmin;
	vvalmin = vvalmax;
	vvalmax = tmp;
    }

    // If vmax < vvalmin or vmin > vvalmax we are done
    if (vmax < vvalmin) {
	rad = radius_ + vmax * tan(cone_angle_);
	loc = location_ + vmax * z_axis_;
	circle = Circle(rad, loc, z_axis_, x_axis_);
	circle.setParamBounds(umin, umax);
	circle.closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
	clo_v = vmax;
	return;
    }
    if (vmin > vvalmax) {
	rad = radius_ + vmin * tan(cone_angle_);
	loc = location_ + vmin * z_axis_;
	circle = Circle(rad, loc, z_axis_, x_axis_);
	circle.setParamBounds(umin, umax);
	circle.closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
	clo_v = vmin;
	return;
    }

    // If we get here, then we effectively have a bounded patch on the
    // cone and the closest point lies on one of the edge curves!

    // If the vvalmin/max bounds are in the domain, we adopt these for
    // our bounds.
    vmin = (vmin > vvalmin ? vmin : vvalmin);
    vmax = (vmax < vvalmax ? vmax : vvalmax);

    // Examine the edge curves...

    double tmp_clo_u, tmp_clo_v, tmp_clo_dist;
    Point tmp_clo_pt;

    // Bottom - a circle
    rad = radius_ + vmin * tan(cone_angle_);
    loc = location_ + vmin * z_axis_;
    circle = Circle(rad, loc, z_axis_, x_axis_);
    circle.setParamBounds(umin, umax);
    circle.closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
    clo_v = vmin;
    if (clo_dist < epsilon)
	return;

    // Top - a circle
    rad = radius_ + vmax * tan(cone_angle_);
    loc = location_ + vmax * z_axis_;
    circle = Circle(rad, loc, z_axis_, x_axis_);
    circle.setParamBounds(umin, umax);
    circle.closestPoint(pt, umin, umax, tmp_clo_u, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_v = vmax;
    if (tmp_clo_dist < clo_dist) {
	clo_u = tmp_clo_u;
	clo_v = tmp_clo_v;
	clo_pt = tmp_clo_pt;
	clo_dist = tmp_clo_dist;
	if (clo_dist < epsilon)
	    return;
    }

    // Are there more edges?
    if (fabs(umax - umin - 2.0 * M_PI) < epsilon)
	return;

    // Left - a line
    line = getLine(umin);
    line->closestPoint(pt, vmin, vmax, tmp_clo_v, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_u = umin;
    if (tmp_clo_dist < clo_dist) {
	clo_u = tmp_clo_u;
	clo_v = tmp_clo_v;
	clo_pt = tmp_clo_pt;
	clo_dist = tmp_clo_dist;
	if (clo_dist < epsilon)
	    return;
    }

    // Right - a line
    line = getLine(umax);
    line->closestPoint(pt, vmin, vmax, tmp_clo_v, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_u = umax;
    if (tmp_clo_dist < clo_dist) {
	clo_u = tmp_clo_u;
	clo_v = tmp_clo_v;
	clo_pt = tmp_clo_pt;
	clo_dist = tmp_clo_dist;
    }

    return;
}


//===========================================================================
void Cone::closestBoundaryPoint(const Point& pt,
				double&        clo_u,
				double&        clo_v, 
				Point&       clo_pt,
				double&        clo_dist,
				double epsilon,
				const RectDomain* rd,
				double *seed) const
//===========================================================================
{
    // This is a bit like cheating...

    SplineSurface* sf = geometrySurface();
    sf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
			     rd, seed);
    delete sf;
}


//===========================================================================
void Cone::getBoundaryInfo(Point& pt1, Point& pt2,
			    double epsilon, SplineCurve*& cv,
			    SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
void Cone::turnOrientation()
//===========================================================================
{
    swapParameterDirection();
}



//===========================================================================
void Cone::swapParameterDirection()
//===========================================================================
{
    MESSAGE("swapParameterDirection() not implemented.");
}


//===========================================================================
void Cone::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    MESSAGE("reverseParameterDirection() not implemented.");
}


//===========================================================================
bool Cone::isDegenerate(bool& b, bool& r,
			bool& t, bool& l, double tolerance) const
//===========================================================================
{
    bool res = false;
    b = false;
    r = false;
    t = false;
    l = false;
    double kmin = radius_ + parameterDomain().vmin() * tan(cone_angle_);
    if (kmin == 0.0) {
	b = true;
	res = true;
    }
    double kmax = radius_ + parameterDomain().vmax() * tan(cone_angle_);
    if (kmax == 0.0) {
	t = true;
	res = true;
    }
    return res;
}


//===========================================================================
void Cone::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Not this kind of degeneracy
    return;
}

//===========================================================================
void Cone::setCoordinateAxes()
//===========================================================================
{
    // The x-, y- and z-axes defines a right-handed coordinate system.

    z_axis_.normalize();
    Point tmp = x_axis_ - (x_axis_ * z_axis_) * z_axis_;
    if (tmp.length() == 0.0)
	THROW("X-axis parallel to Z-axis.");

    x_axis_ = tmp;
    y_axis_ = z_axis_.cross(x_axis_);
    x_axis_.normalize();
    y_axis_.normalize();

}


//===========================================================================
void Cone::setParameterBounds(double from_upar, double from_vpar,
			      double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
	THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
	THROW("First v-parameter must be strictly less than second.");
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
	THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
	THROW("(to_upar - from_upar) must not exceed 2pi.");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
bool Cone::isBounded() const
//===========================================================================
{
    // It is enough to check the v-direction, since the u-direction is
    // always bounded.

    return domain_.vmin() > -numeric_limits<double>::infinity() &&
	domain_.vmax() < numeric_limits<double>::infinity();

}


//===========================================================================
SplineSurface* Cone::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Cone::createSplineSurface() const
//===========================================================================
{
    // Based on SISL functions s1021 and s1022.

    // Not properly tested... @jbt

    // First handle the case if not bounded
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    if (!isBounded()) {
	double max = 1.0e8; // "Large" number...
	if (vmin == -numeric_limits<double>::infinity())
	    vmin = -max;
	if (vmax == numeric_limits<double>::infinity())
	    vmax = max;
    }

    // Knot vector around the cone.
    double et1[12]; 
    for (int ki = 0; ki < 12; ki++) {
	if (ki == 0 || ki == 1 || ki == 2)
	    et1[ki] = (double)0.0;
	else if (ki == 3 || ki == 4)
	    et1[ki] = 0.5 * M_PI;
	else if (ki == 5 || ki == 6)
	    et1[ki] = M_PI;
	else if (ki == 7 || ki == 8)
	    et1[ki] = 1.5 * M_PI;
	else if (ki == 9 || ki == 10 || ki == 11)
	    et1[ki] = 2.0 * M_PI;
    }
    // Knot vector along the cone.
    double et2[4]; 
    for (int ki = 0; ki < 4; ki++) {
	if (ki == 0 || ki == 1)
	    et2[ki] = vmin;
	else if (ki == 2 || ki == 3)
	    et2[ki] = vmax;
    }


    // Coefficients
    double rcoef[72];
    double weight = 1.0/sqrt(2.0);
    Point bottom_pos = location_ + vmin * z_axis_;
    Point top_pos = location_ + vmax * z_axis_;
    double brad = radius_ + vmin * tan(cone_angle_);
    double trad = radius_ + vmax * tan(cone_angle_);
    Point b1_axis = brad * x_axis_;
    Point b2_axis = brad * y_axis_;
    Point t1_axis = trad * x_axis_;
    Point t2_axis = trad * y_axis_;
    for (int ki = 0; ki < 3; ki++){
	rcoef[ki] = bottom_pos[ki] + b1_axis[ki];
	rcoef[4 + ki] = weight*(bottom_pos[ki] + b1_axis[ki] + b2_axis[ki]);
	rcoef[8 + ki] = bottom_pos[ki] + b2_axis[ki];
	rcoef[12 + ki] = weight*(bottom_pos[ki] - b1_axis[ki] + b2_axis[ki]);
	rcoef[16 + ki] = bottom_pos[ki] - b1_axis[ki];
	rcoef[20 + ki] = weight*(bottom_pos[ki] - b1_axis[ki] - b2_axis[ki]);
	rcoef[24 + ki] = bottom_pos[ki] - b2_axis[ki];
	rcoef[28 + ki] = weight*(bottom_pos[ki] + b1_axis[ki] - b2_axis[ki]);
	rcoef[32 + ki] = rcoef[ki];

	rcoef[36 + ki] = top_pos[ki] + t1_axis[ki];
	rcoef[40 + ki] = weight*(top_pos[ki] + t1_axis[ki] + t2_axis[ki]);
	rcoef[44 + ki] = top_pos[ki] + t2_axis[ki];
	rcoef[48 + ki] = weight*(top_pos[ki] - t1_axis[ki] + t2_axis[ki]);
	rcoef[52 + ki] = top_pos[ki] - t1_axis[ki];
	rcoef[56 + ki] = weight*(top_pos[ki] - t1_axis[ki] - t2_axis[ki]);
	rcoef[60 + ki] = top_pos[ki] - t2_axis[ki];
	rcoef[64 + ki] = weight*(top_pos[ki] + t1_axis[ki] - t2_axis[ki]);
	rcoef[68 + ki] = rcoef[36 + ki];
    }
    // The rational weights
    for (int ki = 3; ki < 72; ki += 4) {
	if (ki == 3 || ki == 11 || ki == 19 || ki == 27 || ki == 35
	    || ki == 39 || ki == 47 || ki == 55 || ki == 63 || ki == 71)
	    rcoef[ki] = 1.0;
	else
	    rcoef[ki] = weight;
    }

    int ncoefs1 = 9;
    int ncoefs2 = 2;
    int order1 = 3;
    int order2 = 2;
    int dim = 3;
    bool rational = true;
    SplineSurface surface(ncoefs1, ncoefs2, order1, order2,
			  et1, et2, rcoef, dim, rational);

    // Extract subpatch. We need all this because 'surface' is not
    // arc-length parametrized in the u-direction. We also need to
    // avoid the singularity of the cone...
    double umin = domain_.umin();
    double umax = domain_.umax();
    Point pt, tmppt;
    point(pt, umax - umin, 0.0);
    double tmpu, tmpdist;
    Circle circle(radius_, location_, z_axis_, x_axis_);
    SplineCurve* scircle = circle.geometryCurve();
    scircle->closestPoint(pt, 0.0, 2.0 * M_PI, tmpu, tmppt, tmpdist);
    double epsilon = 1.0e-10;
    if (tmpu < epsilon && umax - umin == 2.0 * M_PI) {
	tmpu = 2.0 * M_PI;
    }
    SplineSurface* subpatch = surface.subSurface(0.0, vmin, tmpu, vmax);
    subpatch->basis_u().rescale(umin, umax);
    translateSplineSurf(-location_, *subpatch);
    rotateSplineSurf(z_axis_, umin, *subpatch);
    translateSplineSurf(location_, *subpatch);
    delete scircle;

    return subpatch;
}


//===========================================================================
shared_ptr<Line> Cone::getLine(double upar) const 
//===========================================================================
{
    Point cossin = cos(upar) * x_axis_ + sin(upar) * y_axis_;
    Point loc = location_ + radius_ * cossin;
    Point dir = tan(cone_angle_) * cossin + z_axis_;
    shared_ptr<Line> line(new Line(loc, dir));
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    line->setParamBounds(vmin, vmax);
    return line;
}


//===========================================================================


} // namespace Go
