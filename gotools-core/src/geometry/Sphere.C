//==========================================================================
//                                                                          
// File: Sphere.C                                                            
//                                                                          
// Created: Tue Nov 18 16:14:22 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Sphere.C,v 1.6 2009-03-06 16:10:08 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>
#include <limits>


using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::streamsize;


namespace Go
{


// Constructor.
//===========================================================================
Sphere::Sphere(double radius,
	       Point location, Point z_axis, Point x_axis)
    : radius_(radius),
      location_(location), z_axis_(z_axis), x_axis_(x_axis)
//===========================================================================
{
    // The use of the z- and x-axis in this constructor - and in the
    // definition of a sphere - is based on the use of the
    // 'axis2_placement_3d' entity in STEP.

    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }
    setCoordinateAxes();

    Array<double, 2> ll(0.0, -0.5 * M_PI);
    Array<double, 2> ur(2.0 * M_PI, 0.5 * M_PI);
    domain_ = RectDomain(ll, ur);
}


// Destructor
//===========================================================================
Sphere::~Sphere()
//===========================================================================
{
}

//===========================================================================
void Sphere::read (std::istream& is)
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
       >> x_axis_;

    setCoordinateAxes();
    Array<double, 2> ll(0.0, -0.5 * M_PI);
    Array<double, 2> ur(2.0 * M_PI, 0.5 * M_PI);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Sphere::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;
    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
int Sphere::dimension() const
//===========================================================================
{
    return location_.dimension();
}

//===========================================================================
ClassType Sphere::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
BoundingBox Sphere::boundingBox() const
//===========================================================================
{
    // A rather unefficient hack...
    Sphere* sph = const_cast<Sphere*>(this);
    SplineSurface* tmp = sph->geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
const RectDomain& Sphere::parameterDomain() const
//===========================================================================
{
    return domain_;
}


//===========================================================================
std::vector<CurveLoop> 
Sphere::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Sphere::normalCone() const
//===========================================================================
{
    // A rather unefficient hack...
    Sphere* sph = const_cast<Sphere*>(this);
    SplineSurface* tmp = sph->geometrySurface();
    DirectionCone dc = tmp->normalCone();
    delete tmp;
    return dc;
}


//===========================================================================
DirectionCone Sphere::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    // A rather unefficient hack...
    Sphere* sph = const_cast<Sphere*>(this);
    SplineSurface* tmp = sph->geometrySurface();
    DirectionCone dc = tmp->tangentCone(pardir_is_u);
    delete tmp;
    return dc;
}


//===========================================================================
void Sphere::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    pt = location_ 
	+ radius_ * (cos(vpar) * (cos(upar) * x_axis_ + sin(upar) * y_axis_)
		     + sin(vpar) * z_axis_);
}


//===========================================================================
void Sphere::point(std::vector<Point>& pts, 
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
    pts[1] = radius_ * cos(vpar) * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[2] = radius_ * (-sin(vpar) * (cos(upar) * x_axis_ 
				      + sin(upar) * y_axis_) 
			+ cos(vpar) * z_axis_);
    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void Sphere::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    n = cos(vpar) * (cos(upar) * x_axis_ + sin(upar) * y_axis_)
	+ sin(vpar) * z_axis_;;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Sphere::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
Sphere* Sphere::subSurface(double from_upar, double from_vpar,
			   double to_upar, double to_vpar,
			   double fuzzy) const
//===========================================================================
{
    Sphere* sphere = clone();
    sphere->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return sphere;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
Sphere::subSurfaces(double from_upar, double from_vpar,
		      double to_upar, double to_vpar,
		      double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Sphere> sphere(subSurface(from_upar, from_vpar,
					 to_upar, to_vpar));
    res.push_back(sphere);
    return res;
}


//===========================================================================
double 
Sphere::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    MESSAGE("Does not make sense. Return arbitrarily zero.");
    return 0.0;
}


//===========================================================================
void Sphere::closestPoint(const Point& pt,
			  double&        clo_u,
			  double&        clo_v, 
			  Point&         clo_pt,
			  double&        clo_dist,
			  double         epsilon,
			  const RectDomain* domain_of_interest,
			  double   *seed) const
//===========================================================================
{
    // Find relevant domain of interest
    RectDomain curr_domain_of_interest = domain_;
    if (domain_of_interest != NULL) {
	curr_domain_of_interest.intersectWith(*domain_of_interest);
    }

    // Algorithm:
    // 1) Find closest point on latitudinal circle (u-direction), including
    //    the bounds on u given by the domain of interest
    // 2) Find closest point on longitudinal circle (v-direction), using the
    //    u value from 1) and including the bounds on v given by the
    //    domain of interest

    double umin = curr_domain_of_interest.umin();
    double umax = curr_domain_of_interest.umax();
    double vmin = curr_domain_of_interest.vmin();
    double vmax = curr_domain_of_interest.vmax();

    // Find closest point on latitudinal circle (u-direction)
    double vmean = 0.5 * (vmin + vmax);
    shared_ptr<Circle> latitudinal_circle = getLatitudinalCircle(vmean);
    latitudinal_circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);

    // Find closest point on longitudinal circle (v-direction)
    shared_ptr<Circle> longitudinal_circle = getLongitudinalCircle(clo_u);
    longitudinal_circle->closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);

    // We have what we need
    return;


// #define SPHERE_CLOSEST_POINT_DEBUG
// #ifdef SPHERE_CLOSEST_POINT_DEBUG

//     // Get the SplineSurface
//     SplineSurface* sf = geometrySurface();

//     // Use the given seed to find values
//     Point tmppt(3);
//     double tmpu, tmpv, tmpdist;
//     double curr_seed[2];
//     if (seed == NULL) {
// 	curr_seed[0] = 0.5 * (umin + umax);
// 	curr_seed[1] = 0.5 * (vmin + vmax);
//     }
//     else {
// 	curr_seed[0] = seed[0];
// 	curr_seed[1] = seed[1];
//     }
//     sf->closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon,
// 		     &curr_domain_of_interest, curr_seed);
//     if (tmpdist < clo_dist - epsilon) {
// 	MESSAGE("*** Sphere::closestPoint() failed! ***");
// 	clo_u = tmpu;
// 	clo_v = tmpv;
// 	clo_pt = tmppt;
// 	clo_dist = tmpdist;
//     }

//     // Try to reseed
//     double seeds[4][2];
//     seeds[0][0] = umin;
//     seeds[0][1] = vmin;
//     seeds[1][0] = umax;
//     seeds[1][1] = vmin;
//     seeds[2][0] = umin;
//     seeds[2][1] = vmax;
//     seeds[3][0] = umax;
//     seeds[3][1] = vmax;
//     for (int i = 0; i < 4; ++i) {
// 	sf->closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon,
// 			 &curr_domain_of_interest, seeds[i]);
// 	if (tmpdist < clo_dist - epsilon) {
// 	    MESSAGE("*** Sphere::closestPoint() failed! ***");
// 	    clo_u = tmpu;
// 	    clo_v = tmpv;
// 	    clo_pt = tmppt;
// 	    clo_dist = tmpdist;
// 	}
//     }

//     delete sf;

// #endif // SPHERE_CLOSEST_POINT_DEBUG
}


//===========================================================================
void Sphere::closestBoundaryPoint(const Point& pt,
				  double&        clo_u,
				  double&        clo_v, 
				  Point&         clo_pt,
				  double&        clo_dist,
				  double epsilon,
				  const RectDomain* rd,
				  double *seed) const
//===========================================================================
{
    MESSAGE("May provide incorrect result - use with caution!");

    // This is a bit like cheating...

    SplineSurface* sf = geometrySurface();
    sf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
			     rd, seed);
    delete sf;
}


//===========================================================================
void Sphere::getBoundaryInfo(Point& pt1, Point& pt2,
			    double epsilon, SplineCurve*& cv,
			    SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
void Sphere::turnOrientation()
//===========================================================================
{
    swapParameterDirection();
}



//===========================================================================
void Sphere::swapParameterDirection()
//===========================================================================
{
    MESSAGE("swapParameterDirection() not implemented.");
}


//===========================================================================
void Sphere::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    MESSAGE("reverseParameterDirection() not implemented.");
}


//===========================================================================
bool Sphere::isDegenerate(bool& b, bool& r,
			  bool& t, bool& l, double tolerance) const
//===========================================================================
{
    bool res = false;
    b = false;
    r = false;
    t = false;
    l = false;
    if (parameterDomain().vmin() == -0.5 * M_PI) {
	b = true;
	res = true;
    }
    if (parameterDomain().vmax() == 0.5 * M_PI) {
	t = true;
	res = true;
    }
    return res;
}


//===========================================================================
bool Sphere::isBounded() const
//===========================================================================
{
  return true;
}


//===========================================================================
void Sphere::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Not this kind of degeneracy (I think)
    return;
}


//===========================================================================
void Sphere::setCoordinateAxes()
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
void Sphere::setParameterBounds(double from_upar, double from_vpar,
				double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
	THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
	THROW("First v-parameter must be strictly less than second.");
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
	THROW("u-parameters must be in [0, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
	THROW("(to_upar - from_upar) must not exceed 2pi.");
    if (from_vpar < -0.5 * M_PI || to_vpar > 0.5 * M_PI)
	THROW("v-parameters must be in [-pi/2, pi/2].");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* Sphere::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Sphere::createSplineSurface() const
//===========================================================================
{
    // Based on SISL function s1023. The SISL function creates a
    // surface with swapped parameter directions, so we need to swap
    // at the end.

    // Knot vector in latitude direction
    double et1[8];
    for (int ki = 0; ki < 3; ki++)
	et1[ki] = -0.5 * M_PI;
    for (int ki = 3; ki < 5; ki++)
	et1[ki] = 0.0;
    for (int ki = 5; ki < 8; ki++)
	et1[ki] = 0.5 * M_PI;

    // Knot vector in longitude direction
    double et2[12];
    for (int ki = 0; ki < 3; ki++)
	et2[ki] = (double)0.;
    for (int ki = 0; ki < 4; ki++) {
	et2[3 + 2*ki]     = (ki + 1)*0.5*M_PI;
	et2[3 + 2*ki + 1] = (ki + 1)*0.5*M_PI;
    }
    et2[11] = 2.0 * M_PI;


    // Coefficients
    double rcoef[180];
    double weight = 1.0 / sqrt(2.0);
    double w1, w2;
    double x_comp, y_comp, z_comp;
    Point x_axis = radius_ * x_axis_;
    Point y_axis = radius_ * y_axis_;
    Point z_axis = radius_ * z_axis_;
    for (int ki = 0; ki < 9; ki++) {
	if (ki == 1 || ki == 3 || ki == 5 || ki == 7)
	    w2 = weight;
	else
	    w2 = (double)1.;

	if (ki == 0 || ki == 1 || ki == 7 || ki == 8)
	    x_comp = (double)1.;
	else if (ki == 3 || ki == 4 || ki == 5)
	    x_comp = - (double)1.;
	else
	    x_comp = (double)0.;

	if (ki == 1 || ki == 2 || ki == 3)
	    y_comp = (double)1.;
	else if (ki == 5 || ki == 6 || ki == 7)
	    y_comp = - (double)1.;
	else
	    y_comp = (double)0.;

	for (int kj = 0; kj < 5; kj++) {
	    if (kj == 1 || kj == 3)
		w1 = weight;
	    else
		w1 = (double)1.;

	    if (kj == 0 || kj == 1)
		z_comp = - (double)1.; // Note: s1023 uses +1.0 here
	    else if (kj == 3 || kj == 4)
		z_comp = (double)1.; // Note: s1023 uses -1.0 here
	    else
		z_comp = (double)0.;

	    w1 *= w2;

	    if (kj == 0 || kj == 4) {
		for (int kl = 0; kl < 3; kl++) {
		    rcoef[4*(ki*5 + kj) + kl] 
			= w1*(location_[kl] + z_comp*z_axis[kl]);
		}
	    }
	    else { 
		for (int kl = 0; kl < 3; kl++) {
		    rcoef[4*(ki*5 + kj) + kl] 
			=  w1*(location_[kl] + x_comp*x_axis[kl]
			       + y_comp*y_axis[kl] + z_comp*z_axis[kl]);
		}
	    }
	    rcoef[4*(ki*5 + kj) + 3] = w1;
	}
    }

    int ncoefs1 = 5;
    int ncoefs2 = 9;
    int order1 = 3;
    int order2 = 3;
    int dim = 3;
    bool rational = true;
    SplineSurface surface(ncoefs1, ncoefs2, order1, order2,
			  et1, et2, rcoef, dim, rational);
    surface.swapParameterDirection();

    // Extract subpatch. We need all this because 'surface' is not
    // arc-length parametrized neither in the u- or v-directions.
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    Point llpt, urpt, tmppt;
    point(llpt, 0.0, vmin);
    point(urpt, umax - umin, vmax);
    double llu, llv, uru, urv, tmpdist;
    double epsilon = 1.0e-10;
    double seed[2];
    seed[0] = 0.0;
    seed[1] = vmin;  
    surface.closestPoint(llpt, llu, llv, tmppt, tmpdist, epsilon, NULL, seed);
    llu = 0.0;
    seed[0] = umax - umin;
    seed[1] = vmax; 
    surface.closestPoint(urpt, uru, urv, tmppt, tmpdist, epsilon, NULL, seed);

    if (0.5 * M_PI - urv < epsilon && vmin > -0.5 * M_PI + epsilon) {
	// Possible trouble with u parameter - trying lr instead - but
	// only if lower edge is not degenerate
	point(urpt, umax - umin, vmin);
	seed[0] = umax - umin;
	seed[1] = vmin;
	surface.closestPoint(urpt, uru, urv, tmppt, tmpdist, epsilon, NULL, seed);
	urv = 0.5 * M_PI;
    }
    if (uru < epsilon && umax - umin == 2.0 * M_PI) {
	uru = 2.0 * M_PI;
    }
    SplineSurface* subpatch = surface.subSurface(llu, llv, uru, urv);
    subpatch->basis_u().rescale(umin, umax);
    subpatch->basis_v().rescale(vmin, vmax);
    GeometryTools::translateSplineSurf(-location_, *subpatch);
    GeometryTools::rotateSplineSurf(z_axis_, umin, *subpatch);
    GeometryTools::translateSplineSurf(location_, *subpatch);

    return subpatch;
}


//===========================================================================
shared_ptr<Circle> Sphere::getLatitudinalCircle(double vpar) const
//===========================================================================
{
    Point centre = location_ + radius_ * sin(vpar) * z_axis_;
    double radius = radius_ * cos(vpar);
    shared_ptr<Circle> circle(new Circle(radius, centre, z_axis_, x_axis_));
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->setParamBounds(umin, umax);
    return circle;
}


//===========================================================================
shared_ptr<Circle> Sphere::getLongitudinalCircle(double upar) const
//===========================================================================
{
    Point udir = cos(upar) * x_axis_ + sin(upar) * y_axis_;
    Point newz = udir.cross(z_axis_);
    shared_ptr<Circle> circle(new Circle(radius_, location_,
					 newz, udir));
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    circle->setParamBounds(vmin, vmax);
    return circle;
}


//===========================================================================


} // namespace Go
