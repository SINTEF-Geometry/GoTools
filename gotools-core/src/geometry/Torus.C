//==========================================================================
//                                                                          
// File: Torus.C                                                             
//                                                                          
// Created: Wed Feb 25 13:22:43 2009                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Torus.C,v 1.3 2009-05-13 07:30:53 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include <vector>
#include <limits>


using Go::SweepSurfaceCreator;
using std::vector;
using std::cout;
using std::endl;
using std::streamsize;
using std::swap;


namespace Go
{


// Constructor.
//===========================================================================
Torus::Torus(double major_radius, double minor_radius,
	     Point location, Point z_axis, Point x_axis,
	     bool select_outer, bool isSwapped)
    : major_radius_(major_radius), minor_radius_(minor_radius),
    location_(location), z_axis_(z_axis), x_axis_(x_axis),
      is_degenerate_torus_(false), select_outer_(select_outer), phi_(-1.0)
//===========================================================================
{
    // The use of the z- and x-axis in this constructor - and in the
    // definition of a torus - is based on the use of the
    // 'axis2_placement_3d' entity in STEP.

    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }
    setCoordinateAxes();
    setDegenerateInfo();
    setDefaultDomain();

    if (isSwapped)
        swapParameterDirection();
}


// Destructor
//===========================================================================
Torus::~Torus()
//===========================================================================
{
}

//===========================================================================
void Torus::read (std::istream& is)
//===========================================================================
{
    // Note on the data format: Wether the torus is degenerate or not
    // depends automatically on the minor radius being greater than
    // the major radius. The select_outer_ flag - which makes sense
    // only for degenerate tori - is included in the format for
    // reasons of convenience. A non-standard parameter domain is not
    // "supported" within the data format at this time, however.

    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim, select_outer_int;
    is >> dim;
    if (dim != 3)
	THROW("Dimension must be 3.");
    location_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> major_radius_
       >> minor_radius_
       >> location_
       >> z_axis_
       >> x_axis_
       >> select_outer_int;
    if (select_outer_int == 0)
	select_outer_ = false;
    else if (select_outer_int == 1)
	select_outer_ = true;
    else
	THROW("Unknown input for select_outer - must be 0 or 1");

    setCoordinateAxes();
    setDegenerateInfo();

    // "Reset" swapping
    isSwapped_ = false;

    // NB: Mind the sequence of parameters!
    double from_upar, from_vpar, to_upar, to_vpar;
    is >> from_upar >> to_upar
        >> from_vpar >> to_vpar;

    // Need to take care of rounding errors: If upars are "roughly"
    // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
    const double pareps = 1.0e-4; // This is admittedly arbitrary...
    if (fabs(from_upar) < pareps && fabs(to_upar - 2.0*M_PI) < pareps) {
        from_upar = 0.0;
        to_upar = 2.0 * M_PI;
    }
    if (fabs(from_vpar) < pareps && fabs(to_vpar - 2.0*M_PI) < pareps) {
        from_vpar = 0.0;
        to_vpar = 2.0 * M_PI;
    }
    if (is_degenerate_torus_) {
        if (select_outer_) {
            if (fabs(from_vpar + phi_) < pareps)
                from_vpar = -phi_;
            if (fabs(to_vpar - phi_) < pareps)
                to_vpar = phi_;
        }
        else {
            if (fabs(from_vpar - phi_) < pareps)
                from_vpar = phi_;
            if (fabs(to_vpar - 2.0 * M_PI + phi_) < pareps)
                to_vpar = 2.0 * M_PI - phi_;
        }
    }

    setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);

    // Swapped flag
    int isSwapped; // 0 or 1
    is >> isSwapped;
    if (isSwapped == 0) {
        // Do nothing
    }
    else if (isSwapped == 1) {
        swapParameterDirection();
    }
    else {
        THROW("Swapped flag must be 0 or 1");
    }
}


//===========================================================================
void Torus::write(std::ostream& os) const
//===========================================================================
{
    // Note on the data format: See comments for read().

    streamsize prev = os.precision(15);
    os << dimension() << endl
       << major_radius_ << endl
       << minor_radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;
    if (select_outer_)
	os << "1" << endl;
    else
	os << "0" << endl;

    // NB: Mind the parameter sequence!
    os << domain_.umin() << " " << domain_.umax() << endl
       << domain_.vmin() << " " << domain_.vmax() << endl;

    if (!isSwapped()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
    }

    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
int Torus::dimension() const
//===========================================================================
{
    return location_.dimension();
}

//===========================================================================
ClassType Torus::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
BoundingBox Torus::boundingBox() const
//===========================================================================
{
    // A rather unefficient hack...
    SplineSurface* tmp = geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
Torus* Torus::clone() const
//===========================================================================
{
    Torus* torus = new Torus(major_radius_, minor_radius_,
		       location_, z_axis_, x_axis_, 
                       select_outer_, isSwapped_);
    torus->domain_ = domain_;
    return torus;
}    
    
//===========================================================================
const RectDomain& Torus::parameterDomain() const
//===========================================================================
{
    if (!isSwapped())
        return domain_;

    // If parameters are swapped, we must make a swapped domain
    Array<double, 2> ll, ur;
    ll[0] = domain_.vmin();
    ll[1] = domain_.umin();
    ur[0] = domain_.vmax();
    ur[1] = domain_.umax();
    orientedDomain_ = RectDomain(ll, ur);
    return orientedDomain_;
}


//===========================================================================
std::vector<CurveLoop> 
Torus::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("allBoundaryLoops() not implemented. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Torus::normalCone() const
//===========================================================================
{
    RectDomain domain = parameterDomain();
    double umin = domain.umin();
    double umax = domain.umax();
    double vmin = domain.vmin();
    double vmax = domain.vmax();
    Point dir;
    normal(dir, 0.5*(umin+umax), 0.5*(vmin+vmax));

    Point ll, ur;
    normal(ll, umin, vmin);
    normal(ur, umax, vmax);
    double angle = std::max(dir.angle(ll), dir.angle(ur));

    return DirectionCone(dir, angle);
}


//===========================================================================
DirectionCone Torus::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    if (isSwapped())
        pardir_is_u = !pardir_is_u;

    if (pardir_is_u) {
	shared_ptr<Circle> circle = getMajorCircle(0.0);
	return circle->directionCone();
    }
    else {
	DirectionCone normals = normalCone();
	double angle = normals.angle();
        RectDomain domain = parameterDomain();
	double umid = 0.5 * (domain.umax() + domain.umin());
	double vmid = 0.5 * (domain.vmax() + domain.vmin());
	vector<Point> pts(3);
	point(pts, umid, vmid, 1);
	return DirectionCone(pts[2], angle);
    }
}


//===========================================================================
void Torus::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar); // In case of swapped
    pt = location_
	+ (major_radius_ + minor_radius_ * cos(vpar)) * (cos(upar) * x_axis_ 
							 + sin(upar) * y_axis_)
	+ minor_radius_ * sin(vpar) * z_axis_;
}


//===========================================================================
void Torus::point(std::vector<Point>& pts, 
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

    // Swap parameters, if needed
    getOrientedParameters(upar, vpar);
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First derivatives
    double cosu = cos(upar);
    double sinu = sin(upar);
    double cosv = cos(vpar);
    double sinv = sin(vpar);
    pts[ind1] = (major_radius_ + minor_radius_ * cosv)
	* (-sinu * x_axis_ + cosu * y_axis_);
    pts[ind2] = minor_radius_ 
	* (-sinv * (cosu * x_axis_ + sinu * y_axis_)
	   + cosv * z_axis_);
    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void Torus::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    // This formula holds for both regular and degenerate tori.
    getOrientedParameters(upar, vpar);
    n = cos(vpar) * (cos(upar) * x_axis_ + sin(upar) * y_axis_)
	+ sin(vpar) * z_axis_;
    if (isSwapped())
        n *= -1.0;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Torus::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
Torus* Torus::subSurface(double from_upar, double from_vpar,
			 double to_upar, double to_vpar,
			 double fuzzy) const
//===========================================================================
{
    Torus* torus = clone();
    torus->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return torus;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
Torus::subSurfaces(double from_upar, double from_vpar,
		  double to_upar, double to_vpar,
		  double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Torus> torus(subSurface(from_upar, from_vpar,
				       to_upar, to_vpar));
    res.push_back(torus);
    return res;
}


//===========================================================================
double 
Torus::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    MESSAGE("nextSegmentVal() doesn't make sense. Returning arbitrarily 0.0.");
    return 0.0;
}


//===========================================================================
void Torus::closestPoint(const Point& pt,
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
    RectDomain curr_domain_of_interest = parameterDomain();
    if (domain_of_interest != NULL) {
	curr_domain_of_interest.intersectWith(*domain_of_interest);
    }

    // Algorithm:
    // 1) Find closest point on major circle (u-direction), including
    //    the bounds on u given by the domain of interest
    // 2) Find closest point on minor circle (v-direction), using the
    //    u value from 1) and including the bounds on v given by the
    //    domain of interest

    double umin = curr_domain_of_interest.umin();
    double umax = curr_domain_of_interest.umax();
    double vmin = curr_domain_of_interest.vmin();
    double vmax = curr_domain_of_interest.vmax();

    getOrientedParameters(umin, vmin);
    getOrientedParameters(umax, vmax);

    // Find closest point on major circle (u-direction)
    double vmean = 0.5 * (vmin + vmax);
    shared_ptr<Circle> major_circle = getMajorCircle(vmean);
    major_circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);

    // Find closest point on minor circle (v-direction)
    shared_ptr<Circle> minor_circle = getMinorCircle(clo_u);
    minor_circle->closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);

    // We have what we need
    getOrientedParameters(clo_u, clo_v);
    return;


// #define TORUS_CLOSEST_POINT_DEBUG
// #ifdef TORUS_CLOSEST_POINT_DEBUG

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
// 	MESSAGE("*** Torus::closestPoint() failed! ***");
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
// 	    MESSAGE("*** Torus::closestPoint() failed! ***");
// 	    clo_u = tmpu;
// 	    clo_v = tmpv;
// 	    clo_pt = tmppt;
// 	    clo_dist = tmpdist;
// 	}
//     }

//     delete sf;

// #endif // TORUS_CLOSEST_POINT_DEBUG

}


//===========================================================================
void Torus::closestBoundaryPoint(const Point& pt,
				double&        clo_u,
				double&        clo_v, 
				Point&       clo_pt,
				double&        clo_dist,
				double epsilon,
				const RectDomain* rd,
				double *seed) const
//===========================================================================
{
    // Algorithm:
    // 1) Find closest point overall
    // 2) If on boundary, return
    // 3) Else, snap to boundary

    // Find closest point overall
    closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);

    // Check if on boundary
    RectDomain dom = parameterDomain();
    double umin = dom.umin();
    double umax = dom.umax();
    double vmin = dom.vmin();
    double vmax = dom.vmax();
    if (fabs(clo_u - umin) < epsilon ||
        fabs(clo_u - umax) < epsilon ||
        fabs(clo_v - vmin) < epsilon ||
        fabs(clo_v - vmax) < epsilon) {
            return;
    }

    // Overall closest point is in the inner of the patch
    double clo_u_inner = clo_u;
    double clo_v_inner = clo_v;

    // Fist on boundary
    clo_u = umin;
    point(clo_pt, umin, clo_v_inner);
    clo_dist = pt.dist(clo_pt);

    // Second
    Point clo_pt_tmp;
    point(clo_pt_tmp, umax, clo_v_inner);
    double clo_dist_tmp = pt.dist(clo_pt_tmp);
    if (clo_dist_tmp < clo_dist) {
        clo_pt = clo_pt_tmp;
        clo_u = umax;
        clo_v = clo_v_inner;
        clo_dist = clo_dist_tmp;
    }

    // Third
    point(clo_pt_tmp, clo_u_inner, vmin);
    clo_dist_tmp = pt.dist(clo_pt_tmp);
    if (clo_dist_tmp < clo_dist) {
        clo_pt = clo_pt_tmp;
        clo_u = clo_u_inner;
        clo_v = vmin;
        clo_dist = clo_dist_tmp;
    }

    // Fourth
    point(clo_pt_tmp, clo_u_inner, vmax);
    clo_dist_tmp = pt.dist(clo_pt_tmp);
    if (clo_dist_tmp < clo_dist) {
        clo_pt = clo_pt_tmp;
        clo_u = clo_u_inner;
        clo_v = vmax;
        clo_dist = clo_dist_tmp;
    }

    return;
}


//===========================================================================
void Torus::getBoundaryInfo(Point& pt1, Point& pt2,
			    double epsilon, SplineCurve*& cv,
			    SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
bool Torus::isDegenerate(bool& b, bool& r,
			bool& t, bool& l, double tolerance) const
//===========================================================================
{
    bool res = false;
    b = false;
    r = false;
    t = false;
    l = false;

    if (!is_degenerate_torus_)
	return res;

    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    if (select_outer_) {
	if (fabs(vmin + phi_) < tolerance) {
	    b = true;
	    res = true;
	}
	if (fabs(vmax - phi_) < tolerance) {
	    t = true;
	    res = true;
	}
    }
    else {
	if (fabs(vmin - phi_) < tolerance) {
	    b = true;
	    res = true;
	}
	if (fabs(vmax - 2.0 * M_PI + phi_) < tolerance) {
	    t = true;
	    res = true;
	}
    }
    if (isSwapped()) {
        swap(b, l);
        swap(t, r);
    }
    return res;
}


//===========================================================================
bool Torus::isBounded() const
//===========================================================================
{
  return true;
}


//===========================================================================
bool Torus::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    closed_dir_u = (domain_.umax() - domain_.umin() == 2.0*M_PI);
    closed_dir_v = (domain_.vmax() - domain_.vmin() == 2.0*M_PI);
    if (isSwapped())
        swap(closed_dir_u, closed_dir_v);
    return (closed_dir_u || closed_dir_v);
}


//===========================================================================
void Torus::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Not this kind of degeneracy
    return;
}


//===========================================================================
void Torus::setSelectOuter(bool select_outer)
//===========================================================================
{
    if (!is_degenerate_torus_)
	return;

    if (select_outer_ == select_outer)
	return;

    select_outer_ = select_outer;
    setDefaultDomain();
}


//===========================================================================
void Torus::setCoordinateAxes()
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
void Torus::setDegenerateInfo()
//===========================================================================
{
    is_degenerate_torus_ = (minor_radius_ > major_radius_);
    if (is_degenerate_torus_) {
	phi_ = acos(-major_radius_ / minor_radius_);
    }
}


//===========================================================================
void Torus::setDefaultDomain()
//===========================================================================
{
    double umin = 0.0;
    double umax = 2.0 * M_PI;
    double vmin = 0.0;
    double vmax = 2.0 * M_PI;
    if (is_degenerate_torus_) {
	if (select_outer_) {
	    vmin = -phi_;
	    vmax = phi_;
	}
	else {
	    vmin = phi_;
	    vmax = 2.0 * M_PI - phi_;
	}
    }
    Array<double, 2> ll(umin, vmin);
    Array<double, 2> ur(umax, vmax);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Torus::setParameterBounds(double from_upar, double from_vpar,
			       double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
	THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
	THROW("First v-parameter must be strictly less than second.");

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
	THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
	THROW("(to_upar - from_upar) must not exceed 2pi.");
    if (is_degenerate_torus_) {
	if (select_outer_) {
	    if (from_vpar < -phi_)
		THROW("First v-parameter must be >= -phi");
	    if (to_vpar > phi_)
		THROW("Second v-parameter must be <= phi");
	}
	else {
	    if (from_vpar < phi_)
		THROW("First v-parameter must be >= phi");
	    if (to_vpar > 2.0 * M_PI - phi_)
		THROW("Second v-parameter must be <= 2pi - phi");
	}
    }
    else {
	if (from_vpar < -2.0 * M_PI || to_vpar > 2.0 * M_PI)
	    THROW("v-parameters must be in [-2pi, 2pi].");
	if (to_vpar - from_vpar > 2.0 * M_PI)
	    THROW("(to_vpar - from_vpar) must not exceed 2pi.");
    }

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* Torus::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Torus::createSplineSurface() const
//===========================================================================
{
    double umin = domain_.umin();
    double umax = domain_.umax();

    shared_ptr<Circle> circle = getMinorCircle(umin);
    shared_ptr<SplineCurve> sccircle(circle->geometryCurve());
    double angle = umax - umin;

    SplineSurface* sstorus 
	= SweepSurfaceCreator::rotationalSweptSurface(*sccircle, angle,
						      location_, z_axis_);
    sstorus->basis_u().rescale(umin, umax);

    if (isSwapped())
        sstorus->swapParameterDirection();

    return sstorus;
}


//===========================================================================
shared_ptr<Circle> Torus::getMajorCircle(double vpar) const
//===========================================================================
{
    Point centre = location_ + minor_radius_ * sin(vpar) * z_axis_;
    double radius = major_radius_ + minor_radius_ * cos(vpar);
    shared_ptr<Circle> circle(new Circle(radius, centre, z_axis_, x_axis_));
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->setParamBounds(umin, umax);
    return circle;
}


//===========================================================================
shared_ptr<Circle> Torus::getMinorCircle(double upar) const
//===========================================================================
{
    Point udir = cos(upar) * x_axis_ + sin(upar) * y_axis_;
    Point centre = location_ + major_radius_ * udir;
    Point newz = udir.cross(z_axis_);
    shared_ptr<Circle> circle(new Circle(minor_radius_, centre,
					 newz, udir));
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    circle->setParamBounds(vmin, vmax);
    return circle;
}


//===========================================================================


} // namespace Go
