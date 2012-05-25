//===========================================================================
//                                                                           
// File: Plane.C                                                  
//                                                                           
// Created: Wed Mar 14 17:32:59 2001                                         
//                                                                           
// Author: Vibeke Skytt, Jan Thomassen
//                                                                           
// Revision: $Id: Plane.C,v 1.13 2009-02-17 13:10:03 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/SplineSurface.h"
#include <vector>
#include <limits>


using std::vector;
using std::endl;
using std::numeric_limits;
using std::streamsize;


namespace Go
{


// Constructor. Input is location and normal
//===========================================================================
Plane::Plane(Point location, Point normal)
    : location_(location), normal_(normal), vec1_(1.0, 0.0, 0.0)
//===========================================================================
{
    if (location.dimension() != 3)
        return;
    setSpanningVectors();

    setParameterBounds(-numeric_limits<double>::infinity(),
                       -numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity());
}

//===========================================================================
Plane::Plane(Point location, Point normal, Point x_axis)
    : location_(location), normal_(normal), vec1_(x_axis)
//===========================================================================
{
    if (location.dimension() != 3)
        return;
    setSpanningVectors();

    setParameterBounds(-numeric_limits<double>::infinity(),
                       -numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity());
}

// Constructor. Input is coefficients of implicit equation
//===========================================================================
Plane::Plane(double a, double b, double c, double d)
//===========================================================================
{
    double tol = 1.0e-8;
    // double frac = d/sqrt(a*a + b*b + c*c);
    // double ta = (fabs(a) < tol) ? 0.0 : frac/a;
    // double tb = (fabs(b) < tol) ? 0.0 : frac/b;
    // double tc = (fabs(c) < tol) ? 0.0 : frac/c;
    int fac = 3;
    if (fabs(a) < tol)
      fac--;
    if (fabs(b) < tol)
      fac--;
    if (fabs(c) < tol)
      fac--;
    double ta =  (fabs(a) < tol) ? 0.0 : d/(fac*a);
    double tb =  (fabs(b) < tol) ? 0.0 : d/(fac*b);
    double tc =  (fabs(c) < tol) ? 0.0 : d/(fac*c);

    location_ = Point(ta, tb, tc);
    normal_ = Point(a, b, c);
    vec1_ = Point(1.0, 0.0, 0.0);
    setSpanningVectors();

    setParameterBounds(-numeric_limits<double>::infinity(),
                       -numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity(),
                       numeric_limits<double>::infinity());
}

// Destructor
//===========================================================================
Plane::~Plane()
//===========================================================================
{
}

//===========================================================================
void Plane::read (std::istream& is)
//===========================================================================
{
    // NB: Parameter sequence in the g2 file format is different
    // than the argument list of the setParameterBounds() function!

    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    location_.resize(dim);
    normal_.resize(dim);
    vec1_.resize(dim);
    is >> location_
       >> normal_
       >> vec1_;
    setSpanningVectors();

    int isBounded; 
    is >> isBounded;
    if (isBounded == 0) {
        // Unbounded
        setParameterBounds(-numeric_limits<double>::infinity(),
                           -numeric_limits<double>::infinity(),
                           numeric_limits<double>::infinity(),
                           numeric_limits<double>::infinity());
    }
    else if (isBounded == 1) {
        // NB: See comment on parameter sequence above!
        double from_upar, from_vpar, to_upar, to_vpar;
        is >> from_upar >> to_upar
           >> from_vpar >> to_vpar;
        setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    }
    else {
        THROW("Bounded flag must be 0 or 1");
    }

}


//===========================================================================
void Plane::write(std::ostream& os) const
//===========================================================================
{
    // NB: Parameter sequence in the g2 file format is different
    // than the argument list of the setParameterBounds() function!

    streamsize prev = os.precision(15);
    os << dimension() << endl
       << location_ << endl
       << normal_ << endl
       << vec1_ << endl;

    if (!isBounded()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
        os << domain_.umin() << " " << domain_.umax() << endl
           << domain_.vmin() << " " << domain_.vmax() << endl;
    }
    os.precision(prev);   // Reset precision to it's previous value
}

//===========================================================================
int Plane::dimension() const
//===========================================================================
{
    return location_.dimension();
}

//===========================================================================
ClassType Plane::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
BoundingBox Plane::boundingBox() const
//===========================================================================
{
    BoundingBox box(3);

    // If the plane is unbounded, return an empty box.
    if (!isBounded())
        return box;

    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    vector<Point> pts;
    Point pt;
    point(pt, umin, vmin);
    pts.push_back(pt);
    point(pt, umax, vmin);
    pts.push_back(pt);
    point(pt, umin, vmax);
    pts.push_back(pt);
    point(pt, umax, vmax);
    pts.push_back(pt);
    box.setFromPoints(pts);

    return box;
}

//===========================================================================
const Domain& Plane::parameterDomain() const
//===========================================================================
{
    // Return the domain. This domain is either set by
    // setParameterBounds(), or it is uninitialized.
    return domain_;
}


//===========================================================================
std::vector<CurveLoop> 
Plane::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    // Does not make sense. Returns empty vector.
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Plane::normalCone() const
//===========================================================================
{
    return DirectionCone(normal_);
}


//===========================================================================
DirectionCone Plane::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    if (pardir_is_u)
        return DirectionCone(vec1_);
    else
        return DirectionCone(vec2_);
}


//===========================================================================
void Plane::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    pt = location_ + upar * vec1_ + vpar * vec2_;
}


//===========================================================================
void Plane::point(std::vector<Point>& pts, 
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

    point(pts[0], upar, vpar);
    if (derivs == 0)
        return;

    // Derivatives are just the spanning vectors
    pts[1] = vec1_;
    pts[2] = vec2_;

    // Second order and higher derivatives vanish. They are already
    // set to zero, so we return.
    return;

}


//===========================================================================
void Plane::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    n = normal_;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Plane::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    // Not yet implemented
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
Plane* Plane::subSurface(double from_upar, double from_vpar,
                         double to_upar, double to_vpar,
                         double fuzzy) const
//===========================================================================
{
    Plane* plane = clone();
    plane->setParameterBounds(from_upar, from_vpar,
                              to_upar, to_vpar);
    return plane;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
Plane::subSurfaces(double from_upar, double from_vpar,
                   double to_upar, double to_vpar,
                   double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Plane> plane(subSurface(from_upar, from_vpar,
                                       to_upar, to_vpar));
    res.push_back(plane);
    return res;
}


//===========================================================================
double 
Plane::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    // Does not make sense. Return arbitrarily zero.
    return 0.0;
}


//===========================================================================
void Plane::closestPoint(const Point& pt,
                         double&        clo_u,
                         double&        clo_v, 
                         Point&       clo_pt,
                         double&        clo_dist,
                         double         epsilon,
                         const RectDomain* domain_of_interest,
                         double   *seed) const
//===========================================================================
{
    clo_pt = projectPoint(pt);
    clo_u = (clo_pt - location_) * vec1_;
    clo_v = (clo_pt - location_) * vec2_;
    if (clo_u < domain_.umin())
        clo_u = domain_.umin();
    if (clo_u > domain_.umax())
        clo_u = domain_.umax();
    if (clo_v < domain_.vmin())
        clo_v = domain_.vmin();
    if (clo_v > domain_.vmax())
        clo_v = domain_.vmax();
    point(clo_pt, clo_u, clo_v);
    clo_dist = (clo_pt - pt).length();
}


//===========================================================================
void Plane::closestBoundaryPoint(const Point& pt,
                                 double&        clo_u,
                                 double&        clo_v, 
                                 Point&       clo_pt,
                                 double&        clo_dist,
                                 double epsilon,
                                 const RectDomain* rd,
                                 double *seed) const
//===========================================================================
{
  // Does not make sense if the plane is unbounded in all four directions.

  Point proj_pt = projectPoint(pt);
  double proj_u = (proj_pt - location_) * vec1_;
  double proj_v = (proj_pt - location_) * vec2_;
  if (proj_u < domain_.umin())
    proj_u = domain_.umin();
  if (proj_u > domain_.umax())
    proj_u = domain_.umax();
  if (proj_v < domain_.vmin())
    proj_v = domain_.vmin();
  if (proj_v > domain_.vmax())
    proj_v = domain_.vmax();

  bool best_found = false;
  double best_dist;
  double inf = numeric_limits<double>::infinity();
  clo_u = proj_u;
  clo_v = proj_v;

  if (domain_.umin() > -inf && (!best_found || best_dist > proj_u - domain_.umin()))
    {
      best_found = true;
      clo_u = domain_.umin();
      best_dist = proj_u - clo_u;
    }
  if (domain_.umax() < inf && (!best_found || best_dist > domain_.umax() - proj_u))
    {
      best_found = true;
      clo_u = domain_.umax();
      best_dist = clo_u - proj_u;
    }
  if (domain_.vmin() > -inf && (!best_found || best_dist > proj_v - domain_.vmin()))
    {
      best_found = true;
      clo_v = domain_.vmin();
      best_dist = proj_v - clo_v;
    }
  if (domain_.vmax() < inf && (!best_found || best_dist > domain_.vmax() - proj_v))
    {
      best_found = true;
      clo_v = domain_.vmax();
      best_dist = clo_v - proj_v;
    }

  if (!best_found)
    THROW("Can not find closestBoundaryPoint(), plane has no boundary");

  point(clo_pt, clo_u, clo_v);
  clo_dist = (clo_pt - pt).length();
}


//===========================================================================
void Plane::getBoundaryInfo(Point& pt1, Point& pt2,
                            double epsilon, SplineCurve*& cv,
                            SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    // Does not make sense. Leave input/output unchanged.
    MESSAGE("getBoundaryInfo() not yet implemented.");
}


//===========================================================================
void Plane::turnOrientation()
//===========================================================================
{
    // Calling swapParameterDirection() may cause problems with finite
    // parameter bounds.
    swapParameterDirection();
}



//===========================================================================
void Plane::swapParameterDirection()
//===========================================================================
{
    // A Plane has a canonical parametrization. Therefore, for a bounded
    // Plane it doesn't make sense to swap parameter directions. How do we 
    // handle this? Assuming an "unbounded swap" will work for now... @jbt

    // Spanning vectors
    Point tmp = vec1_;
    vec1_ = vec2_;
    vec2_ = tmp;
    normal_ = -1.0 * normal_;

    if (isBounded()) {
        MESSAGE("Not properly implemented - check parameter bounds");
    }
}


//===========================================================================
void Plane::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    if (direction_is_u)
        vec1_ = -1.0 * vec1_;
    else
        vec2_ = -1.0 * vec2_;

    normal_ = -1.0 * normal_;

    MESSAGE("Not properly implemented - check parameter bounds");
}


//===========================================================================
bool Plane::isDegenerate(bool& b, bool& r,
                         bool& t, bool& l, double tolerance) const
//===========================================================================
{
    // The canonical parametrization is not degenerate
    b = false;
    r = false;
    t = false;
    l = false;
    return false;
}


//===========================================================================
void Plane::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Not this kind of degeneracy
    return;
}

//===========================================================================
Point Plane::projectPoint(const Point& pnt) const
//===========================================================================
{
    Point vec = pnt - location_;
    Point projected = pnt - (vec*normal_)*normal_;
    return projected;
}


//===========================================================================
double Plane::distance(const Point& pnt) const
//===========================================================================
{
    return pnt.dist(projectPoint(pnt));
    MESSAGE("Computing distance to projected "
            "- not necessarily closest - point.");
}
//===========================================================================
void Plane::setSpanningVectors()
//===========================================================================
{
    // The vectors vec1_, vec2_, and normal_ define a right-handed
    // coordinate system.

    Point tmp = vec1_ - (vec1_ * normal_) * normal_;
    if (tmp.length() == 0.0) {
        // normal_ and vec1_ are parallel. Try x-axis along coordinate
        // x-direction.
        vec1_ = Point(1.0, 0.0, 0.0);
        tmp = vec1_ - (vec1_ * normal_) * normal_;
        if (tmp.length() == 0.0) {
            // Still parallel...? Try the y-direction.
            vec1_ = Point(0.0, 1.0, 0.0);
            tmp = vec1_ - (vec1_ * normal_) * normal_;
            ASSERT(tmp.length() != 0.0);
        }
    }

    vec1_ = tmp;
    vec2_ = normal_.cross(vec1_);
    normal_.normalize();
    vec1_.normalize();
    vec2_.normalize();

}


//===========================================================================
void Plane::setParameterBounds(double from_upar, double from_vpar,
                               double to_upar, double to_vpar)
//===========================================================================
{
    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* Plane::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Plane::createSplineSurface() const
//===========================================================================
{
    int ncoefsu = 2;
    int ncoefsv = 2;
    int ordu = 2;
    int ordv = 2;

    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();

    // Handle the case if not bounded
    if (!isBounded()) {
        double max = 1.0e8; // "Large" number...
        if (umin == -numeric_limits<double>::infinity())
            umin = -max;
        if (umax == numeric_limits<double>::infinity())
            umax = max;
        if (vmin == -numeric_limits<double>::infinity())
            vmin = -max;
        if (vmax == numeric_limits<double>::infinity())
            vmax = max;
    }

    vector<double> knotsu(4);
    knotsu[0] = umin;
    knotsu[1] = umin;
    knotsu[2] = umax;
    knotsu[3] = umax;
    vector<double> knotsv(4);
    knotsv[0] = vmin;
    knotsv[1] = vmin;
    knotsv[2] = vmax;
    knotsv[3] = vmax;
    vector<double> coefs(12);
    Point c0, c1, c2, c3;
    point(c0, umin, vmin);
    point(c1, umax, vmin);
    point(c2, umin, vmax);
    point(c3, umax, vmax);
    for (int d = 0; d < 3; ++d) {
        coefs[d] = c0[d];
        coefs[3+d] = c1[d];
        coefs[6+d] = c2[d];
        coefs[9+d] = c3[d];
    }
    int dim = 3;

    return new SplineSurface(ncoefsu, ncoefsv, ordu, ordv,
                             knotsu.begin(), knotsv.begin(), 
                             coefs.begin(), dim);
}


//===========================================================================
bool Plane::isBounded() const
//===========================================================================
{
    double inf = numeric_limits<double>::infinity();

    return domain_.umin() > -inf && domain_.vmin() > -inf
        && domain_.umax() < inf && domain_.vmax() < inf;

}



//===========================================================================
Plane* Plane::intersect(const RotatedBox& bd_box) const
//===========================================================================
{
    Plane* int_plane = NULL;
//     if (!isBounded()) {
// 	return int_plane; // Well, we could of course also treat the
// 			  // bounded intersection.
//     }

    // We check if the plane is above or below bd_box.
    Point high = bd_box.high_rot();
    Point dir1 = high - location_;
    double angle1 = dir1.angle(normal_);
    bool above = (angle1 < 0.5*M_PI);
    if (above)
        return int_plane; // Empty intersection.

    Point low = bd_box.low_rot();
    Point dir2 = high - location_;
    double angle2 = dir2.angle(normal_);
    bool below = (angle2 > 0.5*M_PI);
    if (below)
        return int_plane; // Empty intersection.

    // As the rotated box is using the same coordinate system as the
    // plane, we do not need to rotate the box when intersecting. We
    // need only find the umin, umax, vmin & vmax.

    // We project the high and low into the plane.
    double clo_u_high, clo_v_high, clo_u_low, clo_v_low, clo_dist;
    double epsdummy = 1.0e-06; // Well, not used ...
    Point clo_pt_high, clo_pt_low;
    closestPoint(high, clo_u_high, clo_v_high, clo_pt_high, clo_dist, epsdummy);
    closestPoint(low, clo_u_low, clo_v_low, clo_pt_low, clo_dist, epsdummy);

    if (isBounded()) {
        // We must pick the part of the box that coincides with the
        // bounded plane, if any.
        double global_u_low = std::max(clo_u_low, domain_.umin());
        double global_v_low = std::max(clo_v_low, domain_.vmin());
        double global_u_high = std::min(clo_u_high, domain_.umax());
        double global_v_high = std::min(clo_v_high, domain_.vmax());
        if (global_u_low > global_u_high || global_v_low > global_v_high)
            return int_plane; // Empty intersection. 
        else {
            int_plane = clone();
            int_plane->setParameterBounds(global_u_low, global_v_low,
                                          global_u_high, global_v_high);
        }
    } else {
        int_plane = clone();
        int_plane->setParameterBounds(clo_u_low, clo_v_low,
                                      clo_u_high, clo_v_high);
    }

    return int_plane;
}

//===========================================================================
bool Plane::isPlanar(Point& normal, double tol)
//===========================================================================
{
  normal = normal_;

  // This surface is a plane
  return true;
}

} // namespace Go
