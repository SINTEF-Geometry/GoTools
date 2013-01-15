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
using std::swap;


namespace Go
{


// Constructor. Input is location and normal
//===========================================================================
Plane::Plane(Point location, Point normal,
             bool isSwapped)
    : location_(location), normal_(normal), vec1_(1.0, 0.0, 0.0)
//===========================================================================
{
    if (location.dimension() != 3)
        return;
    setSpanningVectors();

    double inf = numeric_limits<double>::infinity();
    setParameterBounds(-inf, -inf, inf, inf);

    if (isSwapped)
        swapParameterDirection();
}

//===========================================================================
Plane::Plane(Point location, Point normal, Point x_axis,
             bool isSwapped)
    : location_(location), normal_(normal), vec1_(x_axis)
//===========================================================================
{
    if (location.dimension() != 3)
        return;
    setSpanningVectors();

    double inf = numeric_limits<double>::infinity();
    setParameterBounds(-inf, -inf, inf, inf);

    if (isSwapped)
        swapParameterDirection();
}

// Constructor. Input is coefficients of implicit equation
//===========================================================================
Plane::Plane(double a, double b, double c, double d,
             bool isSwapped)
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

    double inf = numeric_limits<double>::infinity();
    setParameterBounds(-inf, -inf, inf, inf);

    if (isSwapped)
        swapParameterDirection();
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

    // "Reset" swapping
    isSwapped_ = false;

    // Bounded flag
    int isBounded; 
    is >> isBounded;
    if (isBounded == 0) {
        // Unbounded
        double inf = numeric_limits<double>::infinity();
        setParameterBounds(-inf, -inf, inf, inf);
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

    if (!isSwapped()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
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
    // First handle the case if not bounded
    if (!isBounded()) {
        // Create a SplineSurface with a large but finite BoundingBox
        SplineSurface* tmp = geometrySurface();
        BoundingBox box = tmp->boundingBox();
        delete tmp;
        return box;
    }

    BoundingBox box(3);

    // Call parameterDomain() to get the possibly swapped domain
    RectDomain domain = parameterDomain();
    double umin = domain.umin();
    double umax = domain.umax();
    double vmin = domain.vmin();
    double vmax = domain.vmax();
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
Plane* Plane::clone() const
//===========================================================================
{
    Plane* plane = new Plane(location_, normal_, isSwapped_);
    plane->domain_ = domain_;
    return plane;
}
    
    
//===========================================================================
const RectDomain& Plane::parameterDomain() const
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
Plane::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    // Does not make sense. Returns empty vector.
    MESSAGE("allBoundaryLoops() not implemented. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Plane::normalCone() const
//===========================================================================
{
    Point n;
    normal(n, 0.0, 0.0); // Evaluates arbitrarily in (0,0)
    return DirectionCone(n);
}


//===========================================================================
DirectionCone Plane::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    if (isSwapped())
        pardir_is_u = !pardir_is_u;

    if (pardir_is_u)
        return DirectionCone(vec1_);
    else
        return DirectionCone(vec2_);
}


//===========================================================================
void Plane::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar); // In case of swapped
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
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);
    pts[ind1] = vec1_;
    pts[ind2] = vec2_;

    // Second order and higher derivatives vanish. They are already
    // set to zero, so we return.
    return;

}


//===========================================================================
void Plane::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    n = normal_;
    if (isSwapped())
        n *= -1.0;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Plane::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented.");
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
    MESSAGE("nextSegmentVal() doesn't make sense. Returning arbitrarily 0.0.");
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
    getOrientedParameters(clo_u, clo_v);
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
    if (!isBounded())
        THROW("Unbounded plane - no closest boundary point exist.");

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
  double best_dist = 100000.0; // "Large" number
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

  getOrientedParameters(clo_u, clo_v);
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
    MESSAGE("Computing distance to projected "
            "- not necessarily closest - point.");
    return pnt.dist(projectPoint(pnt));
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
    if (from_upar >= to_upar )
        THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

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

    RectDomain dom = parameterDomain();
    double umin = dom.umin();
    double umax = dom.umax();
    double vmin = dom.vmin();
    double vmax = dom.vmax();

    // Handle the case if not bounded
    if (!isBounded()) {
        double max = 1.0e8; // "Large" number...
        double inf = numeric_limits<double>::infinity();
        if (umin == -inf)
            umin = -max;
        if (umax == inf)
            umax = max;
        if (vmin == -inf)
            vmin = -max;
        if (vmax == inf)
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

    SplineSurface* surf = new SplineSurface(ncoefsu, ncoefsv, ordu, ordv,
                             knotsu.begin(), knotsv.begin(), 
                             coefs.begin(), dim);
    return surf;
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
bool Plane::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    closed_dir_u = false;
    closed_dir_v = false;
    return false;
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
        RectDomain domain = parameterDomain();
        double global_u_low = std::max(clo_u_low, domain.umin());
        double global_v_low = std::max(clo_v_low, domain.vmin());
        double global_u_high = std::min(clo_u_high, domain.umax());
        double global_v_high = std::min(clo_v_high, domain.vmax());
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
  if (isSwapped())
      normal *= -1.0;

  // This surface is a plane
  return true;
}

//===========================================================================
bool Plane::isLinear(Point& dir1, Point& dir2, double tol)
//===========================================================================
{
  if (isSwapped())
    {
      dir1 = vec2_;
      dir2 = vec1_;
    }
  else
    {
      dir1 = vec1_;
      dir2 = vec2_;
    }
  return true;
}

} // namespace Go
