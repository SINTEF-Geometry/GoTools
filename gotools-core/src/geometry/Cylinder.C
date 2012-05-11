//==========================================================================
//                                                                          
// File: Cylinder.C                                                          
//                                                                          
// Created: Mon Oct 20 14:18:14 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Cylinder.C,v 1.10 2009-03-04 15:41:34 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Line.h"
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
Cylinder::Cylinder(double radius,
                   Point location, Point z_axis, Point x_axis)
    : radius_(radius),
      location_(location), z_axis_(z_axis), x_axis_(x_axis)
//===========================================================================
{
    // The use of the z- and x-axis in this constructor - and in the
    // definition of a cylinder - is based on the use of the
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
Cylinder::~Cylinder()
//===========================================================================
{
}

//===========================================================================
void Cylinder::read (std::istream& is)
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

    // Must be put inside a BoundedSurface if cylinder is to be
    // bounded.
    Array<double, 2> ll(0.0, -numeric_limits<double>::infinity());
    Array<double, 2> ur(2.0 * M_PI, numeric_limits<double>::infinity());
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Cylinder::write(std::ostream& os) const
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
int Cylinder::dimension() const
//===========================================================================
{
    return location_.dimension();
}

//===========================================================================
ClassType Cylinder::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
BoundingBox Cylinder::boundingBox() const
//===========================================================================
{
    // First handle the case if not bounded
    if (!isBounded()) {
        // Does not make sense - return a default constructed
        // BoundingBox, which is invalid. @jbt
        return BoundingBox();
    }

    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    vector<Point> points;

    // Find bounding box of circle at vmin
    shared_ptr<Circle> vmin_circle = getCircle(vmin);
    BoundingBox vmin_bbox = vmin_circle->boundingBox();
    points.push_back(vmin_bbox.low());
    points.push_back(vmin_bbox.high());

    // Find bounding box of circle at vmax
    shared_ptr<Circle> vmax_circle = getCircle(vmax);
    BoundingBox vmax_bbox = vmax_circle->boundingBox();
    points.push_back(vmax_bbox.low());
    points.push_back(vmax_bbox.high());

    BoundingBox box;
    box.setFromPoints(points);
    return box;
}

//===========================================================================
const Domain& Cylinder::parameterDomain() const
//===========================================================================
{
    return domain_;
}


//===========================================================================
std::vector<CurveLoop> 
Cylinder::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Cylinder::normalCone() const
//===========================================================================
{
    double umin = domain_.umin();
    double umax = domain_.umax();
    Point dir;
    normal(dir, 0.5*(umin+umax), 0.0);

    return DirectionCone(dir, 0.5*(umax-umin));
}


//===========================================================================
DirectionCone Cylinder::tangentCone(bool pardir_is_u) const
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
void Cylinder::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    pt = location_ 
        + radius_ * (cos(upar) * x_axis_ + sin(upar) * y_axis_)
        + vpar * z_axis_;
}


//===========================================================================
void Cylinder::point(std::vector<Point>& pts, 
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
    pts[1] = radius_ * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[2] = z_axis_;
    if (derivs == 1)
        return;

    // Second order and higher derivatives.
    for (int i = 2; i <= derivs; ++i) {
        int index = i*(i+1)/2;
        pts[index] = radius_ * (cos(upar + i*0.5*M_PI) * x_axis_
                                + sin(upar + i*0.5*M_PI) * y_axis_);
    }

}


//===========================================================================
void Cylinder::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    n = cos(upar) * x_axis_ + sin(upar) * y_axis_;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Cylinder::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
Cylinder* Cylinder::subSurface(double from_upar, double from_vpar,
                               double to_upar, double to_vpar,
                               double fuzzy) const
//===========================================================================
{
    Cylinder* cylinder = clone();
    cylinder->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return cylinder;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
Cylinder::subSurfaces(double from_upar, double from_vpar,
                      double to_upar, double to_vpar,
                      double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Cylinder> cylinder(subSurface(from_upar, from_vpar,
                                             to_upar, to_vpar));
    res.push_back(cylinder);
    return res;
}


//===========================================================================
double 
Cylinder::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    MESSAGE("Does not make sense. Return arbitrarily zero.");
    return 0.0;
}


//===========================================================================
void Cylinder::closestPoint(const Point& pt,
                            double&        clo_u,
                            double&        clo_v, 
                            Point&         clo_pt,
                            double&        clo_dist,
                            double         epsilon,
                            const RectDomain* domain_of_interest,
                            double   *seed) const
//===========================================================================
{
    // Find closest point on central axis
    Line axis(location_, z_axis_);
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    axis.closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);

    // Find closest point on the circle at the found v parameter
    shared_ptr<Circle> circle = getCircle(clo_v);
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);

    // That's it - we have what we need...
}


//===========================================================================
void Cylinder::closestBoundaryPoint(const Point& pt,
                                    double&        clo_u,
                                    double&        clo_v, 
                                    Point&         clo_pt,
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
void Cylinder::getBoundaryInfo(Point& pt1, Point& pt2,
                               double epsilon, SplineCurve*& cv,
                               SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
void Cylinder::turnOrientation()
//===========================================================================
{
    MESSAGE("'Turn orientation' is ambigous - did you \n"
            "mean 'reverse parameter direction n'?");
}



//===========================================================================
void Cylinder::swapParameterDirection()
//===========================================================================
{
    MESSAGE("Does not make sense for Cylinders.");
}


//===========================================================================
void Cylinder::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    // This operation does not make sense, because the coordinate
    // system is always right-handed.

    MESSAGE("reverseParameterDirection() not implemented.");
}


//===========================================================================
bool Cylinder::isDegenerate(bool& b, bool& r,
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
void Cylinder::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Not this kind of degeneracy
    return;
}

//===========================================================================
void Cylinder::setCoordinateAxes()
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
void Cylinder::setParameterBounds(double from_upar, double from_vpar,
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
void Cylinder::setParamBoundsU(double from_upar, double to_upar)
//===========================================================================
{
    if (from_upar >= to_upar )
        THROW("First u-parameter must be strictly less than second.");
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
        THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
        THROW("(to_upar - from_upar) must not exceed 2pi.");

    Array<double, 2> ll(from_upar, domain_.vmin());
    Array<double, 2> ur(to_upar, domain_.vmax());
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Cylinder::setParamBoundsV(double from_vpar, double to_vpar)
//===========================================================================
{
    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    Array<double, 2> ll(domain_.umin(), from_vpar);
    Array<double, 2> ur(domain_.umax(), to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
bool Cylinder::isBounded() const
//===========================================================================
{
    // It is enough to check the v-direction, since the u-direction is
    // always bounded.

    return domain_.vmin() > -numeric_limits<double>::infinity() &&
        domain_.vmax() < numeric_limits<double>::infinity();

}


//===========================================================================
SplineSurface* Cylinder::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Cylinder::createSplineSurface() const
//===========================================================================
{
    // Based on SISL functions s1021 and s1022.

    // First handle the case if not bounded
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    if (!isBounded()) {
        double max = 1.0e6;//8; // "Large" number...
        if (vmin == -numeric_limits<double>::infinity())
            vmin = -max;
        if (vmax == numeric_limits<double>::infinity())
            vmax = max;
    }

    // Knot vector around the cylinder.
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
    // Knot vector along the cylinder.
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
    Point axis1 = radius_ * x_axis_;
    Point axis2 = radius_ * y_axis_;
    for (int ki = 0; ki < 3; ki++){
        rcoef[ki] = bottom_pos[ki] + axis1[ki];
        rcoef[4 + ki] = weight*(bottom_pos[ki] + axis1[ki] + axis2[ki]);
        rcoef[8 + ki] = bottom_pos[ki] + axis2[ki];
        rcoef[12 + ki] = weight*(bottom_pos[ki] - axis1[ki] + axis2[ki]);
        rcoef[16 + ki] = bottom_pos[ki] - axis1[ki];
        rcoef[20 + ki] = weight*(bottom_pos[ki] - axis1[ki] - axis2[ki]);
        rcoef[24 + ki] = bottom_pos[ki] - axis2[ki];
        rcoef[28 + ki] = weight*(bottom_pos[ki] + axis1[ki] - axis2[ki]);
        rcoef[32 + ki] = rcoef[ki];

        rcoef[36 + ki] = top_pos[ki] + axis1[ki];
        rcoef[40 + ki] = weight*(top_pos[ki] + axis1[ki] + axis2[ki]);
        rcoef[44 + ki] = top_pos[ki] + axis2[ki];
        rcoef[48 + ki] = weight*(top_pos[ki] - axis1[ki] + axis2[ki]);
        rcoef[52 + ki] = top_pos[ki] - axis1[ki];
        rcoef[56 + ki] = weight*(top_pos[ki] - axis1[ki] - axis2[ki]);
        rcoef[60 + ki] = top_pos[ki] - axis2[ki];
        rcoef[64 + ki] = weight*(top_pos[ki] + axis1[ki] - axis2[ki]);
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
    // arc-length parametrized in the u-direction.
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmid = 0.5 * (vmin + vmax);
    Point pt, tmppt;
    point(pt, umax - umin, vmid);
    double tmpu, tmpv, tmpdist;
    double epsilon = 1.0e-10;
    surface.closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon);
    if (tmpu < epsilon && umax - umin == 2.0 * M_PI) {
        tmpu = 2.0 * M_PI;
    }
    SplineSurface* subpatch = surface.subSurface(0.0, vmin, tmpu, vmax);
    subpatch->basis_u().rescale(umin, umax);
    GeometryTools::translateSplineSurf(-location_, *subpatch);
    GeometryTools::rotateSplineSurf(z_axis_, umin, *subpatch);
    GeometryTools::translateSplineSurf(location_, *subpatch);

    return subpatch;
}


//===========================================================================
shared_ptr<Circle> Cylinder::getCircle(double vpar) const
//===========================================================================
{
    Point centre = location_ + vpar * z_axis_;
    shared_ptr<Circle> circle(new Circle(radius_, centre, z_axis_, x_axis_));
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->setParamBounds(umin, umax);
    return circle;
}


//===========================================================================


} // namespace Go
