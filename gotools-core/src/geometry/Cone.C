/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
using std::streamsize;
using std::swap;


namespace Go
{


// Constructor.
//===========================================================================
Cone::Cone(double radius,
           Point location, Point z_axis, Point x_axis,
           double cone_angle, bool isSwapped)
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
    double inf = numeric_limits<double>::infinity();
    setParameterBounds(0.0, -inf, 2.0 * M_PI, inf);

    if (isSwapped)
        swapParameterDirection();
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
    // NB: Parameter sequence in the g2 file format is different
    // than the argument list of the setParameterBounds() function!

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

    // "Reset" swapping
    isSwapped_ = false;

    int isBounded; 
    is >> isBounded;
    if (isBounded == 0) {
        // Unbounded in v direction

        // NB: See comment on parameter sequence above!
        double from_upar, to_upar;
        is >> from_upar >> to_upar;

        // Need to take care of rounding errors: If upars are "roughly"
        // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
        const double pareps = 1.0e-4; // This is admittedly arbitrary...
        if (fabs(from_upar) < pareps && fabs(to_upar - 2.0*M_PI) < pareps) {
            from_upar = 0.0;
            to_upar = 2.0 * M_PI;
        }
        double inf = numeric_limits<double>::infinity();
        setParameterBounds(from_upar, -inf, to_upar, inf);
    }
    else if (isBounded == 1) {
        // NB: See comment on parameter sequence above!
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
void Cone::write(std::ostream& os) const
//===========================================================================
{
    // NB: Parameter sequence in the g2 file format is different
    // than the argument list of the setParameterBounds() function!

    streamsize prev = os.precision(15);
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl
       << cone_angle_ << endl;

    // NB: Mind the parameter sequence!
    if (!isBounded()) {
        os << "0" << endl
           << domain_.umin() << " " << domain_.umax() << endl;
    }
    else {
        os << "1" << endl
           << domain_.umin() << " " << domain_.umax() << endl
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
    SplineSurface* tmp = geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
Cone* Cone::clone() const
//===========================================================================
{
    Cone* cone = new Cone(radius_, location_, z_axis_, x_axis_,
        cone_angle_, isSwapped_);
    cone->domain_ = domain_;
    return cone;
}

//===========================================================================
const RectDomain& Cone::parameterDomain() const
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
Cone::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("allBoundaryLoops() not implemented. Returns an empty vector.");
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
    double u = 0.5*(umin+umax);
    double v = 0.0;
    if (isSwapped_)
        swap(u, v);
    normal(dir, u, v);

    return DirectionCone(dir, 0.5*(umax-umin));
}


//===========================================================================
DirectionCone Cone::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    if (isSwapped())
        pardir_is_u = !pardir_is_u;

    DirectionCone normals = normalCone();
    Point dir = normals.centre();
    if (pardir_is_u) {
        if (isSwapped())
            dir *= -1.0;
        Point tandir = z_axis_.cross(dir);
        double angle = normals.angle();
        return DirectionCone(tandir, angle);
    }
    else {
        Point tmpdir = z_axis_.cross(dir);
        Point tandir = dir.cross(tmpdir);
        double angle = normals.angle();
        return DirectionCone(tandir, angle);
    }
}


//===========================================================================
void Cone::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar); // In case of swapped
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

    // Swap parameters, if needed
    getOrientedParameters(upar, vpar);
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First derivatives. TESTME
    pts[ind1] = (radius_ + vpar * tan(cone_angle_)) 
        * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[ind2] = tan(cone_angle_) * (cos(upar) * x_axis_ 
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
    getOrientedParameters(upar, vpar);
    double tana = tan(cone_angle_);
    double tana2 = tana * tana;
    n = (cos(upar) * x_axis_ + sin(upar) * y_axis_ - tana * z_axis_)
        / sqrt(1.0 + tana2);
    if (isSwapped())
        n *= -1.0;

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
  //MESSAGE("constParamCurves() not yet implemented");

    // If domain is unbounded int the const par dir there is nothing we can do.

    bool cone_pardir_is_u = (isSwapped()) ? !pardir_is_u : pardir_is_u;
    if (isSwapped())
    {
        MESSAGE("Not yet tested this function with swapped cone!");
    }
    vector<shared_ptr<ParamCurve> > res;
    if (cone_pardir_is_u)
    {
        shared_ptr<ParamCurve> circle = getCircle(parameter);
        res.push_back(circle);
    }
    else
    {
        if (!isBounded())
        {
            MESSAGE("constParamCurves() not supported for unbounded cone in linear direction!");
        }
        else
        {
            double vmin = domain_.vmin();
            double vmax = domain_.vmax();
            Point cv_min = ParamSurface::point(parameter, vmin);
            Point cv_max = ParamSurface::point(parameter, vmax);
            shared_ptr<Line> line(new Line(cv_min, cv_max, vmin, vmax));
            res.push_back(line);
        }
    }

    return res;
}


//===========================================================================
shared_ptr<ParamCurve>
Cone::constParamCurve(double iso_par, bool pardir_is_u,
			  double from, double to) const
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > res;
    bool real_pardir_is_u = (isSwapped()) ? !pardir_is_u : pardir_is_u;
    if (real_pardir_is_u)
    {
	shared_ptr<ParamCurve> circle = getCircle(iso_par);
	shared_ptr<ParamCurve> sub_circle(circle->subCurve(from, to));
	return sub_circle;
    }
    else
    {
	Point par_from(iso_par, from);
	getOrientedParameters(par_from[0], par_from[1]);
	Point par_to(iso_par, to);
	getOrientedParameters(par_to[0], par_to[1]);

	Point cv_min = ParamSurface::point(par_from[0], par_from[1]);
	Point cv_max = ParamSurface::point(par_to[0], par_to[1]);
	shared_ptr<ParamCurve> line(new Line(cv_min, cv_max, from, to));
	return line;
    }
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
    MESSAGE("nextSegmentVal() doesn't make sense. Returning arbitrarily 0.0.");
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
        if (isSwapped_) {
            umin = std::max(umin, domain_of_interest->vmin());
            umax = std::min(umax, domain_of_interest->vmax());
            vmin = std::max(vmin, domain_of_interest->umin());
            vmax = std::min(vmax, domain_of_interest->umax());
        } else {
            umin = std::max(umin, domain_of_interest->umin());
            umax = std::min(umax, domain_of_interest->umax());
            vmin = std::max(vmin, domain_of_interest->vmin());
            vmax = std::min(vmax, domain_of_interest->vmax());
        }
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
        getOrientedParameters(clo_u, clo_v);
        return;
    }
    if (vmin > vvalmax) {
        rad = radius_ + vmin * tan(cone_angle_);
        loc = location_ + vmin * z_axis_;
        circle = Circle(rad, loc, z_axis_, x_axis_);
        circle.setParamBounds(umin, umax);
        circle.closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
        clo_v = vmin;
        getOrientedParameters(clo_u, clo_v);
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
    if (clo_dist < epsilon) {
        getOrientedParameters(clo_u, clo_v);
        return;
    }

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
        if (clo_dist < epsilon) {
            getOrientedParameters(clo_u, clo_v);
            return;
        }
    }

    // Are there more edges?
    if (fabs(umax - umin - 2.0 * M_PI) < epsilon) {
        getOrientedParameters(clo_u, clo_v);
        return;
    }

    // Left - a line
    line = getLine(umin);
    line->closestPoint(pt, vmin, vmax, tmp_clo_v, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_u = umin;
    if (tmp_clo_dist < clo_dist) {
        clo_u = tmp_clo_u;
        clo_v = tmp_clo_v;
        clo_pt = tmp_clo_pt;
        clo_dist = tmp_clo_dist;
        if (clo_dist < epsilon) {
            getOrientedParameters(clo_u, clo_v);
            return;
        }
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

    getOrientedParameters(clo_u, clo_v);
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
    //// In progress... @jbt
    //double inf = numeric_limits<double>::infinity();
    //const double large_number = 1.0e8;

    //double umin = domain_.umin();
    //double umax = domain_.umax();
    //double vmin = domain_.vmin();
    //if (vmin == -inf)
    //    vmin = -large_number;
    //double vmax = domain_.vmax();
    //if (vmax == inf)
    //    vmax = large_number;

    //// Checking the four bounding curves.

    //double clo_u_tmp, clo_v_tmp, clo_dist_tmp;
    //Point clo_pt_tmp;

    //// First
    //shared_ptr<Line> line = getLine(umin);
    //clo_u = umin;
    //line->closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);
    //// Second
    //line = getLine(umax);
    //line->closestPoint(pt, vmin, vmax, clo_v_tmp, clo_pt_tmp, clo_dist_tmp);
    //if (clo_dist_tmp < clo_dist) {
    //    clo_u = clo_u_tmp;
    //    clo_pt = clo_pt_tmp;
    //    clo_dist = clo_dist_tmp;
    //}
    //// Third
    //shared_ptr<Circle> circle = getCircle(vmin);


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
bool Cone::isDegenerate(bool& b, bool& r,
                        bool& t, bool& l, double tolerance) const
//===========================================================================
{
    bool res = false;
    b = false;
    r = false;
    t = false;
    l = false;
    double kmin = radius_ + domain_.vmin() * tan(cone_angle_);
    if (kmin == 0.0) {
        b = true;
        res = true;
        if (isSwapped())
            swap(b, l);
    }
    double kmax = radius_ + domain_.vmax() * tan(cone_angle_);
    if (kmax == 0.0) {
        t = true;
        res = true;
        if (isSwapped())
            swap(t, r);
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

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    double tol = 1.0e-13;
    if (from_upar > -2.0 * M_PI - tol && from_upar < -2.0 * M_PI)
      from_upar = -2.0 * M_PI;
    if (to_upar < 2.0 * M_PI + tol && to_upar >2.0 * M_PI)
      to_upar = 2.0 * M_PI;
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
        THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
        THROW("(to_upar - from_upar) must not exceed 2pi.");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Cone::setParamBoundsU(double from_upar, double to_upar)
//===========================================================================
{
    RectDomain tmp_domain = parameterDomain();
    double from_vpar = tmp_domain.vmin();
    double to_vpar = tmp_domain.vmax();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    double tol = 1.0e-13;
    if (from_upar > -2.0 * M_PI - tol && from_upar < -2.0 * M_PI)
      from_upar = -2.0 * M_PI;
    if (to_upar < 2.0 * M_PI + tol && to_upar >2.0 * M_PI)
      to_upar = 2.0 * M_PI;
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
void Cone::setParamBoundsV(double from_vpar, double to_vpar)
//===========================================================================
{
    RectDomain tmp_domain = parameterDomain();
    double from_upar = tmp_domain.umin();
    double to_upar = tmp_domain.umax();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    Array<double, 2> ll(domain_.umin(), from_vpar);
    Array<double, 2> ur(domain_.umax(), to_vpar);
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
shared_ptr<Circle> Cone::getCircle(double par) const
//===========================================================================
{
    Point centre = location_ + par * z_axis_;
    shared_ptr<Circle> circle(new Circle(radius_, centre, z_axis_, x_axis_));
    // Note: We are using domain_ on purpose, because domain_'s
    // u-direction is always the angular direction, no matter what
    // isSwapped_ is.
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->setParamBounds(umin, umax);
    return circle;
}

    
//===========================================================================
bool Cone::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    closed_dir_u = (domain_.umax() - domain_.umin() == 2.0*M_PI);
    closed_dir_v = false;
    if (isSwapped())
        swap(closed_dir_u, closed_dir_v);
    return (closed_dir_u || closed_dir_v);
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
    double tmpu = umax - umin;
    double tmpv = 0.0;
    getOrientedParameters(tmpu, tmpv);
    point(pt, tmpu, tmpv);
    double tmpdist;
    Circle circle(radius_, location_, z_axis_, x_axis_);
    SplineCurve* scircle = circle.geometryCurve();
    scircle->closestPoint(pt, 0.0, 2.0 * M_PI, tmpu, tmppt, tmpdist);
    double epsilon = 1.0e-10;
    if (tmpu < epsilon && umax - umin == 2.0 * M_PI) {
        tmpu = 2.0 * M_PI;
    }
    SplineSurface* subpatch = surface.subSurface(0.0, vmin, tmpu, vmax);
    subpatch->basis_u().rescale(umin, umax);
    GeometryTools::translateSplineSurf(-location_, *subpatch);
    GeometryTools::rotateSplineSurf(z_axis_, umin, *subpatch);
    GeometryTools::translateSplineSurf(location_, *subpatch);
    delete scircle;

    if (isSwapped())
        subpatch->swapParameterDirection();

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
shared_ptr<ElementaryCurve> 
Cone::getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
			      const Point* start_par_pt, const Point* end_par_pt) const 
//===========================================================================
{
  // Default is not simple elementary parameter curve exists
  shared_ptr<ElementaryCurve> dummy;
  
  // Bookkeeping related to swapped parameters
  int ind1 = 0;
  int ind2 = 1;
  if (isSwapped())
      swap(ind1, ind2);

  double t1, t2;
  int idx;
  bool closed = false;
  if (space_crv->instanceType() == Class_Line)
    {
      if (!((Line*)(space_crv))->isBounded())
	return dummy;   // Project endpoints onto the surface
      t1 = space_crv->startparam();
      t2 = space_crv->endparam();
      idx = ind1;
    }
  else if (space_crv->instanceType() == Class_Circle)
    {
      t1 = space_crv->startparam();
      closed = ((Circle*)(space_crv))->isClosed();
      t2 = (closed) ? 0.5*(t1 + space_crv->endparam()) :
	space_crv->endparam();
      idx = ind2;
    }
  else
    return dummy;

  Point par1(2), par2(2);
  if ((start_par_pt != NULL) && (end_par_pt != NULL))
  {
      par1 = *start_par_pt;
      par2 = *end_par_pt;
  }
  else
  {
      double parval1[2], parval2[2];
      double d1, d2;
      Point close1, close2;
      Point pos1 = space_crv->ParamCurve::point(t1);
      Point pos2 = space_crv->ParamCurve::point(t2);
      closestPoint(pos1, parval1[0], parval1[1], close1, d1, tol);
      closestPoint(pos2, parval2[0], parval2[1], close2, d2, tol);
      if (d1 > tol || d2 > tol)
          return dummy;

      par1[idx] = par2[idx] = 0.5*(parval1[idx] + parval2[idx]);
      par1[1-idx] = parval1[1-idx];
      par2[1-idx] = parval2[1-idx];
      Point mid = this->ParamSurface::point(0.5*(par1[0]+par2[0]), 0.5*(par1[1]+par2[1]));
      Point cv_mid = space_crv->ParamCurve::point(0.5*(t1+t2));
      if (mid.dist(cv_mid) > tol)
      {
          bool dummy_u, dummy_v;
          if (isClosed(dummy_u, dummy_v))
          {
              // Extra check at the seam
              double ptol = 1.0e-4;
              if (parval1[ind1] < ptol)
                  parval1[ind1] = 2.0*M_PI;
              else if (par1[ind1] > 2.0*M_PI - ptol)
                  parval1[ind1] = 0.0;
              else if (par2[ind1] < ptol)
                  parval2[ind1] = 2.0*M_PI;
              else if (par2[ind1] > 2.0*M_PI - ptol)
                  parval2[ind1] = 0.0;
              par1[idx] = par2[idx] = 0.5*(parval1[idx] + parval2[idx]);
              par1[1-idx] = parval1[1-idx];
              par2[1-idx] = parval2[1-idx];
              mid = this->ParamSurface::point(0.5*(par1[0]+par2[0]), 0.5*(par1[1]+par2[1]));
              if (mid.dist(cv_mid) > tol)
                  return dummy;
          }
          else
              return dummy;  // Linear parameter curve not close enough
      }

      if (closed)
          par2[ind1] = par1[ind1] + 2.0*M_PI;

      if (start_par_pt != NULL)
      {
          // MESSAGE("Avoid computing par1.");
          par1 = *start_par_pt;
      }
      if (end_par_pt != NULL)
      {
          // MESSAGE("Avoid computing par2.");
          par2 = *end_par_pt;
      }
  }

  shared_ptr<Line> param_cv(new Line(par1, par2, 
				     space_crv->startparam(), space_crv->endparam()));

  // TEST
  Point p1 = param_cv->ParamCurve::point(param_cv->startparam());
  Point p2 = param_cv->ParamCurve::point(param_cv->endparam());
  
  return param_cv;
}


//===========================================================================
bool Cone::isAxisRotational(Point& centre, Point& axis, Point& vec,
				double& angle)
//===========================================================================
{
  centre = location_;
  axis = z_axis_;
  if (domain_.umin() == 0.0)
    vec = x_axis_;
  else
    {
      Point pt;
      RectDomain domain = parameterDomain();
      point(pt, domain.umin(), domain.vmin());
      vec = pt - location_;
      vec.normalize();
    }
  angle = domain_.umax() - domain_.umin();

  return true;
}

} // namespace Go
