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
    setParameterDomain(0.0, 2.0 * M_PI, -inf, inf);

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
    bool has_param_int = (isBounded >= 10);
    isBounded = isBounded % 10;
    if (isBounded == 0) {
        // Unbounded in v direction

        // NB: See comment on parameter sequence above!
        double from_upar, to_upar;
        is >> from_upar >> to_upar;
	double start_u = from_upar, end_u = to_upar;
	if (has_param_int)
	  {
	    is >> start_u >> end_u;
	  }

        // Need to take care of rounding errors: If upars are "roughly"
        // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
        if (fabs(from_upar) < ptol_ && fabs(to_upar - 2.0*M_PI) < ptol_) {
            from_upar = 0.0;
            to_upar = 2.0 * M_PI;
        }
        double inf = numeric_limits<double>::infinity();
        setParameterBounds(from_upar, -inf, to_upar, inf);
	setParameterDomain(start_u, end_u, -inf, inf);
    }
    else if (isBounded == 1) {
        // NB: See comment on parameter sequence above!
        double from_upar, from_vpar, to_upar, to_vpar;
        is >> from_upar >> to_upar
            >> from_vpar >> to_vpar;
	double start_u = from_upar, end_u = to_upar, start_v = from_vpar, end_v = to_vpar;
	if (has_param_int)
	  {
	    is >> start_u >> end_u >> start_v >> end_v;
	  }

        // Need to take care of rounding errors: If upars are "roughly"
        // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
        if (fabs(from_upar) < ptol_ && fabs(to_upar - 2.0*M_PI) < ptol_) {
            from_upar = 0.0;
            to_upar = 2.0 * M_PI;
        }

        setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
	setParameterDomain(start_u, end_u, start_v, end_v);
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
        os << "10" << endl
           << parbound_.umin() << " " << parbound_.umax() << endl
           << domain_.umin() << " " << domain_.umax() << endl;
    }
    else {
        os << "11" << endl
           << parbound_.umin() << " " << parbound_.umax() << endl
           << parbound_.vmin() << " " << parbound_.vmax() << endl
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
    cone->parbound_ = parbound_;
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
DirectionCone Cone::normalCone() const
//===========================================================================
{
    double umin = parbound_.umin();
    double umax = parbound_.umax();
    Point dir;
    double u = 0.5*(domain_.umin()+domain_.umax());
    double v = (isBounded()) ? domain_.vmin() : 0.0;
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
	  dir *= -1.0;   // This statement lacks in cylinder, what is correct?
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
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    if (isBounded())
      vpar = parbound_.vmin() + 
	(vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
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
    double fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double fac2 = 1.0;
    upar = parbound_.umin() + fac1*(upar-domain_.umin());
    if (isBounded())
      {
	fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	vpar = parbound_.vmin() + fac2*(vpar-domain_.vmin());
      }

    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First order derivatives
    pts[ind1] = fac1*(radius_ + vpar * tan(cone_angle_)) 
      * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[ind2] = fac2*(tan(cone_angle_) * (cos(upar) * x_axis_ 
					  + sin(upar) * y_axis_)  
		      + z_axis_);

    // Second order derivatives
    if (derivs > 1)
      {
	ind1 = 3;
	ind2 = 5;
	if (isSwapped())
	  swap(ind1, ind2);
	pts[ind1] = fac1*fac1*(radius_ + vpar * tan(cone_angle_)) 
	  * (-cos(upar) * x_axis_ - sin(upar) * y_axis_);
	pts[4] = fac1*fac2*tan(cone_angle_)*(-sin(upar)*x_axis_ + cos(upar)*y_axis_);
	pts[ind2].setValue(0.0);
      }

    if (derivs <= 2)
        return;

    // Third order and higher derivatives.
    MESSAGE("Third order or higher derivatives not yet implemented.");

}


//===========================================================================
void Cone::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar);
    double fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double fac2 = 1.0;
	upar = parbound_.umin() + fac1*(upar-domain_.umin());
    if (isBounded())
      {
	fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	vpar = parbound_.vmin() + fac2*(vpar-domain_.vmin());
      }
    double tana = tan(cone_angle_);
    double tana2 = tana * tana;
    n = fac1*fac2*(cos(upar) * x_axis_ + sin(upar) * y_axis_ - tana * z_axis_)
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
    bool cone_pardir_is_u = (isSwapped()) ? !pardir_is_u : pardir_is_u;
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
            shared_ptr<Line> line = getLine(parameter);
            res.push_back(line);
        }
        else
        {
            double vmin = domain_.vmin();
            Point par_from(parameter, vmin);
            getOrientedParameters(par_from[0], par_from[1]);
            Point cv_min = ParamSurface::point(par_from[0], par_from[1]);
            double vmax = domain_.vmax();
            Point par_to(parameter, vmax);
            getOrientedParameters(par_to[0], par_to[1]);
            Point cv_max = ParamSurface::point(par_to[0], par_to[1]);
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
    double fac1 = 
      (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double fac2 = 1.0;
    if (isBounded())
      fac2 =
	(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    if (isSwapped())
      {
	double bound3 = from_upar, bound4 = to_upar;
    	double bound1 = parbound_.umin() + fac1*(from_vpar-domain_.umin());
    	double bound2 = parbound_.umin() + fac1*(to_vpar-domain_.umin());
	if (isBounded())
	  {
	    bound3 = parbound_.vmin() + fac2*(from_upar-domain_.vmin());
	    bound4 = parbound_.vmin() + fac2*(to_upar-domain_.vmin());
	  }
    	cone->setParameterBounds(bound3, bound1, bound4, bound2);
      }
    else
      {
	double bound3 = from_vpar, bound4 = to_vpar;
	double bound1 = parbound_.umin() + fac1*(from_upar-domain_.umin());
	double bound2 = parbound_.umin() + fac1*(to_upar-domain_.umin());
	if (isBounded())
	  {
	    bound3 = parbound_.vmin() + fac2*(from_vpar-domain_.vmin());
	    bound4 = parbound_.vmin() + fac2*(to_vpar-domain_.vmin());
	  }
	cone->setParameterBounds(bound1, bound3, bound2, bound4);
      }
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

    // Operate in the circle and line parameterization of the cone
    double fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double fac2 = 1.0;
    umin = parbound_.umin() + fac1*(umin-domain_.umin());
    umax = parbound_.umin() + fac1*(umax-domain_.umin());
    if (isBounded())
      {
	fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	vmin = parbound_.vmin() + fac2*(vmin-domain_.vmin());
	vmax = parbound_.vmin() + fac2*(vmax-domain_.vmin());
      }

    // Identify the two values of the v-parameter where an unbounded
    // cone is orthogonal to the cone-to-point vector.
    double rad = radius_;
    if (radius_ < epsilon)
        rad = 1.0;
    Point loc = location_;
    Circle circle(rad, loc, z_axis_, x_axis_);
    circle.closestPoint(pt, 0.0, 2.0*M_PI, clo_u, clo_pt, clo_dist);

    Point cossin = cos(clo_u) * x_axis_ + sin(clo_u) * y_axis_;
    loc = location_ + radius_ * cossin;
    Point dir = tan(cone_angle_) * cossin + z_axis_;
    shared_ptr<Line> line(new Line(loc, dir));
    double vvalmin, vvalmax, tmp;
    line->closestPoint(pt, vmin, vmax, vvalmin, clo_pt, clo_dist);
    double clo_u2 = clo_u - M_PI;
    if (clo_u2 < 0.0)
    {
        clo_u2 += 2.0*M_PI;
    }

    cossin = cos(clo_u2) * x_axis_ + sin(clo_u2) * y_axis_;
    loc = location_ + radius_ * cossin;
    dir = tan(cone_angle_) * cossin + z_axis_;
    line = shared_ptr<Line>(new Line(loc, dir));
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

	clo_u = domain_.umin() + 
	  (clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
	if (isBounded())
	  clo_v = domain_.vmin() + 
	    (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
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

	clo_u = domain_.umin() + 
	  (clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
	if (isBounded())
	  clo_v = domain_.vmin() + 
	    (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
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
    circle = Circle(fabs(rad), loc, z_axis_, x_axis_);
    circle.setParamBounds(umin, umax);
    circle.closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
    clo_v = vmin;
    if (rad < 0)
      {
	clo_u += M_PI;
	if (clo_u > umax)
	  clo_u -= (2.0*M_PI);
      }
    if (clo_dist < epsilon) {
      clo_u = domain_.umin() + 
	(clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
      if (isBounded())
	clo_v = domain_.vmin() + 
	  (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
      getOrientedParameters(clo_u, clo_v);
      return;
    }

    // Top - a circle
    rad = radius_ + vmax * tan(cone_angle_);
    loc = location_ + vmax * z_axis_;
    circle = Circle(fabs(rad), loc, z_axis_, x_axis_);
    circle.setParamBounds(umin, umax);
    circle.closestPoint(pt, umin, umax, tmp_clo_u, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_v = vmax;
    if (rad < 0)
      {
	tmp_clo_u += M_PI;
	if (tmp_clo_u > umax)
	  tmp_clo_u -= (2.0*M_PI);
      }
    if (tmp_clo_dist < clo_dist) {
        clo_u = tmp_clo_u;
        clo_v = tmp_clo_v;
        clo_pt = tmp_clo_pt;
        clo_dist = tmp_clo_dist;
        if (clo_dist < epsilon) {
	  clo_u = domain_.umin() + 
	    (clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
	  if (isBounded())
	    clo_v = domain_.vmin() + 
	      (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
             getOrientedParameters(clo_u, clo_v);
            return;
        }
    }

    // Are there more edges?
    if (fabs(umax - umin - 2.0 * M_PI) < epsilon) {
      clo_u = domain_.umin() + 
	(clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
      if (isBounded())
	clo_v = domain_.vmin() + 
	  (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
      getOrientedParameters(clo_u, clo_v);
      return;
    }

    // Left - a line
    cossin = cos(umin) * x_axis_ + sin(umin) * y_axis_;
    loc = location_ + radius_ * cossin;
    dir = tan(cone_angle_) * cossin + z_axis_;
    line = shared_ptr<Line>(new Line(loc, dir));
    line->closestPoint(pt, vmin, vmax, tmp_clo_v, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_u = umin;
    if (tmp_clo_dist < clo_dist) {
        clo_u = tmp_clo_u;
        clo_v = tmp_clo_v;
        clo_pt = tmp_clo_pt;
        clo_dist = tmp_clo_dist;
        if (clo_dist < epsilon) {
	  clo_u = domain_.umin() + 
	    (clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
	  if (isBounded())
	    clo_v = domain_.vmin() + 
	      (clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
	  getOrientedParameters(clo_u, clo_v);
	  return;
        }
    }

    // Right - a line
    cossin = cos(umax) * x_axis_ + sin(umax) * y_axis_;
    loc = location_ + radius_ * cossin;
    dir = tan(cone_angle_) * cossin + z_axis_;
    line = shared_ptr<Line>(new Line(loc, dir));
    line->closestPoint(pt, vmin, vmax, tmp_clo_v, tmp_clo_pt, tmp_clo_dist);
    tmp_clo_u = umax;
    if (tmp_clo_dist < clo_dist) {
        clo_u = tmp_clo_u;
        clo_v = tmp_clo_v;
        clo_pt = tmp_clo_pt;
        clo_dist = tmp_clo_dist;
    }


    clo_u = domain_.umin() + 
      (clo_u - parbound_.umin())*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
    if (isBounded())
      clo_v = domain_.vmin() + 
	(clo_v - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
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
    double kmin = radius_ + parbound_.vmin() * tan(cone_angle_);
    if (kmin == 0.0) {
        b = true;
        res = true;
        if (isSwapped())
            swap(b, l);
    }
    double kmax = radius_ + parbound_.vmax() * tan(cone_angle_);
    if (kmax == 0.0) {
        t = true;
        res = true;
        if (isSwapped())
            swap(t, r);
    }
    return res;
}

//===========================================================================
void Cone::getDegenerateParam(double& par, int& dir) const
//===========================================================================
{
  double eps = 1.0e-10;
  par = 0.0;
  double ang_tan = tan(cone_angle_);
  if (fabs(ang_tan) < eps)
    dir = 0;
  else
    {
      par = -radius_/ang_tan;
      par = domain_.vmin() +
	(par - parbound_.vmin())*(domain_.vmax()-domain_.vmin())/(parbound_.vmax()-parbound_.vmin());
      dir = (isSwapped()) ? 1 : 2;
    }
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

    bool bounded = isBounded();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar > -2.0 * M_PI - ptol_ && from_upar < -2.0 * M_PI)
      from_upar = -2.0 * M_PI;
    if (to_upar < 2.0 * M_PI + ptol_ && to_upar >2.0 * M_PI)
      to_upar = 2.0 * M_PI;
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
        THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
        THROW("(to_upar - from_upar) must not exceed 2pi.");

    double fac1 = 
      (domain_.umax()-domain_.umin())/(parbound_.umax() - parbound_.umin());
    double fac2 = 1.0;
    double start_u = domain_.umin() + fac1*(from_upar-parbound_.umin());
    double end_u = domain_.umin() + fac1*(to_upar-parbound_.umin());
    double start_v = from_vpar, end_v = to_vpar;
    if (bounded)
      {
      fac2 = (domain_.vmax()-domain_.vmin())/(parbound_.vmax() - parbound_.vmin());
      start_v = domain_.vmin() + fac2*(from_vpar-parbound_.vmin());
      end_v = domain_.vmin() + fac2*(to_vpar-parbound_.vmin());
      }

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(start_u, start_v);
    Array<double, 2> ur2(end_u, end_v);
    domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Cone::setParamBoundsU(double from_upar, double to_upar)
//===========================================================================
{
  RectDomain tmp_domain = parbound_;
  double from_vpar = tmp_domain.vmin();
  double to_vpar = tmp_domain.vmax();
  getOrientedParameters(from_upar, from_vpar);
  getOrientedParameters(to_upar, to_vpar);

  if (from_upar > -2.0 * M_PI - ptol_ && from_upar < -2.0 * M_PI)
    from_upar = -2.0 * M_PI;
  if (to_upar < 2.0 * M_PI + ptol_ && to_upar >2.0 * M_PI)
    to_upar = 2.0 * M_PI;
  if (from_upar >= to_upar )
    THROW("First u-parameter must be strictly less than second.");
  if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
    THROW("u-parameters must be in [-2pi, 2pi].");
  if (to_upar - from_upar > 2.0 * M_PI)
    THROW("(to_upar - from_upar) must not exceed 2pi.");

  double fac1 = 
    (domain_.umax()-domain_.umin())/(parbound_.umax() - parbound_.umin());
  double start_u = domain_.umin() + fac1*(from_upar-parbound_.umin());
  double end_u =  domain_.umin() + fac1*(to_upar-parbound_.umin());

  Array<double, 2> ll(from_upar, from_vpar);
  Array<double, 2> ur(to_upar, to_vpar);
  parbound_ = RectDomain(ll, ur);
   
  Array<double, 2> ll2(start_u, domain_.vmin());
  Array<double, 2> ur2(end_u, domain_.vmax());
  domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Cone::setParamBoundsV(double from_vpar, double to_vpar)
//===========================================================================
{
  RectDomain tmp_domain = parbound_;
    double from_upar = tmp_domain.umin();
    double to_upar = tmp_domain.umax();
    bool bounded = isBounded();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    double fac2 = (bounded) ?
      (domain_.vmax()-domain_.vmin())/(parbound_.vmax() - parbound_.vmin()) : 1.0;
    double start_v = (bounded) ?
      domain_.vmin() + fac2*(from_vpar-parbound_.vmin()) : from_vpar;
    double end_v = (bounded) ?
      domain_.vmin() + fac2*(to_vpar-parbound_.vmin()) : to_vpar;

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(domain_.umin(), start_v);
    Array<double, 2> ur2(domain_.umax(), end_v);
    domain_ = RectDomain(ll2, ur2);
}

//===========================================================================
void Cone::setParameterDomain(double startpar_u, double endpar_u, 
			      double startpar_v, double endpar_v)
//===========================================================================
{
  getOrientedParameters(startpar_u, startpar_v);
  getOrientedParameters(endpar_u, endpar_v);
  Array<double, 2> ll(startpar_u, startpar_v);
  Array<double, 2> ur(endpar_u, endpar_v);
  domain_ = RectDomain(ll, ur);
}

//===========================================================================
void Cone::restrictParameterDomain(double startpar_u, double endpar_u, 
				   double startpar_v, double endpar_v)
//===========================================================================
{
  getOrientedParameters(startpar_u, startpar_v);
  getOrientedParameters(endpar_u, endpar_v);
  startpar_u = parbound_.umin() + 
    (startpar_u-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
  endpar_u = parbound_.umin() + 
    (endpar_u-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
  if (isBounded())
    {
      startpar_v = parbound_.vmin() + 
	(startpar_v-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
      endpar_v = parbound_.vmin() + 
	(endpar_v-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    }

  if (isSwapped())
    {
      std::swap(startpar_u, startpar_v);
      std::swap(endpar_u, endpar_v);
    }
  
  setParameterBounds(startpar_u, startpar_v, endpar_u, endpar_v);
}

//===========================================================================
bool Cone::isBounded() const
//===========================================================================
{
    // It is enough to check the v-direction, since the u-direction is
    // always bounded.

    return parbound_.vmin() > -numeric_limits<double>::infinity() &&
        parbound_.vmax() < numeric_limits<double>::infinity();

}

    
//===========================================================================
shared_ptr<Circle> Cone::getCircle(double par) const
//===========================================================================
{
  if (isBounded())
    par = parbound_.vmin() + 
      (par-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
  Point centre = location_ + par * z_axis_;
  const double radius = radius_ + par * tan(cone_angle_);
  shared_ptr<Circle> circle(new Circle(radius, centre, z_axis_, x_axis_));
  // Note: We are using domain_ on purpose, because domain_'s
  // u-direction is always the angular direction, no matter what
  // isSwapped_ is.
  double umin = parbound_.umin();
  double umax = parbound_.umax();
  circle->setParamBounds(umin, umax);
  circle->setParameterInterval(domain_.umin(), domain_.umax());
  return circle;
}

    
//===========================================================================
bool Cone::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
  closed_dir_u = (parbound_.umax() - parbound_.umin() >= 2.0*M_PI - ptol_ &&
		  parbound_.umax() - parbound_.umin() <= 2.0*M_PI + ptol_);
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
    double vmin = parbound_.vmin();
    double vmax = parbound_.vmax();
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
    double umin = parbound_.umin();
    double umax = parbound_.umax();
    double tmpv = 0.5 * (vmin + vmax);
    Point pt, tmppt;
    double tmpu = umax - umin;
    double rad = radius_ + tmpv * tan(cone_angle_);
    Point loc = location_ + tmpv*z_axis_;
    pt = loc + rad*(cos(tmpu)*x_axis_ + sin(tmpu)*y_axis_);
    double tmpdist;
    Circle circle(rad, loc, z_axis_, x_axis_);
    SplineCurve* scircle = circle.geometryCurve();
    double guesspar = tmpu;
    scircle->closestPoint(pt, 0.0, 2.0 * M_PI, tmpu, tmppt, tmpdist,
			  &guesspar);
    double epsilon = 1.0e-10;
    if (tmpu < epsilon && umax - umin == 2.0 * M_PI) {
        tmpu = 2.0 * M_PI;
    }
    SplineSurface* subpatch = surface.subSurface(0.0, vmin, tmpu, vmax);
    subpatch->basis_u().rescale(domain_.umin(), domain_.umax());
    if (isBounded())
      subpatch->basis_v().rescale(domain_.vmin(), domain_.vmax());
    GeometryTools::translateSplineSurf(-location_, *subpatch);
    GeometryTools::rotateSplineSurf(z_axis_, umin, *subpatch);
    GeometryTools::translateSplineSurf(location_, *subpatch);
    delete scircle;

    if (isSwapped())
        subpatch->swapParameterDirection();

    return subpatch;
}

//===========================================================================
SplineSurface* Cone::createNonRationalSpline(double eps) const
//===========================================================================
{
  if (!isBounded())
    {
      MESSAGE("createNonRationalSpline is not implemented in the unbounded case");
      return NULL;
    }

  // First fetch the circular boundary curves
  shared_ptr<Circle> circ1 = getCircle(domain_.vmin());
  shared_ptr<Circle> circ2 = getCircle(domain_.vmax());

  // Get Spline approximation
  shared_ptr<SplineCurve> crv1(circ1->createNonRationalSpline(eps));
  shared_ptr<SplineCurve> crv2(circ2->createNonRationalSpline(eps));

  // Interpolate curves
  // This is a hack because GoTools lacks proper lofting functionality
  // for surfaces
  double knot_diff_tol = 1.0e-8; 
  vector<shared_ptr<SplineCurve> > bd_cvs(2);
  bd_cvs[0] = crv1;
  bd_cvs[1] = crv2;
  GeometryTools::unifyCurveSplineSpace(bd_cvs, knot_diff_tol);
  vector<double> coefs;
  coefs.insert(coefs.end(), crv1->coefs_begin(), crv1->coefs_end());
  coefs.insert(coefs.end(), crv2->coefs_begin(), crv2->coefs_end());
  vector<double> knots2(4);
  knots2[0] = knots2[1] = domain_.vmin();
  knots2[2] = knots2[3] = domain_.vmax();
  SplineSurface *surf = new SplineSurface(crv1->numCoefs(), 2,
					  crv1->order(), 2,
					  crv1->knotsBegin(),
					  knots2.begin(), 
					  coefs.begin(),
					  crv1->dimension());

  if (isSwapped())
    surf->swapParameterDirection();

  return surf;
}

//===========================================================================
shared_ptr<Line> Cone::getLine(double upar) const 
//===========================================================================
{
    // if (upar < 0.0 || upar > 2.0*M_PI)
    // {
    //     cout << "Swapped domain for Cone, calling getLine(). upar: " << upar << ", umin: " << domain_.umin() <<
    //         ", umax: " << domain_.umax() << ", values should be inside [0, 2*M_PI)." << 
    //         ", vmin: " << domain_.vmin() << ", vmax: " << domain_.vmax() << endl;
    // }
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    Point cossin = cos(upar) * x_axis_ + sin(upar) * y_axis_;
    Point loc = location_ + radius_ * cossin;
    Point dir = tan(cone_angle_) * cossin + z_axis_;
    shared_ptr<Line> line(new Line(loc, dir));
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    line->setParamBounds(parbound_.vmin(), parbound_.vmax());
    line->setParameterInterval(vmin, vmax);
    return line;
}


//===========================================================================
shared_ptr<ElementaryCurve> 
Cone::getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
			      const Point* start_par_pt, const Point* end_par_pt) const 
//===========================================================================
{
  double eps = 1.0e-9;
  
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
    double fac = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
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
    double u1 = parbound_.umin() + fac*(parval1[ind1] - domain_.umin());
    double u2 = parbound_.umin() + fac*(parval2[ind1] - domain_.umin());
    double u3 = parbound_.umin() + fac*(par1[ind1] - domain_.umin());
    double u4 = parbound_.umin() + fac*(par2[ind1] - domain_.umin());
    if (mid.dist(cv_mid) > tol)
      {
	bool dummy_u, dummy_v;
	if (isClosed(dummy_u, dummy_v))
          {
	    // Extra check at the seam
	    double ptol = (tol < 0.001) ? 1.0e-4 : 1.0e-2;  
	    if (u1 < ptol)
	      u1 = 2.0*M_PI;
	    else if (u3 > 2.0*M_PI - ptol)
	      u1 = 0.0;
	    else if (u4 < ptol)
	      u2 = 2.0*M_PI;
	    else if (u4 > 2.0*M_PI - ptol)
	      u2 = 0.0;
	    parval1[ind1] = domain_.umin()+(u1 - parbound_.umin())/fac;
	    parval2[ind1] = domain_.umin()+(u2 - parbound_.umin())/fac;
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

    u3 = parbound_.umin() + fac*(par1[ind1] - domain_.umin());
    u4 = parbound_.umin() + fac*(par2[ind1] - domain_.umin());
    if (closed)
      {
	u4 = u3 + 2.0*M_PI;
	par2[ind1] = domain_.umin()+(u4 - parbound_.umin())/fac;
      }

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

    double domain_max = (idx == 1) ? parbound_.umax() : parbound_.vmax();
    bool isreversed = false;
    if (closed && fabs(domain_max-par1[1-idx]) < eps)
      {
	double domain_min = (idx == 1) ? parbound_.umin() : parbound_.vmin();
	par1[1-idx] = domain_min;
	par2[1-idx] = domain_max;
	isreversed = true;
      }

  shared_ptr<Line> param_cv(new Line(par1, par2, 
				     space_crv->startparam(), space_crv->endparam(),
				     isreversed));

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

//===========================================================================
double Cone::radius(double u, double v) const
//===========================================================================
{
    getOrientedParameters(u, v); // In case of swapped
    if (isBounded())
      v = parbound_.vmin() + 
	(v-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    double rad = radius_ + v * tan(cone_angle_);
    return rad;
}

//===========================================================================
  void Cone::enlarge(double len1, double len2, double len3, double len4)
//===========================================================================
{
  // Distances are given in geometry space, compute corresponding parameter
  // distance in the rotational direction
  double u1, u2, v1, v2;
  double h = radius_/tan(cone_angle_);
  if (isSwapped())
    {
      double alpha1 = len3/radius_;
      double alpha2 = len4/radius_;
      double p1 = len1*cos(cone_angle_);
      double p2 = len2*cos(cone_angle_);
      u1 = parbound_.vmin() - p1;
      u2 = parbound_.vmax() + p2;
      v1 = std::max(parbound_.umin() - alpha1, 
		    std::max(-2.0*M_PI, parbound_.umax()-2.0*M_PI));
      v2 = std::min(parbound_.umax() + alpha2, 
		    std::min(2.0*M_PI, parbound_.umin()+2.0*M_PI));
      
      u1 = std::max(u1, -h);
      if (v2 - v1 > 2.0*M_PI)
	{
	  double vdel = v2 - v1 - 2.0*M_PI;
	  v2 -= 0.5*vdel;
	  v1 += 0.5*vdel;
	}
     }
  else
    {
      double alpha1 = len1/radius_;
      double alpha2 = len2/radius_;
      double p1 = len3*cos(cone_angle_);
      double p2 = len4*cos(cone_angle_);
      u1 = std::max(parbound_.umin() - alpha1, 
		    std::max(-2.0*M_PI, parbound_.umax()-2.0*M_PI));
      u2 = std::min(parbound_.umax() + alpha2, 
		    std::min(2.0*M_PI, parbound_.umin()+2.0*M_PI));
      v1 = parbound_.vmin() - p1;
      v2 = parbound_.vmax() + p2;
      if (u2 - u1 > 2.0*M_PI)
	{
	  double udel = u2 - u1 - 2.0*M_PI;
	  u2 -= 0.5*udel;
	  u1 += 0.5*udel;
	}
      v1 = std::max(v1, -h);
    }
  setParameterBounds(u1, v1, u2, v2);
}

//===========================================================================
  void Cone::translate(const Point& vec)
//===========================================================================
{
  location_ += vec;
}

} // namespace Go
