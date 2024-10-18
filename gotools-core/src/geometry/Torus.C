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

#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sisl.h"
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
    if (fabs(from_upar) < ptol_ && fabs(to_upar - 2.0*M_PI) < ptol_) {
        from_upar = 0.0;
        to_upar = 2.0 * M_PI;
    }
    if (fabs(from_vpar) < ptol_ && fabs(to_vpar - 2.0*M_PI) < ptol_) {
        from_vpar = 0.0;
        to_vpar = 2.0 * M_PI;
    }
    if (is_degenerate_torus_) {
        if (select_outer_) {
            if (fabs(from_vpar + phi_) < ptol_)
                from_vpar = -phi_;
            if (fabs(to_vpar - phi_) < ptol_)
                to_vpar = phi_;
        }
        else {
            if (fabs(from_vpar - phi_) < ptol_)
                from_vpar = phi_;
            if (fabs(to_vpar - 2.0 * M_PI + phi_) < ptol_)
                to_vpar = 2.0 * M_PI - phi_;
        }
    }

    setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);

    // Swapped flag
    int isSwapped; // 0 or 1
    is >> isSwapped;
    bool has_param_int = (isSwapped >= 10);
    double start_u = from_upar, end_u = to_upar, start_v = from_vpar, end_v = to_vpar;
    if (has_param_int)
      {
	is >> start_u >> end_u >> start_v >> end_v;
      }
    setParameterDomain(start_u, end_u, start_v, end_v);
    isSwapped = isSwapped % 10;
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
    os << parbound_.umin() << " " << parbound_.umax() << endl
       << parbound_.vmin() << " " << parbound_.vmax() << endl;

    if (!isSwapped()) {
        os << "10" << endl;
    }
    else {
        os << "11" << endl;
    }
    os << domain_.umin() << " " << domain_.umax() << endl
       << domain_.vmin() << " " << domain_.vmax() << endl;

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
    torus->parbound_ = parbound_;
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
      double par = domain_.umin() -
	parbound_.umin()*(domain_.umax()-domain_.umin())/(parbound_.umax()-parbound_.umin());
	shared_ptr<Circle> circle = getMajorCircle(par);
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
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    vpar = parbound_.vmin() + 
      (vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
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
   double fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    upar = parbound_.umin() + fac1*(upar-domain_.umin());
    vpar = parbound_.vmin() + fac2*(vpar-domain_.vmin());

    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First derivatives
    double cosu = cos(upar);
    double sinu = sin(upar);
    double cosv = cos(vpar);
    double sinv = sin(vpar);
    pts[ind1] = fac1*(major_radius_ + minor_radius_ * cosv)
      * (-sinu * x_axis_ + cosu * y_axis_);
    pts[ind2] = fac2*minor_radius_ 
	* (-sinv * (cosu * x_axis_ + sinu * y_axis_)
	   + cosv * z_axis_);

    // Second order derivatives
    if (derivs > 1)
      {
	ind1 = 3;
	ind2 = 5;
	if (isSwapped())
	  swap(ind1, ind2);
	pts[ind1] = -fac1*fac1*(major_radius_ + 
			       minor_radius_ * cosv)*(cosu*x_axis_ +
						      sinu*y_axis_);
	pts[4] = -fac1*fac2*minor_radius_*sinv*(-sinu*x_axis_ + cosu*y_axis_);
	pts[ind2] = -fac2*fac2*minor_radius_*(cosv*(cosu*x_axis_+sinu*y_axis_) +
					      sinv*z_axis_);
      }

    if (derivs <= 2)
	return;

    // Second order and higher derivatives.
    MESSAGE("Third order or higher derivatives not yet implemented.");

}


//===========================================================================
void Torus::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    // This formula holds for both regular and degenerate tori.
    getOrientedParameters(upar, vpar);
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    vpar = parbound_.vmin() + 
      (vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
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
    bool torus_pardir_is_u = (isSwapped()) ? !pardir_is_u : pardir_is_u;

    // Major circle has v as constant parameter, i.e. the circle is parametrized in the u dir.
    shared_ptr<Circle> circle = (torus_pardir_is_u) ? getMajorCircle(parameter) : getMinorCircle(parameter);
    vector<shared_ptr<ParamCurve> > res;
    res.push_back(circle);

    return res;
}


//===========================================================================
Torus* Torus::subSurface(double from_upar, double from_vpar,
			 double to_upar, double to_vpar,
			 double fuzzy) const
//===========================================================================
{
    Torus* torus = clone();
    if (isSwapped())
      {
	double bound1 = parbound_.umin() + 
	  (from_vpar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
	double bound2 = parbound_.umin() + 
	  (to_vpar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
	double bound3 = parbound_.vmin() + 
	  (from_upar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	double bound4 = parbound_.vmin() + 
	  (to_upar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	torus->setParameterBounds(bound3, bound1, bound4, bound2);
      }
    else
      {
	double bound1 = parbound_.umin() + 
	  (from_upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
	double bound2 = parbound_.umin() + 
	  (to_upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
	double bound3 = parbound_.vmin() + 
	  (from_vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	double bound4 = parbound_.vmin() + 
	  (to_vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	torus->setParameterBounds(bound1, bound3, bound2, bound4);
      }
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
    // Algorithm:
    // 1) Find closest point on major circle (u-direction), including
    //    the bounds on u given by the domain of interest
    // 2) Find closest point on minor circle (v-direction), using the
    //    u value from 1) and including the bounds on v given by the
    //    domain of interest

    // Find relevant domain of interest
    RectDomain curr_domain_of_interest = parameterDomain();
    if (domain_of_interest != NULL) {
	curr_domain_of_interest.intersectWith(*domain_of_interest);
    }

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

    double vmin = parbound_.vmin();
    double vmax = parbound_.vmax();
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
#if 0
    closed_dir_u = (parbound_.umax() - parbound_.umin() == 2.0*M_PI);
    closed_dir_v = (parbound_.vmax() - parbound_.vmin() == 2.0*M_PI);
#else
    closed_dir_u = (fabs(parbound_.umax() - parbound_.umin() - 2.0*M_PI) < ptol_);
    closed_dir_v = (fabs(parbound_.vmax() - parbound_.vmin() - 2.0*M_PI) < ptol_);
#endif

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
    parbound_ = RectDomain(ll, ur);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void Torus::setParameterBounds(double from_upar, double from_vpar,
			       double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar)
	THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
	THROW("First v-parameter must be strictly less than second.");

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    if (fabs(from_upar) < ptol_)
      from_upar = 0.0;
    else if (fabs(2.0*M_PI-from_upar) < ptol_)
      from_upar = 2.0*M_PI;
    if (fabs(from_vpar) < ptol_)
      from_vpar = 0.0;
    else if (fabs(2.0*M_PI-from_vpar) < ptol_)
      from_vpar = 2.0*M_PI;
    
    if (fabs(to_upar) < ptol_)
      to_upar = 0.0;
    else if (fabs(2.0*M_PI-to_upar) < ptol_)
      to_upar = 2.0*M_PI;
    if (fabs(to_vpar) < ptol_)
      to_vpar = 0.0;
    else if (fabs(2.0*M_PI-to_vpar) < ptol_)
      to_vpar = 2.0*M_PI;
    
     // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar > -2.0 * M_PI - ptol_ && from_upar < -2.0 * M_PI)
      from_upar = -2.0 * M_PI;
    if (to_upar < 2.0 * M_PI + ptol_ && to_upar > 2.0 * M_PI)
      to_upar = 2.0 * M_PI;
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
	THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
	THROW("(to_upar - from_upar) must not exceed 2pi.");
    if (is_degenerate_torus_) {
	if (select_outer_) {
	  if (from_vpar > -phi_ - ptol_ && from_vpar < -phi_)
	    from_vpar = -phi_;
	  if (to_vpar < phi_ + ptol_ && to_vpar > phi_)
	    to_vpar = phi_;
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
      if (from_vpar > -2.0 * M_PI - ptol_ && from_vpar < -2.0 * M_PI)
	from_vpar = -2.0 * M_PI;
      if (to_vpar < 2.0 * M_PI + ptol_ && to_vpar > 2.0 * M_PI)
	to_vpar = 2.0 * M_PI;
	if (from_vpar < -2.0 * M_PI || to_vpar > 2.0 * M_PI)
	    THROW("v-parameters must be in [-2pi, 2pi].");
	if (to_vpar - from_vpar > 2.0 * M_PI)
	    THROW("(to_vpar - from_vpar) must not exceed 2pi.");
    }

    double start_u = parbound_.umin() + 
      (from_upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double end_u = parbound_.umin() + 
      (to_upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    double start_v = parbound_.vmin() + 
      (from_vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    double end_v = parbound_.vmin() + 
      (to_vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(start_u, start_v);
    Array<double, 2> ur2(end_u, end_v);
    domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Torus::setParameterDomain(double startpar_u, double endpar_u, 
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
    double angle = parbound_.umax() - parbound_.umin();

    SplineSurface* sstorus 
	= SweepSurfaceCreator::rotationalSweptSurface(*sccircle, angle,
						      location_, z_axis_);
    sstorus->basis_u().rescale(umin, umax);

    if (isSwapped())
        sstorus->swapParameterDirection();

    return sstorus;
}


//===========================================================================
SplineSurface* Torus::createNonRationalSpline(double eps) const
//===========================================================================
{
  // First fetch the first circular boundary curve in the minor
  // direction
  shared_ptr<Circle> circ = getMinorCircle(domain_.vmin());

  // Feth non-rational spline approximation
  shared_ptr<SplineCurve> crv(circ->createNonRationalSpline(0.5*eps));

  // Rotate this circle the valid angle around the main axis to 
  // create the non-rational spline surface
  // Note that the result will be rational if the tolerance is equal to zero
  int status;
  SISLCurve *qc = Curve2SISL(*crv);
  double *point = const_cast<double*>(location_.begin());
  double *axis =  const_cast<double*>(z_axis_.begin());
  SISLSurf *qs = NULL;
  s1302(qc, 0.5*eps, parbound_.umax()-parbound_.umin(), point, axis,
	&qs, &status);
  if (status < 0 || qs == NULL)
    return NULL;  // Approximation failed

  SplineSurface *surf = SISLSurf2Go(qs);
  surf->setParameterDomain(domain_.umin(), domain_.umax(),
			   domain_.vmin(), domain_.vmax());
  if (isSwapped())
    surf->swapParameterDirection();

  freeSurf(qs);
  return surf;
}

//===========================================================================
shared_ptr<Circle> Torus::getMajorCircle(double vpar) const
//===========================================================================
{
    vpar = parbound_.vmin() + 
      (vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    Point centre = location_ + minor_radius_ * sin(vpar) * z_axis_;
    double radius = major_radius_ + minor_radius_ * cos(vpar);
    shared_ptr<Circle> circle(new Circle(radius, centre, z_axis_, x_axis_));
    double umin = parbound_.umin();
    double umax = parbound_.umax();
    circle->setParamBounds(umin, umax);
    circle->setParameterInterval(domain_.umin(), domain_.umax());
    return circle;
}


//===========================================================================
shared_ptr<Circle> Torus::getMinorCircle(double upar) const
//===========================================================================
{
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    Point udir = cos(upar) * x_axis_ + sin(upar) * y_axis_;
    Point centre = location_ + major_radius_ * udir;
    Point newz = udir.cross(z_axis_);
    shared_ptr<Circle> circle(new Circle(minor_radius_, centre,
					 newz, udir));
    double vmin = parbound_.vmin();
    double vmax = parbound_.vmax();
    circle->setParamBounds(vmin, vmax);
   circle->setParameterInterval(domain_.vmin(), domain_.vmax());
    return circle;
}


//===========================================================================

//===========================================================================
  void Torus::enlarge(double len1, double len2, double len3, double len4)
//===========================================================================
{
  // Distances are given in geometry space, compute corresponding parameter
  // distances
  double alpha1 = len1/major_radius_;
  double alpha2 = len2/major_radius_;
  double alpha3 = len3/minor_radius_;
  double alpha4 = len4/minor_radius_;
  double lim1 = is_degenerate_torus_ ? (select_outer_ ? -phi_ : phi_) : -2.0*M_PI;
  double lim2 = is_degenerate_torus_ ? (select_outer_ ? phi_ : 2.0*M_PI-phi_) : 
    2.0*M_PI;
  double u1, u2, v1, v2;
  if (isSwapped())
    {
      u1 = std::max(parbound_.vmin() - alpha1, 
		    std::max(lim1, parbound_.vmax()-2.0*M_PI));
      u2 = std::min(parbound_.vmax() + alpha2, 
		    std::min(lim2, parbound_.vmin()+2.0*M_PI));
      v1 = std::max(parbound_.umin() - alpha3, 
		    std::max(lim1, parbound_.umax()-2.0*M_PI));
      v2 = std::min(parbound_.umax() + alpha4, 
		    std::min(lim2, parbound_.umin()+2.0*M_PI));
    }
  else
    {
      u1 = std::max(parbound_.umin() - alpha1, 
		    std::max(lim1, parbound_.umax()-2.0*M_PI));
      u2 = std::min(parbound_.umax() + alpha2, 
		    std::min(lim2, parbound_.umin()+2.0*M_PI));
      v1 = std::max(parbound_.vmin() - alpha3, 
		    std::max(lim1, parbound_.vmax()-2.0*M_PI));
      v2 = std::min(parbound_.vmax() + alpha4, 
		    std::min(lim2, parbound_.vmin()+2.0*M_PI));
    }
  if (u2 - u1 > 2.0*M_PI)
    {
      double udel = u2 - u1 - 2.0*M_PI;
      u2 -= 0.5*udel;
      u1 += 0.5*udel;
    }
  if (v2 - v1 > 2.0*M_PI)
    {
      double vdel = v2 - v1 - 2.0*M_PI;
      v2 -= 0.5*vdel;
      v1 += 0.5*vdel;
    }
  setParameterBounds(u1, v1, u2, v2);
}

//===========================================================================
  void Torus::translate(const Point& vec)
//===========================================================================
{
  location_ += vec;
}


} // namespace Go
