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

#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Circle.h"
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
Cylinder::Cylinder(double radius,
                   Point location, Point z_axis, Point x_axis,
                   bool isSwapped)
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
    double inf = numeric_limits<double>::infinity();
    setParameterBounds(0.0, -inf, 2.0 * M_PI, inf);
    setParameterDomain(0.0, 2.0*M_PI, -inf, inf);

    if (isSwapped)
        swapParameterDirection();
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
	try
	{
	    setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
	    setParameterDomain(start_u, end_u, start_v, end_v);
	}
	catch (...)
	{ // We want the read routine to continue reading data, not a good strategy to throw before all object data is parsed.
	    MESSAGE("Failed setting parameter bounds.");
	}
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
void Cylinder::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;

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
        // Create a SplineSurface with a large but finite BoundingBox
        SplineSurface* tmp = geometrySurface();
        BoundingBox box = tmp->boundingBox();
        delete tmp;
        return box;
    }

    // Note: We are using domain_ on purpose, because domain_'s
    // v-direction is always the "linear" direction, no matter what
    // isSwapped_ is.
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
Cylinder* Cylinder::clone() const
//===========================================================================
{
    Cylinder* cyl = new Cylinder(radius_, location_, z_axis_, x_axis_, 
        isSwapped_);
    cyl->parbound_ = parbound_;
    cyl->domain_ = domain_;
    return cyl;
}
    
    
//===========================================================================
const RectDomain& Cylinder::parameterDomain() const
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
DirectionCone Cylinder::normalCone() const
//===========================================================================
{
    // Using domain_ on purpose, disregarding isSwapped_.
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
DirectionCone Cylinder::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    // This function needs testing...

    if (isSwapped())
        pardir_is_u = !pardir_is_u;

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
    getOrientedParameters(upar, vpar); // In case of swapped
    upar = parbound_.umin() + 
      (upar-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    if (isBounded())
      {
	vpar = parbound_.vmin() + 
	  (vpar-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
      }
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
    point(pts[0], upar, vpar); // The method will swap if needed.
    if (derivs == 0)
        return;

    // Swap parameters, if needed
    getOrientedParameters(upar, vpar);
    double fac1 = 1.0, fac2 = 1.0;
    fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    upar = parbound_.umin() + fac1*(upar-domain_.umin());
    if (isBounded())
      {
	fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	vpar = parbound_.vmin() + fac2*(vpar-domain_.vmin());
      }

    // First derivatives. TESTME
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);
    pts[ind1] = fac1*radius_ * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[ind2] = fac2*z_axis_;
    if (derivs == 1)
        return;

    // Second order and higher derivatives. TESTME!
    for (int i = 2; i <= derivs; ++i) {
      fac1 *= (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
      int index = i*(i+1)/2;
      if (isSwapped())
	index = (i+1)*(i+2)/2 - 1;
      pts[index] = fac1*radius_ * (cos(upar + i*0.5*M_PI) * x_axis_
				   + sin(upar + i*0.5*M_PI) * y_axis_);
    }

}


//===========================================================================
void Cylinder::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar);
    double fac1 = (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    upar = parbound_.umin() + fac1*(upar-domain_.umin());
    if (isBounded())
      {
	double fac2 = (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
	vpar = parbound_.vmin() + fac2*(vpar-domain_.vmin());
	n = fac1*fac2*cos(upar) * x_axis_ + sin(upar) * y_axis_;
      }
    else
      n = fac1*cos(upar) * x_axis_ + sin(upar) * y_axis_;
    if (isSwapped())
        n *= -1.0;
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
Cylinder::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    bool real_pardir_is_u = (isSwapped()) ? !pardir_is_u : pardir_is_u;
    vector<shared_ptr<ParamCurve> > res;
    if (real_pardir_is_u)
    {
        shared_ptr<ParamCurve> circle = getCircle(parameter);
        res.push_back(circle);
    }
    else
    {
        if (!isBounded())
        {
            // Since the point() function will swap the input if parameter directions is swapped, we must do the same.
            Point par_zero(parameter, 0.0);
            getOrientedParameters(par_zero[0], par_zero[1]);
            Point pt_zero = ParamSurface::point(par_zero[0], par_zero[1]);
            // The z_axis_ is the direction of the line.
            shared_ptr<Line> line(new Line(pt_zero, z_axis_));
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
Cylinder::constParamCurve(double iso_par, bool pardir_is_u,
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
Cylinder* Cylinder::subSurface(double from_upar, double from_vpar,
                               double to_upar, double to_vpar,
                               double fuzzy) const
//===========================================================================
{
    Cylinder* cylinder = clone();
    bool bounded = isBounded();
    double fac1 = 1.0, fac2 = 1.0;
    fac1 = 
      (parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());
    if (bounded)
      {
	fac2 =
	  (parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
      }
    if (isSwapped())
      {
	double bound1 = from_vpar, bound2 = to_vpar;
	double bound3 = from_upar, bound4 = to_upar;
	bound1 = parbound_.umin() + fac1*(from_vpar-domain_.umin());
	bound2 = parbound_.umin() + fac1*(to_vpar-domain_.umin());
	if (bounded)
	  {
	    bound3 = parbound_.vmin() + fac2*(from_upar-domain_.vmin());
	    bound4 = parbound_.vmin() + fac2*(to_upar-domain_.vmin());
	  }
    	cylinder->setParameterBounds(bound3, bound1, bound4, bound2);
      }
    else
      {
	double bound1 = from_upar, bound2 = to_upar;
	double bound3 = from_vpar, bound4 = to_vpar;
	bound1 = parbound_.umin() + fac1*(from_upar-domain_.umin());
	bound2 = parbound_.umin() + fac1*(to_upar-domain_.umin());
	if (bounded)
	  {
	    bound3 = parbound_.vmin() + fac2*(from_vpar-domain_.vmin());
	    bound4 = parbound_.vmin() + fac2*(to_vpar-domain_.vmin());
	  }
	cylinder->setParameterBounds(bound1, bound3, bound2, bound4);
      }
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
    axis.setParamBounds(parbound_.vmin(), parbound_.vmax());
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    axis.setParameterInterval(vmin, vmax);
    // No point in using seed for closest point on a line.
    axis.closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);

    // Find closest point on the circle at the found v parameter
    shared_ptr<Circle> circle = getCircle(clo_v);
    double umin = domain_.umin();
    double umax = domain_.umax();
    const int circle_ind = isSwapped() ? 1 : 0;
    const double* circle_seed = (seed) ? &seed[circle_ind] : NULL;
#if 0 // Seems best to let the circle handle this case by actually using the seed.
    // If seed is at the seam we make sure we de not flip over the seam.
    if (circle_seed)
    {
      double seed =  parbound_.umin() + 
	((*circle_seed)-domain_.umin())*(parbound_.umax()-parbound_.umin())/
	(domain_.umax()-domain_.umin());
      double ustart = parbound_.umin() + 
	(umin-domain_.umin())*(parbound_.umax()-parbound_.umin())/
	(domain_.umax()-domain_.umin());
      double uend = parbound_.umin() + 
	(umax-domain_.umin())*(parbound_.umax()-parbound_.umin())/
	(domain_.umax()-domain_.umin());

      const double domain_fraction = 0.1;
      if ((fabs(seed) < epsilon) && fabs(2*M_PI - uend) < epsilon)
	{
	  uend = domain_fraction*2.0*M_PI;
	}
      if ((fabs(seed) < epsilon) && fabs(ustart) < epsilon)
	{
	  uend = domain_fraction*2.0*M_PI;
	}
	umax = domain_.umin() + 
	  (uend-parbound_.umin())*(domain_.umax()-domain_.umin())/
	  (parbound_.umax()-parbound_.umin());
    }
#endif
    circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist, circle_seed);

    // Take care of swapped parameters
    getOrientedParameters(clo_u, clo_v);
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
  if (!isBounded())
    {
      MESSAGE("closestBoundaryPoints does not make sense for unbounded surfaces");
      clo_dist = -1.0;
      return;
    }
  
    RectDomain domain = containingDomain();
    if (!rd)
	rd = &domain;
    shared_ptr<ParamCurve> bdcrv;
    clo_dist = 1.0e10;  // Initialize with a large number
    Point cpt;
    double cdist, cpar;

    // Checking closest point on the bottom boundary
    bdcrv = constParamCurve(rd->vmin(), true, rd->umin(), rd->umax());
    bdcrv->closestPoint(pt, rd->umin(), rd->umax(),
			clo_u, clo_pt, clo_dist, seed);
    clo_v = rd->vmin();

    // Checking the right boundary
    bdcrv = constParamCurve(rd->umax(), false, rd->vmin(), rd->vmax());
    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
			(seed == 0) ? seed : seed+1);
    if (cdist < clo_dist)
      {
	clo_pt = cpt;
	clo_u = rd->umax();
	clo_v = cpar;
	clo_dist = cdist;
      }

   // Checking the upper boundary
    bdcrv = constParamCurve(rd->vmax(), true, rd->umin(), rd->umax());
    bdcrv->closestPoint(pt, rd->umin(), rd->umax(), cpar, cpt, cdist,
			seed);
    if (cdist < clo_dist)
      {
	clo_pt = cpt;
	clo_u = cpar;
	clo_v = rd->vmax();
	clo_dist = cdist;
      }
    
   // Checking the left boundary
    bdcrv = constParamCurve(rd->umin(), false, rd->vmin(), rd->vmax());
    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
			(seed == 0) ? seed : seed+1);
    if (cdist < clo_dist)
      {
	clo_pt = cpt;
	clo_u = rd->umin();
	clo_v = cpar;
	clo_dist = cdist;
      }
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
shared_ptr<ElementaryCurve> 
Cylinder::getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
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
      idx = ind1; // 0
    }
  else if (space_crv->instanceType() == Class_Circle)
    {
      t1 = space_crv->startparam();
      closed = ((Circle*)(space_crv))->isClosed();
      t2 = (closed) ? 0.5*(t1 + space_crv->endparam()) :
	space_crv->endparam();
      idx = ind2; // 1
    }
  else
    return dummy;
      
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

   Point par1(2), par2(2);
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
	  double ptol = 1.0e-4;
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
	  if (mid.dist(cv_mid) > tol) {
              // We try one more thing: Shifting the angular parameters -2pi
              double tmppar1 = u1 - 2.0 * M_PI;
              double tmppar2 = u2 - 2.0 * M_PI;
              if (fabs(u2 - tmppar1) <= 2.0 * M_PI) {
                  u1 = tmppar1;
		  parval1[ind1] = domain_.umin()+(u1 - parbound_.umin())/fac;
              }
              else if (fabs(tmppar2 - u1) <= 2.0 * M_PI) {
                  u2 = tmppar2;
		  parval2[ind1] =domain_.umin()+(u2 - parbound_.umin())/fac;
              }
              else {
                  return dummy;
              }
              par1[idx] = par2[idx] = 0.5*(parval1[idx] + parval2[idx]);
              par1[1-idx] = parval1[1-idx];
              par2[1-idx] = parval2[1-idx];
              mid = this->ParamSurface::point(0.5*(par1[0]+par2[0]), 0.5*(par1[1]+par2[1]));
              if (mid.dist(cv_mid) > tol) {
                  return dummy;
              }
          }
	}
      else
	return dummy;  // Linear parameter curve not close enough
    }

  bool pt1_in_dom = domain_.isInDomain(Vector2D(par1[ind1], par1[ind2]), tol);
  bool pt2_in_dom = domain_.isInDomain(Vector2D(par2[ind1], par2[ind2]), tol);
  if (!(pt1_in_dom && pt2_in_dom))
  {
      MESSAGE("End pt(s) not in domain! Suspecting that the seam must be moved.");
  }

  // bool pt1_at_seam = domain_.isInDomain(Vector2D(par1[ind1], par1[ind2]), tol);
  // bool pt2_at_seam = domain_.isInDomain(Vector2D(par2[ind1], par2[ind2]), tol);

  u3 = parbound_.umin() + fac*(par1[ind1] - domain_.umin());
  u4 = parbound_.umin() + fac*(par2[ind1] - domain_.umin());
  if (closed)
  {
    double sign = (u3 > M_PI) ? -1.0 : 1.0;
    u4 = u3 + sign*2.0*M_PI;
    par2[ind1] = domain_.umin()+(u4 - parbound_.umin())/fac;
  }
  bool pt1_at_seam = std::min(fabs(u3 - parbound_.umin()), fabs(parbound_.umax()) - u3) < tol;
  bool pt2_at_seam = std::min(fabs(u4 - parbound_.umin()), fabs(parbound_.umax()) - u4) < tol;
  if (start_par_pt != NULL)
  {
      //MESSAGE("Avoid computing par1.");
      par1 = *start_par_pt;
      if (pt1_at_seam && pt2_at_seam && (end_par_pt == NULL))
      {
	  // @@sbr201506 We make sure that a seam point does not flip as that is most likely the correct approach.
	  // This should be handled more robustly by first projecting curves that do not follow the seam.
	  // There should also be a marching approach to ensure our choice is correct.
	  par2[ind1] = par1[ind1];
      }
  }
  if (end_par_pt != NULL)
  {
      //MESSAGE("Avoid computing par2.");
      par2 = *end_par_pt;
      if (pt1_at_seam && pt2_at_seam && (start_par_pt == NULL))
      {
	  par1[ind1] = par2[ind1];
      }
  }
  if (space_crv->isReversed())
    std::swap(par1, par2);
  shared_ptr<Line> param_cv(new Line(par1, par2, 
				     space_crv->startparam(), 
				     space_crv->endparam(),
				     space_crv->isReversed()));

#ifndef NDEBUG
  {
      // TEST
      Point p1 = param_cv->ParamCurve::point(param_cv->startparam());
      Point p2 = param_cv->ParamCurve::point(param_cv->endparam());
      int stop_break = 1;
  }
#endif

  return param_cv;
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
    if (to_upar - from_upar > 2.0 * M_PI+ptol_)
        THROW("(to_upar - from_upar) must not exceed 2pi.");
    if (to_upar - from_upar > 2.0 * M_PI)
      to_upar = 2.0 * M_PI - from_upar;

    double fac1 = 
      (domain_.umax()-domain_.umin())/(parbound_.umax() - parbound_.umin());
    double fac2 = 1.0;
    if (bounded)
      fac2 = 
	(domain_.vmax()-domain_.vmin())/(parbound_.vmax() - parbound_.vmin());
    double start_u = domain_.umin() + fac1*(from_upar-parbound_.umin());
    double end_u = domain_.umin() + fac1*(to_upar-parbound_.umin());
    double start_v = (bounded) ?
      domain_.vmin() + fac2*(from_vpar-parbound_.vmin()) : from_vpar;
    double end_v = (bounded) ?
      domain_.vmin() + fac2*(to_vpar-parbound_.vmin()) : to_vpar;

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(start_u, start_v);
    Array<double, 2> ur2(end_u, end_v);
    domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Cylinder::setParamBoundsU(double from_upar, double to_upar)
//===========================================================================
{
  RectDomain tmp_domain = parbound_;
    double from_vpar = tmp_domain.vmin();
    double to_vpar = tmp_domain.vmax();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    if (from_upar > -2.0 * M_PI - ptol_ && from_upar < -2.0 * M_PI)
      from_upar = -2.0 * M_PI;
    if (to_upar < 2.0 * M_PI + ptol_ && to_upar > 2.0 * M_PI)
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
    double end_u = domain_.umin() + fac1*(to_upar-parbound_.umin());

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(start_u, domain_.vmin());
    Array<double, 2> ur2(end_u, domain_.vmax());
    domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Cylinder::setParamBoundsV(double from_vpar, double to_vpar)
//===========================================================================
{
  RectDomain tmp_domain = parbound_;
    double from_upar = tmp_domain.umin();
    double to_upar = tmp_domain.umax();
    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    double fac2 = 1.0;
    double start_v = from_vpar;
    double end_v = to_vpar;
    if (isBounded())
      {
	fac2 = 
	  (domain_.vmax()-domain_.vmin())/(parbound_.vmax() - parbound_.vmin());
	start_v = domain_.vmin() + fac2*(from_vpar-parbound_.vmin());
	end_v = domain_.vmin() + fac2*(to_vpar-parbound_.vmin());
      }

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    parbound_ = RectDomain(ll, ur);

    Array<double, 2> ll2(domain_.umin(), start_v);
    Array<double, 2> ur2(domain_.umax(), end_v);
    domain_ = RectDomain(ll2, ur2);
}


//===========================================================================
void Cylinder::setParameterDomain(double startpar_u, double endpar_u, 
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
bool Cylinder::isBounded() const
//===========================================================================
{
    // It is enough to check the v-direction, since the u-direction is
    // always bounded.

    return parbound_.vmin() > -numeric_limits<double>::infinity() &&
        parbound_.vmax() < numeric_limits<double>::infinity();
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
    double vmin = parbound_.vmin();
    double vmax = parbound_.vmax();
    if (!isBounded()) {
        double max = 1.0e8; // "Large" number...
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
    double umin = parbound_.umin();
    double umax = parbound_.umax();
    double tmpu = umax - umin;
    Array<double, 2> ll(umin, vmin);
    Array<double, 2> ur(umax, vmax);
    RectDomain rd(ll, ur);
    if (umin < 0.0)
      {
	// Make sure to get a large enough resulting surface
	rd = surface.containingDomain();
      }
    if (tmpu < 2.0 * M_PI) {
        double vmid = 0.5 * (vmin + vmax);
        Point pt, tmppt;
        double tmpv = vmid;
	pt = location_ 
	  + radius_ * (cos(tmpu) * x_axis_ + sin(tmpu) * y_axis_)
	  + tmpv * z_axis_;
        double tmpdist;
        double epsilon = 1.0e-10;
	double seed[2];
	seed[0] = tmpu;
	seed[1] = tmpv;
        surface.closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon, 
			     NULL, seed);
        if (tmpu < epsilon) {
            tmpu = 2.0 * M_PI;
        }
    }
    SplineSurface* subpatch = surface.subSurface(0.0, vmin, tmpu, vmax);
    subpatch->basis_u().rescale(domain_.umin(), domain_.umax());
    if (isBounded())
      subpatch->basis_v().rescale(domain_.vmin(), domain_.vmax());
    GeometryTools::translateSplineSurf(-location_, *subpatch);
    GeometryTools::rotateSplineSurf(z_axis_, umin, *subpatch);
    GeometryTools::translateSplineSurf(location_, *subpatch);

    if (isSwapped())
        subpatch->swapParameterDirection();

    return subpatch;
}


//===========================================================================
SplineSurface* Cylinder::createNonRationalSpline(double eps) const
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
shared_ptr<Circle> Cylinder::getCircle(double par) const
//===========================================================================
{
  if (isBounded())
    par = parbound_.vmin() + 
      (par-domain_.vmin())*(parbound_.vmax()-parbound_.vmin())/(domain_.vmax()-domain_.vmin());
    Point centre = location_ + par * z_axis_;
    shared_ptr<Circle> circle(new Circle(radius_, centre, z_axis_, x_axis_));
    // Note: We are using domain_ on purpose, because parbound_'s
    // u-direction is always the angular direction, no matter what
    // isSwapped_ is.
    double umin = parbound_.umin();
    double umax = parbound_.umax();
    circle->setParamBounds(umin, umax);
    circle->setParameterInterval(domain_.umin(), domain_.umax());
    return circle;
}


//===========================================================================
bool Cylinder::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
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
bool Cylinder::isAxisRotational(Point& centre, Point& axis, Point& vec,
				double& angle)
//===========================================================================
{
  centre = location_;
  axis = z_axis_;
  if (parbound_.umin() == 0.0)
    vec = x_axis_;
  else
    {
      Point pt;
      point(pt, parbound_.umin(), parbound_.vmin());
      vec = pt - location_;
      vec.normalize();
    }
  angle = parbound_.umax() - parbound_.umin();

  return true;
}

//===========================================================================
bool Cylinder::isLinear(Point& dir1, Point& dir2, double tol)
//===========================================================================
{
  dir1 = z_axis_;
  dir2.resize(0);
  return true;
}


//===========================================================================
void Cylinder::rotate(double rot_ang_rad)
//===========================================================================
{
    GeometryTools::rotatePoint(z_axis_, rot_ang_rad, x_axis_);
    GeometryTools::rotatePoint(z_axis_, rot_ang_rad, y_axis_);
}

//===========================================================================
  void Cylinder::enlarge(double len1, double len2, double len3, double len4)
//===========================================================================
{
  // Distances are given in geometry space, compute corresponding parameter
  // distance in the rotational direction
  double u1, u2, v1, v2;
  if (isSwapped())
    {
      double alpha1 = len3/radius_;
      double alpha2 = len4/radius_;
      u1 = parbound_.vmin() - len1;
      u2 = parbound_.vmax() + len2;
      v1 = std::max(parbound_.umin() - alpha1, 
		    std::max(-2.0*M_PI, parbound_.umax()-2.0*M_PI));
      v2 = std::min(parbound_.umax() + alpha2, 
		    std::min(2.0*M_PI, parbound_.umin()+2.0*M_PI));
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
      u1 = std::max(parbound_.umin() - alpha1, 
		    std::max(-2.0*M_PI, parbound_.umax()-2.0*M_PI));
      u2 = std::min(parbound_.umax() + alpha2, 
		    std::min(2.0*M_PI, parbound_.umin()+2.0*M_PI));
      v1 = parbound_.vmin() - len3;
      v2 = parbound_.vmax() + len4;
      if (u2 - u1 > 2.0*M_PI)
	{
	  double a_tol = 1.0e-10;
	  double udel = u2 - u1 - 2.0*M_PI + a_tol;
	  u2 -= 0.5*udel;
	  u1 += 0.5*udel;
	}
     }
  setParameterBounds(u1, v1, u2, v2);
}

//===========================================================================
  void Cylinder::translate(const Point& vec)
//===========================================================================
{
  location_ += vec;
}

//===========================================================================
  bool Cylinder::atSeam(int dir, double parval) const
//===========================================================================
{
  int dir2 = dir - 1;
  if (isSwapped())
    dir2 = 1 - dir2;
  if (dir2 != 0)
    return false;  // Not the rotational direction

  double upar = parbound_.umin() + 
      (parval-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());

  if (fabs(parval) < ptol_)
    return true;
  if (fabs(2.0*M_PI - parval) < ptol_)
    return true;

  return false;
}

//===========================================================================
  bool Cylinder::fullPeriod(int dir, double parval1, double parval2) const
//===========================================================================
{
  int dir2 = dir - 1;
  if (isSwapped())
    dir2 = 1 - dir2;
  if (dir2 != 0)
    return false;  // Not the rotational direction

  double upar1 = parbound_.umin() + 
      (parval1-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());

  double upar2 = parbound_.umin() + 
      (parval2-domain_.umin())*(parbound_.umax()-parbound_.umin())/(domain_.umax()-domain_.umin());

  if (fabs(upar2 - upar1) > 2*M_PI-ptol_ && 
      fabs(upar2 - upar1) < 2*M_PI+ptol_)
    return true;

  return false;
}

} // namespace Go
