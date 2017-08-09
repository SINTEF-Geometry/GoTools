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

#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Line.h"
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
Sphere::Sphere(double radius,
	       Point location, Point z_axis, Point x_axis,
               bool isSwapped)
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
    setParameterBounds(0.0, -0.5 * M_PI, 2.0 * M_PI, 0.5 * M_PI);

    if (isSwapped)
        swapParameterDirection();
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

    // "Reset" swapping
    isSwapped_ = false;

    // NB: Mind the sequence of parameters!
    double from_upar, from_vpar, to_upar, to_vpar;
    is >> from_upar >> to_upar
        >> from_vpar >> to_vpar;

    const double numtol = 1.0e-14; // To handle roundoff errors.
    if (fabs(from_upar) < numtol) {
        from_upar = 0.0;
    }
    if (fabs(to_upar - 2.0*M_PI) < numtol) {
        to_upar = 2.0 * M_PI;
    }
    if (fabs(from_vpar + 0.5*M_PI) < numtol) {
        from_vpar = -0.5 * M_PI;
    }
    if (fabs(to_vpar - 0.5*M_PI) < numtol) {
        to_vpar = 0.5 * M_PI;
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
void Sphere::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;

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
    SplineSurface* tmp = geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
Sphere* Sphere::clone() const
//===========================================================================
{
    Sphere* sph = new Sphere(radius_, location_, z_axis_, x_axis_, 
        isSwapped_);
    sph->domain_ = domain_;
    return sph;
}
    
//===========================================================================
const RectDomain& Sphere::parameterDomain() const
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
Sphere::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("allBoundaryLoops() not implemented. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone Sphere::normalCone() const
//===========================================================================
{
    // A rather unefficient hack...
    SplineSurface* tmp = geometrySurface();
    DirectionCone dc = tmp->normalCone();
    delete tmp;
    return dc;
}


//===========================================================================
DirectionCone Sphere::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    // A rather unefficient hack...
    SplineSurface* tmp = geometrySurface();
    DirectionCone dc = tmp->tangentCone(pardir_is_u);
    delete tmp;
    return dc;
}


//===========================================================================
void Sphere::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar); // In case of swapped
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

    // Swap parameters, if needed
    getOrientedParameters(upar, vpar);
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First derivatives
    pts[ind1] = radius_ * cos(vpar) * (-sin(upar) * x_axis_ + cos(upar) * y_axis_);
    pts[ind2] = radius_ * (-sin(vpar) * (cos(upar) * x_axis_ 
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
    getOrientedParameters(upar, vpar);
    n = cos(vpar) * (cos(upar) * x_axis_ + sin(upar) * y_axis_)
	+ sin(vpar) * z_axis_;;
    if (isSwapped())
        n *= -1.0;
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
    MESSAGE("nextSegmentVal() doesn't make sense. Returning arbitrarily 0.0.");
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
    RectDomain curr_domain_of_interest = parameterDomain();;
    if (domain_of_interest != NULL) {
        if (isSwapped_) {
            MESSAGE("Missing handling of swapped domain!");
        }
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

    getOrientedParameters(umin, vmin);
    getOrientedParameters(umax, vmax);

    // Find closest point on latitudinal circle (u-direction)
    double vmean = 0.5 * (vmin + vmax);
    shared_ptr<Circle> latitudinal_circle = getLatitudinalCircle(vmean);
    latitudinal_circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);

    // Find closest point on longitudinal circle (v-direction)
    shared_ptr<Circle> longitudinal_circle = getLongitudinalCircle(clo_u);
    longitudinal_circle->closestPoint(pt, vmin, vmax, clo_v, clo_pt, clo_dist);

    // We have what we need
    getOrientedParameters(clo_u, clo_v);
    return;
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
void Sphere::getBoundaryInfo(Point& pt1, Point& pt2,
			    double epsilon, SplineCurve*& cv,
			    SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
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
    if (domain_.vmin() == -0.5 * M_PI) {
	b = true;
	res = true;
        if (isSwapped())
            swap(b, l);
    }
    if (domain_.vmax() == 0.5 * M_PI) {
	t = true;
	res = true;
        if (isSwapped())
            swap(t, r);
    }
    return res;
}


//===========================================================================
shared_ptr<ElementaryCurve> 
Sphere::getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
				const Point* start_par_pt, const Point* end_par_pt) const 
//===========================================================================
{
  double angtol = 0.01;

  // Default is not simple elementary parameter curve exists
  shared_ptr<ElementaryCurve> dummy;

  // We search for a linear parameter curve. This can only occur if the
  // space curve is a circle with its centre at the sphere centre or the
  // centre at the z-axis of the sphere with where the normal of the plane
  // in which the circle lives is parallel with the z-axis
  Circle *circle = NULL;
  if (space_crv->instanceType() == Class_Circle)
    {
      circle = dynamic_cast<Circle*>(space_crv);
    }

  if (!circle)
    return dummy;

  // Check if the parameter curve corresponding to the space curve is a line
  Point centre = circle->getCentre();
  Point normal = circle->getNormal();
  if (normal.dimension() != 3)
    return dummy;  // No plane specified

  double dist = location_.dist(centre);
  double ang = z_axis_.angle(normal);
  ang = std::min(ang, fabs(M_PI-ang));
  Point pt = location_ + ((centre - location_)*z_axis_)*z_axis_;
  double dist2 = centre.dist(pt);

  // Bookkeeping related to swapped parameters
  int ind1 = 0;
  int ind2 = 1;
  if (isSwapped())
      swap(ind1, ind2);

  if (dist < tol || (dist2 < tol && ang < angtol))
    {
      // Constant parameter curve
      bool closed = circle->isClosed();
      double t1 = circle->startparam();
      double t2 = (closed) ? 0.5*(t1 + circle->endparam()) : circle->endparam();
      int idx;

      if (dist < tol && ang > angtol)
	idx = ind1; // 0
      else
	idx = ind2; // 1

      // Project endpoints (or startpoint and midpoint) of the circle onto
      // the sphere
      double parval1[2], parval2[2];
      double d1, d2;
      Point close1, close2;
      Point pos1 = space_crv->ParamCurve::point(t1);
      Point pos2 = space_crv->ParamCurve::point(t2);
      closestPoint(pos1, parval1[0], parval1[1], close1, d1, tol);
      closestPoint(pos2, parval2[0], parval2[1], close2, d2, tol);
      if (d1 > tol || d2 > tol)
	return dummy;

      if (dist < tol && ang > angtol)
	{
	  // It might be necessary to adjust the parameter value at a
	  // degenerate point to get the correct circle
	  double degdist1 = std::min(fabs(parval1[ind2]-0.5*M_PI), 
				     fabs(parval1[ind2]+0.5*M_PI));
	  double degdist2 = std::min(fabs(parval2[ind2]-0.5*M_PI), 
				     fabs(parval2[ind2]+0.5*M_PI));
	  if (degdist1 < degdist2)
	    parval1[ind1] = parval2[ind1];
	  else
	    parval2[ind1] = parval1[ind1];

          // If both values are at a degenerate point, we need a "tiebreaker"
          if (degdist1 < tol && degdist2 < tol) {
              Point midpt = space_crv->ParamCurve::point(0.5*(t1+t2));
              double midpar[2];
              Point midclo;
              double middist;
              closestPoint(midpt, midpar[0], midpar[1], midclo, middist, tol);
              parval1[ind1] = midpar[ind1];
              parval2[ind1] = midpar[ind1];
          }
	}
	  
      // Check mid point
      Point par1(2), par2(2);
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
//	  MESSAGE("Avoid computing par1.");
	  par1 = *start_par_pt;
      }
      if (end_par_pt != NULL)
      {
//	  MESSAGE("Avoid computing par2.");
	  par2 = *end_par_pt;
      }
      shared_ptr<Line> param_cv(new Line(par1, par2, 
					 space_crv->startparam(), space_crv->endparam()));
      
      // TEST
      Point p1 = param_cv->ParamCurve::point(param_cv->startparam());
      Point p2 = param_cv->ParamCurve::point(param_cv->endparam());
      
      return param_cv;
    }
  else
    return dummy;
}

//===========================================================================
bool Sphere::isBounded() const
//===========================================================================
{
  return true;
}


//===========================================================================
bool Sphere::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    closed_dir_u = (domain_.umax() - domain_.umin() == 2.0*M_PI);
    closed_dir_v = false;
    if (isSwapped())
        swap(closed_dir_u, closed_dir_v);
    return (closed_dir_u || closed_dir_v);
}


///===========================================================================
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

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar < 0.0 || to_upar > 2.0 * M_PI)
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
    double tmpu = 0.0;
    double tmpv = vmin;
    getOrientedParameters(tmpu, tmpv);
    point(llpt, tmpu, tmpv);
    tmpu = umax - umin;
    tmpv = vmax;
    getOrientedParameters(tmpu, tmpv);
    point(urpt, tmpu, tmpv);
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

    if (isSwapped())
        subpatch->swapParameterDirection();

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
bool Sphere::isAxisRotational(Point& centre, Point& axis, Point& vec,
				double& angle)
//===========================================================================
{
  // @@@ VSK. This test is not general enough for a sphere. There are more 
  // freedom in the axis direction and vector direction
  centre = location_;
  axis = z_axis_;
  if (domain_.umin() == 0.0)
    vec = x_axis_;
  else
    {
      Point pt;
      point(pt, domain_.umin(), domain_.vmin());
      vec = pt - location_;
      vec.normalize();
    }
  angle = domain_.umax() - domain_.umin();

  return true;
}

//===========================================================================

} // namespace Go
