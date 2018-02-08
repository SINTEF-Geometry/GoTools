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

#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/SplineSurface.h"


using std::vector;
using std::max;
using std::min;
using std::endl;
using std::streamsize;
using std::swap;


namespace Go
{

//===========================================================================
Disc::Disc(Point centre, double radius, Point x_axis, Point normal,
    bool isSwapped) :
    centre_(centre), radius_(radius), x_axis_(x_axis), z_axis_(normal),
    centre_degen_(true)
//===========================================================================
{
    degen_angles_[0] = 0.0;
    degen_angles_[1] = 0.5 * M_PI;
    degen_angles_[2] = 1.0 * M_PI;
    degen_angles_[3] = 1.5 * M_PI;
    
    setCoordinateAxes();
    setDefaultDomain();

    if (isSwapped)
        swapParameterDirection();
}


//===========================================================================
void Disc::read (std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
        THROW("Invalid geometry file!");
    }

    int dim, centre_degen_int;
    is >> dim;
    centre_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> centre_
       >> radius_
       >> z_axis_
       >> x_axis_
       >> centre_degen_int;
    for (int i = 0; i < 4; ++i)
      is >> degen_angles_[i];

    if (centre_degen_int == 0)
      centre_degen_ = false;
    else if (centre_degen_int == 1)
      centre_degen_ = true;
    else
      THROW("Unknown input for centre_degen - must be 0 or 1");

    setCoordinateAxes();

    // "Reset" swapping
    isSwapped_ = false;

    // Parameter bounds. NOTE: Mind the order of the parameters!
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
  void Disc::write(std::ostream& os) const
  //===========================================================================
  {
    streamsize prev = os.precision(15);
    os << dimension() << endl
       << centre_ << endl
       << radius_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;

    if (centre_degen_)
      os << "1" << endl;
    else
      os << "0" << endl;
    for (int i = 0; i < 4; ++i)
      os << degen_angles_[i] << endl;

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
  int Disc::dimension() const
  //===========================================================================
  {
    return centre_.dimension();
  }

  //===========================================================================
  ClassType Disc::instanceType() const
  //===========================================================================
  {
    return classType();
  }

  //===========================================================================
  BoundingBox Disc::boundingBox() const
  //===========================================================================
  {
    return boundaryCircle().boundingBox();
  }

  //===========================================================================
  Disc* Disc::clone() const
  //===========================================================================
  {
    Disc* newDisc = new Disc(centre_, radius_, x_axis_, z_axis_, isSwapped_);
    newDisc->centre_degen_ = centre_degen_;
    for (int i = 0; i < 4; ++i)
      newDisc->degen_angles_[i] = degen_angles_[i];
    newDisc->domain_ = domain_;
    return newDisc;
  }

  //===========================================================================
  const RectDomain& Disc::parameterDomain() const
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
  vector<CurveLoop> Disc::allBoundaryLoops(double degenerate_epsilon) const
  //===========================================================================
  {
    MESSAGE("allBoundaryLoops() not implemented. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
  }

  //===========================================================================
  DirectionCone Disc::normalCone() const
  //===========================================================================
  {
    Point normal = z_axis_;
    if (isSwapped())
        normal *= -1.0;
    return DirectionCone(normal);
  }

  //===========================================================================
  DirectionCone Disc::tangentCone(bool pardir_is_u) const
  //===========================================================================
  {
    if (isSwapped())
        pardir_is_u = !pardir_is_u;

    double vmin = domain_.vmin();
    double vmax = domain_.vmax();

    vector<Point> pts(3);
    double u = 1.0;
    double v = 0.5*(vmin+vmax);
    if (isSwapped())
        swap(u, v);
    point(pts, u, v, 1);
    if (pardir_is_u)
      return DirectionCone(pts[1], 0.5*(vmax-vmin));
    else
      return DirectionCone(pts[2], 0.5*(vmax-vmin));
  }

  //===========================================================================
  void Disc::point(Point& pt, double upar, double vpar) const
  //===========================================================================
  {
    getOrientedParameters(upar, vpar); // In case of swapped
    pt = centre_
      + upar * (cos(vpar) * x_axis_
		+ sin(vpar) * y_axis_);
  }


  //===========================================================================
  void Disc::point(std::vector<Point>& pts, 
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
    for (int i = 0; i < totpts; ++i)
      {
	if (pts[i].dimension() != dim)
	  pts[i].resize(dim);
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
    double cosv = cos(vpar);
    double sinv = sin(vpar);

    pts[ind1] = cosv * x_axis_ + sinv * y_axis_;
    pts[ind2] = upar * (-sinv * x_axis_ + cosv * y_axis_);

    if (derivs == 1)
      return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");
  }

  //===========================================================================
  void Disc::normal(Point& n, double upar, double vpar) const
  //===========================================================================
  {
    n = z_axis_;
    if (isSwapped())
        n *= -1.0;
  }


  //===========================================================================
  vector<shared_ptr<ParamCurve> >
  Disc::constParamCurves(double parameter, bool pardir_is_u) const
  //===========================================================================
  {
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
  }


  //===========================================================================
  Disc* Disc::subSurface(double from_upar, double from_vpar,
			 double to_upar, double to_vpar,
			 double fuzzy) const
  //===========================================================================
  {
    Disc* newDisc = clone();
    newDisc->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return newDisc;
  }

  //===========================================================================
  vector<shared_ptr<ParamSurface> >
  Disc::subSurfaces(double from_upar, double from_vpar,
		    double to_upar, double to_vpar,
		    double fuzzy) const
  //===========================================================================
  {
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<Disc> newDisc(subSurface(from_upar, from_vpar,
					to_upar, to_vpar));
    res.push_back(newDisc);
    return res;
  }

  //===========================================================================
  double 
  Disc::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() doesn't make sense. Returning arbitrarily 0.0.");
    return 0.0;
  }


  //===========================================================================
  void Disc::closestPoint(const Point& pt,
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
    double umin = curr_domain_of_interest.umin();
    double umax = curr_domain_of_interest.umax();
    double vmin = curr_domain_of_interest.vmin();
    double vmax = curr_domain_of_interest.vmax();

    Point vec = pt - centre_;
    Point projected = pt - (vec*z_axis_)*z_axis_;
    double radius = centre_.dist(projected);
    double angle = x_axis_.angle(projected - centre_);
    double eps = 1.0e-4;
    clo_u = radius;
    clo_v = angle;
    getOrientedParameters(clo_u, clo_v);
    Point testpt;
    point(testpt, clo_u, clo_v);
    if (testpt.dist(projected) > eps) {
        if (!isSwapped()) {
            clo_v = 2.0*M_PI - angle;
        }
        else {
            clo_u = 2.0*M_PI - angle;
        }
    }
    point(testpt, clo_u, clo_v);
    if (testpt.dist(projected) > eps)
        THROW("This should never happen!");

    // Adjust for parameter domain
    if (clo_u < umin)
        clo_u = umin;
    if (clo_u > umax)
        clo_u = umax;
    if (clo_v < vmin)
        clo_v = vmin;
    if (clo_v > vmax)
        clo_v = vmax;

    point(clo_pt, clo_u, clo_v);
    clo_dist = clo_pt.dist(pt);

  }


  //===========================================================================
  void Disc::closestBoundaryPoint(const Point& pt,
				  double&        clo_u,
				  double&        clo_v, 
				  Point&       clo_pt,
				  double&        clo_dist,
				  double epsilon,
				  const RectDomain* rd,
				  double *seed) const
  //===========================================================================
  {
    MESSAGE("closestBoundaryPoint() not yet implemented");
  }


  //===========================================================================
  void Disc::getBoundaryInfo(Point& pt1, Point& pt2,
			     double epsilon, SplineCurve*& cv,
			     SplineCurve*& crosscv, double knot_tol) const
  //===========================================================================
  {
    MESSAGE("getBoundaryInfo() not yet implemented");
  }


  //===========================================================================
  bool Disc::isDegenerate(bool& b, bool& r,
			  bool& t, bool& l, double tolerance) const
  //===========================================================================
  {
    b = false;
    r = false;
    t = false;
    l = true;
    if (isSwapped()) {
        swap(b, l);
        swap(t, r);
    }
    return true;
  }


  //===========================================================================
  void Disc::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
  //===========================================================================
  {
      deg_corners.clear();
      if (domain_.umin() > 0.0)
          return;
      deg_corners.push_back(centre_);
  }


//===========================================================================
bool Disc::isBounded() const
//===========================================================================
{
  return true;
}


//===========================================================================
bool Disc::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    closed_dir_u = (domain_.umax() - domain_.umin() == 2.0*M_PI);
    closed_dir_v = false;
    if (isSwapped())
        swap(closed_dir_u, closed_dir_v);
    return (closed_dir_u || closed_dir_v);
}


//===========================================================================
  SplineSurface* Disc::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* Disc::createSplineSurface() const
//===========================================================================
{
    if (centre_degen_) {
        Circle boundary = boundaryCircle();
        SplineCurve* boundspl = boundary.createSplineCurve();
        int dim = dimension();
        bool rational = true;
        int numu = 2;
        int numv = boundspl->numCoefs();
        int ordu = 2;
        int ordv = boundspl->order();
        vector<double> knotsu(2, 0.0);
        knotsu.push_back(radius_);
        knotsu.push_back(radius_); // (0, 0, r, r)
        vector<double> knotsv(boundspl->knotsBegin(), boundspl->knotsEnd());
        vector<double> rcoefs(boundspl->rcoefs_begin(), boundspl->rcoefs_end());
        vector<double> degpt(centre_.begin(), centre_.begin()+dim);
        degpt.push_back(1.0); // w = 1
        for (int i = 0; i < numv; ++i) {
            rcoefs.insert(rcoefs.begin() + i*2*(dim+1), degpt.begin(), degpt.end());
        }
        SplineSurface* tmpsurf = new SplineSurface(numu, numv, ordu, ordv,
            knotsu.begin(), knotsv.begin(), rcoefs.begin(), dim, rational);
        
        if (isSwapped())
            tmpsurf->swapParameterDirection();

        RectDomain dom = parameterDomain();
        SplineSurface* surf = tmpsurf->subSurface(dom.umin(), dom.vmin(),
            dom.umax(), dom.vmax());
        return surf;
    }
    else {
        MESSAGE("createSplineSurface() not implemented for degenerate corners.");
        return NULL;
    }
    
}

//===========================================================================
void Disc::setParameterBounds(double from_upar, double from_vpar,
			      double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
        THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
        THROW("First v-parameter must be strictly less than second.");

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    const double num_tol = 1.0e-12;
    if ((from_upar < 0.0) && (from_upar > -num_tol))
        from_upar = 0.0;
    if ((from_vpar < -2.0 * M_PI) && (from_vpar > -2.0 * M_PI - num_tol))
        from_vpar = -2.0 * M_PI;
    if ((to_vpar > 2.0 * M_PI) && (to_vpar < 2.0 * M_PI + num_tol))
        to_vpar = 2.0 * M_PI;
    if ((to_vpar - from_vpar > 2.0 * M_PI) && (to_vpar - from_vpar < 2.0 * M_PI - num_tol))
        to_vpar = from_vpar + 2.0*M_PI;

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar < 0.0)
        THROW("from_upar must be >=  0.0");
    if (from_vpar < -2.0 * M_PI || to_vpar > 2.0 * M_PI)
        THROW("v-parameters must be in [-2pi, 2pi].");
    if (to_vpar - from_vpar > 2.0 * M_PI)
        THROW("(to_vpar - from_vpar) must not exceed 2pi.");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


  //===========================================================================
  void Disc::setCoordinateAxes()
  //===========================================================================
  {
    // The x- and y-axes defines a right-handed coordinate system for dimension 2.
    // The x-, y- and z-axes defines a right-handed coordinate system for dimension 3.

    if (dimension() == 2)
      {
	y_axis_ = Point(-x_axis_[1], x_axis_[0]);
	x_axis_.normalize();
	y_axis_.normalize();
      }
    else
      {
	z_axis_.normalize();
	Point tmp = x_axis_ - (x_axis_ * z_axis_) * z_axis_;
	if (tmp.length() == 0.0)
	  THROW("X-axis parallel to Z-axis.");

	x_axis_ = tmp;
	y_axis_ = z_axis_.cross(x_axis_);
	x_axis_.normalize();
	y_axis_.normalize();
      }
  }


  //===========================================================================
  void Disc::setDefaultDomain()
  //===========================================================================
  {
    setParameterBounds(0.0, 0.0, radius_, 2.0 * M_PI);
  }


  //===========================================================================
  Circle Disc::boundaryCircle() const
  //===========================================================================
  {
      // Circle of radius radius_. Parameter domain may have smaller umax.
      Circle c(radius_, centre_, z_axis_, x_axis_);
      c.setParamBounds(domain_.vmin(), domain_.vmax());
      return c;
  }


//===========================================================================
  void Disc::enlarge(double len1, double len2, double len3, double len4)
//===========================================================================
{
  // Distances are given in geometry space, compute corresponding parameter
  // distance in the rotational direction
  double u1, u2, v1, v2;
  if (isSwapped())
    {
      double alpha1 = len1/radius_;
      double alpha2 = len2/radius_;
      u1 = std::max(domain_.vmin() - alpha1, 
		    std::max(-2.0*M_PI, domain_.vmax()-2.0*M_PI));
      u2 = std::min(domain_.vmax() + alpha2, 
		    std::min(2.0*M_PI, domain_.vmin()+2.0*M_PI));
      v1 = domain_.umin() - len4;
      v2 = domain_.umax() + len3;
      if (u2 - u1 > 2.0*M_PI)
	{
	  double udel = u2 - u1 - 2.0*M_PI;
	  u2 -= 0.5*udel;
	  u1 += 0.5*udel;
	}
    }
  else
    {
      double alpha1 = len3/radius_;
      double alpha2 = len4/radius_;
      u1 = domain_.umin() - len1;
      u2 = domain_.umax() + len2;
      v1 = std::max(domain_.vmin() - alpha1, 
		    std::max(-2.0*M_PI, domain_.vmax()-2.0*M_PI));
      v2 = std::min(domain_.vmax() + alpha2, 
		    std::min(2.0*M_PI, domain_.vmin()+2.0*M_PI));
      if (v2 - v1 > 2.0*M_PI)
	{
	  double vdel = v2 - v1 - 2.0*M_PI;
	  v2 -= 0.5*vdel;
	  v1 += 0.5*vdel;
	}
     }
  setParameterBounds(u1, v1, u2, v2);
}


} // namespace Go
