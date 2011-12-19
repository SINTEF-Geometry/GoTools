//===========================================================================
//
// File : Disc.C
//
// Created: Thu Nov  5 12:27:16 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/SplineSurface.h"


using std::vector;
using std::max;
using std::min;
using std::endl;


namespace Go
{

//===========================================================================
Disc::Disc(Point centre, double radius, Point x_axis, Point normal) :
    centre_(centre), radius_(radius), x_axis_(x_axis), z_axis_(normal),
    centre_degen_(false)
//===========================================================================
{
    degen_angles_[0] = 0.0;
    degen_angles_[1] = 0.5 * M_PI;
    degen_angles_[2] = 1.0 * M_PI;
    degen_angles_[3] = 1.5 * M_PI;
    
    setCoordinateAxes();
    setDefaultDomain();
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
    setDefaultDomain();
  }


  //===========================================================================
  void Disc::write(std::ostream& os) const
  //===========================================================================
  {
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
    Disc* newDisc = new Disc(centre_, radius_, x_axis_, z_axis_);
    newDisc->centre_degen_ = centre_degen_;
    for (int i = 0; i < 4; ++i)
      newDisc->degen_angles_[i] = degen_angles_[i];
    return newDisc;
  }

  //===========================================================================
  const RectDomain& Disc::parameterDomain() const
  //===========================================================================
  {
    return domain_;
  }

  //===========================================================================
  CurveLoop Disc::outerBoundaryLoop(double degenerate_epsilon) const
  //===========================================================================
  {
    MESSAGE("Not implememnted. Returns an empty loop.");
    CurveLoop loop;
    return loop;
  }


  //===========================================================================
  vector<CurveLoop> Disc::allBoundaryLoops(double degenerate_epsilon) const
  //===========================================================================
  {
    MESSAGE("Not implememnted. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
  }

  //===========================================================================
  DirectionCone Disc::normalCone() const
  //===========================================================================
  {
    return DirectionCone(z_axis_);
  }

  //===========================================================================
  DirectionCone Disc::tangentCone(bool pardir_is_u) const
  //===========================================================================
  {
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();

    vector<Point> pts;
    point(pts, 1.0, 0.5*(vmin+vmax), 1);
    if (pardir_is_u)
      return DirectionCone(pts[1], 0.5*(vmax-vmin));
    else
      return DirectionCone(pts[2], 0.5*(vmax-vmin));
  }

  //===========================================================================
  void Disc::point(Point& pt, double upar, double vpar) const
  //===========================================================================
  {
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

    // First derivatives
    double cosv = cos(vpar);
    double sinv = sin(vpar);

    pts[1] = cosv * x_axis_ + sinv * y_axis_;
    pts[2] = upar * (-sinv * x_axis_ + cosv * y_axis_);

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
    newDisc->setParameterDomain(from_upar, from_vpar, to_upar, to_vpar);
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
    MESSAGE("Does not make sense. Return arbitrarily zero.");
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
    MESSAGE("closesPoint() not yet implemented");
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
    MESSAGE("closesBoundaryPoint() not yet implemented");
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
  void Disc::turnOrientation()
  //===========================================================================
  {
    swapParameterDirection();
  }



  //===========================================================================
  void Disc::swapParameterDirection()
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not implemented.");
  }


  //===========================================================================
  void Disc::reverseParameterDirection(bool direction_is_u)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not implemented.");
  }

  //===========================================================================
  bool Disc::isDegenerate(bool& b, bool& r,
			  bool& t, bool& l, double tolerance) const
  //===========================================================================
  {
    b = false;
    r = false;
    t = true;
    l = false;
    return true;
  }


  //===========================================================================
  void Disc::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
  //===========================================================================
  {
    MESSAGE("getDegenerateCorners() not implemented.");
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
    MESSAGE("Not implemented.");
    return NULL;
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
  void Disc::setParameterDomain(double from_upar, double from_vpar,
				double to_upar, double to_vpar)
  //===========================================================================
  {
    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
  }


  //===========================================================================
  void Disc::setDefaultDomain()
  //===========================================================================
  {
    setParameterDomain(0.0, 0.0, radius_, 2.0 * M_PI);
  }


  //===========================================================================
  Circle Disc::boundaryCircle() const
  //===========================================================================
  {
    return Circle(radius_, centre_, z_axis_, x_axis_);
  }




} // namespace Go
