//===========================================================================
//
// File : TorusVolume.C
//
// Created: Wed Nov 11 08:26:22 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#include "GoTools/trivariate/TorusVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <limits>


using std::vector;
using std::endl;
using std::shared_ptr;


namespace Go
{


// Constructor
//===========================================================================
TorusVolume::TorusVolume(double major_radius, double minor_radius,
        Point location, Point z_axis, Point x_axis) :
    location_(location),
    major_radius_(major_radius),
    minor_radius_min_(0.0), minor_radius_max_(minor_radius),
    rev_angle_min_(0.0), rev_angle_max_(2 * M_PI),
    z_axis_(z_axis), x_axis_(x_axis)
//===========================================================================
{
    if (location_.dimension() != 3) {
        THROW("Dimension must be 3.");
        return;
    }

    degen_angles_[0] = 0.0;
    degen_angles_[1] = 0.5 * M_PI;
    degen_angles_[2] = 1.0 * M_PI;
    degen_angles_[3] = 1.5 * M_PI;

    setCoordinateAxes();
}

  // Destructor
  //===========================================================================
  TorusVolume::~TorusVolume()
  //===========================================================================
  {
  }

  //===========================================================================
  void TorusVolume::read(std::istream& is)
  //===========================================================================
  {
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim, centre_degen_int;
    is >> dim;
    location_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> location_
       >> z_axis_
       >> x_axis_
       >> major_radius_
       >> minor_radius_min_
       >> minor_radius_max_
       >> rev_angle_min_
       >> rev_angle_max_;

    is >> centre_degen_int;
    for (int i = 0; i < 4; ++i)
      is >> degen_angles_[i];

    if (centre_degen_int == 0)
      centre_degen_ = false;
    else if (centre_degen_int == 1)
      centre_degen_ = true;
    else
      THROW("Unknown input for centre_degen - must be 0 or 1");

    setCoordinateAxes();
  }

  //===========================================================================
  void TorusVolume::write(std::ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl
       << major_radius_ << endl
       << minor_radius_min_ << endl
       << minor_radius_max_ << endl
       << rev_angle_min_ << endl
       << rev_angle_max_ << endl;

    if (centre_degen_)
      os << "1" << endl;
    else
      os << "0" << endl;
    for (int i = 0; i < 4; ++i)
      os << degen_angles_[i] << endl;
  }

  //===========================================================================
  int TorusVolume::dimension() const
  //===========================================================================
  {
    return location_.dimension();
  }

  //===========================================================================
  ClassType TorusVolume::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox TorusVolume::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }


  //===========================================================================
  TorusVolume* TorusVolume::clone() const
  //===========================================================================
  {
    TorusVolume* tor = new TorusVolume(major_radius_, minor_radius_max_,
				       location_, z_axis_, x_axis_);
    tor->minor_radius_min_ = minor_radius_min_;
    tor->rev_angle_min_ = rev_angle_min_;
    tor->rev_angle_max_ = rev_angle_max_;
    tor->centre_degen_ = centre_degen_;
    for (int i = 0; i < 4; ++i)
      tor->degen_angles_[i] = degen_angles_[i];
    return tor;
  }


  //===========================================================================
  DirectionCone TorusVolume::tangentCone(int pardir) const
  //===========================================================================
  {
    MESSAGE("tangentCone() not yet implemented");
    DirectionCone bb;
    return bb;
  }

  //===========================================================================
  const Array<double,6> TorusVolume::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;

    // First direction: Minor radius
    pSpan[0] = minor_radius_min_;
    pSpan[3] = minor_radius_max_;

    // Second direction: Rotation angle around z-axis
    pSpan[1] = rev_angle_min_;
    pSpan[4] = rev_angle_max_;

    // Third direction: Rotation angle around circle
    pSpan[2] = 0.0;
    pSpan[5] = 2 * M_PI;

    return pSpan;
  }

  //===========================================================================
  void TorusVolume::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = location_
      + (major_radius_ + upar * cos(wpar)) * (cos(vpar) * x_axis_ 
					      + sin(vpar) * y_axis_)
      + upar * sin(wpar) * z_axis_;
  }

  //===========================================================================
  void TorusVolume::point(vector<Point>& pts, 
			     double upar, double vpar, double wpar,
			     int derivs,
			     bool u_from_right,
			     bool v_from_right,
			     bool w_from_right,
			     double resolution ) const
  //===========================================================================
  {
    DEBUG_ERROR_IF(derivs < 0,
		   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
    int ptsz;
    ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz< totpts,
		   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
        }
	pts[i].setValue(0.0);
    }

    point(pts[0], upar, vpar, wpar);
    if (derivs == 0)
        return;

    // Derivatives
    double cosv = cos(vpar);
    double sinv = sin(vpar);
    double cosw = cos(wpar);
    double sinw = sin(wpar);

    pts[1] = cosw * (cosv * x_axis_ + sinv * y_axis_) + sinw * z_axis_;
    pts[2] = (major_radius_ + upar * cosw) * (-sinv * x_axis_ + cosv * y_axis_);
    pts[3] = upar * (-sinw * (cosv * x_axis_ + sinv * y_axis_)
		     + cosw * z_axis_);

    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");
  }

  //===========================================================================
  double TorusVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void TorusVolume::closestPoint(const Point& pt,
				    double&        clo_u,
				    double&        clo_v,
				    double&        clo_w,
				    Point&         clo_pt,
				    double&        clo_dist,
				    double         epsilon,
				    double   *seed) const
  //===========================================================================
  {
    MESSAGE("closestPoint() not yet implemented");
  }

  //===========================================================================
  void TorusVolume::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not implemented.");
  }

  //===========================================================================
  void TorusVolume::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not implemented.");
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > TorusVolume::getAllBoundarySurfaces() const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void TorusVolume::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    location_ += vec;
  }

  //===========================================================================
  SplineVolume* TorusVolume::geometryVolume() const
  //===========================================================================
  {
    double cos_v_min = cos(rev_angle_min_);
    double sin_v_min = sin(rev_angle_min_);
    Point rev_x_axis = cos_v_min * x_axis_ + sin_v_min * y_axis_;
    Point rev_y_axis = cos_v_min * y_axis_ - sin_v_min * x_axis_;

    // Start with line from inner to outer radius in tube section (when v = v_min, w = 0.0)
    SplineCurve line(location_ + (major_radius_ + minor_radius_min_) * rev_x_axis,
		     minor_radius_min_,
		     location_ + (major_radius_ + minor_radius_max_) * rev_x_axis,
		     minor_radius_max_);

    // Rotatet to create entire tube section (when v = v_min)
    shared_ptr<SplineSurface> section(SweepSurfaceCreator::rotationalSweptSurface(line,
										  2 * M_PI,
										  location_ + major_radius_ * rev_x_axis,
										  -rev_y_axis));
    SplineSurface* section_ptr = section.get();

    // Rotate again to create entire spline volume
    SplineVolume* result = SweepVolumeCreator::rotationalSweptVolume(*section_ptr, rev_angle_max_ - rev_angle_min_, location_, z_axis_);

    // Reorder parameters and parameter domains
    result->swapParameterDirection(1, 2);
    result->swapParameterDirection(0, 1);
    result->setParameterDomain(result->startparam(0), result->endparam(0),
			       rev_angle_min_, rev_angle_max_,
			       result->startparam(2), result->endparam(2));

    return result;
  }


  //===========================================================================
  void TorusVolume::setParameters(double from_par, double to_par, int pardir)
  //===========================================================================
  {
    if (pardir == 0)
      {
	minor_radius_min_ = from_par;
	minor_radius_max_ = to_par;
	if (from_par > 0.0)
	  centre_degen_ = true;
      }
    else if (pardir == 1)
      {
	rev_angle_min_ = from_par;
	rev_angle_max_ = to_par;
      }
  }


  //===========================================================================
  void TorusVolume::setCoordinateAxes()
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


} // namespace Go
