//===========================================================================
//
// File : CylinderVolume.C
//
// Created: Tue Nov 10 08:11:24 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================




#include "GoTools/trivariate/CylinderVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <limits>


using std::shared_ptr;
using std::min;
using std::max;
using std::numeric_limits;
using std::istream;
using std::ostream;
using std::vector;
using std::endl;

namespace Go
{

  // Constructor
  //===========================================================================
  CylinderVolume::CylinderVolume(Point centre, double radius, Point normal, Point x_axis) :
    centre_(centre),
    radius_min_(0.0),
    radius_max_(radius),
    height_min_(-numeric_limits<double>::infinity()),
    height_max_(numeric_limits<double>::infinity()),
    x_axis_(x_axis),
    z_axis_(normal),
    centre_degen_(false)
  //===========================================================================
  {
    if (centre_.dimension() != 3) {
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
  CylinderVolume::~CylinderVolume()
  //===========================================================================
  {
  }

  //===========================================================================
  void CylinderVolume::read(istream& is)
  //===========================================================================
  {
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim, is_infinite_int, centre_degen_int;
    is >> dim;
    centre_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> centre_
       >> z_axis_
       >> x_axis_
       >> radius_min_
       >> radius_max_;

    is >> is_infinite_int;
    if (is_infinite_int == 1)
      height_min_ = -numeric_limits<double>::infinity();
    else
      is >> height_min_;

    is >> is_infinite_int;
    if (is_infinite_int == 1)
      height_min_ = -numeric_limits<double>::infinity();
    else
      is >> height_min_;

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
  void CylinderVolume::write(ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << centre_ << endl
       << z_axis_ << endl
       << x_axis_ << endl
       << radius_min_ << endl
       << radius_max_ << endl;

    if (height_min_ == -numeric_limits<double>::infinity())
      os << "1" << endl;
    else
      os << "0 " << height_min_ << endl;

    if (height_max_ == numeric_limits<double>::infinity())
      os << "1" << endl;
    else
      os << "0 " << height_max_ << endl;

    if (centre_degen_)
      os << "1" << endl;
    else
      os << "0" << endl;
    for (int i = 0; i < 4; ++i)
      os << degen_angles_[i] << endl;
  }

  //===========================================================================
  int CylinderVolume::dimension() const
  //===========================================================================
  {
    return centre_.dimension();
  }

  //===========================================================================
  ClassType CylinderVolume::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox CylinderVolume::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }


  //===========================================================================
  CylinderVolume* CylinderVolume::clone() const
  //===========================================================================
  {
    CylinderVolume* cyl = new CylinderVolume(centre_, radius_max_, z_axis_, x_axis_);
    cyl->radius_min_ = radius_min_;
    cyl->height_min_ = height_min_;
    cyl->height_max_ = height_max_;
    cyl->centre_degen_ = centre_degen_;
    for (int i = 0; i < 4; ++i)
      cyl->degen_angles_[i] = degen_angles_[i];
    return cyl;
  }


  //===========================================================================
  DirectionCone CylinderVolume::tangentCone(int pardir) const
  //===========================================================================
  {
    MESSAGE("tangentCone() not yet implemented");
    DirectionCone bb;
    return bb;
  }

  //===========================================================================
  const Array<double,6> CylinderVolume::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;

    // First direction: Radius
    pSpan[0] = radius_min_;
    pSpan[3] = radius_max_;

    // Second direction: Angle
    pSpan[1] = 0.0;
    pSpan[4] = 2 * M_PI;

    // Third direction: Height
    pSpan[2] = height_min_;
    pSpan[5] = height_max_;

    return pSpan;
  }

  //===========================================================================
  void CylinderVolume::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = centre_ + upar * (cos(vpar) * x_axis_ + sin(vpar) * y_axis_) + wpar * z_axis_;
  }

  //===========================================================================
  void CylinderVolume::point(vector<Point>& pts, 
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

    point(pts[0], upar, vpar, wpar);
    if (derivs == 0)
        return;

    // Derivatives
    pts[1] = cos(vpar) * x_axis_ + sin(vpar) * y_axis_;
    pts[2] = upar * (-sin(vpar) * x_axis_ + cos(vpar) * y_axis_);
    pts[3] = z_axis_;

    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");
  }

  //===========================================================================
  double CylinderVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void CylinderVolume::closestPoint(const Point& pt,
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
  void CylinderVolume::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not implemented.");
  }

  //===========================================================================
  void CylinderVolume::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not implemented.");
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > CylinderVolume::getAllBoundarySurfaces() const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void CylinderVolume::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    centre_ += vec;
  }

  //===========================================================================
  SplineVolume* CylinderVolume::geometryVolume() const
  //===========================================================================
  {
    // Create volume by starting with line on bottom face from inner to outer surface,
    // then rotate to get bottom disk, then use linear sweep along height to get entire volume
    Point inner_bottom = centre_ + radius_min_ * x_axis_ + height_min_ * z_axis_;
    Point outer_bottom = inner_bottom + (radius_max_ - radius_min_) * x_axis_;
    Point outer_top = outer_bottom + (height_max_ - height_min_) * z_axis_;

    SplineCurve bottom_line(inner_bottom, radius_min_, outer_bottom, radius_max_);
    SplineCurve outer_line(outer_bottom, height_min_, outer_top, height_max_);

    shared_ptr<SplineSurface> bottom_disk(SweepSurfaceCreator::rotationalSweptSurface(bottom_line, 2 * M_PI, centre_, z_axis_));
    SplineSurface* bottom_disk_ptr = bottom_disk.get();
    SplineVolume* result = SweepVolumeCreator::linearSweptVolume(*bottom_disk_ptr, outer_line, outer_bottom);

    result->swapParameterDirection(0, 1);

    return result;
  }


  //===========================================================================
  void CylinderVolume::setParameters(double from_par, double to_par, int pardir)
  //===========================================================================
  {
    if (pardir == 0)
      {
	radius_min_ = from_par;
	radius_max_ = to_par;
	if (from_par > 0.0)
	  centre_degen_ = true;
      }
    else if (pardir == 2)
      {
	height_min_ = from_par;
	height_max_ = to_par;
      }
  }


  //===========================================================================
  void CylinderVolume::setCoordinateAxes()
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
