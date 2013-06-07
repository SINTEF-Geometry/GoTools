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

#include "GoTools/trivariate/ConeVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <limits>


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
  ConeVolume::ConeVolume(double radius, Point location,
			 Point z_axis, Point x_axis,
			 double cone_angle) :
    radius_(radius),
    location_(location), z_axis_(z_axis), x_axis_(x_axis),
    cone_angle_(cone_angle),
    height_min_(-numeric_limits<double>::infinity()),
    height_max_(numeric_limits<double>::infinity()),
    centre_degen_(false)
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
  ConeVolume::~ConeVolume()
  //===========================================================================
  {
  }

  //===========================================================================
  void ConeVolume::read(istream& is)
  //===========================================================================
  {
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim, is_infinite_int, centre_degen_int;
    is >> dim;
    location_.resize(dim);
    z_axis_.resize(dim);
    x_axis_.resize(dim);
    is >> radius_
       >> location_
       >> z_axis_
       >> x_axis_
       >> cone_angle_
       >> height_min_;

    is >> is_infinite_int;
    if (is_infinite_int == 1)
      height_max_ = -numeric_limits<double>::infinity();
    else
      is >> height_max_;

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
  void ConeVolume::write(ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl
       << cone_angle_ << endl
       << height_min_ << endl;

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
  int ConeVolume::dimension() const
  //===========================================================================
  {
    return location_.dimension();
  }

  //===========================================================================
  ClassType ConeVolume::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox ConeVolume::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }


  //===========================================================================
  ConeVolume* ConeVolume::clone() const
  //===========================================================================
  {
    ConeVolume* newCone = new ConeVolume(radius_, location_, z_axis_, x_axis_, cone_angle_);
    newCone->height_min_ = height_min_;
    newCone->height_max_ = height_max_;
    newCone->centre_degen_ = centre_degen_;
    for (int i = 0; i < 4; ++i)
      newCone->degen_angles_[i] = degen_angles_[i];
    return newCone;
  }


  //===========================================================================
  DirectionCone ConeVolume::tangentCone(int pardir) const
  //===========================================================================
  {
    MESSAGE("tangentCone() not yet implemented");
    DirectionCone bb;
    return bb;
  }

  //===========================================================================
  const Array<double,6> ConeVolume::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;

    // First direction: Radius
    pSpan[0] = 0.0;
    pSpan[3] = radius_;

    // Second direction: Angle
    pSpan[1] = 0.0;
    pSpan[4] = 2 * M_PI;

    // Third direction: Height
    pSpan[2] = height_min_;
    pSpan[5] = height_max_;

    return pSpan;
  }

  //===========================================================================
  void ConeVolume::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = location_
      + upar * (radius_ + wpar * tan(cone_angle_)) * (cos(vpar) * x_axis_ 
							+ sin(vpar) * y_axis_)
      + wpar * z_axis_;
  }

  //===========================================================================
  void ConeVolume::point(vector<Point>& pts, 
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
    pts[1] = (radius_ + wpar * tan(cone_angle_))
      * (cos(vpar) * x_axis_ + sin(vpar) * y_axis_);
    pts[2] = upar * (radius_ + wpar * tan(cone_angle_))
      * (-sin(vpar) * x_axis_ + cos(vpar) * y_axis_);
    pts[3] = upar * tan(cone_angle_)
      * (cos(vpar) * x_axis_ + sin(vpar) * y_axis_)
      + z_axis_;

    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");
  }

  //===========================================================================
  double ConeVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void ConeVolume::closestPoint(const Point& pt,
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
  void ConeVolume::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not implemented.");
  }

  //===========================================================================
  void ConeVolume::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not implemented.");
  }

  //===========================================================================
  vector<shared_ptr<ParamSurface> > ConeVolume::getAllBoundarySurfaces() const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void ConeVolume::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    location_ += vec;
  }

  //===========================================================================
  SplineVolume* ConeVolume::geometryVolume() const
  //===========================================================================
  {
    // First create bilinear surface when v = 0
    vector<Point> corner(4);
    corner[0] = location_ + height_min_ * z_axis_;
    corner[1] = location_
      + (radius_ + height_min_ * tan(cone_angle_)) * (x_axis_ )
      + height_min_ * z_axis_;
    corner[2] = location_ + height_max_ * z_axis_;
    corner[3] = location_
      + (radius_ + height_max_ * tan(cone_angle_)) * (x_axis_ )
      + height_max_ * z_axis_;
    vector<double> coefs_section(12);
    for (int i = 0, pos = 0; i < 4; ++i)
      for (int j = 0; j < 3; ++j, ++pos)
	coefs_section[pos] = corner[i][j];
    vector<double> knots_u(4), knots_w(4);
    knots_u[0] = knots_u[1] = 0.0;
    knots_u[2] = knots_u[3] = 1.0;
    knots_w[0] = knots_w[1] = height_min_;
    knots_w[2] = knots_w[3] = height_max_;
    SplineSurface section(2, 2, 2, 2, knots_u.begin(), knots_w.begin(), coefs_section.begin(), 3);

    // Then rotate
    SplineVolume* result = SweepVolumeCreator::rotationalSweptVolume(section, 2 * M_PI, location_, z_axis_);
    result->swapParameterDirection(0, 1);

    return result;
  }


  //===========================================================================
  void ConeVolume::setParameters(double from_par, double to_par, int pardir)
  //===========================================================================
  {
    if (pardir == 2)
      {
	height_min_ = from_par;
	height_max_ = to_par;
      }
  }


  //===========================================================================
  void ConeVolume::setCoordinateAxes()
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
