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

#include "GoTools/trivariate/SphereVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"


using std::vector;
using std::endl;


namespace Go
{


  // Constructor
  //===========================================================================
  SphereVolume::SphereVolume(double radius, Point location, Point z_axis, Point x_axis) :
    radius_(radius),
    location_(location), z_axis_(z_axis), x_axis_(x_axis)
  //===========================================================================
  {
    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }
    setCoordinateAxes();
  }

  // Destructor
  //===========================================================================
  SphereVolume::~SphereVolume()
  //===========================================================================
  {
  }

  //===========================================================================
  void SphereVolume::read(std::istream& is)
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
  }

  //===========================================================================
  void SphereVolume::write(std::ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << radius_ << endl
       << location_ << endl
       << z_axis_ << endl
       << x_axis_ << endl;
  }

  //===========================================================================
  int SphereVolume::dimension() const
  //===========================================================================
  {
    return location_.dimension();
  }

  //===========================================================================
  ClassType SphereVolume::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox SphereVolume::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }

  //===========================================================================
  DirectionCone SphereVolume::tangentCone(int pardir) const
  //===========================================================================
  {
    MESSAGE("tangentCone() not yet implemented");
    DirectionCone bb;
    return bb;
  }

  //===========================================================================
  const Array<double,6> SphereVolume::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;

    // First direction: Distance from centre
    pSpan[0] = 0.0;
    pSpan[3] = radius_;

    // Second direction: Azimuth angle
    pSpan[1] = 0.0;
    pSpan[4] = 2 * M_PI;

    // Third direction: Elevation angle
    pSpan[2] = -0.5 * M_PI;
    pSpan[5] = -0.5 * M_PI;

    return pSpan;
  }

  //===========================================================================
  void SphereVolume::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = location_ 
	+ upar * (cos(wpar) * (cos(vpar) * x_axis_ + sin(vpar) * y_axis_)
		     + sin(wpar) * z_axis_);
  }

  //===========================================================================
  void SphereVolume::point(vector<Point>& pts, 
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
    pts[1] = cos(wpar) * (cos(vpar) * x_axis_ + sin(vpar) * y_axis_)
             + sin(wpar) * z_axis_;
    pts[2] = upar * cos(wpar) * (-sin(vpar) * x_axis_ + cos(vpar) * y_axis_);
    pts[3] = upar * (-sin(wpar) * (cos(vpar) * x_axis_ 
				   + sin(vpar) * y_axis_) 
		     + cos(wpar) * z_axis_);

    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");
  }

  //===========================================================================
  double SphereVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void SphereVolume::closestPoint(const Point& pt,
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
  void SphereVolume::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    MESSAGE("reverseParameterDirection() not implemented.");
  }

  //===========================================================================
  void SphereVolume::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    MESSAGE("swapParameterDirection() not implemented.");
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > 
  SphereVolume::getAllBoundarySurfaces(bool do_clear) const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void SphereVolume::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    location_ += vec;
  }

  //===========================================================================
  SplineVolume* SphereVolume::geometryVolume() const
  //===========================================================================
  {
    // Start with line from centre to north pole (when v = 0.0, w = -PI/2)
    SplineCurve line(location_, 0.0, location_ - radius_ * z_axis_, radius_);

    // Rotatet to create half disk (when v = 0.0)
    shared_ptr<SplineSurface> section(SweepSurfaceCreator::rotationalSweptSurface(line, M_PI, location_, -y_axis_));
    SplineSurface* section_ptr = section.get();

    // Rotate again to create entire spline volume
    SplineVolume* result = SweepVolumeCreator::rotationalSweptVolume(*section_ptr, 2 * M_PI, location_, z_axis_);

    // Reorder parameters and parameter domains
    result->swapParameterDirection(1, 2);
    result->swapParameterDirection(0, 1);
    result->setParameterDomain(result->startparam(0), result->endparam(0),
			       result->startparam(2), result->endparam(2),
			       -M_PI * 0.5, M_PI * 0.5);

    return result;
  }


  //===========================================================================
  void SphereVolume::setCoordinateAxes()
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
