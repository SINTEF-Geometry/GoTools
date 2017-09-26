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

#include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/trivariate/SplineVolume.h"


using std::vector;
using std::endl;


namespace Go
{


  // Constructor
  //===========================================================================
  Parallelepiped::Parallelepiped(Point corner,
				 Point dir_u, Point dir_v, Point dir_w,
				 double len_u, double len_v, double len_w) :
    corner_(corner),
    dir_u_(dir_u), dir_v_(dir_v), dir_w_(dir_w),
    length_u_(len_u), length_v_(len_v), length_w_(len_w)
  //===========================================================================
  {
    if (corner_.dimension() != 3 ||
	dir_u_.dimension() != 3 ||
	dir_v_.dimension() != 3 ||
	dir_w_.dimension() != 3)
      {
	THROW("Dimension must be 3.");
	return;
      }

    dir_u_.normalize();
    dir_v_.normalize();
    dir_w_.normalize();
  }

  // Destructor
  //===========================================================================
  Parallelepiped::~Parallelepiped()
  //===========================================================================
  {
  }

  //===========================================================================
  void Parallelepiped::read(std::istream& is)
  //===========================================================================
  {
    if (!is.good()) {
	THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    if (dim != 3)
	THROW("Dimension must be 3.");
    corner_.resize(dim);
    dir_u_.resize(dim);
    dir_v_.resize(dim);
    dir_w_.resize(dim);

    is >> corner_
       >> dir_u_
       >> length_u_
       >> dir_v_
       >> length_v_
       >> dir_w_
       >> length_w_;
  }

  //===========================================================================
  void Parallelepiped::write(std::ostream& os) const
  //===========================================================================
  {
    os << dimension() << endl
       << corner_ << endl
       << dir_u_ << endl
       << length_u_ << endl
       << dir_v_ << endl
       << length_v_ << endl
       << dir_w_ << endl
       << length_w_ << endl;
  }

  //===========================================================================
  int Parallelepiped::dimension() const
  //===========================================================================
  {
    return corner_.dimension();
  }

  //===========================================================================
  ClassType Parallelepiped::instanceType() const
  //===========================================================================
  {
    return classType();
  }


  //===========================================================================
  BoundingBox Parallelepiped::boundingBox() const
  //===========================================================================
  {
    MESSAGE("boundingBox() not yet implemented");
    BoundingBox bb;
    return bb;
  }

  //===========================================================================
  DirectionCone Parallelepiped::tangentCone(int pardir) const
  //===========================================================================
  {
    if (pardir == 0)
      return DirectionCone(dir_u_);
    else if (pardir == 1)
      return DirectionCone(dir_v_);
    else
      return DirectionCone(dir_w_);
  }

  //===========================================================================
  const Array<double,6> Parallelepiped::parameterSpan() const
  //===========================================================================
  {
    Array<double,6> pSpan;
    pSpan[0] = pSpan[2] = pSpan[4] = 0.0;
    pSpan[1] = length_u_;
    pSpan[3] = length_v_;
    pSpan[5] = length_w_;
    return pSpan;
  }

  //===========================================================================
  void Parallelepiped::point(Point& pt, double upar, double vpar, double wpar) const
  //===========================================================================
  {
    pt = corner_ + upar * dir_u_ + vpar * dir_v_ + wpar * dir_w_;
  }

  //===========================================================================
  void Parallelepiped::point(vector<Point>& pts, 
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
    pts[1] = length_u_ * dir_u_;
    pts[2] = length_v_ * dir_v_;
    pts[3] = length_w_ * dir_w_;

    // Second order and higher derivatives vanish. They are already
    // set to zero, so we return.
  }

  //===========================================================================
  double Parallelepiped::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    MESSAGE("nextSegmentVal() not yet implemented");
    return 0.0;
  }

  //===========================================================================
  void Parallelepiped::closestPoint(const Point& pt,
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
  void Parallelepiped::reverseParameterDirection(int pardir)
  //===========================================================================
  {
    if (pardir == 0)
      {
	corner_ += length_u_ * dir_u_;
	dir_u_ = -dir_u_;
      }
    else if (pardir == 1)
      {
	corner_ += length_v_ * dir_v_;
	dir_v_ = -dir_v_;
      }
    else if (pardir == 2)
      {
	corner_ += length_w_ * dir_w_;
	dir_w_ = -dir_w_;
      }
  }

  //===========================================================================
  void Parallelepiped::swapParameterDirection(int pardir1, int pardir2)
  //===========================================================================
  {
    if ((pardir1 == 0 && pardir2 == 1) || (pardir1 == 1 && pardir2 == 0))
      {
	double tmp_l = length_u_;
	length_u_ = length_v_;
	length_v_ = tmp_l;
	Point tmp_d = dir_u_;
	dir_u_ = dir_v_;
	dir_v_ = tmp_d;
      }
    else if ((pardir1 == 0 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 0))
      {
	double tmp_l = length_u_;
	length_u_ = length_w_;
	length_w_ = tmp_l;
	Point tmp_d = dir_u_;
	dir_u_ = dir_w_;
	dir_w_ = tmp_d;
      }
    else if ((pardir1 == 1 && pardir2 == 2) || (pardir1 == 2 && pardir2 == 1))
      {
	double tmp_l = length_v_;
	length_v_ = length_w_;
	length_w_ = tmp_l;
	Point tmp_d = dir_v_;
	dir_v_ = dir_w_;
	dir_w_ = tmp_d;
      }
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > 
  Parallelepiped::getAllBoundarySurfaces(bool do_clear) const
  //===========================================================================
  {
    MESSAGE("getAllBoundarySurfaces() not implemented.");
    vector<shared_ptr<ParamSurface> > bound_surf;
    return bound_surf;
  }

  //===========================================================================
  void Parallelepiped::translate(const Point& vec)
  //===========================================================================
  {
    ALWAYS_ERROR_IF(dimension() != vec.dimension(), "Volume and translation vector of different dimension");
    corner_ += vec;
  }

  //===========================================================================
  SplineVolume* Parallelepiped::geometryVolume() const
  //===========================================================================
  {
    vector<Point> allCorners(8);
    allCorners[0] = corner_;
    allCorners[1] = corner_ + dir_u_ * length_u_;
    Point v_edge = dir_v_ * length_v_;
    for (int i = 0; i < 2; ++i)
      allCorners[2+i] = allCorners[i] + v_edge;
    Point w_edge = dir_w_ * length_w_;
    for (int i = 0; i < 4; ++i)
      allCorners[4+i] = allCorners[i] + w_edge;

    vector<double> coefs(24);
    for (int i = 0, pos = 0; i < 8; ++i)
      for (int j = 0; j< 3; ++j, ++pos)
	coefs[pos] = allCorners[i][j];

    vector<double> knots_u(4), knots_v(4), knots_w(4);
    knots_u[0] = knots_u[1] = knots_v[0] = knots_v[1] = knots_w[0] = knots_w[1] = 0.0;
    knots_u[2] = knots_u[3] = length_u_;
    knots_v[2] = knots_v[3] = length_v_;
    knots_w[2] = knots_w[3] = length_w_;

    return new SplineVolume(2, 2, 2, 2, 2, 2, knots_u.begin(), knots_v.begin(), knots_w.begin(), coefs.begin(), 3);
  }



} // namespace Go
