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

#include "GoTools/intersections/PlaneInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>


using std::vector;
using std::string;


namespace Go {


//===========================================================================
PlaneInt::PlaneInt()
  : AlgObj3DInt(1)
//===========================================================================
{
}


//===========================================================================
PlaneInt::PlaneInt(Point point, Point normal)
    : AlgObj3DInt(1), point_(point), normal_(normal)
//===========================================================================
{
    // All values in terms_ initialized to 0.0.
    double a, b, c, d;
    computeConstants(point_, normal_, a, b, c, d);
    terms_.push_back(Alg3DElem(a, 1, 0, 0));
    terms_.push_back(Alg3DElem(b, 0, 1, 0));
    terms_.push_back(Alg3DElem(c, 0, 0, 1));
    terms_.push_back(Alg3DElem(d, 0, 0, 0));
    power_basis_ = true;
}


//===========================================================================
PlaneInt::PlaneInt(double a, double b, double c, double d)
    : AlgObj3DInt(1)
//===========================================================================
{
    // All values in terms_ initialized to 0.0.
    terms_.push_back(Alg3DElem(a, 1, 0, 0));
    terms_.push_back(Alg3DElem(b, 0, 1, 0));
    terms_.push_back(Alg3DElem(c, 0, 0, 1));
    terms_.push_back(Alg3DElem(d, 0, 0, 0));
    computePoints();
    power_basis_ = true;
}


//===========================================================================
PlaneInt::~PlaneInt()
//===========================================================================
{
}


//===========================================================================
void PlaneInt::read(std::istream& is)
//===========================================================================
{
  // Expecting input to be on form
  // plane_pt
  // normal
  // @@sbr Possibly separator characters ... As well as comments ...
  Point plane_pt(3), normal(3);
  char character;
  is >> character;
  plane_pt.read(is);
  is >> character;
  normal.read(is);

  point_ = plane_pt;
  normal_ = normal;

    double a, b, c, d;
    computeConstants(point_, normal_, a, b, c, d);
    terms_.push_back(Alg3DElem(a, 1, 0, 0));
    terms_.push_back(Alg3DElem(b, 0, 1, 0));
    terms_.push_back(Alg3DElem(c, 0, 0, 1));
    terms_.push_back(Alg3DElem(d, 0, 0, 0));
}


//===========================================================================
double PlaneInt::a() const
//===========================================================================
{
    double a = 0.0;
    int sum_order;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	sum_order = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1] + terms_[ki].degrees_[2];
	ASSERT(sum_order <= 1);
	if (terms_[ki].degrees_[0] == 1)
	    a += terms_[ki].factor_;
    }

    return a;
}


//===========================================================================
double PlaneInt::b() const
//===========================================================================
{
    double b = 0.0;
    int sum_order;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	sum_order = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1] + terms_[ki].degrees_[2];
	ASSERT(sum_order <= 1);
	if (terms_[ki].degrees_[1] == 1)
	    b += terms_[ki].factor_;
    }

    return b;
}


//===========================================================================
double PlaneInt::c() const
//===========================================================================
{
    double c = 0.0;
    int sum_order;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	sum_order = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1] + terms_[ki].degrees_[2];
	ASSERT(sum_order <= 1);
	if (terms_[ki].degrees_[2] == 1)
	    c += terms_[ki].factor_;
    }

    return c;
}


//===========================================================================
double PlaneInt::d() const
//===========================================================================
{
    double d = 0.0;
    int sum_order;
    for (size_t ki = 0; ki < terms_.size(); ++ki) {
	sum_order = terms_[ki].degrees_[0] +
	    terms_[ki].degrees_[1] + terms_[ki].degrees_[2];
	ASSERT(sum_order <= 1);
	if (sum_order == 0)
	    d += terms_[ki].factor_;
    }

    return d;
}


//===========================================================================
shared_ptr<SplineSurface>
PlaneInt::surface(Point mid_pt,
		  double length_x, double length_y) const
//===========================================================================
{
    // We create a plane with center in origo and normal the z-axis.
    //   int order = 2;
    vector<double> knots(4, 0.0);
    for (int ki = 2; ki < 4; ++ki)
	knots[ki] = knots[ki] = 1.0;
    int dim = 3;
    vector<double> coefs(4*dim, 0.0);
    for (int ki = 0; ki < 2; ++ki) {
	for (int kj = 0; kj < 2; ++kj) {
	    coefs[dim*(2*ki+kj)] = (ki < 1) ? -0.5*length_x : 0.5*length_x;
	    coefs[dim*(2*ki+kj)+1] = (kj < 1) ? -0.5*length_y : 0.5*length_y;
	    coefs[dim*(2*ki+kj)+2] = 0.0;
	}
    }

    shared_ptr<SplineSurface> go_plane
	(new SplineSurface(2, 2, 2, 2,
			   knots.begin(), knots.begin(),
			   coefs.begin(), 3));

    // As mid_pt may lie outside plane we project.
    Point proj_mid_pt = projectPoint(mid_pt);
    // We then rotate and translate to coincide with values in
    // proj_mid_pt & normal_.
    Point normal(0.0, 0.0, 1.0);
    Point rot_axis = normal%normal_;
    double zero = 1e-14;
    if (rot_axis.length() > zero) { // Otherwise the sf needs no rotation.
	// We must compute the angle (in radians) between normals.
	double angle = normal.angle(normal_);
	GeometryTools::rotateSplineSurf(rot_axis, angle, *go_plane);
    }
    GeometryTools::translateSplineSurf(proj_mid_pt, *go_plane);

    return go_plane;
}


//===========================================================================
void PlaneInt::computePoints()
//===========================================================================
{
    double a_fac = a();
    double b_fac = b();
    double c_fac = c();
    double d_fac = d();
    // Setting x = y = 0.0 we get our ref point.
    point_ = (c_fac == 0.0)
	? Point(0.0, 0.0, 0.0)
	: Point(0.0, 0.0, -d_fac/c_fac);

    normal_ = Point(a_fac, b_fac, c_fac);
    normal_.normalize();
}


//===========================================================================
void PlaneInt::computeConstants(Point point, Point normal,
				double& a, double& b,
				double& c, double& d) const
//===========================================================================
{
    // n*(x - x_0, y - y_0, z - z_0) = 0.
    a = normal_[0];
    b = normal_[1];
    c = normal_[2];
    d = -(normal_*point_);
}



//===========================================================================
Point PlaneInt::projectPoint(const Point& space_pt) const
//===========================================================================
{
    double a, b, c, d;
    computeConstants(point_, normal_, a, b, c, d);
    double signed_dist
	= (a*space_pt[0] + b*space_pt[1] + c*space_pt[2] + d)
	/ (normal_.length());

    Point proj_space_pt = space_pt - signed_dist*normal_;
    return proj_space_pt;
}


} // namespace Go
