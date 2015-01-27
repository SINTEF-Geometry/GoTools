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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ElementarySurface.h"
#include <algorithm>
#include <iomanip>
#include <fstream>


using std::vector;
using std::streamsize;
using std::endl;
using std::pair;
using std::make_pair;


namespace Go
{


//===========================================================================
SplineSurface::~SplineSurface()

//===========================================================================
{
}

//===========================================================================
void SplineSurface::read (std::istream& is)
//===========================================================================
{
    // We verify that the object is valid.
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    // Canonical data
    is >> dim_;
    is >> rational_;
    is >> basis_u_;
    is >> basis_v_;
    int nc = basis_u_.numCoefs()*basis_v_.numCoefs();
    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    if (rational_) {
	int n = nc * (dim_ + 1);
	rcoefs_.resize(n);
	for (int i = 0; i < n; ++i)
	    is >> rcoefs_[i];
	coefs_.resize(nc*dim_);
	updateCoefsFromRcoefs();
    } else {
	int n = nc*dim_;
	coefs_.resize(n);
	for (int i = 0; i < n; ++i)
	    is >> coefs_[i];
    }

    is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
}


//===========================================================================
void SplineSurface::write (std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    os << dim_ << ' ' << rational_ << '\n';
    os << basis_u_;
    os << basis_v_;
    int n = basis_u_.numCoefs()*basis_v_.numCoefs();
    int kdim = dim_ + (rational_ ? 1 : 0);
    const vector<double>& co = rational_ ? rcoefs_ : coefs_;
    for (int i = 0; i < n; ++i) {
	os << co[i*kdim];
	for (int d = 1; d < kdim; ++d) {
	    os << ' ' << co[i*kdim + d];
	}
	os << '\n';
    }
    os << endl;
    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
BoundingBox SplineSurface::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    box.setFromArray(&coefs_[0], &coefs_[0] + coefs_.size(), dim_);
    return box;
}

//===========================================================================
CompositeBox SplineSurface::compositeBox() const
//===========================================================================
{
    CompositeBox box(&coefs_[0], dim_, numCoefs_u(), numCoefs_v());
    return box;
}


//===========================================================================
DirectionCone SplineSurface::normalCone(NormalConeMethod method) const
//===========================================================================
{
    ALWAYS_ERROR_IF(dim_ != 3, "Normal only defined in 3D");

    // bool b, r, t, l;
    // double epsilon = 1.0e-8;
    // if (isDegenerate(b, r, t, l, epsilon))
    // 	THROW("Surface is degenerate - no normal cone defined.");

    // Several implementations are made, controlled with the value of
    // the enum 'method':
    //
    // SederbergMeyers: Default. Builds derivative cones in the u and
    // v direction and contructs the normal cone as the smallest cone
    // containing all cross products of vectors from each cone. May
    // give very large cones.
    //
    // SMCornersFirst: Similar, but with corners processed first for
    // the u and v cones.

    if (method == SederbergMeyers) {
	// Method based on T.W. Sederberg and R.J. Meyers, GAGD 5 (1988) 161
	int nu = numCoefs_u();
	int nv = numCoefs_v();
	const double* start = &coefs_[0];
	Point null(0.0, 0.0, 0.0);

	// First make cone in u-direction
	vector<double> coefs_u;
	Point tmp(dim_);
	int i;
	for (i = 0; i < nu-1; ++i) {
	    for (int j = 0; j < nv; ++j) {
		int a = dim_ * (nu*j + i); // (i, j)
		int b = dim_ * (nu*j + i + 1); // (i+1, j)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_u.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_u = (int)coefs_u.size();
	DirectionCone cone_u;
	cone_u.setFromArray(&coefs_u[0], &coefs_u[0] + size_u, dim_);
	if (cone_u.greaterThanPi())
	    return DirectionCone(null, 4.0);

	// Then in v-direction
	vector<double> coefs_v;
	for (i = 0; i < nu; ++i) {
	    for (int j = 0; j < nv-1; ++j) {
		int a = dim_ * (nu*j + i); // (i, j)
		int b = dim_ * (nu*(j+1) + i); // (i, j+1)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_v.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_v = (int)coefs_v.size();
	DirectionCone cone_v;
	cone_v.setFromArray(&coefs_v[0], &coefs_v[0] + size_v, dim_);
	if (cone_v.greaterThanPi())
	    return DirectionCone(null, 4.0);

	// Calculate centre and angle
	Point s = cone_u.centre();
	Point t = cone_v.centre();
	double cos_beta = s.cosAngle(t);
	double beta = acos(cos_beta);
	double theta_s = cone_u.angle();
	double theta_t = cone_v.angle();
	if (theta_s + theta_t >= beta
	    || theta_s + theta_t + beta >= M_PI)
	    return DirectionCone(null, 4.0);

	Point centre = s % t;

	double sin_s = sin(theta_s);
	double sin_t = sin(theta_t);
	double angle = sqrt(sin_s * sin_s + 2.0 * sin_s * sin_t * cos_beta
			    + sin_t * sin_t);
	double sin_beta = sin(beta);
	angle = asin(angle / sin_beta);

	DirectionCone normal_cone(centre, angle);
	return normal_cone;

    } else if (method == SMCornersFirst) {
	// Sederberg-Meyers with corners processed first
	int nu = numCoefs_u();
	int nv = numCoefs_v();
	const double* start = &coefs_[0];
	Point null(0.0, 0.0, 0.0);

	// First make cone in u-direction - push_back the corners first
	vector<double> coefs_u;
	int a, b;
	Point tmp(dim_);
	// Corner (0, 0)
	a = 0;
	b = dim_;
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_u.push_back(tmp[dd]);
	    }
	}
	// Corner (nu-2, 0)
	a = dim_ * (nu-2);
	b = dim_ * (nu-1);
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_u.push_back(tmp[dd]);
	    }
	}
	// Corner (0, nv-1)
	a = dim_ * nu * (nv-1);
	b = dim_ * (1 + nu * (nv-1));
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_u.push_back(tmp[dd]);
	    }
	}
	// Corner (nu-2, nv-1)
	a = dim_ * (nu-2 + nu * (nv-1));
	b = dim_ * (nu-1 + nu * (nv-1));
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_u.push_back(tmp[dd]);
	    }
	}

	int i;
	for (i = 0; i < nu-1; ++i) {
	    for (int j = 0; j < nv; ++j) {
		a = dim_ * (i + nu*j); // (i, j)
		b = dim_ * (i+1 + nu*j); // (i+1, j)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_u.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_u = (int)coefs_u.size();
	DirectionCone cone_u;
	cone_u.setFromArray(&coefs_u[0], &coefs_u[0] + size_u, dim_);
	if (cone_u.greaterThanPi())
	    return DirectionCone(null, 4.0);

	// Then in v-direction
	vector<double> coefs_v;
	// Corner (0, 0)
	a = 0;
	b = dim_ * nu;
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_v.push_back(tmp[dd]);
	    }
	}
	// Corner (nu-1, 0)
	a = dim_ * (nu-1);
	b = dim_ * (nu-1 + nu);
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_v.push_back(tmp[dd]);
	    }
	}
	// Corner (0, nv-2)
	a = dim_ * nu * (nv-2);
	b = dim_ * nu * (nv-1);
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_v.push_back(tmp[dd]);
	    }
	}
	// Corner (nu-1, nv-2)
	a = dim_ * (nu-1 + nu * (nv-2));
	b = dim_ * (nu-1 + nu * (nv-1));
	for (int dd = 0; dd < dim_; ++dd)
	    tmp[dd] = *(start + b + dd) - *(start + a + dd);
	if (tmp.length() > 0.0) { // Ignore vectors of zero length
	    for (int dd = 0; dd < dim_; ++dd) {
		coefs_v.push_back(tmp[dd]);
	    }
	}

	for (i = 0; i < nu; ++i) {
	    for (int j = 0; j < nv-1; ++j) {
		a = dim_ * (i + nu*j); // (i, j)
		b = dim_ * (i + nu*(j+1)); // (i, j+1)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_v.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_v = (int)coefs_v.size();
	DirectionCone cone_v;
	cone_v.setFromArray(&coefs_v[0], &coefs_v[0] + size_v, dim_);
	if (cone_v.greaterThanPi())
	    return DirectionCone(null, 4.0);

	// Calculate centre and angle
	Point s = cone_u.centre();
	Point t = cone_v.centre();
	double cos_beta = s.cosAngle(t);
	double beta = acos(cos_beta);
	double theta_s = cone_u.angle();
	double theta_t = cone_v.angle();
	if (theta_s + theta_t >= beta
	    || theta_s + theta_t + beta >= M_PI)
	    return DirectionCone(null, 4.0);

	Point centre = s % t;

	double sin_s = sin(theta_s);
	double sin_t = sin(theta_t);
	double angle = sqrt(sin_s * sin_s + 2.0 * sin_s * sin_t * cos_beta
			    + sin_t * sin_t);
	double sin_beta = sin(beta);
	angle = asin(angle / sin_beta);

	DirectionCone normal_cone(centre, angle);
	return normal_cone;
    }
    else if (method == sislBased)
      {
	// We are making a cone surrounding the orientating surface on the
	// unit sphere. The cone is representated with centre coordinates
	// and an angle. The orientation is computed from aproximation of
	// the normal to the surface.  Based on the sisl function s1990.

	double tol = 0.1;
	Point null(0.0, 0.0, 0.0);
	Point axis(dim_);
	double angle = 0.0;
	int in1 = numCoefs_u();
	int in2 = numCoefs_v();

	Point corner[4];  // The coefficients making the corner of
	// each patch
	Point diff[4];    // Difference vector between corner
	// coefficients
	Point norm[4];    // Estimated surface normal (cross product
	// between difference vectors)
	int kver, khor;   // The index to the vertice in the upper
	// left corner to the patch to treat.
	vector<double>::const_iterator it1;
	int ki, kj;

	// Here we are treating each patch in the control polygon
	// separately.
	bool first = true;
	for (it1=coefs_begin(), kver=0; kver < (in2-1); kver++, it1+=dim_)
	  for (khor=0; khor < (in1-1); khor++, it1+=dim_)
	    {
	      // Here we make the tangents in each corner of the
	      // patch, and in direction with the clock. The first
	      // and the last vector contains both the first
	      // tangent.
	      corner[0].resize(dim_);
	      corner[0].setValue(&it1[0]);
	      corner[1].resize(dim_);
	      corner[1].setValue(&it1[dim_]);
	      corner[2].resize(dim_);
	      corner[2].setValue(&it1[(in1+1)*dim_]);
	      corner[3].resize(dim_);
	      corner[3].setValue(&it1[in1*dim_]);
	      for (ki=0; ki<4; ki++)
		{
		  kj = ((ki+1) % 4);
		  diff[ki] = corner[kj] - corner[ki];
		}
	  
	      // Here we makes the normales in each corner of the
	      // patch.  We are using a cross product between two
	      // tangents.  The normals ar also normalized.
	      int count = 0;
	      for (ki=0; ki<4; ki++)
		{
		  kj = (ki == 0) ? 3 : ki-1;
		  norm[ki] = diff[kj].cross(diff[ki]);
		  double len = norm[ki].normalize_checked();
		  if (len == 0.0)
		    count++;
		}

	      if (count == 4)
		continue;  // Degenerate control polygon patch. No contribution to the cone

	      if (first)
		{
		  // Computing the center coordinate of the cone
		  for (kj=0; kj<dim_; ++kj)
		    {
		      double tmin = 1.0;
		      double tmax = -1.0;
		      for (ki=0; ki<4; ++ki)
			{
			  tmin = std::min(tmin, norm[ki][kj]);
			  tmax = std::max(tmax, norm[ki][kj]);
			}
		      axis[kj] = 0.5*(tmin + tmax);
		    }
		  double len = axis.normalize_checked();
		  if (len > 0)
		    {
		      // Computing the angle of the cone
		      for (ki=0; ki<4; ++ki)
			{
			  double ang = axis.angle(norm[ki]);
			  angle = std::max(angle, ang);
			}
		      first = false;
		    }
		}
	      else
		{
		  // Computing the new center and angle of the cone
		  for (ki=0; ki<4; ++ki)
		    {
		      double ang = axis.angle(norm[ki]);
		      if (angle + ang > M_PI+tol)
			{
			  // The angle is too large
			  return DirectionCone(null, 4.0);
			}
		      else if (ang > angle)
			{
			  // The normal is not inside the cone. We have to 
			  // compute a new cone
			  double sin_tang = sin(ang);      
			  double delta    = (ang - angle)/2.0;

			  double t1 = sin(delta)/sin_tang;                   
			  double t2 = sin(ang - delta)/sin_tang;  
			  axis  *= t2;
			  axis += norm[ki]*t1;

			  angle = 0.5*(angle + ang);
			}
		    }
		}
	    }
    
	DirectionCone normal_cone(axis, angle);
	return normal_cone;
      }

    else {
      THROW("This can't happen in SplineSurface::normalCone()!");
      return DirectionCone();
    }

    return DirectionCone();

}


//===========================================================================
DirectionCone SplineSurface::normalCone() const
//===========================================================================
{
  return normalCone(sislBased);
}


//===========================================================================
DirectionCone SplineSurface::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    ALWAYS_ERROR_IF(dim_ != 3, "Normal only defined in 3D");

    int nu = numCoefs_u();
    int nv = numCoefs_v();
    const double* start = &coefs_[0];
    Point null(0.0, 0.0, 0.0);
    Point tmp(dim_);
    int i;

    if (pardir_is_u)
    {
	// Make cone in u-direction
	vector<double> coefs_u;
	for (i = 0; i < nu-1; ++i) {
	    for (int j = 0; j < nv; ++j) {
		int a = dim_ * (nu*j + i); // (i, j)
		int b = dim_ * (nu*j + i + 1); // (i+1, j)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_u.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_u = (int)coefs_u.size();
	DirectionCone cone_u;
	cone_u.setFromArray(&coefs_u[0], &coefs_u[0] + size_u, dim_);
	return cone_u;
    }
    else
    {
	// Then in v-direction
	vector<double> coefs_v;
	for (i = 0; i < nu; ++i) {
	    for (int j = 0; j < nv-1; ++j) {
		int a = dim_ * (nu*j + i); // (i, j)
		int b = dim_ * (nu*(j+1) + i); // (i, j+1)
		for (int dd = 0; dd < dim_; ++dd) {
		    tmp[dd] = *(start + b + dd) - *(start + a + dd);
		}
		if (tmp.length() > 0.0) { // Ignore vectors of zero length
		    for (int dd = 0; dd < dim_; ++dd) {
			coefs_v.push_back(tmp[dd]);
		    }
		}
	    }
	}
	int size_v = (int)coefs_v.size();
	DirectionCone cone_v;
	cone_v.setFromArray(&coefs_v[0], &coefs_v[0] + size_v, dim_);
	return cone_v;
    }
}

//===========================================================================
int SplineSurface::dimension() const
//===========================================================================
{
    return dim_;
}


//===========================================================================
ClassType SplineSurface::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
double SplineSurface::startparam_u() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_u_.begin();
    return knot[basis_u_.order() - 1];
}


//===========================================================================
double SplineSurface::startparam_v() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_v_.begin();
    return knot[basis_v_.order() - 1];
}


//===========================================================================
double SplineSurface::endparam_u() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_u_.begin();
    return knot[basis_u_.numCoefs()];
}
//===========================================================================
double SplineSurface::endparam_v() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_v_.begin();
    return knot[basis_v_.numCoefs()];
}


//===========================================================================
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
// const Domain& SplineSurface::parameterDomain() const
// #else
// const RectDomain& SplineSurface::parameterDomain() const
// #endif // _MSC_VER < 1300
// #else
const RectDomain& SplineSurface::parameterDomain() const
// #endif
//===========================================================================
{
    Vector2D ll(basis_u_.startparam(), basis_v_.startparam());
    Vector2D ur(basis_u_.endparam(), basis_v_.endparam());
    domain_ = RectDomain(ll, ur);
    return domain_;
}


//===========================================================================
RectDomain SplineSurface::containingDomain() const
//===========================================================================
{
  // The rectangular domain containing the domain of the surface is
  // the same as the domain of the splinesurface. This is because the
  // domain of a splinesurface is rectangular.
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//   RectDomain rd = dynamic_cast<const RectDomain&>(parameterDomain());
//   return rd;
// #else
//   return parameterDomain();
// #endif // _MSC_VER < 1300
// #else
  return parameterDomain();
// #endif
}

//===========================================================================
bool SplineSurface::inDomain(double u, double v) const
//===========================================================================
{
    if (u < startparam_u() || u > endparam_u())
	return false;
    if (v < startparam_v() || v > endparam_v())
	return false;

    return true;
}

//===========================================================================
Point SplineSurface::closestInDomain(double u, double v) const
//===========================================================================
{
  double u1 = std::min(std::max(u, startparam_u()), endparam_u());
  double v1 = std::min(std::max(v, startparam_v()), endparam_v());
  return Point(u1, v1);
}

//===========================================================================
CurveLoop SplineSurface::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    // Test for degeneracy.
    bool deg[4];
    if (degenerate_epsilon < 0.0)
      deg[0] = deg[1] = deg[2] = deg[3] = false; // All curves are wanted
    else
      isDegenerate(deg[0], deg[1], deg[2], deg[3], degenerate_epsilon);
    std::vector< shared_ptr< ParamCurve > >  vec;
    for (int edgenum = 0; edgenum < 4; ++edgenum) {
	if (!deg[edgenum]) {
	    shared_ptr<ParamCurve> edgecurve (edgeCurve(edgenum));
	    if (edgenum == 2 || edgenum == 3)
		edgecurve->reverseParameterDirection();
	    vec.push_back(edgecurve);
	}
    }

    return CurveLoop(vec, (degenerate_epsilon < 0.0) ? DEFAULT_SPACE_EPSILON :
		     degenerate_epsilon);
}


//===========================================================================
std::vector<CurveLoop> SplineSurface::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
  // There is only one boundary loop...
  std::vector<CurveLoop> cvloopvec;
  cvloopvec.push_back(outerBoundaryLoop(degenerate_epsilon));
  return cvloopvec;
}

//===========================================================================
SplineSurface* SplineSurface::mirrorSurface(const Point& pos, 
					    const Point& norm) const
//===========================================================================
{
  Point normal = norm;
  normal.normalize();
  vector<double>::const_iterator sc = rational_ ? rcoefs_begin() : coefs_begin();
  int ncoefs = numCoefs_u()*numCoefs_v();
  int kdim = dim_ + (rational_ == true);
  vector<double> sc2(ncoefs*kdim);
  int ki, kj;
  for (ki=0; ki<ncoefs; ++ki, sc+=kdim)
    {
      Point tmp(sc, sc+dim_);
      Point tmp2 = ((tmp - pos)*normal) * normal;
      for (kj=0; kj<dim_; ++kj)
	sc2[ki*kdim+kj] = tmp[kj] - 2.0*tmp2[kj];
      if (rational_)
	sc2[ki*kdim+kdim] = sc[kdim];
    }

  SplineSurface* mirrored = new SplineSurface(basis_u_, basis_v_,
					      &sc2[0], dim_, rational_);

  return mirrored;
}


//===========================================================================
void SplineSurface::interpolate(Interpolator& interpolator1,
				  Interpolator& interpolator2,
				  int num_points1,
				  int num_points2,
				  int dim,
				  const double* param1_start,
				  const double* param2_start,
				  const double* data_start)
//===========================================================================
{
    
    std::vector<double> stage1coefs;

    // Interpolate in the second parameter direction
    interpolator2.interpolate(num_points2, 
			      num_points1*dim,
			      param2_start,
			      data_start,
			      stage1coefs);

    basis_v_ = interpolator2.basis();
    int num_coefs2 = basis_v_.numCoefs();

    // Transpose the stage1coefs array
    SplineUtils::transpose_array(dim, num_coefs2, num_points1, &stage1coefs[0]);

    // Interpolate in the first parameter direction
    interpolator1.interpolate(num_points1,
			      dim*num_coefs2,
			      param1_start,
			      &stage1coefs[0],
			      coefs_);

    basis_u_ = interpolator1.basis();
    int num_coefs1 = basis_u_.numCoefs();

    SplineUtils::transpose_array(dim, num_coefs1, num_coefs2, &coefs_[0]);

    dim_ = dim;
}

//===========================================================================
double SplineSurface::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    double startpar = (dir == 0) ? startparam_u() : startparam_v();
    double endpar = (dir == 0) ? endparam_u() : endparam_v();
    const BsplineBasis& basis = (dir == 0) ? basis_u() : basis_v();

  if (!forward && par <= startpar)
    return startpar;
  
  if (forward && par >= endpar)
    return endpar;

  std::vector<double>::const_iterator knot;
  if (forward) {
    par += fabs(tol);
    knot = std::upper_bound(basis.begin(),basis.end(),par);
    if (knot == basis.end())
      return endpar;
    else
      return *knot;
  }
  else {
    par -= fabs(tol);
    for (knot=basis.end()-1; knot>basis.begin(); --knot) {
      if (*knot < par)
	return *knot;
    }
    return *(basis.begin());
  }
}

//===========================================================================
void SplineSurface::replaceCoefficient(int ix, Point coef)
//===========================================================================
{
  ASSERT(dim_ == coef.dimension());
  vector<double>::iterator c1 = coefs_begin() + ix*dim_;
  for (int ki=0; ki<dim_; ++ki)
    c1[ki] = coef[ki];

  if (rational_)
    {
      vector<double>::iterator c2 = rcoefs_begin() + ix*dim_;
      for (int ki=0; ki<dim_; ++ki)
	c2[ki] = c1[ki]*c2[dim_];
    }
}

//===========================================================================
void SplineSurface::getWeights(std::vector<double>& weights) const
//===========================================================================
{
    int ncoefs = basis_u_.numCoefs()*basis_v_.numCoefs();
    if ((int)weights.size() != ncoefs)
	weights.resize(ncoefs);
    if (rational_)
    {
	int ki, ki2;
	for (ki=0, ki2=0; ki<ncoefs; ++ki, ki2+=(dim_+1))
	    weights[ki] = rcoefs_[ki2+dim_];
    }
    else
	std::fill(weights.begin(), weights.end(), 1.0);
    
}


//===========================================================================
bool SplineSurface::isDegenerate(bool& bottom, bool& right, bool& top, 
				   bool& left, double epsilon) const
//===========================================================================
{
  if (degen_.is_set_ && fabs(degen_.tol_ - epsilon) < 1.0e-15)
    {
      bottom = degen_.b_;
      right = degen_.r_;
      top = degen_.t_;
      left = degen_.l_;
    }
  else if (basis_u_.isKreg() && basis_v_.isKreg())
    {
    int i0 = numCoefs_u();
    int i1 = numCoefs_v();

    // If there is only one point (i.e. order is one in this direction...)
    // the surface is not considered degenerate
    bottom = true;
    int idx[4] = {0, (i0-1)*dim_, (i1-1)*i0*dim_, 0};
    int delta[4] = {dim_, dim_*i0, dim_, dim_*i0};
    bool b[4] = {true, true, true, true};
    for (int side = 0; side < 4; ++side) {
	int imax = side==0 || side==2 ? i0 : i1;
	Point cumuldist(dim_);
	std::fill(cumuldist.begin(), cumuldist.end(), 0.0);
	for (int i = 1; i < imax && b[side]; ++i) {
	    for (int j = 0; j < dim_; ++j) {
		cumuldist[j] += fabs(coefs_[idx[side] + i*delta[side] + j]
			    - coefs_[idx[side] + (i-1)*delta[side] + j]);
		if (cumuldist[j] > epsilon) {
		    b[side] = false;
		    break;
		}
	    }
	}
    }

    bottom = b[0];
    right = b[1];
    top = b[2];
    left = b[3];

    degen_.is_set_ = true;
    degen_.tol_ = epsilon;
    degen_.b_ = bottom;
    degen_.l_ = left;
    degen_.t_ = top;
    degen_.r_ = right;
    }
  else
    {
// #ifdef _MSC_VER
//       SplineSurface *sfkreg = dynamic_cast<SplineSurface*>(subSurface(startparam_u(),
// 					   startparam_v(),
// 					   endparam_u(),
// 					   endparam_v()));
// #else
      SplineSurface *sfkreg = subSurface(startparam_u(),
					   startparam_v(),
					   endparam_u(),
					   endparam_v());
// #endif
      sfkreg->isDegenerate(bottom, right, top, left, epsilon);
      degen_.is_set_ = true;
      degen_.tol_ = epsilon;
      degen_.b_ = bottom;
      degen_.l_ = left;
      degen_.t_ = top;
      degen_.r_ = right;
      delete sfkreg;
    }


    return left || right || top || bottom;
}

//===========================================================================
void SplineSurface::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    // Parameter values in corners
    double param[8];
    param[0] = param[4] = startparam_u();
    param[2] = param[6] = endparam_u();
    param[1] = param[3] = startparam_v();
    param[5] = param[7] = endparam_v();

    // For all corners
    vector<Point> derivs(3);
    double ang;
    for (int ki=0; ki<4; ki++)
    {
	point(derivs, param[2*ki], param[2*ki+1], 1);
	ang = derivs[1].angle(derivs[2]);
	if (fabs(ang) < tol || fabs(M_PI-ang) < tol)
	    deg_corners.push_back(Point(param[2*ki], param[2*ki+1]));
    }
}

//===========================================================================
void SplineSurface::getCornerPoints(vector<pair<Point,Point> >& corners) const
//===========================================================================
{
  corners.resize(4);

  // Parameter values in corners
  double param[8];
  param[0] = param[6] = startparam_u();
  param[2] = param[4] = endparam_u();
  param[1] = param[3] = startparam_v();
  param[5] = param[7] = endparam_v();

  // For all corners
  Point pos;
  for (int ki=0; ki<4; ki++)
    {
      point(pos, param[2*ki], param[2*ki+1]);
      corners[ki] = make_pair(pos, Point(param[2*ki], param[2*ki+1]));
    }
}

 //===========================================================================
void SplineSurface::turnOrientation()
//===========================================================================
{
  swapParameterDirection();
}

//===========================================================================
void SplineSurface::swapParameterDirection()
//===========================================================================
{
    if (rational_) {
	SplineUtils::transpose_array(dim_+1, numCoefs_v(), numCoefs_u(),
			&(activeCoefs()[0]));
	updateCoefsFromRcoefs();
    } else {
	SplineUtils::transpose_array(dim_, numCoefs_v(), numCoefs_u(),
			&(activeCoefs()[0]));
    }
    basis_u_.swap(basis_v_);
    if (degen_.is_set_) {
	std::swap(degen_.b_, degen_.l_);
	std::swap(degen_.t_, degen_.r_);
    }
    // @@ Maybe we should do something about the spatial boundary?
}

//===========================================================================
void SplineSurface::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    if (direction_is_u) {
	// This could be done more rapidly on-the-spot, but for the moment,
	// the current implementation will do....
	swapParameterDirection();
	reverseParameterDirection(false);
	swapParameterDirection();
    } else {
	// we will treat the surface coefficients as curve coefficients 
	// on a spline curve in high dimensionality
	int kdim = dim_ + (rational_ ? 1 : 0) ;
	int n = numCoefs_v();
	int high_dimension = numCoefs_u() * kdim;
	std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
	for (int i = 0; i < n/2; ++i) {
	    for (int dd = 0; dd < high_dimension; ++dd) {
		std::swap(co[i * high_dimension + dd], co[(n - 1 - i) * high_dimension + dd]);
	    }
	}
	basis_v_.reverseParameterDirection();
	if (rational_) {
	    updateCoefsFromRcoefs();
	}
	if (degen_.is_set_) {
	    std::swap(degen_.b_, degen_.t_);
	}
    }
}

//===========================================================================
void SplineSurface::setParameterDomain(double u1, double u2, 
					 double v1, double v2)
//===========================================================================
{
  basis_u_.rescale(u1, u2);
  basis_v_.rescale(v1, v2);
  Vector2D ll(basis_u_.startparam(), basis_v_.startparam());
  Vector2D ur(basis_u_.endparam(), basis_v_.endparam());
  domain_ = RectDomain(ll, ur);
} 

//===========================================================================
void SplineSurface::removeKnot_u(double upar)
//===========================================================================
{
    // We write sf as spline curve, remove knot from cv, transfer back to sf.
    swapParameterDirection();
    removeKnot_v(upar);
    swapParameterDirection();
}

//===========================================================================
void SplineSurface::removeKnot_v(double vpar)
//===========================================================================
{
    // We write sf as spline curve, remove knot from cv, transfer back to sf.
    int kdim = rational_ ? dim_+1 : dim_;
    // Make a hypercurve from this surface
    SplineCurve cv(numCoefs_v(), order_v(), basis_v_.begin(),
		   activeCoefs().begin(), kdim*numCoefs_u(), false);
    // Insert the knot in the curve
    cv.removeKnot(vpar);
    // Insert knot into basis
    basis_v_.removeKnot(vpar);
    // Copy back the coefficients
    if (!rational_) {
	coefs_.resize(numCoefs_u()*numCoefs_v()*dim_);
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    } else {
	rcoefs_.resize(numCoefs_u()*numCoefs_v()*(dim_+1));
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
	updateCoefsFromRcoefs();
    }

}

//===========================================================================
void SplineSurface:: makeSurfaceKRegular()
//===========================================================================
{
  // Check if the surface is k-regular already
  if (basis_u_.endMultiplicity(true) == basis_u_.order() &&
      basis_u_.endMultiplicity(false) == basis_u_.order() &&
      basis_v_.endMultiplicity(true) == basis_v_.order() &&
      basis_v_.endMultiplicity(false) == basis_v_.order())
    return;  // No knot insertion necessary

  // Pick surface inside parameter domain
  shared_ptr<SplineSurface> kreg = 
    shared_ptr<SplineSurface>(subSurface(startparam_u(),
					 startparam_v(),
					 endparam_u(),
					 endparam_v()));
  swap(*(kreg.get()));
}

//===========================================================================
SplineCurve* SplineSurface::edgeCurve(int ccw_edge_number) const
//===========================================================================
{
    const BsplineBasis& bas
	= (ccw_edge_number == 0 || ccw_edge_number == 2) ? basis_u_ : basis_v_;

    const BsplineBasis& opposite_bas
	= (ccw_edge_number == 0 || ccw_edge_number == 2) ? basis_v_ : basis_u_;


    // Check if the basis has k-tupple knots in the current direction
    bool start = (ccw_edge_number == 0 || ccw_edge_number == 3) ? true : false;
    int mult = opposite_bas.endMultiplicity(start);
    if (mult < opposite_bas.order())
      {
	double par = (start) ? opposite_bas.startparam() : opposite_bas.endparam();
	return constParamCurve(par, (ccw_edge_number == 0 || ccw_edge_number == 2));
      }

    // Handle the rational case
    int kdim = dim_ + (rational_ ? 1 : 0);
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    SplineCurve* sc;

    if (ccw_edge_number == 0) {
	sc = new SplineCurve(bas.numCoefs(), bas.order(),
			       bas.begin(), co.begin(),
			       dim_, rational_);
    } else {
	std::vector<double> vec (bas.numCoefs()*kdim);
	for (int i = 0; i < bas.numCoefs(); ++i) {
	    for (int dd = 0; dd < kdim; ++dd) {
		double c;
		switch (ccw_edge_number) {
		case 1:
		    {
			c = co[(numCoefs_u()*(i+1)-1)*kdim + dd];
			break;
		    }
		case 2:
		    {
			c = co[(numCoefs_u()*(numCoefs_v()-1) + i)*kdim + dd];
			break;
		    }
		case 3:
		    {
			c = co[(numCoefs_u()*i)*kdim + dd];
			break;
		    }
		default:
		    {
			THROW("Wrong edge number: " << ccw_edge_number);
			break;
		    }
		}
		vec[i*kdim + dd] = c;
	    }
	}
	sc = new SplineCurve(bas.numCoefs(), bas.order(),
			       bas.begin(), vec.begin(),
			       dim_, rational_);
    }
    return sc;
}



//===========================================================================
double SplineSurface::knotSpan(int pardir, int iknot) const
//===========================================================================
{
  vector<double>::const_iterator knots = basis(pardir).begin();
  int nmbc = (pardir == 0) ? numCoefs_u() : numCoefs_v();
  int order = (pardir == 0) ? order_u() : order_v();
  if (iknot < 0 || iknot >= nmbc + order - 1)
    return 0.0;
  return knots[iknot+1] - knots[iknot];
}

//===========================================================================
SplineCurve*
SplineSurface::constParamCurve (double parameter,
				bool pardir_is_u) const
//===========================================================================
{
    // Handle the rational case
    int kdim = dim_ + (rational_ ? 1 : 0);
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    std::vector<double> huge_curve_coefs;
    if (!pardir_is_u) {
	// We must flip the surface
	int nu = numCoefs_u();
	int nv = numCoefs_v();
	huge_curve_coefs.reserve(nu*nv*kdim);
	for (int i = 0; i < nu; ++i)
	    for (int j = 0; j < nv; ++j)
		for (int k = 0; k < kdim; ++k)
		    huge_curve_coefs.push_back(co[(j*nu+i)*kdim + k]);
    } else {
	huge_curve_coefs = co;
    }

    const BsplineBasis& bas = pardir_is_u ? basis_u_ : basis_v_;
    const BsplineBasis& other_bas = pardir_is_u ? basis_v_ : basis_u_;

    int num = other_bas.numCoefs();
    int order = other_bas.order();
    std::vector<double>::const_iterator knotstart = other_bas.begin();
    std::vector<double>::const_iterator coefstart = huge_curve_coefs.begin();
    int dim = bas.numCoefs() * kdim;
    SplineCurve huge_curve(num, order, knotstart, coefstart, dim, false);
    std::vector<double> coefs_wanted(dim);
    Point p(&coefs_wanted[0], &coefs_wanted[0] + dim, false);
    huge_curve.point(p, parameter);
    num = bas.numCoefs();
    order = bas.order();
    knotstart = bas.begin();
    SplineCurve* sc = new SplineCurve (num, order,
				       knotstart, coefs_wanted.begin(),
				       dim_, rational_);
    return sc;
}



//===========================================================================
void 
SplineSurface::constParamCurve (double parameter,
				  bool pardir_is_u,
				  SplineCurve*& cv,
				  SplineCurve*& crosscv) const
//===========================================================================
{
    // Handle the rational case
    int kdim = dim_ + (rational_ ? 1 : 0);
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    std::vector<double> huge_curve_coefs;
    if (!pardir_is_u) {
	// We must flip the surface
	int nu = numCoefs_u();
	int nv = numCoefs_v();
	huge_curve_coefs.reserve(nu*nv*kdim);
	for (int i = 0; i < nu; ++i)
	    for (int j = 0; j < nv; ++j)
		for (int k = 0; k < kdim; ++k)
		    huge_curve_coefs.push_back(co[(j*nu+i)*kdim + k]);
    } else {
	huge_curve_coefs = co;
    }

    const BsplineBasis& bas = pardir_is_u ? basis_u_ : basis_v_;
    const BsplineBasis& other_bas = pardir_is_u ? basis_v_ : basis_u_;

    int num = other_bas.numCoefs();
    int order = other_bas.order();
    std::vector<double>::const_iterator knotstart = other_bas.begin();
    std::vector<double>::const_iterator coefstart = huge_curve_coefs.begin();
    int dim = bas.numCoefs() * kdim;
    SplineCurve huge_curve (num, order,
			      knotstart, coefstart,
			      dim, false);

    std::vector<Point> p;
    Point p1(dim), p2(dim);
    p.push_back(p1);
    p.push_back(p2);
    huge_curve.point(p, parameter, 1);
    num = bas.numCoefs();
    order = bas.order();
    knotstart = bas.begin();
    cv = new SplineCurve (num, order, knotstart, p[0].begin(),
			    dim_, rational_);
    if (rational_)
      {
	//@@@VSK. Not implemented yet
	MESSAGE("Cross tangent of rational curve is approximated");

	// Evaluate the cross tangents in the Greville parameters 
	// of the surface
	int kn = bas.numCoefs();
	vector<double> pts(kn*dim_);
	vector<double> params(kn);
	int ki, kj;
	vector<Point> der(3);
	double sf_par[2];
	int idx1 = (pardir_is_u) ? 0 : 1;
	int idx2 = 1 - idx1;
	sf_par[idx1] = parameter;
	for (ki=0; ki<kn; ki++)
	  {
	    sf_par[idx2] = bas.grevilleParameter(ki);
	    params[ki] = sf_par[idx2];
	    point(der, sf_par[0], sf_par[1], 1);
	    for (kj=0; kj<dim_; kj++)
	      pts[ki*dim_+kj] = der[1+idx1][kj];
	  }

	// Interpolate to get cross parameter curve
	vector<double> coefs;
	vector<int> tang_idx;
	vector<double> tang_pts;
	BsplineBasis basis = bas;
	SplineInterpolator interpolate;
	interpolate.setBasis(basis);
	interpolate.interpolate(params, pts, tang_idx, tang_pts, coefs);

	// Create cross tangent curve
	crosscv = new SplineCurve (num, order, knotstart, coefs.begin(),
				     dim_, false);
      }
    else
      crosscv = new SplineCurve (num, order, knotstart, p[1].begin(),
			    dim_, rational_);

}


//===========================================================================
vector<shared_ptr<ParamCurve> >
SplineSurface::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > return_cvs;
    return_cvs.push_back(shared_ptr<ParamCurve>(constParamCurve(parameter,
								pardir_is_u)));

    return return_cvs;
}


//===========================================================================
void SplineSurface::getConstParamCurves(const std::vector<double>& params_u,
					const std::vector<double>& params_v,
					std::vector<shared_ptr<SplineCurve> >& curves_u,
					std::vector<shared_ptr<SplineCurve> >& curves_v)
//===========================================================================
{
  curves_u.resize(0);
  curves_v.resize(0);

  for (size_t i = 0; i < params_u.size(); ++i)
    curves_u.push_back(shared_ptr<SplineCurve> (constParamCurve(params_u[i], false)));
  for (size_t i = 0; i < params_v.size(); ++i)
    curves_v.push_back(shared_ptr<SplineCurve> (constParamCurve(params_v[i], true)));
}



//===========================================================================
void
SplineSurface::getBoundaryIdx(Point& pt1, Point& pt2, 
			      double epsilon, int& bdindex,
			      double& par1, double& par2, double knot_tol) const
//--------------------------------------------------------------------------
  // Given two points on the surface boundary, find the number of the
  // corresponding boundary and the curve parameter of the closest points
  // on this surface boundary.
  //
  // Ordering of boundaries:
  //                       1
  //           ----------------------
  //           |                    |
  //         2 |                    | 3
  //      v    |                    |
  //      ^    ----------------------
  //      |-> u            0
//===========================================================================
{
  // Find parameter value between which the boundary curve passes.
  double u1, v1, u2, v2;
  Point cl1, cl2;
  double d1, d2;
  double tol = 1.0e-7;  // Tolerance in the parameter domain.
  closestBoundaryPoint(pt1, u1, v1, cl1, d1, tol);
  closestBoundaryPoint(pt2, u2, v2, cl2, d2, tol);

  bdindex = -1;
  if (d1 > epsilon || d2 > epsilon)
    return;        // Point not on surface

  // As we are seeking a boundary curve, we know that for both points at
  // least one parameter must be an endpoint in a u- or v-knot vector.
  // We use a value of 1e-05, as the basis' knotIntervalFuzzy functions
  // default 1e-12 tolerance is too strict.
  // @@ The value may be given as a parameter?
  basis_u().knotIntervalFuzzy(u1, knot_tol);
  basis_v().knotIntervalFuzzy(v1, knot_tol);
  basis_u().knotIntervalFuzzy(u2, knot_tol);
  basis_v().knotIntervalFuzzy(v2, knot_tol);

  double startu = startparam_u();
  double endu = endparam_u();
  double startv = startparam_v();
  double endv = endparam_v();
  if (fabs(u1-u2) < fabs(v1-v2) && fabs(u1-u2) < 0.01*(endu-startu))
    {
      double umid = 0.5*(u1+u2);
      par1 = v1;
      par2 = v2;
      bdindex = (fabs(umid - startu) < fabs(endu - umid))
	? 2 : 3;
    }
  else if (fabs(v1-v2) < fabs(u1-u2) && fabs(v1-v2) < 0.01*(endv-startv))
    {
      double vmid = 0.5*(v1+v2);
      par1 = u1;
      par2 = u2;
      bdindex = (fabs(vmid - startv) < fabs(endv - vmid))
	? 0 : 1;
    }
  else
    {
      // No clear boundary is found. Is it possible to save the
      // situation? Check degeneracy. May assume there is only one degenerate edge.
      // In such a scenario we switch a parameter of the closest point.
	if (degen_.is_set_) {
	    if (degen_.b_) {
		if (v1 < v2)
		    u1 = u2;
		else
		    u2 = u1;
	    }
	    else if (degen_.t_) {
		if (v1 > v2)
		    u1 = u2;
		else
		    u2 = u1;
	    }
	    else if (degen_.l_) {
		if (u1 < u2)
		    v1 = v2;
		else
		    v2 = v1;
	    }
	    else if (degen_.r_) {
		if (u1 > u2)
		    v1 = v2;
		else
		    v2 = v1;
	    }

	    // We check whether new points are close enough.
	    Point new_pt1 = ParamSurface::point(u1, v1);
	    Point new_pt2 = ParamSurface::point(u2, v2);
	    double new_u1, new_u2, new_v1, new_v2;
	    closestBoundaryPoint(new_pt1, new_u1, new_v1, cl1, d1, tol);
	    closestBoundaryPoint(new_pt2, new_u2, new_v2, cl2, d2, tol);
	    if (d1 > epsilon || d2 > epsilon)
		return;        // Point not on surface
	    
	    if (fabs(u1-u2) < fabs(v1-v2)) {
		double umid = 0.5*(u1+u2);
		par1 = v1;
		par2 = v2;
		bdindex = (fabs(umid - startu) < fabs(endu - umid))
		    ? 2 : 3;
	    }
	    else if (fabs(v1-v2) < fabs(u1-u2)) {
		double vmid = 0.5*(v1+v2);
		par1 = u1;
		par2 = u2;
		bdindex = (fabs(vmid - startv) < fabs(endv - vmid))
		    ? 0 : 1;
	    }
	}
	else {
	    return;
	}
	
    }
      
  return;
}

//===========================================================================
void
SplineSurface::getBoundaryIdx(Point& pt1, double epsilon, int& bdindex,
			      double knot_tol) const
//--------------------------------------------------------------------------
  // Given one point on the surface boundary, find the number of the
  // corresponding boundary.
  //
  // Ordering of boundaries:
  //                       1
  //           ----------------------
  //           |                    |
  //         2 |                    | 3
  //      v    |                    |
  //      ^    ----------------------
  //      |-> u            0
//===========================================================================
{
  // Find parameter value between which the boundary curve passes.
  double u1, v1;
  Point cl1;
  double d1;
  double tol = 1.0e-7;  // Tolerance in the parameter domain.
  closestBoundaryPoint(pt1, u1, v1, cl1, d1, tol);

  bdindex = -1;
  if (d1 > epsilon)
    return;        // Point not on surface

  basis_u().knotIntervalFuzzy(u1, knot_tol);
  basis_v().knotIntervalFuzzy(v1, knot_tol);

  double startu = startparam_u();
  double endu = endparam_u();
  double startv = startparam_v();
  double endv = endparam_v();

  if (std::min(u1-startu, endu-u1) < std::min(v1-startv, endv-v1))
    bdindex = (u1-startu < endu-u1) ? 2 : 3;
  else
    bdindex = (v1-startv < endv-v1) ? 0 : 1;
}

//===========================================================================
void SplineSurface::appendSurface(ParamSurface* sf, int join_dir, bool repar)
//===========================================================================
{
    int cont = 1;
    double dist_dummy = 0;
    appendSurface(sf, join_dir, cont, dist_dummy, repar);
}

//===========================================================================
double SplineSurface::appendSurface(ParamSurface* sf, int join_dir,
				  int cont, double& dist, bool repar)
//===========================================================================
{
    ASSERT(sf->instanceType() == Class_SplineSurface);
 
   // Describe the sfs as curves in the given parameter
    // direction
    vector<shared_ptr<SplineSurface> > sfs;
    sfs.push_back(shared_ptr<SplineSurface>(clone()));
    sfs.push_back(shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(sf->clone())));
    // Make sure that the surfaces are described in the same spline space
    // in the join direction
    // Should be included, but must be fixed first to avoid new knots
    // due to numeric noice
    // @@@ VSK. I am not sure if unifySurfaceSplineSpace always does the
    // right thing. Check it up!
    double eps = 1.0e-6;  // A rather arbitrary choice 
    GeometryTools::unifySurfaceSplineSpace(sfs, eps, 2-join_dir+1);
    bool make_rational = false;
    if (sfs[0]->rational())// && (!sfs[1]->rational()))
      make_rational = true;
    if (/*(!sfs[0]->rational()) && */sfs[1]->rational())
      make_rational = true;
    //sfs.push_back(shared_ptr<SplineSurface>(sf->clone()));
    vector<shared_ptr<SplineCurve> > curves;
    double max_wgt_diff = 0.0;
    for (size_t ki = 0; ki < sfs.size(); ++ki) {
      if (make_rational)
	sfs[ki]->representAsRational();

      if (sfs[ki]->rational())
	{
	  double wgt_diff =
	    sfs[ki]->setAvBdWeight(1.0, join_dir-1, (ki == 1));
	  max_wgt_diff = std::max(max_wgt_diff, wgt_diff);
	}

	shared_ptr<SplineCurve> cv;
	cv = GeometryTools::representSurfaceAsCurve(*sfs[ki], join_dir);
	curves.push_back(cv);
    }
    // curves[0]->appendCurve(curves[1].get(), (make_rational) ? 0 : cont, 
    // 			   dist, repar);
    curves[0]->appendCurve(curves[1].get(), (rational()) ? 0 : cont, 
    			   dist, repar);

    // Represent the curve as a surface
    const BsplineBasis& common_bas = sfs[0]->basis(2-join_dir); 
    shared_ptr<SplineSurface> joined_sf;
    joined_sf = GeometryTools::representCurveAsSurface(*curves[0], join_dir, common_bas, 
					rational() || make_rational);

    *this = *joined_sf;
    return 0.5*max_wgt_diff;
}

//===========================================================================
void
SplineSurface::getBoundaryInfo(Point& pt1, Point& pt2,
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
  double par1, par2;
  int bdidx;
  getBoundaryIdx(pt1, pt2, epsilon, bdidx, par1, par2, knot_tol);
  if (bdidx < 0)
    return;

  getBoundaryInfo(par1, par2, bdidx, cv, crosscv, knot_tol);
  return;
}

//===========================================================================
void
SplineSurface::getBoundaryInfo(double par1, double par2,
				 int bdindex, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
  // Set no output
  cv = crosscv = 0;
  // Following intuition, the curves go from pt1 to pt2.
  bool turn_curves = false;

  if (par1 > par2)
    turn_curves = true;
  double bdpar;
  switch (bdindex)
    {
    case 0:
      bdpar = startparam_v();
      break;
    case 1:
      bdpar = endparam_v();
      break;
    case 2:
      bdpar = startparam_u();
      break;
    case 3:
      bdpar = endparam_u();
      break;
    default:
      return;
    }

  // Ftech constant parameter curve in par. dir.
  SplineCurve *c1=0, *c2=0;
  constParamCurve(bdpar, bdindex<=1, c1, c2);
// #ifdef _MSC_VER
//   cv = dynamic_cast<SplineCurve*>(c1->subCurve(std::min(par1,par2),
// 						     std::max(par1,par2), knot_tol));
//   crosscv = dynamic_cast<SplineCurve*>(c2->subCurve(std::min(par1,par2),
// 							  std::max(par1,par2), knot_tol));
// #else
  cv = c1->subCurve(std::min(par1,par2), std::max(par1,par2), knot_tol);
  crosscv = c2->subCurve(std::min(par1,par2), std::max(par1,par2), knot_tol);
// #endif
  if (bdindex == 0 || bdindex == 2)
    // We must turn cross-curve so that it points outwards.
    for (std::vector<double>::iterator iter = crosscv->coefs_begin();
	 iter != crosscv->coefs_end(); ++iter)
      iter[0] *= -1.0;
  delete c1;
  delete c2;

  if (turn_curves) {
      cv->reverseParameterDirection();
      crosscv->reverseParameterDirection();
  }
  // else not a boundary curve. Return no curves.



}

//===========================================================================
int
SplineSurface::boundaryIndex(Point& param_pt1, Point& param_pt2) const
//===========================================================================
{
  if (param_pt1.dimension() != 2 || param_pt2.dimension() != 2)
    return -1;   // No identification of boundary is possible

  double umin = startparam_u();
  double umax = endparam_u();
  double vmin = startparam_v();
  double vmax = endparam_v();
  double tol1 = 0.1*(umax-umin);
  double tol2 = 0.1*(vmax-vmin);
  if (fabs(param_pt1[0]-param_pt2[0]) < fabs(param_pt1[1]-param_pt2[1]))
    {
      // A boundary in second parameter direction is expected.
      if (fabs(param_pt1[0]-umin) < tol1 && fabs(param_pt2[0]-umin) < tol1)
	return 3;
      else if (fabs(param_pt1[0]-umax) < tol1 && fabs(param_pt2[0]-umax) < tol1)
	return 1;
    }
  else
    {
      // A boundary in first parameter direction is expected.
      if (fabs(param_pt1[1]-vmin) < tol2 && fabs(param_pt2[1]-vmin) < tol2)
	return 0;
      else if (fabs(param_pt1[1]-vmax) < tol2 && fabs(param_pt2[1]-vmax) < tol2)
	return 2;
    }
  return -1;
}

//===========================================================================
void SplineSurface::swap(SplineSurface& other)
//===========================================================================
{
    std::swap(dim_, other.dim_);
    std::swap(rational_, other.rational_);
    basis_u_.swap(other.basis_u_);
    basis_v_.swap(other.basis_v_);
    coefs_.swap(other.coefs_);
    rcoefs_.swap(other.rcoefs_);
    std::swap(domain_, other.domain_);
    spatial_boundary_.swap(other.spatial_boundary_);
    std::swap(degen_, other.degen_);
}

//===========================================================================
bool SplineSurface::replaceBoundaryCurve(int bd_nmb, 
					 shared_ptr<SplineCurve> bd_crv,
					 bool unify)
//===========================================================================
{
  if ((rational_ && !bd_crv->rational()) ||
      (!rational_ && bd_crv->rational()))
    return false;
  if (dim_ != bd_crv->dimension())
    return false;

  // Make sure that the curve has the same parameter domain as the surface
  shared_ptr<SplineCurve> crv = shared_ptr<SplineCurve>(bd_crv->clone());

  BsplineBasis* basis = (bd_nmb == 0 || bd_nmb == 1) ? &basis_v_ : &basis_u_;

  basis->check();

  crv->setParameterInterval(basis->startparam(), basis->endparam());

  // Make sure that the curve contains all surface knot values
  if (unify)
    {
      vector<double> diff1;
      diff1 = crv->basis().missingKnots(*basis);
      //       std::set_difference(basis->begin(), basis->end(),
      // 			  crv->knotsBegin(), crv->knotsEnd(),
      // 			  std::back_inserter(diff1));
      // Make sure that the surface contain all curve knots
      vector<double> diff2;
      diff2 = basis->missingKnots(crv->basis());
      //       std::set_difference(crv->knotsBegin(), crv->knotsEnd(),
      // 			  basis->begin(), basis->end(),
      // 			  std::back_inserter(diff2));

      // To avoid extremely dense knots
      size_t ki;
      int left = 0;
      for (ki=0; ki<diff1.size(); ++ki)
	left = crv->basis().knotIntervalFuzzy(diff1[ki]);
	
      if (diff1.size() > 0)
	crv->insertKnot(diff1);

      for (ki=0; ki<diff2.size(); ++ki)
	left = basis->knotIntervalFuzzy(diff2[ki]);

      if (bd_nmb == 0 || bd_nmb == 1)
	{
	  if (diff2.size() > 0)
	    insertKnot_v(diff2);
	  basis_v_.check();
	}
      else
	{
	  if (diff2.size() > 0)
	    insertKnot_u(diff2);
	  basis_u_.check();
	}
    }

  // Replace coefficients
  int nmb_coef = crv->numCoefs();
  int dim = crv->dimension() + rational_;
  vector<double>::iterator cv_coef = 
    (rational_) ? crv->rcoefs_begin() : crv->coefs_begin();
  vector<double>::iterator sf_coef =  (rational_) ? rcoefs_begin() : coefs_begin();

  int stop = nmb_coef*dim;
  int ki, kj;
  if (bd_nmb == 3)
    sf_coef += basis_u_.numCoefs()*(basis_v_.numCoefs()-1)*dim;
  if (bd_nmb == 2 || bd_nmb == 3)
    {
      for (ki=0; ki<stop; ++ki, ++cv_coef, ++sf_coef)
	*sf_coef = *cv_coef;
    }
  if (bd_nmb == 1)
    sf_coef += (basis_u_.numCoefs()-1)*dim;
  if (bd_nmb == 0 || bd_nmb == 1)
    {
      for (ki=0; ki<nmb_coef; ++ki)
	{
	  for (kj=0; kj<dim; ++kj, ++cv_coef, ++sf_coef)
	    *sf_coef = *cv_coef;
	  sf_coef += (basis_u_.numCoefs()-1)*dim;
	}
    }

  if (rational_)
    updateCoefsFromRcoefs();

  return true;
}


// Added by KMO for ICADA usage.
//===========================================================================
void SplineSurface::deform(const std::vector<double>& vec, int vdim)
//===========================================================================
{
  int i, j;
  vector<double>::iterator it;
  if (vdim == 0) vdim = dim_;
  ALWAYS_ERROR_IF(vec.size()*dim_ < coefs_.size()*vdim, "Deformation vector is too short");

  // Change rcoefs_ if rational
  if (rational_)
    for (it = rcoefs_.begin(), j = 0; it != rcoefs_.end(); ++it)
      {
	double w = it[dim_];
	for (i = 0; i < dim_ && i < vdim; ++i)
	  it[i] += vec[j+i] * w;
	it += dim_;
	j += vdim;
      }

  // Change coefs, both when rational and not rational
  for (it = coefs_.begin(), j = 0; it != coefs_.end(); )
    {
      for (i = 0; i < dim_ && i < vdim; ++i)
	it[i] += vec[j+i];
      it += dim_;
      j += vdim;
    }
}


//===========================================================================
void SplineSurface::add(const SplineSurface* other, double tol)
//===========================================================================
{
  int ord_u = basis_u_.order();
  int ord_v = basis_v_.order();
  int ncoefs_u = basis_u_.numCoefs();
  int ncoefs_v = basis_v_.numCoefs();
  int ncoefs = ncoefs_u * ncoefs_v;
  int kdim = rational_ ? dim_+1 : dim_;
  ALWAYS_ERROR_IF(ord_u != other->basis_u_.order(), "Surfaces have different order in first parameter direction");
  ALWAYS_ERROR_IF(ord_v != other->basis_v_.order(), "Surfaces have different order in second parameter direction");
  ALWAYS_ERROR_IF(ncoefs_u != other->basis_u_.numCoefs(), "Surfaces have different number of coefficients in first parameter direction");
  ALWAYS_ERROR_IF(ncoefs_v != other->basis_v_.numCoefs(), "Surfaces have different number of coefficients in second parameter direction");
  ALWAYS_ERROR_IF(dim_ != other->dim_, "Surfaces has different geometry space dimension");
  ALWAYS_ERROR_IF(rational_ != other->rational_, "Can not add rational to non-rational surface");

  vector<double>::const_iterator knots, knots_other;
  knots = basis_u_.begin();
  knots_other = other->basis_u_.begin();
  for (int i = 0; i < ncoefs_u + ord_u; ++i)
    ALWAYS_ERROR_IF(knots[i] != knots_other[i], "Surfaces have different knot vector in first parameter direction");
  knots = basis_v_.begin();
  knots_other = other->basis_v_.begin();
  for (int i = 0; i < ncoefs_v + ord_v; ++i, ++knots, ++knots_other)
    ALWAYS_ERROR_IF((*knots) != (*knots_other), "Surfaces have different knot vector in second parameter direction");

  if (rational_)
    {
      vector<double>::const_iterator weights = rcoefs_.begin() + dim_;
      vector<double>::const_iterator weights_other = other->rcoefs_.begin() + dim_;
      for (int i = 0; i < ncoefs; ++i, weights += kdim, weights_other += kdim)
	ALWAYS_ERROR_IF(fabs ((*weights) - (*weights_other)) >= tol, "Surfaces have different weights");
    }

  // All tests passed - we can now add coefficients

  vector<double>::iterator coefs_it = coefs_.begin();
  vector<double>::const_iterator other_coefs_it = other->coefs_.begin();
  for (int i = 0; i < ncoefs * kdim; ++i, ++coefs_it, ++other_coefs_it)
    (*coefs_it) += (*other_coefs_it);

  // Update rcoefs_ in rational case

  if (rational_)
    {
      coefs_it = coefs_.begin();
      vector<double>::iterator rcoefs_it = rcoefs_.begin();
      for (int i = 0; i < ncoefs; ++i, coefs_it += dim_, rcoefs_it += kdim)
	for (int k = 0; k < dim_; ++k)
	  rcoefs_it[k] = rcoefs_it[dim_] * coefs_it[k];
    }

}


//===========================================================================
void SplineSurface::representAsRational()
//===========================================================================
{
  if (rational_)
    return;   // This surface is already rational

  int in = basis_u_.numCoefs() * basis_v_.numCoefs();
  rcoefs_.resize(in*(dim_+1));
  int ki, kj;
  for (ki=0; ki<in; ++ki)
    {
      for (kj=0; kj<dim_; ++kj)
	rcoefs_[ki*(dim_+1)+kj] = coefs_[ki*dim_+kj];
      rcoefs_[ki*(dim_+1)+dim_] = 1.0;
    }
  rational_ = true;
}


//===========================================================================
double SplineSurface::setAvBdWeight(double wgt, int pardir, bool at_start)
//===========================================================================
{
  if (!rational_)
    return 0.0;   // This surface is not rational

  int kdim = dimension() + 1;
  vector<double>::iterator c1 = rcoefs_begin();
  int offset = 0;
  int kn1 = numCoefs_u();
  int kn2 = numCoefs_v();
  double avwgt = 0.0;
  double maxwgt = 0.0;
  double minwgt = HUGE;
  int ki;
  if (pardir == 0)
    {
      // Compute average weight
      int ndel = kdim*kn1;
       if (!at_start)
	{
	  offset += kdim*(kn1-1);
	}
       for (ki=0; ki<kn2; offset+=ndel, ++ki) {
	 double weight = c1[offset + dim_];
	 avwgt += weight;
	 maxwgt = std::max(maxwgt, weight);
	 minwgt = std::min(minwgt, weight);
       }
      avwgt /= (double)kn2;
    }
  else
    {
      // Compute average weight
      int ndel = kdim;
       if (!at_start)
	{
	  c1 += kdim*kn1*(kn2-1);
	}
       for (ki=0; ki<kn1; offset+=ndel, ++ki) {
	 double weight = c1[offset + dim_];
	 avwgt += weight;
	 maxwgt = std::max(maxwgt, weight);
	 minwgt = std::min(minwgt, weight);
       }
      avwgt /= (double)kn1;
     }

  // Set new weights
  double fac = wgt/avwgt;  // Weights are supposed to be positive and not zero
  for (size_t ki=0; ki<rcoefs_.size(); ++ki)
    rcoefs_[ki] *= fac;

  return maxwgt - minwgt;
}

//===========================================================================
bool SplineSurface::isAxisRotational(Point& centre, Point& axis, Point& vec,
				     double& angle)
//===========================================================================
{
  if (elementary_surface_.get())
    return elementary_surface_->isAxisRotational(centre, axis, vec, angle);
  else
    return false;
}

//===========================================================================
bool SplineSurface::isPlanar(Point& normal, double tol)
//===========================================================================
{
  if (elementary_surface_.get())
    return elementary_surface_->isPlanar(normal, tol);
  else
    {
      return ParamSurface::isPlanar(normal, tol);
    }
}

//===========================================================================
shared_ptr<ElementarySurface> SplineSurface::getElementarySurface()
//===========================================================================
{
    if (is_elementary_surface_)
        return elementary_surface_;
    return shared_ptr<ElementarySurface>();
}


//===========================================================================
void SplineSurface::setElementarySurface(shared_ptr<ElementarySurface> elsurf)
//===========================================================================
{
    elementary_surface_ = elsurf;
    is_elementary_surface_ = true;
}


//===========================================================================
bool SplineSurface::checkElementarySurface()
//===========================================================================
{
    MESSAGE("Not yet implemented");
    return true;
}


//===========================================================================
//
//                 Private helper functions
//
//===========================================================================



//===========================================================================
void SplineSurface::updateCoefsFromRcoefs()
//===========================================================================
{
    coefs_.resize(numCoefs_u()*numCoefs_v()*dim_);
    SplineUtils::make_coef_array_from_rational_coefs(&rcoefs_[0],
					&coefs_[0],
					numCoefs_u()*numCoefs_v(),
					dim_);
}




} // namespace Go
