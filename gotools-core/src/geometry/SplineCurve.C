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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/ElementaryCurve.h"

#include <iomanip>

using namespace std;

namespace Go
{

//===========================================================================
SplineCurve::SplineCurve(const Point& pnt1, const Point& pnt2)
  // Make a linear spline curve interpolating given end points.
//===========================================================================
  : dim_(pnt1.dimension()), rational_(false), is_elementary_curve_(false)
{
    ALWAYS_ERROR_IF(pnt1.dimension() != pnt2.dimension(),"Parameter mismatch.");

  // Chord length parameterization (but not a very small knot
  // interval).
  double knotdist = pnt1.dist(pnt2);
  double knots[4];
  knots[0] = knots[1] = 0.0;
  knots[2] = knots[3] = std::max(knotdist, 0.01);

  // Linear spline space
  BsplineBasis b2(2, 2, knots);
  basis_ = b2;

  // Make coefficients
  int ki;
  for (ki=0; ki<dim_; ki++)
    coefs_.push_back(pnt1[ki]);
  for (ki=0; ki<dim_; ki++)
    coefs_.push_back(pnt2[ki]);
}


//===========================================================================
SplineCurve::SplineCurve(const Point& pnt1, double startpar,
			 const Point& pnt2, double endpar)
  // Make a linear spline curve interpolating given end points.
//===========================================================================
  : dim_(pnt1.dimension()), rational_(false), is_elementary_curve_(false)
{
    ALWAYS_ERROR_IF(pnt1.dimension() != pnt2.dimension(),"Parameter mismatch.");

  // Set knot vector
  //  double knotdist = pnt1.dist(pnt2);
  double knots[4];
  knots[0] = knots[1] = startpar;
  knots[2] = knots[3] = endpar;

  // Linear spline space
  BsplineBasis b2(2, 2, knots);
  basis_ = b2;

  // Make coefficients
  int ki;
  for (ki=0; ki<dim_; ki++)
    coefs_.push_back(pnt1[ki]);
  for (ki=0; ki<dim_; ki++)
    coefs_.push_back(pnt2[ki]);
}


//===========================================================================
SplineCurve::~SplineCurve()
//===========================================================================
{

}

//===========================================================================
void SplineCurve::read (std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }
    is >> dim_;
    is >> rational_;
    is >> basis_;
    int nc = basis_.numCoefs();
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
void SplineCurve::write (std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);

    os << dim_ << ' ' << rational_ << '\n';
    os << basis_;
    int n = basis_.numCoefs();
    int kdim = dim_ + (rational_ ? 1 : 0);
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    for (int i = 0; i < n; ++i) {
	os << co[i*kdim];
	for (int d = 1; d < kdim; ++d) {
	    os << ' ' << co[i*kdim + d];
	}
	os << '\n';
    }
    os << std::endl;
    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
BoundingBox SplineCurve::boundingBox() const
//===========================================================================
{
    BoundingBox box;
    box.setFromArray(&coefs_[0], &coefs_[0] + coefs_.size(), dim_);

    return box;
}


//===========================================================================
CompositeBox SplineCurve::compositeBox() const
//===========================================================================
{
    CompositeBox box(&coefs_[0], dim_, numCoefs(), 1);
    return box;
}


//===========================================================================
DirectionCone SplineCurve::directionCone() const
//===========================================================================
{
    shared_ptr<SplineCurve> dc(derivCurve(1));

    DirectionCone cone;
    cone.setFromArray(&(dc->coefs_[0]), 
		      &(dc->coefs_[0]) + (dc->coefs_).size(), dim_);
    return cone;
}


//===========================================================================
int SplineCurve::dimension() const
//===========================================================================
{
    return dim_;
}


//===========================================================================
ClassType SplineCurve::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
double SplineCurve::startparam() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_.begin();
    return knot[basis_.order() - 1];
}


//===========================================================================
double SplineCurve::endparam() const
//===========================================================================
{
    std::vector<double>::const_iterator knot = basis_.begin();
    return knot[basis_.numCoefs()];
}


//===========================================================================
SplineCurve* SplineCurve::geometryCurve()
//===========================================================================
{
//   return this;
    // We return a copy of this object, to avoid differing memory
    // handling depending on curve type.
    return clone();
}


//===========================================================================
void SplineCurve::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    int kdim = dim_ + (rational_ ? 1 : 0);
    int n = numCoefs();
    int i;
    Point tmp(kdim);
    std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    for (i = 0; i < n/2; ++i) {
	for (int dd = 0; dd < kdim; ++dd) {
	    tmp[dd] = co[i*kdim + dd];
	    co[i*kdim + dd] = co[(n-1-i)*kdim + dd];
	    co[(n-1-i)*kdim + dd] = tmp[dd];
	}
    }
    basis_.reverseParameterDirection();
    if (rational_) {
	updateCoefsFromRcoefs();
    }

    double xtmp;
    if (dim_ == 2 && switchparam)
      {
	// Switch the x- and y-coefficients of the curve.
	for (i=0; i<n; i++)
	  {
	    xtmp = co[i*kdim];
	    co[i*kdim] = co[i*kdim+1]; 
	    co[i*kdim+1] = xtmp;
	  }
	if (rational_) 
	  updateCoefsFromRcoefs();
      }
}


//===========================================================================
bool SplineCurve::isDegenerate(double degenerate_epsilon)
//===========================================================================
{
    bool deg = true;
    double cumuldist = 0.0;
    Point last(&coefs_[0], &coefs_[dim_], false);
    Point current(dim_);
    for (int i = 1; i < numCoefs(); ++i) {
	for (int dd = 0; dd < dim_; ++dd) {
	    current[dd] = coefs_[i*dim_ + dd];
	}
	cumuldist += last.dist(current);
	last = current;	   
	if (cumuldist >= degenerate_epsilon) {
	    deg = false;
	    break;
	}
    }
    return deg;
}


//===========================================================================
void SplineCurve::getWeights(std::vector<double>& weights) const
//===========================================================================
{
  int ncoefs = basis_.numCoefs();
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
void SplineCurve::interpolate(Interpolator& interpolator,
				int num_points,
				int dim,
				const double* param_start,
				const double* data_start)
//===========================================================================
{
    interpolator.interpolate(num_points, dim, param_start, data_start,
			     coefs_);
    basis_ = interpolator.basis();
    dim_ = dim;
    rational_ = false;
}


//===========================================================================
void SplineCurve::setParameterInterval(double t1, double t2)
//===========================================================================
{
    basis_.rescale(t1, t2);
}


//===========================================================================
double SplineCurve::nextSegmentVal(double par, bool forward, double tol) const
//===========================================================================
{
  if (!forward && par <= startparam())
    return startparam();
  
  if (forward && par >= endparam())
    return endparam();

  std::vector<double>::const_iterator knot;
  if (forward) {
    par+= fabs(tol);
    knot = std::upper_bound(basis().begin(),basis().end(),par);
    if (knot == basis().end())
      return endparam();
    else
      return *knot;
  }
  else {
    par -= fabs(tol);
    for (knot=basis().end()-1; knot>basis().begin(); --knot) { 
      if (*knot < par)
	return *knot;
    }
    return *(basis().begin());
  }
}

//===========================================================================
void SplineCurve::swap(SplineCurve& other)
//===========================================================================
{
    std::swap(dim_, other.dim_);
    std::swap(rational_, other.rational_);
    basis_.swap(other.basis_);
    coefs_.swap(other.coefs_);
    rcoefs_.swap(other.rcoefs_);
}


// Added by KMO for ICADA usage.
//===========================================================================
void SplineCurve::deform(const std::vector<double>& vec, int vdim)
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
  void SplineCurve::equalBdWeights(bool at_start)
//===========================================================================
{
  if (!rational_)
    return;  // Non-rational, all weights are equal to one. Nothing to do

  // REFERENCES : "Moebius Reparametrizations of Rational B-splines",
  //              E.T.Y. Lee and Miriam L. Lucian, CAGD 8 (1991) 213-215
  // The formula for the knots in the reparametrized curve is as follows:
  //
  // psi(t) = (alpha*t + beta)/(gamma*t + delta) =
  //          (w_0(t_1-t_0)(t_2-t)khi_0 + (t-t_0)(t_2(w_1-w_0)+w_0*t_1-w_1*t_0)khi_2)/
  //          (w_0(t_1-t_0)(t_2-t) + (t-t_0)(t_2(w_1-w_0)+w_0*t_1-w_1*t_0)),
  //
  //          where w_0 is end-weight, w_1 is neighbour-weight, t_0 & t_1 corr params,
  //          t_2 param in other end, khi_0 & khi_2 new interval (we're actually going
  //          to use input interval).
  //
  // The unknowns in psi are thus given as follows (setting khi_0 = t_0 & khi_2 = t_2):
  //
  // alpha = (t_2(w_1-w_0)+w_0*t_1-w_1*t_0)*t_2 - w_0(t_1-t_0)t_0
  // beta = t_0*t_2(t_2-t_0)(w_0-w_1)
  // gamma = (t_2-t_0)(w_1-w_0)
  // delta = w_0(t_1*t_2 - t_0*t_1) + w_1(t_0*t_0 - t_0*t_2)
  //
  // The new weights are given by the the formula:
  //
  // w_i' = w_i/(mu(t_(i+1))*(mu(t_(i+2))*...*(mu(t_(i+k-1)))),
  //
  // where mu(t) = gamma*t + delta
 
  int kdim = dim_ + 1;
  int in = basis_.numCoefs();
  int ik = basis_.order();

  double w_0 = (at_start) ? rcoefs_[dim_] : rcoefs_[kdim*(in-1)+dim_];
  double w_1 = (at_start) ? rcoefs_[kdim+dim_] : rcoefs_[kdim*(in-2)+dim_];

  // Check if the weights are already equal
  double eps = 1.0e-10;
  if (fabs(w_1 - w_0) < eps)
    return;

  vector<double>::iterator st = basis_.begin();
  double t_0 = (at_start) ? st[ik-1] : st[in];
  double t_1 = (at_start) ? st[ik] : st[in-1];
  double t_2 = (at_start) ? st[in] : st[ik-1];

  double alpha = (t_2*(w_1 - w_0) + w_0*t_1 - w_1*t_0)*t_2 - w_0*(t_1 - t_0)*t_0;
  double beta = t_0*t_2*(t_2 - t_0)*(w_0 - w_1);
  double gamma = (t_2 - t_0)*(w_1 - w_0);
  double delta = w_0*(t_1*t_2 - t_0*t_1) + w_1*(t_0*t_0 - t_0*t_2);

  vector<double> mu(ik+in);

  // We run through knots, reparametrizing.
  int ki, kj;
  for (ki = 0; ki < ik + in; ++ki)
    {
      mu[ki] = gamma*st[ki] + delta;
      if (fabs(mu[ki]) < eps)
	THROW("Not possible to perform Moebious reparametrization");
      double psi = (alpha*st[ki] + beta)/(mu[ki]);
      st[ki] = psi;
    }

  // The mu vector may contain large values. As mult of all weights with the
  // same factor does not alter the curve, we make sure we do not end up
  // with too high values.
  double min_mu = mu[0];
  double max_mu = mu[0];
  for (ki = 1; ki < (int)mu.size(); ++ki)
    {
      if (mu[ki] < min_mu)
	min_mu = mu[ki];
      if (mu[ki] > max_mu)
	max_mu = mu[ki];
    }
  double frac_mu = 2.0/(min_mu + max_mu);
  for (ki = 0; ki < (int)mu.size(); ++ki)
    {
      mu[ki] *= frac_mu;
    }

  // We then compute our new weights (and update homogeneous coordinates).
  for (ki = 0; ki < in; ++ki)
    {
      double prod = 1.0;
      for (kj = 1; kj < ik; ++kj)
	{
	  prod *= mu[ki+kj];
	}
      rcoefs_[ki*kdim+dim_] /= prod;
      for (kj = 0; kj < dim_; ++kj)
	{
	  rcoefs_[ki*kdim+kj] = coefs_[ki*dim_+kj]*rcoefs_[ki*kdim+dim_];
	}
    }

  // Update coefs_
  updateCoefsFromRcoefs();
 }

//===========================================================================
  void SplineCurve::representAsRational()
//===========================================================================
{
  if (rational_)
    return;   // This curve is already rational

  int in = basis_.numCoefs();
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
  void SplineCurve::setBdWeight(double wgt, bool at_start)
//===========================================================================
{
  if (!rational_)
    return;   // No weights 

  double wgt2 = (at_start) ? rcoefs_[dim_] : rcoefs_[rcoefs_.size()-1];
  double fac = wgt/wgt2;  // Weights are supposed to be positive and not zero
  for (size_t ki=0; ki<rcoefs_.size(); ++ki)
    rcoefs_[ki] *= fac;
}

//===========================================================================
  void SplineCurve::replaceEndPoint(Point pnt, bool at_start)
//===========================================================================
{
  if (at_start)
    makeKnotStartRegular();
  else
    makeKnotEndRegular();

  if (rational_)
    {
      vector<double>::iterator coefs = rcoefs_begin();
      int idx = (at_start) ? 0 : (numCoefs()-1)*(dim_+1);
      for (int ki=0; ki<dim_; ki++)
	coefs[idx+ki] = pnt[ki]*coefs[idx+dim_];
      updateCoefsFromRcoefs();
    }
  else
    {
      vector<double>::iterator coefs = coefs_begin();
      int idx = (at_start) ? 0 : (numCoefs()-1)*dim_;
      for (int ki=0; ki<dim_; ki++)
	coefs[idx+ki] = pnt[ki];
    }
}

//===========================================================================
void SplineCurve::translateCurve(const Point& dir)
//===========================================================================
{
  vector<double>::iterator c1 = 
    (rational_) ? rcoefs_begin() : coefs_begin();
  vector<double>::iterator c2 = 
    (rational_) ? rcoefs_end() : coefs_end();
  int dim = dim_ + rational_;
  for (; c1<c2; c1+=dim)
    {
      if (rational_)
	{
	  double wgt = c1[2];
	  c1[0] = (c1[0]/wgt + dir[0])*wgt;
	  c1[1] = (c1[1]/wgt + dir[1])*wgt;
	}
      else
	{
	  c1[0] += dir[0];
	  c1[1] += dir[1];
	}
    }
  if (rational_)
    updateCoefsFromRcoefs();

}

//===========================================================================
  void SplineCurve::translateSwapCurve(const Point& dir, double sgn, int pdir)
//===========================================================================
{
  vector<double>::iterator c1 = 
    (rational_) ? rcoefs_begin() : coefs_begin();
  vector<double>::iterator c2 = 
    (rational_) ? rcoefs_end() : coefs_end();
  int dim = dim_ + rational_;
  for (; c1<c2; c1+=dim)
    {
      if (rational_)
	{
	  double wgt = c1[2];
	  if (pdir == 1 || pdir == 3)
	    c1[0] = (sgn*c1[0]/wgt + dir[0])*wgt;
	  if (pdir == 2 || pdir == 3)
	    c1[1] = (sgn*c1[1]/wgt + dir[1])*wgt;
	}
      else
	{
	  if (pdir == 1 || pdir == 3)
	    c1[0] = (sgn*c1[0]+dir[0]);
	  if (pdir == 2 || pdir == 3)
	    c1[1] = (sgn*c1[1]+dir[1]);
	}
    }
  if (rational_)
    updateCoefsFromRcoefs();

}
//===========================================================================
shared_ptr<ElementaryCurve> SplineCurve::getElementaryCurve()
//===========================================================================
{
    if (is_elementary_curve_)
        return elementary_curve_;
    return shared_ptr<ElementaryCurve>();
}


//===========================================================================
void SplineCurve::setElementaryCurve(shared_ptr<ElementaryCurve> elcurve)
//===========================================================================
{
    elementary_curve_ = elcurve;
    is_elementary_curve_ = true;
}


//===========================================================================
bool SplineCurve::checkElementaryCurve()
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
void SplineCurve::updateCoefsFromRcoefs()
//===========================================================================
{
    int num_coefs = rcoefs_.size() / (dim_+1);
    coefs_.resize(num_coefs*dim_);
    SplineUtils::make_coef_array_from_rational_coefs(&rcoefs_[0],
					&coefs_[0],
					num_coefs,
					dim_);
}

//===========================================================================
bool SplineCurve::isAxisRotational(Point& centre, Point& axis, Point& vec,
				   double& angle)
//===========================================================================
{
  if (elementary_curve_.get())
    return elementary_curve_->isAxisRotational(centre, axis, vec, angle);
  else
    return false;
}

//===========================================================================
bool SplineCurve::isLinear(Point& dir, double tol)
//===========================================================================
{
  if (elementary_curve_.get())
    return elementary_curve_->isLinear(dir, tol);
  else
    {
      DirectionCone cone = directionCone();
      dir = cone.centre();

      // Check the cone angle for planarity
      return (cone.angle() < tol);
    }
}

//===========================================================================
bool SplineCurve::isInPlane(const Point& loc, const Point& axis,
			    double eps, Point& normal) const
//===========================================================================
{
  if (dim_ != loc.dimension())
    return false;

  // Select the vertex that gives the most significant distance to the location
  vector<double>::const_iterator cf = coefs_.begin();
  vector<double>::const_iterator cf2 = coefs_.end();

  double max_dist = 0.0;
  Point pos;
  for (; cf<cf2; cf+=dim_)
    {
      Point tmp(cf, cf+dim_);
      double dist = loc.dist(tmp);
      if (dist > max_dist)
	{
	  max_dist = dist;
	  pos = tmp;
	}
    }
  if (max_dist < eps)
    return false;  // This is not really true as the curve is degenerate, but
  // it is not possible to define the normal

  // Define plane normal
  normal = (pos - loc).cross(axis);
  normal.normalize();

  // Check if all coefficients lie in this plane
  cf = coefs_.begin();
  for (; cf<cf2; cf+=dim_)
    {
      Point tmp(cf, cf+dim_);
      double dist = (tmp - loc)*normal;
      if (dist > eps)
	return false;  // Not planar
    }
  return true;  // Planar
}

//===========================================================================
bool SplineCurve::isInPlane(const Point& norm,
			   double eps, Point& pos) const
//===========================================================================
{
  if (elementary_curve_.get())
    return elementary_curve_->isInPlane(norm, eps, pos);
  else
    {
      if (dim_ != norm.dimension())
	return false;
      
      // For each control vertex, check if the vector between this
      // vertex and the start vertex is perpendicular to the given
      // plane normal
      vector<double>::const_iterator cf = coefs_.begin();
      vector<double>::const_iterator cf2 = coefs_.end();

      pos = Point(cf, cf+dim_);
      for (cf+=dim_; cf<cf2; cf+=dim_)
	{
	  Point tmp(cf, cf+dim_);
	  Point vec = tmp - pos;
	  if (vec.length() < eps)
	    continue;

	  double ang = vec.angle(norm);
	  if (fabs(0.5*M_PI-ang) >= eps)
	    return false;
	}
      return true;
    }
	  
}

} // namespace Go








