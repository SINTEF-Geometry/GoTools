//===========================================================================
//                                                                           
// File: Spline2FunctionInt.C                                                
//                                                                           
// Created: Mon Sep 27 16:14:46 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Spline2FunctionInt.C,v 1.30 2006-12-19 13:49:02 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Spline2FunctionInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/SurfaceCreators.h"


using std::vector;
using std::max;

namespace Go
{


//===========================================================================
Spline2FunctionInt::Spline2FunctionInt(shared_ptr<SplineSurface> surf)
    : Param2FunctionInt(surf), spsf_(surf)
//===========================================================================
{
    ASSERT(spsf_.get() != NULL);
}


//===========================================================================
Spline2FunctionInt::Spline2FunctionInt(shared_ptr<SplineSurface> surf,
				       Param2FunctionInt *parent)
    : Param2FunctionInt(surf, parent), spsf_(surf)
//===========================================================================
{
    ASSERT(spsf_.get() != NULL);
}


//===========================================================================
bool Spline2FunctionInt::hasInnerKnots(int pardir) const
//===========================================================================
{
    return (pardir == 1) ? (spsf_->numCoefs_v() > spsf_->order_v()) :
	(spsf_->numCoefs_u() > spsf_->order_u());
}


//===========================================================================
bool Spline2FunctionInt::hasCriticalValsOrKnots(int pardir) const
//===========================================================================
{
    return (hasCriticalVals(pardir) || hasInnerKnots(pardir));
}


//===========================================================================
struct sort_distance
//===========================================================================
{
    double mid;
    sort_distance(double start, double end)
    { mid = 0.5*(start+end); }

    bool operator()(double a, double b) const
    {
	return (fabs(a-mid) < fabs(b-mid));
    }
};


//===========================================================================
vector<double> Spline2FunctionInt::getInnerKnotVals(int pardir, bool sort) const 
//===========================================================================
{
  // First fetch all inner knots
  vector<double> vals;
  const BsplineBasis basis = (pardir == 1) ? spsf_->basis_v() : spsf_->basis_u();
  if (basis.numCoefs() == basis.order())
      return vals;

  std::vector<double>::const_iterator et = basis.begin();
  int kk = basis.order();
  int kn = basis.numCoefs();
  vals.push_back(et[kk]);
  int ki;
  for (ki=kk+1; ki<kn; ki++)
    if (et[ki] > vals[vals.size()-1])
      vals.push_back(et[ki]);

  if (sort)
    {
      // Sort knot vector with respect to the distance from the midpoint
      // of the current parameter interval
      sort_distance compare(basis.startparam(), basis.endparam());
      std::sort(vals.begin(), vals.end(), compare);
    }

  return vals;
}

//===========================================================================
vector<double> Spline2FunctionInt::getCriticalValsAndKnots(int pardir) const 
//===========================================================================
{
    vector<double> critical = getCriticalVals(pardir);
    vector<double> knots = getInnerKnotVals(pardir, false);
    vector<double> vals;
    double ta = startParam(pardir);
    double tb = endParam(pardir);
    vals.push_back(ta);
    int ki, kj;
    for (ki=0, kj=0; ;)
    {
	if (ki >= int(critical.size()) && kj >= int(knots.size()))
	    break;
	else if (ki >= int(critical.size()))
	{
	    if (knots[kj] > vals[vals.size()-1])
		vals.push_back(knots[kj]);
	    kj++;
	}
	else if (kj >= int(knots.size()))
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
	else if (critical[ki] < knots[kj])
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
	else
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
    }
    if (tb > vals[vals.size()-1])
	vals.push_back(tb);

    return vals;
}


//===========================================================================
double Spline2FunctionInt::startParam(int pardir) const
//===========================================================================
{
    return (pardir == 1) ? spsf_->startparam_v() : spsf_->startparam_u();
}


//===========================================================================
double Spline2FunctionInt::endParam(int pardir) const
//===========================================================================
{
    return (pardir == 1) ? spsf_->endparam_v() : spsf_->endparam_u();
}


//===========================================================================
shared_ptr<Param2FunctionInt>
Spline2FunctionInt::makeIntFunction(shared_ptr<ParamSurface> surf)
//===========================================================================
{
    shared_ptr<SplineSurface> spline_surf
	= dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
    shared_ptr<Spline2FunctionInt> surf_int =
	shared_ptr<Spline2FunctionInt>
	(new Spline2FunctionInt(spline_surf, this));
    return surf_int;
}


//===========================================================================
int Spline2FunctionInt::
knotIntervalFuzzy(int pardir, double& t, double tol) const
//===========================================================================
{
  const BsplineBasis basis
      = (pardir == 1) ? spsf_->basis_v() : spsf_->basis_u();
  int i = basis.knotIntervalFuzzy(t, tol);
  return i;
}


//===========================================================================
double Spline2FunctionInt::
nextSegmentVal(int pardir, double par, bool forward) const
//===========================================================================
{
  if (!forward && par <= startParam(pardir))
    return startParam(pardir);
  
  if (forward && par >= endParam(pardir))
    return endParam(pardir);

  /*
  std::vector<double>::const_iterator knot;
  if (forward) {
    knot = std::upper_bound(spsf_->basis().begin(),spsf_->basis().end(),par);
    return *knot;
  }
  else {
    for (knot=spsf_->basis().end()-1; knot>spsf_->basis().begin(); --knot) { 
      if (*knot < par)
	return *knot;
    }
    return *spsf_->basis().begin();
  }
  */

  const BsplineBasis basis = (pardir == 1) ? spsf_->basis_v() : spsf_->basis_u();

  int i = basis.knotInterval(par);
  if (forward)
      return basis.begin()[i+1];
  else
      return basis.begin()[i];
}


//===========================================================================
void  
Spline2FunctionInt::getBoundaryObjects(vector<shared_ptr<BoundaryFunctionInt> >& bd_objs)
//===========================================================================
{
    // We need to extract the four bd cvs (i.e. bd functions).
    int ki;
    for (ki = 0; ki < 4; ++ki) {
	int pardir;
	double tpar;
	shared_ptr<SplineCurve> bd_cv(extractBdCurve(*spsf_, ki, pardir, tpar));
	// The boundary curve should then be put in a Spline1FunctionInt.
	shared_ptr<Spline1FunctionInt> spline1_int
	    (new Spline1FunctionInt(bd_cv));

	shared_ptr<BoundaryFunctionInt> bd_obj(new BoundaryFunctionInt(spline1_int,
								       pardir, tpar));

	boundary_obj_.push_back(bd_obj);
	bd_objs.push_back(bd_obj);
    }
}


//===========================================================================
shared_ptr<SplineSurface> Spline2FunctionInt::surface3D()
//===========================================================================
{
    shared_ptr<SplineSurface> sf_3d = SurfaceCreators::insertParamDomain(*spsf_);

    return sf_3d;
}


//===========================================================================
bool Spline2FunctionInt::monotone(Point& dir, double tol) const
//===========================================================================
{
    ASSERT(spsf_->dimension() == 1); // Well, probably checked already.
    // We start by constructing the derivative of our spline function.
    // (dS/du, dS/dv)
    shared_ptr<SplineSurface> grad(createGradSurface());

    // Since we're only interested in the sign of the coefficients,
    // and all the weights are assumed to be positive, we may look at
    // the spatial coefs only.  We then run through the coefs,
    // checking if they are all above or below 0.
    //int ki;
    // We do not want values too close to 0.0 (without being 0.0).

    // Calculate maximal value of the coefficients. This will be used
    // to set the threshold, which is a suitable fraction of the
    // maximal value.
    typedef vector<double>::iterator iter;
    double maxval = 0.0;
    for (iter it = grad->coefs_begin(); it != grad->coefs_end(); ++it) {
	if (maxval < fabs(*it)) {
	    maxval = fabs(*it);
	}
    }

    //double threshold = 1.0e-15;
    //double threshold = tol;
    double threshold = 0.001 * maxval;
    vector<double>::iterator it = grad->coefs_begin();
    // The direction cone expects all elements to be non-zero.
    // We're accept non-strict monotonicity, as long as it is not
    // total (along all the sf).
    vector<double> non_zero_coefs;
    while (it != grad->coefs_end()) {
	if (fabs(it[0]) > threshold || fabs(it[1]) > threshold) {
	    non_zero_coefs.push_back(it[0]);
	    non_zero_coefs.push_back(it[1]);
	}
	it += 2;
    }

    // We then compute the directioncone for the gradient surface.
    // S(u,v) is monotone iff there exist a direction (dir_u, dir_v)
    // such that S(u'+t*dir_u, v'+t*dir_v) is monotone.  This is
    // equivalent to the gradient surface lying in a halfplane,
    // i.e. au + bv > 0.0.  Which again is equivalent to the box of
    // the polar coordinates spanning less than pi.
    int dim = grad->dimension();
    DirectionCone cone;
//     double cone_zero_tol = 1e-10;
//     cone.setZeroTol(cone_zero_tol);
    int nmb = (int)(non_zero_coefs.end() - non_zero_coefs.begin());
    if (nmb == 0) {
	// we cannot set cone from empty array
	MESSAGE("Spline2FunctionInt::monotone() -> encountered zero array, "
		"could not set DirectionCone.");
	dir = Point(0.0, 0.0);
	return true; // We assert monotonity as all coefs were below threshold.
    }
    cone.setFromArray(&(non_zero_coefs[0]), &(non_zero_coefs[0]) + nmb, dim);
    bool monotone = (!cone.greaterThanPi());
    dir = cone.centre();

    return monotone;
}


//===========================================================================
shared_ptr<SplineSurface> Spline2FunctionInt::createGradSurface() const
//===========================================================================
{
    // We differentiate the surface in both parameter directions.
    shared_ptr<SplineSurface> S_u(spsf_->derivSurface(1, 0));
    shared_ptr<SplineSurface> S_v(spsf_->derivSurface(0, 1));

#ifdef INTERSECTIONS_DEBUG
    if (0) { // We write to file (u, v, S_u) & (u, v, S_v).
	vector<double>::const_iterator max_iter =
	    max_element(S_u->coefs_begin(), S_u->coefs_end());
	vector<double>::const_iterator min_iter =
	    min_element(S_u->coefs_begin(), S_u->coefs_end());
	double max_int = max(S_u->endparam_u() - S_u->startparam_u(),
			     S_u->endparam_v() - S_u->startparam_v());
	double max_abs_val = max(-min_iter[0], max_iter[0]);
	if (max_abs_val < 2.0*max_int) {
	    double factor = 2.0*max_int/max_abs_val;
	    vector<double>::iterator iter = S_u->coefs_begin();
	    while (iter != S_u->coefs_end()) {
		iter[0] *= factor;
		++iter;
	    }
	    if (S_u->rational()) {
		iter = S_u->rcoefs_begin();
		while (iter != S_u->rcoefs_end()) {
		    iter[0] *= factor;
		    iter += 2;
		}
	    }
	}
	// Then the next sf
	max_iter = max_element(S_v->coefs_begin(), S_v->coefs_end());
	min_iter = min_element(S_v->coefs_begin(), S_v->coefs_end());
	max_int = max(S_v->endparam_u() - S_v->startparam_u(),
		      S_v->endparam_v() - S_v->startparam_v());
	max_abs_val = max(-min_iter[0], max_iter[0]);
	if (max_abs_val < 2.0*max_int) {
	    double factor = 2.0*max_int/max_abs_val;
	    vector<double>::iterator iter = S_v->coefs_begin();
	    while (iter != S_v->coefs_end()) {
		iter[0] *= factor;
		++iter;
	    }
	    if (S_v->rational()) {
		iter = S_v->rcoefs_begin();
		while (iter != S_v->rcoefs_end()) {
		    iter[0] *= factor;
		    iter += 2;
		}
	    }
	}

	shared_ptr<SplineSurface> S_u_3d = SurfaceCreators::insertParamDomain(*S_u);
	shared_ptr<SplineSurface> S_v_3d = SurfaceCreators::insertParamDomain(*S_v);
	// We also create the separating 0.0-sf.
	vector<double> lin_knots_u(4, S_u->startparam_u());
	lin_knots_u[2] = lin_knots_u[3] = S_u->endparam_u();
	vector<double> lin_knots_v(4, S_u->startparam_v());
	lin_knots_v[2] = lin_knots_v[3] = S_u->endparam_v();
	vector<double> coefs(4, 0.0);
	shared_ptr<SplineSurface> sf_0(new SplineSurface(2, 2, 2, 2,
							 lin_knots_u.begin(), lin_knots_v.begin(),
							 coefs.begin(), 1));
	shared_ptr<SplineSurface> sf_0_3d = SurfaceCreators::insertParamDomain(*sf_0);
	std::ofstream debug("tmp/debug.g2");
	S_u_3d->writeStandardHeader(debug);
	S_u_3d->write(debug);
	S_v_3d->writeStandardHeader(debug);
	S_v_3d->write(debug);
	sf_0_3d->writeStandardHeader(debug);
	sf_0_3d->write(debug);
    }
#endif // INTERSECTIONS_DEBUG

    // We next make sure that they share order.
    int order1 = max(S_u->order_u(), S_v->order_u());
    int order2 = max(S_u->order_v(), S_v->order_v());
    S_u->raiseOrder(order1 - S_u->order_u(), order2 - S_u->order_v());
    S_v->raiseOrder(order1 - S_v->order_u(), order2 - S_v->order_v());

    // We merge S_u & S_v into a 2D-surface.
    int kdim = 1 + spsf_->rational();
    int nmb_coefs = S_u->numCoefs_u()*S_u->numCoefs_v();
    int new_kdim = 2 + spsf_->rational();
    vector<double> new_coefs(nmb_coefs*new_kdim);
    vector<double>::const_iterator iter_u = (spsf_->rational()) ?
	S_u->rcoefs_begin() : S_u->coefs_begin();
    vector<double>::const_iterator iter_v = (spsf_->rational()) ?
	S_v->rcoefs_begin() : S_v->coefs_begin();
    for (int ki = 0; ki < nmb_coefs; ++ki) {
	new_coefs[ki*new_kdim] = iter_u[ki*kdim];
	new_coefs[ki*new_kdim+1] = iter_v[ki*kdim];
	if (spsf_->rational()) {
	    new_coefs[ki*new_kdim+2] = iter_u[ki*kdim+1];
	}
    }

    // Finally we create the gradient surface using our new vector of 2D coefs.
    shared_ptr<SplineSurface> grad_sf(new SplineSurface(S_u->numCoefs_u(),
							S_u->numCoefs_v(),
							S_u->order_u(),
							S_u->order_v(),
							S_u->basis_u().begin(),
							S_u->basis_v().begin(),
							new_coefs.begin(), 2,
							spsf_->rational()));
				      
    return grad_sf;
}


//===========================================================================


} // namespace Go

