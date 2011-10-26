//===========================================================================
//                                                                           
// File: SplineCurveInt.C 
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: SplineCurveInt.C,v 1.21 2006-07-07 14:46:50 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/RotatedBox.h"


using std::vector;
using std::shared_ptr;


namespace Go {


//===========================================================================
SplineCurveInt::SplineCurveInt(shared_ptr<ParamCurve> curve,
			       ParamGeomInt *parent)
  : ParamCurveInt(curve, parent)
//===========================================================================
{
  spcv_= std::dynamic_pointer_cast<SplineCurve, ParamCurve>(curve);
  if (checkPeriodicity(0) >= 1)
  {
      spcv_ = (shared_ptr<SplineCurve>)(spcv_->subCurve(spcv_->startparam(),
							spcv_->endparam()));
      curve_ = (shared_ptr<ParamCurve>)(spcv_);
  }
}


//===========================================================================
SplineCurveInt::SplineCurveInt(shared_ptr<ParamCurve> curve)
  : ParamCurveInt(curve)
//===========================================================================
{
  spcv_= std::dynamic_pointer_cast<SplineCurve, ParamCurve>(curve);
}


//===========================================================================
int 
SplineCurveInt::checkPeriodicity(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0);
    int per = analyzePeriodicity(*(spcv_.get()));
    return per;
}


//===========================================================================
bool SplineCurveInt::hasInnerKnots(int pardir) const
//===========================================================================
{
  return (spcv_->numCoefs() > spcv_->order());
}


//===========================================================================
struct sort_distance
//===========================================================================
{
    // Functor used to compare elements of a vector for sorting in
    // getInnerKnotVals()

    double mid;
    sort_distance(double start, double end)
    { mid = 0.5*(start+end); }

    bool operator()(double a, double b) const
    {
	return (fabs(a-mid) < fabs(b-mid));
    }
};


//===========================================================================
vector<double> SplineCurveInt::getInnerKnotVals(int pardir, bool sort) const 
//===========================================================================
{
    // First fetch all inner knots
    vector<double> vals;
    if (spcv_->numCoefs() == spcv_->order())
	return vals;

    vector<double>::const_iterator et = spcv_->basis().begin();
    int kk = spcv_->order();
    int kn = spcv_->numCoefs();
    vals.push_back(et[kk]);
    for (int ki = kk+1; ki < kn; ki++) {
	if (et[ki] > vals.back()) {
	    vals.push_back(et[ki]);
	}
    }

    if (sort) {
	// Sort knot vector with respect to the distance from the midpoint
	// of the current parameter interval
	sort_distance compare(curve_->startparam(), curve_->endparam());
	std::sort(vals.begin(), vals.end(), compare);
    }

    return vals;
}


//===========================================================================
bool SplineCurveInt::hasCriticalValsOrKnots(int pardir) const 
//===========================================================================
{
  ASSERT(pardir == 0 || pardir == 1);
  bool critical = hasCriticalVals(pardir);
  return (hasInnerKnots(pardir) || critical);
}


//===========================================================================
vector<double> SplineCurveInt::getCriticalValsAndKnots(int pardir) const 
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
std::shared_ptr<ParamCurveInt> 
   SplineCurveInt::makeIntObject(shared_ptr<ParamCurve> curve)
//===========================================================================
{
  shared_ptr<SplineCurveInt> curve_int =
    shared_ptr<SplineCurveInt>(new SplineCurveInt(curve, this));
  return curve_int;
}


//===========================================================================
int SplineCurveInt::getMeshSize(int dir)
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return ParamCurveInt::getMeshSize(dir);
  else
      return (dir == 0) ? spcv_->numCoefs() : 1;
}


//===========================================================================
vector<double>::iterator SplineCurveInt::getMesh()
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return ParamCurveInt::getMesh();
  else
      return spcv_->coefs_begin();
}


//===========================================================================
double SplineCurveInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return ParamCurveInt::paramFromMesh(dir, idx);
  else
      return (dir == 0) ? spcv_->basis().grevilleParameter(idx) : 0.0;
}


//===========================================================================
bool 
SplineCurveInt::isSpline()
//===========================================================================
{
    return true;
}


//===========================================================================
double
SplineCurveInt::getOptimizedConeAngle(Point& axis1, Point& axis2)
 //===========================================================================
{
    // We are making a plane spanned by the centre axis of given
    // cones. Then we project each tangent to this plane and compute
    // the angle beetween these projections and the centre of the
    // cone. Based on the sisl function s1796.

    int in = spcv_->numCoefs();
    vector<double>::iterator coefs = spcv_->coefs_begin();
    double angle = 0.0;

    vector<double>::iterator it1;
    int ki;
    Point endpt[2];   // End coefficients of the current segment
    Point diff;       // Difference vector between coefficients
    for (it1=coefs, ki=0; ki<in-1; ki++, it1+=dim_)
    {
	// Here we make an aproximative tangent to the curve
	// using the control polygon. The tangent is normalized	    
	endpt[0].resize(dim_);
	endpt[0].setValue(&it1[0]);
	endpt[1].resize(dim_);
	endpt[1].setValue(&it1[dim_]);
	diff = endpt[1] - endpt[0];
	double len = diff.normalize_checked();
	if (len == 0.0)
	    diff = axis1;

	double t1 = axis1*diff;
	double t2 = axis2*diff;
	if (t2 <= 0.0)
	    continue;

	double tang = t1/sqrt(t1*t1 + t2*t2);
	tang = std::min(tang, 1.0);
	tang = std::max(tang, -1.0);
	tang = acos(tang);
	angle = std::max(angle, tang);
    }

    return angle;
}


//===========================================================================
int SplineCurveInt::knotIntervalFuzzy(double& t, double tol) const
//===========================================================================
{
  int i = spcv_->basis().knotIntervalFuzzy(t, tol);
  return i;
}


//===========================================================================
RotatedBox SplineCurveInt::getRotatedBox(std::vector<Point>& axis) const
//===========================================================================
{
    RotatedBox box(spcv_->coefs_begin(), dimension(), spcv_->numCoefs(), 
		   1, &axis[0]);
    return box;
}


//===========================================================================
double SplineCurveInt::
nextSegmentVal(double par, bool forward, double tol) const
//===========================================================================
{
    return curve_->nextSegmentVal(par, forward, tol);
}


//===========================================================================


} // namespace Go

