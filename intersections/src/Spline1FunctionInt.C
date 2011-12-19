//===========================================================================
//                                                                           
// File: Spline1FunctionInt.C                                                
//                                                                           
// Created: Fri Sep 24 15:53:31 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Spline1FunctionInt.C,v 1.16 2006-09-01 12:25:34 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"


using std::vector;


namespace Go {


//===========================================================================
Spline1FunctionInt::Spline1FunctionInt(shared_ptr<ParamCurve> curve)
    : Param1FunctionInt(curve)
//===========================================================================
{
    spcv_= dynamic_pointer_cast<SplineCurve, ParamCurve>(curve);
}


//===========================================================================
Spline1FunctionInt::Spline1FunctionInt(shared_ptr<ParamCurve> curve,
				       ParamFunctionInt *parent)
    : Param1FunctionInt(curve, parent)
//===========================================================================
{
    spcv_= dynamic_pointer_cast<SplineCurve, ParamCurve>(curve);
}


//===========================================================================
bool Spline1FunctionInt::hasInnerKnots(int pardir) const
//===========================================================================
{
    return (spcv_->numCoefs() > spcv_->order());
}

//===========================================================================
bool Spline1FunctionInt::hasCriticalValsOrKnots(int pardir) const
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
vector<double> Spline1FunctionInt::
getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
  // First fetch all inner knots
  vector<double> vals;
  if (spcv_->numCoefs() == spcv_->order())
    return vals;

  std::vector<double>::const_iterator et = spcv_->basis().begin();
  int kk = spcv_->order();
  int kn = spcv_->numCoefs();
  vals.push_back(et[kk]);
  int ki;
  for (ki=kk+1; ki<kn; ki++)
    if (et[ki] > vals[vals.size()-1])
      vals.push_back(et[ki]);

  if (sort)
    {
      // Sort knot vector with respect to the distance from the midpoint
      // of the current parameter interval
      sort_distance compare(curve_->startparam(), curve_->endparam());
      std::sort(vals.begin(), vals.end(), compare);
    }

  return vals;
}


//===========================================================================
vector<double> Spline1FunctionInt::getCriticalValsAndKnots(int pardir) const
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
int Spline1FunctionInt::getMeshSize(int dir)
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return Param1FunctionInt::getMeshSize(dir);
  else
      return (dir == 0) ? spcv_->numCoefs() : 1;
}

//===========================================================================
vector<double>::iterator Spline1FunctionInt::getMesh()
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return Param1FunctionInt::getMesh();
  else
      return spcv_->coefs_begin();
}

//===========================================================================
double Spline1FunctionInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
  int meshsize = 5;
  if (spcv_->numCoefs() < meshsize)
      return Param1FunctionInt::paramFromMesh(dir, idx);
  else
      return (dir == 0) ? spcv_->basis().grevilleParameter(idx) : 0.0;
}

//===========================================================================
shared_ptr<Param1FunctionInt> 
Spline1FunctionInt::makeIntFunction(shared_ptr<ParamCurve> curve)
//===========================================================================
{
  shared_ptr<Spline1FunctionInt> curve_int =
    shared_ptr<Spline1FunctionInt>(new Spline1FunctionInt(curve, this));
  return curve_int;
}


//===========================================================================
int Spline1FunctionInt::knotIntervalFuzzy(double& t, double tol) const
//===========================================================================
{
  int i = spcv_->basis().knotIntervalFuzzy(t, tol);
  return i;
}


//===========================================================================
double Spline1FunctionInt::nextSegmentVal(double par, bool forward) const
//===========================================================================
{

  if (!forward && par <= curve_->startparam())
    return curve_->startparam();
  
  if (forward && par >= curve_->endparam())
    return curve_->endparam();

  /*
  std::vector<double>::const_iterator knot;
  if (forward) {
    knot = std::upper_bound(spcv_->basis().begin(),spcv_->basis().end(),par);
    return *knot;
  }
  else {
    for (knot=spcv_->basis().end()-1; knot>spcv_->basis().begin(); --knot) { 
      if (*knot < par)
	return *knot;
    }
    return *spcv_->basis().begin();
  }
  */

  int i = spcv_->basis().knotInterval(par);
  if (forward)
    return spcv_->basis().begin()[i+1];
  else
    return spcv_->basis().begin()[i];
}


//===========================================================================
bool Spline1FunctionInt::monotone(Point& dir, double tol) const
//===========================================================================
{
    // We start by constructing the derivative cv of our spline cv.
    shared_ptr<SplineCurve> der_cv(spcv_->derivCurve(1));

    // Even though the curve should be rational we are still
    // guaranteed that the curve is monotone if all coefs have the
    // same sign.

    // We then run through the coefs, checking if they are all above
    // or below 0.
    double threshold = 1.0e-08; //0.0; @@sbr Suitable value?
    threshold = std::max(threshold, tol);
    vector<double>::const_iterator scoef = der_cv->coefs_begin();
    // By using a threshold we set a tolerance belt around 0.0.
    bool above = (scoef[0] > -threshold);
    bool below = (scoef[0] < threshold);
    int ki;
    int kn = der_cv->numCoefs();
    for (ki = 1; ki < kn; ++ki) {
	if (scoef[ki] < -threshold) // We do not want 0.0 derivative ...
	    above = false;
	if (scoef[ki] > threshold)
	    below = false;
    }

    // dir = Point(1.0);
    dir.resize(1); // is it perhaps this that is intended by the
		   // above, commented-out line?
    dir[0] = 1.0; 

    return (above || below);
}


//===========================================================================


} // namespace Go

