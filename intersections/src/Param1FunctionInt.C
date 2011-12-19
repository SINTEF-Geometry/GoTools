//===========================================================================
//                                                                           
// File: Param1FunctionInt.C
//                                                                           
// Created: Tue Sep 21 12:25:30 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Param1FunctionInt.C,v 1.23 2006-11-03 14:15:12 jbt Exp $
//                                                                           
// Description: More or less taken directly from ParamCurveInt.C
//                                                                           
//===========================================================================


#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"


using std::vector;
using std::min;
using std::pair;
using std::sort;


namespace Go {


//===========================================================================
Param1FunctionInt::Param1FunctionInt(shared_ptr<ParamCurve> curve)
  : curve_(curve), deg_tol_(-1.0)
//===========================================================================
{
//     ASSERT(curve_->dimension() == 1);
    dim_ = 1; //dimension(); //1;
    parentcurve_ = 0;
}


//===========================================================================
Param1FunctionInt::Param1FunctionInt(shared_ptr<ParamCurve> curve,
				     ParamFunctionInt *parent)
  : curve_(curve), deg_tol_(-1.0)
//===========================================================================
{
    dim_ = 1; //curve_->dimension();
    parentcurve_ = dynamic_cast<Param1FunctionInt*>(parent);
    ASSERT(parent != 0);
}


//===========================================================================
Param1FunctionInt::~Param1FunctionInt()
//===========================================================================
{
}


//===========================================================================
Param1FunctionInt* Param1FunctionInt::getParam1FunctionInt()
//===========================================================================
{ 
    return this;
}


//===========================================================================
shared_ptr<ParamCurve> Param1FunctionInt::getParamCurve()
//===========================================================================
{
    return curve_;
}


//===========================================================================
shared_ptr<const ParamCurve> Param1FunctionInt::getParamCurve() const
//===========================================================================
{
    return curve_;
}


//===========================================================================
shared_ptr<ParamCurve>
Param1FunctionInt::getParentParamCurve(double& start, double& end)
//===========================================================================
{
  start = startparam();
  end = endparam();
  if (parentcurve_)
    return parentcurve_->getParentParamCurve();
  else
    return curve_;
}


//===========================================================================
shared_ptr<ParamCurve> Param1FunctionInt::getParentParamCurve()
//===========================================================================
{
    if (parentcurve_)
	return parentcurve_->getParentParamCurve();
    else
	return curve_;
}


//===========================================================================
int Param1FunctionInt::numParams() const
//===========================================================================
{
  return 1;
}


//===========================================================================
void Param1FunctionInt::getLengthAndWiggle(double *length, double *wiggle)
//===========================================================================
{
  int nsample = 5;
  int ki;
  double tl, tw;
  Point pt1, pt2, vec1, vec2, pt1_2d, pt2_2d;
  double tpar, tint;
  tint = (endparam() - startparam())/(double)(nsample-1);
  tpar = startparam();
  curve_->point(pt1, tpar);
  pt1_2d = Point(tpar, pt1[0]);
  for (ki=1, tpar+=tint, tl=0.0, tw=0.0; ki<nsample; ki++, tpar+=tint)
    {
      curve_->point(pt2, tpar);
      pt2_2d = Point(tpar, pt2[0]);
      tl += pt1_2d.dist(pt2_2d);
      vec2 = pt2_2d - pt1_2d;
      if (ki>1)
	tw += vec1.angle(vec2);
      pt1_2d = pt2_2d;
      vec1 = vec2;
    }

  *length = tl;
  *wiggle = tw;
}


//===========================================================================
bool Param1FunctionInt::hasInnerKnots(int pardir) const
//===========================================================================
{
  return false;
}


//===========================================================================
bool Param1FunctionInt::hasCriticalVals(int pardir) const
//===========================================================================
{
  return (segment_.size() > 0);
}

//===========================================================================
bool Param1FunctionInt::hasCriticalValsOrKnots(int pardir) const
//===========================================================================
{
  return (segment_.size() > 0);
}

//===========================================================================
bool Param1FunctionInt::canDivide(int pardir)
//===========================================================================
{
    // @@sbr Stop dividing if almost (tol-equal) flat?
    bool is_deg = isDegenerate(deg_tol_, 0, NULL);
    if (is_deg)
	return false;

    // This test must check whether the parameter interval of the curve
    // is long enough and probably also some other critera
    // @@@ VSK, for the time being ...
    double ta = startparam();
    double tb = endparam();
    double tc = 0.5*(ta+tb);

    double rel_par_res = 1e-10; // = epsge_->getRelParRes(); @@sbr
    if (min(fabs(tc-ta), fabs(tb-tc)) < rel_par_res)
	return false;
    else
	return true;
}


//===========================================================================
static bool compare_prio(std::pair<double, int> seg1, 
			 std::pair<double, int> seg2)
//===========================================================================
{
    // Comparison function used by getCriticalVals()

    if (seg1.second < seg2.second)
	return false;
    else
	return true;
}


//===========================================================================
vector<double> Param1FunctionInt::getCriticalVals(int pardir) const
//===========================================================================
{
  // Some values might be more critical than others (according to the
  // priority in segment_). Other might not be really critical. Should
  // all values be returned?
  // What about sorting? According to priority? According to the distance
  // from the middle of this curve? At least it exists more information
  // valuable for sorting here than in the Intersectors.

  // @@@ VSK. For the time being I sort according to priorites, but not
  // accoring to the parameter. I return all values.

  vector<double> vals;
  if (segment_.size() == 0)
    return vals;

  // Sort the segment array.
  vector<pair<double, int> > segment_copy = segment_;
  sort(segment_copy.begin(), segment_copy.end(), compare_prio);

  int ki;
  for (ki=0; ki<int(segment_copy.size()); ki++)
    vals.push_back(segment_copy[ki].first);
  
  return vals;
}

//===========================================================================
vector<double> Param1FunctionInt::getCriticalValsAndKnots(int pardir) const 
//===========================================================================
{
    vector<double> vals = getCriticalVals(pardir);
    return vals;
}

//===========================================================================
int Param1FunctionInt::getMeshSize(int dir)
//===========================================================================
{
  int meshsize = 5;
  if (dir != 0)
    return 1;
  else
    {
      if (mesh_.size() == 0)
	makeMesh(meshsize);
      return (int)mesh_.size()/dim_;
    }
}

//===========================================================================
vector<double>::iterator Param1FunctionInt::getMesh()
//===========================================================================
{
  int meshsize = 5;
  if (mesh_.size() == 0)
    makeMesh(meshsize);
  return mesh_.begin();
}

//===========================================================================
double Param1FunctionInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
  if (dir != 0 || idx < 0 || idx >= int(mesh_.size()))
    return 0.5*(startparam() + endparam());
  else
    return ((double)(idx)*endparam() + 
	    (double)(mesh_.size()-idx)*startparam())/(double)(mesh_.size());
}


//===========================================================================
vector<double> Param1FunctionInt::getInnerKnotVals(int pardir, bool sort) const 
//===========================================================================
{
  // Nothing for a general parametric curve
  vector<double> vals;
  return vals;
}


//===========================================================================
bool Param1FunctionInt::boundaryPoint(const double* par, double tol) const 
//===========================================================================
{
  ASSERT(par != NULL);

  double d1 = fabs(*par - startParam(1));
  double d2 = fabs(*par - endParam(1));

  return min(d1, d2) < tol;
}


//===========================================================================
void Param1FunctionInt::
subdivide(int pardir, double par, 
	  vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
	  vector<shared_ptr<ParamFunctionInt> >& bd_objs)
//===========================================================================
{
  double start = startparam();
  double end = endparam();

  // Perform subdivision
  // @@@ VSK, This calls should be made more effective, but can do for
  // the time being
  shared_ptr<ParamCurve> curve1, curve2;
  curve1 = shared_ptr<ParamCurve>(curve_->subCurve(start, par));
  curve2 = shared_ptr<ParamCurve>(curve_->subCurve(par, end));
  DEBUG_ERROR_IF(curve1.get()==0 || curve2.get()==0, "Error in subdivide");

  // Get the subdivision point
  // @@@ VSK, Can also be done more effective by simply fetching the
  // boundary coefficients
  shared_ptr<Point> subdivpt = shared_ptr<Point>(new Point(dim_));
  point(*subdivpt, &par); // We can not evaluate directly on the curve_.
//   curve1->point(*(subdivpt.get()), par);

  // Make intersection objects
  subdiv_objs.push_back(makeIntFunction(curve1));
  subdiv_objs.push_back(makeIntFunction(curve2));
  bd_objs.push_back(shared_ptr<ParamFunctionInt>
		    (new Param0FunctionInt((*subdivpt)[0], this)));
}


//===========================================================================
shared_ptr<Param1FunctionInt> 
   Param1FunctionInt::makeIntFunction(shared_ptr<ParamCurve> curve)
//===========================================================================
{
  shared_ptr<Param1FunctionInt> curve_int =
      shared_ptr<Param1FunctionInt>(new Param1FunctionInt(curve, this));
  return curve_int;
}


//===========================================================================
CompositeBox Param1FunctionInt::compositeBox() const
//===========================================================================
{
  return curve_->compositeBox();
}


//===========================================================================
bool Param1FunctionInt::monotone(Point& dir, double tol) const
//===========================================================================
{
    // @@sbr It seems we need more info to decide upon monotonicity.
    // Implement using sample approach!
    // Possibly store pts (and corr params).
    return false;
}


//===========================================================================
void Param1FunctionInt::
getBoundaryObjects(vector<shared_ptr<BoundaryFunctionInt> >& bd_objs)
//===========================================================================
{
  // First endpoint
  shared_ptr<Point> point1 = shared_ptr<Point>(new Point(dim_));
  double tmin = startparam();
  point(*point1, &tmin);
  shared_ptr<Param0FunctionInt> startpoint =
      shared_ptr<Param0FunctionInt>(new Param0FunctionInt((*point1)[0]));
  shared_ptr<BoundaryFunctionInt> bd1 =
      shared_ptr<BoundaryFunctionInt>
      (new BoundaryFunctionInt(startpoint, 0, startparam()));
  boundary_obj_.push_back(bd1);
  bd_objs.push_back(bd1);

  // Second endpoint
  shared_ptr<Point> point2 = shared_ptr<Point>(new Point(dim_));
  double tmax = endparam();
  point(*point2, &tmax);
  shared_ptr<Param0FunctionInt> endpoint =
      shared_ptr<Param0FunctionInt>(new Param0FunctionInt((*point2)[0]));
  shared_ptr<BoundaryFunctionInt> bd2 =
      shared_ptr<BoundaryFunctionInt>
      (new BoundaryFunctionInt(endpoint, 0, endparam()));
  boundary_obj_.push_back(bd2);
  bd_objs.push_back(bd2);
}


//===========================================================================
void Param1FunctionInt::assureInRange(double& t)
//===========================================================================
{
  double tmin = curve_->startparam();
  double tmax = curve_->endparam();
  if (t < tmin)
    t = tmin;
  else if (t > tmax)
    t = tmax;
}


//===========================================================================
int Param1FunctionInt::knotIntervalFuzzy(double& t, double tol) const
//===========================================================================
{
  double tmin = curve_->startparam();
  double tmax = curve_->endparam();
  if (t < tmin || fabs(t-tmin)<tol) {
    t = tmin;
    return 0;
  }
  else if (t > tmax || fabs(t-tmax)<tol) {
    t = tmax;
    return 1;
  }
  else
    return 0;

}


//===========================================================================
double Param1FunctionInt::nextSegmentVal(double par, bool forward) const
//===========================================================================
{
    // No inner knots expected to exist.
    if (forward) {
	return curve_->endparam();
    } else {
	return curve_->startparam();
    }
}


//===========================================================================
bool Param1FunctionInt::isDegenerate(double epsge, int dir, double *par)
//===========================================================================
{
    ASSERT(dir==0);

    if (deg_tol_ < 0.0 || deg_tol_ != epsge)
    {
	// Check degeneracy
	deg_tol_ = epsge;
    }

    CompositeBox box = compositeBox();
    Point low = box.low(0.0, 0.0);
    Point high = box.high(0.0, 0.0);
    int dim = box.dimension();
    for (int ki=0; ki<dim; ki++)
	if (high[ki] - low[ki] > epsge)
	    return false;
   
    return true;
}


//===========================================================================
void Param1FunctionInt::makeMesh(int size)
//===========================================================================
{
  mesh_.clear();
  mesh_.reserve(size*dimension());
  int ki;
  double ta = startparam();
  double tb = endparam();
  double tint = (tb - ta)/(double)(size-1);
  double par;
  Point pt;
  for (ki=0, par=ta; ki<size; ki++, par+=tint)
    {
      curve_->point(pt, par);
      mesh_.insert(mesh_.end(), pt.begin(), pt.end());
    }
      
}


//===========================================================================


} // namespace Go

