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

#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/utils/RotatedBox.h"


using std::min;
using std::pair;
using std::make_pair;
using std::vector;


namespace Go {


// //==========================================================================
// ParamCurveInt::ParamCurveInt(shared_ptr<ParamCurve> curve)
//   : curve_(curve)
// //==========================================================================
// {
//   dim_ = curve_->dimension();
// }

//===========================================================================
ParamCurveInt::ParamCurveInt(shared_ptr<ParamCurve> curve,
			     ParamGeomInt* parent)
    : ParamGeomInt(parent), curve_(curve), lw_set_(false)
//===========================================================================
{
  dim_ = curve_->dimension();
}

//===========================================================================
int ParamCurveInt::numParams() const
//===========================================================================
{
  return 1;
}

//===========================================================================
ParamCurveInt* ParamCurveInt::getParamCurveInt()
//===========================================================================
{
  return this;
}

//===========================================================================
shared_ptr<ParamCurve> ParamCurveInt::getParamCurve()
//===========================================================================
{
  return curve_;
}

//===========================================================================
shared_ptr<const ParamCurve> ParamCurveInt::getParamCurve() const
//===========================================================================
{
  return curve_;
}

//===========================================================================
shared_ptr<ParamCurve>
ParamCurveInt::getParentParamCurve(double& start, double& end)
//===========================================================================
{
    start = startparam();
    end = endparam();
    if (parent_ && parent_->numParams() == 1)
	{
	    ParamCurveInt *parentcv
		= dynamic_cast<ParamGeomInt*>(parent_)->getParamCurveInt();
	    if (parentcv == 0)
		return curve_;
	    else
		return parentcv->getParentParamCurve();
	}
    else
	return curve_;
}

//===========================================================================
shared_ptr<ParamCurve>
ParamCurveInt::getParentParamCurve()
//===========================================================================
{
  if (parent_ && parent_->numParams() == 1)
    {
	ParamCurveInt *parentcv
	    = dynamic_cast<ParamGeomInt*>(parent_)->getParamCurveInt();
      if (parentcv == 0)
	return curve_;
      else
	return parentcv->getParentParamCurve();
    }
  else
    return curve_;
}

//===========================================================================
bool ParamCurveInt::isDegenerate(double epsge, int dir, double *par)
//
//===========================================================================
{
    ASSERT(dir==0);
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
void ParamCurveInt::getLengthAndWiggle(double *length, double *wiggle)
//===========================================================================
{
    if (lw_set_)
    {
	length[0] = length_;
	wiggle[0] = wiggle_;
    }
    else
    {

	int nsample = 5;
	int ki;
	double tl, tw;
	double d1=0.0, d2=0.0;
	Point pt1, pt2, vec1, vec2;
	double tpar, tint;
	tint = (endparam() - startparam())/(double)(nsample-1);
	tpar = startparam();
	curve_->point(pt1, tpar);
	for (ki=1, tpar+=tint, tl=0.0, tw=0.0; ki<nsample; ki++, tpar+=tint)
	{
	    curve_->point(pt2, tpar);
	    d2 = pt1.dist(pt2);
	    tl += d2;
	    vec2 = pt2 - pt1;
	    if (ki>1 && d1 > 0.0 && d2 > 0.0)
		tw += vec1.angle(vec2);
	    pt1 = pt2;
	    vec1 = vec2;
	    d1 = d2;
	}

	length_ = *length = tl;
	wiggle_ = *wiggle = tw;
	lw_set_ = true;
    }
}

//===========================================================================
bool ParamCurveInt::hasInnerKnots(int pardir) const
//===========================================================================
{
  return false;
}

//===========================================================================
bool ParamCurveInt::hasCriticalVals(int pardir) const
//===========================================================================
{
  return (segment_.size() > 0);
}

//===========================================================================
void ParamCurveInt::setCriticalVal(int pardir, double par) 
//===========================================================================
{
    if (segment_.size() == 0 || par > segment_[segment_.size()-1].first)
	segment_.push_back(make_pair(par,1));
    else
    {
	double prev = startParam(0);
	for (int ki=0; ki<(int)(segment_.size()); ki++)
	{
	    if (prev < par && par < segment_[ki].first)
	    {
		segment_.insert(segment_.begin()+ki, make_pair(par,1));
		break;
	    }
	    prev = segment_[ki].first;
	}
    }
}

//===========================================================================
bool ParamCurveInt::hasCriticalValsOrKnots(int pardir) const 
//===========================================================================
{
  return (segment_.size() > 0 || checkPeriodicity(pardir) >= 1);
}

//===========================================================================
bool ParamCurveInt::canDivide(int pardir)
//===========================================================================
{
  // This test must check whether the parameter interval of the curve
  // is long enough and probably also some other critera
  // @@@ VSK, for the time being ...
  // @@@ VSK, This function could be moved one step up in the class
  // hierarchy if it continues to be similar the the surface function
  double ta = startparam();
  double tb = endparam();
  double tc = 0.5*(ta+tb);
  
  //return (tc != ta && tc != tb);
    double rel_par_res = 1e-10; // = epsge_->getRelParRes();
    if (min(fabs(tc-ta), fabs(tb-tc)) < rel_par_res)
	return false;
    else
	return true;
}

static bool compare_prio(std::pair<double, int> seg1, 
			 std::pair<double, int> seg2)
{
  if (seg1.second < seg2.second)
    return false;
  else
    return true;
}

//===========================================================================
vector<double> ParamCurveInt::getCriticalVals(int pardir) const
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
  if (segment_.size() > 0)
  {
      // Sort the segment array.
      vector<pair<double, int> > segment_copy = segment_;
      std::sort(segment_copy.begin(), segment_copy.end(), compare_prio);

      int ki;
      for (ki=0; ki<int(segment_copy.size()); ki++)
	  vals.push_back(segment_[ki].first);
  }
  
  // The seem of a periodic curve is also a critical parameter
  int per = checkPeriodicity(pardir);
  if (per >= 1)
      vals.insert(vals.begin(), startParam(pardir));

  return vals;
}

//===========================================================================
vector<double> ParamCurveInt::getCriticalValsAndKnots(int pardir) const
//===========================================================================
{
    vector<double> vals = getCriticalVals(pardir);
    return vals;
}

//===========================================================================
vector<double> ParamCurveInt::getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
  // Nothing for a general parametric curve
  vector<double> vals;
  return vals;
}


//===========================================================================
bool ParamCurveInt::boundaryPoint(const double* par, double tol) const
//===========================================================================
{
  ASSERT(par != NULL);

  double d1 = fabs(*par - startParam(0));
  double d2 = fabs(*par - endParam(0));

  return min(d1, d2) < tol;
}

//===========================================================================
void ParamCurveInt::subdivide(int pardir, double par, 
			      vector<shared_ptr<ParamGeomInt> >& subdiv_objs,
			      vector<shared_ptr<ParamGeomInt> >& bd_objs)
//===========================================================================
{
    double start = startparam();
    double end = endparam();
    double p_interval = end - start;

    // Perform subdivision
    // @@@ VSK, This calls should be made more effective, but can do for
    // the time being
    shared_ptr<ParamCurve> curve1, curve2;
    shared_ptr<ParamCurve> crv;
    crv = getParentParamCurve();
    int per = -1; // checkPeriodicity(0);
    if (per >= 1) {
	// Periodic curve. Make one non-periodic sub curve
	while (par >= end)
	    par -= p_interval;
	curve1 = shared_ptr<ParamCurve>(crv->subCurve(par, par+p_interval));
	DEBUG_ERROR_IF(curve1.get()==0, "Error in subdivide");
    } else {
	curve1 = shared_ptr<ParamCurve>(crv->subCurve(start, par));
	curve2 = shared_ptr<ParamCurve>(crv->subCurve(par, end));
	DEBUG_ERROR_IF(curve1.get()==0 || curve2.get()==0,
		       "Error in subdivide");
    }
    // Get the subdivision point
    // @@@ VSK, Can also be done more effective by simply fetching the
    // boundary coefficients
    shared_ptr<Point> subdivpt = shared_ptr<Point>(new Point(dim_));
    curve1->point(*(subdivpt.get()), par);

    // Make intersection objects
    subdiv_objs.push_back(makeIntObject(curve1));
    if (curve2.get() != 0)
	subdiv_objs.push_back(makeIntObject(curve2));
    bd_objs.push_back(shared_ptr<ParamPointInt>
		      (new ParamPointInt(subdivpt, this)));
}

//===========================================================================
shared_ptr<ParamCurveInt> 
   ParamCurveInt::makeIntObject(shared_ptr<ParamCurve> curve)
//===========================================================================
{
  shared_ptr<ParamCurveInt> curve_int =
      shared_ptr<ParamCurveInt>(new ParamCurveInt(curve, this));
  return curve_int;
}

//===========================================================================
CompositeBox  ParamCurveInt::compositeBox() const 
//===========================================================================
{
  return curve_->compositeBox();
}

//===========================================================================
DirectionCone  ParamCurveInt::directionCone() const
//===========================================================================
{
  return curve_->directionCone();
}

//===========================================================================
int 
ParamCurveInt::checkPeriodicity(int pardir) const
//===========================================================================
{
    // Don't know anything about periodicity
    return -1;  // Open, non-periodic curve
}

//===========================================================================
void  
ParamCurveInt::
getBoundaryObjects(std::vector<shared_ptr<BoundaryGeomInt> >& bd_objs)
//===========================================================================
{
    // Check periodicity
    int per = checkPeriodicity();
    int nmb_bd = (per < 1) ? 2 : 0;

  if (int(boundary_obj_.size()) == nmb_bd)
    bd_objs.insert(bd_objs.begin(), boundary_obj_.begin(), 
		   boundary_obj_.end());
  else
    {
      boundary_obj_.clear();

      if (per < 1)
      {
	  // First endpoint
	  shared_ptr<Point> point1 = shared_ptr<Point>(new Point(dim_));
	  curve_->point(*(point1.get()), startparam());
	  shared_ptr<ParamPointInt> startpoint = 
	      shared_ptr<ParamPointInt>(new ParamPointInt(point1, this));
	  shared_ptr<BoundaryGeomInt> bd1 = 
	      shared_ptr<BoundaryGeomInt>
	      (new BoundaryGeomInt(startpoint, 0, startparam()));
	  boundary_obj_.push_back(bd1);
	  bd_objs.push_back(bd1);

	  // Second endpoint
	  shared_ptr<Point> point2 = shared_ptr<Point>(new Point(dim_));
	  curve_->point(*(point2.get()), endparam());
	  shared_ptr<ParamPointInt> endpoint = 
	      shared_ptr<ParamPointInt>(new ParamPointInt(point2, this));
	  shared_ptr<BoundaryGeomInt> bd2 = 
	      shared_ptr<BoundaryGeomInt>
	      (new BoundaryGeomInt(endpoint, 0, endparam()));
	  boundary_obj_.push_back(bd2);
	  bd_objs.push_back(bd2);
      }
    }
}

//===========================================================================
int ParamCurveInt::getMeshSize(int dir)
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
vector<double>::iterator ParamCurveInt::getMesh()
//===========================================================================
{
  int meshsize = 5;
  if (mesh_.size() == 0)
    makeMesh(meshsize);
  return mesh_.begin();
}

//===========================================================================
double ParamCurveInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
    int nmesh = ((int)(mesh_.size()))/dim_;
    if (dir != 0 || idx < 0 || idx >= nmesh)
	return 0.5*(startparam() + endparam());
    else
	return ((double)(idx)*endparam() + 
	    (double)(nmesh-idx)*startparam())/(double)(nmesh);
}

//===========================================================================
void ParamCurveInt::makeMesh(int size) const
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
bool 
ParamCurveInt::isSpline()
//===========================================================================
{
    return false;
}

//===========================================================================
double
ParamCurveInt::getOptimizedConeAngle(Point& axis1, Point& axis2)
//===========================================================================
{
    return M_PI;
}

//===========================================================================
void ParamCurveInt::assureInRange(double& t)
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
int ParamCurveInt::knotIntervalFuzzy(double& t, double tol) const
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
RotatedBox ParamCurveInt::getRotatedBox(std::vector<Point>& axis) const
//===========================================================================
{
    int size = 100;
    makeMesh(size);
    RotatedBox box(&mesh_[0], dimension(), (int)mesh_.size(), 1, &axis[0]);
    return box;
}

//===========================================================================
void ParamCurveInt::axisFromEndpts(Point& axis) const
//===========================================================================
{
    double start = startParam(0);
    double end = endParam(0);
    Point pt1, pt2;
    point(pt1, &start);
    point(pt2, &end);
    axis = pt2 - pt1;
}

//===========================================================================
double ParamCurveInt::
nextSegmentVal(double par, bool forward, double tol) const
//===========================================================================
{
    return curve_->nextSegmentVal(par, forward, tol);
}


//===========================================================================


} // namespace Go
