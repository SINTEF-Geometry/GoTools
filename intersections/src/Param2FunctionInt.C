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

#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"


using std::vector;
using std::cerr;
using std::endl;
using std::min;
using std::pair;
using std::sort;


namespace Go {


//===========================================================================
Param2FunctionInt::Param2FunctionInt(shared_ptr<ParamSurface> surf)
    : surf_(surf), deg_tol_(-1.0),
      domain_(surf_->containingDomain()), temp_point_array_(3) // Only 1 sf.
//===========================================================================
{
  ASSERT(surf_->dimension() == 1);
  dim_ = 1;
  parentfunction_ = 0;
}


//===========================================================================
Param2FunctionInt::Param2FunctionInt(shared_ptr<ParamSurface> surf,
				     Param2FunctionInt *parent)
  : surf_(surf), deg_tol_(-1.0),
    domain_(surf_->containingDomain()), temp_point_array_(3) // Only 1 sf.
//===========================================================================
{
  dim_ = surf_->dimension();
  parentfunction_ = parent;
}


//===========================================================================
Param2FunctionInt* Param2FunctionInt::getParam2FunctionInt()
//===========================================================================
{
    return this;
}


//===========================================================================
shared_ptr<ParamCurve> 
Param2FunctionInt::getIsoCurve(double param_start, 
			       double param_end, 
			       double isoval, 
			       bool pardir_is_u) const
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > const_curves = 
	surf_->constParamCurves(isoval, pardir_is_u);
    if (const_curves.size() != 1) {
	cerr << "Design error!!\n"
	     << "Code does not support fragmented isocurves, such as\n"
	     << "might arise from trimmed surfaces with holes.\n"
	     << "You now continue completely at your own responsibility!"
	     << endl;
    }
    return shared_ptr<ParamCurve>(const_curves[0]->subCurve(param_start,
							    param_end));
}


//===========================================================================
shared_ptr<ParamCurve> 
Param2FunctionInt::getConstantParameterCurve(int pardir, double par)
//===========================================================================
{
    // Get the current parameter domain
    vector<double> mima = getMima();

    // End-parameter values of the curve
    double param1[2], param2[2];
    int idx2 = 1-pardir;
    param1[pardir] = param2[pardir] = par;
    param1[idx2] = mima[2*idx2];
    param2[idx2] = mima[2*idx2+1];

    // Make constant parameter curve as a curve-on-surface curve
    Point pp1(param1[0], param1[1]), pp2(param2[0], param2[1]);
    shared_ptr<SplineCurve> pcrv =
	shared_ptr<SplineCurve>(new SplineCurve(pp1, param1[idx2],
						pp2, param2[idx2]));
    shared_ptr<ParamCurve> constcrv =
	shared_ptr<ParamCurve>(new CurveOnSurface(getParamSurface(),
						  pcrv, true));

    return constcrv;
}


//===========================================================================
shared_ptr<ParamSurface>
Param2FunctionInt::getParentParamSurface(RectDomain& domain)
//===========================================================================
{
  domain = domain_;
  if (parent_ && parent_->numParams() == 2)
  {
      Param2FunctionInt *parentsf
	  = dynamic_cast<ParamFunctionInt*>(parent_)->getParam2FunctionInt();
      if (parentsf == 0)
	  return surf_;
      else
	  return parentsf->getParentParamSurface();
  }
  else
    return surf_;
}


//===========================================================================
shared_ptr<ParamSurface> 
Param2FunctionInt::getParentParamSurface()
//===========================================================================
{
  if (parent_ && parent_->numParams() == 2)
  {
      Param2FunctionInt *parentsf
	  = dynamic_cast<ParamFunctionInt*>(parent_)->getParam2FunctionInt();
      if (parentsf == 0)
	  return surf_;
      else
	  return parentsf->getParentParamSurface();
  }
  else
    return surf_;
}


//===========================================================================
shared_ptr<ParamSurface> Param2FunctionInt::getParamSurface()
//===========================================================================
{
    return surf_;
}


//===========================================================================
shared_ptr<const ParamSurface> Param2FunctionInt::getParamSurface() const
//===========================================================================
{
    return surf_;
}


//===========================================================================
int Param2FunctionInt::numParams() const
//===========================================================================
{
  return 2;
}


//===========================================================================
vector<double> Param2FunctionInt::getMima() const
//===========================================================================
{
    vector<double> mima(4);
    mima[0] = domain_.umin();
    mima[1] = domain_.umax();
    mima[2] = domain_.vmin();
    mima[3] = domain_.vmax();
    return mima;
}


//===========================================================================
void Param2FunctionInt::getLengthAndWiggle(double *length, double *wiggle)
//===========================================================================
{
    int nsample_u = 5; // We sample in a nsample*nsample grid.
    int nsample_v = 5;
    int ki, kj;
    double tl1 = 0.0, tw1 = 0.0, tl2 = 0.0, tw2 = 0.0;
    Point pt, vec1_1, vec1_2, vec2_1, vec2_2;
    double umin = startParam(0);
    double ustep = (endParam(0) - umin)/(nsample_u - 1);
    double vmin = startParam(1);
    double vstep = (endParam(1) - vmin)/(nsample_v - 1);
    vector<Point> pts; // We store our sampled pts in a vector for
                       // easy access.  We look upon the function as
                       // embedded in 3D: (u, v, f(u, v)).
    for (ki = 0; ki < nsample_v; ++ki) {
	double vpar = vmin + ki*vstep;
	for (kj = 0; kj < nsample_u; ++kj) {
	    double upar = umin + kj*ustep;
	    surf_->point(pt, upar, vpar);
	    pts.push_back(Point(upar, vpar, pt[0]));
	}
    }

    // We then compute the total dist and wiggle in both parameter directions.
    for (ki = 0; ki < nsample_u; ++ki) {
	for (kj = 0; kj < nsample_v; ++kj) {
	    if (ki < nsample_u - 1) {
		tl1 += pts[kj*nsample_u+ki].dist(pts[kj*nsample_u+ki+1]);
		if (ki > 0) {
		    vec1_1 = pts[kj*nsample_u+ki] - pts[kj*nsample_u+ki-1];
		    vec1_2 = pts[kj*nsample_u+ki+1] - pts[kj*nsample_u+ki];
		    if (vec1_1.length2() != 0.0 && vec1_2.length2() != 0.0)
			tw1 += vec1_1.angle(vec1_2);
		}
	    }
	    if (kj < nsample_v - 1) {
		tl2 += pts[kj*nsample_u+ki].dist(pts[(kj+1)*nsample_u+ki]);
		if (kj > 0) {
		    vec2_1 = pts[kj*nsample_u+ki] - pts[(kj-1)*nsample_u+ki];
		    vec2_2 = pts[(kj+1)*nsample_u+ki] - pts[kj*nsample_u+ki];
		    if (vec2_1.length2() != 0.0 && vec2_2.length2() != 0.0)
			tw2 += vec2_1.angle(vec2_2);
		}
	    }
	}
    }

    // To minimize dependency on # samples we divide by nsample (in case
    // function is used as a measure of sf characteristic).
    length[0] = tl1/nsample_v;
    wiggle[0] = tw1/nsample_v;
    length[1] = tl2/nsample_u;
    wiggle[1] = tw2/nsample_u;
}


//===========================================================================
bool Param2FunctionInt::hasInnerKnots(int pardir) const
//===========================================================================
{
  return false;
}


//===========================================================================
bool Param2FunctionInt::hasCriticalVals(int pardir) const
//===========================================================================
{
  ASSERT(pardir == 0 || pardir == 1);
  return (segment_[pardir].size() > 0);
//   return (pardir == 2) ? (segment_v_.size() > 0) : (segment_u_.size() > 0);
}

//===========================================================================
bool Param2FunctionInt::hasCriticalValsOrKnots(int pardir) const
//===========================================================================
{
    return hasCriticalVals(pardir);
}

//===========================================================================
bool Param2FunctionInt::canDivide(int pardir)
//===========================================================================
{
    // Flatness should be handled by sortParameterDirections() in intersector.
//     // @@sbr Stop dividing if almost (tol-equal) flat?
//     // Maybe more corresponding to shouldDivide() ... ?
//     bool is_deg = isDegenerate(deg_tol_, pardir);
//     if (is_deg)
// 	return false;

    // This test must check whether the parameter interval of the surf
    // is long enough and probably also some other critera
    // @@@ VSK, for the time being ...
    double ta = startParam(pardir);
    double tb = endParam(pardir);
    double tc = 0.5*(ta+tb);

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
vector<double> Param2FunctionInt::getCriticalVals(int pardir) const
//===========================================================================
{
  // Some values might be more critical than others (according to the
  // priority in segment_). Other might not be really critical. Should
  // all values be returned?
  // What about sorting? According to priority? According to the distance
  // from the middle of this surf? At least it exists more information
  // valuable for sorting here than in the Intersectors.

  // @@@ VSK. For the time being I sort according to priorites, but not
  // accoring to the parameter. I return all values.
  int ki;
  vector<double> vals;
  if (pardir == 1) {
      if (segment_[1].size() == 0)
	  return vals;

      // Sort the segment array.
      vector<pair<double, int> > segment_copy = segment_[1];
      sort(segment_copy.begin(), segment_copy.end(), compare_prio);

      for (ki=0; ki<int(segment_[1].size()); ki++)
	  vals.push_back(segment_copy[ki].first);
  } else {
      if (segment_[0].size() == 0)
	  return vals;

      // Sort the segment array.
      vector<pair<double, int> > segment_copy = segment_[0];
      std::sort(segment_copy.begin(), segment_copy.end(), compare_prio);

      for (ki=0; ki<int(segment_[0].size()); ki++)
	  vals.push_back(segment_copy[ki].first);
  }
  
  return vals;
}

//===========================================================================
vector<double> Param2FunctionInt::getCriticalValsAndKnots(int pardir) const
//===========================================================================
{
  vector<double> vals = getCriticalVals(pardir);
  return vals;
}

//===========================================================================
int Param2FunctionInt::getMeshSize(int dir)
//===========================================================================
{
    THROW("Not yet implemented!");
}

//===========================================================================
double Param2FunctionInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
    THROW("Not yet implemented!");
}

//===========================================================================
std::vector<double>::iterator Param2FunctionInt::getMesh()
//===========================================================================
{
    THROW("Not yet implemented!");
}

//===========================================================================
vector<double> Param2FunctionInt::getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
  // Nothing for a general parametric surf
  vector<double> vals;
  return vals;
}


//===========================================================================
double Param2FunctionInt::startParam(int pardir) const
//===========================================================================
{
    RectDomain cont_dom = surf_->containingDomain();
    return (pardir == 0) ? cont_dom.umin() : cont_dom.vmin();
}


//===========================================================================
double Param2FunctionInt::endParam(int pardir) const
//===========================================================================
{
    RectDomain cont_dom = surf_->containingDomain();
    return (pardir == 0) ? cont_dom.umax() : cont_dom.vmax();
}


//===========================================================================
bool Param2FunctionInt::boundaryPoint(const double* par, double eps) const
//===========================================================================
{
  const Domain& dom = getSurface()->parameterDomain();

  Vector2D par_pt(par[0], par[1]);
  return dom.isInDomain(par_pt, eps);
}


//===========================================================================
void Param2FunctionInt::
subdivide(int pardir, double par,
	  vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
	  vector<shared_ptr<ParamFunctionInt> >& bd_objs)
//===========================================================================
{
  double start = startParam(pardir);
  double end = endParam(pardir);

  int other_pardir = (pardir == 1) ? 0 : 1;
  double other_start = startParam(other_pardir);
  double other_end = endParam(other_pardir);

  shared_ptr<SplineSurface> spline_sf =
      dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);
  ASSERT(spline_sf.get() != 0); // @@sbr Does not cope well with a
				// trimmed sf ...

  // Perform subdivision
  // @@@ VSK, This calls should be made more effective, but can do for
  // the time being
  shared_ptr<ParamSurface> sf1, sf2;
  // @@sbr But what if sf is trimmed ...
  double fuzzy_tol = DEFAULT_PARAMETER_EPSILON; //1e-08;
  try {
      if (pardir == 1) {
	  sf1 = shared_ptr<ParamSurface>(spline_sf->
					 subSurface(other_start, start,
						    other_end, par,
						    fuzzy_tol));
	  sf2 = shared_ptr<ParamSurface>(spline_sf->
					 subSurface(other_start, par,
						    other_end, end,
						    fuzzy_tol));
      } else {
	  sf1 = shared_ptr<ParamSurface>(spline_sf->
					 subSurface(start, other_start,
						    par, other_end,
						    fuzzy_tol));
	  sf2 = shared_ptr<ParamSurface>(spline_sf->
					 subSurface(par, other_start,
						    end, other_end,
						    fuzzy_tol));
      }
  } catch (...) {
      THROW("Failed extracting subsurface!");
  }
  DEBUG_ERROR_IF(sf1.get()==0 || sf2.get()==0, "Error in subdivide");

  // Get the subdivision curve
  // pardir == 1 => par is a v-par.
  int bdidx = (pardir == 0) ? 3 : 1; // @@sbr Is this correct?
  int sub_pardir;
  double sub_tpar;
  shared_ptr<SplineSurface> spline_sf1 =
      dynamic_pointer_cast<SplineSurface, ParamSurface>(sf1);
  shared_ptr<ParamCurve> subsf_cv(extractBdCurve(*spline_sf1, bdidx,
						 sub_pardir, sub_tpar));
  // Make intersection objects
  subdiv_objs.push_back(makeIntFunction(sf1));
  subdiv_objs.push_back(makeIntFunction(sf2));
  bd_objs.push_back(shared_ptr<ParamFunctionInt>
		    (new Spline1FunctionInt(subsf_cv, this))); // @@sbr
							       // Instead
							       // of
							       // Param1F...
}

//===========================================================================
shared_ptr<Param2FunctionInt> 
   Param2FunctionInt::makeIntFunction(shared_ptr<ParamSurface> surf)
//===========================================================================
{
  shared_ptr<Param2FunctionInt> surf_int =
      shared_ptr<Param2FunctionInt>(new Param2FunctionInt(surf, this));
  return surf_int;
}


//===========================================================================
CompositeBox Param2FunctionInt::compositeBox() const
//===========================================================================
{
  return surf_->compositeBox();
}


//===========================================================================
bool Param2FunctionInt::monotone(Point& dir, double tol) const
//===========================================================================
{
    // @@sbr It seems we need more info to decide upon monotonicity.
    // Implement using sample approach!
    // Possibly store pts (and corr params).
    return false;
}


//===========================================================================
void Param2FunctionInt::
getBoundaryObjects(vector<shared_ptr<BoundaryFunctionInt> >& bd_objs)
//===========================================================================
{
    // @@sbr Here we need to do something ...
    THROW("Should be overrun by inherited method!");
}


//===========================================================================
double Param2FunctionInt::isolateDegPar(int dir, int deg_edge, 
					double threshold)
//===========================================================================
{
    // Specify boundary parameters
    int  deg = (deg_edge == 3) ? 1 : deg_edge;
    int bd_idx = 2*dir + deg - 1;
    double bd_par = (deg == 1) ? startParam(dir) : endParam(dir);
    
    if (deg == 0 || deg_tol_ < 0 || !bd_deg_[bd_idx])
	return bd_par;  // No degeneracy

    // Compute the end parameter of the domain corresponding to a 
    // piece of the surface of size fac*epsge from the degenerate edge
    double fac = 10.0;
    int nsample = 5;
    double param[2], delta;
    int ki;

    param[dir] = bd_par;
    param[1-dir] = startParam(1-dir);
    delta = (endParam(1-dir) - param[1-dir])/(double)(nsample-1);

    Point der[2];
    double length;
    double delta_par = 0.0;
    for (ki=0; ki<nsample; ki++, param[1-dir] += delta)
    {
	derivs(param[0], param[1], der[0], der[1]);
	length = der[dir].length();
	delta_par = std::max(delta_par, fac*deg_tol_/length);
    }

    return (deg == 1) ? bd_par + delta_par : bd_par - delta_par;
}


//===========================================================================
void Param2FunctionInt::assureInRange(int pardir, double& t)
//===========================================================================
{
  double tmin = startParam(pardir);
  double tmax = endParam(pardir);
  if (t < tmin)
    t = tmin;
  else if (t > tmax)
    t = tmax;
}


//===========================================================================
int Param2FunctionInt::
knotIntervalFuzzy(int pardir, double& t, double tol) const
//===========================================================================
{
  double tmin = startParam(pardir);
  double tmax = endParam(pardir);
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
double Param2FunctionInt::
nextSegmentVal(int pardir, double par, bool forward) const
//===========================================================================
{
    // No inner knots expected to exist.
    if (forward) {
	return endParam(pardir);
    } else {
	return startParam(pardir);
    }
}


//===========================================================================
bool Param2FunctionInt::isDegenerate(double epsge, int pardir)
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);

    if (deg_tol_ < 0.0 || deg_tol_ != epsge)
    {
	// Check degeneracy
	deg_tol_ = epsge;
    }

    shared_ptr<SplineSurface> spline_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);
    if (spline_sf.get() != 0) {
	// We find diff between max and min elem along lines in input dir.
	int nmb_cvs = (pardir == 0)
	    ? spline_sf->numCoefs_v() : spline_sf->numCoefs_u();
	int nmb_samples = (pardir == 0)
	    ? spline_sf->numCoefs_u() : spline_sf->numCoefs_v();
	double max_diff = -1.0;
	int cv_step = (pardir == 0) ? spline_sf->numCoefs_u() : 1;
	int pt_step = (pardir == 0) ? 1 : spline_sf->numCoefs_u();
	for (int ki = 0; ki < nmb_cvs; ++ki) {
	    double min_val = std::numeric_limits<double>::max();
	    double max_val = std::numeric_limits<double>::lowest();
	    for (int kj = 0; kj < nmb_samples; ++kj) {
		double val = spline_sf->coefs_begin()[ki*cv_step+kj*pt_step];
		if (val < min_val)
		    min_val = val;
		if (val > max_val)
		    max_val = val;
	    }
	    double diff = max_val - min_val;
	    if (diff > max_diff)
		max_diff = diff;
	}
	if (max_diff < deg_tol_)
	    return true;
    }

//     CompositeBox box = compositeBox();
//     Point low = box.low(0.0, 0.0);
//     Point high = box.high(0.0, 0.0);
//     int dim = box.dimension();
//     for (int ki=0; ki<dim; ki++)
// 	if (high[ki] - low[ki] > epsge)
// 	    return false;
   
//     return true;

    return false;
}


//===========================================================================
bool Param2FunctionInt::isDegenerate(double epsge, int pardir, double* par)
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    if (par != 0) { // @@sbr Will be replaced by extraction of iso-cv.
                    // Currently looking at all corr coefs in spec dir.
	MESSAGE("Currently not treating iso-cvs, but doing something.");
    }

    return isDegenerate(epsge, pardir);
}


//===========================================================================
SplineCurve*
Param2FunctionInt::extractBdCurve(const SplineSurface& spline_sf, int bdidx,
				  int& pardir, double& tpar) const
//===========================================================================
{
    double tmin = (bdidx < 2) 
	? spline_sf.startparam_u() : spline_sf.startparam_v();
    double tmax = (bdidx < 2)
	? spline_sf.endparam_u() : spline_sf.endparam_v();
    if (bdidx == 0)
	tpar = spline_sf.startparam_v();
    else if (bdidx == 1)
	tpar = spline_sf.endparam_v();
    else if (bdidx == 2)
	tpar = spline_sf.startparam_u();
    else if (bdidx == 3)
	tpar = spline_sf.endparam_u();
    else
	THROW("Unexpected bd idx.");

    SplineCurve *bd_cv = NULL;
    SplineCurve *cross_cv = NULL;
    double rel_par_res = 1e-12; // = epsge_->getRelParRes();
    spline_sf.getBoundaryInfo(tmin, tmax, bdidx, bd_cv, cross_cv,
			      rel_par_res);

    delete cross_cv;

    pardir = (bdidx < 2) ? 1 : 0; // pardir is dir in which tpar lies.
    return bd_cv;
}


//===========================================================================


} // namespace Go

