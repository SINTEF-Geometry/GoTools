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

#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/RotatedBox.h"


using std::vector;
using std::cerr;
using std::endl;
using std::min;
using std::pair;


namespace Go {


// //==========================================================================
// ParamSurfaceInt::ParamSurfaceInt(shared_ptr<ParamSurface> surf)
//   : ParamGeomInt(0), 
//     surf_(surf), 
//     deg_tol_(-1.0), 
//     domain_(surf_->containingDomain()), 
//     temp_point_array_(6)
// //==========================================================================
// {
//   dim_ = surf_->dimension();
//   nmesh_[0] = nmesh_[1] = -1;
// }


//===========================================================================
ParamSurfaceInt::ParamSurfaceInt(shared_ptr<ParamSurface> surf,
				 ParamGeomInt* parent)
    : ParamGeomInt(parent), surf_(surf), deg_tol_(-1.0), deg_triang_(false),
      domain_(surf_->containingDomain()), lw_set_(false),
      temp_point_array_(6),
      implicit_tol_(-1.0), impl_deg_(3)
//===========================================================================
{
    dim_ = surf_->dimension();
    nmesh_[0] = nmesh_[1] = -1;
    bd_deg_[0] = bd_deg_[1] = bd_deg_[2] = bd_deg_[3] = false;
    deg_domain_ = domain_;
}


//===========================================================================
ParamSurfaceInt::~ParamSurfaceInt()
//===========================================================================
{
}
   

 //===========================================================================
shared_ptr<ParamSurfaceInt> 
ParamSurfaceInt::makeIntObject(shared_ptr<ParamSurface> surf)
//===========================================================================
{
    shared_ptr<ParamSurfaceInt> surf_int =
	shared_ptr<ParamSurfaceInt>(new ParamSurfaceInt(surf, this));
    return surf_int;
}


 //===========================================================================
shared_ptr<ParamCurveInt> 
   ParamSurfaceInt::makeIntCurve(shared_ptr<ParamCurve> crv,  
				 ParamGeomInt *parent)
//===========================================================================
{
    shared_ptr<ParamCurveInt> curve_int =
	shared_ptr<ParamCurveInt>(new ParamCurveInt(crv, parent));
    return curve_int;
}


//===========================================================================
ParamSurfaceInt* ParamSurfaceInt::getParamSurfaceInt()
//===========================================================================
{
    return this;
}


//===========================================================================
shared_ptr<ParamSurface> 
ParamSurfaceInt::getParentParamSurface(RectDomain& domain)
//===========================================================================
{
    domain = domain_;
    if (parent_ && parent_->numParams() == 2) {
	ParamSurfaceInt *parentsf
	    = dynamic_cast<ParamSurfaceInt*>(parent_)->getParamSurfaceInt();
	if (parentsf == 0) {
	    return surf_;
	} else {
	    return parentsf->getParentParamSurface();
	}
    } else {
	return surf_;
    }
}


//===========================================================================
shared_ptr<ParamCurve> 
ParamSurfaceInt::getIsoCurve(double param_start, 
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
ParamSurfaceInt::getConstantParameterCurve(int dir, double par)
//===========================================================================
{
    // Get the current parameter domain
    vector<double> mima = getMima();

    // End-parameter values of the curve
    double param1[2], param2[2];
    int opp_dir = 1-dir;
    param1[dir] = param2[dir] = par;
    param1[opp_dir] = mima[2*opp_dir];
    param2[opp_dir] = mima[2*opp_dir+1];

    // Make constant parameter curve as a curve-on-surface curve
    Point pp1(param1[0], param1[1]), pp2(param2[0], param2[1]);
    shared_ptr<SplineCurve> pcrv =
	shared_ptr<SplineCurve>(new SplineCurve(pp1, param1[opp_dir],
						pp2, param2[opp_dir]));
    shared_ptr<ParamCurve> constcrv = 
	shared_ptr<ParamCurve>(new CurveOnSurface(getParamSurface(),
						  pcrv, true));

    return constcrv;
}

//===========================================================================
shared_ptr<ParamCurve> 
ParamSurfaceInt::getConstantParameterCurve(int dir, double par, double tmin, 
					   double tmax)
//===========================================================================
{
    // Get the current parameter domain
    vector<double> mima = getMima();

    // End-parameter values of the curve
    double param1[2], param2[2];
    int opp_dir = 1-dir;
    param1[dir] = param2[dir] = par;
    param1[opp_dir] = std::max(tmin, mima[2*opp_dir]);
    param2[opp_dir] = std::min(tmax, mima[2*opp_dir+1]);

    // Make constant parameter curve as a curve-on-surface curve
    Point pp1(param1[0], param1[1]), pp2(param2[0], param2[1]);
    shared_ptr<SplineCurve> pcrv =
	shared_ptr<SplineCurve>(new SplineCurve(pp1, param1[opp_dir],
						pp2, param2[opp_dir]));
    shared_ptr<ParamCurve> constcrv = 
	shared_ptr<ParamCurve>(new CurveOnSurface(getParamSurface(),
						  pcrv, true));

    return constcrv;
}

//===========================================================================
void 
ParamSurfaceInt::minimumAlongCurve(int dir, double par, double tmin, double tmax,
				   Point& pnt, double& minpar, Point& minval, 
				   double& mindist)
//===========================================================================
{
    // First describe a constant parameter curve as a curve-on-surface
    shared_ptr<ParamCurve> curve = getConstantParameterCurve(dir, par, tmin, tmax);

    // Iterate to the cloest point
    curve->closestPoint(pnt, tmin, tmax, minpar, minval, mindist);
}


//===========================================================================
shared_ptr<ParamSurface> 
ParamSurfaceInt::getParentParamSurface()
//===========================================================================
{
    if (parent_ && parent_->numParams() == 2) {
	ParamSurfaceInt *parentsf = 
	    dynamic_cast<ParamSurfaceInt*>(parent_)->getParamSurfaceInt();
	if (parentsf == 0) {
	    return surf_;
	} else {
	    return parentsf->getParentParamSurface();
	}
    } else {
	return surf_;
    }
}


//===========================================================================
double ParamSurfaceInt::startParam(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return domain_.lowerLeft()[pardir];
}


//===========================================================================
double ParamSurfaceInt::endParam(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return domain_.upperRight()[pardir];
}


//===========================================================================
bool ParamSurfaceInt::boundaryPoint(const double* par, double eps) const 
//===========================================================================
{
    const Domain& dom = getParamSurface()->parameterDomain();

    Vector2D par_pt(par[0], par[1]);
    return dom.isOnBoundary(par_pt, eps);
}


//===========================================================================
vector<double> ParamSurfaceInt::getMima() const
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
bool ParamSurfaceInt::inCorner(const double *par, double epspar) const
//===========================================================================
{
    if (!(fabs(par[0]-domain_.umin()) < epspar || 
	  fabs(par[0]-domain_.umax()) < epspar))
	return false;
    if (!(fabs(par[1]-domain_.vmin()) < epspar
	  || fabs(par[1]-domain_.vmax()) < epspar))
	return false;
    return true;
}

//===========================================================================
bool ParamSurfaceInt::atDegenerateBd(const double *par, double epsge, double epspar) 
//===========================================================================
{
    // First check if the surface is degenerate
    if (!isDegenerate(epsge))
	return false;

    int degen;
    double deg_par;
    for (int dir=0; dir<2; dir++)
    {
	// Check degeneracy in the current direction
	degen = isDegenerate(epsge, dir);
	if (degen == 1 || degen == 3)
	{
	    deg_par = startParam(dir);
	    if (fabs(deg_par - par[dir]) < epspar)
		return true;
	}
	if (degen == 2 || degen == 3)
	{
	    deg_par = endParam(dir);
	    if (fabs(deg_par - par[dir]) < epspar)
		return true;
	}
    }
    return false;
}

//===========================================================================
void ParamSurfaceInt::getLengthAndWiggle(double *length, double *wiggle)
//===========================================================================
{
    if (lw_set_) {
	length[0] = length_[0];
	length[1] = length_[1];
	wiggle[0] = wiggle_[0];
	wiggle[1] = wiggle_[1];
    } else {
	int meshsize = 5;
	if (nmesh_[0] < 3 || nmesh_[1] < 3)
	    makeMesh(meshsize, meshsize);

	int ki, kj, idx, idx2;
	double currl, currw;
	double d1=0.0, d2=0.0;
	Point pt1(dim_), pt2(dim_);
	Point vec1, vec2;

	// 1. parameter direction
	currl = currw = 0.0;
	for (kj=0, idx=0; kj<nmesh_[1]; kj++) {
	    pt1.setValue(&mesh_[idx]);
	    for (ki=1, idx+=dim_; ki<nmesh_[0]; ki++, idx+=dim_) {
		pt2.setValue(&mesh_[idx]);
		d2 = pt1.dist(pt2);
		currl += d2;
		vec2 = pt2 - pt1;
		if (ki>1 && d1 > 0.0 && d2 > 0.0)
		    currw += vec1.angle(vec2);
		pt1 = pt2;
		vec1 = vec2;
		d1 = d2;
	    }
	}
	currl /= (double)(nmesh_[1]);
	currw /= (double)(nmesh_[1]);
	length_[0] = length[0] = currl;
	wiggle_[0] = wiggle[0] = currw;
      
	// 2. parameter direction
	currl = currw = 0.0;
	for (ki=0, idx=0; ki<nmesh_[0]; ki++, idx+=dim_) {
	    pt1.setValue(&mesh_[idx]);
	    for (kj=1, idx2=idx+dim_*nmesh_[0]; kj<nmesh_[1]; 
		 kj++, idx2+=dim_*nmesh_[0]) {
		pt2.setValue(&mesh_[idx2]);
		currl += pt1.dist(pt2);
		vec2 = pt2 - pt1;
		if (kj>1)
		    currw += vec1.angle(vec2);
		pt1 = pt2;
		vec1 = vec2;
	    }
	}
	currl /= (double)(nmesh_[0]);
	currw /= (double)(nmesh_[0]);
	length_[1] = length[1] = currl;
	wiggle_[1] = wiggle[1] = currw;
	lw_set_ = true;
    }
}


//===========================================================================
bool ParamSurfaceInt::hasInnerKnots(int pardir) const
//===========================================================================
{
  return false;
}


//===========================================================================
bool ParamSurfaceInt::hasCriticalVals(int pardir) const 
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return (segment_[pardir].size() > 0 || checkPeriodicity(pardir) >= 1);
}


//===========================================================================
bool ParamSurfaceInt::hasCriticalValsOrKnots(int pardir) const 
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return (segment_[pardir].size() > 0);
}


//===========================================================================
bool ParamSurfaceInt::canDivideTinyTriang(int pardir)
//===========================================================================
{
    // This test must check whether the parameter interval of the
    // surface in this direction is long enough and probably also some
    // other critera
    // @@@ VSK, for the time being ...
    // @@@ VSK, This function could be moved one step up in the class
    // hierarchy if it continues to be similar the the curve function
    // Check for degeneracy and degenerate triangle
    if (deg_tol_ > 0.0) {
	if (deg_triang_) {
	    if (bd_deg_[2*pardir] || bd_deg_[2*pardir+1]) {
		return false;
	    }
	} else {
	    if (bd_deg_[2*(1-pardir)] || bd_deg_[2*(1-pardir)+1])
		return false;
	}
    }
    return true;
}

//===========================================================================
bool ParamSurfaceInt::canDivide(int pardir)
//===========================================================================
{
    // This test must check whether the parameter interval of the
    // surface in this direction is long enough and probably also some
    // other critera
    // @@@ VSK, for the time being ...
    // @@@ VSK, This function could be moved one step up in the class
    // hierarchy if it continues to be similar the the curve function


    // Check for size of parameter domain (we need some parameter
    // tolerance here)
    double ta = startParam(pardir);
    double tb = endParam(pardir);
    double tc = 0.5*(ta+tb);
    
    //return (tc != ta && tc != tb);
    double rel_par_res = 1e-10; // = epsge_->getRelParRes();
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
    // Comparison function used for sorting in getCriticalVals()

    if (seg1.second < seg2.second)
	return false;
    else
	return true;
}


//===========================================================================
vector<double> ParamSurfaceInt::getCriticalVals(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    // Some values might be more critical than others (according to
    // the priority in segment_). Other might not be really
    // critical. Should all values be returned?
    // What about sorting? According to priority? According to the
    // distance from the middle of this curve? At least it exists more
    // information valuable for sorting here than in the Intersectors.
    // @@@ VSK. For the time being I sort according to priorites, but
    // not accoring to the parameter. I return all values.

    vector<double> vals;
    int ki;
    if (segment_[pardir].size() > 0) {
	// Sort the segment array.
	vector<pair<double, int> > segment_copy = segment_[pardir];
	std::sort(segment_copy.begin(), segment_copy.end(), compare_prio);

	for (ki=0; ki<int(segment_copy.size()); ki++) {
	    vals.push_back(segment_copy[ki].first);
	}
    }
    // The seem of a periodic surface is also a critical parameter
    int per = checkPeriodicity(pardir);
    if (per >= 1) {
	vals.insert(vals.begin(), startParam(pardir));
    }

    return vals;
}


//===========================================================================
vector<double> ParamSurfaceInt::getCriticalValsAndKnots(int pardir) const
//===========================================================================
{
    vector<double> vals = getCriticalVals(pardir);
    return vals;
}


 //===========================================================================
bool ParamSurfaceInt::isSpline()
//===========================================================================
{
    return false;
}


//===========================================================================
double
ParamSurfaceInt::getOptimizedConeAngle(Point& axis1, Point& axis2)
//===========================================================================
{
    return M_PI;
}


//===========================================================================
void ParamSurfaceInt::
knotIntervalFuzzy(double& u, double&v, double utol, double vtol) const
//===========================================================================
{
    const double& umin = domain_.umin();
    const double& umax = domain_.umax();    
    const double& vmin = domain_.vmin();
    const double& vmax = domain_.vmax();
    u = (u < umin + utol) ? umin :
	(u > umax - utol) ? umax : u;
    v = (v < vmin + vtol) ? vmin :
	(v > vmax - vtol) ? vmax : v;
}


//===========================================================================
double ParamSurfaceInt::
nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    if (forward) {
	return endParam(dir);
    } else {
	return startParam(dir);
    }
}


//===========================================================================
vector<double> ParamSurfaceInt::getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
    // Nothing for a general parametric surface
    vector<double> vals;
    return vals;
}


//===========================================================================
int ParamSurfaceInt::getMeshSize(int dir)
//===========================================================================
{
    int meshsize = 5;
    if (dir != 0 && dir != 1) {
	return 1;
    } else {
	if (int(mesh_.size()) != nmesh_[0]*nmesh_[1]) {
	    makeMesh(meshsize, meshsize);
	}
	return nmesh_[dir];
    }
}


//===========================================================================
CompositeBox ParamSurfaceInt::compositeBox() const 
//===========================================================================
{
  return surf_->compositeBox();
}


//===========================================================================
DirectionCone ParamSurfaceInt::directionCone() const
//===========================================================================
{
    if (cone_.greaterThanPi() < 0)
	cone_ = surf_->normalCone();
    return cone_;
}


//===========================================================================
int ParamSurfaceInt::checkPeriodicity(int pardir) const
//===========================================================================
{
    // Don't know anything about periodicity
    return -1;  // Open, non-periodic curve
}


//===========================================================================
void ParamSurfaceInt::
getBoundaryObjects(std::vector<shared_ptr<BoundaryGeomInt> >& bd_objs)
//===========================================================================
{
    if (boundary_obj_.size() != 0) {
	// Assumes that the boundary curves are fetched correctly
	bd_objs.insert(bd_objs.begin(), boundary_obj_.begin(), 
		       boundary_obj_.end());
    } else {
	boundary_obj_.clear();

	// This is a very general implementation that should be
	// specialized in the subclasses
	int ki, kj;
	std::vector<CurveLoop> bounds = surf_->allBoundaryLoops();
	for (ki=0; ki<int(bounds.size()); ki++) {
	    for (kj=0; kj<int(bounds[ki].size()); kj++) {
		shared_ptr<ParamCurve> curr_cv = bounds[ki][kj];
		shared_ptr<ParamCurveInt> curr_cv_int = 
		    shared_ptr<ParamCurveInt>
		    (new ParamCurveInt(curr_cv, this));
		shared_ptr<BoundaryGeomInt> bd = 
		    shared_ptr<BoundaryGeomInt>
		    (new BoundaryGeomInt(curr_cv_int, -1, 0.0));
		// @@@ VSK. How can we check if this curve is a
		// constant parameter curve? How can we specify a
		// trimming curve
		boundary_obj_.push_back(bd);
		bd_objs.push_back(bd);
	    }
	}
    }
}


//===========================================================================
vector<shared_ptr<ParamSurfaceInt> >
ParamSurfaceInt::subSurfaces(double from_upar, double from_vpar,
			     double to_upar, double to_vpar,
			     double fuzzy) 
//===========================================================================
{
    shared_ptr<ParamSurface> srf = getParamSurface();
    vector<shared_ptr<ParamSurface> > sub_sfs = 
	srf->subSurfaces(from_upar, from_vpar, to_upar, to_vpar, fuzzy);
    
    vector<shared_ptr<ParamSurfaceInt> > sub_sfs_int;
    for (size_t ki=0; ki<sub_sfs.size(); ki++)
	sub_sfs_int.push_back(makeIntObject(sub_sfs[ki]));
    
    return sub_sfs_int;
}


//===========================================================================
void ParamSurfaceInt::
subdivide(int pardir, double par, 
	  vector<shared_ptr<ParamGeomInt> >& subdiv_objs,
	  vector<shared_ptr<ParamGeomInt> >& bd_objs)
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    double ta1 = domain_.umin();
    double ta2 = domain_.vmin();
    double tb1 = domain_.umax();
    double tb2 = domain_.vmax();

    vector<shared_ptr<ParamSurface> > sub1, sub2;
    shared_ptr<ParamSurface> srf;
//     srf = getParentParamSurface();
    srf = getParamSurface();
    if (pardir == 0) {
	double p_interval = tb1 - ta1;
	int per = -1; // checkPeriodicity(pardir)
	if (per >= 1) {
	    while (par >= tb1) {
		par -= p_interval;
	    }
	    sub1 = srf->subSurfaces(par, ta2, par+p_interval, tb2);
	} else {
	    sub1 = srf->subSurfaces(ta1, ta2, par, tb2);
	    sub2 = srf->subSurfaces(par, ta2, tb1, tb2);
	}
    } else {
 	double p_interval = tb2 - ta2;
	int per = -1; // checkPeriodicity(pardir)
	if (per >= 1) {
	    while (par >= tb2) {
		par -= p_interval;
	    }
	    sub1 = srf->subSurfaces(ta1, par, tb1, par+p_interval);
	} else {
	    sub1 = srf->subSurfaces(ta1, ta2, tb1, par);
	    sub2 = srf->subSurfaces(ta1, par, tb1, tb2);
	}
    }
    for (size_t ki = 0; ki < sub1.size(); ki++)
	subdiv_objs.push_back(makeIntObject(sub1[ki]));
    for (size_t ki = 0; ki < sub2.size(); ki++)
	subdiv_objs.push_back(makeIntObject(sub2[ki]));

    if (getDegTriang()) {
	for (size_t ki = 0; ki < subdiv_objs.size(); ki++) {
	    ParamSurfaceInt *tmp = (ParamSurfaceInt*)(subdiv_objs[ki].get());
	    tmp->setDegTriang();
	}
    }

    vector<shared_ptr<ParamCurve> > div_crvs;
    bool u_dir = (pardir == 1);
    for (size_t ki = 0; ki < sub1.size(); ki++) {
	div_crvs.clear();
	div_crvs  = surf_->constParamCurves(par, u_dir);
	for (size_t kj = 0; kj < div_crvs.size(); kj++) {
	    shared_ptr<ParamCurveInt> curve_int =
		makeIntCurve(div_crvs[kj], this);
	    bd_objs.push_back(curve_int);
	}
    }
}


//===========================================================================
vector<double>::iterator ParamSurfaceInt::getMesh()
//===========================================================================
{
    int meshsize = 5;
    if (int(mesh_.size()) != nmesh_[0]*nmesh_[1])
	makeMesh(meshsize, meshsize);
    return mesh_.begin();
}


//===========================================================================
double ParamSurfaceInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
    if (dir != 0 && dir != 1) {
	return 0.0;  // Arbitrary
    } else {
	if (idx < 0 || idx >= nmesh_[dir]) {
	    return 0.5*(domain_.lowerLeft()[dir] + domain_.upperRight()[dir]);
	} else {
	    return ((double)(nmesh_[dir]-idx)*domain_.lowerLeft()[dir] + 
		    (double)(idx)*domain_.upperRight()[dir])/
		(double)(nmesh_[dir]);
	}
    }
}


//===========================================================================
void ParamSurfaceInt::makeMesh(int size1, int size2) const
//===========================================================================
{
    mesh_.clear();
    mesh_.reserve(size1*size2*dim_);
    nmesh_[0] = size1;
    nmesh_[1] = size2;
    double ta1 = domain_.umin();
    double ta2 = domain_.vmin();
    double tb1 = domain_.umax();
    double tb2 = domain_.vmax();
    double tint1 = (tb1 - ta1)/(double)(size1-1);
    double tint2 = (tb2 - ta2)/(double)(size2-1);
    int ki, kj;
    double par1, par2;
    Point pt;
    for (kj = 0, par2 = ta2; kj < size2; kj++, par2 += tint2) {
	for (ki = 0, par1 = ta1; ki < size1; ki++, par1 += tint1) {
	    surf_->point(pt, par1, par2);
	    mesh_.insert(mesh_.end(), pt.begin(), pt.end());
	}
    }
      
}


//===========================================================================
bool ParamSurfaceInt::isSimple()
//===========================================================================
{
    // Estimates if the current surface is simple enough for a
    // singularity iteration. Checks the span of the normal cone and
    // the size of the surface

    if (cone_.greaterThanPi() < 0)
      {
	try {
	  cone_ = surf_->normalCone();
	}
	catch (...)
	  {
	    return false;  // Probably degenerate somehow
	  }
      }
    if (cone_.angle() < M_PI/3.0)
	return true;
    else
	return false;
}


 //===========================================================================
bool ParamSurfaceInt::isDegenerate(double epsge)
//===========================================================================
{
    if (deg_tol_ < 0.0 || deg_tol_ != epsge) {
	// Check degeneracy
	deg_tol_ = epsge;
	bool is_deg;
	is_deg = surf_->isDegenerate(bd_deg_[2], bd_deg_[1], 
				   bd_deg_[3], bd_deg_[0], epsge);
	computeDegDomain(epsge);
	return is_deg;
    } else {
	return (bd_deg_[0] || bd_deg_[1] || bd_deg_[2] || bd_deg_[3]);
    }
}


 //===========================================================================
const RectDomain&  ParamSurfaceInt::getDegDomain(double epsge)
//===========================================================================
{
    if (deg_tol_ < 0.0 || deg_tol_ != epsge) 
    {
	isDegenerate(epsge);
    }
    return deg_domain_;
}

 //===========================================================================
void ParamSurfaceInt::computeDegDomain(double aepsge)
//===========================================================================
{
    if (deg_tol_ < 0.0 || deg_tol_ != aepsge) 
    {
	isDegenerate(aepsge);
    }
    else
    {
	double deg_fac = 0.6;
	double deg_par;
	Point ll(domain_.umin(), domain_.vmin());
	Point ur(domain_.umax(), domain_.vmax());
	if (bd_deg_[0])
	{
	    deg_par = isolateDegPar(0, 1, 0.0, &deg_fac);
	    ll.setValue(deg_par, ll[1]);
	}
	if (bd_deg_[1])
	{
	    deg_par = isolateDegPar(0, 2, 0.0, &deg_fac);
	    ur.setValue(deg_par, ur[1]);
	}
	if (bd_deg_[2])
	{
	    deg_par = isolateDegPar(1, 1, 0.0, &deg_fac);
	    ll.setValue(ll[0],deg_par);
	}
	if (bd_deg_[3])
	{
	    deg_par = isolateDegPar(1, 2, 0.0, &deg_fac);
	    ur.setValue(ur[0], deg_par);
	}

	Vector2D tmp1(ll[0],ll[1]), tmp2(ur[0],ur[1]);
	RectDomain tmp(tmp1, tmp2);
	deg_domain_ = tmp;
    }
}

//===========================================================================
int ParamSurfaceInt::isDegenerate(double epsge, int dir)
//===========================================================================
{
    ASSERT(dir==0 || dir==1);

    int degen = 0;
    if (deg_tol_ < 0.0 || deg_tol_ != epsge) {
	// Check degeneracy
	bool is_deg;
	is_deg = isDegenerate(epsge);
    }
    if (bd_deg_[2*dir])
	degen += 1;
    if (bd_deg_[2*dir+1])
	degen += 2;

    return degen;
}


//===========================================================================
bool ParamSurfaceInt::isDegenerate(double epsge, int dir, double *par)
//===========================================================================
{
    // NB! Assumes that a surface is not degenerate in the inner

    ASSERT(dir==0 || dir==1);

    bool degen = false;
    if (deg_tol_ < 0.0 || deg_tol_ != epsge)
    {
	// Check degeneracy
	bool is_deg;
	is_deg = isDegenerate(epsge);
    }

    if (dir == 0)
    {
	// Looks for degeneracy in the first parameter direction
	// Get parameter interval in the second parameter direction
	double ta = startParam(1);
	double tb = endParam(1);
	double pdel = 0.001*(tb - ta);
	pdel = std::max(pdel, 1.0e-10);
	if (par[1]-ta < pdel && bd_deg_[2])
	    degen = true;
	else if (tb-par[1] && bd_deg_[3])
	    degen = true;
    }
    else
    {
	// Looks for degeneracy in the second parameter direction
	// Get parameter interval in the first parameter direction
	double ta = startParam(0);
	double tb = endParam(0);
	double pdel = 0.001*(tb - ta);
	pdel = std::max(pdel, 1.0e-10);
	if (par[0]-ta < pdel && bd_deg_[0])
	    degen = true;
	else if (tb-par[0] && bd_deg_[1])
	    degen = true;
    }

    return degen;
}


//===========================================================================
double ParamSurfaceInt::isolateDegPar(int dir, int deg_edge, 
				      double threshold, double *deg_factor)
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
    double fac = (deg_factor == NULL) ? 100.0 : (*deg_factor);
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
double ParamSurfaceInt::getParOffBd(int dir, bool atstart, double tol) const
//===========================================================================
{
    // Estimate the parameter value of a surface a specified distance
    // from a given edge
    // Specify boundary parameters
    double bd_par = (atstart) ? startParam(dir) : endParam(dir);
    
    // Compute the end parameter of the domain corresponding to a 
    // piece of the surface of size tol from the specified edge
    int nsample = 10;
    double param[2], delta;
    int ki;

    param[dir] = bd_par;
    param[1-dir] = startParam(1-dir);
    delta = (endParam(1-dir) - param[1-dir])/(double)(nsample-1);

    Point der[2];
    double length;
    double delta_par = endParam(dir) - startParam(dir);
    for (ki=0; ki<nsample; ki++, param[1-dir] += delta)
    {
	derivs(param[0], param[1], der[0], der[1]);
	length = der[dir].length();
	delta_par = std::min(delta_par, tol/length);
    }

    return (atstart) ? bd_par + delta_par : bd_par - delta_par;
}


//===========================================================================
RotatedBox ParamSurfaceInt::getRotatedBox(std::vector<Point>& axis) const
//===========================================================================
{
    int size = 100;
    makeMesh(size, size);
    RotatedBox box(&mesh_[0], dimension(), nmesh_[0], nmesh_[1], &axis[0]);
    return box;
}


//===========================================================================
void ParamSurfaceInt::axisFromCorners(Point& axis1, Point& axis2) const
//===========================================================================
{
    double cpar[8];
    cpar[0] = cpar[6] = startParam(0);
    cpar[2] = cpar[4] = endParam(0);
    cpar[1] = cpar[3] = startParam(1);
    cpar[5] = cpar[7] = endParam(1);

    // Evaluate the surface in the corner
    vector<Point> corner(4);
    int ki;
    for (ki=0; ki<4; ki++)
	point(corner[ki], cpar+2*ki);

//    axis1 = 0.5*((corner[1] - corner[0]) + (corner[2] - corner[3]));
//    axis2 = 0.5*((corner[3] - corner[0]) + (corner[2] - corner[1]));
    axis1 = corner[2] - corner[0];
    axis2 = corner[3] - corner[1];
 }


//===========================================================================
void ParamSurfaceInt::splitAtG0(double angtol,
				vector<shared_ptr<ParamSurfaceInt> >& subG1)
//===========================================================================
{
    // Split surface at G1 discontinuities.  Currently only relevant
    // for spline surfaces
    return;
}


//===========================================================================
void ParamSurfaceInt::getSingularity(double eps, double sing_par[], 
				     Point& sing_pt, double& sing_val, 
				     double *seed)
//===========================================================================
{
    // Iterate for a surface singularity
    surf_->singularity(sing_par[0], sing_par[1], sing_pt, sing_val,
		       eps, &domain_, seed);
}


//===========================================================================
bool ParamSurfaceInt::canSelfIntersect(double epsge) const
//===========================================================================
{
    // Check if the current surface may selfintersect. Uses normal and
    // tangent cones.

    // First fetch normal cone
    if (cone_.greaterThanPi() < 0)
	cone_ = directionCone();

    if (cone_.greaterThanPi())
	return true;

    // Check tangent cones
    DirectionCone tan_u = surf_->tangentCone(true);
    double dot_nu = cone_.centre()*tan_u.centre();
    double angle_u = tan_u.angle();
    if (fabs(dot_nu) >= cos(angle_u/2.0))
	return true;

    DirectionCone tan_v = surf_->tangentCone(false);
    double dot_nv = cone_.centre()*tan_v.centre();
    double angle_v = tan_v.angle();
    if (fabs(dot_nv) >= cos(angle_v/2.0))
	return true;

    return false;  // No possibility for selfintersection
}


// //==========================================================================
// bool ParamSurfaceInt::markDegenerate(double epsge, double eps_threshold)
// //==========================================================================
// {
//     // @@@ VSK. Must first check if a boundary is degenerate. For
//     // spline surfaces this is implemented in the geometry
//     // library. What about other parametric surfaces, trimmed surfaces
//     // and ...? I think it need to be a possible request for all
//     // parametric surfaces. What does it mean for a trimmed surface to
//     // be degenerate? The trimming curve is degenerate or the
//     // underlying surface is degenerate? I don't think a degenerate
//     // trimming curve is serious, and it will not necessarily occur
//     // even if the underlying surface is degenerate.
//     bool isdegen = isDegenerate(epsge);

//     if (!isdegen)
// 	return false;  // No degenerate edges
// }


//===========================================================================


} // namespace Go
