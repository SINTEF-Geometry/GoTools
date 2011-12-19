//===========================================================================
//                                                                           
// File: Par1FuncIntersector.C                                               
//                                                                           
// Created: Tue Sep 21 14:06:12 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par1FuncIntersector.C,v 1.15 2006-11-03 14:15:12 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Par1FuncIntersector.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/Par0FuncIntersector.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"


using std::vector;
using std::setprecision;


namespace Go {


//===========================================================================
Par1FuncIntersector::Par1FuncIntersector(shared_ptr<ParamFunctionInt> func,
					 shared_ptr<ParamFunctionInt> C,
					 shared_ptr<GeoTol> epsge,
					 Intersector* prev,
					 int eliminated_parameter,
					 double eliminated_value)
    : IntersectorFuncConst(func, C, epsge, prev,
			   eliminated_parameter, eliminated_value)
//===========================================================================
{
  
}


//===========================================================================
Par1FuncIntersector::~Par1FuncIntersector()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
Par1FuncIntersector::lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
					   shared_ptr<ParamFunctionInt> obj2,
					   Intersector* prev,
					   int eliminated_parameter,
					   double eliminated_value)
//===========================================================================
{
    // Both objects should be of the type ParamPointInt.
    shared_ptr<Intersector> curr_inter; 
    curr_inter = shared_ptr<Intersector>(new Par0FuncIntersector(obj1, obj2, epsge_, prev,
								 eliminated_parameter,
								 eliminated_value));
    return curr_inter;
}


//===========================================================================
int Par1FuncIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between a point and a function can occur only if the
  // function is degenerate. So ...

  return 0;
}


//===========================================================================
//
// Purpose : Both curve & pt lie within the same epsion ball. It is not a
//           simple case, and there is intersection.
//           This is an uwanted situation. We have to make a result that
//           is as consistent as possible since we cannot subdivide anymore.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004.
//===========================================================================
void Par1FuncIntersector::microCase()
//===========================================================================
{
//     GO_ERROR("Not yet implemented!", InputError());

    // Fetch all intersection points belonging to the two curves
    // The intersection points should be sorted according to the parameter
    // of one of the curves
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getSortedIntersections(int_pts);

    int ki;
    if (int_pts.size() == 0) {
	// No intersection point exist. Construct one internal to the two
	// curves
	// It is probably good enough to take the midpoint since the 
	// curve is so small
	double mid_par1 = 0.5*(func_int_->startParam(0) +
			       func_int_->endParam(0));
	double C_par = C_->startParam(0); // Well, that would be 0.0.
      
	int_results_->addIntersectionPoint(func_int_, C_, getTolerance(),
					   &mid_par1, &C_par);
    } else if (int_pts.size() == 1) {
	// One intersection point. Leave it like this.
	;
    } else {
	// More than one intersection point. Connect.
	for (ki=1; ki<int(int_pts.size()); ki++)
	    int_pts[ki-1]->connectTo(int_pts[ki], MICRO_PAR1FUNC);
    }
}


//===========================================================================
//
// Purpose : Given two parametric curve in a simple case situation, iterate
//           to the intersection point, if any. 
//
// Written by : Sverre Briseid, SINTEF, Sep 2004
//
//===========================================================================
int Par1FuncIntersector::updateIntersections()
//===========================================================================
{
  // Fetch already existing intersection points
  vector<shared_ptr<IntersectionPoint> > int_pts;
  int_results_->getIntersectionPoints(int_pts);

  if (int_pts.size() > 0)
    {
      // At most one intersection point is expected. One or more points
      // exist already. Leave it like this.
	// @@sbr Do we expect object to be monotone?
      return 0;
    }

  // Iterate to an intersection point between the two curves
  double par, dist, tmin, tmax;
  doIterate(par, dist, tmin, tmax);

  if ((par != tmin) && (par != tmax) && (dist <= epsge_->getEpsge()))
    {
      // An intersection point is found. Represent it in the data
      // structure
      // @@@ vsk, In sisl there is a check if the found intersection
      // point is very close to the endpoints of the curve. In that case
      // it is dismissed since the intersections at the endpoints are found
      // already. Here we have already checked if there exist any
      // intersection points. Thus, we know that there does not exist any
      // intersection point at the endpoint.
      // @@sbr But for tangential intersections an endpoint may not be discovered ...
      // As the int_pts.size() == 0 and an intersection point nevertheless exist
      // this should be such a (or similar) situation.
      // @@sbr Another approach could be to avoid subdividing if the obj (2d)
      // is monotone in that direction.
      int_results_->addIntersectionPoint(func_int_, 
					 C_,
					 getTolerance(),
					 &par,
					 NULL);
      return 1; // A new point is found
    }

  //#endif // 0

  return 0;
}


//===========================================================================
//
// Purpose : Given a parametric curve and no simple case situation, 
//           subdivide the curve to produce more simple sub problems.
//           Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004
//
//===========================================================================
int Par1FuncIntersector::doSubdivide()
//===========================================================================
{
//     MESSAGE("Under construction!");

//     int perm[2];  // Two curves means two parameter directions that need sorting

    // Intersection objects belonging to the next recursion level
    vector<shared_ptr<ParamFunctionInt> > sub_functions;

    // Sort the curves according to importance of subdivision according to
    // properties of the curves and already computed intersections
    double deg_tol = epsge_->getEpsge();
    func_int_->isDegenerate(deg_tol, 0, NULL); // @@sbr Just to set deg_tol ...
    bool can_subdiv = func_int_->canDivide(0);
    if (!can_subdiv)
	return 0;   // Not possible to subdivide any of the curves

    // For each parameter direction in prioritized order, fetch an
    // appropriate subdivision parameter, and perform subdivision
    double subdiv_par;
    vector<shared_ptr<ParamFunctionInt> > func_sub;
    vector<shared_ptr<ParamFunctionInt> > subdivpt;  // @@@ vsk, A parametric point?

    bool found = (getSubdivisionParameter(0, subdiv_par) != 0);
    if (!found) {
	// No parameter value is found. Move this parameter direction
	// to the end of the permutation array, and decrease the
	// number of parameter directions where it is possible to 
	// subdivide.
	return 0;
    }

    // Subdivide the current curve
    try {
	func_int_->subdivide(0, subdiv_par, func_sub, subdivpt);
    } catch (...) {
	func_sub.clear();
	subdivpt.clear();
    }
    if (subdivpt.size() < 1 || func_sub.size() == 0)
	return 0;  // No new objects 

    // Intersect the subdivision points with the other object
    shared_ptr<Intersector> subdiv_intersector = 
	lowerOrderIntersector(subdivpt[0], C_,
			      this, 0, subdiv_par);

    // Is it here relevant to fetch existing intersection points
    // and/or insert points into intersection curves before computing
    // intersection points with the subdivision points. Normally, this
    // is mostly of interest when surfaces are involved, but intersection
    // intervals might exist?
    // These computations do anyway involve the intersection pool, but
    // they might be trigged somehow. The parameter direction and value
    // are required information

    subdiv_intersector->compute(false);

    // Check quality of intersection points
    // If the quality is not sufficient, find a new subdivision parameter
    // and repeat the process.
    // Otherwise
    int_results_->includeReducedInts(subdiv_intersector->getIntPool());

    sub_functions.insert(sub_functions.end(), func_sub.begin(), func_sub.end());

    for (size_t ki = 0; ki < sub_functions.size(); ++ki) {
	shared_ptr<Intersector> intersector = 
	    shared_ptr<Intersector>(new Par1FuncIntersector(sub_functions[ki],
							    C_, epsge_, this));
// 	intersector->getIntPool()->setPoolInfo(int_results_);
	sub_intersectors_.push_back(intersector);
    }

    return 1;
}


//===========================================================================
void Par1FuncIntersector::doIterate(double& par1, double& dist,
				    double& tmin, double& tmax)
//===========================================================================
{
  // Iterate to an intersection point. We know that we have got one curve
  // and one point. Fetch the data instances
  Point pt;
  //double tpar = 0.0;
//   double *par=0;

  // We need create a good initial guess parameter.
  // As the curve is not guaranteed to be monotone we otherwise risk finding
  // a local minimum which is not global.

  Param1FunctionInt* par1_func_int = func_int_->getParam1FunctionInt();
  ASSERT(par1_func_int != 0);
  shared_ptr<ParamCurve> curve = par1_func_int->getParentParamCurve(tmin, tmax);
//   par = &clo_par;

  Point c_pt;
  C_->point(c_pt, NULL);

  // Perform iteration
  Point clo_pt;
  double seed[1];
  getIterationSeed(seed);
  curve->closestPoint(c_pt, tmin, tmax, par1, clo_pt, dist, seed);

  // @@sbr We check if we iterated to a pt too close to tmax-param. If so
  // move to tmax-param and check dist.
  double ptol = 1000.0*epsge_->getRelParRes();
  if ((par1 != tmin) && (fabs(par1 - tmin) < ptol)) {
      par1 = tmin;
      clo_pt = curve->point(par1);
      dist = c_pt.dist(clo_pt);
  } else if ((par1 != tmax) && (fabs(tmax - par1) < ptol)) {
      par1 = tmax;
      clo_pt = curve->point(par1);
      dist = c_pt.dist(clo_pt);
  }
}


//===========================================================================
int Par1FuncIntersector::getSubdivisionParameter(int dir, double& par)
//===========================================================================
{
    // Purpose : Return (what seems to be) the most suitable
    // subdivision param.

  int ki, kj, sgn;
//   int nmbdir1 = func_int_->numParams();
//   int nmbdir2 = C_->numParams();
  double ptol = 100.0*epsge_->getRelParRes();
  double gtol = 100.0*epsge_->getEpsge(); //100.0*epsge_->getEpsge();

  // Set pointer to the intersection object corresponding to the parameter
  // direction
  ParamFunctionInt *obj = func_int_.get();
  int pdir = 0;
  double ta = obj->startParam(pdir);
  double tb = obj->endParam(pdir);
  double fac = 0.1;
  double del = fac*(tb - ta);

  // Get critical parameters
  // @@@ Critical parameters of different priority? Sorting?
  vector<double> critical_pars = obj->getCriticalVals(pdir);

  int size = (int)critical_pars.size();
  int is_critical;

  if (size > 0)
    {
      // Check suitability of the critical parameters
      for (ki=0; ki<size; ki++)
	{
//	  is_critical = int_results_->inInfluenceArea(pdir, critical_pars[ki]);
//	  if (is_critical == 0 || is_critical == 2)
	    par = critical_pars[ki];
	    if (par >= ta + del && par <= tb - del)
	    {
		return 1;
	    }
	}
    }

  // Look for a suitable knot in which to subdivide. First fetch
  // the knots sorted according to the distance to the mid-parameter
  vector<double> knot_vals = obj->getInnerKnotVals(pdir, true);
  size = (int)knot_vals.size();
   if (size > 0)
    {
      // Check suitability of the knots
      for (ki=0; ki<size; ki++)
	{
	    par = knot_vals[ki];
//	  is_critical = int_results_->inInfluenceArea(pdir, knot_vals[ki]);
//	  if (is_critical == 0 || is_critical == 2)
	    
	    if (par >= ta + del && par <= tb - del)
	    {
		return 2;
	    }
	}
    }
 
   // Check for inner intersection points in which to subdivide
   vector<double> inner_ints = int_results_->getSortedInnerInts(pdir);
   size = (int)inner_ints.size();
  if (size > 0)
    {
      // Check suitability of the intersection points
      for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
	   ki+=sgn*kj, kj++, sgn*=-1)
	{
	    par = inner_ints[ki];
// 	  is_critical = int_results_->inInfluenceArea(pdir, inner_ints[ki]);
// 	  if (is_critical == 0 || is_critical == 2)
	    if (par >= ta + del && par <= tb - del)
	    {
		return 3;
	    }
	}
    }

  // If a singularity exist we use that as a subdivision point.

  // Iterate for a closest point in which to subdivide
  if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_ = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }
  
  // Check if a closest point exist already
  double dist = 2*gtol; // At least larger than gtol.
  if (singularity_info_->hasPoint())
  {
      dist = 0.0;
      par = singularity_info_->getParam(dir);
  }
  else
  {
      // Iterate for a closest point
      double tmin, tmax;
      doIterate(par, dist, tmin, tmax);
      if (dist < gtol) {

	  // We also check if the found parameter lies close to the iso-cv of an
	  // existing int_pt in prev_intersector_.
	  shared_ptr<IntersectionPool> prev_pool = prev_intersector_->getIntPool();
	  vector<shared_ptr<IntersectionPoint> > int_pts;
	  prev_pool->getSortedIntersections(int_pts);
	  Point c_pt;
	  C_->point(c_pt, NULL);
	  double max_knot_dist = epsge_->getEpsge(); //epsge_->getRelParRes(); @@sbr Too large value?
	  for (size_t ki = 0; ki < int_pts.size(); ++ki) {
	      const double* prev_par = int_pts[ki]->getPar1();
	      if (fabs(prev_par[dir] - par) < max_knot_dist) {
		  // We check if par[dir] is close enough, i.e. within epsge_.
		  Point func_pt;
		  func_int_->point(func_pt, &prev_par[dir]);
		  if (c_pt.dist(func_pt) < epsge_->getEpsge()) {
		      par = prev_par[dir];
		  }
	      }
	  }

	  singularity_info_->setSingularPoint(&par, 1);
	  std::cout << setprecision(15) << "Singularity: " << par << std::endl;
      }
  }

  if (dist < gtol && par > ta+ptol && par < tb-ptol)
  {
      // Check the parameter value
      vector<shared_ptr<IntersectionPoint> > int_pts;
      is_critical = int_results_->inInfluenceArea(dir, par);
      if (is_critical == 0 || is_critical == 2)
	  return 4;     
  }


//   // Subdivide in the middle
//   par = 0.5*(ta+tb);
//   return true;
  
    // Subdivide at a suitable parameter
  double divpar = 0.5*(ta+tb);
  del = 0.1*(tb-ta);
  double tint = del;
  for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=del, sgn*=-1)
    {
      is_critical = int_results_->inInfluenceArea(dir, divpar);
      if (is_critical == 0 || is_critical == 2)
	{
	  par = divpar;
	  return 5;
	}
    }

  return 0;  // No subdivision parameter found
}


void Par1FuncIntersector::getIterationSeed(double seed[])
{
    // We do it the easy way, computing closest point to the mesh.
//     int mesh_size = func_int_->getMeshSize(0);
    int mesh_size = func_int_->getMeshSize(0);
    ASSERT(mesh_size > 0);

    vector<double>::const_iterator mesh = func_int_->getMesh();

    Param0FunctionInt* c_func = C_->getParam0FunctionInt();
    double C = c_func->getValue();
    int clo_ind = 0;
    double clo_dist = fabs(C - mesh[0]);
    for (int ki = 1; ki < mesh_size; ++ki) {
	double dist = fabs(C - mesh[ki]);
	if (dist < clo_dist) {
	    clo_ind = ki;
	    clo_dist = dist;
	}
    }

    // We must now calculate the corresponding parameter value.
    // Assuming the grid is regular (evenly divided).
    double tmin = func_int_->startParam(0);
    double tmax = func_int_->endParam(0);
    double tstep = (tmax - tmin)/(mesh_size - 1);
    double clo_par = tmin + clo_ind*tstep;

    seed[0] = clo_par;
}


} // namespace Go
