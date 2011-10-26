//===========================================================================
//                                                                           
// File: CvCvIntersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: CvPtIntersector.C,v 1.32 2006-12-15 13:34:57 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/CvPtIntersector.h"
#include "GoTools/intersections/PtPtIntersector.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionPool.h"

using std::vector;
using std::cout;
using std::shared_ptr;


namespace Go
{


//===========================================================================
CvPtIntersector::CvPtIntersector(shared_ptr<ParamGeomInt> obj1,
				 shared_ptr<ParamGeomInt> obj2,
				 shared_ptr<GeoTol> epsge, 
				 Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)

  : Intersector2Obj(obj1, obj2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  // Set array indices to the object array. We know that one object is
  // a curve an the other is a point, but need to find out which one
  // is which
  ParamCurveInt *curve = obj1->getParamCurveInt();
  if (curve == 0)
    {
      cv_idx_ = 1;
      pt_idx_ = 0;
    }
  else
    {
      cv_idx_ = 0;
      pt_idx_ = 1;
    }
}


//===========================================================================
CvPtIntersector::CvPtIntersector(shared_ptr<ParamGeomInt> obj1,
				 shared_ptr<ParamGeomInt> obj2,
				 double epsge, 
				 Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)

  : Intersector2Obj(obj1, obj2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  // Set array indices to the object array. We know that one object is
  // a curve an the other is a point, but need to find out which one
  // is which
  ParamCurveInt *curve = obj1->getParamCurveInt();
  if (curve == 0)
    {
      cv_idx_ = 1;
      pt_idx_ = 0;
    }
  else
    {
      cv_idx_ = 0;
      pt_idx_ = 1;
    }
}


//===========================================================================
CvPtIntersector::~CvPtIntersector()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
CvPtIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* prev,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
  // Set up intersection problem between two points
  shared_ptr<PtPtIntersector> curr_inter; 
    curr_inter = shared_ptr<PtPtIntersector>(new PtPtIntersector(obj1, 
								 obj2, 
								 epsge_, 
								 prev,
								 eliminated_parameter,
								 eliminated_value));

  return curr_inter;
}


//===========================================================================
int CvPtIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between a point and a curve can occur only if the
  // curve is degenerate. 
    return 0;
}


void CvPtIntersector::microCase()
{

  // If no intersection point exist, make one
    vector<shared_ptr<IntersectionPoint> > ints;
    int_results_->getIntersectionPoints(ints);

    if (ints.size() == 0)
    {
	// The point is in the middle of the surface, not too exact, but
	// everything is very small.
	// iterate down to an intersection
	double start, end;
	Point pt, clo_pt;
	double clo_par, clo_dist;
	double tpar = double(0);
	ParamCurveInt* cv_int = obj_int_[cv_idx_]->getParamCurveInt();
	shared_ptr<ParamCurve> curve = cv_int->getParentParamCurve(start, end);
	DEBUG_ERROR_IF(curve.get() == 0, "Error in data structure.");
	obj_int_[pt_idx_]->point(pt, &tpar); // parameter value irrelevant
	double mid;
	int_results_->getMidParameter(&mid);
 	curve->closestPoint(pt, start, end, clo_par, clo_pt, clo_dist, &mid);
 	if (clo_dist < epsge_->getEpsge()) {
	    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], epsge_, &clo_par, &clo_par);
	} 
    }
  else if (ints.size() == 1)
    {
      // One intersection point. Leave it like this.
      ;
    }
    else
    {
      // More than one intersection point. Connect.
      for (int ki=1; ki<int(ints.size()); ki++)
	ints[ki-1]->connectTo(ints[ki], COINCIDENCE_CVPT);
    }
	
}


//===========================================================================
int CvPtIntersector::updateIntersections()
//===========================================================================
{
  // Fetch already existing intersection points
  vector<shared_ptr<IntersectionPoint> > int_pts;
  int_results_->getIntersectionPoints(int_pts);

  if (int_pts.size() > 0)
    {
      // At most one intersection point is expected. One or more points
      // exist already. Leave it like this.
      return 0;
    }

  // Iterate to an intersection point. We now that we have got one curve
  // and one point. Fetch the data instances
  double start, end;
  Point pt;
  double tpar = 0.0;
  double clo_par, clo_dist;
  double *par1=0, *par2=0;

  ParamCurveInt* cv_int = obj_int_[cv_idx_]->getParamCurveInt();
  shared_ptr<ParamCurve> curve = 
      cv_int->getParentParamCurve(start, end);
  DEBUG_ERROR_IF(curve.get() == 0, "Error in data structure");
  obj_int_[pt_idx_]->point(pt, &tpar);  // The parameter value is not relevant 
  if (cv_idx_ == 0)
    par1 = &clo_par;
  else
    par2 = &clo_par;

  double seed;
  if (hasSingularityInfo())
  {
      shared_ptr<SingularityInfo> sing = getSingularityInfo();
      if (sing->hasPoint())
	  seed = sing->getParam(0);
      else
	  getSeedIteration(&seed);
  }
  else
      getSeedIteration(&seed);

  // Perform iteration
  double ptol = 100.0*epsge_->getRelParRes();
  Point clo_pt;
  curve->closestPoint(pt, start, end, clo_par, clo_pt, clo_dist, &seed);
    
  if (clo_dist <= epsge_->getEpsge())
    {
	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << "CvPt. Int. pt. found " << clo_par << ", dist: " << clo_dist;
	}
      // An intersection point is found. Represent it in the data
      // structure
      // @@@ vsk, In sisl there is a check if the found intersection
      // point is very close to the endpoints of the curve. In that case
      // it is dismissed since the intersections at the endpoints are found
      // already. Here we have already checked if there exist any
      // intersection points. Thus, we know that there does not exist any
      // intersection point at the endpoint.
	// VSK, 0906 But maybe the endpoint already is in intersection, but that
	// point is represented elsewhere at a more exact position. Make a check.

	if (clo_par-obj_int_[cv_idx_]->startParam(0) < ptol || 
	    obj_int_[cv_idx_]->endParam(0)-clo_par < ptol)
	{
	    if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
		cout << " Point dismissed" << std::endl;
	    }
	    return 0;
	}

	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << " Point registered" << std::endl;
	    }
	shared_ptr<IntersectionPoint> tmp = 
	    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], 
					       epsge_, par1, par2);
      return 1; // A new point is found
    }

  return 0;
}


//===========================================================================
//
// Purpose : Given one linear parametric curve and one point, compute the intersection.
//
// Written by : Vibeke Skytt, 1204
//
//===========================================================================
int CvPtIntersector::linearCase()
//===========================================================================
{
  return 0;
}

///////////////////////////////////////////////////////////////////
//===========================================================================
//
// Purpose : Given one parametric curve, one point and no simple case 
//           situation, subdivide the surface to produce more simple sub 
//           problems. Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Vibeke Skytt, 1104
//
//===========================================================================
int CvPtIntersector::doSubdivide()
//===========================================================================
{
  int nmb_subdiv;  // Number of parameter directions in which to subdivide
  //int nmbdir1 = obj_int_[0]->numParams();
  int ki, kj, kr;

  // Sort the parameter directions according to importance of subdivision 
  // according to properties of the objects and already computed intersections
  nmb_subdiv = sortParameterDirections();

  if (nmb_subdiv == 0)
    return 0;   // Not possible to subdivide the curve 

  // Intersection objects belonging to the next recursion level
  vector<shared_ptr<ParamGeomInt> > sub_objects1;
  vector<shared_ptr<ParamGeomInt> > sub_objects2;

  // For each parameter direction in prioritized order, fetch an
  // appropriate subdivision parameter, and perform subdivision
  double subdiv_par;
  SubdivisionClassification found;
  vector<shared_ptr<ParamGeomInt> > obj_sub;
  vector<shared_ptr<ParamGeomInt> > subdivobj;

  for (ki=0; ki<nmb_subdiv; ki++)
    {
      found = getSubdivisionParameter(subdiv_par);
      if (found == CANNOT_SUBDIVIDE)
	{
	  // No parameter value is found. Move this parameter direction
	  // to the end of the permutation array, and decrease the
	  // number of parameter directions where it is possible to 
	  // subdivide.
	}

      // Subdivide the current curve
      obj_sub.clear();
      subdivobj.clear();
      try {
      obj_int_[cv_idx_]->subdivide(0, subdiv_par, obj_sub, subdivobj);
      } catch (...)
      {
	  obj_sub.clear();
	  subdivobj.clear();
      }
	  
	  if (obj_sub.size() < 1 || subdivobj.size() == 0)
	    continue;  // No new objects 

	  for (kr = 0; kr<int(subdivobj.size()); kr++)
	    {
	      // Intersect the subdivision object with the other object
	      // @@@ VSK Faktorenes orden er ikke likegyldig!!
		shared_ptr<Intersector> subdiv_intersector
		    (new PtPtIntersector(cv_idx_==0 ? subdivobj[kr] : obj_int_[0],
					 cv_idx_==1 ? subdivobj[kr] : obj_int_[1],
					 epsge_, this, 0, subdiv_par));
							   
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
	    }

	  if (cv_idx_ == 0)
	      sub_objects1.insert(sub_objects1.end(), obj_sub.begin(), 
				  obj_sub.end());
	  else
	      sub_objects2.insert(sub_objects2.end(), obj_sub.begin(), 
				  obj_sub.end());
    }

  // Create new intersector objects
  if (sub_objects1.size() == 0)
      sub_objects1.push_back(obj_int_[0]);
  if (sub_objects2.size() == 0)
      sub_objects2.push_back(obj_int_[1]);
  for (ki=0; ki<int(sub_objects1.size()); ki++)
    for (kj=0; kj < int(sub_objects2.size()); kj++)
      {
	shared_ptr<Intersector> intersector = 
	  (shared_ptr<Intersector>)(new CvPtIntersector(sub_objects1[ki],
							sub_objects2[kj],
							epsge_, this));
	sub_intersectors_.push_back(intersector);
      }

  return 1;
}

//===========================================================================
//
// Purpose : Fetch information related to the objects in order to decide
//           which curves that can be subdivided and which parameter direction
//           is most important to subdivide.
//
// Written by : Vibeke Skytt, 1004
//
//===========================================================================
int CvPtIntersector::sortParameterDirections()
//===========================================================================
{
    // Subdivide the curve if possible
    bool can_divide = obj_int_[cv_idx_]->canDivide(0);
    if (can_divide)
	return 1; 
    else 
	return 0;
}


//===========================================================================
//
// Purpose :
//
// Written by : Vibeke Skytt, 1004
//
//===========================================================================
SubdivisionClassification CvPtIntersector::getSubdivisionParameter(double& par)
//===========================================================================
{
  int ki, kj, sgn;
  //int nmbdir1 = obj_int_[0]->numParams();
  //int nmbdir2 = obj_int_[1]->numParams();

  // Set pointer to the intersection object corresponding to the parameter
  // direction
  ParamGeomInt *obj = obj_int_[cv_idx_].get();
  int pdir = 0;
  double ta = obj->startParam(pdir);
  double tb = obj->endParam(pdir);
  double fac = 0.1;
  double del = fac*(tb - ta);

  // Get critical parameters
  // @@@ Critical parameters of different priority? Sorting?
  vector<double> critical_pars = obj->getCriticalVals(pdir);

  int size = (int)critical_pars.size();
  //int is_critical;
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
	      return DIVIDE_CRITICAL;
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
	      return DIVIDE_KNOT;
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
	      return DIVIDE_INT;
	    }
	}
    }

  // Subdivide in the middle
  par = 0.5*(ta+tb);
  return DIVIDE_PAR;

//   double divpar = 0.5*(ta+tb);
//   double tdel = 0.1*(tb-ta);
//   double tint = del;
//   for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=tdel, sgn*=-1)
//     {
// 	par = divpar;
// //       is_critical = int_results_->inInfluenceArea(pdir, divpar);
// //       if (is_critical == 0 || is_critical == 2)
// 	if (par >= ta + del && par <= tb - del)
// 	{
// 	  return true;
// 	}
//     }

//   return false;  // No subdivision parameter found
}



} // namespace Go
