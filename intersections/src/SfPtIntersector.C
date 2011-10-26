//===========================================================================
//                                                                           
// File: CvCvIntersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: SfPtIntersector.C,v 1.35 2007-11-01 14:31:45 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/SfPtIntersector.h"
#include "GoTools/intersections/CvPtIntersector.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/utils/RotatedBox.h"


using std::vector;
using std::cout;
using std::shared_ptr;


namespace Go
{


//===========================================================================
SfPtIntersector::SfPtIntersector(shared_ptr<ParamGeomInt> obj1,
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
  ParamSurfaceInt *surf = obj1->getParamSurfaceInt();
  if (surf == 0)
    {
      sf_idx_ = 1;
      pt_idx_ = 0;
    }
  else
    {
      sf_idx_ = 0;
      pt_idx_ = 1;
    }
}


//===========================================================================
SfPtIntersector::~SfPtIntersector()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
SfPtIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* prev,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
  // Set up intersection problem between two points
  shared_ptr<CvPtIntersector> curr_inter; 
    curr_inter = shared_ptr<CvPtIntersector>(new CvPtIntersector(obj1, 
								 obj2, 
								 epsge_, 
								 prev,
								 eliminated_parameter,
								 eliminated_value));

  return curr_inter;
}


//===========================================================================
int SfPtIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between a point and a surface can occur only if the
  // surface has a degenerate boundary. Make sure not to connect in the 
    // case of a closedsurface. Check midpoint.
    int ki;
    vector<shared_ptr<IntersectionPoint> > ints;
    int_results_->getIntersectionPoints(ints);

    // Remove points not lying at a degenerate boundary. First fetch surface.
    ParamSurfaceInt *surf = obj_int_[sf_idx_]->getParamSurfaceInt();
    double epsge = epsge_->getEpsge();
    double epspar = epsge_->getRelParRes();
    bool at_deg;
    for (ki=0; ki< int(ints.size()); ki++)
    {
	at_deg = surf->atDegenerateBd((sf_idx_== 0) ? ints[ki]->getPar1() : ints[ki]->getPar2(),
					   epsge, epspar);
	if (!at_deg && !selfint_case_)
	{
	    ints.erase(ints.begin()+ki);
	    ki--;
	}
    }

    //static bool do_connect = false;
    if (ints.size() != 2 || (!at_deg && !selfint_case_))
	return 0;

    // Check if the partial derivitives in the surface in the intersection
    // points indicates a coincidence, i.e. the angle bewteen the partial
    // derivatives is small
    Point pt1 = ints[0]->getPoint();
    Point pt2 = ints[1]->getPoint();
     double ang_tol = 0.1*epsge_->getAngleTol();
    vector<Point> der1(3), der2(3);
    surf->point(der1, (sf_idx_ == 0) ? ints[0]->getPar1() : ints[0]->getPar2(), 1);
    surf->point(der2, (sf_idx_ == 0) ? ints[1]->getPar1() : ints[1]->getPar2(), 1);
    double ang1 = der1[1].angle(der1[2]);
    double ang2 = der2[1].angle(der2[2]);
    if (ang1 < ang_tol || fabs(M_PI-ang1) < ang_tol);
    else
	return 0; 
    if (ang2 < ang_tol || fabs(M_PI-ang2) < ang_tol);
    else
	return 0; 
   double param[2], clo_par[2], dist;
    for (ki=0; ki<2; ki++)
	param[ki] = 0.5*(ints[0]->getPar(ki) + ints[1]->getPar(ki));

    doIterate(clo_par, dist, param);

    // Check if we end up in one of the initial points
    double ptol = epsge_->getRelParRes();
    for (ki=0; ki<2; ki++)
    {
	if (fabs(clo_par[ki] - ints[0]->getPar(ki)) > ptol)
	    break;
    }
    if (ki == 2)
	return 0;

    for (ki=0; ki<2; ki++)
    {
	if (fabs(clo_par[ki] - ints[1]->getPar(ki)) > ptol)
	    break;
    }
    if (ki == 2)
	return 0;

    // An extra check on whether the found points is close to the initial ones
    double fac = 1.0e-4;
    Point par1(ints[0]->getPar(0),ints[0]->getPar(1));
    Point par2(ints[1]->getPar(0),ints[1]->getPar(1));
    double dist_pnt = par1.dist(par2);
    Point close(clo_par[0],clo_par[1]);
    if (close.dist(par1) < fac*dist_pnt || close.dist(par2) < fac*dist_pnt)
	return 0;  // The new point is close enough to one of the initial
                   // ones to try some more subdivision

    // Check if we may have reached a different branch or the zero-curve is very
    // curved. In that case, some more subdivision may pay off
    Point guess(param[0],param[1]);
    // double fac2 = 0.75;
    double fac3 = 0.05;
    if (close.dist(guess) > fac3*close.dist(par1) ||
	close.dist(guess) > fac3*close.dist(par2))
	return 0;

    Point ptmid;
    obj_int_[sf_idx_]->point(ptmid, clo_par);

    if (dist < epsge_->getEpsge() && pt1.dist(ptmid) + ptmid.dist(pt2) < epsge_->getEpsge())
    {
	    
	if (selfint_case_)
	{
	    // Make an extra check on the configuration of the found intersection
	    // points
	    double ta[2], tb[2], tdel[2];
	    double frac = 0.5;
	    int ki;
	    for (ki=0; ki<2; ki++)
	    {
		ta[ki] = obj_int_[sf_idx_]->startParam(ki);
		tb[ki] = obj_int_[sf_idx_]->endParam(ki);
		tdel[ki] = fabs(ints[0]->getPar(ki) - ints[1]->getPar(ki));
	    }
	    if (tdel[0] < frac*(tb[0] - ta[0]) ||
		tdel[1] < frac*(tb[1] - ta[1]))
	    {
		// Only a small part of the surface is close to the 
		// intersection points. Check constant parameter curves for 
		// simplicity of surface.
		ParamSurfaceInt *sf = obj_int_[sf_idx_]->getParamSurfaceInt();
		shared_ptr<ParamSurface> psf = sf->getParamSurface();
		for (ki=0; ki<2; ki++)
		{
		    double par = 0.5*(ints[0]->getPar(ki) + 
				      ints[1]->getPar(ki));
		    vector<shared_ptr<ParamCurve> > ccrv = 
			psf->constParamCurves(par, (ki==1));
		    if (ccrv.size() != 1)
			return 0;
		    try {
			DirectionCone cone = ccrv[0]->directionCone();
			if (cone.greaterThanPi())
			    return 0;  // Complex surface
		    } catch (...)
		    {
			;  // Degenerate curve, not complex
		    }
		    
		}
	    }
	}

    	// Make sure that the points are connected
	ints[0]->connectTo(ints[1], COINCIDENCE_SFPT);
		
	return 1;
    }
    else 
	return 0;
}
 

//===========================================================================
void SfPtIntersector::microCase()
//===========================================================================
{
  // If no intersection point exist, make one
    vector<shared_ptr<IntersectionPoint> > ints;
    int_results_->getIntersectionPoints(ints);

    if (ints.size() == 0)
    {
	// The point is in the middle of the surface, not too exact, but
	// everything is very small.
	// The alternative would be to run an iteration.
	double mid[2];
	double dist;
	doIterate(mid, dist);
	if (dist < epsge_->getEpsge()) {
	    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
					       epsge_, mid, mid);
	}
    }
  else if (ints.size() == 1)
    {
      // One intersection point. Leave it like this.
      ;
    }
    else if (ints.size() == 2)
    {
	// Make sure the intersection points are connected
	ints[0]->connectTo(ints[1], MICRO_SFPT);
    }
}


//===========================================================================
int SfPtIntersector::updateIntersections()
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

  // Iterate to the intersection point
  double clo_par[2], clo_dist;
  double ptol = 100.0*epsge_->getRelParRes();
  doIterate(clo_par, clo_dist);
  if (clo_dist <= epsge_->getEpsge())
    {
 	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << "SfPt. Int. pt. found " << clo_par[0] << " " << clo_par[1];
	    cout << ", dist: " << clo_dist;
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
	int idx = (obj_int_[0]->numParams() == 2) ? 0 : 1;
	if (clo_par[0]-obj_int_[idx]->startParam(0) < ptol || 
	    obj_int_[idx]->endParam(0)-clo_par[0] < ptol ||
	    clo_par[1]-obj_int_[idx]->startParam(1) < ptol || 
	    obj_int_[idx]->endParam(1)-clo_par[1] < ptol)
	{
	    if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
		cout << " Point dismissed" << std::endl;
	    }
	    return 0;
	}

	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << " Point registered" << std::endl;
	    }
      int_results_->addIntersectionPoint(obj_int_[0], 
					 obj_int_[1],
					 getTolerance(),
					 clo_par, clo_par); 

      return 1; // A new point is found
    }

  return 0;
}


void SfPtIntersector::doIterate(double clo_par[2], double& clo_dist, double *guess)
{
  // Iterate to an intersection point. We know that we have got one curve
  // and one point. Fetch the data instances
  Point pt;
  double tpar = 0.0;
  RectDomain domain;

  ParamSurfaceInt* sf_int = obj_int_[sf_idx_]->getParamSurfaceInt();
  shared_ptr<ParamSurface> surf = 
      sf_int->getParentParamSurface(domain);
  DEBUG_ERROR_IF(surf.get() == 0, "Error in data structure");
  obj_int_[pt_idx_]->point(pt, &tpar);  // The parameter value is not relevant 

  // Compute seed
  double seed[2];
  if (guess)
  {
      seed[0] = guess[0];
      seed[1] = guess[1];
  }
  else if (hasSingularityInfo())
  {
      shared_ptr<SingularityInfo> sing = getSingularityInfo();
      if (sing->hasPoint())
      {
	  seed[0] = sing->getParam(0);
	  seed[1] = sing->getParam(1);
      }
      else
	  getSeedIteration(seed);
  }
  else
      getSeedIteration(seed);

  // Perform iteration
  Point clo_pt;
  surf->closestPoint(pt, clo_par[0], clo_par[1], clo_pt, clo_dist,
		     epsge_->getEpsge(), &domain, seed);
}
    
///////////////////////////////////////////////////////////////////
//
//  Purpose    : Perform a rotated box test between a curve and a 
//               surface
//
//  Written by : Vibeke Skytt, 0105
//
/////////////////////////////////////////////////////////////////

int SfPtIntersector::performRotatedBoxTest(double eps1, double eps2)
{
    //double ang_tol = epsge_->getAngleTol();
  ParamPointInt *point = obj_int_[pt_idx_]->getParamPointInt();
  ParamSurfaceInt *surf = obj_int_[sf_idx_]->getParamSurfaceInt();
  DEBUG_ERROR_IF(point == 0 || surf == 0,
	   "Inconsistence in the data structure");

    //Only done for simple surfaces
    if (!surf->isSimple())
	return 1;

    // Define coordinate system for the rotation. First check if an
    // intersection point between the geometry objects is found
    vector<Point> axis(2), der(3);
    Point norm;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    //int nmb_rec = nmbRecursions();
    if (int_pts.size() > 0)
    {
	// Make coordinate system from the partial derivatives of the surface
	// in the (first) intersection point
	const double *param = (sf_idx_ == 0) ? int_pts[0]->getPar1() : 
	    int_pts[0]->getPar2();
	obj_int_[sf_idx_]->point(der, param, 1);
	axis[0] = der[1];
	norm = der[1].cross(der[2]);
//	axis[1] = norm.cross(axis[0]);
	axis[1] = norm;
	if (getenv("DEBUG_ROTATEDBOX") && *(getenv("DEBUG_ROTATEDBOX")) == '1')
	    std::cout << "Intersection point found " << std::endl;
	
    }
    else
    {
	surf->axisFromCorners(der[0], der[1]);
	norm = der[0].cross(der[1]);
//	axis[1] = norm.cross(axis[0]);
	axis[0] = der[0];
	axis[1] = norm;
	if (getenv("DEBUG_ROTATEDBOX") && *(getenv("DEBUG_ROTATEDBOX")) == '1')
	    std::cout << "No intersection point " << std::endl;
	
    }
  
    // Make rotated boxes
    if (axis[0].length() < epsge_->getNumericalTol() ||
	axis[1].length() < epsge_->getNumericalTol()) 
	return 1;
    axis[0].normalize();
    axis[1].normalize();
    if (axis[0].angle_smallest(axis[1]) < epsge_->getNumericalTol())
	return 1;

    RotatedBox box1 = surf->getRotatedBox(axis);
    RotatedBox box2 = point->getRotatedBox(axis);
    bool overlap;
    overlap = box1.overlaps(box2, eps1, 0.0);
    if (overlap)
    {
	overlap = box1.overlaps(box2, eps1, eps2);
	if (foundIntersectionNearBoundary())
	    return 1;
    }
    if (overlap)
    {
	std::swap(axis[0], axis[1]);
	box1 = surf->getRotatedBox(axis);
	box2 = point->getRotatedBox(axis);
	overlap = box1.overlaps(box2, eps1, 0.0);
	if (overlap)
	{
	    overlap = box1.overlaps(box2, eps1, eps2);
	    if (foundIntersectionNearBoundary())
		return 1;
	}
    }
    return (overlap) ? 1 : 0;
}

///////////////////////////////////////////////////////////////////
//===========================================================================
//
// Purpose : Given one linear parametric surface and one point, compute the intersection.
//
// Written by : Vibeke Skytt, 1204
//
//===========================================================================
int SfPtIntersector::linearCase()
//===========================================================================
{
  return 0;
}

///////////////////////////////////////////////////////////////////
//===========================================================================
//
// Purpose : Given one parametric surface, one point and no simple case 
//           situation, subdivide the surface to produce more simple sub 
//           problems. Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Vibeke Skytt, 1104
//
//===========================================================================
int SfPtIntersector::doSubdivide()
//===========================================================================
{
  int perm[2];  // Two parameter directions that need sorting
  int nmb_subdiv;  // Number of parameter directions in which to subdivide
  int nmbdir1 = obj_int_[0]->numParams();
  int ki, kj, kr;

  // Sort the parameter directions according to importance of subdivision 
  // according to properties of the objects and already computed intersections
  nmb_subdiv = sortParameterDirections(perm);

  if (nmb_subdiv == 0)
    return 0;   // Not possible to subdivide neither the curve nor
                // the surface

  // Intersection objects belonging to the next recursion level
  // First they point to the objects at this level. Then the current objects
  // are replaced as subdivision is taking place.
  int nbobj[2];
  vector<shared_ptr<ParamGeomInt> > sub_objects;
  for (ki=0; ki<2; ki++)
    {
      sub_objects.push_back(obj_int_[ki]);
    }
  nbobj[0] = nbobj[1] = 1;    // Currently one instance of each intersection
                              // object at this level

  // For each parameter direction in prioritized order, fetch an
  // appropriate subdivision parameter, and perform subdivision
  double subdiv_par;
  SubdivisionClassification found;
  vector<shared_ptr<ParamGeomInt> > obj_sub;
  vector<shared_ptr<ParamGeomInt> > subdivobj;

  for (ki=0; ki<nmb_subdiv; ki++)
    {
      found = getSubdivisionParameter(perm[ki], subdiv_par);
      if (found == CANNOT_SUBDIVIDE)
	{
	  // No parameter value is found. Move this parameter direction
	  // to the end of the permutation array, and decrease the
	  // number of parameter directions where it is possible to 
	  // subdivide.
	}

      // Subdivide the current object (or objects. If the current object
      // is the surface then it can have been subdivided in one parameter
      // direction already and we have two sub surfaces to subdivide)
      int idxobj = (perm[ki] >= nmbdir1);
      int idx = idxobj*nbobj[0];
      for (kj=0; kj<nbobj[idxobj]; kj++)
	{
	    obj_sub.clear();
	    subdivobj.clear();
	    try {
		sub_objects[idx+kj]->subdivide(perm[ki]-idxobj*nmbdir1, 
					       subdiv_par, obj_sub, subdivobj);
	    } catch (...)
	    {
		obj_sub.clear();
		subdivobj.clear();
	    }
	  if (obj_sub.size() < 1 || subdivobj.size() == 0)
	    continue;  // No new objects 

	  for (kr=0; kr<int(obj_sub.size()); kr++)
	      obj_sub[kr]->setParent(obj_int_[idxobj].get());
      
	  for (kr = 0; kr<int(subdivobj.size()); kr++)
	    {
	      // Intersect the subdivision object with the other object
	      // @@@ VSK Faktorenes orden er ikke likegyldig!!
	      shared_ptr<Intersector> subdiv_intersector = 
		lowerOrderIntersector((idxobj==0) ? subdivobj[kr] : 
				      obj_int_[0], 
				      (idxobj==0) ? obj_int_[1] :
				      subdivobj[kr], 
				      this, perm[ki], subdiv_par);

	      // Is it here relevant to fetch existing intersection points
	      // and/or insert points into intersection curves before computing
	      // intersection points with the subdivision points. Normally, this
	      // is mostly of interest when surfaces are involved, but intersection
	      // intervals might exist?
	      // These computations do anyway involve the intersection pool, but
	      // they might be trigged somehow. The parameter direction and value
	      // are required information

	      if (found == DIVIDE_SING && hasSingularityInfo())
	      {
		  // Transfer singularity information to the child
		  // process
		  subdiv_intersector->setSingularityInfo(getSingularityInfo(), 
							 perm[ki]);
	      }

	      subdiv_intersector->compute(false);

	  // Check quality of intersection points
	  // If the quality is not sufficient, find a new subdivision parameter
	  // and repeat the process.
	  // Otherwise
	      int_results_->includeReducedInts(subdiv_intersector->getIntPool());
	    }

	  // Replace pointers to current objects with pointers to the
	  // sub objects after subdivision. The new objects can be the 
	  // final sub objects or they can be subdivided once more depending
	  // on the number of parameter directions elected for subdivision
	  // in the current object.
	  sub_objects[idx+kj] = obj_sub[0];
	  sub_objects.insert(sub_objects.begin()+idx+kj+1, obj_sub.begin()+1,
			      obj_sub.end());
	  nbobj[idxobj] += ((int)obj_sub.size()-1);
	  kj += ((int)obj_sub.size()-1);

	}
    }

  // Create new intersector objects
  for (ki=0; ki<nbobj[0]; ki++)
    for (kj=0; kj<nbobj[1]; kj++)
      {
	shared_ptr<Intersector> intersector = 
	  shared_ptr<Intersector>(new SfPtIntersector(sub_objects[ki],
						  sub_objects[nbobj[0]+kj],
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
int SfPtIntersector::sortParameterDirections(int perm[])
//===========================================================================
{
  // @@@ VSK. This routine is extremely similar to the one for curve-curve
  // intersections. Should it be moved one level up in the hierarchy

  double length[2], wiggle[2];  // Estimated length and wiggliness
  bool inner_knots[2], critical_val[2], can_divide[2];
  bool has_inner_ints[2];
  double rel_wiggle = 0.1;
  double rel_length = 0.1;
  int nmbdir[2];
  int ki, kj, kr;

  nmbdir[0] = obj_int_[0]->numParams();
  nmbdir[1] = obj_int_[1]->numParams();

  // Fetch information from the curve and surface
  if (nmbdir[0] > 0)
      obj_int_[0]->getLengthAndWiggle(length, wiggle);
  if (nmbdir[1] > 0)
  obj_int_[1]->getLengthAndWiggle(length+nmbdir[0], wiggle+nmbdir[0]);

  for (ki=0, kr=0; ki<2; ki++)
    for (kj=0; kj<nmbdir[ki]; kj++, kr++)
      {
	inner_knots[kr] = obj_int_[ki]->hasInnerKnots(kj);

	critical_val[kr] = obj_int_[ki]->hasCriticalVals(kj);

	can_divide[kr] = obj_int_[ki]->canDivide(kj);

	has_inner_ints[kr] = int_results_->hasPointsInInner(kr);
      }

  double med_len = 0.0, med_wiggle = 0.0;
  for (ki=0; ki<2; ki++)
  {
      med_len += length[ki];
      med_wiggle += wiggle[ki];
  }
  med_len /= (double)2;
  med_wiggle /= (double)2;

  double min_length = std::max(0.5*med_len, epsge_->getEpsge());
  double min_wiggle = std::max(0.5*med_wiggle, 0.02);

  // Fetch information from the intersection pool according to inner 
  // intersection points

  // Number of parameter directions is two.
  int size = 2;

  int curr = 0;
  int min_nmb = 0;

  // Initiate permutation array
  for (ki=0; ki<size; ki++)
    perm[ki] = ki;

  // Sort according to the given values
  // First sort out the directions where subdivsion is impossible
  for (ki=0; ki<size; ki++)
    {
      if (!can_divide[perm[ki]])
	{
	  if (perm[ki] < size-1)
	    std::swap(perm[ki], perm[size-1]);
	  size--;
	}
    }

  // First priority is parameter directions with critical values
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (critical_val[perm[kj]] && !critical_val[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (critical_val[perm[ki]])
	  curr++;

  // Next criterium is inner knots
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (inner_knots[perm[kj]] && !inner_knots[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (inner_knots[perm[ki]])
	  curr++;


  // Existing intersection points internal to the parameter interval
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (has_inner_ints[perm[kj]] && !has_inner_ints[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (has_inner_ints[perm[ki]])
	  curr++;

  // Finally length and curvature
  min_nmb = curr;
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if ((length[perm[kj]] > length[perm[ki]] && 
	   wiggle[perm[kj]] > wiggle[perm[ki]]) ||
	  (0.2*length[perm[kj]] > length[perm[ki]] && 
	   (wiggle[perm[kj]] > rel_wiggle*wiggle[perm[ki]] ||
	       wiggle[perm[ki]] < min_wiggle)) ||
	  (0.2*wiggle[perm[kj]] > wiggle[perm[ki]] && 
	   (length[perm[kj]] > rel_length*length[perm[ki]] ||
	    length[perm[ki]] < min_length)))
	{
	  std::swap(perm[ki], perm[kj]);
	}

  // Check that the minimum conditions for subdivision is satisfied
  for (ki=size-1; ki>=min_nmb; ki--)
      if (length[perm[ki]] < min_length && wiggle[perm[ki]] < min_wiggle)
	  size--;

  return size;
}


//===========================================================================
SubdivisionClassification
SfPtIntersector::getSubdivisionParameter(int dir, double& par)
//===========================================================================
{
    // @@@ VSK. This routine is extremely similar to the one for
    // curve-curve intersections. Should it be moved one level up in the
    // hierarchy

    int ki, kj, sgn;
//     double ptol = 100.0*epsge_->getRelParRes();
    double gtol = 100.0*epsge_->getEpsge();
    int nmbdir1 = obj_int_[0]->numParams();
    //int nmbdir2 = obj_int_[1]->numParams();

    // Set pointer to the intersection object corresponding to the parameter
    // direction
    ParamGeomInt *obj = 
	(dir < nmbdir1) ? obj_int_[0].get() : obj_int_[1].get();
    int pdir = (dir < nmbdir1) ? dir : dir-nmbdir1;
    double ta = obj->startParam(pdir);
    double tb = obj->endParam(pdir);
    double frac = 0.1*(tb - ta);

    // Get critical parameters
    // @@@ Critical parameters of different priority? Sorting?
    vector<double> critical_pars = obj->getCriticalVals(pdir);

    int size = (int)critical_pars.size();
    int is_critical = 0;
    if (size > 0) {
	// Check suitability of the critical parameters
	for (ki=0; ki<size; ki++) {
// 	    is_critical
// 		= int_results_->inInfluenceArea(dir, critical_pars[ki]);
	    if (is_critical == 0 || is_critical == 2)
		{
		    par = critical_pars[ki];
		    return DIVIDE_CRITICAL;
		}
	}
    }


    // Check for inner intersection points in which to subdivide
    vector<double> inner_ints = int_results_->getSortedInnerInts(dir);
    size = (int)inner_ints.size();
    if (size > 0){
	// Check suitability of the intersection points
	for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
	     ki+=sgn*kj, kj++, sgn*=-1) {
// 	    is_critical = int_results_->inInfluenceArea(dir, inner_ints[ki]);
	    par = inner_ints[ki];
	    if (par > ta + frac && par < tb - frac) 
	    {
		return DIVIDE_INT;
	    }
	}
    }

    // Iterate for a closest point in which to subdivide
    if (singularity_info_.get() == 0) {
	// Make empty singularity info instance
	singularity_info_
	    = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }
  
    // Check if a closest point exist already
    double param[2], dist=0.0;
    if (singularity_info_->hasPoint()) {
	par = singularity_info_->getParam(dir);
    } else {
	// Iterate for a closest point in which to subdivide
	doIterate(param, dist);
	if (false) {
// 	if (dist < gtol && param[dir] > ta+frac && param[dir] < tb-frac) {
	    // Doesn't work
	    // Look for a singular point nearby
	    Point sing_pt, sf_pt, pt;
	    double sing_par[2], sing_val, sing_dist;
	    double tptpar = 0.0;
	    obj_int_[pt_idx_]->point(pt, &tptpar);
	    ParamSurfaceInt *sf = obj_int_[sf_idx_]->getParamSurfaceInt();
	    sf->getSingularity(epsge_->getEpsge(), sing_par, sing_pt, 
			       sing_val, param);
	    sf->point(sf_pt, sing_par);
	    sing_dist = sf_pt.dist(pt);
	    if (sing_val < gtol && sing_dist <= dist) {
		param[0] = sing_par[0];
		param[1] = sing_par[1];
	    }
	}
	    
	singularity_info_->setSingularPoint(param, 2);
	par = param[dir];
    }
    if (dist < gtol && par > ta+frac && par < tb-frac) {
	return DIVIDE_SING;
    }
    
    // Look for a suitable knot in which to subdivide. First fetch the
    // knots sorted according to the distance to the mid-parameter
    vector<double> knot_vals = obj->getInnerKnotVals(pdir, true);
    size = (int)knot_vals.size();
    if (size > 0) {
	// Check suitability of the knots
	for (ki=0; ki<size; ki++) {
// 	    is_critical = int_results_->inInfluenceArea(dir, knot_vals[ki]);
	    if (is_critical == 0 || is_critical == 2) {
		par = knot_vals[ki];
		return DIVIDE_KNOT;
	    }
	}
    }

    // Subdivide at a suitable parameter
    double divpar = 0.5*(ta+tb);
    double del = 0.1*(tb-ta);
    double tint = del;
    for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=del, sgn*=-1) {
// 	is_critical = int_results_->inInfluenceArea(dir, divpar);
	if (is_critical == 0 || is_critical == 2) {
	    par = divpar;
	    return DIVIDE_PAR;
	}
    }

    return CANNOT_SUBDIVIDE;  // No subdivision parameter found
}


//===========================================================================


} // namespace Go
