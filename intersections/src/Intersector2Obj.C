//===========================================================================
//                                                                           
// File: Intersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: Intersector2Obj.C,v 1.80 2007-06-22 16:33:32 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Intersector2Obj.h"
#include "GoTools/utils/CompositeBox.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/intersections/ParamGeomInt.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/Values.h"
#include <vector>

#include "GoTools/geometry/ObjectHeader.h" // for debugging
#include <fstream> // For debugging
#include <stdio.h> // for debugging
#include <iostream>
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
//#include <iostream> // @@debug purposes


using namespace Go;
using std::vector;
using std::cout;


//===========================================================================
Intersector2Obj::Intersector2Obj(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 shared_ptr<GeoTol> epsge,
				 Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector(epsge, prev)
//===========================================================================
{
  obj_int_[0] = obj1;
  obj_int_[1] = obj2;
  shared_ptr<IntersectionPool> parent_pool;
  if (prev) {
    parent_pool = prev->getIntPool();
  }
  int_results_ = 
    shared_ptr<IntersectionPool>(new IntersectionPool(obj1, 
						      obj2, 
						      parent_pool, 
						      eliminated_parameter,
						      eliminated_value));
  selfint_case_ = (prev) ? prev->isSelfintCase() : 0;
  if (prev && eliminated_parameter < 0)
  {
      singularity_info_ = (shared_ptr<SingularityInfo>)
	  (new SingularityInfo(prev->getSingularityInfo()));
  }
  if (hasSingularityInfo())
  {
      int ki, kj;
      int nmb = obj1->numParams() + obj2->numParams();
      if (nmb > 0)
      {
	  vector<double> start(nmb);
	  vector<double> end(nmb);
	  for (ki=0, kj=0; ki<obj1->numParams(); ki++, kj++)
	  {
	      start[kj] = obj1->startParam(ki);
	      end[kj] = obj1->endParam(ki);
	  }
	  for (ki=0; ki<obj2->numParams(); ki++, kj++)
	  {
	      start[kj] = obj2->startParam(ki);
	      end[kj] = obj2->endParam(ki);
	  }
      singularity_info_->cleanUp(start, end, epsge_->getRelParRes());
      }
  }
}


//===========================================================================
Intersector2Obj::Intersector2Obj(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 double epsge, Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector(epsge, prev)
//===========================================================================
{
  obj_int_[0] = obj1;
  obj_int_[1] = obj2;
  shared_ptr<IntersectionPool> parent_pool;
  if (prev) {
    parent_pool = prev->getIntPool();
  }
  int_results_ = 
    shared_ptr<IntersectionPool>(new IntersectionPool(obj1, 
						      obj2, 
						      parent_pool, 
						      eliminated_parameter,
						      eliminated_value));
  selfint_case_ = (prev) ? prev->isSelfintCase() : 0;
}


//===========================================================================
Intersector2Obj::~Intersector2Obj()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
void Intersector2Obj::print_objs()
//===========================================================================
{
}


//===========================================================================
int Intersector2Obj::performInterception()
//===========================================================================
{
    // Purpose: Do interception test between two parametric objects.

    int do_intercept;

    // Get the boxes belonging to the objects
    CompositeBox box1 = obj_int_[0]->compositeBox();
    CompositeBox box2 = obj_int_[1]->compositeBox();

    // Check overlap
    double eps = 0.5*epsge_->getEpsge();
    double overlap = 0.0;
//     if (box1.getOverlap(box2, overlap, eps, -0.5*eps))
    if (box1.overlaps(box2, eps, 0.0)) {

// 	// TESTING
// 	static int nmb_sep = 0;
// 	static int nmb_sep2 = 0;
// 	int nmbpar = numParams();
// 	int rec = nmbRecursions();
// 	bool sep = false;
   // 	if (nmbpar == 3) {
// 	    if (!box1.overlaps(box2, eps, -eps)) {
// 		nmb_sep++;
// 		sep = true;
// 	    }
// 	    if  (!box1.overlaps(box2, 0.0, 0.0)) {
// 		nmb_sep2++;
// 		sep = true;
// 	    }
// 	    if (sep) {
// 		std::cout << " Box test separation " << nmb_sep;
// 		std::cout << ", no fat: " << nmb_sep2 << ", recursion: ";
// 		std::cout << rec << std::endl;
// 	    }
// 	}
	
	// Check if the boxes both lie in an area of size equal to the
	// tolerance

	// Some check for this should be included in the bounding box
	Point low1 = box1.low(0.0, 0.0);
	Point high1 = box1.high(0.0, 0.0);
	Point low2 = box2.low(0.0, 0.0);
	Point high2 = box2.high(0.0, 0.0);
	int dim = box1.dimension();
	int ki;
	for (ki=0; ki<dim; ki++) {
	    if (std::max(high1[ki],high2[ki]) - std::min(low1[ki],low2[ki]) > 
		epsge_->getEpsge()) {
		break;
	    }
	}
	if (ki == dim) {
	    do_intercept = 2;
	} else {
	    // Remember statistical information
	    if (!hasComplexityInfo()) {
		complexity_info_
		    = shared_ptr<ComplexityInfo>(new ComplexityInfo());
	    }
	    complexity_info_->setBoxOverlap(overlap);

	    // Check with boxes reduced at the boundary
	    if (!box1.overlaps(box2, eps, -eps) &&
		!foundIntersectionNearBoundary()) {
		do_intercept = 0;
// 		if (sep)
// 		    std::cout << "Interception by boundary reduced"
// 			      << std::endl;
	    } else {
		// Check if a test on rotated boxes should be called
		int numpt = int_results_->numIntersectionPoints();
		if (numpt == 2) {
		    // Check if the two points are connected along a
		    // constant parameter direction
		    vector<shared_ptr<IntersectionPoint> > intpts;
		    int_results_->getIntersectionPoints(intpts);
		    if (intpts[0]->isConnectedTo(intpts[1])) {
			shared_ptr<IntersectionLink> link =
			    intpts[0]->getIntersectionLink(intpts[1].get());
			if (link->isIsoparametric())
			    numpt = 1; // In practice one intersection instance
		    }
		}

		do_intercept = 1;  // No interception yet
		if (numpt <= 1) {
// 		if (numpt <= 2) {

		    // At most one intersection point is found at the
		    // boundaries. Try rotated box tests
		    do_intercept = performRotatedBoxTest(eps, -eps);
	    
		    if (getenv("DEBUG_ROTATEDBOX") 
			&& *(getenv("DEBUG_ROTATEDBOX")) == '1')
			std::cout << "Rotated box called, result: "
				  << do_intercept << std::endl;
	    
		}
		
	    }
	}
    } else {
	do_intercept = 0;
    }
    if (getenv("DEBUG_BOX") && *getenv("DEBUG_BOX") == '1') {
	int npar = obj_int_[0]->numParams()
	    + obj_int_[1]->numParams();
	if (npar == 4 || (npar == 2 && getenv("SUBDIV_CVCV")
			  && *getenv("SUBDIV_CVCV") == '1')) {
	    bool touch = box1.overlaps(box2, 0.0);
	    std::cout << "Intercept " << do_intercept <<", almost "
		      << touch << std::endl;
	    if (do_intercept == 1 && !touch) {
		writeOut();
	    }
	}
    }

    return do_intercept;
}



//===========================================================================
int Intersector2Obj::complexIntercept()
//===========================================================================
{
    // Purpose: Do interception test between two parametric objects using
    //          methods related to approximate implicitization

    int not_intercept = 1;

    // Make a call to a virtual function that
    // initially only has a proper implementation in
    // SfSfIntersector. The other Intersectors just
    // returns no interception. This can be extended
    // later.

    // If both surfaces are simple (or at least one,
    // we can do some tuning on this, at least one
    // surface can be implicitized and finds an
    // implicite surface within the tolerance:

    // - Implicitize the "simplest surface". It must
    // be Bezier. If both are Bezier, use the one with
    // smalles cone angle

    // - Use, at least initially, a tolerance of
    // - 0.5*epsge_->getEpsge()

    // - Later we should also have a test on how
    // constant the size of the gradient of the
    // implicite surface is.

    // - The 1D tolerance must be made according to an
    // estimate of the gradient in the relevant area.
    // eps1D = eps3D * | grad H(F(u,v)| if H is the
    // implicit approximation of one surface and F is
    // the other surface Talk to Tor.

    // - Make the expression H(F(u,v), i.e. put the
    // coefficients of the other surface into the
    // implicit one. This is probably done in the
    // algebraic class. What do you do for planes and
    // spheres?

    // - Make the functional intersector

    // - Call performInterception in the functional
    // intersector to check if it possibly can be an
    // intersection. Thus, this function must be
    // accessible.

    // We create the implicit representation of the
    // objects.  I.e. we implizitice one of the sfs,
    // then plug the other into the impl sf.

    // If both objects are spline-sfs one of the sfs
    // will be implicitized.

    not_intercept = performInterceptionByImplicitization(); //1;
    //not_intercept = 1;
		    
    if (not_intercept == 1) {

	// @@@ TODO. This is a different implicitization
	// and a different virtual function in
	// SfSfIntersector.

	// The situation: An intersection curve is found
	// at the boundary between two surfaces. None of
	// the surfaces has a very large cone, probably
	// the cone angle is less than PI/3. The angle
	// between the cone centres is not very
	// large. Less than PI/6 or PI/4?  Jan has made a
	// new implicitization method that makes an
	// implicit surface that follows the common
	// boundary curve between the two surfaces
	// (represented as a spline curve) and is ruled
	// with respect to a vector which is computed as
	// the average of the cone centres. Class
	// ImplicitizeCurveAndVectorAlgo in
	// implicitization library.

	// The test applies only to spline surfaces (at
	// least for the time being).

	// It is tested if the spline coefficients of the
	// two surfaces, respectively, lie on different
	// sides of the implicit surface. If so, no more
	// intersections are possible.

	// I think much of the work is done in the
	// SfSfIntersector.

	// Test on whether there is an intersection curve
	// at the boundaries is done by the
	// IntersectionPool. Other pieces of code can be
	// placed in for instance SplineSurfaceInt or the
	// algebraic surface class, but I am not sure yet.

	// THIS IMPLICITIZATION IS AFTER THE OTHER ONE.

	not_intercept = interceptionBySeparationSurface();
    }
	    
    return not_intercept ? 0 : 1;
}

//===========================================================================
int Intersector2Obj::performRotatedBoxTest(double eps1, double eps2)
//===========================================================================
{
    // Purpose: Default implementation for rotated box test. Do
    // nothing

    return 1;   // Default result is overlap
}


//===========================================================================
int Intersector2Obj::simpleCase()
//===========================================================================
{
    //  Purpose: Do simple case test between two parametric objects.
    //  Testing between directional cones is performed


    // Get the direction cones.  Should the complete cones be computed
    // if they are larger than pi?
    DirectionCone cone1, cone2;
    try {
	cone1 = obj_int_[0]->directionCone();
    } catch (...) {

	// Check if the object is degenerate
	CompositeBox box = obj_int_[0]->compositeBox();
	Point low = box.low(0.0, 0.0);
	Point high = box.high(0.0, 0.0);
	int dim = box.dimension();
	int ki;
	for (ki=0; ki<dim; ki++) {
	    if (high[ki] - low[ki] > epsge_->getEpsge()) {
		// Can't accept this kind of degeneracy as simple case
		return 0;
	    }
	}
	Point centre(dim);
	centre.setValue(0.0);
	DirectionCone tmp(centre, 0.0);
	cone1 = tmp;
    }
    try {
	cone2 = obj_int_[1]->directionCone();
    } catch (...) {

	// Check if the object is degenerate
	CompositeBox box = obj_int_[1]->compositeBox();
	Point low = box.low(0.0, 0.0);
	Point high = box.high(0.0, 0.0);
	int dim = box.dimension();
	int ki;
	for (ki=0; ki<dim; ki++) {
	    if (high[ki] - low[ki] > epsge_->getEpsge()) {
		// Can't accept this kind of degeneracy as simple case
		return 0;
	    }
	}
	Point centre(dim);
	centre.setValue(0.0);
	DirectionCone tmp(centre, 0.0);
	cone2 = tmp;
    }

    if (getenv("DEBUG_CONE") && *getenv("DEBUG_CONE") == '1') {
	int npar = obj_int_[0]->numParams()
	    + obj_int_[1]->numParams();
	if (npar == 4 || (npar == 2 && getenv("SUBDIV_CVCV") &&
			  *getenv("SUBDIV_CVCV") == '1')) {
	    cout << "Simple case, angle1  " << cone1.angle();
	    cout << ", angle 2 " << cone2.angle() << " overlap ";
	    cout << cone1.overlaps(cone2) << std::endl;
	}
    }

    // Check if any of the cones are greater than pi. In that case,
    // it is no simple case.
    // NB! It might be that this test should be replaced by something
    // more strict in the future.
    if (cone1.greaterThanPi() || cone2.greaterThanPi()) {
// 	int simple_case = simpleCaseByImplicitization();
// 	if (simple_case == 1)
// 	    return 1;
	return 0;
    }

    // Check overlap.  @@@jbt - Sinha's Theorem.
    // Note that this test must be overruled for curve-surface
    // intersection. In that case the cones are expected to be
    // perpendicular
    if (cone1.overlaps(cone2)) {
	double overlap;
	Point axis1 = cone1.centre();
	Point axis2 = cone2.centre();
	double angle = axis1.angle_smallest(axis2);
	overlap = std::max(angle + cone1.angle() + cone2.angle() - M_PI,
			   cone1.angle() + cone2.angle() - angle);

	// Remember statistical information
	if (!hasComplexityInfo())
	    complexity_info_
		= shared_ptr<ComplexityInfo>(new ComplexityInfo());
	complexity_info_->setConeOverlap(overlap);

	// Extra test if there is a small overlap between the cones. 
	if (angle < M_PI - epsge_->getAngleTol() &&
	    angle > epsge_->getAngleTol() &&
	    cone1.angle() <= 1.3*angle &&
	    cone2.angle() <= 1.3*angle) {
	    int simple_case =  simpleCase2(axis1, axis2);
	    if (simple_case == 0 && numParams() == 4) {
		// Not simple. Make yet another try.
		bool has_bd_ints[8];  // Maximum number of boundaries
		for (int i = 0; i < 8; ++i) {
		    // We initialize in order to avoid that Valgrind
		    // reports this as an error. @jbt.
		    has_bd_ints[i] = false;
		}
		for (int ki = 0; ki < 2; ki++) {
		    int kr = 0;
		    for (int kj=0; kj<obj_int_[ki]->numParams(); kj++) {
			double par = obj_int_[ki]->startParam(kj);
			has_bd_ints[kr] = 
			    int_results_->hasIntersectionPoints(2*ki+kj, par);
			par = obj_int_[ki]->endParam(kj);
			has_bd_ints[kr+1] = 
			    int_results_->hasIntersectionPoints(2*ki+kj, par);
			kr += 2;
		    }
		}

		double eps = epsge_->getEpsge();
		DirectionCone red_cone1 = 
		    obj_int_[0]->reducedDirectionCone(has_bd_ints, eps);
		DirectionCone red_cone2 = 
		    obj_int_[1]->reducedDirectionCone
		    (has_bd_ints + 2*obj_int_[0]->numParams(), 
		     eps);

		if (red_cone1.greaterThanPi() < 0 
		    || red_cone2.greaterThanPi() < 0) {
		    // Failed to make a direction cone.
		    // This is a complex case
		    complexity_info_->setComplex();
		    return 0;
		}
		    
		if (red_cone1.overlaps(red_cone2)) {
		    simple_case = 0;
		} else {
		    simple_case = 1;
		}
	    }

	    if (simple_case == 1) {
		return 1;
	    }
	}
	return 0;   // Overlap and too large for further investigations
    }

    // Otherwise it is a simple case (Sinha's Theorem)
    return 1;
}


//===========================================================================
int Intersector2Obj::complexSimpleCase()
//===========================================================================
{
    //  Purpose: Do simple case test between two parametric objects.
    //  Methods of approximate implicitization is used.

    // Make a call to a virtual function that initially only has a
    // proper implementation in SfSfIntersector. The other
    // Intersectors just returns no simple case. This can be
    // extended later.
    // Most of the actions is the same as in
    // performInterception. Only relatively simple surfaces are
    // implicitized.
    // Monotone trend in the functional intersector is called.
    // Since the same implicit approximation of the same surface
    // can be used in different contexts, it is important to store
    // it to save some time. A natural place for storage is
    // ParamSurfaceInt.

    int simple_case = simpleCaseByImplicitization();
    //int simple_case = 0; // Simple case is not confirmed

    return simple_case;
}

//===========================================================================
int Intersector2Obj::simpleCase2(Point& axis1, Point& axis2)
//===========================================================================
{
    //  Purpose: Do extra simple case test between two parametric
    //  objects. Default implementation if no particular action is
    //  taken.

    return 0;
}


//===========================================================================
bool Intersector2Obj::isLinear()
//===========================================================================
{
    //  Purpose: Check if both objects involved in the intersection
    //  are linear.

    // For each object
    for (int ki=0; ki<2; ki++)
    {
	if (!obj_int_[ki]->isLinear(epsge_->getEpsge()))
	    return false;
    }
    return true;   // Both objects are linear
}


//===========================================================================
int Intersector2Obj::getBoundaryIntersections()
//===========================================================================
{

    int ki, kj;
    bool only_second = false;  // Compute boundary intersections only
                               // from the second object
    static bool do_post_iter = true;

    // Check if a previous intersector at the same level exist.  In
    // that case the boundary intersections are computed already, and
    // it is nothing to do.
    if (prev_intersector_ != 0
	&& prev_intersector_->numParams() == numParams()) {
	return 0;  // Not necessary to compute
    }

    // Check if a previous intersector of a higher level exist.  In
    // that case the same pair of boundary conditions may have been
    // computed already
    if (prev_intersector_ != 0) {

	// An intersector exist, and the number of parameter
	// directions is different.

	// Check if number two of the current objects belong to the
	// second array of boundary objects in the previous
	// intersector. In that case, the intersections on the
	// boundaries of the first of the current objects are treated
	// already.

	// NB! This depends on the sequencing on the boundary
	// computation operations. Always treat the first array of
	// boundary objects first.

	int nmb_bd = prev_intersector_->nmbBdObj(1);
	for (kj=0; kj<nmb_bd; kj++) {
	    ParamObjectInt *bd_obj =  
		prev_intersector_->getBoundaryObject(1,kj)->getObject().get();
	    if (obj_int_[1].get() == bd_obj) {
		only_second = true;  // Not necessary to compute the
	                             // boundary intersections, they
	                             // are computed already
	    }
	}
    }

    if (isSelfintCase() == 1)
	only_second = false;

    // The boundary intersections must be computed.
    // Make the appropriate intersectors
    // Get the boundary objects
    // Keep the boundary intersectors until this function goes out of
    // scope to be able to search through already defined intersector
    // for relevant intersection results.

    shared_ptr<Intersector> bd_intersector;
    std::vector<shared_ptr<BoundaryGeomInt> > bd_obj1;
    std::vector<shared_ptr<BoundaryGeomInt> > bd_obj2;
    if (!only_second)
	obj_int_[0]->getBoundaryObjects(bd_obj1);
    obj_int_[1]->getBoundaryObjects(bd_obj2);
    
    int nmbpar1, nmbpar2;
    nmbpar1 = obj_int_[0]->numParams();
    nmbpar2 = obj_int_[1]->numParams();
    for (ki=0; ki < int(bd_obj1.size()); ki++) {
	bd_intersector = lowerOrderIntersector(bd_obj1[ki]->getObject(), 
					       obj_int_[1], 
					       this,
					       bd_obj1[ki]->getDir(),
					       bd_obj1[ki]->getPar());
	bd_intersector->compute();
	int nmb_orig = int_results_->numIntersectionPoints();
	int_results_->includeReducedInts(bd_intersector->getIntPool());
	// Post iterate new points
// 	if (nmb_orig < int_results_->numIntersectionPoints())
// 	    postIterate(nmb_orig, bd_obj1[ki]->getDir());

	// Avoid false connections in the degenerate case
	removeDegenerateConnections();

	// Remove inner points in constant parameter intersection links
	// (vsk, 0609) and isolated points identical to the existing ones
	int_results_->cleanUpPool(nmb_orig, epsge_->getEpsge());
    }
    if (do_post_iter)
    {
    for (ki=0; ki<numParams(); ki++)
      postIterate(0, ki);
    }
  
    for (ki=0; ki < int(bd_obj2.size()); ki++) {
	bd_intersector = lowerOrderIntersector(obj_int_[0], 
					       bd_obj2[ki]->getObject(),
					       this,
					       bd_obj2[ki]->getDir()+nmbpar1,
					       bd_obj2[ki]->getPar());
	bd_intersector->compute();
	int nmb_orig = int_results_->numIntersectionPoints();
	int_results_->includeReducedInts(bd_intersector->getIntPool());
// 	if (nmb_orig < int_results_->numIntersectionPoints())
// 	    postIterate(nmb_orig, bd_obj1[ki]->getDir());

	// Avoid false connections in the degenerate case
	removeDegenerateConnections();

	// Remove inner points in constant parameter intersection
	// links
	// (vsk, 0609) and isolated points identical to the existing ones
	int_results_->cleanUpPool(nmb_orig, epsge_->getEpsge());
    }
    if (do_post_iter)
    {
    for (ki=0; ki<numParams(); ki++) {
      postIterate(0, ki);
      }
    }

    return 1;  // Boundary intersections computed
}


//===========================================================================
void Intersector2Obj::getSeedIteration(double seed[])
//===========================================================================
{
    //  Purpose: Compute seed to iteration routines based on mesh
    //  fetched from the intersecton objects. The mesh can be a
    //  control polygon, but can also be a rough facetted
    //  representation of the object.

    int nmbel1[2], nmbel2[2];  // Maximum of two parameter directions
                               // for each object
    int numpar1 = obj_int_[0]->numParams();
    int numpar2 = obj_int_[1]->numParams();
    int dim = obj_int_[0]->dimension();
    vector<double>::iterator mesh1, mesh2;
    double ta[4], tb[4];

    // Initialize number mesh-points in each parameter directions
    nmbel1[0] = nmbel1[1] = nmbel2[0] = nmbel2[1] = 1;
    for (int ki=0; ki<numpar1; ki++) {
	nmbel1[ki] = obj_int_[0]->getMeshSize(ki);
	ta[ki] = obj_int_[0]->startParam(ki);
	tb[ki] = obj_int_[0]->endParam(ki);
    }
    mesh1 = obj_int_[0]->getMesh();
    for (int ki=0; ki<numpar2; ki++) {
	nmbel2[ki] = obj_int_[1]->getMeshSize(ki);
	ta[numpar1+ki] = obj_int_[1]->startParam(ki);
	tb[numpar1+ki] = obj_int_[1]->endParam(ki);
    }
    mesh2 = obj_int_[1]->getMesh();

    vector<double>::iterator pt1 = mesh1; 
    vector<double>::iterator  pt2 = mesh2;

    // Find the mesh points with smallest distance
    double dist;
    int min_idx1[2], min_idx2[2];
    min_idx1[0] = min_idx1[1] = min_idx2[0] = min_idx2[1] = 0;
    double mindist = Utils::distance_squared(&pt1[0], &pt1[0] + dim, &pt2[0]);
    for (int kj=0; kj<nmbel1[1]; kj++) {
	for (int ki=0; ki<nmbel1[0]; ki++) {
	    for (int kr=0; kr<nmbel2[1]; kr++) {
		pt2 = mesh2;
		for (int kh=0; kh<nmbel2[0]; kh++) {
		    dist = Utils::distance_squared(&pt1[0], &pt1[0] + dim, &pt2[0]);
		    if (dist < mindist) {
			mindist = dist;
			min_idx1[0] = ki;
			min_idx1[1] = kj;
			min_idx2[0] = kh;
			min_idx2[1] = kr;
		    }
		    pt2 += dim;
		}
	    }
	    pt1 += dim;
	}
    }
  
    // Fetch parameter values associated with the minimum distance
    // mesh points
    int ind = 0;
    for (int ki=0; ki<numpar1; ki++) {
	seed[ind++] = obj_int_[0]->paramFromMesh(ki, min_idx1[ki]);
// 	seed[kr++] = 0.5 * (obj_int_[0]->startParam(ki)
// 			    + obj_int_[0]->endParam(ki)); 
    }

    for (int ki=0; ki<numpar2; ki++){
	seed[ind++] = obj_int_[1]->paramFromMesh(ki, min_idx2[ki]);
// 	seed[kr++] = 0.5 * (obj_int_[1]->startParam(ki)
// 			    + obj_int_[1]->endParam(ki));
    }

    // Don't want the seed in one endpoint
    double ptol = epsge_->getRelParRes();
    for (int kr=0; kr<numpar1+numpar2; kr++) {
	if (seed[kr]-ta[kr] < ptol) {
	    seed[kr] = ta[kr] + 0.25*(tb[kr]-ta[kr]);
	}
	if (tb[kr] - seed[kr] < ptol) {
	    seed[kr] = tb[kr] - 0.25*(tb[kr]-ta[kr]);
	}
    }
}


//===========================================================================
bool Intersector2Obj::atBoundary(const double *par1, vector<bool>& boundaries)
//===========================================================================
{
    //  Purpose: Check if an intersection point lies at any boundaries
    // The boundaries are arranged start parameter, end parameter for each 
    // parameter direction
    int total_num_par = obj_int_[0]->numParams() + obj_int_[1]->numParams();

    double rel_par_res = epsge_->getRelParRes();

    boundaries.resize(total_num_par*2);
    int ki;
    bool at_bd = false;
    for (ki=0; ki<total_num_par; ++ki)
    {
	boundaries[2*ki] = boundaries[2*ki+1]  = false;
	double min_par_val;
	double max_par_val;
	if (ki < obj_int_[0]->numParams()) {
	    min_par_val = obj_int_[0]->startParam(ki);
	    max_par_val = obj_int_[0]->endParam(ki);
	} else {
	    min_par_val = obj_int_[1]->startParam(ki - obj_int_[0]->numParams());
	    max_par_val = obj_int_[1]->endParam(ki - obj_int_[0]->numParams());
	}

	if (fabs(min_par_val - par1[ki]) < rel_par_res)
	{
	    boundaries[2*ki] = true;
	    at_bd = true;
	}
	    
	if (fabs(max_par_val - par1[ki]) < rel_par_res)
	{
	    boundaries[2*ki+1] = true;
	    at_bd = true;
	}
	    
    }

    return at_bd;
}

//===========================================================================
bool Intersector2Obj::atSameBoundary(const double *par1, const double *par2)
//===========================================================================
{
    //  Purpose: Check if two intersection points represented by their parameter
    //           values lies at a common object boundary
    int total_num_par = obj_int_[0]->numParams() + obj_int_[1]->numParams();

    double rel_par_res = epsge_->getRelParRes();

    bool at_bd1[8], at_bd2[8];
    int ki;
    for (ki=0; ki<total_num_par; ++ki)
    {
	at_bd1[2*ki] = at_bd1[2*ki+1] = at_bd2[2*ki] = at_bd2[2*ki+1] = false;
	double min_par_val;
	double max_par_val;
	if (ki < obj_int_[0]->numParams()) {
	    min_par_val = obj_int_[0]->startParam(ki);
	    max_par_val = obj_int_[0]->endParam(ki);
	} else {
	    min_par_val = obj_int_[1]->startParam(ki - obj_int_[0]->numParams());
	    max_par_val = obj_int_[1]->endParam(ki - obj_int_[0]->numParams());
	}

	if (fabs(min_par_val - par1[ki]) < rel_par_res)
	    at_bd1[2*ki] = true;
	    
	if (fabs(max_par_val - par1[ki]) < rel_par_res)
	    at_bd1[2*ki+1] = true;
	    
	if (fabs(min_par_val - par2[ki]) < rel_par_res)
	    at_bd2[2*ki] = true;
	    
	if (fabs(max_par_val - par2[ki]) < rel_par_res)
	    at_bd2[2*ki+1] = true;
    }

    for (ki=0; ki<2*total_num_par; ++ki)
	if (at_bd1[ki] && at_bd2[ki])
	    return true;

    return false;
}

//===========================================================================
void Intersector2Obj::printDebugInfo()
//===========================================================================
{
    int npar1 = obj_int_[0]->numParams();
    int npar2 = obj_int_[1]->numParams();
    int ki, kj;

    if (npar1 + npar2 == 2 && (npar1 == 2 || npar2 == 2))
    {
	// Point-surface
	if (getenv("SUBDIV_SFPT") && (*getenv("SUBDIV_SFPT")) == '1')
	{
	    double ta1[2], ta2[2], tb1[2], tb2[2];
	    std::cout << "================================================"
		      << std::endl;
	    std::cout << "Domain 1: ";
	    for (ki=0; ki<npar1; ki++) {
		ta1[ki] = obj_int_[0]->startParam(ki);
		tb1[ki] = obj_int_[0]->endParam(ki);
		std::cout << ta1[ki] << " ";
		std::cout << tb1[ki] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << "Domain 2: ";
	    for (ki=0; ki<npar2; ki++) {
		ta2[ki] = obj_int_[1]->startParam(ki);
		tb2[ki] = obj_int_[1]->endParam(ki);
		std::cout << ta2[ki] << " ";
		std::cout << tb2[ki] << " ";
	    }
	    std::cout << std::endl;

	    vector<shared_ptr<IntersectionPoint> > ipoint;
	    int_results_->getIntersectionPoints(ipoint);
	    std::cout << "Intersection points : " << ipoint.size()
		      << std::endl;
	    for (ki=0; ki < int(ipoint.size()); ki++) {
		for (kj=0; kj<npar1+npar2; kj++) {
		    std::cout << ipoint[ki]->getPar(kj) << "  ";
		}
		std::cout << ipoint[ki]->numNeighbours() << "  ";
		std::cout << ipoint[ki]->getDist();
		std::cout << std::endl;
	    }
	    int stop_break;
	    stop_break = 1;
	}
    }

    // Curve-surface
    if (npar1 + npar2 == 3) {
	if (getenv("SUBDIV_SFCV") && (*getenv("SUBDIV_SFCV")) == '1') {
	    double ta1[2], ta2[2], tb1[2], tb2[2];
	    std::cout << "================================================"
		      << std::endl;
	    std::cout << "Domain 1: ";
	    for (ki=0; ki<npar1; ki++) {
		ta1[ki] = obj_int_[0]->startParam(ki);
		tb1[ki] = obj_int_[0]->endParam(ki);
		std::cout << ta1[ki] << " ";
		std::cout << tb1[ki] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << "Domain 2: ";
	    for (ki=0; ki<npar2; ki++) {
		ta2[ki] = obj_int_[1]->startParam(ki);
		tb2[ki] = obj_int_[1]->endParam(ki);
		std::cout << ta2[ki] << " ";
		std::cout << tb2[ki] << " ";
	    }
	    std::cout << std::endl;

	    vector<shared_ptr<IntersectionPoint> > ipoint;
	    int_results_->getIntersectionPoints(ipoint);
	    std::cout << "Intersection points : " << ipoint.size()
		      << std::endl;
	    for (ki=0; ki < int(ipoint.size()); ki++) {
		for (kj=0; kj<npar1+npar2; kj++) {
		    std::cout << ipoint[ki]->getPar(kj) << "  ";
		}
		std::cout << ipoint[ki]->numNeighbours() << "  ";
		std::cout << ipoint[ki]->getDist();
		std::cout << std::endl;
	    }

	int stop_break;
	stop_break = 1;
	}
	if (prev_intersector_ == 0
	    || prev_intersector_->numParams() > numParams()) {
	    // Write objects to file
	    std::ofstream debug("geom_cv_sf_out.g2");
	    ParamCurveInt *curve1 = obj_int_[0]->getParamCurveInt();
	    ParamCurveInt *curve2 = obj_int_[1]->getParamCurveInt();
	    ParamSurfaceInt *surf1 = obj_int_[0]->getParamSurfaceInt();
	    ParamSurfaceInt *surf2 = obj_int_[1]->getParamSurfaceInt();
	    if (curve1) {
		SplineCurve *cv = curve1->getParamCurve()->geometryCurve();
		cv->writeStandardHeader(debug);
		cv->write(debug);
	    }
	    if (curve2) {
		SplineCurve *cv = curve2->getParamCurve()->geometryCurve();
		cv->writeStandardHeader(debug);
		cv->write(debug);
	    }
	    if (surf1) {
		shared_ptr<ParamSurface> srf = surf1->getParamSurface();
		srf->writeStandardHeader(debug);
		srf->write(debug);
	    }
	    if (surf2) {
		shared_ptr<ParamSurface> srf = surf2->getParamSurface();
		srf->writeStandardHeader(debug);
		srf->write(debug);
	    }
	}
    }
    
    // Curve-curve
    if (npar1 == 1 && npar2 == 1
	&& (prev_intersector_ == 0
	    || prev_intersector_->numParams() > numParams())) {
	std::ofstream debug("geom_cv_cv_out.g2");
	ParamCurveInt *curve1 = obj_int_[0]->getParamCurveInt();
	ParamCurveInt *curve2 = obj_int_[1]->getParamCurveInt();
	if (curve1) {
	    SplineCurve *cv = curve1->getParamCurve()->geometryCurve();
	    cv->writeStandardHeader(debug);
	    cv->write(debug);
	}
	if (curve2) {
	    SplineCurve *cv = curve2->getParamCurve()->geometryCurve();
	    cv->writeStandardHeader(debug);
	    cv->write(debug);
	}
    }

    // Surface-surface or surface-point (if write_point = true)
    bool write_point = (getenv("DEBUG_POINT")
			&& (*getenv("DEBUG_POINT"))=='1');
    if (npar1 + npar2 == 4
	|| (write_point && ((npar1 == 2 && npar2 == 0)
			    || (npar1 == 0 && npar2 == 2)))) {
	//bool write_sfs = false;
	bool write_sfs = true;
	int nmb_rec = nmbRecursions();
	if (nmb_rec == 1
	    || (prev_intersector_ != 0
		&& prev_intersector_->isSelfIntersection())) {
	    write_sfs = true;
	}
	if (write_sfs) {
	    std::ofstream debug("geom_sf_sf_out.g2");
	    ParamSurfaceInt *surf1 = obj_int_[0]->getParamSurfaceInt();
	    ParamSurfaceInt *surf2 = obj_int_[1]->getParamSurfaceInt();
	    if (surf1) {
		shared_ptr<ParamSurface> srf = surf1->getParamSurface();
		srf->writeStandardHeader(debug);
		srf->write(debug);
	    }
	    if (surf2) {
		shared_ptr<ParamSurface> srf = surf2->getParamSurface();
		srf->writeStandardHeader(debug);
		srf->write(debug);
	    }
	}

	double ta1[2], ta2[2], tb1[2], tb2[2];
	std::cout << "================================================"
		  << std::endl;
	std::cout << "Domain 1: ";
	for (ki=0; ki<npar1; ki++) {
	    ta1[ki] = obj_int_[0]->startParam(ki);
	    tb1[ki] = obj_int_[0]->endParam(ki);
	    std::cout << ta1[ki] << " ";
	    std::cout << tb1[ki] << " ";
	}
	std::cout << std::endl;
	std::cout << "Domain 2: ";
	for (ki=0; ki<npar2; ki++) {
	    ta2[ki] = obj_int_[1]->startParam(ki);
	    tb2[ki] = obj_int_[1]->endParam(ki);
	    std::cout << ta2[ki] << " ";
	    std::cout << tb2[ki] << " ";
	}
	std::cout << std::endl;

	vector<shared_ptr<IntersectionPoint> > ipoint;
	int_results_->getIntersectionPoints(ipoint);
	std::cout << "Intersection points in this pool: "
		  << ipoint.size() << std::endl;
	for (ki=0; ki < int(ipoint.size()); ki++) {
	    for (kj=0; kj<npar1+npar2; kj++) {
		std::cout << ipoint[ki]->getPar(kj) << "  ";
	    }
	    Point p1 = ipoint[ki]->getPoint1();
	    Point p2 = ipoint[ki]->getPoint1();
	    for (kj=0; kj<p1.dimension(); kj++)
		std::cout << p1[kj] << "  ";
	    for (kj=0; kj<p2.dimension(); kj++)
		std::cout << p2[kj] << "  ";
	    std::cout << ipoint[ki]->getDist() << "  ";
	    std::cout << ipoint[ki]->getSingularityType();
	    std::cout << std::endl;
	}

	if (getenv("DEBUG_PAR") && (*getenv("DEBUG_PAR"))=='1') {
	    int_results_->writeDebug();
	}
	int stop_break;
	stop_break = 1;
    }
	
}


//===========================================================================
