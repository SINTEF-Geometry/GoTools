//===========================================================================
//                                                                           
// File: Intersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: Intersector.C,v 1.59 2007-11-01 14:34:46 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Intersector.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/GeoTol.h"


using std::cout;
using std::endl;


namespace Go {


//===========================================================================
Intersector::Intersector(double epsge, Intersector* prev)
    : //int_results_(shared_ptr<IntersectionPool>(new IntersectionPool())),
      prev_intersector_(prev)
//===========================================================================
{
    epsge_ = shared_ptr<GeoTol>(new GeoTol(epsge));
}


//===========================================================================
Intersector::Intersector(shared_ptr<GeoTol> epsge, Intersector *prev)
    : //int_results_(shared_ptr<IntersectionPool>(new IntersectionPool())),
      prev_intersector_(prev)
//===========================================================================
{
    epsge_ = shared_ptr<GeoTol>(new GeoTol(epsge.get()));
}


//===========================================================================
void Intersector::compute(bool compute_at_boundary)
//===========================================================================
{
    // Purpose: Compute the topology of the current intersection

    // Make sure that no "dead intersection points" exist in the pool,
    // i.e. points that have been removed when compute() has been run
    // on sibling subintersectors.
    int_results_->synchronizePool();

    // Make sure that all intersection points at the
    // boundary/boundaries of the current object are already computed
    if (compute_at_boundary)
	getBoundaryIntersections();

    // Remove inner points in constant parameter intersection
    // links
    // (vsk, 0609) and isolated points identical to the existing ones
    int_results_->cleanUpPool();
    int nmb_orig = int_results_->numIntersectionPoints();

    if (getenv("DEBUG") && *(getenv("DEBUG")) == '1') {
	try {
	    printDebugInfo();
	} catch (...) {
	    MESSAGE("Failed printing debug info, continuing.");
	}
    }

    // Check if any intersections are possible in the inner of the
    // objects
    int status_intercept = performInterception();

    // Branch on the outcome of the interseption test
    if (status_intercept == 0) {
	// No intersection is possible
    } else if (status_intercept == 2) {
	// Both objects are too small for further processing.
	// Handle micro case
	microCase();
    } else if (degTriangleSimple()) {
	// This situation is currently relevant only for intersections
	// between two parametric surfaces. It will probably at some
	// stage be relevant for two-parametric functions.
	// All the necessary connections are made
    } else if (checkCoincidence()) {
	// The two objects coincide. The representation is already
	// updated according to this situation
    } else {
	// status_intercept == 1

	// Intersections might exist. Check for simple case. 0 = Maybe
	// simple case; 1 = Confirmed simple case.
	int status_simplecase = simpleCase();

	if (status_simplecase == 1) {
	    // Confirmed simple case.
	    // Compute intersection points or curves according to the
	    // properties of this particular intersection
	    updateIntersections();
	} else if (isLinear()) {
	    // Linearity is a simple case, but it is important to
	    // check for coincidence before trying to find/connect
	    // intersections as the simple case criteria is not
	    // satisfied
	    updateIntersections();
	}
	else if (complexIntercept())
	{
	    // Interception by more complex algorithms is performed
	    // (implicitization). No further intersectsions are found
	    // to be possible
	}
	else if (complexSimpleCase())
	{
	    // Simple case test by more complex algorithms is performed
	    // (implicitization). A simple case is found.
	    updateIntersections();
	} else if (!complexityReduced()) {
	    // For the time being, write documentation of the
	    // situation to a file
	    handleComplexity();
	} else {
	    // It is necessary to subdivide the current objects
	    doSubdivide();
	    
	    int nsubint = int(sub_intersectors_.size());
	    for (int ki = 0; ki < nsubint; ki++) {
		sub_intersectors_[ki]->getIntPool()
		    ->includeCoveredNeighbourPoints();
		sub_intersectors_[ki]->compute();
	    }
	}
    }

//     // Write intersection point diagnostics
//     if (numParams() == 4) {
// 	writeIntersectionPoints();
//     }

    if (prev_intersector_ && prev_intersector_->numParams() > numParams())
    {
	// Remove inner points in constant parameter intersection
	// links
	// (vsk, 0609) and isolated points identical to the existing ones
	int_results_->cleanUpPool(nmb_orig);

	// No more recursion at this level. Post iterate the intersection points
	doPostIterate();
    }

    // Prepare output intersection results
    if (prev_intersector_ == 0 || prev_intersector_->isSelfIntersection())
    {
	/*if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') {
	    cout << "Status after cleaning up pool:" << endl;
	    writeIntersectionPoints();
	    }*/

	// Remove loose ends of intersection links in the inner
	//int_results_->weedOutClutterPoints();
	if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') 
	{
	    cout << "Status after removing clutter points:" << endl;
	    writeIntersectionPoints();
	    int_results_->writeDebug();
	}

	// Remove loose ends of intersection links in the inner
	//int_results_->weedOutClutterPoints();
	int_results_->cleanUpPool(0);

	if (true /*getenv("DO_REPAIR") && *(getenv("DO_REPAIR")) == '1'*/) 
	{
	    if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') 
	    {
		cout << "Starting repair" << endl;
	    }
	    repairIntersections();

	    if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') 
	    {
		cout << "Status after repairing intersections:" << endl;
		writeIntersectionPoints();
	    }
	}
    }

    if (prev_intersector_ == 0) {
	// Top level intersector
	/*if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') {
	    cout << "Status after removing clutter points:" << endl;
	    writeIntersectionPoints();
	    }*/

// 	if (/*true */getenv("DO_REPAIR") && *(getenv("DO_REPAIR")) == '1') {
// 	    repairIntersections();

// 	    if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') {
// 		cout << "Status after repairing intersections:" << endl;
// 		writeIntersectionPoints();
// 	    }
// 	}

	if (getenv("DEBUG_FINISH") && *(getenv("DEBUG_FINISH")) == '1') {
	    int_results_->writeDebug();
	}

	int_results_->makeIntersectionCurves();
    }
}


//===========================================================================
void Intersector::
getResult(std::vector<shared_ptr<IntersectionPoint> >& int_points,
	  std::vector<shared_ptr<IntersectionCurve> >& int_curves)
//===========================================================================
{
    //  Purpose: Fetch intersection results
    
    if (int_results_.get()) {
	int_results_->getResult(int_points, int_curves);
	return;
    }
    MESSAGE("Warning: Null pointer to IntersectionPool in Intersector"
	    " detected in Intersector::getResult");
}


//===========================================================================
bool Intersector::validateSiblingPools()
//===========================================================================
{
    if (prev_intersector_ == 0) {
	return int_results_->validate();
    } else {
	int nsiblings = (int)prev_intersector_->sub_intersectors_.size();
	for (int i = 0; i < nsiblings; ++i) {
	    if (!prev_intersector_->sub_intersectors_[i]->int_results_->
		validate()) {
		return false;
	    }
	}
    }
    return true;

}


//===========================================================================
void Intersector::setHighPriSing(double* par)
//===========================================================================
{
    // Purpose: Instruct the intersector about known singular points

    int nmbpar = numParams();
    if (!hasSingularityInfo()) {
	if (prev_intersector_
	    && prev_intersector_->hasSingularityInfo()
	    && (prev_intersector_->numParams() == nmbpar)) {
	    singularity_info_ = (shared_ptr<SingularityInfo>)
		(new SingularityInfo(prev_intersector_
				     ->getSingularityInfo()));
	}
	else {
	    // Make empty singularity info instance
	    singularity_info_ = (shared_ptr<SingularityInfo>)
		(new SingularityInfo());
	}
    }
    singularity_info_->setHighPriSing(par, nmbpar);
}


//===========================================================================
void Intersector::writeIntersectionPoints() const
//===========================================================================
{
    cout << "****************************************" << endl;
    int_results_->writeIntersectionPoints();
    int_results_->writeIntersectionLinks();
    cout << "****************************************" << endl;

    return;
}


//===========================================================================


} // namespace Go
