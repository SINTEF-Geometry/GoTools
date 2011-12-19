//===========================================================================
//                                                                           
// File: Par2FuncIntersector.C                                               
//                                                                           
// Created: Mon Sep 27 14:49:22 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par2FuncIntersector.C,v 1.25 2006-11-03 14:43:22 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Par2FuncIntersector.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/closestPtCurves.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/PtPtIntersector.h"
#include "GoTools/intersections/Par1FuncIntersector.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/extremalPtSurfSurf.h"
#include "GoTools/intersections/IntersectionPool.h"


using std::vector;
using std::pair;
using std::make_pair;
using std::swap;


namespace Go
{


//===========================================================================
Par2FuncIntersector::Par2FuncIntersector(shared_ptr<ParamFunctionInt> func,
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
Par2FuncIntersector::~Par2FuncIntersector()
//===========================================================================
{
    // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
Par2FuncIntersector::lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
					   shared_ptr<ParamFunctionInt> obj2,
					   Intersector* prev,
					   int eliminated_parameter,
					   double eliminated_value)
//===========================================================================
{
    // Expecting objects to be of type Param1FunctionInt & Param0FunctionInt.
    shared_ptr<Intersector> curr_inter
	(new Par1FuncIntersector(obj1, obj2, epsge_, prev,
				 eliminated_parameter,
				 eliminated_value));

    // We should also include existing intersection points with input
    // parameter (eliminated_value).

    return curr_inter;
}


//===========================================================================
int Par2FuncIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between a point and a function can occur only if the
  // function is degenerate. So ...

  return 0;
}


//===========================================================================
//
// Purpose :
//
// Written by : Sverre Briseid, SINTEF, Sep 2004.
//===========================================================================
void Par2FuncIntersector::microCase()
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
	int_results_->getMidParameter(mid);

	shared_ptr<IntersectionPoint> tmp
	    = int_results_->addIntersectionPoint(func_int_, C_,
						 epsge_, mid, mid);
    }
  else if (ints.size() == 1)
    {
      // One intersection point. Leave it like this.
      ;
    }
    else if (ints.size() == 2)
    {
	// Make sure the intersection points are connected
	ints[0]->connectTo(ints[1], MICRO_PAR2FUNC);
    }
}


//===========================================================================
//
// Purpose : Given two parametric surfaces in a simple case situation, 
//           connect intersection points at the surface boundaries into
//           tracks to represent the topology of the intersection problem
//
// Written by : Vibeke Skytt, 1104
//
//===========================================================================
int
Par2FuncIntersector::updateIntersections()
//===========================================================================
{
//     // Fetch ParamSurface information
     Param2FunctionInt* func_int = func_int_->getParam2FunctionInt();
     Param0FunctionInt* C_int = C_->getParam0FunctionInt();
//     ASSERT(func_int != 0); // && surf2 != 0);

//     ParamSurface *srf1 = surf1->getParamSurface().get();
//     ParamSurface *srf2 = surf2->getParamSurface().get();

//     // Simple case is detected by checking overlap in the normal cones of the
//     // surface. Sort the intersection points at the boundaries of this surface
//     // with respect to the  normal vector of the plane spanned by the
//     // cone axis.
//     DirectionCone cone1 = obj_int_[0]->directionCone();
//     DirectionCone cone2 = obj_int_[1]->directionCone();
//     Point axis1 = cone1.centre();
//     Point axis2 = cone2.centre();
//     Point vec = axis1 % axis2;

     shared_ptr<ParamSurface> func_sf = func_int->getSurface();
     double C = C_int->getValue();

     // Fetch the intersection points at the boundary
     // sorted?
     // @@sbr We must somehow sort the intersection points along the bd.
     // It makes sense to use the parametric intersection directions.
     Param2FunctionInt* func2_int = func_int_->getParam2FunctionInt();
     Point mon_dir;
     func2_int->monotone(mon_dir);
     Point sorting_dir = mon_dir; //(2);
//      sorting_dir[0] = -mon_dir[1];
//      sorting_dir[1] = mon_dir[0];
     sorting_dir.normalize();
     vector<shared_ptr<IntersectionPoint> > bd_ints;
     int_results_->getSortedBdInts(sorting_dir, bd_ints);

    // Classify the intersection points
    int ki;
    int nmb_nottouch = 0, nmb_ints;
    int nmb_in = 0, nmb_out = 0, nmb_perpendicular = 0;
    int nmb_branch;
    IntPtClassification type, type2;
    vector<pair<shared_ptr<IntersectionPoint>,IntPtClassification> > bd_ints_typed;
    nmb_ints = (int)bd_ints.size();
    if (nmb_ints == 0)
	return 1;  // No points to connect. Finished.

    for (ki=0; ki<nmb_ints; ki++)
    {
	type = bd_ints[ki]->getClassification(func_sf.get(), C);

	// Check if another intersection branch exist
	nmb_branch = bd_ints[ki]->numBranches();
	if (nmb_branch == 2)
	{
	    MESSAGE("Currently implementing support for multiple branches!");
	    type2 = bd_ints[ki]->getClassification(func_sf.get(), 1);
	    if (type == DIR_TOUCH)
		type = type2;
	    else if (type2 != DIR_TOUCH)
	    {
		bd_ints_typed.push_back(make_pair(bd_ints[ki],type2));
		if (type2 != DIR_TOUCH)
		    nmb_nottouch++;
		if (type2 == DIR_IN)
		    nmb_in++;
		if (type2 == DIR_OUT)
		    nmb_out++;
	    }
	}
	bd_ints_typed.push_back(make_pair(bd_ints[ki],type));
	if (type != DIR_TOUCH)
	    nmb_nottouch++;
	if (type == DIR_IN)
	    nmb_in++;
	if (type == DIR_OUT)
	    nmb_out++;
	if (type == DIR_PERPENDICULAR)
	    nmb_perpendicular++;
    }

    // Move corner points to the end of the array
    if (nmb_nottouch < int(bd_ints_typed.size()))
    {
	for (ki=0; ki<nmb_nottouch; ki++)
	    if (bd_ints_typed[ki].second == DIR_TOUCH)
	    {
		pair<shared_ptr<IntersectionPoint>, IntPtClassification> curr
		    = bd_ints_typed[ki];
		bd_ints_typed.erase(bd_ints_typed.begin()+ki, 
				    bd_ints_typed.begin()+ki+1);
		bd_ints_typed.push_back(curr);
		ki--;
	    }
    }

    // Perform connection of intersection points into track according to the
    // current situation
    // First check if more than one relevant point exist
    if (nmb_ints == 1 || 
	(nmb_nottouch == 1 && 
	 !(nmb_ints == 2 && (nmb_in == 1 || nmb_out == 1 || nmb_perpendicular == 1))))	
	return 1;   // Finished, not necessary to connect

    // Then check if the intersection points are connected already
    if (isConnected(bd_ints_typed, nmb_nottouch))
	return 1;   // Nothing more to do

   // Try to make connections between the intersection points according
    // to their sequence. The classification must be consistent with the
    // sequence
    // @@sbr As a quick fix we move the first in pt to the end if it is DIR_OUT.
    size_t nmb_moved = 0;
    while (nmb_moved < bd_ints_typed.size()) {
	if (bd_ints_typed[0].second == DIR_OUT) {
	    bd_ints_typed.push_back(bd_ints_typed[0]);
	    bd_ints_typed.erase(bd_ints_typed.begin());
	    ++nmb_moved;
	} else {
	    break;
	}
    }
    if (connectDirected(bd_ints_typed, nmb_nottouch))
	return 1;    // Connections performed

    // If only two relevant (not corner touching) points exist and one points in     
    // or one points out, the two points must be connected. A corner touch
    // is overruled if one point points in or out
    if (nmb_ints == 2 &&
	((nmb_in + nmb_out == 2) || nmb_in == 1 || nmb_out == 1))
    {
	// Connect the two non-touching points
	bd_ints_typed[0].first->connectTo(bd_ints_typed[1].first,
					  SIMPLE_TWO_POINTS); //nmb_nottouch].first);
	return 1;
    }

    if (nmb_ints == 2 && 
	((nmb_in + nmb_out  + nmb_perpendicular == 2) || nmb_in == 1 || nmb_out == 1
	    || nmb_perpendicular == 1))
    {
	// Check if a connection between the two non-touch points is correct
	if (canConnect(bd_ints_typed[0].first, bd_ints_typed[1].first)) //nmb_nottouch].first))
	{
	    bd_ints_typed[0].first->connectTo(bd_ints_typed[1].first,
					      SIMPLE_TWO_POINTS); //nmb_nottouch].first);
	}
	return 1;
    }



    // This situation is not handled. Write debug output.
    // This situation is not handled. Write debug output.
    if (getenv("DEBUG") && (*getenv("DEBUG"))=='1')
    {
	std::cout << "Inconsistent conditions for connections found " << std::endl;
	writeDebugConnect(bd_ints_typed);
    }
    return 0;
}


//===========================================================================
bool Par2FuncIntersector::
isConnected(vector<shared_ptr<IntersectionPoint> > bd_ints,
	    int nmbbd)
//===========================================================================
{
    // Purpose : Check if all the points in a set of intersection
    // points at the surface boundaries is connected to some other
    // point in the set

    vector<bool> found_connection(nmbbd, false);
    int ki, kj;
    for (ki=0; ki<nmbbd; ki++)
    {
	for (kj=ki+1; kj<nmbbd; kj++)
	{
	    if (found_connection[ki] && found_connection[kj])
		continue;
	    if (bd_ints[ki]->isConnectedTo(bd_ints[kj]))
		found_connection[ki] = found_connection[kj] = true;
	}
    }
    for (ki=0; ki<nmbbd; ki++)
	if (!found_connection[ki])
	    return false;

    return true;
}


//===========================================================================
bool Par2FuncIntersector::
isConnected(vector<pair<shared_ptr<IntersectionPoint>,
	    IntPtClassification> >& bd_ints,
	    int nmbbd)
//===========================================================================
{
    // Purpose : Check if all the points in a set of intersection
    // points at the surface boundaries is connected to some other
    // point in the set

    vector<bool> found_connection(nmbbd, false);
    int ki, kj;
    for (ki=0; ki<nmbbd; ki++)
    {
	for (kj=ki+1; kj<nmbbd; kj++)
	{
	    if (found_connection[ki] && found_connection[kj])
		continue;
	    if (bd_ints[ki].first == bd_ints[kj].first)
		continue;  // The same branch point

	    if (bd_ints[ki].first->isConnectedTo(bd_ints[kj].first))
		found_connection[ki] = found_connection[kj] = true;
	}
    }
    for (ki=0; ki<nmbbd; ki++)
	if (!found_connection[ki])
	    return false;

    return true;
}


//===========================================================================
bool Par2FuncIntersector::
connectDirected(vector<pair<shared_ptr<IntersectionPoint>,
		IntPtClassification> >& bd_ints,
		int nmbbd)
//===========================================================================
{
    // Purpose : Connect boundary intersection in correct sequence if
    // their type makes it possible

    if (nmbbd == 1)
	return false;  // No rule

    // @@@ VSK. Should we test on distance or influence area. Connect anyway
    // if the points are close enough
    // First test if connection is possible
    int ki;
    bool has_connected = false;
    for (ki=1; ki<nmbbd; ki++)
    {
	// If the same branch point
	if (bd_ints[ki-1].first == bd_ints[ki].first);
	// If already connected pass to the next point
	else if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first));
	else if (bd_ints[ki-1].second == DIR_OUT)
	    return false;   // The intersection curve cannot start outwards
	else if (bd_ints[ki-1].second == DIR_IN &&
		 bd_ints[ki].second == DIR_IN)
	    return false;   // Inconsistent directions
	else if ((bd_ints[ki-1].second == DIR_IN || bd_ints[ki-1].second == DIR_PERPENDICULAR) &&
		 (bd_ints[ki].second == DIR_OUT || bd_ints[ki].second == DIR_PERPENDICULAR))
	    ki++;  // OK
	else if ((bd_ints[ki-1].second == DIR_IN || bd_ints[ki-1].second == DIR_PERPENDICULAR) &&
		 bd_ints[ki].second == DIR_PARALLEL);  // OK
	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
		 (bd_ints[ki].second == DIR_OUT || bd_ints[ki].second == DIR_PERPENDICULAR)); 
	else
	    return false;  // No rule
    }
    if (ki == nmbbd-1 &&
	(bd_ints[ki].second == DIR_IN || bd_ints[ki].second == DIR_OUT ||
	    bd_ints[ki].second == DIR_PERPENDICULAR))
	return false;    // Hanging boundary point not possible to connect

    // Perform connections
    for (ki=1; ki<nmbbd; ki++)
    {
	// If the same branch point
	if (bd_ints[ki-1].first == bd_ints[ki].first);
	// If already connected pass to the next point
	else if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first));
	else if ((bd_ints[ki-1].second == DIR_IN  || bd_ints[ki-1].second == DIR_PERPENDICULAR) &&
		 (bd_ints[ki].second == DIR_OUT || bd_ints[ki].second == DIR_PERPENDICULAR))
	{
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
	    has_connected = true;
	    ki++; 
	}
	else if ((bd_ints[ki-1].second == DIR_IN || bd_ints[ki-1].second == DIR_PERPENDICULAR) &&
		 bd_ints[ki].second == DIR_PARALLEL)
	{
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
	    has_connected = true;
	}
	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
		 (bd_ints[ki].second == DIR_OUT || bd_ints[ki].second == DIR_PERPENDICULAR))
	{
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
	    has_connected = true;
	}
    }
    return has_connected;
}
//     int ki;
//     for (ki=1; ki<nmbbd; ki++)
//     {
// 	// If already connected pass to the next point
// 	if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first));
// 	else if (bd_ints[ki-1].second == DIR_OUT)
// 	    return false;   // The intersection curve cannot start outwards
// 	else if (bd_ints[ki-1].second == DIR_IN &&
// 		 bd_ints[ki].second == DIR_IN)
// 	    return false;   // Inconsistent directions
// 	else if (bd_ints[ki-1].second == DIR_IN &&
// 		 bd_ints[ki].second == DIR_OUT)
// 	    ki++;  // OK
// 	else if (bd_ints[ki-1].second == DIR_IN &&
// 		 bd_ints[ki].second == DIR_PARALLEL);  // OK
// 	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
// 		 bd_ints[ki].second == DIR_OUT);   // OK
// 	else
// 	    return false;  // No rule
//     }
//     if (ki == nmbbd-1 &&
// 	(bd_ints[ki].second == DIR_IN || bd_ints[ki].second == DIR_OUT))
// 	return false;    // Hanging boundary point not possible to connect

//     // Perform connections
//     for (ki=1; ki<nmbbd; ki++)
//     {
// 	// If already connected pass to the next point
// 	if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first));
// 	else if (bd_ints[ki-1].second == DIR_IN &&
// 		 bd_ints[ki].second == DIR_OUT)
// 	{
// 	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
// 	    ki++; 
// 	}
// 	else if (bd_ints[ki-1].second == DIR_IN &&
// 		 bd_ints[ki].second == DIR_PARALLEL)
// 	{
// 	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
// 	}
// 	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
// 		 bd_ints[ki].second == DIR_OUT)
// 	{
// 	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first, LINK_UNDEFINED);
// 	}
//     }
//     return true;
// }


//===========================================================================
bool
Par2FuncIntersector::canConnect(shared_ptr<IntersectionPoint> pt1,
				shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    // Purpose : Check if a connection between two intersection points
    // at the surface boundaries is possible

    // Try to iterate to an intersection point at a constant parameter curve
    // between the two points. First define the curve.
    vector<double> par1 = pt1->getPar();
    vector<double> par2 = pt2->getPar();
    double maxdist = 0.0, pdist;
    int max_ind=-1, ki;
    int size = (int)(par1.size());
    for (ki=0; ki<size; ki++)
    {
	pdist = fabs(par1[ki]-par2[ki]);
	if (pdist > maxdist)
	{
	    maxdist = pdist;
	    max_ind = ki;
	}
    }
    if (max_ind < 0)
	return false;   // Do not connect

    Param2FunctionInt *func_int = func_int_->getParam2FunctionInt();
    Param0FunctionInt *C_int = C_->getParam0FunctionInt();
    shared_ptr<ParamCurve> crv;
    shared_ptr<ParamSurface> srf;
    RectDomain domain;  
//     int cv_idx, sf_idx;
//    if (max_ind < 2)
//     {
	crv = func_int->getConstantParameterCurve(max_ind,
						  0.5*(par1[max_ind]+par2[max_ind]));
// 	cv_idx = 0;
// 	sf_idx = 1;
//     }
//    else
//     {
// 	crv = func_int->getConstantParameterCurve(max_ind-2,
// 					       0.5*(par1[max_ind]+par2[max_ind]));
// // 	srf = surf1->getParentParamSurface(domain);
// // 	cv_idx = 2;
// // 	sf_idx = 0;
//     }
	
   // Make seed
   int other_ind = (max_ind == 0) ? 1 : 0;
   double seed = 0.5*(par1[other_ind]+par2[other_ind]);

   // Iterate to a closest point
   double clo_par, clo_dist;
   Point c_pt(1);
   c_pt[0] = C_int->getValue();
   Point clo_pt;
   crv->closestPoint(c_pt, crv->startparam(), crv->endparam(), clo_par, clo_pt, clo_dist, &seed);
   if (clo_dist < epsge_->getEpsge())
       return true;  // An intersection point is found
   else
       return false;
}

//===========================================================================
//
// Purpose : Given a 2-parametric function and no simple case situation, 
//           subdivide the function to produce more simple sub problems.
//           Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004
//
//===========================================================================
int Par2FuncIntersector::doSubdivide()
//===========================================================================
{
    int ki;

    // Intersection objects belonging to the next recursion level

    // Sort the curves according to importance of subdivision according to
    // properties of the curves and already computed intersections
    double deg_tol = epsge_->getEpsge();
    func_int_->isDegenerate(deg_tol, 0, NULL); // @@sbr Just to set deg_tol ...
    bool can_subdiv0 = func_int_->canDivide(0);
    bool can_subdiv1 = func_int_->canDivide(1);
    if (!can_subdiv0 && !can_subdiv1)
 	return 0;   // Not possible to subdivide any of the curves

    if (getenv("SUBDIV_FUNC2") && *(getenv("SUBDIV_FUNC2"))=='1') {
	std::cout << "================================================" << std::endl;
	std::cout << "Domain 1: ";
	for (ki=0; ki<2; ki++) {
	    std::cout << func_int_->startParam(ki) << " ";
	    std::cout << func_int_->endParam(ki) << " ";
	}
	std::cout << std::endl;
    }

    int perm[2]; // perm[0] currently the most suitable direction for subdivision.
    //     int deg_edge[4];
    int nmb_subdiv = sortParameterDirections(perm); //, deg_edge);
    if (nmb_subdiv == 0)
	return 0;   // Not possible to subdivide neither the curve nor
                    // the surface


    vector<shared_ptr<ParamFunctionInt> > sub_functions;
    sub_functions.push_back(func_int_);
    //int nbobj = 1;

    // For each parameter direction in prioritized order, fetch an
    // appropriate subdivision parameter, and perform subdivision
    int kj;
    double subdiv_par;
    vector<shared_ptr<ParamFunctionInt> > func_sub;
    vector<shared_ptr<ParamFunctionInt> > subdiv_func;
    int found;

    for (ki=0; ki<nmb_subdiv; ki++)
      {
	int pardir = perm[ki];
	found = getSubdivisionParameter(perm[ki], subdiv_par);
	if (!found) {
	    // No parameter value is found. Move this parameter direction
	    // to the end of the permutation array, and decrease the
	    // number of parameter directions where it is possible to 
	    // subdivide.
	    swap(perm[0], perm[1]);
// 	    found = getSubdivisionParameter(perm[ki], subdiv_par);
// 	    if (!found) {
// 		// We're done, nothing more to do.
// 		return 0;
// 	    }
	    nmb_subdiv--;
	    ki--;
	    continue;
	}

	if (getenv("SUBDIV_FUNC2") && (*getenv("SUBDIV_FUNC2")) == '1')
	    {
		std::cout << "Subdivide dir = " << perm[ki] << " par = " << subdiv_par;
		std::cout << " criterium = " << found << std::endl;
	    }

	// Subdivide the current object (or objects. If the current object
	// is the surface then it can have been subdivided in one parameter
	// direction already and we have two sub surfaces to subdivide)
// 	int idx = 0; //ki; //(subdiv_func.size() == 1) ? ki : 
	for (kj = 0; kj < int(sub_functions.size()); ++kj) {
	    func_sub.clear();
	    subdiv_func.clear();
	    // Subdivide the current surf
	    try {
		sub_functions[kj]->subdivide(pardir, subdiv_par, func_sub, subdiv_func);
	    } catch (...) {
		func_sub.clear();
		subdiv_func.clear();
	    }
	    if (func_sub.size() <= 1 || subdiv_func.size() == 0)
		continue;  // No new objects (both questions are really the
	    // same test)
	    for (size_t kr=0; kr<func_sub.size(); kr++)
		// We're interested in the func_int_, not the temp obj's.
		func_sub[kr]->setParent(func_int_.get()); //sub_functions[kj].get());
      
	    for (size_t kr = 0; kr < subdiv_func.size(); kr++) {
		// Intersect the subdivision object with the other object
		// @@@ VSK Faktorenes orden er ikke likegyldig!!
		shared_ptr<Intersector> subdiv_intersector = 
		    lowerOrderIntersector(subdiv_func[kr], C_,
					  this, perm[ki], subdiv_par);

		// Is it here relevant to fetch existing intersection points
		// and/or insert points into intersection curves before computing
		// intersection points with the subdivision points. Normally, this
		// is mostly of interest when surfaces are involved, but intersection
		// intervals might exist?
		// These computations do anyway involve the intersection pool, but
		// they might be trigged somehow. The parameter direction and value
		// are required information

		// @@sbr But for the case with intersection between plane and tube
		// this leads to a problem as only half of the intersection points
		// along the tangential intersection are picked up ...
		subdiv_intersector->compute(false);

		// Check quality of intersection points
		// If the quality is not sufficient, find a new subdivision parameter
		// and repeat the process.
		// Otherwise
		int_results_->includeReducedInts(subdiv_intersector->getIntPool());
	    }
	    sub_functions[kj] = func_sub[0];
	    sub_functions.insert(sub_functions.begin()+kj+1, func_sub.begin()+1,
				 func_sub.end());
// 	    idx += func_sub.size() - 1;
	    kj += (int)func_sub.size() - 1;
	    // 	  sub_functions.insert(sub_functions.end(), func_sub.begin(), func_sub.end());
	}
    }

    for (ki = 0; ki < int(sub_functions.size()); ++ki) {
	shared_ptr<Intersector> intersector = 
	    shared_ptr<Intersector>(new Par2FuncIntersector(sub_functions[ki],
							    C_, epsge_, this));
	// 	intersector->getIntPool()->setPoolInfo(int_results_);
	sub_intersectors_.push_back(intersector);
    }

    return 1;
}


//===========================================================================
//
// Purpose : Fetch information related to the 2par-function in order to decide
//           if we can subdivide and in which direction.
//
// Written by : Vibeke Skytt, 0804
//
//===========================================================================
int Par2FuncIntersector::sortParameterDirections(int perm[]) //, int deg_edge[])
//===========================================================================
{
  double length[2], wiggle[2];  // Estimated length and wiggliness
  bool inner_knots[2], critical_val[2], can_divide[2], is_deg[2];
  bool has_inner_ints[2];
  double rel_wiggle = 0.1;
  double rel_length = 0.1;
  double deg_tol = epsge_->getEpsge();

  // Fetch information from the sfs.
  // @@sbr Hmm, the wiggle for a 1d-surface ... Makes sense if looked upon as (u, v, f) in R^3.
  func_int_->getLengthAndWiggle(length, wiggle);

  double med_len = 0.0, med_wiggle = 0.0;
  for (int ki=0; ki<2; ki++)
  {
      med_len += length[ki];
      med_wiggle += wiggle[ki];
  }
  med_len /= (double)2;
  med_wiggle /= (double)2;

  double min_length = std::max(0.5*med_len, epsge_->getEpsge());
  double min_wiggle = std::max(0.5*med_wiggle, 0.02);

  for (int ki = 0; ki < 2; ++ki) {
      inner_knots[ki] = func_int_->hasInnerKnots(ki);

      critical_val[ki] = func_int_->hasCriticalVals(ki);

      can_divide[ki] = func_int_->canDivide(ki);

      is_deg[ki] = func_int_->isDegenerate(deg_tol, ki, NULL); // Corr to flatness for a function.
      // pardir == 1 equals v-dir?
      has_inner_ints[ki] = int_results_->hasPointsInInner(ki);
  }

  // Fetch information from the intersection pool according to inner 
  // intersection points
  // Number of parameter directions is two.
  int size = 2;
  int curr = 0;

  int min_nmb = 0;

  // Initiate permutation array
  for (int ki=0; ki<size; ki++)
    perm[ki] = ki;

  // Sort according to the given values
  // First sort out the directions where subdivsion is impossible
  for (int ki=0; ki<size; ki++)
    {
      if (!can_divide[perm[ki]])
	{
	  if (perm[ki] < size-1)
	    std::swap(perm[ki], perm[size-1]);
	  size--;
	}
    }

  // First sort out the directions where subdivsion is impossible
  for (int ki=0; ki<size; ki++)
    {
      if (is_deg[perm[ki]])
	{
	  if (perm[ki] < size-1)
	    std::swap(perm[ki], perm[size-1]);
	  size--;
	}
    }

  // First priority is parameter directions with critical values
  for (int ki=curr; ki<size; ki++)
    for (int kj=ki+1; kj<size; kj++)
      if (critical_val[perm[kj]] && !critical_val[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	  curr++;
	  break;
	}

  // Next criterium is inner knots
  for (int ki=curr; ki<size; ki++)
    for (int kj=ki+1; kj<size; kj++)
      if (inner_knots[perm[kj]] && !inner_knots[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	  curr++;
	  break;
	}

  // Existing intersection points internal to the parameter interval
  for (int ki=curr; ki<size; ki++)
    for (int kj=ki+1; kj<size; kj++)
      if (has_inner_ints[perm[kj]] && !has_inner_ints[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	  curr++;
	  break;
	}

  // Finally length and curvature
  min_nmb = curr;
  for (int ki=curr; ki<size; ki++)
    for (int kj=ki+1; kj<size; kj++)
      if ((length[perm[kj]] > length[perm[ki]] && 
	   wiggle[perm[kj]] > wiggle[perm[ki]]) ||
	  (length[perm[kj]] > length[perm[ki]] && 
	   wiggle[perm[kj]] > rel_wiggle*wiggle[perm[ki]]) ||
	  (wiggle[perm[kj]] > wiggle[perm[ki]] && 
	   length[perm[kj]] > rel_length*length[perm[ki]]))
	{
	  std::swap(perm[ki], perm[kj]);
	  curr++;
	  break;
	}

  // Check that the minimum conditions for subdivision is satisfied
  for (int ki=size-1; ki>=min_nmb; ki--, size--)
      if (length[perm[ki]] >= min_length || wiggle[perm[ki]] >= min_wiggle)
	  break;

  return size;
}


//===========================================================================
//
// Purpose : Return (what seems to be) the most suitable subdivision param.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004
//
//===========================================================================
int Par2FuncIntersector::getSubdivisionParameter(int dir, double& par)
//===========================================================================
{
  int ki, sgn;
  int nmbdir1 = func_int_->numParams();
  //int nmbdir2 = C_->numParams();

  // Set pointer to the intersection object corresponding to the parameter
  // direction
  ParamFunctionInt *obj = func_int_.get();
  //Param2FunctionInt *func = obj->getParam2FunctionInt();

  int pdir = (dir < nmbdir1) ? dir : dir-nmbdir1;
  double ta = obj->startParam(pdir);
  double tb = obj->endParam(pdir);
  double frac = 0.1;

  // Get critical parameters
  // @@@ Critical parameters of different priority? Sorting?
  vector<double> critical_pars = obj->getCriticalVals(pdir);

  int size = (int)critical_pars.size();
  int is_critical = 0;
  if (size > 0)
    {
      // Check suitability of the critical parameters
      for (ki=0; ki<size; ki++)
	{
	    vector<shared_ptr<IntersectionPoint> > int_pts;
	  is_critical = int_results_->inInfluenceArea(dir, critical_pars[ki], 
						      int_pts);
	  if (int_pts.size() > 0)
	  {
	      // Check if the subdivision position is OK with respect to the
	      // angle between the subdivision curve and the intersection curve
	      is_critical = checkSubdivParam(dir, critical_pars[ki], ta, tb, int_pts);
	  }
	  if (is_critical == 0 || is_critical == 2)
	    {
	      par = critical_pars[ki];
// 	      splitIntResults(int_pts, dir, par, ta, tb);
	      return 2;
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
	  vector<shared_ptr<IntersectionPoint> > int_pts;
	  is_critical = int_results_->inInfluenceArea(dir, knot_vals[ki], 
						      int_pts);
	  if (int_pts.size() > 0)
	  {
	      // Check if the subdivision position is OK with respect to the
	      // angle between the subdivision curve and the intersection curve
	      is_critical = checkSubdivParam(dir, knot_vals[ki], ta, tb, int_pts);
	  }
	  if (is_critical == 0 || is_critical == 2)
	    {
	      par = knot_vals[ki];
// 	      splitIntResults(int_pts, dir, par, ta, tb);
	      return 2;
	    }
// //	  is_critical = int_results_->inInfluenceArea(dir, knot_vals[ki]);
// 	  if (is_critical == 0 || is_critical == 2)
// 	    {
// 	      par = knot_vals[ki];
// 	      return true;
// 	    }
	}
    }

   // Check for intersection points internal to a boundary in which to subdivide    
    vector<shared_ptr<IntersectionPoint > > ip;
    int_results_->getIntersectionPoints(ip);
    size = (int)ip.size();
    if (size > 0)
      {
      // Check suitability of the intersection points
//        for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
//  	   ki+=sgn*kj, kj++, sgn*=-1)
        for (ki=0; ki<size; ki++)
	  {
	    par = ip[ki]->getPar(dir);
//	    par = inner_ints[ki];
	    if (par < ta+frac*(tb-ta) ||
		par > tb-frac*(tb-ta))
		continue;

	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = int_results_->inInfluenceArea(dir, par,
						      int_pts);
	    if (int_pts.size() > 0)
	      {
		// Check if the subdivision position is OK with respect to the
		// angle between the subdivision curve and the intersection curve
		is_critical = checkSubdivParam(dir, par, ta, tb, int_pts);
	      }
	    if (is_critical == 0 || is_critical == 2)
	      {
// 	        splitIntResults(int_pts, dir, par, ta, tb);
		return 5;
	    }
	}
    }
//    // Check for inner intersection points in which to subdivide
//    vector<double> inner_ints = int_results_->getSortedInnerInts(dir);
//    size = inner_ints.size();
//   if (size > 0)
//     {
//       // Check suitability of the intersection points
//       for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
// 	   ki+=sgn*kj, kj++, sgn*=-1)
// 	{
// //	  is_critical = int_results_->inInfluenceArea(dir, inner_ints[ki]);
// 	  if (is_critical == 0 || is_critical == 2)
// 	    {
// 	      par = inner_ints[ki];
// 	      return true;
// 	    }
// 	}
//     }

  
  // Subdivide at a suitable parameter
  double divpar = 0.5*(ta+tb);
  double del = 0.1*(tb-ta);
  double tint = del;
  for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=del, sgn*=-1)
    {
      vector<shared_ptr<IntersectionPoint> > int_pts;
      is_critical = int_results_->inInfluenceArea(dir, divpar, int_pts);
      if (int_pts.size() > 0)
      {
	  // Check if the subdivision position is OK with respect to the
	  // angle between the subdivision curve and the intersection curve
	  is_critical = checkSubdivParam(dir, divpar, ta, tb, int_pts);
      }
      if (is_critical == 0 || is_critical == 2)
	{
	  par = divpar;
// 	  splitIntResults(int_pts, dir, par, ta, tb);
	  return 6;
// 	  par = divpar;
// 	  return true;
	}
    }

  return 0;  // No subdivision parameter found
}

//===========================================================================
//
// Purpose :
//
// Written by : Vibeke Skytt, 1104
//
//===========================================================================
int
Par2FuncIntersector::checkSubdivParam(int dir, double par, double ta, double tb,
				  vector<shared_ptr<IntersectionPoint> >& int_pts)
//
// Return value : = 1 : Not a good parameter
//                = 2 : Parameter accepted
//===========================================================================
{
    int pdir = (dir < 2) ? dir : dir-2;
    Param2FunctionInt* func_int = func_int_->getParam2FunctionInt();
    int ki;
    double pval;
    double frac = 0.2;
    double ang_tol = std::min(10.0*epsge_->getRefAng(), 0.1);
    double eps = epsge_->getEpsge();

    // Skip the intersection points at a safe distance from the subdivision value
    int npts = (int)int_pts.size();
    for (ki=0; ki<npts; ki++)
    {
	pval = int_pts[ki]->getPar(dir);
	if (fabs(par - pval) >= frac*(tb - ta) &&
	    !(int_pts[ki]->inInfluenceAreaBracket(dir, pval)))
	{
	    std::swap(int_pts[ki], int_pts[npts-1]);
	    npts--;
	    ki--;
	}
   }

//     if (npts > 2)
// 	return 1;   // Not a good subdivision parameter

    if (npts > 1)
    {
	// Fetch the first and the last point
	int idx = (dir >= 2)*2 + (1-pdir);
	vector<shared_ptr<IntersectionPoint> > pts(2);
	pts[0] = pts[1] = int_pts[0];
	for (ki=1; ki<npts; ki++)
	{
	    if (int_pts[ki]->getPar(idx) < pts[0]->getPar(idx))
		pts[0] = int_pts[ki];
	    if (int_pts[ki]->getPar(idx) > pts[1]->getPar(idx))
		pts[1] = int_pts[ki];
	}

	// Check if the complete subdivision curve is in intersection.
	// First check if the two intersection points lie at the same boundary.
	// In that case, the subdivision parameter cannot be accepted
// 	if (int_results_->atSameBoundary(pts[0], pts[1]))
// 	    return 1;

// 	// Check the sequence of the intersectionpoint
// 	if (int_pts[0]->getPar(dir) > int_pts[1]->getPar(dir))
// 	    std::swap(int_pts[0], int_pts[1]);

	int coincide = checkIsoCurve(pdir, (pdir == dir), par, pts);
	if (coincide == 2)
	    return 2;    // Constant parameter curve OK
    }

    // For all (one or two) intersection point, check that the angle between
    // the direction of the candidate subdivision curve and the intersection
    // curve is not too acute.
    for (ki=0; ki<npts; ki++)
    {
	pval = int_pts[ki]->getPar(dir);
	if (fabs(par - pval) >= 0.5*frac*(tb - ta))
	    continue; 
	if (int_pts[ki]->getSingularityType() == HIGHER_ORDER_POINT ||
	    int_pts[ki]->getSingularityType() == ISOLATED_POINT)
	    continue;
	Point pos1 = int_pts[ki]->getPar1Point();
	Point tang1 = int_pts[ki]->getTangent();
	const double *param = (dir == pdir) ? int_pts[ki]->getPar1() 
	    : int_pts[ki]->getPar2();
	Point pos2(param[0], param[1]);
// 	func_int->point(pos2, param);
	Point der_u, der_v;
	func_int->derivs(param[0], param[1], der_u, der_v);
	double alpha = tang1.angle((pdir == 0) ?
				   Point(0.0, 1.0) : Point(1.0, 0.0)); //der_v : der_u);
// 	double alpha = tang1.angle((pdir == 0) ? der_v : der_u);
	if (pos1.dist(pos2) > eps)
	{
	    Point vec = pos1 - pos2;
	    alpha += vec.angle_smallest((pdir == 0) ?
					Point(0.0, 1.0) : Point(1.0, 0.0)); //der_v : der_u);
	}
	if (alpha < ang_tol || fabs(M_PI-alpha) < ang_tol)
	    return 1;

	if (int_pts[ki]->getSingularityType() == BRANCH_POINT)
	{
	    tang1 = int_pts[ki]->getTangent(true);
	    alpha = tang1.angle((pdir == 0) ?
				Point(0.0, 1.0) : Point(1.0, 0.0)); //der_v : der_u);
	    if (pos1.dist(pos2) > eps)
	    {
		Point vec = pos1 - pos2;
		alpha += vec.angle_smallest((pdir == 0) ?
					    Point(0.0, 1.0) : Point(1.0, 0.0)); //der_v : der_u);
	    }
	    if (alpha < ang_tol || fabs(M_PI-alpha) < ang_tol)
		return 1;
	}
    }

    return 2;

//     int pdir = (dir < 2) ? dir : dir-2;
//     int ki;

//     Param2FunctionInt* func_int = func_int_->getParam2FunctionInt();

//     if (int_pts.size() > 2)
// 	return 1;   // Not a good subdivision parameter

//     if (int_pts.size() == 2)
//     {
// 	// Check if the complete subdivision curve is in intersection.
// 	// First check if the two intersection points lie at the same boundary.
// 	// In that case, the subdivision parameter cannot be accepted
// 	if (int_results_->atSameBoundary(int_pts[0], int_pts[1]))
// 	    return 1;

// 	// Check the sequence of the intersectionpoint
// 	if (int_pts[0]->getPar(dir) > int_pts[1]->getPar(dir))
// 	    std::swap(int_pts[0], int_pts[1]);

// 	int coincide = checkIsoCurve(pdir, (pdir == dir), par, int_pts);
// 	if (coincide == 2)
// 	    return 2;    // Constant parameter curve OK
//     }

//     // For all (one or two) intersection point, check that the angle between
//     // the direction of the candidate subdivision curve and the intersection
//     // curve is not too acute.
//     for (ki=0; ki<int(int_pts.size()); ki++)
//     {
// 	Point tang1 = int_pts[ki]->getTangent();
// 	const double *param = (dir == pdir) ? int_pts[ki]->getPar1() 
// 	    : int_pts[ki]->getPar2();
// 	Point der_u, der_v;
// 	func_int->derivs(param[0], param[1], der_u, der_v);
// 	double alpha = tang1.angle((pdir == 0) ? der_v : der_u);
// 	if (alpha < epsge_->getRefAng())
// 	    return 1;
//     }

//     return 2;
}

//===========================================================================
//
// Purpose :
//
// Written by : Vibeke Skytt, 1104
//
//===========================================================================
int 
Par2FuncIntersector::checkIsoCurve(int pdir, bool first, double par,
				   vector<shared_ptr<IntersectionPoint> > int_pts)
//===========================================================================
{
    Param2FunctionInt* func_int = func_int_->getParam2FunctionInt();
    shared_ptr<ParamSurface> surf = func_int->getParamSurface();
    // Represent the candidate subdivision curve as curve-in-surface
    // curve
    Point parval[2];
    for (int ki=0; ki<2; ki++) {
       if (first) {
	   parval[ki].setValue(int_pts[ki]->getPar1()[0], 
			       int_pts[ki]->getPar1()[1]);
       } else {
	   parval[ki].setValue(int_pts[ki]->getPar2()[0], 
			       int_pts[ki]->getPar2()[1]);
       }
   }

    double ta = parval[0][1-pdir];
    double tb = parval[1][1-pdir];

    // We should also extract the corresponding space-curve from the surface.
    shared_ptr<SplineCurve> space_sub_cv;
    if (surf->instanceType() == Class_SplineSurface) {
	shared_ptr<SplineSurface> spline_sf =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
	bool pardir_is_u = (pdir == 1);
	shared_ptr<SplineCurve> space_cv(spline_sf->constParamCurve(par, pardir_is_u));
	space_sub_cv = shared_ptr<SplineCurve>(space_cv->subCurve(ta, tb));
    }
    BoundingBox box((space_sub_cv.get() != 0) ? space_sub_cv->boundingBox() : surf->boundingBox());
    double epsge = epsge_->getEpsge();
    Param0FunctionInt* C_int = C_->getParam0FunctionInt();
    ASSERT(C_int != 0);
    double C = C_int->getValue();
    bool coincide = ((fabs(box.low()[0] - C) < epsge) && (fabs(box.high()[0] - C) < epsge));

    if (coincide)
	return 2;

    return 0;
}

//===========================================================================
bool 
Par2FuncIntersector::getSubdivAtSing(int dir, double ta, double tb,
				     double& par)
//===========================================================================
{
    // Purpose : Search for a subdivision parameter in a singular
    // intersection point in a surface, the point is either already
    // found at the surface boundary or an iteration will be performed
    // to look for a singularity.

    Param2FunctionInt* func2_int = func_int_->getParam2FunctionInt();

    // Fetch singularity information from the previous intersector if
    // necessary and if such an information exist
    if (singularity_info_.get() == 0 && prev_intersector_ &&
	prev_intersector_->numParams() == 2 && 
	prev_intersector_->hasSingularityInfo()) {
	singularity_info_ = (shared_ptr<SingularityInfo>)
	    (new SingularityInfo(prev_intersector_->getSingularityInfo()));
    }
    else if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_
	    = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }

    // First fetch all singular intersection points (branch points) at
    // the surface boundaries
    bool is_critical;
    int subdiv_ok = -1;
    vector<shared_ptr<IntersectionPoint> > bd_sing;
    int_results_->getBranchPoints(bd_sing);
    for (size_t ki=0; ki<bd_sing.size(); ki++) {
	// Check the current singularity to see if it is a suitable
	// subdivision point
	par = bd_sing[ki]->getPar(dir);
	vector<shared_ptr<IntersectionPoint> > int_pts;
	is_critical = (int_results_->inInfluenceArea(dir, par, int_pts) != 0);
	if (int_pts.size() > 0) {
	    // Check if the subdivision position is OK with respect to
	    // the angle between the subdivision curve and the
	    // intersection curve
	    subdiv_ok = checkSubdivParam(dir, par, ta, tb, int_pts);
	}
	if (subdiv_ok == 0 || subdiv_ok == 2)
	    return true;
    }

    // No suitable singularity at the boundary is found. Check the
    // singularity information from the intersector.
    if (singularity_info_->hasPoint()) {
	par = singularity_info_->getParam(dir);
	double ta = func2_int->startParam(dir);
 	double tb = func2_int->endParam(dir);
	if (par > ta && par < tb) {
	    // Check the parameter value
	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = (int_results_->inInfluenceArea(dir, par, int_pts) != 0);
	    if (int_pts.size() > 0) {
		// Check if the subdivision position is OK with
		// respect to the angle between the subdivision curve
		// and the intersection curve
		subdiv_ok = checkSubdivParam(dir, par, ta, tb, int_pts);
	    }
	    if (subdiv_ok == 0 || subdiv_ok == 2)
		return true;
	}
    }

    
    // Iterate for a singularity. First check if the current surface
    // is monotone.  Make surface pointers

    Point mon_dir;
    bool monotone = func2_int->monotone(mon_dir);
    int bound_simple = 5;
    // Our second object is not singular as it is a constant,
    // equivalent to a plane for the 2d-par function setting.
    if ((!singularity_info_->iterationDone() && (monotone))
	|| singularity_info_->nmbSimple1() > bound_simple) {
	// Iterate for a singular intersection point
	
 	RectDomain domain1, domain2;
	shared_ptr<ParamSurface> surf = 
	    func2_int->getParentParamSurface(domain1);

	//double angle_tol = epsge_->getAngleTol();
	double seed[4];
	//double result_par[4];
	int constraints[2];
	double constraints_par[2];
	double limit[8];
	limit[0] = domain1.umin();
	limit[1] = domain1.umax();
	limit[2] = domain1.vmin();
	limit[3] = domain1.vmax();
	limit[4] = domain2.umin();
	limit[5] = domain2.umax();
	limit[6] = domain2.vmin();
	limit[7] = domain2.vmax();

	// Get seed. If intersection points at the boundaries exist,
	// use the middle point between the intersections
	if (int_results_->hasIntersectionPoints()) {
	    int_results_->getMidParameter(seed);
	} else {
	    // Use the midpoint
	    for (int ki=0; ki<4; ki++)
		seed[ki] = 0.5*(limit[2*ki] + limit[2*ki+1]);
	}
	constraints[0] = constraints[1] = 0;
	constraints_par[0] = constraints_par[1] = 0.0;

// 	// Iterate
// 	int isfound;
// 	isfound = extremalPtSurfSurf(surf.get(), surf2.get(), constraints, 
// 				     constraints_par, limit, seed, 
// 				     result_par, angle_tol);

// 	// Update singularity information
// 	singularity_info_->addIterationCount();

// 	// Check output
// 	if (isfound) {
// 	    // Check if the found singularity is a branch point
	    
// 	    singularity_info_->setSingularPoint(result_par, 4);

// 	    // Normal quality check for subdivision parameters
// 	    par = result_par[dir];
// 	    vector<shared_ptr<IntersectionPoint> > int_pts;
// 	    is_critical = int_results_->inInfluenceArea(dir, par, int_pts);
// 	    if (int_pts.size() > 0) {
// 		// Check if the subdivision position is OK with
// 		// respect to the angle between the subdivision curve
// 		// and the intersection curve
// 		subdiv_ok = checkSubdivParam(dir, par, int_pts);
// 	    }
// 	    if (subdiv_ok == 0 || subdiv_ok == 2)
// 		return true;
// 	}
    }
    // @@sbr Only relevant for sf-sf-intersection
//     singularity_info_->addSimpleCount(simple1, simple2);
	    
    return false;  // No singularity suitable for subdivision is found
}


//===========================================================================
void Par2FuncIntersector::
writeDebugConnect(vector<pair<shared_ptr<IntersectionPoint>,
		  IntPtClassification> >& bd_ints)
//===========================================================================
{
    // Purpose : Write debug info to a file

    static int number = 0;
    number++;
    char buffer[20];
    sprintf(buffer, "debug/connect.out");
    sprintf(buffer, "%i",number);

    Param2FunctionInt* func2_int = func_int_->getParam2FunctionInt();
    ASSERT(func2_int != 0);
    shared_ptr<ParamSurface> sf = func2_int->getParamSurface();

    // Open debug output

   std::ofstream debug(buffer);
    int ki, kj;
    for (ki=0; ki<int(bd_ints.size()); ki++)
    {
	vector<double> par = bd_ints[ki].first->getPar();
	for (kj=0; kj<int(par.size()); kj++)
	    debug << par[kj] << " ";
	debug << bd_ints[ki].second << std::endl;
    }
    sf->writeStandardHeader(debug);
    sf->write(debug);    
}


//===========================================================================


} // namespace Go
