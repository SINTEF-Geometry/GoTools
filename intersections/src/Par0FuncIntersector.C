//===========================================================================
//                                                                           
// File: Par0FuncIntersector.C                                               
//                                                                           
// Created: Fri Oct  1 10:36:02 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par0FuncIntersector.C,v 1.5 2006-03-09 15:16:56 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/Par0FuncIntersector.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"


using std::vector;


namespace Go {


//===========================================================================
Par0FuncIntersector::Par0FuncIntersector(shared_ptr<ParamFunctionInt> func,
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
Par0FuncIntersector::~Par0FuncIntersector()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
Par0FuncIntersector::lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
					   shared_ptr<ParamFunctionInt> obj2,
					   Intersector* prev,
					   int eliminated_parameter,
					   double eliminated_value)
//===========================================================================
{
    // Both objects should be of the type ParamPointInt.
    shared_ptr<Par0FuncIntersector> curr_inter; 
    return curr_inter;
}


//===========================================================================
int Par0FuncIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between a point and a curve can occur only if the
  // curve is degenerate. So ...
  return 0;
}


//===========================================================================
//
// Purpose : Two curves lie both within the same epsion ball. It is not a
//           simple case, and there is intersection.
//           This is an uwanted situation. We have to make a result that
//           is as consistent as possible since we cannot subdivide anymore.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004.
//===========================================================================
void Par0FuncIntersector::microCase()
//===========================================================================
{
    int nmb_pt = int_results_->numIntersectionPoints();
    if (nmb_pt > 0)
	return;  // Nothing more to do

    // Check for an intersection between the two points. First
    // fetch the point instances
    Point pt1, pt2;
    double tpar = 0.0;
    func_int_->point(pt1, &tpar);  // The parameter value is not relevant 
    // in this case
    C_->point(pt2, &tpar);
    if (pt1.dist(pt2) <= epsge_->getEpsge()) {
	int_results_->addIntersectionPoint(func_int_, C_, getTolerance(), 0, 0);
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
int Par0FuncIntersector::updateIntersections()
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

  Point func_pt, c_pt;
  func_int_->point(func_pt, NULL);
  C_->point(c_pt, NULL);

  double dist = c_pt.dist(func_pt);
  if (dist <= epsge_->getEpsge())
    {
      // An intersection point is found. Represent it in the data
      // structure
      // @@@ vsk, In sisl there is a check if the found intersection
      // point is very close to the endpoints of the curve. In that case
      // it is dismissed since the intersections at the endpoints are found
      // already. Here we have already checked if there exist any
      // intersection points. Thus, we know that there does not exist any
      // intersection point at the endpoint.
      int_results_->addIntersectionPoint(func_int_, 
					 C_,
					 getTolerance(),
					 NULL, // @@sbr How is this handled from above ...
					 NULL);

      return 1; // A new point is found
    }

  return 0;
}


//===========================================================================
//
// Purpose : Given two parametric curve and no simple case situation, 
//           subdivide the curves to produce more simple sub problems.
//           Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Sverre Briseid, SINTEF, Sep 2004
//
//===========================================================================
int Par0FuncIntersector::doSubdivide()
//===========================================================================
{
    return 0;
}


} // namespace Go
