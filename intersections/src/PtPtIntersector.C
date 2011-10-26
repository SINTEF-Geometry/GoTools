//===========================================================================
//                                                                           
// File: CvCvIntersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: PtPtIntersector.C,v 1.12 2006-03-09 15:16:57 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/PtPtIntersector.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"


using std::shared_ptr;


namespace Go {


//===========================================================================
PtPtIntersector::PtPtIntersector(shared_ptr<ParamGeomInt> point1, 
				 shared_ptr<ParamGeomInt> point2,
				 shared_ptr<GeoTol> epsge, 
				 Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector2Obj(point1, point2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  
}

//===========================================================================
PtPtIntersector::PtPtIntersector(shared_ptr<ParamGeomInt> point1, 
				 shared_ptr<ParamGeomInt> point2,
				 double epsge,  
				 Intersector *prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector2Obj(point1, point2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  
}

//===========================================================================
PtPtIntersector::~PtPtIntersector()
//===========================================================================
{
  // Currently empty
}

//===========================================================================
shared_ptr<Intersector> 
PtPtIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* prev,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
  // No lower order intersector exist!
  // It does not make sense to return anything.
  // We should never enter this routine
  shared_ptr<PtPtIntersector> curr_inter; 

  return curr_inter;
}

//===========================================================================
int PtPtIntersector::checkCoincidence()
//===========================================================================
{
  // Coincidence between two points

  return 0;
}

//===========================================================================
void PtPtIntersector::microCase()
//===========================================================================
{
    // At most one intersection points is expected. Check if there exist
    // one already
    int nmb_pt = int_results_->numIntersectionPoints();
    if (nmb_pt > 0)
	return;  // Nothing more to do

  // Check for an intersection between the two points. First
  // fetch the point instances
  Point pt1, pt2;
  double tpar = 0.0;
  obj_int_[0]->point(pt1, &tpar);  // The parameter value is not relevant 
                               // in this case
  obj_int_[1]->point(pt2, &tpar);
  if (pt1.dist(pt2) <= epsge_->getEpsge())
    {
      int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], getTolerance(), 0, 0);
    }
}

//===========================================================================
int PtPtIntersector::updateIntersections()
//===========================================================================
{
  // Interation does not make sense
  return 0;
}

//===========================================================================
//
// Purpose : Given points.
//
// Written by : Vibeke Skytt, 1204
//
//===========================================================================
int PtPtIntersector::linearCase()
//===========================================================================
{
  return 0;
}

//===========================================================================
int PtPtIntersector::doSubdivide()
//===========================================================================
{
  // Subdivision does not make sense in this case
  return 0;
}

} // namespace Go
