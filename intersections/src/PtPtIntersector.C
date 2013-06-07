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

#include "GoTools/intersections/PtPtIntersector.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"



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
