//===========================================================================
//                                                                           
// File: ParamPointInt.C 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: ParamPointInt.C,v 1.14 2006-03-14 09:29:11 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/utils/CompositeBox.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/utils/DirectionCone.h"


using std::vector;


namespace Go {


//===========================================================================
ParamPointInt::ParamPointInt(shared_ptr<Point> point, ParamGeomInt* parent)
    : ParamGeomInt(parent), point_(point),
      coords_(point->begin(), point->end())
//===========================================================================
{
  dim_ = point_->dimension();
}


//===========================================================================
CompositeBox  ParamPointInt::compositeBox() const
//===========================================================================
{
  CompositeBox bbox(*(point_.get()), *(point_.get()));
  return bbox;
}


//===========================================================================
RotatedBox  ParamPointInt::getRotatedBox(std::vector<Point>& axis) const
//===========================================================================
{
    RotatedBox box(point_->begin(), dimension(), 1, 1, &axis[0]);
    return box;
}


//===========================================================================
DirectionCone ParamPointInt::directionCone() const
//===========================================================================
{
  Point centre(dim_);
  centre.setValue(0.0);
  DirectionCone cone(centre, 0.0);
  return cone;
}


//===========================================================================
void ParamPointInt::
getBoundaryObjects(vector<shared_ptr<BoundaryGeomInt> >& bd_objs)
//===========================================================================
{
  // A point has no boundaries
  return;
}


//===========================================================================


} // namespace Go
