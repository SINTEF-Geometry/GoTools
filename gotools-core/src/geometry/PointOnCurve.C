//===========================================================================
//                                                                           
// File: PointOnCurve.C
//                                                                           
// Created: Dec 2008
//                                                                           
// Author: Vibeke Skytt
// 
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/PointOnCurve.h"

using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::vector;

//===========================================================================
PointOnCurve::PointOnCurve()
//===========================================================================
{
  t1_ = t2_ = par_ = 0.0;
}

//===========================================================================
PointOnCurve::PointOnCurve(shared_ptr<ParamCurve> curve, double par)
    : par_(par), crv_(curve)
//===========================================================================
{
    t1_ = crv_->startparam();
    t2_ = crv_->endparam();

    point_ = crv_->point(par_);
}

//===========================================================================
PointOnCurve::PointOnCurve(shared_ptr<ParamCurve> curve, Point pnt)
    : point_(pnt), crv_(curve)
//===========================================================================
{

    t1_ = crv_->startparam();
    t2_ = crv_->endparam();

    double clo_dist;
    Point clo_pt;
    crv_->closestPoint(point_, t1_, t2_, par_, clo_pt, clo_dist);
}

//===========================================================================
PointOnCurve::~PointOnCurve()
//===========================================================================
{
}

//===========================================================================
Point PointOnCurve::getPos() const 
//===========================================================================
{
    return point_;
}

//===========================================================================
void PointOnCurve::evaluate(int der, std::vector<Point>& deriv) const
//===========================================================================
{
  if (crv_.get())
    crv_->point(deriv, par_, der);
}

//===========================================================================
void PointOnCurve::setParInterval(double start, double end)
//===========================================================================
{
    ASSERT(start <= par_ && par_ <= end);
    t1_ = start;
    t2_ = end;
}


