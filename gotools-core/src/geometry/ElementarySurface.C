//==========================================================================
//                                                                          
// File: ElementarySurface.C                                                 
//                                                                          
// Created: Thu Sep  4 16:49:54 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ElementarySurface.C,v 1.3 2008-12-10 14:52:10 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/SplineSurface.h"


using std::swap;


namespace Go
{


// Constructor
//===========================================================================
ElementarySurface::ElementarySurface()
    : isSwapped_(false) //, isReversedU_(false), isReversedV_(false)
//===========================================================================
{
}


// Destructor
//===========================================================================
ElementarySurface::~ElementarySurface()
//===========================================================================
{
}


//===========================================================================
CurveLoop ElementarySurface::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
  if (isBounded())
    {
      // If surface is bounded, use outerBoundaryLoop on geometrySurface as default
      // Might be overridden in subclasses
      shared_ptr<SplineSurface> asSplineSurface = shared_ptr<SplineSurface>(geometrySurface());
      return asSplineSurface->outerBoundaryLoop(degenerate_epsilon);
    }
  else
    {
      // If surface is unbounded, give an empty loop
      CurveLoop loop;
      return loop;
    }
}


//===========================================================================
RectDomain ElementarySurface::containingDomain() const
//===========================================================================
{
    RectDomain dom = dynamic_cast<const RectDomain&>(parameterDomain());
    return dom;
}


//===========================================================================
bool ElementarySurface::inDomain(double u, double v) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    double eps = 1.0e-8;  // A small number
    return parameterDomain().isInDomain(pnt, eps);
}

//===========================================================================
Point ElementarySurface::closestInDomain(double u, double v) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    Array<double, 2> close(0.0, 0.0);
    double eps = 1.0e-8;  // A small number
    parameterDomain().closestInDomain(pnt, close, eps);
    return Point(close.x(), close.y());
}

//===========================================================================
double ElementarySurface::area(double tol) const
//===========================================================================
{
    return geometrySurface()->area(tol);
}

//===========================================================================
void ElementarySurface::getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const
//===========================================================================
{
  // VSK 0910. Maybe not the best solution in the long run
  geometrySurface()->getCornerPoints(corners);
}

//===========================================================================
bool ElementarySurface::isBounded() const
//===========================================================================
{
    // Assume unbounded by default
    return false;
}

//===========================================================================
shared_ptr<ElementaryCurve> 
ElementarySurface::getElementaryParamCurve(ElementaryCurve* space_crv, double tol) const 
//===========================================================================
{
  // Default is not simple elementary parameter curve exists
  shared_ptr<ElementaryCurve> dummy;
  return dummy;
}


////===========================================================================
//void ElementarySurface::reverseParameterDirection(bool direction_is_u)
////===========================================================================
//{
//    // If direction_is_u==true, this function will reverse the direction
//    // that is *currently* u. In other words, it depends on isSwapped_.
//
//    // NOTE: This function probably doesn't make sense for unbounded
//    // surfaces...
//
//    if (isSwapped())
//        direction_is_u = !direction_is_u;
//
//    if (direction_is_u)
//        isReversedU_ = !isReversedU_;
//    else
//        isReversedV_ = !isReversedV_;
//}
    
    
//===========================================================================
void ElementarySurface::swapParameterDirection()
//===========================================================================
{
    isSwapped_ = !isSwapped_;
}    
    
    
////===========================================================================
//bool ElementarySurface::isReversedU() const
////===========================================================================
//{
//    return isReversedU_;
//}    
    
    
////===========================================================================
//bool ElementarySurface::isReversedV() const
////===========================================================================
//{
//    return isReversedV_;
//}    
    
    
//===========================================================================
bool ElementarySurface::isSwapped() const
//===========================================================================
{
    return isSwapped_;
}    
    
    
//===========================================================================
void ElementarySurface::getOrientedParameters(double& u, double& v) const 
//===========================================================================
{
    if (isSwapped_) {
        swap(u, v);
    }
    //if (isReversedU_) {
    //    // Do something
    //}
    //if (isReversedV_) {
    //    // Do someting
    //}
}



} // namespace Go
