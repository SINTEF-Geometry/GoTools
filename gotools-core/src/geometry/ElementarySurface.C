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

#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/RectDomain.h"


using std::swap;
using std::vector;

namespace Go
{


// Constructor
//===========================================================================
ElementarySurface::ElementarySurface()
  : ParamSurface(), ptol_(1.0e-12), isSwapped_(false)
//===========================================================================
{
}


// Destructor
//===========================================================================
ElementarySurface::~ElementarySurface()
//===========================================================================
{
}


// //===========================================================================
// CurveLoop ElementarySurface::outerBoundaryLoop(double degenerate_epsilon) const
// //===========================================================================
// {
//   if (isBounded())
//     {
//       // If surface is bounded, use outerBoundaryLoop on geometrySurface as default
//       // Might be overridden in subclasses
//       shared_ptr<SplineSurface> asSplineSurface = shared_ptr<SplineSurface>(geometrySurface());
//       return asSplineSurface->outerBoundaryLoop(degenerate_epsilon);
//     }
//   else
//     {
//       // If surface is unbounded, give an empty loop
//       CurveLoop loop;
//       return loop;
//     }
// }


//===========================================================================
RectDomain ElementarySurface::containingDomain() const
//===========================================================================
{
    RectDomain dom = dynamic_cast<const RectDomain&>(parameterDomain());
    return dom;
}


//===========================================================================
bool ElementarySurface::inDomain(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isInDomain(pnt, eps);
}

//===========================================================================
int ElementarySurface::inDomain2(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isInDomain2(pnt, eps);
}

//===========================================================================
bool ElementarySurface::onBoundary(double u, double v, double eps) const 
//===========================================================================
{
    Array<double, 2> pnt(u, v);
    return parameterDomain().isOnBoundary(pnt, eps);
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
SplineSurface* ElementarySurface::asSplineSurface()
//===========================================================================
{
    return createSplineSurface();
}

//===========================================================================
bool ElementarySurface::isClosed(bool& closed_dir_u, bool& closed_dir_v) const
//===========================================================================
{
    // Assume not closed by default
    closed_dir_u = false;
    closed_dir_v = false;
    return false;
}

//===========================================================================
shared_ptr<ElementaryCurve> 
ElementarySurface::getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
					   const Point* start_par_pt, const Point* end_par_pt) const 
//===========================================================================
{
  // Default is not simple elementary parameter curve exists
  shared_ptr<ElementaryCurve> dummy;
  return dummy;
}


//===========================================================================
void ElementarySurface::turnOrientation()
//===========================================================================
{
    // Default behaviour is to swap parameter directions.
    swapParameterDirection();
}


//===========================================================================
void ElementarySurface::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    MESSAGE("reverseParameterDirection() not implemented.");
}


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

//===========================================================================
double 
ElementarySurface::nextSegmentVal(int dir, double par, bool forward, 
				  double tol) const
//===========================================================================
{
  RectDomain dom = containingDomain();
  if (dir == 0)
    return (forward) ? dom.umax() : dom.umin();
  else
    return (forward) ? dom.vmax() : dom.vmin();
}

//===========================================================================
CurveLoop ElementarySurface::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
  if (isBounded())
    {
      // Test for degeneracy.
      bool deg[4];
      if (degenerate_epsilon < 0.0)
	deg[0] = deg[1] = deg[2] = deg[3] = false; // All curves are wanted
      else
	isDegenerate(deg[0], deg[1], deg[2], deg[3], degenerate_epsilon);
      RectDomain dom = containingDomain();
      vector<double> par(4);
      par[0] = dom.vmin();
      par[1] = dom.umax();
      par[2] = dom.vmax();
      par[3] = dom.umin();
      bool u_dir = true;
      std::vector< shared_ptr< ParamCurve > >  vec;
      for (int edgenum = 0; edgenum < 4; ++edgenum) {
	if (!deg[edgenum]) {
	  vector<shared_ptr<ParamCurve> > edgecurve =
	    constParamCurves(par[edgenum], u_dir);
	  if (edgecurve.size() > 0)
	    {
	      if (edgenum == 2 || edgenum == 3)
		edgecurve[0]->reverseParameterDirection();
	      vec.push_back(edgecurve[0]);
	    }
	}
	u_dir = !u_dir;
      }

      if (vec.size() > 0)
	return CurveLoop(vec, (degenerate_epsilon < 0.0) ? DEFAULT_SPACE_EPSILON :
			 degenerate_epsilon);
      else
	{
	  CurveLoop loop;
	  return loop;
	}
    }
  else
    {
      // If surface is unbounded, give an empty loop
      CurveLoop loop;
      return loop;
    }

}


//===========================================================================
std::vector<CurveLoop> ElementarySurface::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
  // There is only one boundary loop...
  std::vector<CurveLoop> cvloopvec;
  cvloopvec.push_back(outerBoundaryLoop(degenerate_epsilon));
  return cvloopvec;
}


} // namespace Go
