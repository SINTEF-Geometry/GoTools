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

#include <algorithm>
#include "GoTools/geometry/RectDomain.h"

#ifdef __BORLANDC__
#include <cmath> // For fabs. Should be required by VC++ and GCC as well...
using std::fabs;
#endif

namespace Go {


//===========================================================================
RectDomain::RectDomain(const Array<double, 2>& corner1, 
		       const Array<double, 2>& corner2)
//===========================================================================
{
    ll_[0] = std::min(corner1[0], corner2[0]);
    ur_[0] = std::max(corner1[0], corner2[0]);
    ll_[1] = std::min(corner1[1], corner2[1]);
    ur_[1] = std::max(corner1[1], corner2[1]);
}


//===========================================================================
RectDomain::~RectDomain()
//===========================================================================
{
}


//  //===========================================================================
//  DomainType RectDomain::domainType() const
//  //===========================================================================
//  {
//      return RectDomain;
//  }


//===========================================================================
bool RectDomain::isInDomain(const Array<double, 2>& point, 
			    double tolerance) const
//===========================================================================
{
    if (point[0] < ll_[0] - tolerance)
	return false;
    if (point[0] > ur_[0] + tolerance)
	return false;
    if (point[1] < ll_[1] - tolerance)
	return false;
    if (point[1] > ur_[1] + tolerance)
	return false;
    return true;
}


//===========================================================================
int RectDomain::isInDomain2(const Array<double, 2>& point, 
			    double tolerance) const
//===========================================================================
{
    if (point[0] < ll_[0] - tolerance)
	return 0;
    if (point[0] > ur_[0] + tolerance)
	return 0;
    if (point[1] < ll_[1] - tolerance)
	return 0;
    if (point[1] > ur_[1] + tolerance)
	return 0;
    if (point[0] < ll_[0] + tolerance)
	return 2;
    if (point[0] > ur_[0] - tolerance)
	return 2;
    if (point[1] < ll_[1] + tolerance)
	return 2;
    if (point[1] > ur_[1] - tolerance)
	return 2;
    return 1;
}


//===========================================================================
bool RectDomain::isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const
//===========================================================================
{
    // Being on the boundary implies being in the domain
    if (!isInDomain(point, tolerance))
	return false;

    if (point[0] < ll_[0] + tolerance)
	return true;
    if (point[0] > ur_[0] - tolerance)
	return true;
    if (point[1] < ll_[1] + tolerance)
	return true;
    if (point[1] > ur_[1] - tolerance)
	return true;
    return false;
}

//===========================================================================
bool RectDomain::isOnBoundary(const Array<double, 2>& point, 
			      double tolerance, int& bd, int& bd2) const
//===========================================================================
{
  // Being on the boundary implies being in the domain
  bd = bd2 = -1;
  if (!isInDomain(point, tolerance))
    return false;

  if (point[0] < ll_[0] + tolerance)
    {
      bd = 0;
      if (point[1] < ll_[1] + tolerance)
	bd2 = 2;
      else if (point[1] > ur_[1] - tolerance)
	bd2 = 3;
      return true;
    }
  if (point[0] > ur_[0] - tolerance)
    {
      bd = 1;
      if (point[1] < ll_[1] + tolerance)
	bd2 = 2;
      else if (point[1] > ur_[1] - tolerance)
	bd2 = 3;
      return true;
    }
  if (point[1] < ll_[1] + tolerance)
    {
      bd = 2;
      return true;
    }
  if (point[1] > ur_[1] - tolerance)
    {
      bd = 3;
      return true;
    }
  return false;
}


//===========================================================================
bool RectDomain::isOnCorner(const Array<double, 2>& point, 
			      double tolerance) const
//===========================================================================
{
    // Being on the corner implies being in the domain
    if (!isInDomain(point, tolerance))
	return false;

    if (point[0] < ll_[0] + tolerance && point[1] < ll_[1] + tolerance)
	return true;
    if (point[0] < ll_[0] + tolerance && point[1] > ur_[1] - tolerance)
	return true;
    if (point[0] > ur_[0] - tolerance && point[1] < ll_[1] + tolerance)
	return true;
    if (point[0] > ur_[0] - tolerance && point[1] > ur_[1] - tolerance)
	return true;
    return false;
}



//===========================================================================
int RectDomain::whichBoundary(const Array<double, 2>& point1, 
			      const Array<double, 2>& point2, 
			      double tolerance) const
//===========================================================================
{
    if (fabs(point1[0] - ll_[0]) < tolerance && 
	fabs(point2[0] - ll_[0]) < tolerance)
	return 0;   // umin
    if (fabs(point1[0] - ur_[0]) < tolerance && 
	fabs(point2[0] - ur_[0]) < tolerance)
	return 1;   // umax
    if (fabs(point1[1] - ll_[1]) < tolerance && 
	fabs(point2[1] - ll_[1]) < tolerance)
	return 2;   // vmin
    if (fabs(point1[1] - ur_[1]) < tolerance && 
	fabs(point2[1] - ur_[1]) < tolerance)
	return 3;   // vmax
    return -1;  // Not a common boundary
}

//===========================================================================
void RectDomain::closestInDomain(const Array<double, 2>& point,
				 Array<double, 2>& clo_pt,
				 double tolerance) const
//===========================================================================
{
    clo_pt = point;
    if (clo_pt[0] < ll_[0] - tolerance)
	clo_pt[0] = ll_[0];
    if (clo_pt[0] > ur_[0] + tolerance)
	clo_pt[0] = ur_[0];
    if (clo_pt[1] < ll_[1] - tolerance)
	clo_pt[1] = ll_[1];
    if (clo_pt[1] > ur_[1] + tolerance)
	clo_pt[1] = ur_[1];
}


//===========================================================================
void RectDomain::closestOnBoundary(const Array<double, 2>& point,
				   Array<double, 2>& clo_bd_pt,
				   double tolerance) const
//===========================================================================
{
    clo_bd_pt = point;
    double dist1 = fabs(point[0] - ll_[0]);
    double dist2 = fabs(point[0] - ur_[0]);
    double dist3 = fabs(point[1] - ll_[1]);
    double dist4 = fabs(point[1] - ur_[1]);
    if (std::min(dist1, dist2) < std::min(dist3, dist4)) {
	if (dist1 < dist2) {
	    clo_bd_pt[0] = ll_[0];
	} else {
	    clo_bd_pt[0] = ur_[0];
	}
    } else {
	if (dist3 < dist4) {
	    clo_bd_pt[1] = ll_[1];
	} else {
	    clo_bd_pt[1] = ur_[1];
	}
    }

    return;
}

//===========================================================================
void RectDomain::addUnionWith(const RectDomain& rd)
//===========================================================================
{
    if (rd.ll_[0] < ll_[0])
	ll_[0] = rd.ll_[0];
    if (rd.ll_[1] < ll_[1])
	ll_[1] = rd.ll_[1];
    if (rd.ur_[0] > ur_[0])
	ur_[0] = rd.ur_[0];
    if (rd.ur_[1] > ur_[1])
	ur_[1] = rd.ur_[1];
}

//===========================================================================
void RectDomain::intersectWith(const RectDomain& rd)
//===========================================================================
{
    ll_[0] = std::max(rd.ll_[0], ll_[0]);
    ll_[1] = std::max(rd.ll_[1], ll_[1]);
    ur_[0] = std::min(rd.ur_[0], ur_[0]);
    ur_[1] = std::min(rd.ur_[1], ur_[1]);
}

//===========================================================================
  bool RectDomain::overlap(const RectDomain& rd, double tol)
//===========================================================================
{
  if (ll_[0] > rd.ur_[0]+tol)
    return false;
  if (ur_[0] < rd.ll_[0]-tol)
    return false;
  if (ll_[1] > rd.ur_[1]+tol)
    return false;
  if (ur_[1] < rd.ll_[1]-tol)
    return false;
  return true;
  // if (ll_[0] < rd.ur_[0]+tol && ll_[0] > rd.ll_[0]-tol)
  //   return true;
  // if (ll_[1] < rd.ur_[1]+tol && ll_[1] > rd.ll_[1]-tol)
  //   return true;
  // if (rd.ll_[0] < ur_[0]+tol && rd.ll_[0] > ll_[0]-tol)
  //   return true;
  // if (rd.ll_[1] < ur_[1]+tol && rd.ll_[1] > ll_[1]-tol)
  //   return true;
  // return false;
}

//===========================================================================
bool RectDomain::overlap(const RectDomain& rd, double tol1, double tol2)
//===========================================================================
{
  if (ll_[0] > rd.ur_[0]+tol1)
    return false;
  if (ur_[0] < rd.ll_[0]-tol1)
    return false;
  if (ll_[1] > rd.ur_[1]+tol2)
    return false;
  if (ur_[1] < rd.ll_[1]-tol2)
    return false;
  return true;
}
} // namespace Go
