//===========================================================================
//                                                                           
// File: RectDomain.C                                               
//                                                                           
// Created: Fri Sep  1 15:37:23 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectDomain.C,v 1.17 2009-05-13 07:30:53 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
	return 2;   // umax
    if (fabs(point1[1] - ur_[1]) < tolerance && 
	fabs(point2[1] - ur_[1]) < tolerance)
	return 3;   // umax
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
	clo_pt[0] = ll_[0];
    if (clo_pt[1] > ur_[1] + tolerance)
	clo_pt[0] = ur_[0];
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

} // namespace Go
