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

#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionLink.h"

using std::vector;
using namespace Go;


//===========================================================================
void SfSfIntersector::handleComplexity()
//===========================================================================
{
    // Purpose : Complexity not reduced. Try to connect intersection
    // points at the boundaries.

    // @@@ TODO. Here we do implicitization as the first choice. We are already
    // in the SfSfIntersector, so virtual functions are not necessary.
    // Implicitization and further treatment is performed in a separate function
    // in this file.
    // The implicit surface is made or fetched. Use the same tests before 
    // implicitization is performed as before.
    // The other surface is put into the implicit one, and the intersection
    // in the functional intersector is performed (call to compute). The problems
    // are mainly:
    // 1. This will be a complex situation in 1D
    // 2. To connect intersection results found at the boundary and the
    //    results from 1D
    // 3. Results from 1D lacks one parameter pair. This information must be
    //    generated. Closest point?
    // It is still some open questions. Should we give information about boundary
    // intersections present at this point to the functional intersector? That
    // would simplify the merging of results later. On the other hand may the
    // intersections found at this stage be quite messy, and make the 1D intersection
    // problem more complex than necessary.
    // THIS IS NOT THE PLACE TO START. ONLY IMPLEMENTED AFTER WE HAVE SOME EXPERIENCE
    // WITH THE OTHER IMPLICITIZATION.

    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (getenv("COMPLEX_SFSF") && *(getenv("COMPLEX_SFSF"))=='1')
    {
	std::cout << "Complex situation, " << (int)int_pts.size();
	std::cout << " intersection points" << std::endl;
	for (size_t k1=0; k1<int_pts.size(); k1++)
	{
	    for (int k2=0; k2<4; k2++)
		std::cout << int_pts[k1]->getPar(k2) << "  ";
	    std::cout << int_pts[k1]->getDist() << std::endl;
	}
	writeDebugLinear(int_pts);
    }

    if (int_pts.size() == 0)
	return;    // No intersection points to connect

    // Sort intersections according to existing connections
    vector<shared_ptr<IntersectionGroup> > groups;
    makeIntersectionGroups(int_pts, groups);
    if (groups.size() == 1)
	return;   // All intersection points are somehow connected already

    classifyGroups(groups);

    // Count groups with intersection points having non-touching tangent directions 
    size_t ki;
    int not_touch = 0;
    for (ki=0; ki<groups.size(); ki++)
    {
	if (groups[ki]->main_type != DIR_TOUCH)
	    not_touch++;
    }
   
    if (not_touch <= 1)
	return;   // Nothing to connect

    if (not_touch == 2)
    {
	// Two groups of points. These should either be connected or not. Check.
	shared_ptr<IntersectionGroup> curr[2];
	int idx = 0;
	for (ki=0; ki<groups.size(); ki++)
	    if (groups[ki]->main_type != DIR_TOUCH)
		curr[idx++] = groups[ki];

	tryConnectGroups(curr[0], curr[1]);
    }
    else
    {
	// More than two groups. Try to separate the groups into sets of not
	// more than two.
    }
}

//===========================================================================
void SfSfIntersector::selfintComplex(IntersectionPoint *sing)
//===========================================================================
{

    double ptol = epsge_->getRelParRes();

    // Remove double instances of points
    int_results_->removeDoublePoints();

    // Fetch intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (getenv("COMPLEX_SFSF") && *(getenv("COMPLEX_SFSF"))=='1')
    {
	std::cout << "Complex situation, " << (int)int_pts.size();
	std::cout << " intersection points" << std::endl;
	for (size_t k1=0; k1<int_pts.size(); k1++)
	{
	    for (int k2=0; k2<4; k2++)
		std::cout << int_pts[k1]->getPar(k2) << "  ";
	    std::cout << int_pts[k1]->getDist() << std::endl;
	}
	writeDebugLinear(int_pts);
    }

    // No self intersection singularites except the specified one are expected.
    // Dismiss such points
    // First compute the maximum distance between the parameters of a point in one
    // parameter direction to set a reasonable tolerance for equality in this context
    double fac = 1.0e-4;
    double par_dist[2];
    par_dist[0] = par_dist[1] = 0.0;
    double dp1, dp2;
    size_t ki, kh;
    for (ki=0; ki<int_pts.size(); ki++)
    {
	par_dist[0] = std::max(par_dist[0], 
			       fabs(int_pts[ki]->getPar1()[0] - int_pts[ki]->getPar2()[0]));
	par_dist[1] = std::max(par_dist[1], 
			       fabs(int_pts[ki]->getPar1()[1] - int_pts[ki]->getPar2()[1]));
    }

    for (ki=0; ki<int_pts.size(); ki++)
    {
	if (int_pts[ki].get() == sing)
	    continue;

	if (int_pts[ki]->isConnectedTo(sing))
	    continue;  // Singularity link

	dp1 = fabs(int_pts[ki]->getPar1()[0] - int_pts[ki]->getPar2()[0]);
	dp2 = fabs(int_pts[ki]->getPar1()[1] - int_pts[ki]->getPar2()[1]);
	if ((dp1 < std::max(ptol, fac*par_dist[0]) &&
	     dp2 < std::max(ptol, fac*par_dist[1])) ||
	    int_pts[ki]->getSingularityType() == HIGHER_ORDER_POINT)
	{
	    int_results_->removeIntPoint(int_pts[ki]);
	    int_pts.erase(int_pts.begin()+ki);
	    ki--;
	}
    }

    if (int_pts.size() == 0)
	return;    // No intersection points to connect

    // Set parent point information (twin)
    size_t curr_size = int_pts.size();
    double param[4];
    for (ki=0; ki<curr_size; ki++)
    {
	for (kh=ki+1; kh<curr_size; kh++)
	{
// 	    if (int_pts[ki].get() == int_pts[kh].get())
// 		continue;

	    dp1 = fabs(int_pts[ki]->getPar1()[0] - int_pts[kh]->getPar2()[0]);
	    dp2 = fabs(int_pts[ki]->getPar1()[1] - int_pts[kh]->getPar2()[1]);
	    if (dp1 < std::max(ptol, fac*par_dist[0]) &&
		dp2 < std::max(ptol, fac*par_dist[1]))
	    {
		// Make sure that the doublettes have exactly the same parameter pairs
		if (int_pts[kh]->getDist() < int_pts[ki]->getDist())
		{
		    param[0] = int_pts[kh]->getPar(2);
		    param[1] = int_pts[kh]->getPar(3);
		    param[2] = int_pts[kh]->getPar(0);
		    param[3] = int_pts[kh]->getPar(1);
		    int_pts[ki]->replaceParameter(param);
		}
		else
		{
		    param[0] = int_pts[ki]->getPar(2);
		    param[1] = int_pts[ki]->getPar(3);
		    param[2] = int_pts[ki]->getPar(0);
		    param[3] = int_pts[ki]->getPar(1);
		    int_pts[kh]->replaceParameter(param);
		}
		if (int_pts[ki]->numNeighbours() >= int_pts[kh]->numNeighbours())
		    int_pts[kh]->setParentPoint(int_pts[ki]);
		else
		    int_pts[ki]->setParentPoint(int_pts[kh]);
		break;  // Doublette found
	    }
	}
    }

    // Sort intersections according to existing connections
    vector<shared_ptr<IntersectionGroup> > groups;
    makeIntersectionGroups(int_pts, groups);
    if (groups.size() == 1)
	return;   // All intersection points are somehow connected already

    classifyGroups(groups);

    // Count groups with intersection points having non-touching tangent directions 
    int not_touch = 0;
    for (ki=0; ki<groups.size(); ki++)
    {
	if (groups[ki]->main_type != DIR_TOUCH)
	    not_touch++;
    }
   
    if (not_touch <= 1)
	return;   // Nothing to connect

    // Dismiss the group containing the singularity
    int kj;
    bool found = false;
    if (sing)
    {
	for (kj=0; kj<(int)(groups.size()); kj++)
	{
	    for (ki=0; ki<groups[kj]->points.size(); ki++)
		if (groups[kj]->points[ki].get() == sing)
		{
		    found = true;
		    break;
		}

	    if (groups[kj]->main_type == DIR_TOUCH || ki<groups[kj]->points.size())
	    {
		groups.erase(groups.begin()+kj, groups.begin()+kj+1);
		kj--;
	    }
	}
    
	// Dismiss double instances of intersection groups
	mergeGroups(groups, sing);
    }

    if (found)
    {
	for (ki=0; ki<groups.size(); ki++)
	{
	    if (groups[ki]->main_type == DIR_HIGHLY_SINGULAR)
		continue;   // Probably a new instance of the already known singularity
	    connectToSing(groups[ki], sing);
	}
    }
    else if (groups.size() == 2)
    {
	// Two groups of points. These should either be connected or not. Check.
	tryConnectGroups(groups[0], groups[1]);
    }
    else
    {
	// More than two groups. Try to separate the groups into sets of not
	// more than two.
    }
}

bool compare_ints1(const shared_ptr<IntersectionPoint> l1, 
		   const shared_ptr<IntersectionPoint>  l2)
{ 
    if (l1->getPar(0) < l2->getPar(0))
	return 1;
    else if (l1->getPar(0) > l2->getPar(0))
	return 0;
    else if (l1->getPar(1) < l2->getPar(1))
	return 1;
    else //if (l1->getPar(1) > l2->getPar(1))
	return 0;
}


bool compare_ints2(const shared_ptr<IntersectionPoint> l1, 
		   const shared_ptr<IntersectionPoint>  l2)
{ 
    if (l1->getPar(1) < l2->getPar(1))
	return 1;
    else if (l1->getPar(1) > l2->getPar(1))
	return 0;
    else if (l1->getPar(0) < l2->getPar(0))
	return 1;
    else //if (l1->getPar(0) > l2->getPar(0))
	return 0;
}

//===========================================================================
void SfSfIntersector::selfintComplex2()
//===========================================================================
{
    double ptol = epsge_->getRelParRes();

    // Remove double instances of points
    int_results_->removeDoublePoints();

    // Post iterate
    /*for (int kh=0; kh<4; kh++)
      postIterate2(0, kh, 0.0, false, false);*/

    // Fetch intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (getenv("COMPLEX_SFSF") && *(getenv("COMPLEX_SFSF"))=='1')
    {
	std::cout << "Complex situation 2, " << (int)int_pts.size();
	std::cout << " intersection points" << std::endl;
	for (size_t k1=0; k1<int_pts.size(); k1++)
	{
	    for (int k2=0; k2<4; k2++)
		std::cout << int_pts[k1]->getPar(k2) << "  ";
	    std::cout << int_pts[k1]->getDist() << std::endl;
	}
	writeDebugLinear(int_pts);
    }

    // No self intersection singularites are expected. Dismiss such points
    // First compute the maximum distance between the parameters of a point in one
    // parameter direction to set a reasonable tolerance for equality in this context
    size_t ki, kj, kh;
    double fac = 1.0e-4;
    double par_dist[2];
    par_dist[0] = par_dist[1] = 0.0;
    double dp1, dp2;
    shared_ptr<IntersectionPoint> dummy;
    for (ki=0; ki<int_pts.size(); ki++)
    {
	par_dist[0] = std::max(par_dist[0], 
			       fabs(int_pts[ki]->getPar1()[0] - int_pts[ki]->getPar2()[0]));
	par_dist[1] = std::max(par_dist[1], 
			       fabs(int_pts[ki]->getPar1()[1] - int_pts[ki]->getPar2()[1]));
    }

    for (ki=0; ki<int_pts.size(); ki++)
    {
	dp1 = fabs(int_pts[ki]->getPar1()[0] - int_pts[ki]->getPar2()[0]);
	dp2 = fabs(int_pts[ki]->getPar1()[1] - int_pts[ki]->getPar2()[1]);
	if ((dp1 < std::max(ptol, fac*par_dist[0]) &&
	     dp2 < std::max(ptol, fac*par_dist[1])) ||
	    int_pts[ki]->getSingularityType() == HIGHER_ORDER_POINT)
	{
	    int_results_->removeIntPoint(int_pts[ki]);
	    int_pts.erase(int_pts.begin()+ki);
	    ki--;
	}
    }

    // Make sure that doublettes exist for each remaining point since there is an
    // unsymmetry in the boundary intersections
    size_t curr_size = int_pts.size();
    double param[4];
    size_t k1, k2;
    for (ki=0; ki<curr_size; ki++)
    {
	for (kj=0; kj<curr_size; kj++)
	{
	    if (int_pts[ki].get() == int_pts[kj].get())
		continue;

	    dp1 = fabs(int_pts[ki]->getPar1()[0] - int_pts[kj]->getPar2()[0]);
	    dp2 = fabs(int_pts[ki]->getPar1()[1] - int_pts[kj]->getPar2()[1]);
	    if (dp1 < std::max(ptol, fac*par_dist[0]) &&
		dp2 < std::max(ptol, fac*par_dist[1]))
	    {
		if (kj > ki)
		{
		    // Make sure that the doublettes have exactly the same parameter pairs
		    if (int_pts[kj]->getDist() < int_pts[ki]->getDist())
		    {
			k1 = kj;
			k2 = ki;
		    }
		    else
		    {
			k1 = ki;
			k2 = kj;
		    }
		    param[0] = int_pts[k1]->getPar(2);
		    param[1] = int_pts[k1]->getPar(3);
		    param[2] = int_pts[k1]->getPar(0);
		    param[3] = int_pts[k1]->getPar(1);
		    if (std::min(fabs(param[2]-param[0]),fabs(param[3]-param[1])) >= ptol)
		    {
			// Make sure to place points exactly at the boundaries
			if (fabs(param[2]-param[0]) < fabs(param[3]-param[1]))
			{
			    double tpar;
			    if (fabs(obj_int_[0]->startParam(0)-param[0]) <
				fabs(obj_int_[0]->endParam(0)-param[0]))
				tpar = obj_int_[0]->startParam(0);
			    else
				tpar = obj_int_[0]->endParam(0);
			    param[0] = param[2] = tpar;
			}
			else
			{
			    double tpar;
			    if (fabs(obj_int_[0]->startParam(1)-param[1]) <
				fabs(obj_int_[0]->endParam(1)-param[1]))
				tpar = obj_int_[0]->startParam(1);
			    else
				tpar = obj_int_[0]->endParam(1);
			    param[1] = param[3] = tpar;
			}
			    
		    }
		    int_pts[k2]->replaceParameter(param);
		    std::swap(param[0],param[2]);
		    std::swap(param[1],param[3]);
		    int_pts[k1]->replaceParameter(param);

		    if (int_pts[ki]->numNeighbours() >= int_pts[kj]->numNeighbours())
			int_pts[kj]->setParentPoint(int_pts[ki]);
		    else
			int_pts[ki]->setParentPoint(int_pts[kj]);
		}
		break;  // Doublette found
	    }
	}

	if (kj >= curr_size)
	{
	    // No doublette found. Create it
	    double param[4];
	    for (kh=0; kh<4; kh++)
		param[kh] = int_pts[ki]->getPar((int)kh);

	    // Check if the doublette lies within the current parameter domains
	    for (kh=0; kh<2; kh++)
		if (param[2+kh] < obj_int_[0]->startParam((int)kh)-ptol || 
		    param[2+kh] > obj_int_[0]->endParam((int)kh)+ptol)
		    break;

	    if (kh == 2)
	    {
		for (kh=0; kh<2; kh++)
		    if (param[kh] < obj_int_[1]->startParam((int)kh)-ptol || 
			param[kh] > obj_int_[1]->endParam((int)kh)+ptol)
			break;
	    }
		
	    if (kh >= 2)
	    {
		if (std::min(fabs(param[2]-param[0]),fabs(param[3]-param[1])) >= ptol)
		{
		    // Make sure to place points exactly at the boundaries
		    if (fabs(param[2]-param[0]) < fabs(param[3]-param[1]))
		    {
			double tpar;
			if (fabs(obj_int_[0]->startParam(0)-param[0]) <
			    fabs(obj_int_[0]->endParam(0)-param[0]))
			    tpar = obj_int_[0]->startParam(0);
			else
			    tpar = obj_int_[0]->endParam(0);
			param[0] = param[2] = tpar;
		    }
		    else
		    {
			double tpar;
			if (fabs(obj_int_[0]->startParam(1)-param[1]) <
			    fabs(obj_int_[0]->endParam(1)-param[1]))
			    tpar = obj_int_[0]->startParam(1);
			else
			    tpar = obj_int_[0]->endParam(1);
			param[1] = param[3] = tpar;
		    }
		}
		int_pts[ki]->replaceParameter(param);
		shared_ptr<IntersectionPoint> tmp_pt = 
		    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], epsge_,
						       param+2, param);
		tmp_pt->setParentPoint(int_pts[ki]);
		int_pts.push_back(tmp_pt);
	    }
	}
    }
	
    if (int_pts.size() <= 1)
	return;  // No connection is possible

    // Sort the points around the parameter domain and check if the point
    // IN and OUT of the domain every second time. If two consequtive points
    // point IN or OUT, remove the points with less accuracy
    Point mid(int_pts[0]->getPar(0), int_pts[0]->getPar(1));
    for (ki=1; ki<int_pts.size(); ki++)
    {
	mid[0] += int_pts[ki]->getPar(0);
	mid[1] += int_pts[ki]->getPar(1);
    }
    mid[0] /= (int)(int_pts.size());
    mid[1] /= (int)(int_pts.size());

    Point curr(int_pts[0]->getPar(0), int_pts[0]->getPar(1));
    Point vec = curr - mid;
    Point orth(vec[1], -vec[0]);
    vector<double> ang;
    ang.reserve(int_pts.size());
    ang.push_back(0.0);
    for (ki=1; ki<int_pts.size(); ki++)
    {
	// Compute angle
	curr.setValue(int_pts[ki]->getPar(0), int_pts[ki]->getPar(1));
	Point vec2 = curr - mid;
	double tmp_ang = vec.angle(vec2);
	if (vec2*orth > 0.0)
	    tmp_ang = 2.0*M_PI - tmp_ang;
	ang.push_back(tmp_ang);
    }

    for (ki=0; ki<int_pts.size(); ki++)
    {
	for (kj=ki+1; kj<int_pts.size(); kj++)
	{
	    if (ang[kj] < ang[ki])
	    {
		std::swap(int_pts[ki],int_pts[kj]);
		std::swap(ang[ki], ang[kj]);
	    }
	}
    }

    // Fetch ParamSurface information
    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    ParamSurface *srf1 = surf1->getParamSurface().get();
    ParamSurface *srf2 = surf2->getParamSurface().get();

    // Check classification
     IntPtClassification type1, type2;
//     for (ki=0; ki<int_pts.size(); ki++)
//     {
// 	kj = ki+1;
// 	if (kj >= int_pts.size())
// 	    kj = 0;

// 	type1 = int_pts[ki]->getClassification(srf1, srf2);
// 	type2 = int_pts[kj]->getClassification(srf1, srf2);

// 	if ((type1 == DIR_IN && type2 == DIR_IN) ||
// 	    (type1 == DIR_OUT && type2 == DIR_OUT))
// 	{
// 	    // Inconsistent directions
// 	    if (int_pts[ki]->getDist() > int_pts[kj]->getDist())
// 	    {
// 		int_results_->removeIntPoint(int_pts[ki]);
// 		int_pts.erase(int_pts.begin()+ki);
// 	    }
// 	    else
// 	    {
// 		int_results_->removeIntPoint(int_pts[kj]);
// 		int_pts.erase(int_pts.begin()+kj);
// 	    }

// 	    ki--;
// 	}
//     }

    // Connect
    if (int_pts.size() == 2)
    {
	int_pts[0]->connectTo(int_pts[1],COMPLEX_SFSF);
    }
    else
    {
	for (ki=0; ki<int_pts.size(); ki++)
	{
	    kj = ki+1;
	    if (kj >= int_pts.size())
		kj = 0;

	    type1 = int_pts[ki]->getClassification(srf1, srf2);
	    type2 = int_pts[kj]->getClassification(srf1, srf2);

	    /*if (type1 == DIR_IN && (type2 == DIR_OUT ||
	      type2 == DIR_PARALLEL ||
	      type2 == DIR_PERPENDICULAR ||
	      type2 == DIR_UNDEF))
	      int_pts[ki]->connectTo(int_pts[kj],COMPLEX_SFSF);
	      else if (type2 == DIR_OUT && (type1 == DIR_PARALLEL ||
	      type1 == DIR_PERPENDICULAR ||
	      type1 == DIR_UNDEF))
	      int_pts[ki]->connectTo(int_pts[kj],COMPLEX_SFSF);*/

	    if (type1 == DIR_IN || type1 == DIR_OUT ||
		type2 == DIR_IN || type2 == DIR_OUT)
	    {
		if ((type1 == DIR_IN && type2 == DIR_IN) ||
		    (type1 == DIR_OUT && type2 == DIR_OUT))
		    continue;

		// Possibility for connection. Check if the point is on
		// different edges
		// Check also if the two points are doublettes. These should not be connected
		dp1 = fabs(int_pts[ki]->getPar1()[0] - int_pts[kj]->getPar2()[0]);
		dp2 = fabs(int_pts[ki]->getPar1()[1] - int_pts[kj]->getPar2()[1]);
		bool doublette = (dp1 < std::max(ptol, fac*par_dist[0]) &&
				  dp2 < std::max(ptol, fac*par_dist[1]));
		if (!doublette &&
		    int_results_->atDifferentBoundary(int_pts[ki], int_pts[kj]) &&
		    !int_results_->atSameBoundary(int_pts[ki], int_pts[kj]))
		{
		    // Ensure that both or no points is a twin
		    if (int_pts[ki]->hasParentPoint() && !int_pts[kj]->hasParentPoint())
		    {
			// Check if int_pts[kj] has a twin
			for (kh=0; kh<int_pts.size(); kh++)
			    if (int_pts[kh]->hasParentPoint() &&
				int_pts[kh]->parentPoint().get() == int_pts[kj].get())
				break;
			if (kh < int_pts.size())
			{
			    if (int_pts[ki]->parentPoint()->numNeighbours() < 
				int_pts[kj]->numNeighbours())
			    {
				int_pts[ki]->parentPoint()->setParentPoint(int_pts[ki]);
				int_pts[ki]->setParentPoint(dummy);
			    }
			    else
			    {
				int_pts[kj]->setParentPoint(int_pts[kh]);
				int_pts[kh]->setParentPoint(dummy);
			    }
			}
			else
			{
			    int_pts[ki]->parentPoint()->setParentPoint(int_pts[ki]);
			    int_pts[ki]->setParentPoint(dummy);
			}
		    }
		    else if (int_pts[kj]->hasParentPoint() && !int_pts[ki]->hasParentPoint())
		    {
			// Check if int_pts[kj] has a twin
			for (kh=0; kh<int_pts.size(); kh++)
			    if (int_pts[kh]->hasParentPoint() &&
				int_pts[kh]->parentPoint().get() == int_pts[ki].get())
				break;
			if (kh < int_pts.size())
			{
			    if (int_pts[kj]->parentPoint()->numNeighbours() < 
				int_pts[ki]->numNeighbours())
			    {
				int_pts[kj]->parentPoint()->setParentPoint(int_pts[kj]);
				int_pts[kj]->setParentPoint(dummy);
			    }
			    else
			    {
				int_pts[ki]->setParentPoint(int_pts[kh]);
				int_pts[kh]->setParentPoint(dummy);
			    }
			}
			else
			{
			    int_pts[kj]->parentPoint()->setParentPoint(int_pts[kj]);
			    int_pts[kj]->setParentPoint(dummy);
			}
		    }

		    int_pts[ki]->connectTo(int_pts[kj],COMPLEX_SFSF);
		    ki++;
		}
	    }
	}
    }
	    

   

/*// Sort points according to parameter direction and tangent direction
    std::sort(int_pts.begin(), int_pts.end(), compare_ints1);
    for (ki=1; ki<int_pts.size(); ki+=2)
    {
	Point tan1 = int_pts[ki-1]->getTangent();
	Point tan2 = int_pts[ki]->getTangent();
	if (tan1*tan2 <= 0.0)
	    break;
    }

    if (ki<int_pts.size())
    {
	std::sort(int_pts.begin(), int_pts.end(), compare_ints2);
	for (ki=1; ki<int_pts.size(); ki+=2)
	{
	    Point tan1 = int_pts[ki-1]->getTangent();
	    Point tan2 = int_pts[ki]->getTangent();
	    if (tan1*tan2 <= 0.0)
		break;
	}
    }
    if (ki >= int_pts.size())
    {
	for (ki=1; ki<int_pts.size(); ki+=2)
	int_pts[ki-1]->connectTo(int_pts[ki],COMPLEX_SFSF);
	} */

}


// Intersection points at the boundaries are collected into groups where
// each intersection point is connected to some other point in the group
//===========================================================================
void 
SfSfIntersector::makeIntersectionGroups(vector<shared_ptr<IntersectionPoint> >& pts,
					vector<shared_ptr<IntersectionGroup> >& groups)
//===========================================================================
{
    int size = (int)(pts.size());
    vector<bool> sorted(size, false);
    int ki, kj, kh;
    for (ki=0; ki<size; ki++)
    {
	if (sorted[ki])
	    continue;  // The point is already placed in a group

	vector<shared_ptr<IntersectionPoint> > curr_pts;
	curr_pts.push_back(pts[ki]);
	sorted[ki] = true;

	for (kj=ki+1; kj<size; kj++)
	{
	    if (sorted[kj])
		continue;
	    for (kh=0; kh<(int)(curr_pts.size()); kh++)
	    {
		if (pts[kj]->isConnectedTo(curr_pts[kh]))
		    break;
	    }
	    if (kh < (int)(curr_pts.size()))
	    {
		// A new point can be added to the group
		curr_pts.push_back(pts[kj]);
		sorted[kj] = true;
		kj = ki;   // We have to check the remaining points again
	    }
	}

	// Make a new group
	shared_ptr<IntersectionGroup> curr_group(new IntersectionGroup(curr_pts));
	groups.push_back(curr_group);
    }

    // All points are sorted
}

//===========================================================================
void
SfSfIntersector::classifyGroups(vector<shared_ptr<IntersectionGroup> >& groups)
//===========================================================================
{
     // Fetch ParamSurface information
    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    ASSERT(surf1 != 0 && surf2 != 0);

    ParamSurface *srf1 = surf1->getParamSurface().get();
    ParamSurface *srf2 = surf2->getParamSurface().get();

   size_t size = groups.size();
   size_t ki;
    for (ki=0; ki<size; ki++)
    {
	groups[ki]->classifyPoints(srf1, srf2);
	groups[ki]->classifyCurve();
    }
}

//===========================================================================
void
SfSfIntersector::tryConnectGroups(shared_ptr<IntersectionGroup> group1,
				  shared_ptr<IntersectionGroup> group2)
//===========================================================================
{
    // Given two groups of intersection points that might be interesting
    // for connection. Check if there exist an intersection internal to 
    // the surface between these groups of intersection points.
    // First fetch the closest points (parameter domain) in the two groups.
    size_t ki;
    double dist2, min_dist2 = 1.0e6;  // A large number
    double dp;
    shared_ptr<IntersectionPoint> pt1, pt2, ptmin1, ptmin2;
    for (ki=0; ki<group1->points.size(); ki++)
    {
	pt1 = group1->points[ki];
	pt2 = group2->closestToPoint(pt1.get());
	int ind;
	for (ind=0, dist2=0.0; ind<4; ind++)
	{
	    dp = pt1->getPar(ind) - pt2->getPar(ind);
	    dist2 += dp*dp;
	}
	if (dist2 < min_dist2)
	{
	    min_dist2 = dist2;
	    ptmin1 = pt1;
	    ptmin2 = pt2;
	}
    }
	
    if (ptmin1->getPoint().dist(ptmin2->getPoint()) < epsge_->getEpsge())
    {
	// The distance between the groups is less than the tolerance. Connect.
	ptmin1->connectTo(ptmin2, COMPLEX_SFSF);
	return;
    }

    // Iterate for an intersection point between the two groups
    // It could that we should test closer to the other group if one of them are
    // tangential for stability reasons
    bool found = canConnect(ptmin1, ptmin2);
    if (found)
    {
	// Find appropriate points for connection
	shared_ptr<IntersectionPoint> pnt1, pnt2;
	IntPtClassification type1, type2;
	type1 = group1->main_type;
	type2 = group2->main_type;
	if (group1->tangential)
	    pnt1 = group1->getBranchPoint();
	if (group2->tangential)
	    pnt2 = group2->getBranchPoint();
	if (pnt1.get()==0 && group1->isIso())
	    pnt1 = ptmin1;
	if (pnt2.get()==0 && group2->isIso())
	    pnt2 = ptmin2;
	if (pnt1.get()==0)
	    pnt1 = group1->getMainPoint(type1, type2);
	if (pnt2.get()==0)
	    pnt2 = group2->getMainPoint(type2, type1);
	if (pnt1.get() == 0)
	    pnt1 = ptmin1;
	if (pnt2.get() == 0)
	    pnt2 = ptmin2;

	pnt1->connectTo(pnt2, COMPLEX_SFSF);
    }
}


//===========================================================================
void
SfSfIntersector::mergeGroups(vector<shared_ptr<IntersectionGroup> >& group,
			     IntersectionPoint *sing)
//===========================================================================
{
    if (group.size() == 1)
	return;    // No possibility for double representation

    // Compute medium distance to singularity
    shared_ptr<IntersectionPoint> ptmin, ptmin2;
    double dist, med_dist = 0.0;
    size_t ki, kj;
    int kr;
    for (ki=0; ki<group.size(); ki++)
    {
	ptmin = group[ki]->closestToPoint(sing);
	for (kr=0, dist=0.0; kr<4; kr++)
	    dist += (ptmin->getPar(kr)-sing->getPar(kr))*(ptmin->getPar(kr)-sing->getPar(kr));
	med_dist += dist;
    }
    med_dist /= (double)(group.size());
	
    for (ki=0; ki<group.size(); ki++)
    {
	ptmin = group[ki]->closestToPoint(sing);
	for (kj=ki+1; kj<group.size(); kj++)
	{
	    if (group[ki]->main_type != group[kj]->main_type)
		continue;   // Not a double representation

	    ptmin2 = group[kj]->closestToPoint(sing);
	    for (kr=0, dist=0.0; kr<4; kr++)
		dist += (ptmin->getPar(kr)-ptmin2->getPar(kr))*
		    (ptmin->getPar(kr)-ptmin2->getPar(kr));
	    if (dist > 0.5*med_dist)
		continue;   // Not very close groups

	    // Remove the current group
	    group.erase(group.begin()+ki, group.begin()+ki+1);
	    ki--;

	    break;
	}
    }
}


//===========================================================================
void
SfSfIntersector::connectToSing(shared_ptr<IntersectionGroup> group,
			       IntersectionPoint *sing)
//===========================================================================
{
    shared_ptr<IntersectionPoint> ptmin;
    ptmin = group->closestToPoint(sing);

    // Find appropriate point for connection
    shared_ptr<IntersectionPoint> pnt;
    IntPtClassification type = group->main_type;
    if (group->tangential)
	pnt = group->getBranchPoint();
    if (pnt.get()==0 && group->isIso())
	    pnt = ptmin;
    if (pnt.get()==0)
	pnt = group->getMainPoint(type, DIR_HIGHLY_SINGULAR);
    if (pnt.get() == 0)
	pnt = ptmin;

    sing->connectTo(pnt, COMPLEX_SFSF);
}


//===========================================================================
void SfSfIntersector::IntersectionGroup::classifyPoints(ParamSurface *srf1,
							ParamSurface *srf2)
//===========================================================================
{
    size_t ki;
    for (ki=0; ki<points.size(); ki++)
    {
	IntPtClassification curr_type = points[ki]->getClassification(srf1, srf2);
	type[ki] = curr_type;
	if (points[ki]->numBranches() == 2)
	{
	    curr_type = points[ki]->getClassification(srf1, srf2, 1);
	    points.insert(points.begin()+ki, points[ki]);
	    type.insert(type.begin()+ki, curr_type);
	    ki++;
	}
    }
    for (ki=0; ki<type.size(); ki++)
    {
	if (main_type == DIR_UNDEF)
	    main_type = type[ki];
	else if (main_type == type[ki]);
	else if (main_type == DIR_TOUCH && type[ki] != DIR_UNDEF)
	    main_type = type[ki];
	else if (type[ki] == DIR_IN && main_type != DIR_OUT)
	    main_type = DIR_IN;
	else if (type[ki] == DIR_OUT && main_type != DIR_IN)
	    main_type = DIR_OUT;
	else if (type[ki] == DIR_PERPENDICULAR && 
		 main_type != DIR_IN && main_type != DIR_OUT)
	    main_type = DIR_PERPENDICULAR;
	else
	{
	    main_type = DIR_UNDEF;
	    break;
	}
    }
}

//===========================================================================
void SfSfIntersector::IntersectionGroup::classifyCurve()
//===========================================================================
{
    if (points.size() == 1)
	return;  // No curve to classify

    // If more than one point belong to the group, check if the
    // curve is iso parametric and/or tangential
    tangential = true;
    int dir;
    for (dir=0; dir<4; dir++)
	iso[dir] = true; 

    for (size_t ki=0; ki<points.size(); ki++)
    {
	if (!points[ki]->pointIsSingular())
	    tangential = false;  // Not all points are singular

	vector<IntersectionPoint*> neighbours;
	points[ki]->getNeighbours(neighbours);
	size_t kj, kh;
	for (kj=0; kj<neighbours.size(); kj++)
	{
	    for (kh=0; kh<points.size(); kh++)
		if (points[kh].get() == neighbours[kj])
		    break;
	    if (kh == points.size())
		continue;  // Neighbouring point not in current group
	    shared_ptr<IntersectionLink> lnk
		= points[ki]->getIntersectionLink(neighbours[kj]);
	    for (dir=0; dir<4; dir++)
		if (!(lnk->isIsoparametricIn(dir)))
		    iso[dir] = false;
	}
    }
}

//===========================================================================
shared_ptr<IntersectionPoint>
SfSfIntersector::IntersectionGroup::closestToPoint(IntersectionPoint *pt1)
//===========================================================================
{
    double dist2, min_dist2 = 1.0e6;  // A large number
    double dp;
    shared_ptr<IntersectionPoint> pt2;
    for (size_t ki=0; ki<points.size(); ki++)
    {
	int ind;
	for (ind=0, dist2=0.0; ind<4; ind++)
	{
	    dp = points[ki]->getPar(ind) - pt1->getPar(ind);
	    dist2 += dp*dp;
	}
	if (dist2 < min_dist2)
	{
	    min_dist2 = dist2;
	    pt2 = points[ki];
	}
    }
    return pt2;
}
   
//===========================================================================
shared_ptr<IntersectionPoint>
SfSfIntersector::IntersectionGroup::getMainPoint(IntPtClassification &type_curr, 
						 IntPtClassification type_other)
//===========================================================================
{
    // Try to find an intersection point of good quality which is 
    // compatible with the other type. 
    // Sort points according to quality
    size_t ki, kj, size = points.size();
    vector<double> edist(size);
    vector<size_t> perm(size);
    shared_ptr<IntersectionPoint> found_pt;

    for (ki=0; ki<size; ki++)
    {
	edist[ki] = points[ki]->getDist();
	perm[ki] = ki;
    }
    
    for (ki=0; ki<size; ki++)
	for (kj=ki+1; kj<size; kj++)
	    if (edist[perm[ki]] > edist[perm[kj]])
		std::swap(perm[ki], perm[kj]);

    // Fetch the first point of appropriate type
    switch (type_other)
    {
	case DIR_IN:
	{
	    for (ki=0; ki<size; ki++)
		if (type[perm[ki]] == DIR_OUT || type[perm[ki]] == DIR_PERPENDICULAR)
		    break;

	    if (ki == size)
	    {
		for (ki=0; ki<size; ki++)
		    if (type[perm[ki]] == DIR_PARALLEL || 
			type[perm[ki]] == DIR_HIGHLY_SINGULAR)
			break;
	    }

	    if (ki < size)
	    {
		found_pt = points[perm[ki]];
		type_curr = type[perm[ki]];
	    }
	}
	break;

	case DIR_OUT:
	{
	    for (ki=0; ki<size; ki++)
		if (type[perm[ki]] == DIR_IN || type[perm[ki]] == DIR_PERPENDICULAR)
		    break;

	    if (ki == size)
	    {
		for (ki=0; ki<size; ki++)
		    if (type[perm[ki]] == DIR_PARALLEL || 
			type[perm[ki]] == DIR_HIGHLY_SINGULAR)
			break;
	    }

	    if (ki < size)
	    {
		found_pt = points[perm[ki]];
		type_curr = type[perm[ki]];
	    }
	}
	break;

	case DIR_PERPENDICULAR:
	case DIR_PARALLEL: 
	case DIR_HIGHLY_SINGULAR:
	case DIR_UNDEF:
	case DIR_TOUCH:
	{
	    for (ki=0; ki<size; ki++)
		if (type[perm[ki]] == DIR_IN || type[perm[ki]] == DIR_OUT ||
		    type[perm[ki]] == DIR_PERPENDICULAR)
		    break;

	    if (ki == size)
	    {
		for (ki=0; ki<size; ki++)
		    if (type[perm[ki]] == DIR_PARALLEL || 
			type[perm[ki]] == DIR_HIGHLY_SINGULAR)
			break;
	    }

	    if (ki < size)
	    {
		found_pt = points[perm[ki]];
		type_curr = type[perm[ki]];
	    }
	}
	break;
	
    }
    return found_pt;
}


//===========================================================================
shared_ptr<IntersectionPoint>
SfSfIntersector::IntersectionGroup::getBranchPoint()
//===========================================================================
{
    shared_ptr<IntersectionPoint> found_pt;

    // For the time being
    return found_pt;
}


//===========================================================================
