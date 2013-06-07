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

#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/SurfaceAssembly.h"

#include <vector>

using namespace Go;
using std::vector;
using std::pair;


//===========================================================================
SurfaceAssembly::SurfaceAssembly(shared_ptr<ParamSurfaceInt> surf,
				 vector<pair<double,int> > u_div,
				 vector<pair<double,int> > v_div,
				 vector<RectDomain> sing_domain,
				 double rel_par_res)
//===========================================================================
{
    // Make a copy to prepar for knot insertion (if spline)
    shared_ptr<ParamSurface> srf(surf->getParamSurface()->clone());
    if (surf->isSpline())
	surf_ = shared_ptr<ParamSurfaceInt>(new SplineSurfaceInt(srf,surf.get())); 
    else
	surf_ = shared_ptr<ParamSurfaceInt>(new ParamSurfaceInt(srf,surf.get())); 
    u_div_ = u_div;
    v_div_ = v_div;
    sing_domain_ = sing_domain;
    ptol_ = rel_par_res;
    closed_in_u_ = false;  // @@@ Only for the time being
    closed_in_v_ = false;  // @@@ Only for the time being
    idx_sub_ = 0;
    idx_assembly_ = 0;

    refineSurf();
}


//===========================================================================
SurfaceAssembly::~SurfaceAssembly()
//===========================================================================
{
    // No explicit action
}


//===========================================================================
int SurfaceAssembly::getNmbSubSurface()
//===========================================================================
{
    int kdiv1 = (int)u_div_.size();
    int kdiv2 = (int)v_div_.size();
    return (kdiv1+1)*(kdiv2+1);
}


//===========================================================================
bool SurfaceAssembly::getNextSubSurface(shared_ptr<ParamSurfaceInt>& sub_sf,
					int& idx, int& sing_idx)
//===========================================================================
{

    // Find start and end parameter values
    idx = idx_sub_;
    idx_sub_++;
    int ksize1 = (int)(u_div_.size()) + 1;
    int ksize2 = (int)(v_div_.size()) + 1;
    if (idx >= ksize1*ksize2)
	return false;

    int ind1 = (idx % ksize1);
    int ind2 = idx/ksize1;

    double ta1 = (ind1 == 0) ? surf_->startParam(0) : u_div_[ind1-1].first;
    double tb1 = (ind1 == ksize1-1) ? surf_->endParam(0) : u_div_[ind1].first;
    double ta2 = (ind2 == 0) ? surf_->startParam(1) : v_div_[ind2-1].first;
    double tb2 = (ind2 == ksize2-1) ? surf_->endParam(1) : v_div_[ind2].first;

    vector<shared_ptr<ParamSurface> > subsfs = 
	surf_->getParamSurface()->subSurfaces(ta1, ta2, tb1, tb2, ptol_);
    
    // @@@ Normally only one sub surface will be made. The exception is trimmed
    // surfaces with complex trimming. Since, at the moment, only spline surfaces
    // are handled, we take a simple solution
    if (surf_->isSpline())
	sub_sf = (shared_ptr<ParamSurfaceInt>)(new SplineSurfaceInt(subsfs[0],
								    surf_.get()));
    else
	sub_sf = (shared_ptr<ParamSurfaceInt>)(new ParamSurfaceInt(subsfs[0],
								   surf_.get()));

    sing_idx = -1;
    for (size_t ki=0; ki<sing_domain_.size(); ki++)
    {
	if (ta1 > sing_domain_[ki].umin()-ptol_ && tb1 < sing_domain_[ki].umax()+ptol_ &&
	    ta2 > sing_domain_[ki].vmin()-ptol_ && tb2 < sing_domain_[ki].vmax()+ptol_)
	{
	    sing_idx = (int)ki;
	    break;  // Sub surface contained in singular domain
	}
    }

    return true;
}


//===========================================================================
bool SurfaceAssembly::getNextAssembly(shared_ptr<ParamSurfaceInt>& assembly,
				      int& idx, bool& sing_sub)
//===========================================================================
{

    // Find start  parameter values
    idx = idx_assembly_;
    idx_assembly_++;
    int ksize1 = (int)(u_div_.size()) + 1;
    int ksize2 = (int)(v_div_.size()) + 1;
    if (idx >= ksize1*ksize2)
	return false;

    int ind1 = (idx % ksize1);
    int ind2 = idx/ksize1;

    // Check if the current assembly is in the end of the surface and
    // treated already
    if (ind1 > 0 && ind1 == ksize1-1 && (u_div_[ind1-1].second % 10) == 0)
	return false;
    if (ind2 > 0 && ind2 == ksize2-1 && (v_div_[ind2-1].second % 10) == 0)
	return false;

    double ta1 = (ind1 == 0) ? surf_->startParam(0) : u_div_[ind1-1].first;
    double ta2 = (ind2 == 0) ? surf_->startParam(1) : v_div_[ind2-1].first;

    // Find end parameter values. The size of an assembly depends on the
    // existence of singularities in the assembly boundaries.

    // Initially, the size of an assembly is set as one sub surface.
    // Check if it can be increased
    int indb1 = ind1 + 1;
    int indb2 = ind2 + 1;

    if (indb1 < ksize1 && indb2 < ksize2 && 
	(u_div_[indb1-1].second % 10) == 0 && (v_div_[indb2-1].second % 10) == 0)
    {
	indb1++;
	indb2++;
    }
    else 
    {
	if (indb1 < ksize1 && (u_div_[indb1-1].second % 10) == 0)
	    indb1++;
	if (indb2 < ksize2 && (v_div_[indb2-1].second % 10) == 0)
	    indb2++;
	/*if (indb1 == ind1+1 && indb1 < ksize1 && 
	    (indb2 == ksize2 || (v_div_[indb2-1].second % 10) == 0) &&
	    (ind2 == 0 || (v_div_[ind2-1].second % 10) == 0))
	    indb1++;
	if (indb2 == ind2+1 && indb2 < ksize2 && 
	    (indb1 == ksize1 || (u_div_[indb1-1].second % 10) == 0) &&
	    (ind1 == 0 || (u_div_[ind1-1].second % 10) == 0))
	    indb2++;*/
    }

    // Get end parameters
    double tb1 = (indb1 >= ksize1) ? surf_->endParam(0) : 
	u_div_[indb1-1].first;
    double tb2 = (indb2 >= ksize2) ? surf_->endParam(1) : 
	v_div_[indb2-1].first;

    vector<shared_ptr<ParamSurface> > subsfs = 
	surf_->getParamSurface()->subSurfaces(ta1, ta2, tb1, tb2, ptol_);
    
    // @@@ Normally only one sub surface will be made. The exception is trimmed
    // surfaces with complex trimming. Since, at the moment, only spline surfaces
    // are handled, we take a simple solution
    if (surf_->isSpline())
	assembly = (shared_ptr<ParamSurfaceInt>)(new SplineSurfaceInt(subsfs[0],
								      surf_.get()));
    else
	assembly = (shared_ptr<ParamSurfaceInt>)(new ParamSurfaceInt(subsfs[0],
								     surf_.get()));

    // Set singularity information
    sing_sub = false;
    bool potential_sing = false;
    if (ind1 > 0 && (u_div_[ind1-1].second % 10) == 2)
	potential_sing = true;
    if (indb1 < ksize1 && (u_div_[indb1-1].second % 10) == 2)
	potential_sing = true;
    if (ind2 > 0 && (v_div_[ind2-1].second % 10) == 2)
	potential_sing = true;
    if (indb2 < ksize2 && (v_div_[indb2-1].second % 10) == 2)
	potential_sing = true;
    if (potential_sing)
    {
	// Check if the current sub surface is included in a singular domain
	size_t ki;
	for (ki=0; ki<sing_domain_.size(); ki++)
	{
	    if (ta1 > sing_domain_[ki].umin()-ptol_ && tb1 < sing_domain_[ki].umax()+ptol_ &&
		ta2 > sing_domain_[ki].vmin()-ptol_ && tb2 < sing_domain_[ki].vmax()+ptol_)
		break;  // Sub surface contained in singular domain
	}
	if (ki < sing_domain_.size())
	    sing_sub = true;
    }

    // Check if the index to the next assembly should be increased
    int ind3 = (idx_assembly_ % ksize1);
//     int ind4 = idx_assembly_/ksize1;
    if (indb1 > ind1+1 && ind3 == ksize1-1)
	idx_assembly_++;

    return true;
}


// This function makes a check on whether two sub surfaces are
// neighbours, but if the surfaces meet at a singularity they are NOT
// classified as neighbours

// //===========================================================================
bool SurfaceAssembly::subSfNeighbour(int idx1, int idx2)
//===========================================================================
{
    int ksize1 = (int)(u_div_.size()) + 1;
    // int ksize2 = (int)(v_div_.size()) + 1;
    int ind1_1 = (idx1 % ksize1);
    int ind1_2 = idx1/ksize1;
    int ind2_1 = (idx2 % ksize1);
    int ind2_2 = idx2/ksize1;

    if (abs(ind2_1 - ind1_1) > 1 || abs(ind2_2 - ind1_2) > 1)
	return false;  // Not neighbours (@@@ closed surfaces not considered)

    // Check if the sub surfaces meet at a singularity
    int ki = -1, kj=-1;
    if (abs(ind1_1-ind2_1) == 1  && abs(ind1_2-ind2_2) == 1)
    {
	ki = std::max(ind1_1, ind2_1) - 1;
	kj = std::max(ind1_2, ind2_2) - 1;
	/*if (((u_div_[ki].second % 10) > 0 && (v_div_[kj].second % 10) > 0) ||
	    ((u_div_[ki].second % 10) > 0 && kj > 0 && (v_div_[kj-1].second % 10) > 0 ) ||
	    ((u_div_[ki].second % 10) > 0 && kj < ksize2-1 && (v_div_[kj+1].second % 10) > 0) ||
	    ((v_div_[kj].second % 10) > 0 && ki > 0 && (u_div_[ki-1].second % 10) > 0) ||
	    ((v_div_[kj].second % 10) > 0 && ki < ksize1-1 && (u_div_[ki+1].second % 10) > 0))
	    return false;   // Singularity*/
	if ((u_div_[ki].second % 10) == 0 && (v_div_[kj].second % 10) == 0)
	    return true;
    }
    else if (abs(ind1_1-ind2_1) == 1 && ind1_2 == ind2_2)
    {
	ki = std::min(ind1_1, ind2_1);
	kj = ind1_2;
	/*if ((u_div_[ki].second % 10) > 0 && 
	    ((kj > 0 && (v_div_[kj-1].second % 10) > 0) ||
	     (kj > 1 && (v_div_[kj-2].second % 10) > 0) ||
	     (kj < ksize2-1 && (v_div_[kj].second % 10) > 0) ||
	     (kj < ksize2-2 && (v_div_[kj+1].second % 10) > 0)))
	     return false;*/
	if ((u_div_[ki].second % 10) == 0)
	    return true;
    }
    else if (ind1_1 == ind2_1 && abs(ind1_2-ind2_2) == 1)
    {
	ki = ind1_1;
	kj = std::min(ind1_2, ind2_2);
	/*if ((v_div_[kj].second % 10) > 0 && 
	    ((ki > 0 && (u_div_[ki-1].second % 10) > 0) ||
	     (ki > 1 && (u_div_[ki-2].second % 10) > 0) ||
	     (ki < ksize1-1 && (u_div_[ki].second % 10) > 0) ||
	     (ki < ksize1-2 && (u_div_[ki+1].second % 10) > 0)))
	     return false;*/
	if ((v_div_[kj].second % 10) == 0)
	    return true;
    }

    return false;
/*    return true;  // Neighbours */
}


//===========================================================================
bool SurfaceAssembly::doTouch(int idx1, int idx2)
//===========================================================================
{
    int ksize1 = (int)(u_div_.size()) + 1;
    // int ksize2 = (int)(v_div_.size()) + 1;
    int ind1_1 = (idx1 % ksize1);
    int ind1_2 = idx1/ksize1;
    int ind2_1 = (idx2 % ksize1);
    int ind2_2 = idx2/ksize1;

    if (abs(ind2_1 - ind1_1) > 1 || abs(ind2_2 - ind1_2) > 1)
	return false;  // Not neighbours (@@@ closed surfaces not considered)
    else
	return true;
}

//===========================================================================
bool SurfaceAssembly::touchAtSingularity(int idx1, int idx2)
//===========================================================================
{
    int ksize1 = (int)(u_div_.size()) + 1;
    int ksize2 = (int)(v_div_.size()) + 1;
    int ind1_1 = (idx1 % ksize1);
    int ind1_2 = idx1/ksize1;
    int ind2_1 = (idx2 % ksize1);
    int ind2_2 = idx2/ksize1;

    if (abs(ind2_1 - ind1_1) > 1 || abs(ind2_2 - ind1_2) > 1)
	return false;  // Not neighbours (@@@ closed surfaces not considered)

    // Check if the sub surfaces meet at a singularity
    int ki = -1, kj=-1;
    if (abs(ind1_1-ind2_1) == 1  && abs(ind1_2-ind2_2) == 1)
    {
	ki = std::max(ind1_1, ind2_1) - 1;
	kj = std::max(ind1_2, ind2_2) - 1;
	if (((u_div_[ki].second % 10) > 0 && (v_div_[kj].second % 10) > 0) ||
	    ((u_div_[ki].second % 10) > 0 && kj > 0 && (v_div_[kj-1].second % 10) > 0) ||
	    ((u_div_[ki].second % 10) > 0 && kj < ksize2-1 && (v_div_[kj+1].second % 10) > 0) ||
	    ((v_div_[kj].second % 10) > 0 && ki > 0 && (u_div_[ki-1].second % 10) > 0) ||
	    ((v_div_[kj].second % 10) > 0 && ki < ksize1-1 && (u_div_[ki+1].second % 10) > 0))
	    return true;   // Singularity
    }
    else if (abs(ind1_1-ind2_1) == 1 && ind1_2 == ind2_2)
    {
	ki = std::min(ind1_1, ind2_1);
	kj = ind1_2;
	if ((u_div_[ki].second % 10) > 0 && 
	    ((kj > 0 && (v_div_[kj-1].second % 10) > 0) ||
	     (kj > 1 && (v_div_[kj-2].second % 10) > 0) ||
	     (kj < ksize2-1 && (v_div_[kj].second % 10) > 0) ||
	     (kj < ksize2-2 && (v_div_[kj+1].second % 10) > 0)))
	    return true;
    }
    else if (ind1_1 == ind2_1 && abs(ind1_2-ind2_2) == 1)
    {
	ki = ind1_1;
	kj = std::min(ind1_2, ind2_2);
	if ((v_div_[kj].second % 10) > 0 && 
	    ((ki > 0 && (u_div_[ki-1].second % 10) > 0) ||
	     (ki > 1 && (u_div_[ki-2].second % 10) > 0) ||
	     (ki < ksize1-1 && (u_div_[ki].second % 10) > 0) ||
	     (ki < ksize1-2 && (u_div_[ki+1].second % 10) > 0)))
	    return true;
    }

    return false;  // Neighbours
}

//===========================================================================
bool SurfaceAssembly::isInPrevAssembly(int idx1, int idx2)
//===========================================================================
{
    int ksize1 = (int)(u_div_.size()) + 1;
    int ksize2 = (int)(v_div_.size()) + 1;
    int ind1_1 = (idx1 % ksize1);
    // int ind1_2 = idx1/ksize1;
    int ind2_1 = (idx2 % ksize1);
    // int ind2_2 = idx2/ksize1;

    // Get next index of division value from a previous assembly in u-direction
    int next1 = -1, next2 = -1;
    int ki;
    for (ki=ind1_1; ki<ksize1-1; ki++)
	if (u_div_[ki].second >= 10)
	{
	    next1 = ki;
	    break;
	}

    for (ki=ind2_1; ki<ksize1-1; ki++)
	if (u_div_[ki].second >= 10)
	{
	    next2 = ki;
	    break;
	}

    if (next1 >=0 && next1 == next2)
	return true;  // These sub surfaces are already intersected

    // v-direction
    next1 = next2 = -1;
    for (ki=ind1_1; ki<ksize2-1; ki++)
	if (v_div_[ki].second >= 10)
	{
	    next1 = ki;
	    break;
	}

    for (ki=ind2_1; ki<ksize2-1; ki++)
	if (v_div_[ki].second >= 10)
	{
	    next2 = ki;
	    break;
	}

    if (next1 >=0 && next1 == next2)
	return true;  // These sub surfaces are already intersected

    return false;
}

//===========================================================================
bool SurfaceAssembly::isInFirstAssembly(shared_ptr<ParamSurfaceInt> sub_srf)
//===========================================================================
{
    double ta1 = surf_->startParam(0);
    double ta2 = surf_->startParam(1);
    double tb1 = surf_->endParam(0);
    double tb2 = surf_->endParam(1);
    double tc1 = sub_srf->startParam(0);
    double tc2 = sub_srf->startParam(1);
    double td1 = sub_srf->endParam(0);
    double td2 = sub_srf->endParam(1);
    int ksize1 = (int)(u_div_.size()) + 1;
    int ksize2 = (int)(v_div_.size()) + 1;

    if (tc1 < ta1-ptol_ || tc2 < ta2-ptol_ || td1 > tb1+ptol_ || td2 > tb2+ptol_)
	return false;  // Sub surface not inside

    // 1. parameter direction
    if (ksize1 > 1 && td1 > u_div_[1].first+ptol_)
	return false;

    if (ksize1 > 0 && td1 > u_div_[0].first+ptol_ && u_div_[0].second%10 > 0)
	return false;

    // 2. parameter direction
    if (ksize2 > 1 && td2 > v_div_[1].first+ptol_)
	return false;

    if (ksize2 > 0 && td2 > v_div_[0].first+ptol_ && v_div_[0].second%10 > 0)
	return false;

    return true;
}


//===========================================================================
int SurfaceAssembly::
getSubSurfaceIndex(shared_ptr<ParamSurfaceInt> sub_srf,
		   bool& at_end)
//===========================================================================
{
    double ta1 = surf_->startParam(0);
    double ta2 = surf_->startParam(1);
    double tb1 = surf_->endParam(0);
    double tb2 = surf_->endParam(1);
    double tc1 = sub_srf->startParam(0);
    double tc2 = sub_srf->startParam(1);
    double td1 = sub_srf->endParam(0);
    double td2 = sub_srf->endParam(1);

    if (tc1 < ta1-ptol_ || tc2 < ta2-ptol_ || td1 > tb1+ptol_ || td2 > tb2+ptol_)
	return -1;  // Sub surface not inside

    // 1. parameter direction
    double tstart, tend;
    size_t ki;
    int kn1=-1, kn2=-1;
    for (ki=0, tstart=ta1; ki<=(u_div_.size()); ki++, tstart=tend)
    {
	tend = (ki < u_div_.size()) ? u_div_[ki].first : tb1;
	if (fabs(tstart-tc1) < ptol_ && fabs(tend-td1) < ptol_)
	{
	    kn1 = (int)ki;
	    break;
	}
    }

    // 2. parameter direction
    for (ki=0, tstart=ta2; ki<=(v_div_.size()); ki++, tstart=tend)
    {
	tend = (ki < v_div_.size()) ? v_div_[ki].first : tb2;
	if (fabs(tstart-tc2) < ptol_ && fabs(tend-td2) < ptol_)
	{
	    kn2 = (int)ki;
	    break;
	}
    }

    if (kn1 == -1 || kn2 == -1)
	return -1;  // No index found

    at_end = (kn1 == (int)(u_div_.size()) || kn2 == (int)(v_div_.size()));
    return kn2*((int)u_div_.size()+1) + kn1;
}


//===========================================================================
void SurfaceAssembly::refineSurf()
//===========================================================================
{
    if (!surf_->isSpline()) {
	MESSAGE("Surface to refine was not a SplineSurface.  Aborting "
		"refinement procedure.");
    }
    shared_ptr<SplineSurface> splinesurf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_->getParamSurface());
    ALWAYS_ERROR_IF(!splinesurf.get(), "Failed cast in refineSurf.");
    
    BsplineBasis u_basis = splinesurf->basis(0);
    BsplineBasis v_basis = splinesurf->basis(1);
    const int u_order = u_basis.order();
    const int v_order = v_basis.order();
    
    // determining which knots to insert
    vector<double> new_u_knots, new_v_knots;
    int i;
    const double u_tol = ptol_ * (u_basis.endparam() - u_basis.startparam());
    for (i = 0; i < int(u_div_.size()); ++i) {
	double tmp = u_div_[i].first;
	u_basis.knotIntervalFuzzy(tmp, u_tol);
	int mult = u_order - u_basis.knotMultiplicity(tmp);
	new_u_knots.insert(new_u_knots.end(), mult, tmp);
    }
    const double v_tol = ptol_ * (v_basis.endparam() - v_basis.startparam());
    for (i = 0; i < int(v_div_.size()); ++i) {
	double tmp = v_div_[i].first;
	v_basis.knotIntervalFuzzy(tmp, v_tol);
	int mult = v_order - v_basis.knotMultiplicity(tmp);
	new_v_knots.insert(new_v_knots.end(), mult, tmp);
    }
    // inserting the knots
    splinesurf->insertKnot_u(new_u_knots);
    splinesurf->insertKnot_v(new_v_knots);
}


//===========================================================================
