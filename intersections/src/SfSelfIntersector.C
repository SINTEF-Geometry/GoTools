//===========================================================================
//                                                                           
// File: SfSfIntersector.C 
//                                                                           
// Created: 
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: SfSelfIntersector.C,v 1.31 2007-11-01 14:31:45 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/SfSelfIntersector.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/SfCvIntersector.h"
#include "GoTools/intersections/SfPtIntersector.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/SurfaceAssembly.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/LineCloud.h"

#include <vector>

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;
using std::ofstream;


//===========================================================================
SfSelfIntersector::SfSelfIntersector(shared_ptr<ParamSurfaceInt> surf,
				     double epsge,
				     Intersector *prev)
    : Intersector(epsge, prev)
//===========================================================================
{
    surf_ = surf;
    //epsge_ = shared_ptr<GeoTol>(new GeoTol(epsge));
    shared_ptr<IntersectionPool> parent_pool;
    if (prev) {
    parent_pool = prev->getIntPool();
    }
    int_results_ = 
	shared_ptr<IntersectionPool>(new IntersectionPool(surf, surf, 
							  parent_pool));
    max_rec_ = 4;   // Initial guess

    if (prev && prev->isSelfIntersection())
    {
	// Transfer information of singularity boxes and unions to this instance
	SfSelfIntersector *prev_int = dynamic_cast<SfSelfIntersector*>(prev);
	vector<double> mima = surf->getMima();  // Parameter domain of surface

	for (size_t ki=0; ki<prev_int->sing_union_.size(); ki++)
	    if (prev_int->sing_union_[ki].isInside(&mima[0]))
	    {
		SingUnion curr_union = prev_int->sing_union_[ki];
		for (int kj=0; kj<(int)(curr_union.singbox_idx_.size()); kj++)
		{
		    SingBox curr_box = prev_int->sing_box_[curr_union.singbox_idx_[kj]];
		    sing_box_.push_back(curr_box);
		    curr_union.singbox_idx_[kj] = kj;
		}
		sing_union_.push_back(curr_union);
	    }
    }
}


//===========================================================================
SfSelfIntersector::SfSelfIntersector(shared_ptr<ParamSurfaceInt> surf,
				     shared_ptr<GeoTol>  epsge,
				     Intersector *prev)
    : Intersector(epsge, prev)
//===========================================================================
{
    surf_ = surf;
    //epsge_ = shared_ptr<GeoTol>(new GeoTol(epsge.get()));
    shared_ptr<IntersectionPool> parent_pool;
    if (prev) {
    parent_pool = prev->getIntPool();
    }
    int_results_ = 
	shared_ptr<IntersectionPool>(new IntersectionPool(surf, surf, 
							  parent_pool));
    max_rec_ = 4;   // Initial guess

    if (prev && prev->isSelfIntersection())
    {
	// Transfer information of singularity boxes and unions to this instance
	SfSelfIntersector *prev_int = dynamic_cast<SfSelfIntersector*>(prev);
	vector<double> mima = surf->getMima();  // Parameter domain of surface

	for (size_t ki=0; ki<prev_int->sing_union_.size(); ki++)
	    if (prev_int->sing_union_[ki].isInside(&mima[0]))
	    {
		SingUnion curr_union = prev_int->sing_union_[ki];
		for (int kj=0; kj<(int)(curr_union.singbox_idx_.size()); kj++)
		{
		    SingBox curr_box = prev_int->sing_box_[curr_union.singbox_idx_[kj]];
		    sing_box_.push_back(curr_box);
		    curr_union.singbox_idx_[kj] = kj;
		}
		sing_union_.push_back(curr_union);
	    }
    }
}


//===========================================================================
SfSelfIntersector::~SfSelfIntersector()
//===========================================================================
{
  // Transfer complex domains to previous intersector
    if (prev_intersector_)
    {
	for (size_t ki=0; ki<complex_domain_.size(); ki++)
	    prev_intersector_->addComplexDomain(complex_domain_[ki]);
    }
}

//===========================================================================
void SfSelfIntersector::compute(bool compute_at_boundary)
//===========================================================================
{
    // Purpose: Compute topology of selfintersection results

    // First make a test to check if the current surface can
    // selfintersect at all
    if (!surf_->canSelfIntersect(epsge_->getEpsge()))
    {
	non_selfint_.push_back(surf_);
	return;
    }
    
    // Split the surface at G1 discontinuities and compute
    // selfintersections separate for each part of the surface
    vector<shared_ptr<ParamSurfaceInt> > subG1;
    surf_->splitAtG0(epsge_->getAngleTol(), subG1);

    // For each sub surface, compute selfintersections
    size_t ki, kj, kr, kh;
    vector<vector<shared_ptr<ParamSurfaceInt> > > nonself;
    nonself.resize(subG1.size());
    bool complex = false;   // Whether or not the sub-pieces are complex cases
    for (ki=0; ki<subG1.size(); ki++)
    {
	shared_ptr<SfSelfIntersector> 
	    sub(new SfSelfIntersector(subG1[ki], epsge_,
				      this));
	complex = sub->computeG1();  // Complexity does not occur at this level

	if (getenv("DO_REPAIR") && *(getenv("DO_REPAIR")) == '1') 
	{
	sub->repairIntersections();
	}

 	// Fetch all non-selfinterseting surfaces which are sub
	// surfaces of the sub surfaces. Collect them in the dedicated
	// class vector, but do also keep track on which item that
	// belongs to which sub surface
	vector<shared_ptr<ParamSurfaceInt> > sf_pure = 
	    sub->getNonSelfintersecting();
	for (kj=0; kj<sf_pure.size(); kj++)
	{
	    non_selfint_.push_back(sf_pure[kj]);
	    nonself[ki].push_back(sf_pure[kj]);
	}
   }

    // Compute intersections between sub surfaces
    for (ki=0; ki<subG1.size(); ki++)
	for (kj=ki+1; kj<subG1.size(); kj++)
	{
	    // Intersect only non selfintersectig surface. Thus, these
	    // surfaces have to be fetched from the collection belonging
	    // to each sub surface
	    for (kr=0; kr<nonself[ki].size(); kr++)
	    {
		for (kh=0; kh<nonself[kj].size(); kh++)
		{
		    shared_ptr<Intersector> 
			sfsfint(new SfSfIntersector(nonself[ki][kr],
						    nonself[kj][kh],
						    epsge_,
						    this));

		    // In this case neighbouring surfaces may be
		    // intersected.  Remove intersections at common
		    // boundaries if they do not belong to a curve
		    // sequence that leaves these boundaries. First
		    // check if the current surfaces have a common
		    // boundary

		    // It is also necessary to consider closed
		    // surfaces.  Check for extra intersection results
		    // also at the seem.

		    // Remove artificial intersection results
		    sfsfint->getIntPool()->removeBoundaryIntersections();
		}
	    }
	}

    // Prepare output intersection results
    if (prev_intersector_ == 0)
    {
	if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
	{
	    writeDebugComplex(1, complex_domain_);
	}

	if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
	{
	    int_results_->writeDebug();
	}

	// Top level intersector
	int_results_->cleanUpPool();

	// Remove loose ends of intersection links in the inner
	int_results_->weedOutClutterPoints();

	if (getenv("DO_REPAIR") && *(getenv("DO_REPAIR")) == '1') 
	{
	repairIntersections();
	}

	int_results_->makeIntersectionCurves();
	if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
	{
	    int_results_->writeDebug();
	}
    }
}

//===========================================================================
bool SfSelfIntersector::computeG1()
//===========================================================================
{
    // Purpose: Compute topology of selfintersection results of
    // surfaces which already is know to be G1

    // First make a test to check if the current surface can
    // selfintersect at all
    if (!surf_->canSelfIntersect(epsge_->getEpsge()))
    {
	non_selfint_.push_back(surf_);
	return false;
    }

    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
    {
	std::ofstream debug("selfint_out.g2");
	shared_ptr<ParamSurface> srf = surf_->getParamSurface();
	srf->writeStandardHeader(debug);
	srf->write(debug);
    }

    // Check recursion level
    int nmb_rec = nmbRecursions();

    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
    {
	std::cout << "Recursion level: " << nmb_rec << std::endl;
    }
    
    
    if (nmb_rec >= max_rec_)
    {
	// Connect intersections if possible

	// Store domain
	Vector2D lower(surf_->startParam(0),surf_->startParam(1));
	Vector2D upper(surf_->endParam(0),surf_->endParam(1));
	RectDomain dom(lower, upper);
	addComplexDomain(dom);

	if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
	{
	    std::cout <<"Complex domain: " << lower[0] << " ";
	    std::cout << upper[0] << " " << lower[1] << " " << upper[1];
	    std::cout << std::endl;
	}

	// Stop recursion
	return true;
    }
    
    
    // Define a surface split into a suitable number of pieces. Make
    // sure to split at singularities (points with vanishing surface
    // normal)
    bool intersection_computed = setSubdivision();
    if (intersection_computed)
	return false;  // Finished with the current sub surface

    // For each assembly of sub surface, compute selfintersections
    // Collect surfaces into assemblies, but do not make an assembly
    // that contains a singularity.  Closed surfaces implies an
    // assembly over the seem of the surface
    size_t kj;
    vector<vector<shared_ptr<ParamSurfaceInt> > > nonself;
    nonself.resize(div_sf_->getNmbSubSurface());
    div_sf_->resetAssemblyIndex();
    shared_ptr<ParamSurfaceInt> curr_assembly;
    int curr_idx, sub_sf_idx;
    int curr_sub = 0;
    bool at_end = false;
    bool sing_sub;
    bool complex_case = false;   // Whether or not the sub-pieces are complex cases
    while (div_sf_->getNextAssembly(curr_assembly, curr_idx, sing_sub))
    {

	double ta1 = curr_assembly->startParam(0);
	double tb1 = curr_assembly->endParam(0);
	double ta2 = curr_assembly->startParam(1);
	double tb2 = curr_assembly->endParam(1);
	if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
	{
	    cout << "Assembly: " << curr_idx << ", "
		 << ta1 << " " << tb1 << " "
		 << ta2 << " " << tb2 << endl;
	}

	if (sing_sub)
	{
	    // The surface assembly may have been handled already. 
	    // Check if the assembly contain a singularity
	    //if (hasSingularity(curr_assembly))  Not required any more (06.12)
	    //{
		if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
		{
		    cout << "Singularity block" << std::endl;
		}
		continue;
		//}
	}

	shared_ptr<SfSelfIntersector>
	    sub(new SfSelfIntersector(curr_assembly, epsge_, this));
	bool local_complex_case = sub->computeG1();
	if (local_complex_case)
	    complex_case = true;

	// Fetch all non-selfinterseting surfaces which are sub
	// surfaces of the sub surfaces. Collect them in the dedicated
	// class vector, but do also keep track on which item that
	// belongs to which sub surface
	vector<shared_ptr<ParamSurfaceInt> > sf_pure = 
	    sub->getNonSelfintersecting();

	// Avoid to store the same piece twice for overlapping
	// assemblies, and sort the piece with regard to sub surfaces
	for (kj=0; kj<sf_pure.size(); kj++)
	{
	    sub_sf_idx = div_sf_->getSubSurfaceIndex(sf_pure[kj], at_end);
	    if (sub_sf_idx != curr_sub && !at_end)
		continue;   // To avoid multible storage

	    non_selfint_.push_back(sf_pure[kj]); 
	    if (sub_sf_idx >= 0)
		nonself[sub_sf_idx].push_back(sf_pure[kj]);
	}
	curr_sub++;
	if (at_end)
	    curr_sub++;
    }

    if (complex_case && !hasSingularity(surf_))
    {
	// The maximum recursion level is reached, but self intersection
	// occuring at the current surface may not be computed. Compute the
	// self intersection by computing and connecting boundary intersection
	// at sub surface level
	// VSK, 070925. If it has a singularity, the sub surface has already been treated
	// with respect to self intersections
	shared_ptr<ParamSurfaceInt> curr_sub;
	int idx, sing_idx;
	vector<int> is_handled(4, 0);
	while (div_sf_->getNextSubSurface(curr_sub, idx, sing_idx))
	{
	    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
	    {
		double ta1 = curr_sub->startParam(0);
		double tb1 = curr_sub->endParam(0);
		double ta2 = curr_sub->startParam(1);
		double tb2 = curr_sub->endParam(1);
		cout << "Complex sub domain: " << curr_idx << ", "
		     << ta1 << " " << tb1 << " "
		     << ta2 << " " << tb2 << endl;
	    }

	    int in_complex_domain = isInComplexDomain(curr_sub);
	    if (in_complex_domain < 0)
	    {
		if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
		    std::cout << "Not included in the complex domain " << std::endl;
		continue;
	    }

	    if (div_sf_->isInPrevAssembly(idx, idx))
	    {
		if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
		    std::cout << "Is in prev assembly " << std::endl;

		continue;  // The current sub surface is already handled
	    }

	    // Check boundaries for self intersections and connect if possible
	    computeBoundaryIntersections(curr_sub, 0, &is_handled[0]);
	}
    }

    // Compute intersections between sub surfaces
    double eps2 = 0.5*epsge_->getEpsge();
    shared_ptr<ParamSurfaceInt> curr_sub1, curr_sub2;
    int idx1, idx2;
    int sing_idx1, sing_idx2;
    div_sf_->resetSubIndex();
    while (div_sf_->getNextSubSurface(curr_sub1, idx1, sing_idx1))
    {
	if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
	{
	    std::cout << std::endl;
	}
	div_sf_->resetSubIndex(idx1+1);
	while (div_sf_->getNextSubSurface(curr_sub2, idx2, sing_idx2))
	{
	    // DEBUG inf
	    double ta1 = curr_sub1->startParam(0);
	    double tb1 = curr_sub1->endParam(0);
	    double ta2 = curr_sub1->startParam(1);
	    double tb2 = curr_sub1->endParam(1);
	    double ta3 = curr_sub2->startParam(0);
	    double tb3 = curr_sub2->endParam(0);
	    double ta4 = curr_sub2->startParam(1);
	    double tb4 = curr_sub2->endParam(1);

	    // Check if the two subsurfaces are neighbours
	    // and do not meet at a singularity
	    int in_complex1 = isInComplexDomain(curr_sub1);
	    int in_complex2 = isInComplexDomain(curr_sub2);
	    if ((complex_case && in_complex1 >= 0 && in_complex2 >= 0) || 
		!div_sf_->subSfNeighbour(idx1, idx2))
	    {
		if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
		{
		    cout << "Sub surface: " << idx1 << ", "
			 << ta1 << " " << tb1 << " "
			 << ta2 << " " << tb2 << endl
			 << idx2 << ", "<< ta3 << " " << tb3 << " "
			 << ta4 << " " << tb4 << endl;
		}

		// Check if the two surface pieces are already treated
		// as a singularity block
		if (sing_idx1 >= 0 && sing_idx2 >= 0 && sing_idx1 == sing_idx2)
		{
		    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
		    {
			std::cout << "Sub surfaces in singularity block" << std::endl;
		    }
		    continue;
		}

		// Check if the two sub surfaces are already intersected in a previous
		// recursion
		bool is_in_prev = div_sf_->isInPrevAssembly(idx1, idx2);
		const SfSelfIntersector* prev = dynamic_cast<const SfSelfIntersector*>(prev_intersector_);
		bool is_in_first1=false, is_in_first2=false;
		if (complex_case && prev && prev->div_sf_.get())
		{
		    is_in_first1 = prev->div_sf_->isInFirstAssembly(curr_sub1);
		    is_in_first2 = prev->div_sf_->isInFirstAssembly(curr_sub2);
		}
		//if ( is_in_prev && !(complex_case && (is_in_first1 && is_in_first2)))
		//if ( is_in_prev && !(complex_case))
		if ( is_in_prev && !(complex_case && in_complex1 >= 0 && in_complex2 >= 0))
		    continue;
		    
		// Make an initial box test to intercept if the
		// subsurfaces cannot intersect
		CompositeBox box1 = curr_sub1->compositeBox();
		CompositeBox box2 = curr_sub2->compositeBox();
		if (!box1.overlaps(box2, eps2, eps2))
		    continue;

		bool touch = div_sf_->doTouch(idx1, idx2);

		if (getenv("WRITE_SING") && (*getenv("WRITE_SING"))=='1')
		{
		    std::ofstream debug1("sing_sf1.g2");
		    shared_ptr<ParamSurface> srf = curr_sub1->getParamSurface();
		    srf->writeStandardHeader(debug1);
		    srf->write(debug1);
		    std::ofstream debug2("sing_sf2.g2");
		    srf = curr_sub2->getParamSurface();
		    srf->writeStandardHeader(debug2);
		    srf->write(debug2);
		}

		/*if (sing_touch)
		{
		    // Check if the union of the cone corresponding to
		    // the two subsurfaces indicates that the combined
		    // surface may self intersect
		    DirectionCone cone1 = curr_sub1->directionCone();
		    DirectionCone cone2 = curr_sub2->directionCone();
		    cone1.addUnionWith(cone2);
		    if (!cone1.greaterThanPi())
			continue;  // No need to compute intersections
			}*/

		// Make sure that the choosen surfaces are either not
		// neighbours or contain a singularity
		shared_ptr<SfSfIntersector> 
		    sfsfint(new SfSfIntersector(curr_sub1, curr_sub2,
						epsge_, this));

		// Check if the two sub surfaces meet in a singularity
		double sing[4];
		if (shareSingularity(curr_sub1, curr_sub2, sing))
		    sfsfint->setHighPriSing(sing);

		if (touch)
		  sfsfint->setSelfintCase(2);

		if (getenv("DEBUG_DIV2") && (*getenv("DEBUG_DIV2"))=='1')
		{
		    sfsfint->getIntPool()->writeDebug();
		}

		/*if (sing_touch && ( sing_idx1 >= 0 || sing_idx2 >= 0))
		{
		    vector<pair<double, double> > max_curv1;
		    vector<pair<double, double> > max_curv2;
		    vector<pair<double, double> > max_curv3;
		    vector<pair<double, double> > max_curv4;

		    getMaxCurvatures(curr_sub1, 2, max_curv1, max_curv2);
		    getMaxCurvatures(curr_sub2, 2, max_curv3, max_curv4);

		    int stop_break;  // Debug functionality
		    stop_break = 1;
		    }*/

		// Compute intersections
		sfsfint->compute();

		// Remove loose ends of intersection links in the inner of the surfaces
		//sfsfint->getIntPool()->weedOutClutterPoints();

		if (getenv("DEBUG_DIV2") && (*getenv("DEBUG_DIV2"))=='1')
		{
		    sfsfint->getIntPool()->writeDebug();
		}

		// In this case neighbouring surfaces may be
		// intersected only if they touch a
		// singularity. Remove intersections at common
		// boundaries if they do not belong to a curve
		// sequence that leaves these boundaries.  First check
		// if the current surfaces have a common boundary

		// Remove artificial intersection results
		sfsfint->getIntPool()->removeBoundaryIntersections();

		int_results_->selfIntersectParamReorganise
		    (sfsfint->getIntPool());
		if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
		{
		    cout << "Nmb intpt: "
			 << int_results_->numIntersectionPoints() << endl;
		    if (getenv("DEBUG_DIV2") && (*getenv("DEBUG_DIV2"))=='1')
		    {
			int_results_->writeDebug();
		    }
		}
		if (touch)
		{
		    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
		    {
			std::cout << "Neighbour" << std::endl;
		    }
		}
	    }
	}

	div_sf_->resetSubIndex(idx1+1);
    }

    return false;
}


//===========================================================================
bool SfSelfIntersector::setSubdivision()
//===========================================================================
{
    // Purpose: Define splitting of the surface

    // Search for singularities (parameter values in which the surface
    // normal vanishes
    vector<shared_ptr<IntersectionPoint> > sing;
    bool computed = getSingularities(sing);
    bool performed = false;

    // Check relevance of singularities
    vector<pair<double,int> > u_div;
    vector<pair<double,int> > v_div;
    vector<RectDomain> sing_domain;



    if (computed)
    {
	ofstream debug1("init_single_block.dsp");

	// Store singularity information
	// Get normal surface
	shared_ptr<ParamSurfaceInt> normsf = surf_->getNormalSurface();
	singularity_info_->addIterationCount();
	double param[2];
	for (size_t ki=0; ki<sing.size(); ki++)
	{
	    param[0] = sing[ki]->getPar(0);
	    param[1] = sing[ki]->getPar(1);
	    singularity_info_->addSingularPoint(param, 2);

	    // Estimate size of singularity box
	    vector<pair<double, int> > sing_box(4);   // umin, umax, vmin, vmax
	    estimateSingBox(normsf, param, sing_box);

	    SingBox curr_box(sing_box, sing[ki]);
	    sing_box_.push_back(curr_box);

	    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
	    {
		if (sing_box[2].second == 1)
		    debug1 << "fg: magenta" << endl;
		else if (sing_box[2].second == 2)
		    debug1 << "fg: green" << endl;
		else 
		    debug1 << "fg: orange" << endl;
		debug1 << "lin:" << endl;
		debug1 << sing_box[0].first << " " << sing_box[2].first << endl;
		debug1 << sing_box[1].first << " " << sing_box[2].first << endl;

		if (sing_box[1].second == 1)
		    debug1 << "fg: magenta" << endl;
		else if (sing_box[1].second == 2)
		    debug1 << "fg: green" << endl;
		else 
		    debug1 << "fg: orange" << endl;
		debug1 << "lin:" << endl;
		debug1 << sing_box[1].first << " " << sing_box[2].first << endl;
		debug1 << sing_box[1].first << " " << sing_box[3].first << endl;

		if (sing_box[3].second == 1)
		    debug1 << "fg: magenta" << endl;
		else if (sing_box[3].second == 2)
		    debug1 << "fg: green" << endl;
		else 
		    debug1 << "fg: orange" << endl;
		debug1 << "lin:" << endl;
		debug1 << sing_box[1].first << " " << sing_box[3].first << endl;
		debug1 << sing_box[0].first << " " << sing_box[3].first << endl;

		if (sing_box[0].second == 1)
		    debug1 << "fg: magenta" << endl;
		else if (sing_box[3].second == 2)
		    debug1 << "fg: green" << endl;
		else 
		    debug1 << "fg: orange" << endl;
		debug1 << "lin:" << endl;
		debug1 << sing_box[0].first << " " << sing_box[3].first << endl;
		debug1 << sing_box[0].first << " " << sing_box[2].first << endl;
	    }

	
	}
    }

    makeSingularityUnions();

    // Define division. First check if the domain must be split according to unions
    splitUnions(u_div, v_div);
    
    if (u_div.size() == 0 && v_div.size() == 0)
    {
	if (sing_union_.size() == 1 && sing_box_.size() > 1)
	{
	    divideOneUnion(u_div, v_div);
	}
	else
	{
	    sortSingularities(sing, u_div, v_div);
	}


	// Modify division values to take possible knot line into consideration
	// and avoid subdivision close to these
	modifyDivisionValues(sing, u_div, 0);
	modifyDivisionValues(sing, v_div, 1);

	// Refine subdivision curves where it cross division lines
	refineSingularityLinks(sing, u_div, 0);
	refineSingularityLinks(sing, v_div, 1);

	// Traverse the elements in the given grid, and define topology if required
	performed = defineSingTopology(sing, u_div, v_div, sing_domain);

	if (!performed)
	{
	    // Add already existing division parameters from the previous recursion
	    getPrevDivision(u_div, v_div);

	    // Get the number of divisions in the two parameter directions
	    int nmb_u, nmb_v;
	    getNmbDivision(nmb_u, nmb_v);

	    // Set splitting parameters
	    setDivisionValues(u_div, nmb_u, 0);
	    setDivisionValues(v_div, nmb_v, 1);

	    // Modify division values to take possible knot line into consideration
	    // and avoid subdivision close to these
	    modifyDivisionValues(sing, u_div, 0);
	    modifyDivisionValues(sing, v_div, 1);

	}
    }
    else
    {
	// Modify division values to take possible knot line into consideration
	// and avoid subdivision close to these
	modifyDivisionValues(sing, u_div, 0);
	modifyDivisionValues(sing, v_div, 1);
    }

//    vector<pair<double,double> > max_curv1, max_curv2;
//    getMaxCurvatures(max_curv1, max_curv2);

    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
    {
	std::cout << "Nmb div u: " << u_div.size() << std::endl;
	for (size_t ki=0; ki<u_div.size(); ki++)
	    std::cout << "(" << u_div[ki].first
		      << "," << u_div[ki].second << ") ";
	std::cout << std::endl;
	std::cout << "Nmb div v: " << v_div.size() << std::endl;
	for (size_t ki=0; ki<v_div.size(); ki++)
	    std::cout << "(" << v_div[ki].first
		      << "," << v_div[ki].second << ") ";
	std::cout << std::endl;

	writeDebugComplex(2, sing_domain);
    }

    div_sf_ = (shared_ptr<SurfaceAssembly>)
	(new SurfaceAssembly(surf_, u_div, v_div, sing_domain, epsge_->getRelParRes()));

    return performed;
}


//===========================================================================
bool SfSelfIntersector::
getSingularities(vector<shared_ptr<IntersectionPoint> >& sing)
//===========================================================================
{
    // Purpose: Search for parameter values in which the surface
    // normal vanishes

    // Parameter domain
    vector<double> mima = surf_->getMima();
    //double ptol = epsge_->getRelParRes();

    // Fetch singularity information from the previous intersector if
    // necessary and if such an information exist
    if (singularity_info_.get() == 0 && prev_intersector_ &&
	prev_intersector_->hasSingularityInfo())
    {
	singularity_info_ = (shared_ptr<SingularityInfo>)
	    (new SingularityInfo(prev_intersector_->getSingularityInfo(),
				 true));
    }
    else if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_
	    = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }

    // Check if a singularity search is performed already
    if (singularity_info_->iterationDone())
    {
	// Fetch singularities
	for (size_t kr=0; kr<sing_box_.size(); kr++)
	    sing.push_back(sing_box_[kr].sing_);
// 	// Fetch singularities
// 	if (singularity_info_->hasPoint())
// 	{
// 	    int nmb_point = singularity_info_->getNmbPoints(2);
// 	    for (int kr=0; kr<nmb_point; kr++)
// 	    {
// 		Point sing_pt = singularity_info_->getPoint(kr, 2);
// 		if (sing_pt[0] > mima[0]-ptol && sing_pt[0] < mima[1]+ptol &&
// 		    sing_pt[1] > mima[2]-ptol && sing_pt[1] < mima[2]+ptol)
// 		sing.push_back(sing_pt);
// 	    }
// 	}
	return false;   // No point in searching for singularities more than once
    }



    // Define origo
    int dim = surf_->dimension();
    shared_ptr<Point> origo(new Point(dim));
    origo->setValue(0.0);
    shared_ptr<ParamPointInt> origo_int(new ParamPointInt(origo));

    // Compute normal surface
    shared_ptr<ParamSurfaceInt> normsf = surf_->getNormalSurface();
    if (normsf.get() == 0)
	return false;  // @@@ Not a spline surface
    
    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
    {
	std::ofstream debug("normsf_out.g2");
	shared_ptr<ParamSurface> srf = normsf->getParamSurface();
	srf->writeStandardHeader(debug);
	srf->write(debug);
    }

     // Intersect
    shared_ptr<SfPtIntersector>
	sfptint(new SfPtIntersector(normsf, origo_int, epsge_));
    sfptint->setSelfintCase(1);

    sfptint->compute();
							
    vector<shared_ptr<IntersectionPoint> > intpts;
    vector<shared_ptr<IntersectionCurve> > intcrv;
    sfptint->getResult(intpts, intcrv);

    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
    {
	sfptint->getIntPool()->writeDebug(1);
	std::ofstream outf("sing_pnt.g2");
	int ki, kj;
	std::cout << "Number of singular points: " << intpts.size() << std::endl;
	for (ki=0; ki<int(intpts.size()); ki++)
	{
	    std::vector<IntersectionPoint*> neighbours;
	    intpts[ki]->getNeighbours(neighbours);
	    std::vector<double> par = intpts[ki]->getPar();
	    for (kj=0; kj<int(par.size()); kj++)
		std::cout << par[kj] << " ";
	    std::cout << neighbours.size() << "  ";
	    std::cout << intpts[ki]->getDist() << std::endl;

	    vector<Point> pt1(6);
	    vector<Point> pt2(3);
	    double par1[2];
	    par1[0] = par[0];
	    par1[1] = par[1];
	    surf_->point(pt1, par1, 2);
	    normsf->point(pt2, par1, 1);
	    int j1, j2;
	    for (j1=0; j1<6; j1++)
		for (j2=0; j2<3; j2++)
		    std::cout << pt1[j1][j2] << " ";
	    std::cout << std::endl;
	    for (j1=0; j1<3; j1++)
		for (j2=0; j2<3; j2++)
		    std::cout << pt2[j1][j2] << " ";
	    double angle = pt2[1].angle(pt2[2]);
	    std::cout << angle << std::endl;
	    std::cout << std::endl;

	    outf << "400 1 0 4 255 0 0 255" << std::endl;
	    outf << "1" << std::endl;
	    pt1[0].write(outf);
	    outf << "\n";

	}

	std::ofstream outg("sing_crv.g2");
	for (ki=0; ki<int(intcrv.size()); ki++)
	{
	    std::vector<double> guide_pt;
	    int nguide = intcrv[ki]->numGuidePoints(); 
	    std::cout << "Number of guide points: " << nguide << std::endl;
	    guide_pt.reserve(3*nguide);
	    for (kj=0; kj<nguide; kj++)
	    {
		shared_ptr<IntersectionPoint>
		    currpt = intcrv[ki]->getGuidePoint(kj);
		std::vector<double> par = currpt->getPar();
		vector<Point> pt1(6);
		vector<Point> pt2(3);
		double par1[2];
		par1[0] = par[0];
		par1[1] = par[1];
		surf_->point(pt1, par1, 2);
		normsf->point(pt2, par1, 1);
		int j1, j2;
		for (j1=0; j1<6; j1++)
		    for (j2=0; j2<3; j2++)
			std::cout << pt1[j1][j2] << " ";
		std::cout << std::endl;
		for (j1=0; j1<3; j1++)
		    for (j2=0; j2<3; j2++)
			std::cout << pt2[j1][j2] << " ";
		double angle = pt2[1].angle(pt2[2]);
		std::cout << angle << std::endl;
		std::cout << std::endl;
		guide_pt.insert(guide_pt.end(), pt1[0].begin(), pt1[0].end());
		if (kj>0 && kj<nguide-1)
		    guide_pt.insert(guide_pt.end(), pt1[0].begin(), pt1[0].end()); 
	    }
	    LineCloud line(&guide_pt[0], nguide-1);
	    line.writeStandardHeader(outg);
	    line.write(outg);
	}

    }

    // Transform the singularities to the current intersection pool.
    // First fetch all intersection points
    size_t ki, kj, kr;
    int idx;
    double param[4];
    IntersectionPoint *p1, *p2, *pcurr;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    sfptint->getIntPool()->getIntersectionPoints(int_pts);

    // Simplify singularity result by removing inner points on
    // intersection links that lie close to at least one of its
    // neighbours
    double frac = 0.01;
    double closetol = frac*sqrt((mima[1]-mima[0])*(mima[1]-mima[0]) +
				(mima[3]-mima[2])*(mima[3]-mima[2]));
    double ang_tol = 0.2;
    double pdist, ang;
    for (idx=0; idx<(int)(int_pts.size()); idx++)
    {
	vector<shared_ptr<IntersectionLink> > links;
	int_pts[idx]->getNeighbourLinks(links);
	if (links.size() == 2)
	{
	    // Let pcurr point to the current point, and p1 and p2 to
	    // the two neighbours
	    pcurr = int_pts[idx].get();
	    p1 = links[0]->getOtherPoint(pcurr);
	    p2 = links[1]->getOtherPoint(pcurr);
	    
	    // Remove only if pcurr lies in the rectangle defined by
	    // p1 and p2
	    if (pcurr->getPar(0) < std::min(p1->getPar(0), p2->getPar(0)) ||
		pcurr->getPar(0) > std::max(p1->getPar(0), p2->getPar(0)))
		continue;
	    if (pcurr->getPar(1) < std::min(p1->getPar(1), p2->getPar(1)) ||
		pcurr->getPar(1) > std::max(p1->getPar(1), p2->getPar(1)))
		continue;

	    // Check also angle
	    Point vec1(pcurr->getPar(0)-p1->getPar(0), pcurr->getPar(1)-p1->getPar(1));
	    Point vec2(p2->getPar(0)-pcurr->getPar(0), p2->getPar(1)-pcurr->getPar(1));
	    ang = vec1.angle(vec2);
	    if (ang > ang_tol)
		continue;

	    for (kj=0; kj<links.size(); kj++)
	    {
		// Compute distance to each neighbour
		links[kj]->getIntersectionPoints(p1, p2);
		pdist = sqrt((p1->getPar(0)-p2->getPar(0))*
			     (p1->getPar(0)-p2->getPar(0)) +
			     (p1->getPar(1)-p2->getPar(1))*
			     (p1->getPar(1)-p2->getPar(1)));

		if (pdist < closetol)
		{
		    // Remove the point from the pool to keep
		    // connectivity information
		    sfptint->getIntPool()->removeIntPoint(int_pts[idx]);
		    int_pts.erase(int_pts.begin()+idx,
				  int_pts.begin()+idx+1);
		    break;
		}
	    }
	    if (kj < links.size())
		idx = -1;
	}
    }
		
    //ofstream debug1("init_single_block.dsp");

    
    vector<shared_ptr<IntersectionPoint> > int_pts2(int_pts.size());
    for (ki=0; ki<int_pts.size(); ki++)
    {
	param[0] = param[2] = int_pts[ki]->getPar(0);
	param[1] = param[3] = int_pts[ki]->getPar(1);

	int_pts2[ki]
	    = int_results_->addIntersectionPoint(surf_, surf_, epsge_,
						 param, param+2);
	
	// Check accuracy
	if (int_pts2[ki]->getDist() > epsge_->getEpsge())
	{
	    // Remove the singular point. It is not an intersection point
	    int_results_->removeIntPoint(int_pts2[ki]);
	    continue;
	}
	    
	// Transform intersection links to already stored points
	vector<shared_ptr<IntersectionLink> > links;
	int_pts[ki]->getNeighbourLinks(links);
	for (kj=0; kj<links.size(); kj++)
	{
	    // Get the other point
	    p1 = links[kj]->getOtherPoint(int_pts[ki].get());
	    for (kr=0; kr<ki; kr++)
	    {
		if (p1 == int_pts[kr].get())
		{
		    int_pts2[ki]->connectTo(int_pts2[kr], LINK_UNDEFINED);
		    break;
		}
	    }
	}
	sing.push_back(int_pts2[ki]);

    }

    return true;
}


bool compare_divpar(const pair<double, int>& l1, const pair<double,int>& l2)
{ 
    return (l1.first < l2.first);
}

//===========================================================================
void SfSelfIntersector::
getPrevDivision(vector<pair<double,int> >& divpar_u,
		vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    const SfSelfIntersector* prev = dynamic_cast<const SfSelfIntersector*>(prev_intersector_);
    vector<double> mima = surf_->getMima();

    if (prev == 0 || prev->div_sf_.get() == 0)
	return;  // No previous division lines

    // Get division values or previous recursion
    vector<pair<double,int> > u_prev = prev->div_sf_->getUdiv();
    vector<pair<double,int> > v_prev = prev->div_sf_->getVdiv();

    // Add relevant division values
    for (size_t ki=0; ki<u_prev.size(); ki++)
    {
	if (u_prev[ki].first > mima[0] && u_prev[ki].first < mima[1])
	{
	    divpar_u.push_back(u_prev[ki]);
	    divpar_u[divpar_u.size()-1].second += 10;
	}
    }

    for (size_t ki=0; ki<v_prev.size(); ki++)
    {
	if (v_prev[ki].first > mima[2] && v_prev[ki].first < mima[3])
	{
	    divpar_v.push_back(v_prev[ki]);
	    divpar_v[divpar_v.size()-1].second += 10;
	}
    }

    // Sort
    std::sort(divpar_u.begin(), divpar_u.end(), compare_divpar);
    std::sort(divpar_v.begin(), divpar_v.end(), compare_divpar);
}

//===========================================================================
void SfSelfIntersector::
sortSingularities(vector<shared_ptr<IntersectionPoint> >& sing,
		  vector<pair<double,int> >& divpar_u,
		  vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    // Parameter domain
    vector<double> mima = surf_->getMima();
    double ptol = 0.0000001;
    ptol = std::min(ptol, 0.01*(mima[1]-mima[0]));
    ptol = std::min(ptol, 0.01*(mima[3]-mima[2]));
    static double frac = 0.1;

    // Dismiss not connected, not standalone singularities that do not
    // belong to a list
    //validateSingularities(sing);  // VSK, 0612. Is this still required?
    if (sing.size() == 0)
	return;   // No further singularity discussions are required

    // Make a grid through all singularities
    vector<double> u_sing(sing.size()), v_sing(sing.size());
    for (size_t ki=0; ki<sing.size(); ki++)
    {
	u_sing[ki] = sing[ki]->getPar(0);
	v_sing[ki] = sing[ki]->getPar(1);
    }

    // Sort 
    std::sort(u_sing.begin(), u_sing.end());
    std::sort(v_sing.begin(), v_sing.end());

    // Adapt the size of the parameter tolerance to the size of the
    // singularity area
    double tdel1 = u_sing[u_sing.size()-1] - u_sing[0];
    double tdel2 = v_sing[v_sing.size()-1] - v_sing[0];
    ptol = std::max(ptol, 0.01*std::min(tdel1,tdel2));
    

    // Dismiss too close parameter values
    for (size_t ki=1; ki<u_sing.size(); ki++)
    {
	if (u_sing[ki] - u_sing[ki-1] < ptol)
	    u_sing.erase(u_sing.begin()+ki, u_sing.begin()+ki+1);
    }
    for (size_t ki=1; ki<v_sing.size(); ki++)
    {
	if (v_sing[ki] - v_sing[ki-1] < ptol)
	    v_sing.erase(v_sing.begin()+ki, v_sing.begin()+ki+1);
    }



    if (sing_box_.size() == 1)
    {
	// Use existing box
	if (sing_box_[0].box_limit_[0].first - mima[0] <
	    frac*(sing_box_[0].box_limit_[1].first-sing_box_[0].box_limit_[0].first))
	    divpar_u.push_back(make_pair(mima[0], 2));
	else
	    divpar_u.push_back(make_pair(sing_box_[0].box_limit_[0].first, 2));

	if (mima[1] - sing_box_[0].box_limit_[1].first <
	    frac*(sing_box_[0].box_limit_[1].first-sing_box_[0].box_limit_[0].first))
	    divpar_u.push_back(make_pair(mima[1], 2));
	else
	    divpar_u.push_back(make_pair(sing_box_[0].box_limit_[1].first, 2));

	if (sing_box_[0].box_limit_[2].first - mima[2] <
	    frac*(sing_box_[0].box_limit_[3].first-sing_box_[0].box_limit_[2].first))
	    divpar_v.push_back(make_pair(mima[2], 2));
	else
	    divpar_v.push_back(make_pair(sing_box_[0].box_limit_[2].first, 2));

	if (mima[3] - sing_box_[0].box_limit_[2].first  <
	    frac*(sing_box_[0].box_limit_[3].first-sing_box_[0].box_limit_[2].first))
	    divpar_v.push_back(make_pair(mima[3], 2));
	else
	    divpar_v.push_back(make_pair(sing_box_[0].box_limit_[3].first, 2));
    }
    else
    {

	// Make grid lines between singularities, but stand alone
	// singularities with large distance to the closest singularity
	// may define the position of a grid line.
	makeDivisionLines(sing, u_sing, 0, mima[0], mima[1], divpar_u);
	makeDivisionLines(sing, v_sing, 1, mima[2], mima[3], divpar_v);
    }

    // Check if all division lines are required
    //checkDivisionLines(sing, divpar_u, divpar_v);
}
	

//===========================================================================
bool SfSelfIntersector::defineSingTopology(vector<shared_ptr<IntersectionPoint> >& sing,
					   vector<pair<double,int> >& divpar_u,
					   vector<pair<double,int> >& divpar_v,
					   vector<RectDomain>& sing_domain)
//===========================================================================
{
    // Traverse the elements in the given grid, and check if the
    // topology should be defined at this stage
    double t1, t2;
    double ta1, ta2, tb1, tb2;  // Limits of current sub surface domain
    vector<double> mima = surf_->getMima();  // Parameter domain of surface
    double pres = epsge_->getRelParRes();

    vector<pair<int,int> > sing_cells;

    for (size_t idx2=1; idx2<divpar_v.size(); idx2++)
    {
	ta2 = divpar_v[idx2-1].first;
	tb2 = divpar_v[idx2].first;
	for (size_t idx1=1; idx1<divpar_u.size(); idx1++)
	{
	    vector<shared_ptr<IntersectionPoint> > inside_sing;
	    ta1 = divpar_u[idx1-1].first;
	    tb1 = divpar_u[idx1].first;

	    // Check if the current box contain a singularity
	    for (size_t ki=0; ki<sing.size(); ki++)
	    {
		t1 = sing[ki]->getPar(0);
		t2 = sing[ki]->getPar(1);
		if (t1 > ta1+pres && 
		    t1 < tb1-pres &&
		    t2 > ta2+pres && 
		    t2 < tb2-pres)
		    inside_sing.push_back(sing[ki]);
	    }
	    if (inside_sing.size() == 0)
		continue;   // No singular point

	    // Choose the middle singular point
	    int idx_sing = 0;
	    double min_dist = 1.0e6, dist;
	    Point mid_pt(0.5*(ta1+tb1),0.5*(ta2+tb2));
	    for (size_t ki=0; ki<inside_sing.size(); ki++)
	    {
		t1 = inside_sing[ki]->getPar(0);
		t2 = inside_sing[ki]->getPar(1);
		Point curr(t1, t2);
		dist = mid_pt.dist(curr);
		if (dist < min_dist)
		{
		    idx_sing = (int)(ki);
		    min_dist = dist;
		}
	    }

	    // Check if any of the boundaries are intersected already
	    int is_handled[4];
	    is_handled[0] = is_handled[1] = is_handled[2] = is_handled[3] = 0;

	    // Run through the previous handled singularity cells
	    for (size_t ki=0; ki<sing_cells.size(); ki++)
	    {
		if ((int)idx1 == sing_cells[ki].first &&
		    (int)idx2 == sing_cells[ki].second+1)
		    is_handled[2] = 1;
		if ((int)idx1 == sing_cells[ki].first+1 &&
		    (int)idx2 == sing_cells[ki].second)
		    is_handled[0] = 1;
	    }
		

	    // Compute selfintersection points at the boundaries of
	    // the current subsurface
	    vector<shared_ptr<ParamSurfaceInt> > sub_sfs =
		surf_->subSurfaces(ta1, ta2, tb1, tb2, pres);

	    for (size_t ki=0; ki<sub_sfs.size(); ki++) {
	      computeBoundaryIntersections(sub_sfs[ki], 
					   inside_sing[idx_sing].get(),
					   is_handled);
	    }

	    // Indicate that this cell has been subject to boundary intersection

	    // Store domain
	    Vector2D lower(ta1, ta2);
	    Vector2D upper(tb1, tb2);
	    RectDomain dom(lower, upper);
	    sing_domain.push_back(dom);

	    // Check if the singular domain is equal to the current sub surface.
	    // In that case, we are finished.
	    if (fabs(ta1-mima[0]) < pres && fabs(tb1-mima[1]) < pres &&
		fabs(ta2-mima[2]) < pres && fabs(tb2-mima[3]) < pres)
		return true;
	    
	}	
    }   
    return false;
}

bool compare(const pair<double, double>& l1, const pair<double,double>& l2)
{ 
    if (l1.first < l2.first)
	return true;
    else if (l2.first < l1.first)
	return false;
    else
    return (l1.second < l2.second);
}

//===========================================================================
void SfSelfIntersector::
splitUnions(vector<pair<double,int> >& divpar_u,
	    vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    if (sing_union_.size() == 0)
	return;

    static double frac = 0.5;
    double ptol0 = epsge_->getRelParRes();
    double ptol = 100.0*ptol0;

    // U-direction
    double start = surf_->startParam(0);
    double end = surf_->endParam(0);
    double tlast;

    // Fetch union limits
    vector<pair<double,double> > u_limit;
    double meanlength = 0.0, length;
    for (size_t ki=0; ki<sing_union_.size(); ki++)
    {
	u_limit.push_back(make_pair(sing_union_[ki].limit_[0],sing_union_[ki].limit_[1]));
	length = sing_union_[ki].limit_[1] - sing_union_[ki].limit_[0];
	meanlength += length;
    }
    meanlength /= (double)(sing_union_.size());

    // Sort
    std::sort(u_limit.begin(), u_limit.end(), compare);

    if (u_limit[0].first - start > frac*meanlength)
	divpar_u.push_back(make_pair(u_limit[0].first, 3));
    for (size_t ki=0; ki<u_limit.size()-1; ki++)
    {
	if (u_limit[ki+1].first < u_limit[ki].second && 
	    u_limit[ki].second <= u_limit[ki+1].second+ptol0);  // Overlapping unions in this direction. 
	                                                       //Do not split.
	else if (u_limit[ki].second > u_limit[ki+1].second)
	{
	    std::swap(u_limit[ki], u_limit[ki+1]);
	}
	else if (u_limit[ki+1].first > u_limit[ki].first &&
		 u_limit[ki+1].first < u_limit[ki].second)
	    divpar_u.push_back(make_pair(u_limit[ki].second, 3));

	else if (u_limit[ki+1].first - u_limit[ki].second < frac*meanlength)
	    // Define split at the average position between the singularity unions
	    divpar_u.push_back(make_pair(0.5*(u_limit[ki].second+u_limit[ki+1].first), 3));
	else
	{
	    // Split at both union limits
	    divpar_u.push_back(make_pair(u_limit[ki].second, 3));
	    divpar_u.push_back(make_pair(u_limit[ki+1].first, 3));
	}
    }
    
    tlast = u_limit[0].second;
    for (size_t ki=1; ki<u_limit.size(); ki++)
	tlast = std::max(tlast, u_limit[ki].second);
    if (end - tlast > frac*meanlength)
	divpar_u.push_back(make_pair(tlast, 3));
    
    if (divpar_u.size()>0 && 
	divpar_u[divpar_u.size()-1].first > end-ptol)
	divpar_u.erase(divpar_u.end()-1, divpar_u.end());
    
    
    // V-direction
    start = surf_->startParam(1);
    end = surf_->endParam(1);

    // Fetch union limits
    vector<pair<double,double> > v_limit;
    meanlength = 0.0;
    for (size_t ki=0; ki<sing_union_.size(); ki++)
    {
	v_limit.push_back(make_pair(sing_union_[ki].limit_[2],sing_union_[ki].limit_[3]));
	length = sing_union_[ki].limit_[3] - sing_union_[ki].limit_[2];
	meanlength += length;
    }
    meanlength /= (double)(sing_union_.size());

    // Sort
    std::sort(v_limit.begin(), v_limit.end(), compare);

    if (v_limit[0].first - start > frac*meanlength)
	divpar_v.push_back(make_pair(v_limit[0].first, 3));
    for (size_t ki=0; ki<v_limit.size()-1; ki++)
    {
	if (v_limit[ki+1].first < v_limit[ki].second &&
	    v_limit[ki].second <= v_limit[ki+1].second+ptol);  // Overlapping unions in this direction. 
	                                                       // Do not split.
	else if (v_limit[ki].second > v_limit[ki+1].second)
	{
	    std::swap(v_limit[ki], v_limit[ki+1]);
	}
	else if (v_limit[ki+1].first > v_limit[ki].first &&
		 v_limit[ki+1].first < v_limit[ki].second)
	    divpar_v.push_back(make_pair(v_limit[ki].second, 3));

	else if (v_limit[ki+1].first - v_limit[ki].second < frac*meanlength)
	    // Define split at the average position between the singularity unions
	    divpar_v.push_back(make_pair(0.5*(v_limit[ki].second+v_limit[ki+1].first), 3));
	else
	{
	    // Split at both union limits
	    divpar_v.push_back(make_pair(v_limit[ki].second, 3));
	    divpar_v.push_back(make_pair(v_limit[ki+1].first, 3));
	}
    }
    tlast = v_limit[0].second;
    for (size_t ki=1; ki<v_limit.size(); ki++)
	tlast = std::max(tlast, v_limit[ki].second);
    if (end - tlast > frac*meanlength)
	divpar_v.push_back(make_pair(tlast, 3));
    
    if (divpar_v.size()>0 && 
	divpar_v[divpar_v.size()-1].first > end-ptol)
	divpar_v.erase(divpar_v.end()-1, divpar_v.end());
    
    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
    {
	vector<double> mima = surf_->getMima();
	ofstream debug1("union_div.dsp");
	debug1 << "fg: magenta" << endl;
	for (size_t ki=0; ki<divpar_u.size(); ki++)
	{
	    debug1 << "lin:" << endl;
	    debug1 << divpar_u[ki].first << "  " << mima[2] << endl;
	    debug1 << divpar_u[ki].first << "  " << mima[3] << endl;
	}
	for (size_t ki=0; ki<divpar_v.size(); ki++)
	{
	    debug1 << "lin:" << endl;
	    debug1 << mima[0] << "  " << divpar_v[ki].first << endl;
	    debug1 << mima[1] << "  " << divpar_v[ki].first << endl;
	}
    }
    
}

bool compare_box1(const SingBox& l1, const SingBox& l2)
{ 
    return (l1.sing_->getPar(0) < l2.sing_->getPar(0));
}

bool compare_box2(const SingBox& l1, const SingBox& l2)
{ 
    return (l1.sing_->getPar(1) < l2.sing_->getPar(1));
}

//===========================================================================
void SfSelfIntersector::
divideOneUnion(vector<pair<double,int> >& divpar_u,
	       vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    if (sing_union_.size() != 1)
	return;

    vector<double> mima = surf_->getMima();
    double ptol = 0.0000001;
    double ptol1 = std::max(ptol, 0.01*(mima[1]-mima[0]));
    double ptol2 = std::max(ptol, 0.01*(mima[3]-mima[2]));

    // Sort singularity boxes
    if (mima[1]-mima[0] < mima[3]-mima[2])
	std::sort(sing_box_.begin(), sing_box_.end(), compare_box1);
    else
	std::sort(sing_box_.begin(), sing_box_.end(), compare_box2);

    // Make sure that all singularities are separated by at least one
    // division line
    size_t ki, kj;
    for (ki=0; ki<sing_box_.size(); ki++)
	for (kj=ki+1; kj<sing_box_.size(); kj++)
	{
	    // Check if the singularities are separated already
	    if (isSeparated(sing_box_[ki].sing_, sing_box_[kj].sing_, divpar_u, divpar_v))
		continue;

	    // Define candidate division value. First check overlap in 
	    // both parameter directions
	    double ta1 = sing_box_[ki].box_limit_[0].first;
	    double tb1 = sing_box_[ki].box_limit_[1].first;
	    double tc1 = sing_box_[ki].box_limit_[2].first;
	    double td1 = sing_box_[ki].box_limit_[3].first;
	    double ta2 = sing_box_[kj].box_limit_[0].first;
	    double tb2 = sing_box_[kj].box_limit_[1].first;
	    double tc2 = sing_box_[kj].box_limit_[2].first;
	    double td2 = sing_box_[kj].box_limit_[3].first;
	    double u1 = sing_box_[ki].sing_->getPar(0);
	    double v1 = sing_box_[ki].sing_->getPar(1);
	    double u2 = sing_box_[kj].sing_->getPar(0);
	    double v2 = sing_box_[kj].sing_->getPar(1);
	    double u_val = fabs(u2 - u1);
	    double v_val = fabs(v2 - v1);
	    //double u_val = (ta1 > ta2) ? tb2 - ta1 : tb1 - ta2;
	    //double v_val = (tc1 > tc2) ? td2 - tc1 : td1 - tc2;
	    
	    if (fabs(u2-u1) < ptol1 && fabs(v2-v1) < ptol2)
		continue;  // Too close singularities
	    if (u_val > v_val || fabs(v2-v1) < ptol2)
	    {
		// Make division value in u-direction
		double umin = std::min(u1, u2);
		double umax = std::max(u1, u2);
		double ah = std::max(ta1, ta2);
		double bh = std::min(tb1, tb2);
		double uval = (ah > umin && bh > umin && ah < umax && bh < umax)
		    ? 0.5*(ah+bh) : 0.5*(u1+u2);

		// Get all overlapping singularity boxes
		vector<pair<double, pair<double, double> > > overlap_u = getOverlapBox(umin, umax, 0);

		// Find best position of division value
		double u_div = getSplitValue(overlap_u, umin, umax, uval);
		divpar_u.push_back(make_pair(u_div, 2));
	    }
	    else
	    {
		// Make division value in u-direction
		double vmin = std::min(v1, v2);
		double vmax = std::max(v1, v2);
		double ah = std::max(tc1, tc2);
		double bh = std::min(td1, td2);
		double vval = (ah > vmin && bh > vmin && ah < vmax && bh < vmax)
		    ? 0.5*(ah+bh) : 0.5*(v1+v2);

		// Get all overlapping singularity boxes
		vector<pair<double, pair<double, double> > >  overlap_v = getOverlapBox(vmin, vmax, 1);

		// Find best position of division value
		double v_div = getSplitValue(overlap_v, vmin, vmax, vval);
		divpar_v.push_back(make_pair(v_div, 2));
	    }
	}

    // Add division values outside the singularity boxes if required
    double umin = surf_->endParam(0);
    double umax = surf_->startParam(0);
    double vmin = surf_->endParam(1);
    double vmax = surf_->startParam(1);
    double udel = 0.25*(umin - umax)/(double)(divpar_u.size())+1;
    double vdel = 0.25*(vmin - vmax)/(double)(divpar_v.size())+1;
    for (ki=0; ki<sing_box_.size(); ki++)
    {
	umin = std::min(umin,sing_box_[ki].box_limit_[0].first);
	umax = std::max(umax,sing_box_[ki].box_limit_[1].first);
	vmin = std::min(vmin,sing_box_[ki].box_limit_[2].first);
	vmax = std::max(vmax,sing_box_[ki].box_limit_[3].first);
    }
    if (umin > surf_->startParam(0) + udel)
	divpar_u.push_back(make_pair(umin, 2));
    else
	divpar_u.push_back(make_pair(surf_->startParam(0), 2));
    if (umax < surf_->endParam(0) - udel)
	divpar_u.push_back(make_pair(umax, 2));
    else
	divpar_u.push_back(make_pair(surf_->endParam(0), 2));
    if (vmin > surf_->startParam(1) + vdel)
	divpar_v.push_back(make_pair(vmin, 2));
    else
	divpar_v.push_back(make_pair(surf_->startParam(1), 2));
    if (vmax < surf_->endParam(1) - vdel)
	divpar_v.push_back(make_pair(vmax, 2));
    else
	divpar_v.push_back(make_pair(surf_->endParam(1), 2));

    // Sort
    std::sort(divpar_u.begin(), divpar_u.end(), compare_divpar);
    std::sort(divpar_v.begin(), divpar_v.end(), compare_divpar);

    if (getenv("DEBUG_DIV") && (*getenv("DEBUG_DIV"))=='1')
    {
	ofstream debug1("init_div.dsp");
	debug1 << "fg: red" << endl;
	for (ki=0; ki<divpar_u.size(); ki++)
	{
	    debug1 << "lin:" << endl;
	    debug1 << divpar_u[ki].first << "  " << mima[2] << endl;
	    debug1 << divpar_u[ki].first << "  " << mima[3] << endl;
	}
	for (ki=0; ki<divpar_v.size(); ki++)
	{
	    debug1 << "lin:" << endl;
	    debug1 << mima[0] << "  " << divpar_v[ki].first << endl;
	    debug1 << mima[1] << "  " << divpar_v[ki].first << endl;
	}
    }

}

//===========================================================================
bool SfSelfIntersector::isSeparated(shared_ptr<IntersectionPoint> sing1,
				    shared_ptr<IntersectionPoint> sing2,
				    vector<pair<double,int> >& divpar_u,
				    vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    // u-direction
    for (size_t ki=0; ki<divpar_u.size(); ki++)
    {
	double s1 = sing1->getPar(0);
	double s2 = sing2->getPar(0);
	double dp = divpar_u[ki].first;
	if ((s1-dp)*(s2-dp) < 0.0)
	    return true;   // Division value lies between singularities
    }

    // v-direction
    for (size_t ki=0; ki<divpar_v.size(); ki++)
    {
	double s1 = sing1->getPar(1);
	double s2 = sing2->getPar(1);
	double dp = divpar_v[ki].first;
	if ((s1-dp)*(s2-dp) < 0.0)
	    return true;   // Division value lies between singularities
    }

    return false;  // Not separated
}

//===========================================================================
vector<pair<double, pair<double, double> > >
SfSelfIntersector::getOverlapBox(double minval, double maxval, int dir)
//===========================================================================
{
    vector<pair<double, pair<double, double> > > overlaps;
    for (size_t ki=0; ki<sing_box_.size(); ki++)
    {
	double ta = sing_box_[ki].box_limit_[2*dir].first;
	double tb = sing_box_[ki].box_limit_[2*dir+1].first;
	if (!((ta < minval && tb < minval) ||
	      (ta > maxval && tb > maxval)))
	    overlaps.push_back(make_pair(sing_box_[ki].sing_->getPar(dir),
					 make_pair(ta, tb)));
    }
    return overlaps;
}

bool compare_overlap(const pair<double, pair<double,double> >& l1, 
		     const pair<double, pair<double,double> >& l2)
{ 
    return (l1.first < l2.first);
}

//===========================================================================
double 
SfSelfIntersector::getSplitValue(vector<pair<double, pair<double, double> > >& overlap,
				 double minval, double maxval, double divval)
//===========================================================================
{
    // VSK, 0207. Currently a very simple version. Must be refined

    double frac = 0.1;
    if (divval < minval + frac*(maxval-minval))
	divval = minval + frac*(maxval-minval);
    if (divval > maxval - frac*(maxval-minval))
	divval = maxval - frac*(maxval-minval);

    // Sort overlap boxes
    std::sort(overlap.begin(), overlap.end(), compare_overlap);

    // Find singularities on each side of the candidate division value
    size_t size = overlap.size();
    double val;
    if (overlap[0].first > divval)
    {
	val = overlap[0].second.first;
    }
    else if (overlap[size-1].first <= divval)
    {
	val = overlap[size-1].second.second;
    }
    else
    {
	size_t ki;
	for (ki=1; ki<size; ki++)
	    if (overlap[ki-1].first <= divval && overlap[ki].first > divval)
		break;
	//val = 0.5*(overlap[ki-1].second.second + overlap[ki].second.first);
	val = 0.5*(overlap[ki-1].first + overlap[ki].first);
    }
    if (val < minval || val > maxval)
	val = 0.5*(minval + maxval);
    return val;
}

//===========================================================================
void SfSelfIntersector::
refineSingularityLinks(vector<shared_ptr<IntersectionPoint> >& sing,
		       vector<pair<double,int> >& divpar,
		       int dir)
//===========================================================================
{
    // Refine subdivision curves where it cross division lines
    // First fetch all links
    vector<shared_ptr<IntersectionLink> > links, tmp_links;
    for (size_t ki=0; ki<sing.size(); ki++)
    {
	sing[ki]->getNeighbourLinks(tmp_links);
	
	for (size_t kr=0; kr<tmp_links.size(); kr++)
	{
	    // Check if the link exist already
	    size_t kj;
	    for (kj=0; kj<links.size(); kj++)
		if (links[kj] == tmp_links[kr])
		    break;
	    if (kj == links.size())
		links.push_back(tmp_links[kr]);
	}
    }

    if (links.size() == 0)
	return;   // No links to refine


    // Get normal surface
     shared_ptr<ParamSurfaceInt> normsf = surf_->getNormalSurface();
    if (normsf.get() == 0)
	return;  // @@@ Not a spline surface

    // Define origo
    Point origo(0.0, 0.0, 0.0);
   
    IntersectionPoint *p1, *p2;
    double pmin, pmax, pmin2, pmax2;
    double param[2];
    for (size_t ki=0; ki<links.size(); ki++)
    {
	// Fetch singular intersection points
	links[ki]->getIntersectionPoints(p1, p2);
	if (p2->getPar(dir) < p1->getPar(dir))
	    std::swap(p1, p2);

	pmin = p1->getPar(dir);
	pmax = p2->getPar(dir);
	pmin2 = p1->getPar(1-dir);
	pmax2 = p2->getPar(1-dir);


	// Find crosses in the u_direction
	size_t idx1, idx2;
	for (idx1=0; idx1<divpar.size() && divpar[idx1].first<=pmin; 
	     idx1++);
	for (idx2=idx1; idx2<divpar.size() && divpar[idx2].first<pmax; 
	     idx2++);

	for (size_t kj=idx1; kj<idx2; kj++)
	{
	    // Iterate for a singular point. First fetch constant
	    // parameter curve
	    param[dir] = divpar[kj].first;
	    shared_ptr<ParamCurve> normcrv = 
		normsf->getConstantParameterCurve(dir, param[dir]);

	    // Iterate
	    double seed = ((pmax - param[dir])*pmin2 
			   + (param[dir] - pmin)*pmax2)/(pmax - pmin);
	    double tstart = std::min(pmin2, pmax2) - 0.5*fabs(pmax2 - pmin2);
	    tstart = std::max(normcrv->startparam(), tstart);
	    double tend = std::max(pmin2, pmax2) + 0.5*fabs(pmax2 - pmin2);
	    tend = std::min(normcrv->endparam(), tend);

	    double s_dist;
	    Point s_pt;
	    normcrv->closestPoint(origo, tstart, tend, param[1-dir], s_pt, 
				  s_dist,&seed);

	    if (s_dist > 2.0*epsge_->getEpsge())
		continue;  // Don't know what to do in this case

	    // Insert singular point
	    shared_ptr<IntersectionPoint> tmp = 
		int_results_->addIntersectionPoint(surf_, surf_, 
						   epsge_, param, param);
	    
	    p1->connectTo(tmp, LINK_UNDEFINED);
	    p2->connectTo(tmp, LINK_UNDEFINED);
	    p1->disconnectFrom(p2);

	    sing.push_back(tmp);
	    p1 = tmp.get();
	    pmin = p1->getPar(dir);
	    pmin2 = p1->getPar(1-dir);
	}
	    
    }
}


//===========================================================================
void SfSelfIntersector::
makeDivisionLines(vector<shared_ptr<IntersectionPoint> >& sing,
		  vector<double>& sing_par, int dir,
		  double start, double end,
		  vector<pair<double,int> >& divpar)
//===========================================================================
{

    // Add extra grid lines external to the singularity domain and
    // between grid lines of large distance
    if (false /*sing_par.size() == 1*/)
	divpar.push_back(make_pair(sing_par[0], 1));
    else
    {
	double tb1, tb2, td;
	double tc = 0.01*(end-start);
	double tdel = (sing_par.size() == 1) ? tc : sing_par[sing_par.size()-1] - sing_par[0];
	double td1 = (sing_par.size() == 1) ? tdel : tdel/(double)(2*(sing_par.size()-1));
	double ptol = std::max(std::min(0.00001, 0.01*(end-start)),
			       0.01*tdel);
	if (td1 > tc)
	    td1 = tc;

	// Get knot values
	vector<double> knots = surf_->getInnerKnotVals(dir, false); // Not
	                                                            // sorted
	size_t knot_idx = 0;
	double tknot, tdiv;
	bool found;
	double frac = 0.2;

	if (sing_par[0] > start+ptol)
	{
	    tb1 = (sing_par[0] - td1 < start+0.5*td1)
		? start : sing_par[0] - td1;
	    tdiv = frac*(sing_par[0] - start);
	    found = closeKnot(tb1, knots, knot_idx, tdiv, tknot);
	    if (found && tknot < sing_par[0]-ptol)
		tb1 = tknot;
	    divpar.push_back(make_pair(tb1, 2));
	}
	else
	    divpar.push_back(make_pair(start, 1));

	for (size_t ki=1; ki<sing_par.size(); ki++)
	{
	    td = sing_par[ki]-sing_par[ki-1];
	    tdiv = frac*td;
	    if (td <8*td1)
	    {
		tb1 = sing_par[ki-1]+0.5*td;
		found = closeKnot(tb1, knots, knot_idx, tdiv, tknot);
		if (found)
		    tb1 = tknot;
		divpar.push_back(make_pair(tb1, 2));
	    }
	    else
	    {
		tb1 = sing_par[ki-1]+td1;
		found = closeKnot(tb1, knots, knot_idx, tdiv, tknot);
		if (found)
		    tb1 = tknot;
		tb2 = sing_par[ki]-td1;
		found = closeKnot(tb2, knots, knot_idx, tdiv, tknot);
		if (found)
		    tb2 = tknot;
		if (tb1 > divpar[divpar.size()-1].first + ptol)
		    divpar.push_back(make_pair(tb1, 2));
		divpar.push_back(make_pair(tb2, 2));
	    }
	}

	if (sing_par[sing_par.size()-1] < end-ptol)
	{
	    tb2 = (sing_par[sing_par.size()-1] + td1 > end-0.5*td1) ? 
		end : sing_par[sing_par.size()-1] + td1;
	    tdiv = frac*(end - sing_par[sing_par.size()-1]);
	    found = closeKnot(tb2, knots, knot_idx, tdiv, tknot);
	    if (found && tknot > sing_par[sing_par.size()-1]+ptol)
		tb2 = tknot;
	    divpar.push_back(make_pair(tb2, 2));
	}
	else
	    divpar.push_back(make_pair(end, 1));
    }
}  


//===========================================================================
void 
SfSelfIntersector::modifyDivisionValues(vector<shared_ptr<IntersectionPoint> >& sing,
					vector<pair<double,int> >& divpar, int dir)
{
    // Sort singularities in the current parameter direction. Include end points
    // of surface
    size_t ki, kj;
    vector<double> sing_par(sing.size()+2);
    sing_par[0] = surf_->startParam(dir);
    for (ki=0; ki<sing.size(); ki++)
    {
	sing_par[ki+1] = sing[ki]->getPar(dir);
    }
    sing_par[ki+1] = surf_->endParam(dir);
    std::sort(sing_par.begin(), sing_par.end());

    // Fetch knots
    vector<double> knots = surf_->getInnerKnotVals(dir, false); // Not sorted
    double td, td1, td2;
    double fac = 0.2;
    double curr_knot;
    size_t idx = 0;

    // For each division value, set the range of possible modification and look
    // for knots in this range
    for (ki=0, kj=0; ki<divpar.size(); ki++)
    {
	for (; kj<sing_par.size()-1; kj++)
	{
	    if (sing_par[kj+1] > divpar[ki].first)
		break;
	}
	if (ki == 0)
	    td1 = divpar[ki].first - sing_par[kj];
	else
	    td1 = std::min(divpar[ki].first - sing_par[kj], 
			   divpar[ki].first - divpar[ki-1].first);

	if (ki == divpar.size()-1)
	    td2 = sing_par[kj+1] - divpar[ki].first;
	else
	    td2 = std::min(sing_par[kj+1] - divpar[ki].first, 
			   divpar[ki+1].first - divpar[ki].first);
	td = fac*(std::min(td1, td2));

	if (closeKnot(divpar[ki].first, knots, idx, td, curr_knot))
	    divpar[ki].first = curr_knot;
    }
}

			  
//===========================================================================
bool SfSelfIntersector::
closeKnot(double par, vector<double>& knots, size_t& knot_idx,
	  double par_div, double& curr_knot)
//===========================================================================
{
    size_t kj, kr, kmin;
    double tmin;

    curr_knot = par;

    // Find closest knot within reasonable distance
    for (; knot_idx<knots.size(); knot_idx++)
	if (knots[knot_idx] >= par - par_div)
	    break;

    for (kr=knot_idx; kr<knots.size(); kr++)
	if (knots[kr] > par + par_div)
	    break;

    if (knot_idx >= knots.size() || kr == knot_idx)
	return false;  // No knots in specified interval

    tmin = fabs(knots[knot_idx]-par);
    kmin = knot_idx;
    for (kj=knot_idx+1; kj<kr; kj++)
    {
	double tdum = fabs(knots[kj] - par);
	if (tdum < tmin)
	{
	    tmin = tdum;
	    kmin = kj;
	}
    }

    curr_knot = knots[kmin];
    return true;
}

//===========================================================================
void SfSelfIntersector::
checkDivisionLines(vector<shared_ptr<IntersectionPoint> >& sing,
		   vector<pair<double,int> >& divpar_u,
		   vector<pair<double,int> >& divpar_v)
//===========================================================================
{
    if (divpar_u.size() < 2 || divpar_v.size() < 2)
	return;   // Not much to do

    double tdel1 = divpar_u[divpar_u.size()-1].first - divpar_u[0].first;
    double frac1 = tdel1/(double)(divpar_u.size()-1);
    double tdel2 = divpar_v[divpar_v.size()-1].first - divpar_v[0].first;
    double frac2 = tdel2/(double)(divpar_v.size()-1);

    if (tdel1 < tdel2)
    {
	checkOneDivLine(sing, 0, frac1, divpar_u, divpar_v);
	checkOneDivLine(sing, 1, frac2, divpar_v, divpar_u);
    }
    else
    {
	checkOneDivLine(sing, 1, frac2, divpar_v, divpar_u);
	checkOneDivLine(sing, 0, frac1, divpar_u, divpar_v);
    }


}


//===========================================================================
void SfSelfIntersector::
checkOneDivLine(vector<shared_ptr<IntersectionPoint> >& sing,
		int dir, double frac, vector<pair<double,int> >& div1,
		vector<pair<double,int> >& div2)
//===========================================================================
{
     size_t ki, kj, kr, kh;
     for (ki=1; ki<div1.size()-1; ki++)
     {
	 if (div1[ki].second != 2)
	     continue;    // Not a candidate for removal
	  
	 // Collect singular points on the cells next to the current
	 // division line
	 double tmin = div1[ki-1].first;
	 double tmed = div1[ki].first;
	 double tmax = div1[ki+1].first;
	 if (tmed - tmin > frac && tmax - tmed > frac)
	     continue;   // Not a narrow division interval

	 vector<shared_ptr<IntersectionPoint> > sing1, sing2;
	 for (kj=0; kj<sing.size(); kj++)
	 {
	     double par = sing[kj]->getPar(dir);
	     if (par > tmin && par <= tmed)
		 sing1.push_back(sing[kj]);
	     if (par >= tmed && par < tmax)
		 sing2.push_back(sing[kj]);
	 }

	 double t1, t2;
	 for (kj=0; kj<sing1.size(); kj++)
	 {
	     t1 = sing1[kj]->getPar(1-dir);
	     for (kr=0; kr<sing2.size(); kr++)
	     {
		 t2 = sing2[kr]->getPar(1-dir);
		 if (t2 < t1)
		     std::swap(t1, t2);

		 // Traverse the v-division lines to find a line
		 // dividing the two singularities
		 for (kh=0; kh<div2.size(); kh++)
		 {
		     if (t1 < div2[kh].first && t2 > div2[kh].first)
			 break;
		 }

		 if (kh == div2.size())
		     break;  // No division line
	     }

	     if (kr < sing2.size())
		 break;  // No division line
	 }
	 if (kj == sing1.size())
	 {
	     // Division lines found for all instances
	     // Remove current division line
	     div1.erase(div1.begin()+ki, div1.begin()+ki+1);
	     ki--;
	 }
     }
}


//===========================================================================
void SfSelfIntersector::
validateSingularities(vector<shared_ptr<IntersectionPoint> >& sing)
//===========================================================================
{
    // Purpose: Dismiss singularities that do not belong to a self
    // intersection curve or are stand alone

    // Parameter domain
    vector<double> mima = surf_->getMima();

    // Parameter distances to check if a singularity is stand alone
    double frac = 0.1;
    double del1 = frac*(mima[1] - mima[0]);
    double del2 = frac*(mima[3] - mima[2]);

    // For each singularity
    for (size_t ki=0; ki<sing.size(); ki++)
    {
	// Check if it is connected in a list
	if (sing[ki]->numNeighbours() > 0)
	    continue;    // This point should not be removed

	// Check if there is anothter singularity within stand alone
	// distance
	size_t kj;
	double dmin1 = std::max(mima[0], sing[ki]->getPar(0)-del1);
	double dmax1 = std::min(mima[1], sing[ki]->getPar(0)+del1);
	double dmin2 = std::max(mima[2], sing[ki]->getPar(1)-del2);
	double dmax2 = std::min(mima[3], sing[ki]->getPar(1)+del2);
	bool found = false;
	for (kj=0; kj<sing.size(); kj++)
	{
	    if (ki == kj)
		continue;

	    if (sing[kj]->getPar(0) < dmin1 || sing[kj]->getPar(0) > dmax1 ||
		sing[kj]->getPar(1) < dmin2 || sing[kj]->getPar(1) > dmax2)
		continue;  // Not withing stand alone distance or
			   // closest distance

	    found = true;
	    // Adapt closest distance
	    if (sing[kj]->getPar(0) < sing[ki]->getPar(0))
		dmin1 = sing[kj]->getPar(0);
	    if (sing[kj]->getPar(0) > sing[ki]->getPar(0))
		dmax1 = sing[kj]->getPar(0);
	    if (sing[kj]->getPar(1) < sing[ki]->getPar(1))
		dmin2 = sing[kj]->getPar(1);
	    if (sing[kj]->getPar(1) > sing[ki]->getPar(1))
		dmax2 = sing[kj]->getPar(1);
	}

	if (!found)
	    continue;  // Stand alone singularity. Keep it

	// Fetch a surface piece around the singularity half distance
	// to the next singularity
	dmin1 = 0.5*(dmin1 + sing[ki]->getPar(0));
	dmax1 = 0.5*(dmax1 + sing[ki]->getPar(0));
	dmin2 = 0.5*(dmin2 + sing[ki]->getPar(1));
	dmax2 = 0.5*(dmax2 + sing[ki]->getPar(1));
	vector<shared_ptr<ParamSurfaceInt> > sub_sfs =
	    surf_->subSurfaces(dmin1, dmin2, dmax1, dmax2, 
			       epsge_->getRelParRes());

	for (kj=0; kj<sub_sfs.size(); kj++)
	{
	    if (hasBoundaryIntersections(sub_sfs[kj]))
		break;
	}

	if (kj == sub_sfs.size())
	{
	    // No selfintersections along the sub surface boundary is
	    // found.  Dismiss the singularity
	    int_results_->removeIntPoint(sing[ki]);
	    sing.erase(sing.begin()+ki, sing.begin()+ki+1);
	}
    }
	    
}


//===========================================================================
bool SfSelfIntersector::
hasBoundaryIntersections(shared_ptr<ParamSurfaceInt> sub_sf)
//===========================================================================
{
    // Purpose: Check if a surface intersect itself at the boundaries
    static double fac = 0.1;
    double aeps = fac*epsge_->getEpsge();

    // Fetch boundary curves
    std::vector<shared_ptr<BoundaryGeomInt> > bd_obj;
    sub_sf->getBoundaryObjects(bd_obj);
    size_t kh;
    for (kh=0; kh<bd_obj.size(); kh++) {
	shared_ptr<ParamObjectInt> bd_obj2 = bd_obj[kh]->getObject();
	shared_ptr<ParamCurveInt> bd_cv = 
	    dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(bd_obj2);  


	// Intersect one boundary curve with the surface
	SfCvIntersector sfcvintersect (sub_sf, bd_cv, aeps);
	sfcvintersect.setSelfintCase(1);
	sfcvintersect.compute();

	// Remove the trivial intersection
	sfcvintersect.getIntPool()->removeBoundaryIntersections();

	vector<shared_ptr<IntersectionPoint> > intpts;
	vector<shared_ptr<IntersectionCurve> > intcrv;

	sfcvintersect.getResult(intpts, intcrv);

	if (intpts.size() > 0 || intcrv.size() > 0)
	    return true;
    }
    return false;
}


//===========================================================================
void SfSelfIntersector::
computeBoundaryIntersections(shared_ptr<ParamSurfaceInt> sub_sf,
			     IntersectionPoint *sing,
			     int is_handled[])
//===========================================================================
{
    // Purpose: Check if a surface intersect itself at the boundaries

    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
    {
	std::ofstream debug("self_sub_out.g2");
	shared_ptr<ParamSurface> srf = sub_sf->getParamSurface();
	srf->writeStandardHeader(debug);
	srf->write(debug);
    }

    static double fac = 0.1;
    double aeps = fac*epsge_->getEpsge();
    double ptol = epsge_->getRelParRes();
    vector<double> mima = sub_sf->getMima();  // Parameter domain of surface

    // Make surface-surface intersector for storage of results on the
    // correct format
    SfSfIntersector sfsfintersect(sub_sf, sub_sf, aeps, this);

    // Fetch boundary curves
    std::vector<shared_ptr<BoundaryGeomInt> > bd_obj;
    sub_sf->getBoundaryObjects(bd_obj);
    size_t kh;
    int nmb1, nmb2=0, nmb3; // Debug
    size_t nmb_prev = 0;
    nmb3 = getIntPool()->numIntersectionPoints();
    for (kh=0; kh<bd_obj.size(); kh++) 
    {
	int bd_idx;
	if (bd_obj[kh]->getDir() == 0)
	    bd_idx = (fabs(bd_obj[kh]->getPar()-mima[0])<ptol) ? 0 : 1;
	else
	    bd_idx = (fabs(bd_obj[kh]->getPar()-mima[2])<ptol) ? 2 : 3;

	shared_ptr<ParamObjectInt> bd_obj2 = bd_obj[kh]->getObject();
	shared_ptr<ParamCurveInt> bd_cv = 
	    dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(bd_obj2);  


	// Intersect one boundary curve with the surface
	SfCvIntersector sfcvintersect (sub_sf, bd_cv, aeps, &sfsfintersect,
				       bd_obj[kh]->getDir()+2, 
				       bd_obj[kh]->getPar());
	sfcvintersect.setSelfintCase(1);

	// Existing intersection points may be singular
	vector<shared_ptr<IntersectionPoint> > sfcv_ints;
	sfcvintersect.getIntPool()->getIntersectionPoints(sfcv_ints);
	for (size_t kr=0; kr<sfcv_ints.size(); kr++)
	{
	    if (sfcv_ints[kr]->hasParentPoint() &&
		sfcv_ints[kr]->parentPoint()->getSingularityType() == HIGHER_ORDER_POINT)
	    {
		double param[3];
		for (int kh=0; kh<3; kh++)
		    param[kh] = sfcv_ints[kr]->getPar(kh);
		sfcvintersect.setHighPriSing(param);
		sfcvintersect.getSingularityInfo()->setHighPriSingType(KEEP_SING);
	    }
	}
	
	//if (is_handled[bd_idx] == 0)
	    sfcvintersect.compute();

	    sfcvintersect.postIterateBd();

	// Remove the trivial intersection
	sfcvintersect.getIntPool()->cleanUpPool(0, aeps);
	sfcvintersect.getIntPool()->removeBoundaryIntersections(false);
	sfsfintersect.getIntPool()->includeReducedInts
	(sfcvintersect.getIntPool());

	sfsfintersect.postIterate3(nmb2, bd_obj[kh]->getDir()+2);
	nmb1 = sfcvintersect.getIntPool()->numIntersectionPoints();
	nmb2 = sfsfintersect.getIntPool()->numIntersectionPoints();
    }

    /*for (kh=0; kh<4; kh++)
    {
	sfsfintersect.postIterate3(0, kh);
	}*/

    // Now, other sub surface  is intersected with the the boundaries 
    // of one sub surface. To get all instances of the self intersection
    // points, the result must be mirrored. This is the same as intersecting
    // the boundaries of the sub surface with the sub surface.
    //sfsfintersect.getIntPool()->mirrorIntersectionPoints(nmb_prev);

    // Connect boundary intersections
    if (sing)
	sfsfintersect.selfintComplex(sing);
    else
	sfsfintersect.selfintComplex2();

    int_results_->selfIntersectParamReorganise(sfsfintersect.getIntPool());
    sfsfintersect.getIntPool()->removeDoublePoints();

    // Fetch intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    sfsfintersect.getIntPool()->getIntersectionPoints(int_pts);

    // Reset all twint point information
    shared_ptr<IntersectionPoint> dummy;
    for (kh=0; kh<int_pts.size(); kh++)
	int_pts[kh]->setParentPoint(dummy);

    nmb_prev = sfsfintersect.getIntPool()->numIntersectionPoints();

}


//===========================================================================
void SfSelfIntersector::getNmbDivision(int& div1, int& div2)
//===========================================================================
{
    // Purpose: Decide the number of subdivisions in each parameter
    // direction

    double length[2], wiggle[2];  // Estimated length and wiggliness
    bool inner_knots[2], critical_val[2], can_divide[2];
    int nmbdir = 2;   // Number of parameter directions

    surf_->getLengthAndWiggle(length, wiggle);
    for (int ki=0; ki<nmbdir; ki++) {
	inner_knots[ki] = surf_->hasInnerKnots(ki);

	critical_val[ki] = surf_->hasCriticalVals(ki);

	can_divide[ki] = surf_->canDivide(ki);
    }

    // Initial choice 
    div1 = div2 = 5;

    // Modify according to surface size and properties
    double fac1 = 0.2;
    double fac2 = 0.05;
    if (length[0] < fac1*length[1])
	div1 = 3;
    if (length[0] < fac2*length[1])
	div1 = 0;

    if (length[1] < fac1*length[0])
	div2 = 3;
    if (length[1] < fac2*length[0])
	div2 = 0;

    if (wiggle[0] > fac1*wiggle[1] || critical_val[0] || inner_knots[0])
	div1 = std::max(1, div1);

    if (wiggle[1] > fac1*wiggle[0] || critical_val[1] || inner_knots[1])
	div2 = std::max(1, div2);

    if (!can_divide[0])
	div1 = 0;

    if (!can_divide[1])
	div2 = 0;

}


//===========================================================================
void 
SfSelfIntersector::setDivisionValues(vector<pair<double,int> >& divpar,
				     int& nmbpar, int dir)
//===========================================================================
{
    // Purpose: Decide the subdivisions parameters in the given
    // parameter direction

    int ki;
    double ta = surf_->startParam(dir);
    double tb = surf_->endParam(dir);

    // If the endparameters are included in the division values,
    // remove these
    if (divpar.size()>0 && divpar[0].first < ta+epsge_->getRelParRes())
	divpar.erase(divpar.begin(), divpar.begin()+1);
    if (divpar.size()>0 && 
	divpar[divpar.size()-1].first > tb-epsge_->getRelParRes())
	divpar.erase(divpar.end()-1, divpar.end());
	
    int knsing = (int)(divpar.size());
    if (knsing >= nmbpar)
	return;   // Enough split parameters

    // Make sure to include the singularities in the splitting
    // parameters
    nmbpar = std::max(nmbpar, (int)(divpar.size()));  
    // double frac = (tb - ta)/(double)(5*(nmbpar+1));
    //double frac2 = 0.1*frac;

    // Distribute the remaining division parameters. 
    // Make initial distribution
    vector<double> initpar(nmbpar - knsing);
    double tint = (tb - ta)/(double)(initpar.size()+1);
    double par;
    for (par=ta+tint, ki=0; ki<(int)(initpar.size()); par+=tint, ki++)
	initpar[ki] = par;

    // Get knot values
    vector<double> knots = surf_->getInnerKnotVals(dir, false); // Not
								// sorted

    int kj, kr, kh, km, kk, kn;
    int kmin;
    double tmin;
    double tlast = ta;
    if (knsing == 0) {
	// No singular values to consider
	for (ki=0, kj=0, kr=0; ki<(int)(initpar.size()); ki++) {
	    // Find closest knot within reasonable distance
	    for (kj=kr; kj<(int)(knots.size()); kj++)
		if (knots[kj] >= initpar[ki] - 0.4*tint)
		    break;

	    for (kr=kj; kr<(int)(knots.size()); kr++)
		if (knots[kr] > initpar[ki] + 0.4*tint)
		    break;

	    tmin = (kj<(int)(knots.size())) ? fabs(knots[kj]-initpar[ki]) : 
		1.0e6;
	    kmin = (kr > kj) ? kj : -1;
	    for (kh=kj+1; kh<kr; kh++) {
		double tdum = fabs(knots[kh]-initpar[ki]);
		if (tdum < tmin) {
		    tmin = tdum;
		    kmin = kh;
		}
	    }
	    if (kmin >=0 && kmin < int(knots.size()) && 
		knots[kmin] > tlast + 0.1*tint)
		divpar.push_back(make_pair(knots[kmin], 0));
	    else
		divpar.push_back(make_pair(initpar[ki], 0));
	    tlast = divpar[divpar.size()-1].first;
	}
    } else {
	// For each interval between singular values, distribute a
	// suitable number of subdivision parameter values.  First
	// redistribute values in initpar
	vector<pair<double,int> > tmpdiv;
	tmpdiv.reserve(nmbpar);

	// of division values
	double ta1 = ta;
	double tb1 = divpar[0].first; 
	for (ki=0, kj=0, km=0, kn=0; ki<nmbpar; ki++)
	{
	    // First check how many extra division values to place
	    // between the already defined ones
	    for (kk=kj; kk<(int)(initpar.size()); kk++)
		if (initpar[kk] > tb1)
		    break;
	    double tint2 = (tb1 - ta1)/(double)(kk-kj+1);  // Initial parameter
	    // distance between division values

	    if (tint2 > 0.5*tint)
	    {
		// Redistribute initial division values
		for (par=ta1+tint2, kh=kj; kh<kk; kh++, par+=tint2)
		    initpar[kh] = par;

		for (; kj<(int)(initpar.size()) && initpar[kj]<tb1; kj++)
		{
		    // Adjust division parameter according to knot
		    // First find correct knot
		    for (; kn<(int)(knots.size()); kn++)
			if (knots[kn]>= initpar[kj] - 0.4*tint2)
			    break;
		    for (kh = kn; kh < (int)knots.size()
			     && knots[kh]<initpar[kj]+0.4*tint2; kh++);
		    if (kh-kn > 0)
		    {
			// Division value must be adjusted. Find closest knot
			tmin = fabs(initpar[kj]-knots[kn]);
			kmin = kn;
			for (kk=kn+1; kk<kh; kk++)
			    if (fabs(initpar[kj]-knots[kk]) < tmin)
			    {
				tmin = fabs(initpar[kj]-knots[kk]);
				kmin = kk;
			    }
			initpar[kj] = knots[kmin];
		    }

		    tmpdiv.push_back(make_pair(initpar[kj], 0));
		}
	    }

	    if (tb1 < tb-epsge_->getRelParRes())
	    {
		// Add the singular division value
		tmpdiv.push_back(make_pair(tb1, divpar[km].second));
	    }
	    ta1 = tb1;
	    tb1 = (km >= (int)(divpar.size()-1)) ? tb : divpar[++km].first;
	}
	divpar = tmpdiv;  // Copy to output vector

	/*for (km=0, ki=0, kj=0, kr=0; km<=(int)(divpar.size()); 
	     km++, ki=kk) {
	    if (tb1 < initpar[ki]) {
		kk = ki;
		tlast = tb1;
		km--;
	    } else {

		for (kk=ki; kk<(int)(initpar.size()); kk++)
		    if (initpar[kk] > tb1)
			break;

		double tint2 = (tb1 - ta1)/(double)(kk-ki+1);
		for (kn=ki, par=ta1+tint2; kn<kk; kn++, par+=tint2)
		    initpar[kn] = par;

		for (kn=ki; kn<kk; kn++) {
		    // Adjust division parameters according to knots
		    for (kj=kr; kj<(int)(knots.size()); kj++)
			if (knots[kj] >= initpar[kn] - 0.4*tint2)
			    break;

		    for (kr=kj; kr<(int)(knots.size()); kr++)
			if (knots[kr] > initpar[kn] + 0.4*tint2)
			    break;

		    tmin = (kj<(int)(knots.size())) ? 
			fabs(knots[kj]-initpar[kn]) : 1.0e6;
		    kmin = (kr > kj) ? kj : -1;
		    for (kh=kj+1; kh<kr; kh++) {
			double tdum = fabs(knots[kh]-initpar[kn]);
			if (tdum < tmin) {
			    tmin = tdum;
			    kmin = kh;
			}
		    }
		    if (kmin >= 0 
			&& kmin < (int)(knots.size()) 
			&& knots[kmin] > tlast + 0.1*tint) {
			divpar.insert(divpar.begin()+km, 
				      make_pair(knots[kmin], 0));
		    } else {
			divpar.insert(divpar.begin()+km,
				      make_pair(initpar[kn], 0));
		    }
		    tlast = divpar[km].first;
		    km++;
		}
	    }

	    ta1 = tb1;
// 	    tb1 = (km == (int)(divpar.size())-1) ? tb : divpar[km+1].first;
	    tb1 = (km == (int)(divpar.size())-1) ? tb : divpar[km].first;
	} */
    }

}


//===========================================================================
bool
SfSelfIntersector::hasSingularity(shared_ptr<ParamSurfaceInt> curr_sub)
//===========================================================================
{
    if (!hasSingularityInfo())
	return false;  // No singularities are found

    // Check each singularity for whether it lies in the current
    // surfaces
    double ptol = epsge_->getRelParRes();
    int nmb_sing = singularity_info_->getNmbPoints(2);
    for (int ki=0; ki<nmb_sing; ki++)
    {
	Point curr = singularity_info_->getPoint(ki, 2);
	int kj;
	for (kj=0; kj<2; kj++)
	{
	    if (curr[kj] < curr_sub->startParam(kj) - ptol ||
		curr[kj] > curr_sub->endParam(kj) + ptol)
		break;
	}
	if (kj < 2)
	    continue;  // The singularity does not lie in the surface

	// A singularity is found
	return true;
    }

    // Check higher order singularity points belonging to this sub problem
     vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    for (size_t ki=0; ki<int_pts.size(); ki++)
	if (int_pts[ki]->getSingularityType() == HIGHER_ORDER_POINT)
	    return true;
   
    return false;
}


//===========================================================================
bool
SfSelfIntersector::shareSingularity(shared_ptr<ParamSurfaceInt> sub1, 
				    shared_ptr<ParamSurfaceInt> sub2, 
				    double sing[])
//===========================================================================
{
    if (!hasSingularityInfo())
	return false;  // No singularities are found

    // Check each singularity for whether it lies in the current
    // surfaces
    double ptol = epsge_->getRelParRes();
    int nmb_sing = singularity_info_->getNmbPoints(2);
    for (int ki=0; ki<nmb_sing; ki++)
    {
	Point curr = singularity_info_->getPoint(ki, 2);
	int kj;
	for (kj=0; kj<2; kj++)
	{
	    if (curr[kj] < sub1->startParam(kj) - ptol ||
		curr[kj] > sub1->endParam(kj) + ptol)
		break;
	}
	if (kj < 2)
	    continue;  // The singularity does not lie in surface 1

	for (kj=0; kj<2; kj++)
	{
	    if (curr[kj] < sub2->startParam(kj) - ptol ||
		curr[kj] > sub2->endParam(kj) + ptol)
		break;
	}
	if (kj < 2)
	    continue;  // The singularity does not lie in surface 2

	// The singularity lies in both surfaces
	for (kj=0; kj<2; kj++)
	    sing[kj] = sing[2+kj] = curr[kj];
	return true;
    }
    return false;  // No singularity lies in both surfaces
}


//===========================================================================
void SfSelfIntersector::estimateSingBox(shared_ptr<ParamSurfaceInt> normsf,
					double param[], 
					vector<pair<double, int> >& sing_box)
					
//===========================================================================
{
    double curv_tol = 50.0*epsge_->getEpsge();
    double dist_tol = 500.0*epsge_->getEpsge();

    // Parameter domain
    vector<double> mima = surf_->getMima();

    // In each parameter direction (positive and negative) from the singularity,
    // step until both the minimum curvature radius of the boundary curve is larger
    // than curv_tol and the minimum distance between the potential singularity
    // box boundary and the singularity is larger than dist_tol
    double step_u = 0.001*(mima[1] - mima[0]); 
    double step_v = 0.001*(mima[3] - mima[2]);
    double min_u = param[0] - step_u;
    double max_u = param[0] + step_u;
    double min_v = param[1] - step_v;
    double max_v = param[1] + step_v;
    bool found_min_u = false, found_max_u = false, found_min_v = false, found_max_v = false;
    int first_min_u = 0, first_max_u = 0, first_min_v = 0, first_max_v = 0;

    double minpar[2];
    double curvrad, dist;

    while (true) 
    {
	if (found_min_u);
	else if (min_u < mima[0])
	    {
		found_min_u = true;
		first_min_u = std::max(first_min_u, 1);
	    }
	else
	{
	    curvrad = getMinCurvatureRadAlongCurve(normsf, 0, min_u, min_v-step_v, 
						   max_v+step_v, minpar);
	    dist = getMinDistAlongCurve(param, 0, min_u, min_v, max_v);
	    if ((curvrad < 0 || curvrad >= curv_tol) && dist >= dist_tol)
		found_min_u = true;
	    else
		min_u -= step_u;
	    if (!first_min_u)
	    {
		if (curvrad < 0 || curvrad >= curv_tol)
		    first_min_u = 2;
		else if (dist >= dist_tol)
		    first_min_u = 3;
	    }
	    if (curvrad >= 0 && curvrad < curv_tol)
	    {
		if (minpar[1] < min_v)
		    min_v = minpar[1];
		if (minpar[1] > max_v)
		    max_v = minpar[1];
	    }
	}

	if (found_max_u);
	else if (max_u > mima[1])
	    {
		found_max_u = true;
		first_max_u = std::max(first_max_u, 1);
	    }
	else
	{
	    curvrad = getMinCurvatureRadAlongCurve(normsf, 0, max_u, min_v-step_v, 
						   max_v+step_v, minpar);
	    dist = getMinDistAlongCurve(param, 0, max_u, min_v, max_v);
	    if ((curvrad < 0 || curvrad >= curv_tol) && dist >= dist_tol)
		found_max_u = true;
	    else
		max_u += step_u;
	    if (!first_max_u)
	    {
		if (curvrad < 0 || curvrad >= curv_tol)
		    first_max_u = 2;
		else if (dist >= dist_tol)
		    first_max_u = 3;
	    }
	    if (curvrad >= 0 && curvrad < curv_tol)
	    {
		if (minpar[1] < min_v)
		    min_v = minpar[1];
		if (minpar[1] > max_v)
		    max_v = minpar[1];
	    }
	}

	if (found_min_v);
	else if (min_v < mima[2])
	{
	    found_min_v = true;
	    first_min_v = std::max(first_min_v, 1);
	}
	else
	{
	    curvrad = getMinCurvatureRadAlongCurve(normsf, 1, min_v, min_u-step_u, 
						   max_u+step_u, minpar);
	    dist = getMinDistAlongCurve(param, 1, min_v, min_u, max_u);
	    if ((curvrad < 0 || curvrad >= curv_tol) && dist >= dist_tol)
		found_min_v = true;
	    else
		min_v -= step_v;
	    if (!first_min_v)
	    {
		if (curvrad < 0 || curvrad >= curv_tol)
		    first_min_v = 2;
		else if (dist >= dist_tol)
		    first_min_v = 3;
	    }
	    if (curvrad >= 0 && curvrad < curv_tol)
	    {
		if (minpar[0] < min_u)
		    min_u = minpar[0];
		if (minpar[0] > max_u)
		    max_u = minpar[0];
	    }
	}

	if (found_max_v);
	else if (max_v > mima[3])
	{
	    found_max_v = true;
	    first_max_v = std::max(first_max_v, 1);
	}
	else
	{
	    curvrad = getMinCurvatureRadAlongCurve(normsf, 1, max_v, min_u-step_u, 
						   max_u+step_u, minpar);
	    dist = getMinDistAlongCurve(param, 1, max_v, min_u, max_u);
	    if ((curvrad < 0 || curvrad >= curv_tol) && dist >= dist_tol)
		found_max_v =  true;
	    else
		max_v += step_v;
	    if (!first_max_v)
	    {
		if (curvrad < 0 || curvrad >= curv_tol)
		    first_max_v = 2;
		else if (dist >= dist_tol)
		    first_max_v = 3;
	    }
	    if (curvrad >= 0 && curvrad < curv_tol)
	    {
		if (minpar[0] < min_u)
		    min_u = minpar[0];
		if (minpar[0] > max_u)
		    max_u = minpar[0];
	    }
	}

	if (found_min_u && found_max_u && found_min_v && found_max_v)
	    break;
    }

    sing_box[0] = make_pair(std::max(min_u, mima[0]), first_min_u);
    sing_box[1] = make_pair(std::min(max_u, mima[1]), first_max_u);
    sing_box[2] = make_pair(std::max(min_v, mima[2]), first_min_v);
    sing_box[3] = make_pair(std::min(max_v, mima[3]), first_max_v);
}

//===========================================================================
void SfSelfIntersector::makeSingularityUnions()
//===========================================================================
{
    sing_union_.clear();  // Start from scratch

    if (sing_box_.size() == 0)
	return;  // No singularities

    size_t ki, kj;
    vector<int> idx;
    double limit[4];
    if (sing_box_.size() == 1)
    {
	// One box and one union
	vector<int> idx(1, 0);
	for (int kr=0; kr<4; kr++)
	    limit[kr] = sing_box_[0].box_limit_[kr].first;
	SingUnion curr_union(limit, idx);
	sing_union_.push_back(curr_union);
    }  
    else
    {
	vector<int> perm(sing_box_.size());
	for (ki=0; ki<perm.size(); ki++)
	    perm[ki] = (int)ki;

	for (ki=0; ki<perm.size(); ki++)
	{
	    vector<int> idx;
	    for (int kr=0; kr<4; kr++)
		limit[kr] = sing_box_[perm[ki]].box_limit_[kr].first;
	    idx.push_back(perm[ki]);

	    for (kj=ki+1; kj<perm.size(); kj++)
	    {
		// Check if the current singularity box overlaps the current union
		if (sing_box_[perm[kj]].box_limit_[0].first > limit[1])
		    continue;  // Not an overlap
		if (sing_box_[perm[kj]].box_limit_[1].first < limit[0])
		    continue;  // Not an overlap
		if (sing_box_[perm[kj]].box_limit_[2].first > limit[3])
		    continue;  // Not an overlap
		if (sing_box_[perm[kj]].box_limit_[3].first < limit[2])
		    continue;  // Not an overlap

		// An overlap is found. Increase the union
		limit[0] = std::min(limit[0], sing_box_[perm[kj]].box_limit_[0].first);
		limit[1] = std::max(limit[1], sing_box_[perm[kj]].box_limit_[1].first);
		limit[2] = std::min(limit[2], sing_box_[perm[kj]].box_limit_[2].first);
		limit[3] = std::max(limit[3], sing_box_[perm[kj]].box_limit_[3].first);
		idx.push_back(perm[kj]);

		// Reorder permutation array
		if (kj > ki+1)
		{
		    size_t dum = perm[ki+1];
		    perm[ki+1] = perm[kj];
		    perm[kj] = (int)dum;
		}

		ki++;
		kj = ki;
	    }

	    // Union found
	    SingUnion curr_union(limit, idx);
	    sing_union_.push_back(curr_union);
	}
    }
	
    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
    {
	ofstream debug1("union_block.dsp");
	debug1 << "fg: gray" << endl;
	for (ki=0; ki<sing_union_.size(); ki++)
	{
	    debug1 << "lin:" << endl;
	    debug1 << sing_union_[ki].limit_[0] << " " << sing_union_[ki].limit_[2] << endl;
	    debug1 << sing_union_[ki].limit_[1] << " " << sing_union_[ki].limit_[2] << endl;
	    debug1 << "lin:" << endl;
	    debug1 << sing_union_[ki].limit_[1] << " " << sing_union_[ki].limit_[2] << endl;
	    debug1 << sing_union_[ki].limit_[1] << " " << sing_union_[ki].limit_[3] << endl;
	    debug1 << "lin:" << endl;
	    debug1 << sing_union_[ki].limit_[1] << " " << sing_union_[ki].limit_[3] << endl;
	    debug1 << sing_union_[ki].limit_[0] << " " << sing_union_[ki].limit_[3] << endl;
	    debug1 << "lin:" << endl;
	    debug1 << sing_union_[ki].limit_[0] << " " << sing_union_[ki].limit_[3] << endl;
	    debug1 << sing_union_[ki].limit_[0] << " " << sing_union_[ki].limit_[2] << endl;
	}
    }
}


//===========================================================================
double SfSelfIntersector::getMinCurvatureRadAlongCurve(shared_ptr<ParamSurfaceInt> normsf, 
						       int dir, double par, 
						       double tmin, double tmax,
						       double param[])
//===========================================================================
{
    Point origo(0.0, 0.0, 0.0);
    double minpar, mindist;
    Point minval;

    // Find minium length of surface normal along the specified curve
    normsf->minimumAlongCurve(dir, par, tmin, tmax, origo, minpar, minval, mindist);

    // Compute max- and min curvatures in the found point in the surface corresponding
    // to the normal surface
    param[dir] = par;
    param[1-dir] = minpar;

    shared_ptr<ParamSurface> sf = surf_->getParamSurface();
    double k1, k2; // Minium and maximum curvature 
    Point d1, d2;  // Directions corrsponding to the principal curvatures
    CurvatureAnalysis::principalCurvatures(*(sf.get()), param[0], param[1], k1, d1, k2, d2);

    double max_curv = std::max(fabs(k1), fabs(k2));
    double curv_rad = (max_curv < epsge_->getNumericalTol()) ? -1 : 1.0/max_curv;
    return curv_rad;
}

//===========================================================================
int SfSelfIntersector::repairIntersections()
//===========================================================================
{
    // Purpose: This function may be called at the end of the
    // intersection algorithm to repair incorrect holes, branch
    // points, etc., in the intersection curves.
    SfSfIntersector sfsfintersect(surf_, surf_, epsge_, this);

    sfsfintersect.removeIsolatedPoints();

     sfsfintersect.getIntPool()->removeLooseEnds();

     return 0;
}


//===========================================================================
double SfSelfIntersector::getMinDistAlongCurve(double param[], int dir, double par, 
					       double tmin, double tmax)
//===========================================================================
{
    // Evaluate the current surface in the given parameter value
    Point pnt;
    surf_->point(pnt, param);

    // Find the closest point to the given point along the given curve
    double minpar, mindist;
    Point minval;

    surf_->minimumAlongCurve(dir, par, tmin, tmax, pnt, minpar, minval, mindist);
    return mindist;
}

//===========================================================================
int SfSelfIntersector::isInComplexDomain(shared_ptr<ParamSurfaceInt> subsf)
//===========================================================================
{
    double ta1 = subsf->startParam(0);
    double tb1 = subsf->endParam(0);
    double ta2 = subsf->startParam(1);
    double tb2 = subsf->endParam(1);

    for (int ki=0; ki<(int)complex_domain_.size(); ki++)
    {
	double u1 = complex_domain_[ki].umin();
	double u2 = complex_domain_[ki].umax();
	double v1 = complex_domain_[ki].vmin();
	double v2 = complex_domain_[ki].vmax();

	if (ta1 >= u1 && tb1 <= u2 && ta2 >= v1 && tb2 <= v2)
	    return ki;
    }

    return -1;
}

//===========================================================================
void SfSelfIntersector::writeDebugComplex(int file, vector<RectDomain>& domain)
//===========================================================================
{
    // DEBUG
    if (domain.size() == 0)
	return;

    ofstream debug((file==1) ? "self_complex.dsp" : "single_block.dsp");
    debug << "fg: red" << endl;
    for (size_t ki=0; ki<domain.size(); ki++)
    {
	double u1 = domain[ki].umin();
	double u2 = domain[ki].umax();
	double v1 = domain[ki].vmin();
	double v2 = domain[ki].vmax();
	debug << "lin:" << endl;
	debug << u1 << " " << v1 << endl;
	debug << u2 << " " << v1 << endl;
	debug << u2 << " " << v2 << endl;
	debug << u1 << " " << v2 << endl;
	debug << u1 << " " << v1 << endl;
    }

}

//===========================================================================
