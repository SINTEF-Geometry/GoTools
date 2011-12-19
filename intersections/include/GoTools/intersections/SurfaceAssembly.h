//===========================================================================
//                                                                           
// File: SurfaceAssembly.h
//                                                                           
// Created: 30.05.05
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SURFACEASSEMBLY_H
#define _SURFACEASSEMBLY_H

#include "GoTools/intersections/ParamSurfaceInt.h"

namespace Go
{

class SurfaceAssembly
    {
    public:

	/// Constructor
	SurfaceAssembly(shared_ptr<ParamSurfaceInt> surf,
			std::vector<std::pair<double,int> > u_div,
			std::vector<std::pair<double,int> > v_div,
			std::vector<RectDomain> sing_domain,
			double rel_par_res);

	/// Destructor
	~SurfaceAssembly();

	void resetSubIndex(int idx=0)
	    { idx_sub_ = idx; }

	void resetAssemblyIndex()
	    { idx_assembly_ = 0; }

	int getNmbSubSurface();


	bool getNextSubSurface(shared_ptr<ParamSurfaceInt>& sub_sf,
			       int& idx, int& sing_idx);

	bool getNextAssembly(shared_ptr<ParamSurfaceInt>& assembly,
			     int& idx, bool& potential_sing);

	/// This function makes a check on whether two sub surfaces are
	/// neighbours, but if the surfaces meet at a singularity
	/// they are NOT classified as neighbours
	bool subSfNeighbour(int idx1, int idx2);

	bool doTouch(int idx1, int idx2);

	bool touchAtSingularity(int idx1, int idx2);

	bool isInPrevAssembly(int idx1, int idx2);

	bool isInFirstAssembly(shared_ptr<ParamSurfaceInt> sub_srf);

	int getSubSurfaceIndex(shared_ptr<ParamSurfaceInt> sub_srf,
			       bool& at_end);

	std::vector<std::pair<double,int> > getUdiv()
	    {
		return u_div_;
	    }

	std::vector<std::pair<double,int> > getVdiv()
	    {
		return v_div_;
	    }

    private:
	shared_ptr<ParamSurfaceInt> surf_;
	std::vector<std::pair<double,int> > u_div_;
	std::vector<std::pair<double,int> > v_div_;
	std::vector<RectDomain> sing_domain_;
	double ptol_;
	bool closed_in_u_;
	bool closed_in_v_;
	mutable int idx_sub_;
	mutable int idx_assembly_;  // The index corresponds to the
	// sub surface where the lower left corner of the assembly lies.
	// Thus, the index does not count the total number of assemblies

	// Refine the surface according to the given division parameters
	// if possible (to save time at a later stage)
	void refineSurf();
    };
} // namespace Go

#endif  // _SURFACEASSEMBLY_H



