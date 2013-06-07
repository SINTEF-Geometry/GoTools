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

#ifndef _SURFACEASSEMBLY_H
#define _SURFACEASSEMBLY_H

#include "GoTools/intersections/ParamSurfaceInt.h"

namespace Go
{

  /// Collection of sub surfaces used in computation of surface
  /// self intersections
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



