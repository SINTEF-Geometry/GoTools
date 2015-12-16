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

#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include <fstream>


using std::vector;

namespace Go
{

//===========================================================================
double BoundedSurface::area(double tol) const
//===========================================================================
{
    double fac = 10.0;

    // Get surrounding domain
    RectDomain domain = containingDomain();
    
    // Get smallest surrounding surface
    shared_ptr<ParamSurface> base_sf = surface_;
    while (base_sf->instanceType() == Class_BoundedSurface)
	base_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(base_sf)->underlyingSurface();
	
    
    vector<shared_ptr<ParamSurface> > sfs = 
	base_sf->subSurfaces(domain.umin(), domain.vmin(), domain.umax(), domain.vmax());

    double total_area = 0.0;
    size_t kr;
    for (kr=0; kr<sfs.size(); ++kr)
    {
	shared_ptr<SplineSurface> spline_sf = 
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(sfs[kr]);
	if (spline_sf.get())
	    total_area += spline_sf->area(tol);
    }
    if (isIsoTrimmed(0.1*tol))
    {
	// We are finished
	return total_area;
    }

    std::ofstream out_file("tmp_mini_surf.g2");
    // Otherwise, split the current surface into smaller pieces and compare areas
    int nmb_split = 5;
    double u1 = domain.umin();
    double v1 = domain.vmin();
    double u_del = (domain.umax() - u1)/(double)(nmb_split);
    double v_del = (domain.vmax() - v1)/(double)(nmb_split);
    
    int ki, kj;
    size_t kk;
    double curr_area = 0.0;
    vector<shared_ptr<ParamSurface> > all_sfs;
    for (kj=0; kj<nmb_split; ++kj, v1+=v_del)
	for (ki=0, u1=domain.umin(); ki<nmb_split; ++ki, u1+=u_del)
	{
	    vector<shared_ptr<ParamSurface> > sub_sfs = subSurfaces(u1, v1, u1+u_del, v1+v_del);
	    for (kr=0; kr<sub_sfs.size(); ++kr)
	    {
// 		sub_sfs[ki]->writeStandardHeader(out_file);
// 		sub_sfs[ki]->write(out_file);
		RectDomain dom2 = sub_sfs[kr]->containingDomain();
		vector<shared_ptr<ParamSurface> > tmp_sfs = base_sf->subSurfaces(std::min(u1,dom2.umin()),
										 std::min(v1,dom2.vmin()),
										 std::max(u1+u_del,dom2.umax()),
										 std::max(v1+v_del,dom2.vmax()));
		for (kk=0; kk<tmp_sfs.size(); ++kk)
		    curr_area += tmp_sfs[kk]->area(tol);
	    }
	    all_sfs.insert(all_sfs.end(), sub_sfs.begin(), sub_sfs.end());
	}


    if (total_area/curr_area < 1.0 + fac*tol)
	return curr_area;  // The area is close enough

    // Compute recursively
    curr_area = 0.0;
    for (kr=0; kr<all_sfs.size(); ++kr)
	curr_area += all_sfs[kr]->area(tol);

    return curr_area;
}

} // namspace Go

