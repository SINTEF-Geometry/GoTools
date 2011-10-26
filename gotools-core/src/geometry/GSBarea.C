//===========================================================================
//
// File : GSBarea.C
//
// Created: Dec. 08
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include <fstream>


using std::shared_ptr;
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
	base_sf = std::dynamic_pointer_cast<BoundedSurface, ParamSurface>(base_sf)->underlyingSurface();
	
    
    vector<shared_ptr<ParamSurface> > sfs = 
	base_sf->subSurfaces(domain.umin(), domain.vmin(), domain.umax(), domain.vmax());

    double total_area = 0.0;
    size_t kr;
    for (kr=0; kr<sfs.size(); ++kr)
    {
	shared_ptr<SplineSurface> spline_sf = 
	    std::dynamic_pointer_cast<SplineSurface, ParamSurface>(sfs[kr]);
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
		curr_area += tmp_sfs[0]->area(tol);
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

