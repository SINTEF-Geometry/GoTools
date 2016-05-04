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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/BoundedUtils.h"

#include <fstream>
#include <assert.h>

using namespace Go;
using std::cout;
using std::endl;
using std::vector;

bool fixParCvCrossingCylinderSeem(BoundedSurface* trimmed_cyl);


int main(int argc, char *argv[])
{
    if (argc != 3)
    {
	std::cout << "Usage: sfs_file.g2 repaired_sfs_file.g2" << std::endl;
	return -1;
    }

    std::ifstream filein(argv[1]); // Input bd sfs (may contain other objects).
    std::ofstream fileout(argv[2]); // Fixed bd sfs (and unaltered other objects).

    // For BoundedSurface we may choose to recreate all the boundary parameter curves.
    const bool recreate_par_cvs = true;
    if (recreate_par_cvs)
    {
	cout << "Recreating all parameter curves for CurveOnSurface." << endl;
    }

    // Create the default factory
    GoTools::init();

    ObjectHeader header;
    int num_bd_sfs = 0;
    int num_bd_sfs_fixed = 0;
    int num_bd_sfs_fix_failed = 0;
    int obj_id = 0;
    vector<shared_ptr<GeomObject> > objs;
    while (filein)
    {
//	std::cout << "Object number: " << objs.size() << std::endl;
	try {
	    header.read(filein);
	}
	catch (...)
	{
	    MESSAGE("Failed reading the Header!");
	    break; // Assuming we are either done or the rest of the file is garbage ...
	}

	shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
	try
	{
	    if (geom_obj->instanceType() == Class_BoundedSurface)
	    {
		shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(geom_obj);
		bool fix_trim_cvs = false; // We do not want to fix trim cvs from the read routine.
		bd_sf->read(filein, fix_trim_cvs);

                if (recreate_par_cvs) {
                    std::vector<CurveLoop> all_loops = bd_sf->allBoundaryLoops();
                    for (size_t ki = 0; ki < all_loops.size(); ++ki) {
                        for (size_t kj = 0; kj < all_loops[ki].size(); ++kj) {
                            shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(all_loops[ki][kj]);
                            if (cv_on_sf.get() != NULL) {
                                cv_on_sf->unsetParameterCurve();
                            }
                        }
                    }
                }
	    }
	    else
	    {
		geom_obj->read(filein);
	    }
	    objs.push_back(geom_obj);
	}
	catch (...)
	{
	    MESSAGE("Failed reading the GeomObject!");
	    // We mush back even objects which failed
	    objs.push_back(shared_ptr<GeomObject>(NULL));
	}

	Utils::eatwhite(filein);
    }

    for (int kk = 0; kk < objs.size(); ++kk)
    {
	shared_ptr<GeomObject> geom_obj = objs[kk];
	if (geom_obj.get() == NULL)
	{
	    MESSAGE("Missing object " << kk << "!");
	    continue;
	}

	if (geom_obj->instanceType() == Class_BoundedSurface)
	{
	    ++num_bd_sfs;
	    BoundedSurface* bd_sf = dynamic_cast<BoundedSurface*>(geom_obj.get());
	    double epsgeo = bd_sf->getEpsGeo(); // The smallest for all the loops.
	    int valid_state = 0;
	    bool is_valid = bd_sf->isValid(valid_state);

#ifndef NDEBUG
	    std::ofstream debug("tmp/debug.g2");
	    ParamSurface* under_sf = bd_sf->underlyingSurface().get();
	    under_sf->writeStandardHeader(debug);
	    under_sf->write(debug);
	    for (int ki = 0; ki < bd_sf->numberOfLoops(); ++ki)
	    {
		shared_ptr<CurveLoop> loop = bd_sf->loop(ki);
		for (size_t kj = 0; kj < loop->size(); ++kj)
		{
		    shared_ptr<ParamCurve> cv = (*loop)[kj];
		    if (cv->instanceType() == Class_CurveOnSurface)
		    {
			shared_ptr<CurveOnSurface> cv_on_sf =
			    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
            if (cv_on_sf->parameterCurve().get() != NULL) {
			    shared_ptr<SplineCurve> pcv =
				dynamic_pointer_cast<SplineCurve, ParamCurve>
				(cv_on_sf->parameterCurve());
			    if (pcv.get() != NULL)
				SplineDebugUtils::writeSpaceParamCurve(*pcv, debug, 0.0);
			    else
			    {
				cv_on_sf->parameterCurve()->writeStandardHeader(debug);
				cv_on_sf->parameterCurve()->write(debug);
			    }
			}
            if (cv_on_sf->spaceCurve().get() != NULL)
			{
			    cv_on_sf->spaceCurve()->writeStandardHeader(debug);
			    cv_on_sf->spaceCurve()->write(debug);
			}
		    }
		    else
		    {
			cv->writeStandardHeader(debug);
			cv->write(debug);
		    }
		}
	    }
	    double debug_val = 0.0;
#endif

	    if (valid_state < 0)//== -2)
	    {
#if 0
		if (bd_sf->underlyingSurface()->instanceType() == Class_Cylinder)
		{ // Special treatment checking for par cvs crossing the seem.
		    bool fixed = fixParCvCrossingCylinderSeem(bd_sf);
		    if (fixed)
		    {
			MESSAGE("A parameter curve crossed the seem, fixed!");
			is_valid = bd_sf->isValid(valid_state);
			if (is_valid)
			{
			    cout << "Success! obj_id = " << kk << endl;
			    ++num_bd_sfs_fixed;
			    continue;
			}
		    }
		}
#endif

		bool proj_ok = Go::BoundedUtils::createMissingParCvs(*bd_sf);
		if (proj_ok)
		{
		    ;//MESSAGE("Succeded in projecting the space curves.");
		}
		else
		{
		    MESSAGE("Failed projecting the space curves.");
		}
		bd_sf->analyzeLoops();
		is_valid = bd_sf->isValid(valid_state);
		if (is_valid)
		{
                    ++num_bd_sfs_fixed;
		}
	    }

	    if (!is_valid)
	    {
		// MESSAGE("Trying to fix the BoundedSurface!");

		double max_loop_gap = -1.0;
		double max_tol_mult = 1.0e03;
		bool success = bd_sf->fixInvalidSurface(max_loop_gap, max_tol_mult);
		if (success)
		{
//		    cout << "Success! obj_id = " << kk << endl;
		    ++num_bd_sfs_fixed;
		}
		else
		{
		    is_valid = bd_sf->isValid(valid_state);
		    MESSAGE("Failed fixing bd_sf! Status: " << valid_state <<
			    ". Object id = " << kk << ", max_loop_gap: " << max_loop_gap);
		    ++num_bd_sfs_fix_failed;

		    if (valid_state == -1)
		    { // This means we may need to replace the epsgeo with a larger value.
			// @@sbr201310 Replace epsgeo with larger value!
			double epsgeo = bd_sf->getEpsGeo();
			// We use the same tolerance for all the loops.
			int num_loops = bd_sf->numberOfLoops();
			double global_max_loop_sf_dist = 0.0;
			int num_samples = 100; // @@sbr201310 The number of samples should not be a user input.
			for (int ki = 0; ki < num_loops; ++ki)
			{
			    double max_loop_sf_dist = bd_sf->maxLoopSfDist(ki, num_samples);
			    if (max_loop_sf_dist > global_max_loop_sf_dist)
			    {
				global_max_loop_sf_dist = max_loop_sf_dist;
			    }
			}

			const double epsgeo_ratio = 1000.0;
			double new_epsgeo = epsgeo_ratio*global_max_loop_sf_dist;
			std::cout << "epsgeo: " << epsgeo << ", new_epsgeo: " << new_epsgeo << std::endl;
			for (int ki = 0; ki < num_loops; ++ki)
			{
			    if (bd_sf->loop(ki)->getSpaceEpsilon() < new_epsgeo)
			    {
				bd_sf->loop(ki)->setSpaceEpsilon(new_epsgeo);
			    }
			}
			// /// Get a shared pointer to a specific boundary loop
			// shared_ptr<CurveLoop> loop(int idx)
			// { return boundary_loops_[idx]; }

			// std::vector<CurveLoop> all_bd_loops = bd_sf->absolutelyAllBoundaryLoops();
			// for (size_t ki = 0; ki < all_bd_loops.size(); ++ki)
			// {
			//     if (all_bd_loops[ki].getSpaceEpsilon() < new_epsgeo)
			//     {
			// 	all_bd_loops[ki].setSpaceEpsilon(new_epsgeo);
			//     }
			// }
			// shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
			// bd_sf = shared_ptr<BoundedSurface>
			//     (new BoundedSurface(under_sf, all_bd_loops, new_epsgeo));
			bd_sf->analyzeLoops();
			is_valid = bd_sf->isValid(valid_state);
			std::cout << "valid_state after tolerance change: " << valid_state << std::endl;
			if (is_valid)
			{
			    --num_bd_sfs_fix_failed;
			    ++num_bd_sfs_fixed;
			}
		    }
		    else
		    {
#if 0
			// We write to file the underlying surfaces and the space curves.
			shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
			under_sf->writeStandardHeader(fileout);
			under_sf->write(fileout);
			vector<CurveLoop> bd_loops = bd_sf->allBoundaryLoops();
			for (size_t ki = 0; ki < bd_loops.size(); ++ki)
			{
			    for (size_t kj = 0; kj < bd_loops[ki].size(); ++kj)
			    {
				shared_ptr<ParamCurve> par_cv = bd_loops[ki][kj];
				if (par_cv->instanceType() == Class_CurveOnSurface)
				{
				    CurveOnSurface* cv_on_sf = dynamic_cast<CurveOnSurface*>(par_cv.get());
				    if (cv_on_sf->spaceCurve())
				    {
					cv_on_sf->spaceCurve()->writeStandardHeader(fileout);
					cv_on_sf->spaceCurve()->write(fileout);
				    }
				    else
				    {
					MESSAGE("Missing space curve!");
				    }
				}
				else
				{
				    MESSAGE("Unexpected curve type!");
				}
			    }
			}
#endif
		    }
		}
	    }
	    // else
	    // {
	    // 	bd_sf->writeStandardHeader(fileout);
	    // 	bd_sf->write(fileout);
	    // }
	}
	// else
	// {
	//     geom_obj->writeStandardHeader(fileout);
	//     geom_obj->write(fileout);
	// }

    }


    int num_valid_bd_sf = 0;
    int num_invalid_bd_sf = 0;
    for (int kk = 0; kk < objs.size(); ++kk)
    {
	objs[kk]->writeStandardHeader(fileout);
	objs[kk]->write(fileout);
	if (objs[kk]->instanceType() == Class_BoundedSurface)
	{
	    shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface>(objs[kk]);
	    bd_sf->analyzeLoops();
	    int valid_state = 0;
	    bool is_valid = bd_sf->isValid(valid_state);
	    if (is_valid)
	    {
		std::cout << "Valid object id: " << kk << ", instance type: " <<
                    bd_sf->underlyingSurface()->instanceType() << std::endl;
		++num_valid_bd_sf;
	    }
	    else
	    {
		std::cout << "Invalid object id: " << kk << ", instance type: " <<
                    bd_sf->underlyingSurface()->instanceType() << std::endl;
		++num_invalid_bd_sf;
	    }
	}
    }

    std::cout << "num_bd_sfs: " << num_bd_sfs << ", num_bd_sfs_fixed: " << num_bd_sfs_fixed <<
	", num_bd_sfs_fix_failed: " << num_bd_sfs_fix_failed << ", num_valid: " << num_valid_bd_sf <<
	", num_invalid: " << num_invalid_bd_sf << std::endl;

}


// @@sbr This method should be moved to BoundedSurface. The method is
// not that straight forward when handling other closed sfs as we may
// not assume that the surface is cyclic with a common period. For
// these cases we may have to split the surface up into smaller
// pieces. And then it makes sense to use an external routine
// (i.e. not a BoundedSurface member function).
bool fixParCvCrossingCylinderSeem(BoundedSurface* trimmed_cyl)
{
    MESSAGE("Under construction!");
    if (trimmed_cyl->underlyingSurface()->instanceType() != Class_Cylinder)
    {
	return false;
    }

    Cylinder* cyl = dynamic_cast<Cylinder*>(trimmed_cyl->underlyingSurface().get());
    bool params_swapped = cyl->isSwapped();
    RectDomain rect_dom = cyl->containingDomain();
    // We may not assume that the cylinder is parametrized on [0, 2*M_PI).


    // We run through all the bd_cvs, extracting the convex hull of
    // the coefs of the parameter cvs. Comparing this with the domain
    // of the cylinder we can decide if we should and are able to move
    // the seem.  We only need to look at the outer loop (the first)
    // to test for crossing of seem and possibly the translation
    // vector.
    shared_ptr<CurveLoop> outer_loop = trimmed_cyl->loop(0);
    vector<BoundingBox> par_bd_boxes;
    vector<ParamCurve*> par_cvs(outer_loop->size());

#ifndef NDEBUG
    CurveOnSurface* first_cv_on_sf = dynamic_cast<CurveOnSurface*>((*outer_loop)[0].get());
    assert(first_cv_on_sf != NULL);
    if (first_cv_on_sf->parameterCurve() && first_cv_on_sf->spaceCurve())
    {
	Point space_cv_pt = first_cv_on_sf->spaceCurve()->point(first_cv_on_sf->spaceCurve()->startparam());
	Point par_cv_pt =
	    first_cv_on_sf->parameterCurve()->point(first_cv_on_sf->parameterCurve()->startparam());
	Point sf_pt = cyl->ParamSurface::point(par_cv_pt[0], par_cv_pt[1]);
    }
#endif

    for (size_t ki = 0; ki < outer_loop->size(); ++ki)
    {
	CurveOnSurface* cv_on_sf = dynamic_cast<CurveOnSurface*>((*outer_loop)[ki].get());
	assert(cv_on_sf != 0);
	ParamCurve* par_cv = cv_on_sf->parameterCurve().get();
	par_cvs[ki] = par_cv;
	if (par_cv != 0)
	{
	    BoundingBox par_bd_box(2);
	    SplineCurve* spline_cv = par_cv->geometryCurve();
	    Point lin_dir;
	    double eps = 1e-05;
	    if (spline_cv != NULL)
	    {
		// Not handling rational cases at the moment.
		assert(!spline_cv->rational());
//		vector<double>
		par_bd_box.setFromArray(spline_cv->coefs_begin(), spline_cv->coefs_end(), 2);
	    }
	    else if (par_cv->isLinear(lin_dir, eps))
	    {
		vector<Point> pts(2);
		pts[0] = par_cv->point(par_cv->startparam());
		pts[1] = par_cv->point(par_cv->endparam());
		par_bd_box.setFromPoints(pts);
	    }
	    else
	    {
		MESSAGE("Case not handled!");
		return false;
	    }
	    par_bd_boxes.push_back(par_bd_box);
	}
    }

    // We then run through all the bd_box elements, checking if they line inside the RectDomain of the cylinder.
    // We may assume that the cylinder is closed in the u-direction.
    bool closed_dir_u, closed_dir_v;
    cyl->isClosed(closed_dir_u, closed_dir_v);
    if (params_swapped)
    {
	std::swap(closed_dir_u, closed_dir_v);
    }
    assert(closed_dir_u && (!closed_dir_v));
    double umin = rect_dom.umin();
    double umax = rect_dom.umax();
    double u_low = umax;
    double u_high = umin;
    for (size_t ki = 0; ki < par_bd_boxes.size(); ++ki)
    {
	Point low = par_bd_boxes[ki].low();
	if (low[0] < u_low)
	{
	    u_low = low[0];
	}
	Point high = par_bd_boxes[ki].high();
	if (high[0] > u_high)
	{
	    u_high = high[0];
	}
    }

    if ((u_low < umin) || (u_high > umax))
    {
	MESSAGE("A parameter curve lies outside the domain of the cylinder! u_low = " <<
		u_low << ", u_high = " << u_high);

	// We rotate the cylinder (parametrization) such that the seem does not cross a parameter curve.
//	double transl_u = (u_low < umin) ? umin - u_low : umax - u_high;
	double u_span = u_high - u_low;
	double new_u_low = umax;
	double new_u_high = umin;
	vector<double> transl_u(par_bd_boxes.size(), 0.0);
	if (u_span > umax - umin)
	{ // For some cases the only option is to split the bd_sf into multiple sfs.
	    // We could try to move the segment by 2*M_PI and then move the seem on the other side.
	    for (size_t ki = 0; ki < par_bd_boxes.size(); ++ki)
	    {
		Point low = par_bd_boxes[ki].low();
		Point high = par_bd_boxes[ki].high();
		if (low[0] < (umin - u_low))
		{
		    cout << "ki = " << ki << ", translate = " << 2*M_PI << endl;
		    transl_u[ki] = 2*M_PI;
		    if (low[0] + 2*M_PI < new_u_low)
		    {
			new_u_low = low[0] + 2*M_PI;
		    }
		    if (high[0] + 2*M_PI > new_u_high)
		    {
			new_u_high = high[0] + 2*M_PI;
		    }
		}
		else
		{
		    if (low[0] < new_u_low)
		    {
			new_u_low = low[0];
		    }
		    if (high[0] > new_u_high)
		    {
			new_u_high = high[0];
		    }
		}
	    }
	}

	u_span = new_u_high - new_u_low;
	if (u_span > 2*M_PI)
	{
	    // We next run through all cvs setting the common translation.
	    MESSAGE("Moving of seem is not enough to handle this case, u_span = " << u_span);
	}
	// new_u_low will be the value used for moving the seem.
	for (size_t ki = 0; ki < par_bd_boxes.size(); ++ki)
	{
	    transl_u[ki] -= new_u_low;
	    cout << "transl_u[ki] = " << transl_u[ki] << endl;
	}

	// We then run through the par_cvs translating the u-values of the coefs.
	cout << "DEBUG: Soon we will handle this case!" << endl;
	// We rotate the cylinder by u_low - new_u_low.
	double rot_ang_deg = new_u_low;
	cyl->rotate(rot_ang_deg);
	// Finally we run through all curve segments and perfomr the translation in the u-dir.
	for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	{
	    if (par_cvs[ki] != NULL)
	    {
		if (par_cvs[ki]->instanceType() == Class_SplineCurve)
		{
		    SplineCurve* spline_cv = dynamic_cast<SplineCurve*>(par_cvs[ki]);
		    vector<double>::iterator iter = spline_cv->coefs_begin();
		    while (iter != spline_cv->coefs_end())
		    {
			iter[0] += transl_u[ki];
			iter += 2;
		    }
		}
		else if (par_cvs[ki]->instanceType() == Class_Line)
		{
		    Line* line = dynamic_cast<Line*>(par_cvs[ki]);
		    Point transl_vec(2, 0.0);
		    transl_vec[0] = transl_u[ki];
		    line->translateCurve(transl_vec);
		}
		else
		{
		    THROW("Unsupported curve type: " << par_cvs[ki]->instanceType());
		}
	    }
	}

#ifndef NDEBUG
	std::ofstream debug("tmp/debug.g2");
	cyl->writeStandardHeader(debug);
	cyl->write(debug);
	for (size_t ki = 0; ki < outer_loop->size(); ++ki)
	{
	    shared_ptr<ParamCurve> cv = (*outer_loop)[ki];
	    if (cv->instanceType() == Class_CurveOnSurface)
	    {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
        if (cv_on_sf->parameterCurve().get() != NULL) {
		    shared_ptr<SplineCurve> pcv =
			dynamic_pointer_cast<SplineCurve, ParamCurve>
			(cv_on_sf->parameterCurve());
		    if (pcv.get() != NULL)
			SplineDebugUtils::writeSpaceParamCurve(*pcv, debug, 0.0);
		    else
		    {
			cv_on_sf->parameterCurve()->writeStandardHeader(debug);
			cv_on_sf->parameterCurve()->write(debug);
		    }
		}
        if (cv_on_sf->spaceCurve().get() != NULL)
		{
		    cv_on_sf->spaceCurve()->writeStandardHeader(debug);
		    cv_on_sf->spaceCurve()->write(debug);
		}
	    }
	    else
	    {
		cv->writeStandardHeader(debug);
		cv->write(debug);
	    }
	}
	double debug_val = 0.0;
#endif

#ifndef NDEBUG
	// CurveOnSurface* first_cv_on_sf = dynamic_cast<CurveOnSurface*>((*outer_loop)[0].get());
	// assert(first_cv_on_sf != NULL);
	Point new_space_cv_pt = first_cv_on_sf->spaceCurve()->point(first_cv_on_sf->spaceCurve()->startparam());
	Point new_par_cv_pt =
	    first_cv_on_sf->parameterCurve()->point(first_cv_on_sf->parameterCurve()->startparam());
	Point new_sf_pt = cyl->ParamSurface::point(new_par_cv_pt[0], new_par_cv_pt[1]);
#endif

	trimmed_cyl->analyzeLoops();
	return true;
    }

    return false;
}
