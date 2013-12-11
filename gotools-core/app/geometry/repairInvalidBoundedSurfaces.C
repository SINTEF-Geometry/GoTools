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

#include <fstream>

using namespace Go;


int main(int argc, char *argv[])
{
    if (argc != 3)
    {
	std::cout << "Usage: sfs_file (.g2) repaired_sfs_file (.g2)" << std::endl;
	return -1;
    }

    std::ifstream filein(argv[1]); // Input bd sfs (may contain other objects).
    std::ofstream fileout(argv[2]); // Fixed bd sfs (and unaltered other objects).

    // Create the default factory
    GoTools::init();

    ObjectHeader header;
    int cntr = 0;
    int num_bd_sfs = 0;
    int num_bd_sfs_fixed = 0;
    int num_bd_sfs_fix_failed = 0;
    while (filein)
    {
	std::cout << "Object number: " << cntr << std::endl;
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
	    geom_obj->read(filein);
	}
	catch (...)
	{
	    MESSAGE("Failed reading the GeomObject!");
	}
	if (geom_obj->instanceType() == Class_BoundedSurface)
	{
	    ++num_bd_sfs;
	    BoundedSurface* bd_sf = dynamic_cast<BoundedSurface*>(geom_obj.get());
	    double epsgeo = 1.5e-02;// bd_sf->getEpsGeo(); // The smallest for all the loops.
	    int valid_state = 0;
	    bool is_valid = bd_sf->isValid(valid_state);

	    if (valid_state == -2)
	    {
		std::vector<CurveLoop> bd_loops = bd_sf->allBoundaryLoops();		
		//We must create missing parameter curves project space curves).
		for (size_t ki = 0; ki < bd_loops.size(); ++ki)
		{
		    for (size_t kj = 0; kj < bd_loops[ki].size(); ++kj)
		    {
			try
			{
			    if (bd_loops[ki][kj]->instanceType() == Class_CurveOnSurface)
			    {
				CurveOnSurface* cv_on_sf = dynamic_cast<CurveOnSurface*>(bd_loops[ki][kj].get());
				cv_on_sf->ensureParCrvExistence(epsgeo);
			    }
			}
			catch (...)
			{
			    MESSAGE("Failed projecting space curve!");
			}
		    }
		}
		is_valid = bd_sf->isValid(valid_state);
		if (is_valid)
		{
		    MESSAGE("Success!");
		}
	    }

	    if (!is_valid)
	    {
		MESSAGE("Trying to fix the BoundedSurface!");

		double max_loop_gap = -1.0;
		bool success = bd_sf->fixInvalidSurface(max_loop_gap);
		if (success)
		{
		    MESSAGE("Success!");
		    ++num_bd_sfs_fixed;
		    bd_sf->writeStandardHeader(fileout);
		    bd_sf->write(fileout);
		}
		else
		{
		    is_valid = bd_sf->isValid(valid_state);
		    MESSAGE("Failed fixing the BoundedSurface! Status: " << valid_state <<
			    ". Writing underlying surface and the space curves.");
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

			double new_epsgeo = 1.1*global_max_loop_sf_dist;
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
			bd_sf->writeStandardHeader(fileout);
			bd_sf->write(fileout);
		    }
		    else
		    {
			// We write to file the underlying surfaces and the space curves.
			shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
			under_sf->writeStandardHeader(fileout);
			under_sf->write(fileout);
			std::vector<CurveLoop> bd_loops = bd_sf->allBoundaryLoops();
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
		    }
		}
	    }
	    else
	    {
		bd_sf->writeStandardHeader(fileout);
		bd_sf->write(fileout);
	    }
	}
	// else
	// {
	//     geom_obj->writeStandardHeader(fileout);
	//     geom_obj->write(fileout);
	// }

	++cntr;
    }

    std::cout << "num_bd_sfs: " << num_bd_sfs << ", num_bd_sfs_fixed: " << num_bd_sfs_fixed <<
	", num_bd_sfs_fix_failed: " << num_bd_sfs_fix_failed << std::endl;
}
