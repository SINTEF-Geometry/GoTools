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

#define BOOST_TEST_MODULE gotools-core/testHermiteApprEvalSurf
#include <boost/test/included/unit_test.hpp>


#include <fstream>
#include "GoTools/creators/EvalOffsetSurface.h"
#include "GoTools/creators/HermiteApprEvalSurf.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"


using namespace Go;
using std::vector;
using std::ifstream;


BOOST_AUTO_TEST_CASE(testHermiteApprEvalSurf)
{
    // @@sbr201701 Create vector of filenames to be used for the testing.

//    ifstream infile1("data/square.g2"); // Trivial, unit square for z = 0.0, bilinear.
    ifstream infile1("data/spline_surface_1.g2"); // Bicubic 18x17, almost flat.
//    ifstream infile1("data/test_bezier.g2"); // Tricky case with self-intersections for offset dist of 1.23 (and smaller).
    BOOST_CHECK_MESSAGE(infile1.good(), "Input file not found or file corrupted!");
    
    ObjectHeader header;
    shared_ptr<SplineSurface> sf(new SplineSurface());;
    infile1 >> header;
    sf->read(infile1);

    const double offset = 0.1;//1.23;//0.0;//1.23;
    const double epsgeo = 1e-06;//3;//1e-06;
    EvalOffsetSurface eval_offset_sf(sf, offset, epsgeo);

    // Creating the initial grid.
    HermiteApprEvalSurf appr_eval_sf(&eval_offset_sf, epsgeo*0.1, epsgeo*0.1);
    // Refining until within the required tolerance.
    appr_eval_sf.refineApproximation();
    // Creating the surface from the Bezier patches.
    shared_ptr<SplineSurface> offset_sf = appr_eval_sf.getSurface();
    BOOST_CHECK(offset_sf.get());
//    assert(offset_sf.get() != NULL);

    const int num_samples = 12;
    const double umin = sf->startparam_u();
    const double vmin = sf->startparam_v();
    const double umax = sf->endparam_u();
    const double vmax = sf->endparam_v();
    const double ustep = (umax - umin)/(double)(num_samples - 1);
    const double vstep = (vmax - vmin)/(double)(num_samples - 1);
    double max_error = -1.0;
    for (size_t kj = 0; kj < num_samples; ++kj)
    {
        std::cout << "kj = " << kj << std::endl;
        double vpar = vmin + (double)kj*vstep;
        for (size_t ki = 0; ki < num_samples; ++ki)
        {
            double upar = umin + (double)ki*ustep;
            Point base_pt = sf->ParamSurface::point(upar, vpar);
            Point offset_pt = offset_sf->ParamSurface::point(upar, vpar);
            double dist = base_pt.dist(offset_pt);
            double error = fabs(dist - offset);
            std::cout << "error: " << error << std::endl;
            if (error > max_error)
            {
                max_error = error;
            }
            const bool disable_test = false;//true;
            if (!disable_test)
            {
                BOOST_CHECK_CLOSE(offset, dist, epsgeo);
            }
        }
    }
    std::cout << "max_error: " << max_error << std::endl;

    std::ofstream fileout("tmp/testHermiteApprEvalSurf_result.g2");
    offset_sf->writeStandardHeader(fileout);
    offset_sf->write(fileout);

    std::ofstream fileout2("tmp/testHermiteApprEvalSurf_input.g2");
    sf->writeStandardHeader(fileout2);
    sf->write(fileout2);

}
