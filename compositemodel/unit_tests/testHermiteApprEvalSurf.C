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
#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::vector;
using std::ifstream;
using std::string;

BOOST_AUTO_TEST_CASE(testHermiteApprEvalSurf)
{

    vector<string> filenames;
    vector<double> offset, epsgeo;
#if 0

    filenames.push_back("data/square.g2"); // Trivial, unit square for z = 0.0, bilinear.

    filenames.push_back("data/spline_surface_1.g2"); // Bicubic 18x17, almost flat.

    filenames.push_back("data/Offset/fanta_ro2.g2");

    filenames.push_back("data/TopSolid/sfw1.g2");

    filenames.push_back("data/TopSolid/sfw2.g2");

    filenames.push_back("data/TopSolid/sfw1_sfw2.g2");
    
    filenames.push_back("data/Offset/fanta_ro2_sub2.g2");

    // Two adjacent unit planes.
    filenames.push_back("data/TopSolid/TopSolid_Surf__20170313-174324.189_221.g2");
#endif

#if 0
    // Two orthogonal planes joined by a cylinder segment (w/ radius of curvature -1.38843).
    // With epsgeo 1e-04 there is too much data to handle removal of self intersections.
    filenames.push_back("data/Offset/yta4.g2");
    offset.push_back(-1.5);//1.3);//(0.01);//0.1;//1.23; //0.2;
    epsgeo.push_back(1.0e-04);//3);//6;
#endif

#if 0
    // Tricky case with self-intersections for offset dist of appr 0.3 and larger. Surface set contains
    // a degenerate spline surface with the degenerate in the middle of the surface set edge. Results in
    // a bad offset boundary curve. Fix!
    filenames.push_back("data/Offset/fanta_ro2_sub.g2");
    offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
    epsgeo.push_back(1.0e-04);//3);//6;
#endif

#if 1
    // Self-intersections for offset dist of appr 0.3 and larger.
    filenames.push_back("data/test_bezier.g2");
    offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
    epsgeo.push_back(1.0e-04);//3);//6;
#endif
    
#if 0
    // Self-int: offset=1e-02,eps=1e-03.
    // Success with removing self intersections using smoothing: offset=1e-03,eps=1e-04.
    filenames.push_back("data/TopSolid/TopSolid_Surf__20170404-123801.816.g2");
    offset.push_back(0.001);//(0.01);//0.1;//1.23; //0.2;
    epsgeo.push_back(1.0e-04);//3);//6;
#endif

#if 0
    // Degenerate patch (triangle): Ok w/ offset=1e-02,eps=1e-03.
    filenames.push_back("data/Offset/fanta_ro2_sub2b.g2");
    offset.push_back(0.01);//0.1;//1.23; //0.2;
    epsgeo.push_back(1.0e-03);//6;
#endif
    
    for (size_t ki = 0; ki < filenames.size(); ++ki)
    {
        std::cout << "\nTesting offsetSurfaceSet() for file " << filenames[ki] <<
            ", offset: " << offset[ki] << ", epsgeo: " << epsgeo[ki] << std::endl;
        
        ifstream infile(filenames[ki]);
        if (!(infile.good()))
        {
            BOOST_ERROR("Input file not found or file corrupted!");
            continue;
        }

        vector<shared_ptr<ParamSurface> > sfs;
        while (!infile.eof())
        {
            ObjectHeader header;
            infile >> header;
            if (header.classType() != Class_SplineSurface)
            {
                BOOST_ERROR("Input was not a SplineSurface, not yet supported!");
                continue;
            }

            shared_ptr<ParamSurface> sf(new SplineSurface());
            sf->read(infile);

            sfs.push_back(sf);

            Utils::eatwhite(infile);
        }        

        std::string input_filename("tmp/testHermiteApprEvalSurf_input.g2");
        std::ofstream fileout2(input_filename);
        for (size_t kk = 0; kk < sfs.size(); ++kk)
        {
            sfs[kk]->writeStandardHeader(fileout2);
            sfs[kk]->write(fileout2);
        }

        shared_ptr<SplineSurface> offset_sf;
        OffsetSurfaceStatus status = OffsetSurfaceUtils::offsetSurfaceSet(sfs, offset[ki], epsgeo[ki], offset_sf);

        if (offset_sf.get() == NULL)
        {
            BOOST_ERROR("Offset surface was not created.");
            continue;
        }

        std::string output_filename("tmp/testHermiteApprEvalSurf_result.g2");
        std::ofstream fileout(output_filename);
        offset_sf->writeStandardHeader(fileout);
        offset_sf->write(fileout);

        if (status == OFFSET_OK)
        {
            MESSAGE("Success!");
        }
        else
        {
            MESSAGE("Failure! Returned status (!= 0): " << status);
            continue;
        }

        std::cout << "Input and result written to " << input_filename << " and " << output_filename << std::endl;

        BOOST_CHECK_EQUAL(status, 0);
        if (status != 0)
        {
            continue;
        }
        
        // // We currently support testing for 1 input surface only.
        // BOOST_CHECK_EQUAL(sfs.size(), 1);

        double global_max_error = 0.0;
        double global_max_u = -1.0, global_max_v = -1.0;
        double global_max_clo_u, global_max_clo_v;
        for (size_t kk = 0; kk < sfs.size(); ++kk)
        {
            const RectDomain& rect_dom = sfs[kk]->containingDomain();
            const int num_samples = 73;
            const double umin = rect_dom.umin();
            const double vmin = rect_dom.vmin();
            const double umax = rect_dom.umax();
            const double vmax = rect_dom.vmax();
            const double ustep = (umax - umin)/(double)(num_samples - 1);
            const double vstep = (vmax - vmin)/(double)(num_samples - 1);
            double max_error = -1.0;
            double max_u, max_v;
            double max_clo_u, max_clo_v;
            for (size_t kj = 0; kj < num_samples; ++kj)
            {
                double vpar = vmin + (double)kj*vstep;
                for (size_t kh = 0; kh < num_samples; ++kh)
                {
                    double upar = umin + (double)kh*ustep;
                    Point base_pt = sfs[kk]->point(upar, vpar);
//                    Point offset_pt = offset_sf->ParamSurface::point(upar, vpar);
                    double clo_u, clo_v, clo_dist;
                    Point offset_pt;
                    offset_sf->closestPoint(base_pt,
                                            clo_u, clo_v, offset_pt, clo_dist,
                                            epsgeo[ki]);
                    // @@sbr201704 We do not check the direction, assuming that this is handled correctly
                    // by the method.
                    double dist = base_pt.dist(offset_pt);
                    double error = fabs(dist - offset[ki]);
                    if (error > max_error)
                    {
                        max_error = error;
                        max_u = upar;
                        max_v = vpar;
                        max_clo_u = clo_u;
                        max_clo_v = clo_v;
                    }
                    const bool disable_test = false;
                    if (!disable_test)
                    {
                        BOOST_CHECK_LE(error, epsgeo[ki]); // Checking if error <= epsgeo.
                    }
                }
            }
            if (max_error > global_max_error)
            {
                global_max_error = max_error;
                global_max_u = max_u;
                global_max_v = max_v;
                global_max_clo_u = max_clo_u;
                global_max_clo_v = max_clo_v;
            }
            std::cout << "max_error: " << max_error << ", max_u: " << max_u << ", max_v: " << max_v << std::endl;
        }
        std::cout << "global_max_error: " << global_max_error << ", epsgeo: " << epsgeo[ki] << ", ki: " << ki <<
            ", global_max_u: " << global_max_u << ", global_max_v: " << global_max_v <<
            ", global_max_clo_u: " << global_max_clo_u << ", global_max_clo_v: " << global_max_clo_v << std::endl;
    }
}
