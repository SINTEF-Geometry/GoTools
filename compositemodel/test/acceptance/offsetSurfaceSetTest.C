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

#define BOOST_TEST_MODULE compositemodel/offsetSurfaceSetTest
#include <boost/test/included/unit_test.hpp>

#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"

#include <fstream>
#include <math.h>
#include <vector>
#include <string>

using namespace Go;
using std::vector;
using std::string;
using std::ifstream;

// #define TEST_NEW_CASES

// #define TEST_UNSUPPORTED_CASES

// #define TEST_FAILING_CASES

#define TEST_WORKING_CASES

struct Config {
public:
    Config()
        : input_filename("tmp/offsetSurfaceSetTest_input.g2"),
          output_filename("tmp/offsetSurfaceSetTest_result.g2")

    {

        const std::string data_basedir("../../gotools-private-data/compositemodel/Offset");

#ifdef GOTOOLS_TEST_PRIVATE_DATA

        // NEW CASES!!!
#ifdef TEST_NEW_CASES

#endif

        // CASES NOT SUPPORTED (EARLY EXIT WITH ERROR MESSAGE)!!!        
#ifdef TEST_UNSUPPORTED_CASES

        // Degenerate patch (triangle): Ok w/ offset=1e-02,eps=1e-03. Using 0.01*epsgeo as curvature_tol.
        // Contains kink curve which is not an iso curve. Exits with failure.
        filenames.push_back(data_basedir+"/fanta_ro2_sub2b.g2");
        offset.push_back(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-03);//6;

        // Tricky case with a degenerate surface. The degenerate point is not well defined with the
        // normal varying as the evaluation approaches from different iso lines. Will require extensive
        // smoothing. Currently not handled using SplineSurface (equation system too large). We need
        // LRSplineSurface (and even with that the task may be to tricky). Exits with failure.
        filenames.push_back(data_basedir+"/fanta_ro2_sub.g2");
        offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;
#endif

        // FAILING CASES!!!
#ifdef TEST_FAILING_CASES
        // Added 2017-10-09. Crash 3. Tricky case with a degenerate underlying surface, with a trimmed
        // boundary meeting in a degenerate corner (parallell tangents).
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_BoundedSurf__20170929-115934.549.g2");
        offset.push_back(1.0e-03);
        epsgeo.push_back(1.0e-04);
#endif

        // WORKING CASES!!!
#ifdef TEST_WORKING_CASES
        
        // Added 2017-10-09. Crash 2.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf_OffsetTest_20170929-110347.373_406_434.g2");
        offset.push_back(1.0e-03);
        epsgeo.push_back(1.0e-04); // Initially given epsgeo = 1.0e-03, conflicts with topology (too many twin edges).

        // Added 2017-10-09. Crash 1.
        filenames.push_back(data_basedir+"/TopSolid/Source-TopSolid_SplineSurf_OffsetTest_20170929-103855.922_949.g2");
        offset.push_back(1.0e-03);
        epsgeo.push_back(1.0e-04); // Initially given epsgeo = 1.0e-03, conflicts with topology (too many twin edges).

        filenames.push_back(data_basedir+"/TopSolid/TopSolid_BoundedSurf__20170623-173106.658.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);

        // Added 2017-06-23
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf__20170623-173106.544.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);

        // Added 2017-06-23
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf__20170623-173916.090.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);

        // Added 2017-06-23
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf__20170623-175900.205.g2");
//      filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf__20170627-130342.267.g2"); The same set.
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);

        // 2 bilinear quadrats with z = 0.018.
        // Tricky case which required us to to define a transition zone.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_Surf__20170313-174324.189_221.g2");
        offset.push_back(0.1);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;

        // 4 planes meeting in kinks along linear segments.
        // Tricky case which required us to to define a transition zone.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_SplineSurf__20170504-1332_kinks.g2");
        offset.push_back(3.0e-03);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;

        // A single surface with a bump in the interior. Use a negative offset value to create a self
        // intersection in the bump. Easy when using a positive offset dist, tricky when using a negative
        // offset dist. With a large offset value (< -3.0e-03) the self intersection is global.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_bump.g2");
        offset.push_back(-3.0e-03);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;

        // Self-intersections for offset dist of appr 0.3 and larger.
        filenames.push_back(data_basedir+"/test_bezier.g2");
        offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;

        // Self-int: offset=1e-02,eps=1e-03.
        // Success with removing self intersections using smoothing: offset=1e-03,eps=1e-04.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_Surf__20170404-123801.816.g2");
        offset.push_back(0.001);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;

        // Two orthogonal planes joined by a trimmed cylinder segment (w/ radius of curvature -1.38843).
        // 201706: Success after setting the twist vector to zero and reducing precond omega from 0.1 to 0.01.
        filenames.push_back(data_basedir+"/yta4.g2");
        offset.push_back(-1.5);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;

        // Illegal bounded surface, the trim loop is CW! The offset works but the offset surface is flipped.
        filenames.push_back(data_basedir+"/TopSolid/TopSolid_BoundedSurf__20170623-173916.162.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

#endif

        GoTools::init();

    }

public:
    vector<string> filenames;
    vector<double> offset;
    vector<double> epsgeo;
    ObjectHeader header;
    const std::string input_filename;
    const std::string output_filename;

};


BOOST_FIXTURE_TEST_CASE(offsetSurfaceSet, Config)
{
    int num_success = 0;
    int num_failures = 0;
    int num_exceptions = 0;    
    for (size_t ki = 0; ki < filenames.size(); ++ki)
    {
        // Read input arguments
        std::ifstream infile(filenames[ki].c_str());
        BOOST_CHECK_MESSAGE(infile.good(), "Input file not found or file corrupt");

        vector<shared_ptr<ParamSurface> > sfs;
        while (!infile.eof())
        {
            infile >> header;
            shared_ptr<ParamSurface> sf;
            if (header.classType() == Class_SplineSurface)
            {
                sf = shared_ptr<ParamSurface>(new SplineSurface());
            }
            else if (header.classType() == Class_BoundedSurface)
            {
                sf = shared_ptr<ParamSurface>(new BoundedSurface());
            }
            else
            {
                BOOST_ERROR("Input surface type not yet supported: " << header.classType());
                continue;
            }

            sf->read(infile);

            sfs.push_back(sf);

            Utils::eatwhite(infile);
        }        

        std::ofstream fileout2(input_filename);
        for (size_t kk = 0; kk < sfs.size(); ++kk)
        {
            sfs[kk]->writeStandardHeader(fileout2);
            sfs[kk]->write(fileout2);
        }

        shared_ptr<SplineSurface> offset_sf;
        OffsetSurfaceStatus status;
        try
        {
            status = OffsetSurfaceUtils::offsetSurfaceSet(sfs, offset[ki], epsgeo[ki], offset_sf);
        }
        catch (...)
        {
            std::cout << "offsetSurfaceSet(): Caught exception for filename " << filenames[ki] << std::endl;
            BOOST_ERROR("Exception caught");
            ++num_exceptions;
            continue;
        }

        BOOST_CHECK_MESSAGE(status == OFFSET_OK,
                            "offsetSurfaceSet(): Failure! Returned status (!= 0): " << status);
        if (status == OFFSET_OK)
        {
            ++num_success;
            if (offset_sf.get())
            {
                std::ofstream fileout(output_filename);
                offset_sf->writeStandardHeader(fileout);
                offset_sf->write(fileout);
            }
        }
        else
        {
            ++num_failures;
        }
    }

    BOOST_TEST_MESSAGE("\nnum files: " << filenames.size() << ", num success: " << num_success << ", num failures: " <<
                  num_failures << ", num_exceptions: " << num_exceptions);

}
