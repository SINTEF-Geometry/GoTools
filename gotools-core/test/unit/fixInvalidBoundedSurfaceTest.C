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

#define BOOST_TEST_MODULE gotools-core/fixInvalidBoundedSurfaceTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedUtils.h"

using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {

#ifdef GOTOOLS_TEST_PRIVATE_DATA

        // Path relative to build/gotools-extra/step_reader
        const string datadir_priv = "../../gotools-private-data/step_reader/data3/";

        // Ford models.
        infiles.push_back(datadir_priv + "Ford/Ford_Car_Hood_inner_001_obj_378.g2");
        valid_model.push_back(true);
        infiles.push_back(datadir_priv + "Ford/Ford_Car_Hood_inner_001_obj_1194.g2");
        valid_model.push_back(true);

        // CAxMan models.
        infiles.push_back(datadir_priv + "CaxMan/Mould_Final_Version_1/stock_cavity_stage4_model_2_obj_49_mod.g2");
        valid_model.push_back(false); // The example is invalid (the surface seam must be rotated).

        infiles.push_back(datadir_priv +
                          "CaxMan/Mould_Final_Version_1/Stock_Cavity_AM_stage4_R8_8_Printing_colorised_model_2_obj_12233.g2");
	//"CaxMan/Mould_Final_Version_1/12233.g2");
        valid_model.push_back(true);
#endif

        const string datadir_public = "data/"; // Relative to build/gotools-core

#if 0
        // This test case consists of a sphere with a trim curve crossing a pole. This scenario can be
        // handled, but not by BoundedUtils::fixInvalidBoundedSurface().
        infiles.push_back(datadir_public + "test_bounded_sf_3.g2");
        valid_model.push_back(true);

#endif

        GoTools::init();
    }

public:
    ObjectHeader header;
    vector<string> infiles;
    vector<bool> valid_model;
};


BOOST_FIXTURE_TEST_CASE(BoundedSurfaceTest, Config)
{
    vector<shared_ptr<BoundedSurface> > bounded_surfaces;

    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {
        std::cout << "i: " << i << ", filename: " << infiles[i] << std::endl;
        //string filename = "test_bounded_sf_2.g2";
        //string infile = datadir + filename;
        string infile = infiles[i];
        //string infile = "../step_reader/DemEx6woExtBlends.g2";

	// std::cout << "infile: " << infile << std::endl;
        ifstream in(infile.c_str(), std::ifstream::in);
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(in);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        std::cout << "Status for input: is_valid: " << is_valid << ", valid_state: " << valid_state << std::endl;
        if (!is_valid)
        {
#ifndef NDEBUG
            std::cout << "DEBUG: Surface is not valid!" << std::endl;
#endif
            bool success = Go::BoundedUtils::createMissingParCvs(*bs);
            if (success)
            {
                bs->analyzeLoops();
                is_valid = bs->isValid(valid_state);
            }

            if (!is_valid)
            {   // This function can mess things up, only called if there is something wrong with the
                // input (according to our code).  Preferrably we should not get here.
                const double eps_geo = bs->getEpsGeo();
                // For these test cases we do not allow changing of the tolerance.
                const double max_tol_mult = 1.0;//std::min(1.0e03, (1.0/eps_geo)); // Not allowing a tolerance larger than 1.0.
#ifndef NDEBUG
                std::cout << "DEBUG: Surface is not valid after createMissingParCvs()!" << std::endl;
#endif
                Go::BoundedUtils::fixInvalidBoundedSurface(bs, max_tol_mult);
                is_valid = bs->isValid(valid_state);
            }
        }

        const bool result_valid = (is_valid == valid_model[i]);

        BOOST_CHECK_MESSAGE(result_valid, "BoundedSurface " << i << ", valid state: " << valid_state);

        bounded_surfaces.push_back(bs);
    }

}

