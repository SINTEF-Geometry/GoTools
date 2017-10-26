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

        const string datadir_priv = "../../../gotools-data/step_reader/data3/"; // Relative to build/gotools-extra/step_reader
        
        infiles.push_back(datadir_priv + "Cavity_AM_obj_985.g2");

#endif

        const string datadir_public = "data/"; // Relative to build/gotools-core

#if 0
        // This test case consists of a sphere with a trim curve crossing a pole. This scenario can be
        // handled, but not by BoundedUtils::fixInvalidBoundedSurface().
        infiles.push_back(datadir_public + "test_bounded_sf_3.g2");

#endif

        GoTools::init();
    }

public:
    ObjectHeader header;
    vector<string> infiles;
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

        ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(in);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        if (!is_valid)
        {
            Go::BoundedUtils::createMissingParCvs(*bs);
            is_valid = bs->isValid(valid_state);
            if (!is_valid)
            {   // This function can mess things up, only called if there is something wrong with the
                // input (according to our code).  Preferrably we should not get here.
                const double eps_geo = bs->getEpsGeo();
                const double max_tol_mult = std::min(1.0e03, (1.0/eps_geo)); // Not allowing a tolerance larger than 1.0.
                Go::BoundedUtils::fixInvalidBoundedSurface(bs, max_tol_mult);
                is_valid = bs->isValid(valid_state);
            }
        }

        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface " << i << ", valid state: " << valid_state);

        bounded_surfaces.push_back(bs);
    }

}

