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

#define BOOST_TEST_MODULE gotools-core/testBoundedSurface
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
        datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back("test_bounded_sf_2.g2"); // 0
        numloops.push_back(1);

        infiles.push_back("test_bounded_sf_3.g2"); // 1
        numloops.push_back(1);

        GoTools::init();
    }

public:
    ObjectHeader header;
    string datadir;
    vector<string> infiles;
    vector<int> numloops;
};


BOOST_FIXTURE_TEST_CASE(testBoundedSurface, Config)
{
    vector<shared_ptr<BoundedSurface> > bounded_surfaces;

    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {

        //string filename = "test_bounded_sf_2.g2";
        //string infile = datadir + filename;
        string infile = datadir + infiles[i];
        //string infile = "../step_reader/DemEx6woExtBlends.g2";

        ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<BoundedSurface> bs(new BoundedSurface());
        bs->read(in);

        vector<CurveLoop> loops = bs->allBoundaryLoops();
        int nloops = loops.size();
        BOOST_CHECK_EQUAL(nloops, 1);

        int valid_state = 0;
        bool is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface " << i 
            << " valid state: " << valid_state);

        Go::BoundedUtils::fixInvalidBoundedSurface(bs);
        is_valid = bs->isValid(valid_state);
        BOOST_CHECK_MESSAGE(is_valid, "BoundedSurface valid state after fixing: " 
            << valid_state);

        bounded_surfaces.push_back(bs);
    }


}

