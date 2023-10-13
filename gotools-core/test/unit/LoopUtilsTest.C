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


#define BOOST_TEST_MODULE gotools-core/LoopUtilsTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"


using namespace Go;
using std::vector;
using std::string;
using std::ifstream;


struct Config {
public:
    Config()
    {

#ifdef GOTOOLS_TEST_PRIVATE_DATA

        // Path relative to build/gotools-extra/step_reader
        const string datadir_priv = "../../gotools-private-data/step_reader/data3/Ford/";
        
        infiles.push_back(datadir_priv + "Ford_Hood_Outer_001_sf_303.g2");
        infiles.push_back(datadir_priv + "Ford_Hood_Hinge_Reinf_001_sf_7.g2");
        infiles.push_back(datadir_priv + "Ford_Hood_Outer_001_sf_3.g2");
#endif

        datadir = "data/"; // Relative to build/gotools-core

        //infiles.push_back(datadir + "bd_plane_many_holes.g2");
        infiles.push_back(datadir + "trimmed_sphere_deg_seg.g2");
        infiles.push_back(datadir + "trimmed_sphere_no_deg_seg.g2");

        GoTools::init();
    }

public:
    ObjectHeader header;
    string datadir;
    vector<string> infiles;
    vector<int> numobjects;

};


BOOST_FIXTURE_TEST_CASE(LoopUtilsTest, Config)
{
    for (auto infile : infiles)
    {
        ifstream in1(infile.c_str());
        BOOST_CHECK_MESSAGE(in1.good(), "Input file not found or file corrupt");

        shared_ptr<BoundedSurface> bd_sf(new BoundedSurface());
        header.read(in1);
        bd_sf->read(in1);

        int valid_state = -1;
        bool valid = bd_sf->isValid(valid_state);

        BOOST_CHECK_EQUAL(valid, true);

        // The input surface has 1 outer loop only.
        vector<CurveLoop> bd_loops = bd_sf->allBoundaryLoops();
        int cntr = 0;
        for (auto bd_loop : bd_loops)
        {
            // The loop is ccw.
            const double int_tol = 1.0e-03;
            bool ccw = LoopUtils::loopIsCCW(bd_loop, int_tol);
            std::cout << "LoopUtilsTest.C: ccw: " << ccw << std::endl;
            bool ccw_true = (cntr == 0);
            BOOST_CHECK_EQUAL(ccw, ccw_true);
            ++cntr;
        }
    }
}
