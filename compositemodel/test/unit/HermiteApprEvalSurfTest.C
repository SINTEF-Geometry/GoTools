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

#define BOOST_TEST_MODULE gotools-core/HermiteApprEvalSurfTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using std::vector;
using std::ifstream;
using std::string;

struct Config {
public:
    Config()
    {

        const std::string datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back(datadir + "bounded_surface.g2");
        numloops.push_back(1);

        infiles.push_back(datadir + "test_bounded_sf_2.g2");
        numloops.push_back(1);

        infiles.push_back(datadir + "test_bounded_sf_3.g2");
        numloops.push_back(1);

    }

public:
    ObjectHeader header;
    vector<string> infiles;
    vector<int> numloops;
};


BOOST_AUTO_TEST_CASE(HermiteApprEvalSurfTest)
{

    vector<string> filenames;
    vector<double> offset, epsgeo;

    const std::string data_basedir("../../gotools-data/compositemodel/Offset");

#if 1
    filenames.push_back(data_basedir+"/TopSolid/TopSolid_bump.g2");
    offset.push_back(-3.0e-03);
    epsgeo.push_back(1.0e-05);
#endif

#if 1
    // Self-intersections for offset dist of appr 0.3 and larger.
    filenames.push_back(data_basedir+"/test_bezier.g2");
    offset.push_back(0.1);
    epsgeo.push_back(1.0e-03);
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


        shared_ptr<SplineSurface> offset_sf;
        OffsetSurfaceStatus status = OffsetSurfaceUtils::offsetSurfaceSet(sfs, offset[ki], epsgeo[ki], offset_sf);

        if (offset_sf.get() == NULL)
        {
            BOOST_ERROR("Offset surface was not created.");
            continue;
        }

        if (status == OFFSET_OK)
        {
            MESSAGE("Success!");
        }
        else
        {
            MESSAGE("Failure! Returned status (!= 0): " << status);
            continue;
        }

        // Since we may need to handle self intersections and kinks in the input surface set we are not
        // guaranteed that we are within the offset-dist +/- epsgeo for all points. Assuming that the
        // method is capable of reporting on the status.
        BOOST_CHECK_EQUAL(status, OFFSET_OK);
        if (status != 0)
        {
            continue;
        }
    }
}
