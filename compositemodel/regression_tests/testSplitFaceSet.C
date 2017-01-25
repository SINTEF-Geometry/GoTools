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

#define BOOST_TEST_MODULE compositemodel/splitFaceSet
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"

using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {
        datadir = "data/"; // Relative to build/gotools/compositemodel

        infiles.push_back("b25.g2");
        numfaces.push_back(20);

        infiles.push_back("b26.g2");
        numfaces.push_back(16);

        gap = 0.001; // 0.001;
        neighbour = 0.01; // 0.01;
        kink = 0.01;
        approxtol = 0.01;

    }

public:
    string datadir;
    vector<string> infiles;
    vector<int> numfaces;
    double gap;
    double neighbour;
    double kink;
    double approxtol;
};


BOOST_FIXTURE_TEST_CASE(splitFaceSet, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i) {

        string infile = datadir + infiles[i];

        // Read input arguments
        std::ifstream file1(infile.c_str());
        BOOST_CHECK_MESSAGE(file1.good(), "Input file not found or file corrupt");

        std::ofstream file2("outfile.g2");

        CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

        CompositeModel *model = factory.createFromG2(file1);

        SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

        if (sfmodel)
        {
            std::vector<shared_ptr<ftSurface> > faces = sfmodel->allFaces();

            //RegularizeFaceSet reg(faces, gap, kink, true);
            RegularizeFaceSet reg(faces, gap, kink, false);
            std::vector<shared_ptr<ftSurface> > sub_faces = reg.getRegularFaces();

            int nsubfaces = sub_faces.size();
            BOOST_CHECK_EQUAL(nsubfaces, numfaces[i]);

            for (int ki=0; ki < nsubfaces; ++ki)
            {
                shared_ptr<ParamSurface> surf = sub_faces[ki]->surface();
                surf->writeStandardHeader(file2);
                surf->write(file2);
            }
        }
    }
}

