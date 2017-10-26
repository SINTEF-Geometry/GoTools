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

#define BOOST_TEST_MODULE LRSplineSurfaceTest
#include <boost/test/included/unit_test.hpp>
#include <fstream>

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"


using namespace Go;
using std::vector;
using std::string;
using std::ifstream;


struct Config {
public:
    Config()
    {

        datadir = "data/"; // Relative to build/lrsplines2D

        infiles.push_back(datadir + "unit_square_cubic_lr_3d.g2");

    }

public:
    ObjectHeader header;
    string datadir;
    vector<string> infiles;
    vector<int> numobjects;

};


BOOST_FIXTURE_TEST_CASE(subSurface, Config)
{
    // Assuming all infiles are LRSplineSurface. Otherwise the fixture must be changed.
    for (auto iter = infiles.begin(); iter != infiles.end(); ++iter)
    {
	ifstream in1(iter->c_str());
        BOOST_CHECK_MESSAGE(in1.good(), "Input file not found or file corrupt");

	shared_ptr<LRSplineSurface> lr_sf(new LRSplineSurface());
	header.read(in1);
	lr_sf->read(in1);

	// lr_sf subSurface() and evaluation in corner point.
	const double fuzzy = 1e-10;
	double const umin = lr_sf->startparam_u();
	double const umax = lr_sf->endparam_u();
	double const new_umin = 0.9*umin + 0.1*umax;
	double const vmin = lr_sf->startparam_v();
	double const vmax = lr_sf->endparam_v();
	const Point orig_pt = lr_sf->ParamSurface::point(new_umin, vmin);
	std::cout << "orig_pt: " << orig_pt[0] << " " << orig_pt[1] << " " << orig_pt[2] << std::endl;
	// MESSAGE("Calling subSurface().");
        std::cout << "new_umin: " << new_umin << ", vmin: " << vmin << ", umax: " << umax << ", vmax: " << vmax << std::endl;
	shared_ptr<LRSplineSurface> sub_sf(lr_sf->subSurface(new_umin, vmin, umax, vmax, fuzzy));
        std::cout << "Done calling subSurface()." << std::endl;
	// To provoke memory overwrite we call the function another time.
//	shared_ptr<LRSplineSurface> sub_sf2(lr_sf->subSurface(new_umin, vmin, umax, vmax, fuzzy));
	// MESSAGE("Done calling subSurface().");
	const Point sub_pt = sub_sf->ParamSurface::point(new_umin, vmin);
	std::cout << "sub_pt: " << sub_pt[0] << " " << sub_pt[1] << " " << sub_pt[2] << std::endl;
	const double dist = orig_pt.dist(sub_pt);
	std::cout << "dist: " << dist << std::endl;
	const double tol = 1e-14;
	BOOST_CHECK_LT(dist, tol);
    }
}
