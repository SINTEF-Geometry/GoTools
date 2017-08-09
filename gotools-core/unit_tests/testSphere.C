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

#define BOOST_TEST_MODULE gotools-core/testSphere
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/Sphere.h"


using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {
        datadir = "data/"; // Relative to build/gotools-core

        //GoTools::init();
    }

public:
    string datadir;
    vector<string> infiles;
};


BOOST_FIXTURE_TEST_CASE(testSphere, Config)
{
    // A sphere
    double radius = 1.0;
    Point location(0.0, 0.0, 0.0);
    Point z_axis(0.0, 0.0, 1.0);
    Point x_axis(1.0, 0.0, 0.0);
    bool isSwapped = true;
    Sphere sphere(radius, location, z_axis, x_axis, isSwapped);

    // closestPoint()
    double clo_u, clo_v, clo_dist;
    Point clo_pt;
    double epsilon = 1.0e-6;
    Point pt_top(0.0, 0.0, 1.0);
    sphere.closestPoint(pt_top, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, 0.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);

    Point pt_bottom(0.0, 0.0, -1.0);
    sphere.closestPoint(pt_bottom, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, -0.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);

    Point pt1(0.0, -1.0, 0.0);
    sphere.closestPoint(pt1, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, 0.0, epsilon);
    BOOST_CHECK_CLOSE(clo_v, 1.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);

    // A space curve
    Point normal(-1.0, 0.0, 0.0);
    Point circle_x(0.0, 0.0, 1.0);
    bool isReversed = true;
    Circle circle(radius, location, normal, circle_x, isReversed);
    double pi = 3.14159265358979; // M_PI slightly truncated
    circle.setParamBounds(pi, 2.0*M_PI);

    // getElementaryParamCurve()
    double tol = 1.0e-6;
    shared_ptr<ElementaryCurve> elem_cv 
        = sphere.getElementaryParamCurve(&circle, tol);
    BOOST_CHECK_MESSAGE(elem_cv, "Didn't get elementary curve.");

    normal = Point(0.0, 0.0, 1.0);
    circle_x = Point(0.0, -1.0, 0.0);
    Circle circ2(radius, location, normal, circle_x, isReversed);
    circ2.setParamBounds(-0.5*M_PI, 0.5*M_PI);
    sphere.swapParameterDirection();
    elem_cv = sphere.getElementaryParamCurve(&circ2, tol);
    double t0 = elem_cv->startparam();
    Point pt(2);
    elem_cv->point(pt, t0);
    BOOST_CHECK_CLOSE(pt[0], 2.0*M_PI, tol);

}

