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

#define BOOST_TEST_MODULE testSplineSurface
#include <boost/test/included/unit_test.hpp>

#include "GoTools/geometry/SplineSurface.h"


using namespace Go;
using std::vector;


BOOST_AUTO_TEST_CASE(testSplineSurface)
{
    int n = 1;
    int m = 4;
    //BOOST_CHECK_EQUAL(n+1, m);
    BOOST_CHECK_EQUAL(n+3, m);

    // Data from looped_surface.g2
    int dim = 3;
    int ncoefsu = 5;
    int ncoefsv = 2;
    int orderu = 4;
    int orderv = 2;
    double knotsu[] = { 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0 };
    double knotsv[] = { 0.0, 0.0, 2.0, 2.0 };
    double coefs[] = { 
        -1.0, -1.0, -1.0,
        0.5, -1.0, 0.5,
        0.0, -1.0, 2.0,
        -0.5, -1.0, 0.5,
        1.0, -1.0, -1.0,
        -1.0, 1.0, -1.0,
        0.5, 1.0, 0.5,
        0.0, 1.0, 2.0,
        -0.5, 1.0, 0.5,
        1.0, 1.0, -1.0
    };
    SplineSurface surf(ncoefsu, ncoefsv, orderu, orderv, knotsu, knotsv,
        coefs, dim);

    // Check knot vector
    vector<int> multu, multv;
    surf.basis_u().knotMultiplicities(multu);
    surf.basis_v().knotMultiplicities(multv);
    vector<double> knotvalsu, knotvalsv;
    surf.basis_u().knotsSimple(knotvalsu);
    surf.basis_v().knotsSimple(knotvalsv);
    BOOST_CHECK_EQUAL(multu.size(), 3);
    BOOST_CHECK_EQUAL(multu[0], 4);
    BOOST_CHECK_EQUAL(multu[1], 1);
    BOOST_CHECK_EQUAL(multu[2], 4);
    BOOST_CHECK_EQUAL(multv.size(), 2);
    BOOST_CHECK_EQUAL(multv[0], 2);
    BOOST_CHECK_EQUAL(multv[1], 2);
    BOOST_CHECK_EQUAL(knotvalsu.size(), 3);
    BOOST_CHECK_EQUAL(knotvalsu[0], 0.0);
    BOOST_CHECK_EQUAL(knotvalsu[1], 1.0);
    BOOST_CHECK_EQUAL(knotvalsu[2], 2.0);
    BOOST_CHECK_EQUAL(knotvalsv.size(), 2);
    BOOST_CHECK_EQUAL(knotvalsv[0], 0.0);
    BOOST_CHECK_EQUAL(knotvalsv[1], 2.0);

}
