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


#define BOOST_TEST_MODULE gotools-core/CurveOnSurfaceTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedUtils.h"

using namespace std;
using namespace Go;


#if 0
struct Config {
public:
    Config()
    {

        const std::string datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back(datadir + "spline_cylinder.g2");
        GoTools::init();
    }

public:
    ObjectHeader header;
    vector<string> infiles;
};
#endif


BOOST_AUTO_TEST_CASE(sameOrientation)
{
    // We read a cylinder.
    const std::string infile = "data/spline_cylinder.g2"; // Relative to build/gotools-core

    ifstream in(infile.c_str());
    BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
    ObjectHeader header;
    header.read(in);
    shared_ptr<SplineSurface> cyl(new SplineSurface());
    cyl->read(in);

    // We then extract an iso-circle around the cylinder and create the corresponding
    // parameter curve, store is as a CurveOnSurface.
    const RectDomain& rect_dom = cyl->containingDomain();
    const bool pardir_is_u = true;
    const double iso_par = (cyl->isBounded()) ? rect_dom.vmin() : 0.0;
    std::vector<shared_ptr<ParamCurve> > iso_cvs = cyl->constParamCurves(iso_par, pardir_is_u);

    Point start(rect_dom.umin(), iso_par);
    Point end(rect_dom.umax(), iso_par);

    shared_ptr<Line> par_cv(new Line(start, end, rect_dom.umin(), rect_dom.umax()));

    const bool pref_par = false;
    CurveOnSurface cv_on_sf(cyl, par_cv, iso_cvs[0], pref_par);

    bool same_orientation = cv_on_sf.sameOrientation();
    BOOST_CHECK_EQUAL(same_orientation, true);

    // We then reverse the parameter curve. The sameOrientation() function should now return the opposite
    // value.
    par_cv->reverseParameterDirection();
    same_orientation = cv_on_sf.sameOrientation();
    BOOST_CHECK_EQUAL(same_orientation, false);

}

