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


#define BOOST_TEST_MODULE OffsetSurfaceTest
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/OffsetSurface.h"


using namespace Go;
using std::vector;


struct Config {
public:

    Config()
    {

        const std::string datadir = "data/"; // Relative to build/gotools-core

        infiles.push_back(datadir + "square.g2");

    }

public:
    std::vector<std::string> infiles;
    ObjectHeader header;

};


// We test the constParamCurve() function.
BOOST_FIXTURE_TEST_CASE(closestBoundaryPoint, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i)
    {
        std::string infile = infiles[i];

        std::ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<SplineSurface> spline_sf(new SplineSurface());
        spline_sf->read(in);

        const double offset_dist = 1.0;
        const double epsgeo = 1.0e-03;
        OffsetSurface offset_sf(spline_sf, offset_dist, epsgeo);

        // We compare an offset pt along the boundary with the function call.
        double wgt = 0.659;
        const double upar = wgt*spline_sf->startparam_u() + (1.0 - wgt)*spline_sf->endparam_u();
        const double vpar = spline_sf->startparam_v();
        Point bd_pt = offset_sf.ParamSurface::point(upar, vpar);
        Point clo_pt;
        double clo_dist, clo_u, clo_v;
        const double epsilon = 1.0e-08;
        offset_sf.closestBoundaryPoint(bd_pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);

        BOOST_CHECK_SMALL(clo_dist, epsgeo);
    }
}

// We test the constParamCurve() function.
BOOST_FIXTURE_TEST_CASE(constParamCurve, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i)
    {
        std::string infile = infiles[i];

        std::ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<SplineSurface> spline_sf(new SplineSurface());
        spline_sf->read(in);

        const double offset_dist = 1.0;
        const double epsgeo = 1.0e-03;
        OffsetSurface offset_sf(spline_sf, offset_dist, epsgeo);

        // We extract an iso-curve in both dirs, then test the distance.
        double wgt = 0.347;
        double iso_par_u = wgt*spline_sf->startparam_u() + (1.0 - wgt)*spline_sf->endparam_u();
        std::vector<shared_ptr<ParamCurve> > cvs_v_dir = offset_sf.constParamCurves(iso_par_u, false);
        assert(cvs_v_dir.size() == 1);
        wgt = 0.659;
        double test_par = wgt*spline_sf->startparam_v() + (1.0 - wgt)*spline_sf->endparam_v();
        Point offset_pt = offset_sf.ParamSurface::point(iso_par_u, test_par);
        Point iso_cv_pt = cvs_v_dir[0]->point(test_par);
        double dist = offset_pt.dist(iso_cv_pt);
        BOOST_CHECK_SMALL(dist, epsgeo);

        wgt = 0.758;
        double iso_par_v = wgt*spline_sf->startparam_v() + (1.0 - wgt)*spline_sf->endparam_v();
        std::vector<shared_ptr<ParamCurve> > cvs_u_dir = offset_sf.constParamCurves(iso_par_v, true);
        assert(cvs_u_dir.size() == 1);
        wgt = 0.338;
        test_par = wgt*spline_sf->startparam_u() + (1.0 - wgt)*spline_sf->endparam_u();
        offset_pt = offset_sf.ParamSurface::point(test_par, iso_par_v);
        iso_cv_pt = cvs_u_dir[0]->point(test_par);
        dist = offset_pt.dist(iso_cv_pt);
        BOOST_CHECK_SMALL(dist, epsgeo);
    }
}


// We test the constParamCurve() function.
BOOST_FIXTURE_TEST_CASE(allBoundaryLoops, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i)
    {
        std::string infile = infiles[i];

        std::ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<SplineSurface> spline_sf(new SplineSurface());
        spline_sf->read(in);

        const double offset_dist = 1.0;
        const double epsgeo = 1.0e-03;
        OffsetSurface offset_sf(spline_sf, offset_dist, epsgeo);

        std::vector<CurveLoop> bd_loops = offset_sf.allBoundaryLoops();
        assert(bd_loops.size() == 1);

        // We sample a number of points on the loop and compare with the offset_sf.
        for (auto& loop : bd_loops)
        {
            for (auto cv : loop)
            {
                const double wgt = 0.316;
                const double tpar = wgt*cv->startparam() + (1.0 - wgt)*cv->endparam();
                Point loop_pt = cv->point(tpar);
                Point clo_pt;
                double clo_dist, clo_u, clo_v;
                const double epsilon = 1.0e-08;
                offset_sf.closestBoundaryPoint(loop_pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);

                BOOST_CHECK_SMALL(clo_dist, epsgeo);
            }
        }
    }
}


// We test the constParamCurve() function.
BOOST_FIXTURE_TEST_CASE(asSplineSurface, Config)
{
    int nfiles = infiles.size();
    for (int i = 0; i < nfiles; ++i)
    {
        std::string infile = infiles[i];

        std::ifstream in(infile.c_str());
        BOOST_CHECK_MESSAGE(in.good(), "Input file not found or file corrupt");
        header.read(in);
        shared_ptr<SplineSurface> spline_sf(new SplineSurface());
        spline_sf->read(in);

        const double offset_dist = 1.0;
        const double epsgeo = 1.0e-03;
        OffsetSurface offset_sf(spline_sf, offset_dist, epsgeo);

        shared_ptr<SplineSurface> offset_spline_sf(offset_sf.asSplineSurface());

        double wgt = 0.659;
        const double upar = wgt*spline_sf->startparam_u() + (1.0 - wgt)*spline_sf->endparam_u();
        wgt = 0.328;
        const double vpar = wgt*spline_sf->startparam_v() + (1.0 - wgt)*spline_sf->endparam_v();

        Point offset_pt = offset_sf.ParamSurface::point(upar, vpar);

        Point spline_pt = offset_spline_sf->ParamSurface::point(upar, vpar);

        double dist = offset_pt.dist(spline_pt);

        std::cout << "dist: " << dist << std::endl;

        BOOST_CHECK_SMALL(dist, epsgeo);
    }
}
