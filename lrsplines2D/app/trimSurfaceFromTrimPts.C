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


#include "GoTools/lrsplines2D/TrimCrvUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;


int main(int argc, char* argv[])
{
    if (argc != 4)
    {
	cout << "Usage: " << argv[0] << " sf.g2 trim_pts_3d.g2 bd_sf.g2" << endl;
	return -1;
    }

    
    std::ifstream filein_sf(argv[1]);
    std::ifstream filein_trim_pts(argv[2]);
    std::ofstream fileout_bd_sf(argv[3]);

    GoTools go_tools;
    go_tools.init();
    Registrator<LRSplineSurface> r293;

    ObjectHeader header;
    header.read(filein_sf);
    shared_ptr<ParamSurface> sf;
    if (header.classType() == Class_SplineSurface)
    {
	sf = shared_ptr<SplineSurface>(new SplineSurface());
    }
    else if (header.classType() == Class_LRSplineSurface)
    {
	sf = shared_ptr<LRSplineSurface>(new LRSplineSurface());
    }
    else
    {
	MESSAGE("Input surface file should be a SplineSurface or a LRSplineSurface!");
	return -1;
    }
    sf->read(filein_sf);
    if ((sf->dimension() == 1) && (sf->instanceType() == Class_LRSplineSurface))
    {
	cout << "Lifting input surface from 1D to 3D." << endl;
	shared_ptr<LRSplineSurface> lr_spline_sf = dynamic_pointer_cast<LRSplineSurface>(sf);
	lr_spline_sf->to3D();
    }

#ifndef NDEBUG
    {
	MESSAGE("Creating bounded box for input sf.");
	BoundingBox bd_box = sf->boundingBox();
	BoundingBox bd_box2 = sf->boundingBox();
	MESSAGE("Done creating bounded box #2 for input sf.");
    }
#endif

    double scale_z_factor = 5.0;
    if (scale_z_factor != 1.0)
    {
	cout << "\nRescaling the z dimension of the surface by a factor of " << scale_z_factor << "!\n" << endl;
	TrimCrvUtils::scaleZ(*sf, scale_z_factor);
    }

    const double epsgeo = 20;//10;//5;//1e-04;
    const double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
    Point translate_vec;
    vector<double> pts_2d = 
      TrimCrvUtils::readTrimPoints(filein_trim_pts, translate_vec);
    vector<vector<double> > split_pts_2d = 
      TrimCrvUtils::splitTrimPoints(pts_2d, epsgeo, kink_tol);

    TrimCrvUtils::translateSurfaceDomain(sf.get(), translate_vec);
    TrimCrvUtils::translateObject(*sf, translate_vec);

    // We must also translate the domain as it is used by the input pts.


    const int par_dim = 2;
    const int max_iter = 5;
    vector<shared_ptr<SplineCurve> > par_cvs;
    for (size_t ki = 0; ki < split_pts_2d.size(); ++ki)
    {
	shared_ptr<SplineCurve> spline_cv_appr_2d
	    (TrimCrvUtils::approximateTrimPts(split_pts_2d[ki], par_dim, epsgeo, max_iter));
	par_cvs.push_back(spline_cv_appr_2d);
    }

    // // We expect the trim curve to be a spline curve approximating some data points.
    // shared_ptr<SplineCurve> par_trim_cv(new SplineCurve());;
    // header.read(filein_par_trim_cv);
    // if (header.classType() != Class_SplineCurve)
    // {
    // 	MESSAGE("Input curve file should be a SplineCurve!");
    // 	return -1;
    // }
    // else
    // {
    // 	par_trim_cv->read(filein_par_trim_cv);
    // }

    // The curve should be CCW.
    const double int_tol = 1e-06;
    vector<shared_ptr<ParamCurve> > par_cvs2(par_cvs.begin(), par_cvs.end());
    bool loop_is_ccw = LoopUtils::loopIsCCW(par_cvs2, epsgeo, int_tol);
    if (!loop_is_ccw)
    {
	MESSAGE("We should change direction of the loop cv!");
	for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	{
	    par_cvs[ki]->reverseParameterDirection();
	}
	reverse(par_cvs.begin(), par_cvs.end());
    }

    TrimCrvUtils::moveCurveCoefsInsideSurfDomain(sf.get(), par_cvs);

    bool use_linear_segments = true;
    vector<shared_ptr<ParamCurve> > par_loop;
    if (use_linear_segments)
    {
	for (size_t ki = 0; ki < par_cvs.size(); ++ki)
	{
	    vector<shared_ptr<Line> > line_segments = TrimCrvUtils::approximateCurve(*par_cvs[ki], epsgeo);
	    par_loop.insert(par_loop.end(), line_segments.begin(), line_segments.end());
	}
    }
    else
    {
	par_loop.insert(par_loop.end(), par_cvs.begin(), par_cvs.end());
    }

    vector<shared_ptr<CurveOnSurface> > loop;
    for (size_t ki = 0; ki < par_loop.size(); ++ki)
    {
	shared_ptr<CurveOnSurface> cv_on_sf(new CurveOnSurface(sf, par_loop[ki], true));
	loop.push_back(cv_on_sf);
    }
    const bool fix_trim_cvs = false;
    const double epsgeo_bd_sf = 1e-03;
    BoundedSurface bd_sf(sf, loop, epsgeo_bd_sf, fix_trim_cvs);
    int valid_state = 0;
    bool is_valid = bd_sf.isValid(valid_state);
    if (!is_valid)
    {
	MESSAGE("Created invalid BoundedSurface, valid_state = " << valid_state);
    }
    else
    {
	MESSAGE("Surface is valid!");
    }

    bd_sf.writeStandardHeader(fileout_bd_sf);
    bd_sf.write(fileout_bd_sf);

#ifndef NDEBUG
    // Suspecting BoundingBox fails ...
    try
    {
	MESSAGE("Creating bounding box!");
	BoundingBox bd_box = bd_sf.boundingBox();
	Point low = bd_box.low();
	Point high = bd_box.high();
	MESSAGE("Done creating bounding box!");
	cout << "low: " << low[0] << ", " << low[1] << ", " << low[2] << endl;
	cout << "high: " << high[0] << ", " << high[1] << ", " << high[2] << endl;
	MESSAGE("Creating bounding box 2!");
	BoundingBox bd_box2 = bd_sf.boundingBox();
	MESSAGE("Done creating bounding box 2!");
    }
    catch (...)
    {
	MESSAGE("Failed creating bounding box!");
    }
    double debug_error = 0.0;
#endif

}
