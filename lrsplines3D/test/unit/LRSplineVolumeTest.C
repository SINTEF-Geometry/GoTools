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

#define BOOST_TEST_MODULE LRSplineVolumeTest
#include <boost/test/included/unit_test.hpp>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"


using namespace Go;
using std::vector;
using std::string;
using std::ifstream;


struct Config {
public:
    Config()
    {

        datadir = "data/"; // Relative to build/lrsplines3D

        infiles.push_back(datadir + "vol.g2");

    }

public:
    ObjectHeader header;
    string datadir;
    vector<string> infiles;
    vector<int> numobjects;

};

bool my_equal_function (double i, double j)
{
    double tol = 1e-10;//14;
    return (fabs(i - j) < tol);
}


BOOST_FIXTURE_TEST_CASE(refine, Config)
{
    // Assuming all infiles are SplineVolume. Otherwise the fixture must be changed.
    for (auto iter = infiles.begin(); iter != infiles.end(); ++iter)
    {
	ifstream in1(iter->c_str());
        BOOST_CHECK_MESSAGE(in1.good(), "Input file not found or file corrupt");

	shared_ptr<SplineVolume> spline_vol(new SplineVolume());

        try
        {
            header.read(in1);

            // We expect all the files to contain 1 SplineVolume.
            BOOST_CHECK_EQUAL(header.classType(), Class_SplineVolume);

            spline_vol->read(in1);
        }
        catch (...)
        {
            BOOST_ERROR("Reading of lrspline volume failed.");
            continue;
        }

        const double knot_tol = 1.0e-06;
        LRSplineVolume lr_vol(spline_vol.get(), knot_tol);

        int order_u = spline_vol->order(0);
        int order_v = spline_vol->order(1);
        int order_w = spline_vol->order(2);

        int num_coefs_u = spline_vol->numCoefs(0);
        int num_coefs_v = spline_vol->numCoefs(1);
        int num_coefs_w = spline_vol->numCoefs(2);
        // We must count the number of segments. We do this the easiest way
        // by looking at the number of different knots in each direction,
        // (and subtract one).
        vector<double> all_knots_u(spline_vol->basis(0).begin(), spline_vol->basis(0).end());
        vector<double> all_knots_v(spline_vol->basis(1).begin(), spline_vol->basis(1).end());
        vector<double> all_knots_w(spline_vol->basis(2).begin(), spline_vol->basis(2).end());
        vector<double>::iterator last_uniq_u = unique(all_knots_u.begin(), all_knots_u.end(), my_equal_function);
        vector<double>::iterator last_uniq_v = unique(all_knots_v.begin(), all_knots_v.end(), my_equal_function);
        vector<double>::iterator last_uniq_w = unique(all_knots_w.begin(), all_knots_w.end(), my_equal_function);
        vector<double> unique_knots_u(all_knots_u.begin(), last_uniq_u);
        vector<double> unique_knots_v(all_knots_v.begin(), last_uniq_v);
        vector<double> unique_knots_w(all_knots_w.begin(), last_uniq_w);

        int num_segs_u = unique_knots_u.size() - 1;
        int num_segs_v = unique_knots_v.size() - 1;
        int num_segs_w = unique_knots_w.size() - 1;

        vector<LRSplineVolume::Refinement3D> refs_u, refs_v, refs_w;
        int mult = 1;
        Direction3D dir = XDIR; // XDIRFIXED, refining in the x-dir (u).
        if ((num_segs_v > order_v - 1) && (num_segs_w > order_w - 1))
            for (size_t ki = 0; ki < num_segs_u; ++ki)
            {
                // We bisect the interval.
                double wgt = 0.4; // Value in the range [0.0, 1.0].
                double ref_par = wgt*unique_knots_u[ki] + (1.0 - wgt)*unique_knots_u[ki+1];

                // We insert line segments which cover order_v - 1 intervals.
                for (size_t kk = 1; kk < num_segs_w - order_w; kk+=order_w+1)
                {
                    double tmin2 = (1) ? unique_knots_w[kk] : unique_knots_w[kk+1];
                    double tmax2 = (1) ? unique_knots_w[kk+order_w] : unique_knots_w[kk+order_w+1];
                    for (size_t kj = 1; kj < num_segs_v - order_v; kj+=order_v+1)
                    {
//	      assert(unique_knots_v.size() > 3);
                        // To avoid a too regular pattern we alter start
                        // values. Maybe we should use a random function ...
                        double tmin1 = (1) ? unique_knots_v[kj] : unique_knots_v[kj+1];
                        double tmax1 = (1) ? unique_knots_v[kj+order_v] : unique_knots_v[kj+order_v+1];

                        LRSplineVolume::Refinement3D ref;
                        ref.kval = ref_par;
                        ref.start1 = tmin1;
                        ref.end1 = tmax1;
                        ref.start2 = tmin2;
                        ref.end2 = tmax2;
                        ref.d = dir;
                        ref.multiplicity = mult;
                        refs_u.push_back(ref); // Refinements in the u-dir.
                    }
                }
            }

        dir = YDIR; // XDIRFIXED, refining in the x-dir (u).
        if ((num_segs_u > order_u - 1) && (num_segs_w > order_w - 1))
            for (size_t ki = 0; ki < num_segs_v; ++ki)
            {
                // We bisect the interval.
                double wgt = 0.7; // Value in the range [0.0, 1.0].
                double ref_par = wgt*unique_knots_v[ki] + (1.0 - wgt)*unique_knots_v[ki+1];

                // We insert line segments which cover order_v - 1 intervals.
                for (size_t kk = 1; kk < num_segs_u - order_u; kk+=order_u+1)
                {
                    double tmin2 = (1) ? unique_knots_u[kk] : unique_knots_u[kk+1];
                    double tmax2 = (1) ? unique_knots_u[kk+order_u] : unique_knots_u[kk+order_u+1];
                    for (size_t kj = 1; kj < num_segs_w - order_w; kj+=order_w+1)
                    {
//		assert(unique_knots_w.size() > 3);
                        // To avoid a too regular pattern we alter start
                        // values. Maybe we should use a random function ...
                        double tmin1 = (1) ? unique_knots_w[kj] : unique_knots_w[kj+1];
                        double tmax1 = (1) ? unique_knots_w[kj+order_w] : unique_knots_w[kj+order_w+1];

                        LRSplineVolume::Refinement3D ref;
                        ref.kval = ref_par;
                        ref.start1 = tmin1;
                        ref.end1 = tmax1;
                        ref.start2 = tmin2;
                        ref.end2 = tmax2;
                        ref.d = dir;
                        ref.multiplicity = mult;
                        refs_v.push_back(ref); // Refinements in the u-dir.
                    }
                }
            }

        dir = ZDIR; // Refining in the Z-dir (w).
        if ((num_segs_u > order_u - 1) && (num_segs_v > order_v - 1))
            for (size_t ki = 0; ki < num_segs_w; ++ki)
            {
                // We bisect the interval.
                double wgt = 0.3; // Value in the range [0.0, 1.0].
                double ref_par = wgt*unique_knots_w[ki] + (1.0 - wgt)*unique_knots_w[ki+1];

                // We insert line segments which cover order_w - 1 intervals.
                for (size_t kk = 1; kk < num_segs_v - order_v; kk+=order_v+1)
                {
                    double tmin2 = (1) ? unique_knots_v[kk] : unique_knots_v[kk+1];
                    double tmax2 = (1) ? unique_knots_v[kk+order_v] : unique_knots_v[kk+order_v+1];
                    for (size_t kj = 1; kj < num_segs_u - order_u; kj+=order_u+1)
                    {
//		assert(unique_knots_u.size() > 3);
                        // To avoid a too regular pattern we alter start
                        // values. Maybe we should use a random function ...
                        double tmin1 = (1) ? unique_knots_u[kj] : unique_knots_u[kj+1];
                        double tmax1 = (1) ? unique_knots_u[kj+order_u] : unique_knots_u[kj+order_u+1];

                        LRSplineVolume::Refinement3D ref;
                        ref.kval = ref_par;
                        ref.start1 = tmin1;
                        ref.end1 = tmax1;
                        ref.start2 = tmin2;
                        ref.end2 = tmax2;
                        ref.d = dir;
                        ref.multiplicity = mult;
                        refs_w.push_back(ref); // Refinements in the u-dir.
                    }
                }
            }
     
        vector<LRSplineVolume::Refinement3D> all_refs(refs_u.begin(), refs_u.end());
        all_refs.insert(all_refs.end(), refs_v.begin(), refs_v.end());
        all_refs.insert(all_refs.end(), refs_w.begin(), refs_w.end());

        // Perform refinement
        // @@sbr201301 Remove when stable.
        lr_vol.refine(all_refs, true);

        // We compare the original with the refined sf in a randomly selected point.
        double wgt = 0.783;
        const double upar = wgt*spline_vol->startparam(0) + (1.0 - wgt)*spline_vol->startparam(0);
        const double vpar = wgt*spline_vol->startparam(1) + (1.0 - wgt)*spline_vol->startparam(1);
        const double wpar = wgt*spline_vol->startparam(2) + (1.0 - wgt)*spline_vol->startparam(2);
        Point spline_pt, lrspline_pt;
        spline_vol->point(spline_pt, upar, vpar, wpar);
        lr_vol.point(lrspline_pt, upar, vpar, wpar);
        const double dist = spline_pt.dist(lrspline_pt);
	const double tol = 1e-14;
	BOOST_CHECK_LT(dist, tol);
    }
}


// Include the test when the function is implemented.
#if 0
BOOST_FIXTURE_TEST_CASE(subVolume, Config)
{
    // Assuming all infiles are SplineVolume. Otherwise the fixture must be changed.
    for (auto iter = infiles.begin(); iter != infiles.end(); ++iter)
    {
	ifstream in1(iter->c_str());
        BOOST_CHECK_MESSAGE(in1.good(), "Input file not found or file corrupt");

	shared_ptr<SplineVolume> spline_vol(new SplineVolume());

        // BOOST_ERROR("Test hangs during call to read(). Fix!");
        // continue;

        try
        {
            header.read(in1);

            // We expect all the files to contain 1 SplineVolume.
            BOOST_CHECK_EQUAL(header.classType(), Class_SplineVolume);

            spline_vol->read(in1);
        }
        catch (...)
        {
            BOOST_ERROR("Reading of lrspline volume failed.");
            continue;
        }

        const double knot_tol = 1.0e-06;
        LRSplineVolume lr_vol(spline_vol.get(), knot_tol);

	// lr_vol subVolume() and evaluation in corner point.
	const double fuzzy = 1e-10;
	double const umin = lr_vol.startparam_u();
	double const umax = lr_vol.endparam_u();
	double const new_umin = 0.9*umin + 0.1*umax;
	double const vmin = lr_vol.startparam_v();
	double const vmax = lr_vol.endparam_v();
	double const wmin = lr_vol.startparam_w();
	double const wmax = lr_vol.endparam_w();

        Point orig_pt;
        lr_vol.point(orig_pt, new_umin, vmin, wmin);
	std::cout << "orig_pt: " << orig_pt[0] << " " << orig_pt[1] << " " << orig_pt[2] << std::endl;
	// MESSAGE("Calling subVolume().");
        std::cout << "new_umin: " << new_umin << ", vmin: " << vmin << ", wmin: " << wmin <<
            ", umax: " << umax << ", vmax: " << vmax << ", wmax: " << wmax << std::endl;

        shared_ptr<LRSplineVolume> sub_vol;
        try
        {
            sub_vol = shared_ptr<LRSplineVolume>(lr_vol.subVolume(new_umin, vmin, wmin, umax, vmax, wmax, fuzzy));
        }
        catch (...)
        {
            MESSAGE("Failed extracting subVolume().");
            continue;
        }

        std::cout << "Done calling subVolume()." << std::endl;
        // To provoke memory overwrite we call the function another time.
//	shared_ptr<LRSplineVolume> sub_vol2(lr_vol.subVolume(new_umin, vmin, umax, vmax, fuzzy));
        // MESSAGE("Done calling subVolume().");
        Point sub_pt;
        sub_vol->point(sub_pt, new_umin, vmin, wmin);
        std::cout << "sub_pt: " << sub_pt[0] << " " << sub_pt[1] << " " << sub_pt[2] << std::endl;
        const double dist = orig_pt.dist(sub_pt);
        std::cout << "dist: " << dist << std::endl;
        const double tol = 1e-14;
        BOOST_CHECK_LT(dist, tol);
    }
}
#endif
