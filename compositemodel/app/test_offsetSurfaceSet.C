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

#include "GoTools/compositemodel/OffsetSurfaceUtils.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"

#include <fstream>
#include <math.h>
#include <vector>
#include <string>

using namespace Go;
using std::vector;
using std::string;
using std::ifstream;

int main( int argc, char* argv[] )
{
    vector<string> filenames;
    vector<double> offset, epsgeo;

    GoTools::init();
    
    // Test number of input arguments
    std::string output_filename("tmp/testHermiteApprEvalSurf_result.g2");
    if (argc == 5)
    {
        // Read input arguments
        std::string filename_in(argv[1]);

        double offset_val = atof(argv[2]);
        double epsgeo_val = atof(argv[3]);

        output_filename = std::string(argv[4]);

        filenames.push_back(filename_in);
        offset.push_back(offset_val);
        epsgeo.push_back(epsgeo_val);
    }
    else if (argc == 1) // No argument given, using hardcoded values.
    {

        // NEW CASES!!!

        // WORKING CASES!!!

#if 1
        filenames.push_back("data/TopSolid/TopSolid_BoundedSurf__20170623-173106.658.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

#if 0 // Added 2017-06-23
        filenames.push_back("data/TopSolid/TopSolid_SplineSurf__20170623-173106.544.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

#if 0 // Added 2017-06-23
        filenames.push_back("data/TopSolid/TopSolid_SplineSurf__20170623-173916.090.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

#if 0 // Added 2017-06-23
        filenames.push_back("data/TopSolid/TopSolid_SplineSurf__20170623-175900.205.g2");
//      filenames.push_back("data/TopSolid/TopSolid_SplineSurf__20170627-130342.267.g2"); 2017-06-17: The same data set.
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

#if 0
        // 2 bilinear quadrats with z = 0.018.
        // Tricky case which required us to to define a transition zone.
        filenames.push_back("data/TopSolid/TopSolid_Surf__20170313-174324.189_221.g2");
        offset.push_back(0.1);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;
#endif

#if 0
        // 4 planes meeting in kinks along linear segments.
        // Tricky case which required us to to define a transition zone.
        filenames.push_back("data/TopSolid/TopSolid_SplineSurf__20170504-1332_kinks.g2");
        offset.push_back(3.0e-03);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;
#endif

#if 0
        // A single surface with a bump in the interior. Use a negative offset value to create a self
        // intersection in the bump. Easy when using a positive offset dist, tricky when using a negative
        // offset dist. With a large offset value (< -3.0e-03) the self intersection is global.
        filenames.push_back("data/TopSolid/TopSolid_bump.g2");
        offset.push_back(-3.0e-03);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-05);//3);//6;
#endif

#if 0
        // Self-intersections for offset dist of appr 0.3 and larger.
        filenames.push_back("data/test_bezier.g2");
        offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;
#endif
    
#if 0
        // Self-int: offset=1e-02,eps=1e-03.
        // Success with removing self intersections using smoothing: offset=1e-03,eps=1e-04.
        filenames.push_back("data/TopSolid/TopSolid_Surf__20170404-123801.816.g2");
        offset.push_back(0.001);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;
#endif

#if 0
        // Two orthogonal planes joined by a trimmed cylinder segment (w/ radius of curvature -1.38843).
        // 201706: Success after setting the twist vector to zero and reducing precond omega from 0.1 to 0.01.
        filenames.push_back("data/Offset/yta4.g2");
        offset.push_back(-1.5);//1.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;
#endif

#if 0
        // Illegal bounded surface, the trim loop is CW! The offset works but the offset surface is flipped.
        filenames.push_back("data/TopSolid/TopSolid_BoundedSurf__20170623-173916.162.g2");
        offset.push_back(5.0e-03);
        epsgeo.push_back(1.0e-03);
#endif

        // CASES NOT SUPPORTED (EARLY EXIT WITH ERROR MESSAGE)!!!
        
#if 0
        // Degenerate patch (triangle): Ok w/ offset=1e-02,eps=1e-03. Using 0.01*epsgeo as curvature_tol.
        // Contains kink curve which is not an iso curve. Exits with failure.
        filenames.push_back("data/Offset/fanta_ro2_sub2b.g2");
        offset.push_back(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-03);//6;
#endif

#if 0
        // Tricky case with a degenerate surface. The degenerate point is not well defined with the
        // normal varying as the evaluation approaches from different iso lines. Will require extensive
        // smoothing. Currently not handled using SplineSurface (equation system too large). We need
        // LRSplineSurface (and even with that the task may be to tricky). Exits with failure.
        filenames.push_back("data/Offset/fanta_ro2_sub.g2");
        offset.push_back(0.3);//(0.01);//0.1;//1.23; //0.2;
        epsgeo.push_back(1.0e-04);//3);//6;
#endif

        // FAILING CASES!!!

    }
    else
    {
        std::cout << "Input arguments : filename_sf_set offset_dist epsgeo" << std::endl;

        exit(-1);
    }

    for (size_t ki = 0; ki < filenames.size(); ++ki)
    {
        std::cout << "\nTesting offsetSurfaceSet() for file " << filenames[ki] <<
            ", offset: " << offset[ki] << ", epsgeo: " << epsgeo[ki] << std::endl;
        
        ifstream infile(filenames[ki]);
        if (!(infile.good()))
        {
            std::cout << "Input file '" << filenames[ki] << "' not found or file corrupted!" << std::endl;
            continue;
        }

        vector<shared_ptr<ParamSurface> > sfs;
        while (!infile.eof())
        {
            ObjectHeader header;
            infile >> header;
            shared_ptr<ParamSurface> sf;
            if (header.classType() == Class_SplineSurface)
            {
                sf = shared_ptr<ParamSurface>(new SplineSurface());
            }
            else if (header.classType() == Class_BoundedSurface)
            {
                sf = shared_ptr<ParamSurface>(new BoundedSurface());
            }
            else
            {
                std::cout << "Input surface type not yet supported: " << header.classType() << std::endl;
                continue;
            }

            sf->read(infile);

            sfs.push_back(sf);

            Utils::eatwhite(infile);
        }        

        std::string input_filename("tmp/testHermiteApprEvalSurf_input.g2");
        std::ofstream fileout2(input_filename);
        for (size_t kk = 0; kk < sfs.size(); ++kk)
        {
            sfs[kk]->writeStandardHeader(fileout2);
            sfs[kk]->write(fileout2);
        }

        shared_ptr<SplineSurface> offset_sf;
        OffsetSurfaceStatus status;
        try
        {
            status = OffsetSurfaceUtils::offsetSurfaceSet(sfs, offset[ki], epsgeo[ki], offset_sf);
        }
        catch (...)
        {
            std::cout << "Caught exception for filename " << filenames[ki] << std::endl;
            continue;
        }
        
        if (offset_sf.get() == NULL)
        {
            std::cout << "Offset surface was not created. Return status: " << status;
            continue;
        }

        std::ofstream fileout(output_filename);
        offset_sf->writeStandardHeader(fileout);
        offset_sf->write(fileout);

        if (status == OFFSET_OK)
        {
            std::cout << "Success!" << std::endl;
        }
        else
        {
            std::cout << "Failure! Returned status (!= 0): " << std::endl;
            continue;
        }

        std::cout << "Input and result written to " << input_filename << " and " << output_filename << std::endl;

        if (status != 0)
        {
            std::cout << "Status not 0, something went wrong!" << std::endl;
            continue;
        }
        
        double global_max_error = 0.0;
        double global_max_u = -1.0, global_max_v = -1.0;
        double global_max_clo_u, global_max_clo_v;
        int num_error = 0;
        for (size_t kk = 0; kk < sfs.size(); ++kk)
        {
            const RectDomain& rect_dom = sfs[kk]->containingDomain();
            const int num_samples = 73;
            const double umin = rect_dom.umin();
            const double vmin = rect_dom.vmin();
            const double umax = rect_dom.umax();
            const double vmax = rect_dom.vmax();
            const double ustep = (umax - umin)/(double)(num_samples - 1);
            const double vstep = (vmax - vmin)/(double)(num_samples - 1);
            double max_error = -1.0;
            double max_u, max_v;
            double max_clo_u, max_clo_v;
            for (size_t kj = 0; kj < num_samples; ++kj)
            {
                double vpar = vmin + (double)kj*vstep;
                for (size_t kh = 0; kh < num_samples; ++kh)
                {
                    double upar = umin + (double)kh*ustep;
                    Point base_pt = sfs[kk]->point(upar, vpar);
//                    Point offset_pt = offset_sf->ParamSurface::point(upar, vpar);
                    double clo_u, clo_v, clo_dist;
                    Point offset_pt;
                    offset_sf->closestPoint(base_pt,
                                            clo_u, clo_v, offset_pt, clo_dist,
                                            epsgeo[ki]);
                    // We do not check the direction (i.e. offset sign), assuming that this is handled
                    // correctly by the method.
                    double dist = base_pt.dist(offset_pt);
                    double sign = (offset[ki] > 0.0) ? 1.0 : -1.0;
                    double error = fabs(offset[ki] - sign*dist);
                    if (error > max_error)
                    {
                        max_error = error;
                        max_u = upar;
                        max_v = vpar;
                        max_clo_u = clo_u;
                        max_clo_v = clo_v;
                    }
                    const bool disable_test = false;
                    if (!disable_test)
                    {
                        if (error > epsgeo[ki]) // Checking if error <= epsgeo.
                        {
                            ++num_error;
                            //MESSAGE("Error too large, error: " << error << ", epsgeo: " << epsgeo[ki]);
                        }
                    }
                }
            }
            if (max_error > global_max_error)
            {
                global_max_error = max_error;
                global_max_u = max_u;
                global_max_v = max_v;
                global_max_clo_u = max_clo_u;
                global_max_clo_v = max_clo_v;
            }
            std::cout << "sf # " << kk << ", max_error: " << max_error << ", max_u: " << max_u <<
                ", max_v: " << max_v << std::endl;
        }
        std::cout << "num_error: " << num_error << std::endl;
        std::cout << "ki: " << ki << ", global_max_error: " << global_max_error << ", epsgeo: " << epsgeo[ki] <<
            ", global_max_u: " << global_max_u << ", global_max_v: " << global_max_v <<
            ", global_max_clo_u: " << global_max_clo_u << ", global_max_clo_v: " << global_max_clo_v << std::endl;
    }
}
