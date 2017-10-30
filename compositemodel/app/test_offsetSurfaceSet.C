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

#define TEST_ALL_WORKING_CASES

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
    else
    {
        std::cout << "Input arguments : filename_sf_set offset_dist epsgeo" << std::endl;

        exit(-1);
    }

    int num_success = 0;
    int num_failures = 0;
    int num_exceptions = 0;
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
            std::cout << "offsetSurfaceSet(): Caught exception for filename " << filenames[ki] << std::endl;
            ++num_exceptions;
            continue;
        }
        
        if (status == OFFSET_OK)
        {
            std::cout << "offsetSurfaceSet(): Success!" << std::endl;
            ++num_success;
        }
        else
        {
            std::cout << "offsetSurfaceSet(): Failure! Returned status (!= 0): " << std::endl;
            ++num_failures;
            continue;
        }

        std::ofstream fileout(output_filename);
        offset_sf->writeStandardHeader(fileout);
        offset_sf->write(fileout);

        std::cout << "Input and result written to " << input_filename << " and " << output_filename << std::endl;

        // The offset dist may deviate from the required dist due to removed self intersections and kinks
        // in the surface set. We assume that the method is capable of reporting on any deviations.
        bool check_sample_distance = false;
        if (check_sample_distance)
        {
            double global_max_error = 0.0;
            double global_max_u = -1.0, global_max_v = -1.0;
            double global_max_clo_u, global_max_clo_v;
            int num_error = 0;
            for (size_t kk = 0; kk < sfs.size(); ++kk)
            {
                const RectDomain& rect_dom = sfs[kk]->containingDomain();
                const int num_samples = 7;//73;
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

    std::cout << "\nnum files: " << filenames.size() << ", num success: " << num_success << ", num failures: " <<
        num_failures << ", num_exceptions: " << num_exceptions << std::endl;
}
