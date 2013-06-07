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

#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineApproximator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include <vector>
#include <fstream>

using namespace std;

//===========================================================================
//                                                                           
// File: interpol_curve_free.C       
//                                                                           
/// Description:
///
/// This program reads a point data set from a file, interpolates a spline curve
/// through the points and writes an output file with the spline curve and the
///  input points.
/// The file is written in a format that can be read by the program 'goview'
/// for a graphic display of the curve and the data points.
/// The two file names are given by the user on the command line.
/// The file 'data/curve1.dat' can be used as an input file.
///
/// This example program uses uniform parametrization (0,1,...n-1), but the
/// recommended methods are chord length or centripetal parametrization.
/// The spline's endpoint conditions are set to 'Free', meaning no conditions
/// at all.
/// 
/// Points are computed by interpolating the spline at the same parameter values
/// as the input data, and the maximum distance between the original and the
/// interpolated points is written to the screen.
///
/// Input file format :
///
///  space-dimension  number-of-points
///  x1 y1 z1
///  x2 y2 z2
///     .
///     .
///  xn yn zn
///
/// The space-dimension must be equal to three if you want to use 'goview'
/// for plotting, and the number-of-points must be greater than three.
///
//===========================================================================


int main(int argc, char** argv)
{
    if (argc != 3)
    {
      cout << "\nUsage: " << argv[0]
	   << " point_infile  curve_point_outfile\n"
	   << "The outfile must end with '.g2'\n" << endl;
      exit(-1);
    }

    // Read input file
    ifstream file(argv[1]);  // Open input file.
    if (!file) {
	cerr << "\nFile error. Could not open file: " << argv[1] << endl;
	return 1;
    }
    cout << "\nRunning program " << argv[0]
	 << "\nInfile     = " << argv[1]
	 << "\nOutfile    = " << argv[2]
	 <<  '\n' << endl;

    int dim;    // Space dimension. 
    int numpt;  // Number of input points.
    file >> dim >> numpt;
    ALWAYS_ERROR_IF(dim !=3,
		    "Dimension must be 3 to look at the curve in 'goview'");
    ALWAYS_ERROR_IF(numpt < 4,
		    "At least four points are needed to make a cubic spline curve.");
    std::vector<double> param(numpt);
    std::vector<double> data(numpt * dim);
    for (int i = 0; i < numpt; ++i) {
	param[i] = static_cast<double>(i);  // Set uniform parameter values.
	for (int dd = 0; dd < dim; ++dd) {
	    file >> data[i*dim + dd];       // Read input point coordinates.
	}
    }
    file.close();

    // Create SplineInterpolator and SplineCurve
    Go::SplineInterpolator interpol;  // Create an empty SplineInterpolator.
    interpol.setFreeConditions();     // Set endpoint conditions to 'Free'
                                      // meaning no conditions at all.
                                      // Requires at least four data points.

    Go::SplineCurve cv;       // Create an empty SplineCurve.
    cv.interpolate(           // Use the input data to initialize the
		              // interpolator and the spline.
		   interpol,  // Reference to the SplineInterpolator object
		              // specifying the interpolation method to use.
		   numpt,     // Number of points to interpolate.
		   dim,       // Dimension of the Euclidean space in which the
                              // points lie (usually 2 or 3).
		   &param[0], // Pointer to the start of an array expressing the
		              // parameter values of the given data points.
		   &data[0]); // Pointer to the array where the coordinates of
                              // the data points are stored.

    // Write curve to file
    ofstream fout(argv[2]);
     // Class_SplineCurve=100 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary_data=0
    cv.writeStandardHeader(fout); // write header.
    fout << cv;    // write spline curve data.

    // Write input points to file
    // Class_PointCloud=400 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "400 1 0 4 255 0 0 255" << endl; // Header.
    fout << numpt << endl;
    for (int i = 0; i < numpt; ++i) {
	int ip = i*dim;
	Go::Point inp_point(data[ip], data[ip+1], data[ip+2]);
	fout << inp_point << endl;  // write input point coordinates.
    }
    fout.close();

    // Compute max distance between the interpolated and corresponding
    // input points.
    double max_dist = 0.0;
    for (int i = 0; i < numpt; ++i) {
	int ip = i*dim;
	Go::Point inp_point(data[ip], data[ip+1], data[ip+2]);
	Go::Point spline_point(dim);
	cv.point(spline_point, param[i]);  // Interpolate at param[i].
	double dist = inp_point.dist(spline_point);
	max_dist = std::max(dist, max_dist);
    }
    cout << "\nMaximum distance between an interpolated point and the "
	 << "corresponding input point is " << max_dist << '\n' << endl;

    return 0;
}
