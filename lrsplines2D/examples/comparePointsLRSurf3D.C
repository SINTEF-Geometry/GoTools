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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/FileUtils.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read parameterized point cloud and corresponding surface from file.
/// Compute distances between points and surface. Prepare for visualization.
///
/// The input point cloud is read from
/// gotools/gotools-data/lrsplines2D/examples/data. The path is hardcoded.
/// Note that the example must be run from a build directory directly
/// placed under gotools to find the data. With another location, modify
/// the path.
/// The input surface is computed by the example program
/// approximateParPointsWithLRSurf and expected to be found in data/approx_lrsurf.g2
/// Average distance and largest negative and positive distance between the
/// surface and the point cloud is written to standard output.
/// The points are grouped and assigned a colour according to the distance
/// to the surface. The grouped points are written do data/coloured_points.g2
/// and can be visualized in goview_vol_and_lr.
//                                                                           
//===========================================================================

// rgb colours: Largest distance below, distance zero, Largest distance above
int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main(int argc, char *argv[])
{
  // Prepare for reading the point file
  std::string infile1("../../gotools-data/lrsplines2D/examples/data/param_pointcloud.txt");
  std::string infile2("data/approx_lrsurf.g2");
  std::ifstream pointsin(infile1.c_str());
  std::ifstream sfin(infile2.c_str());
  
  // Prepare for output
  std::string outfile("data/coloured_points.g2");
  std::ofstream of(outfile.c_str());
  
  // Read point file. The file contains u-, and
  // v-parameters and x-, y- and z-values in this sequence. The values
  // are separated by comma or space
  // Read point cloud
  vector<double> data;  // Data points (x, y, z)
  int n_val = 5;   // Number of values given for each data point
  vector<double> extent(2*n_val);   // Limits for points in all coordinates
  int nmb_pts = 0;
  FileUtils::readTxtPointFile(pointsin, n_val, data, nmb_pts, extent);
  
  std::cout << "Parameter domain: [" << extent[0] << ", " << extent[1] << "] x [" << extent[2] << ", " << extent[3] << "]" << std::endl;;

   // Read header specifying the type of geometry entity
  ObjectHeader header;
  try {
    header.read(sfin);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(1);
    }
  
  // Read surface
  // The following assumes that the specified file contains an LR B-spline
  // surface
  // Create empty surface
  shared_ptr<LRSplineSurface> surf(new LRSplineSurface());
  surf->read(sfin);
  if (!surf.get())
    {
      std::cerr << "The file contains no LR B-spline surface" << std::endl;
      exit(1);
    }

  // Dimension of geometry space
  int dim = surf->dimension();
  if (dim != 3)
    exit(1);   // Unexpected geometry space dimension

  int nmb_points; // Number of points
  double max_above, max_below, avdist, avdist_abs;  // Maximum and average 
  // distance between surface and point cloud
  vector<double> pointsdist;  // For each point: x, y, z, distance to surface
  bool close_dist = true;  // Compute closest point in surface to given
  // point and examine the distance

  // Compute distances and statistics
  LRApproxApp::computeDistPointSpline3D(data, surf, max_above, max_below,
					avdist, avdist_abs, nmb_points,
					pointsdist, close_dist);

  std::cout << "Number of points: " << nmb_points << std::endl;
  std::cout << "Maximum distance below the surface: " << max_below << std::endl;
  std::cout << "Maximum distance above the surface: " << max_above << std::endl;
  std::cout << "Average distance: " << avdist << std::endl;
  std::cout << "Average absolute distance: " << avdist_abs << std::endl;

  // Classify points according to distance
  int nmb_level = 4;  // Number of positive distance levels
  double max_level = std::max(fabs(max_below), max_above);
  double del = max_level/(double)(nmb_level+1);

  vector<double> limits(2*nmb_level+1);  // Limit distance between groups
  vector<vector<double> > level_points(2*nmb_level+2); // x, y, z-values
  // for each group

  // Set distance levels
  limits[nmb_level] = 0;
  for (int ki=1; ki<=nmb_level; ++ki)
    {
      limits[nmb_level-ki] = -ki*del;
      limits[nmb_level+ki] = ki*del;
    }

  int nval = 4;  // Number of value for each entry in pointsdist
  for (size_t ki=0; ki<pointsdist.size(); ki+=nval)
    {
      double val = pointsdist[ki+3];
      
      // Find classification
      size_t ka;
      for (ka=0; ka<limits.size(); ++ka)
	if (val < limits[ka])
	  {
	    level_points[ka].push_back(pointsdist[ki]);
	    level_points[ka].push_back(pointsdist[ki+1]);
	    level_points[ka].push_back(pointsdist[ki+2]);
	    break;
	  }
      if (ka == limits.size())
	{
	  level_points[ka].push_back(pointsdist[ki]);
	  level_points[ka].push_back(pointsdist[ki+1]);
	  level_points[ka].push_back(pointsdist[ki+2]);
	}
    }

  // Write point groups to file
  for (size_t ki=0; ki<level_points.size(); ++ki)
    {
      double lim = (ki < level_points.size()-1) ? limits[ki] : max_above;
      std::cout << "Upper level: " << lim << ", no. of pts: " << level_points[ki].size()/3 << std::endl;
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      double cc[3];
      if (ki <= nmb_level)
      	{
      	  cc[0] = ((nmb_level-ki)*colors[0][0] + ki*colors[1][0])/nmb_level;
      	  cc[1] = ((nmb_level-ki)*colors[0][1] + ki*colors[1][1])/nmb_level;
      	  cc[2] = ((nmb_level-ki)*colors[0][2] + ki*colors[1][2])/nmb_level;
      	}
      else
      	{
      	  cc[0] = ((ki-nmb_level-1)*colors[2][0] + 
      		   (2*nmb_level-ki+1)*colors[1][0])/nmb_level;
      	  cc[1] = ((ki-nmb_level-1)*colors[2][1] + 
      		   (2*nmb_level-ki+1)*colors[1][1])/nmb_level;
      	  cc[2] = ((ki-nmb_level-1)*colors[2][2] + 
      		   (2*nmb_level-ki+1)*colors[1][2])/nmb_level;
      	}

      // Write header and colour
      of << "400 1 0 4 " << cc[0] << " " << cc[1];
      of << " " << cc[2] << " 255" << std::endl;

      // Write points
      level_cloud.write(of);
    }
}



  
