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
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read point cloud from file. The given point cloud represents a function
/// with diagonal steep area and a local elevation.
/// Set appropriate parameters and approximate the points with an LR surface
/// in 1D (function). The point cloud is parameterized on its
/// x- and y-parameters.
///
/// The app PointCloud2LR performs the same operations taking a point cloud
/// as input. Several parameters are also given as input (mandatory and
/// optional). Use -h to get an explanation of all parameters.
///
/// Input to the example is read from
/// gotools/gotools-data/lrsplines2D/examples/data. The path is hardcoded.
/// Note that the example must be run from a build directory directly
/// placed under gotools to find the data. With another location, modify
/// the path.
/// The resulting surface and corresponding mesh is written to
/// data/approx_lrfunc.g2 and data/approx_mesh.eps. The point cloud in
/// g2-format (for visualization) is written to data/pointcloud.g2
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Prepare for reading the point file
  std::string infile("../../gotools-data/lrsplines2D/examples/data/pointcloud.txt");
  std::ifstream pointsin(infile.c_str());
  
  // Prepare for output
  std::string outfile0("data/pointcloud.g2");
  std::string outfile1("data/approx_lrfunc.g2");
  std::string outfile2("data/approx_mesh.eps");
  std::ofstream of0(outfile0.c_str());
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());
  
  // Possible types of input files
  // The g2 file format is described in the documentation. See under
  // "Related Pages". The txt file contains x-, y- and z-values separated
  // by comma or space
  char keys[2][8] = {"g2", "txt"};
  int ptstype = FileUtils::fileType((char*)infile.c_str(), keys, 6);
  if (ptstype < 0)
    {
      std::cout << "ERROR: File type not recognized" << std::endl;
      return 1;
    }

  // Read point cloud
  vector<double> data;  // Data points (x, y, z)
  int del = 3;   // Number of values given for each data point
  vector<double> extent(2*del);   // Limits for points in all coordinates
  int nmb_pts = 0;
  if (ptstype == 0)
    {
     // Read g2 file
      ObjectHeader header;
      PointCloud3D points;
      try {
	header.read(pointsin);
	points.read(pointsin);
      }
      catch (...)
	{
	  std::cout << "ERROR: Not a valid point file" << std::endl;
	  return 1;
	}

      // Compute point extent from a box containing all points
      BoundingBox box = points.boundingBox();
      Point low = box.low();
      Point high = box.high();
      nmb_pts = points.numPoints();
      data.insert(data.end(), points.rawData(), points.rawData()+3*nmb_pts);
      for (int ki=0; ki<3; ++ki)
	{
	  extent[2*ki] = low[ki];
	  extent[2*ki+1] = high[ki];
	}
    }
  else
    FileUtils::readTxtPointFile(pointsin, del, data, nmb_pts, extent);
  
  std::cout << "Domain in x and y: [" << extent[0] << ", " << extent[1] << "] x [" << extent[2];
  std::cout << ", " << extent[3] << "]" << std::endl;
  std::cout << "Height range: " << extent[4] << " - " << extent[5] << std::endl;

  // Write point cloud to g2 file. The point cloud can be visualized in
  // goview or goview_vol_and_lr
  PointCloud3D cloud(data.begin(), nmb_pts);
  cloud.writeStandardHeader(of0);
  cloud.write(of0);

  // Move point cloud to origo. This operation stabilizes the approximation
  // operation given large values for the x- and y- coordinates. This
  // frequently occurs when dealing with geographic point clouds
  Point mid(0.5*(extent[2*(del-3)] + extent[2*(del-3)+1]),
	    0.5*(extent[2*(del-2)] + extent[2*(del-2)+1]), 0.0);
  for (int ki=0; ki<nmb_pts; ++ki)
    for (int kj=del-3; kj<del-1; ++kj)
      {
	data[del*ki+kj] -= mid[kj-del+3];
      }

  // The approximation is performed in LRSurfApprox. Set parameters.
  int dim = del-2;   // Geometry dimension of surface (in this case function)
  double eps = 0.005;  // Approximation tolerance. Given the height range
  // and the face that the point cloud is noisy, this is an appropriate
  // value. The tolerance guides refinement and is not guaranteed to be
  // satisfied for all points
  int num_iter = 6;  // Number of iterations performing refinement and
  // point approximation
  int degree = 2;    // A quadratic surface balances smoothness with
  // shape flexibility
  int nmb_coef = 6;  // Initial number of coefficients in both parameter
  // directions. If the extent of the point cloud is far from quadratic,
  // this should be reflected in the number of coefficients in each parameter
  // direction. This is done automatically in PointCloud2LR.
  bool init_mba = false;  // An initial surface is initiated by computing
  // a tensor product spline surface by least squares approximation
  bool close_dist = false;   // Project points to the surface (function) to
  // compute distance between point and surface
  bool repar = false;  // Perform parameter iteration. Not relevant for
  // function and not active

  // Initiate approximation engine
  shared_ptr<LRSurfApprox> approx(new LRSurfApprox(nmb_coef, degree+1,
						   nmb_coef, degree+1,
						   data, dim,
						   eps, init_mba, close_dist,
						   repar));

  int tomba = 2;  // Use least square approximation with a smoothing term (LS) 
  // for the first  iteration, then switch to multilevel B-spline
  // approximation (MBA)
  approx->setSwitchToMBA(tomba);

  // In LS, apply larger smoothing weigth at the surface boundaries
  approx->setSmoothBoundary(true);
  
  // Several approaches exists in refinement. Set strategy
  int refcat = 1;  // Full span
  int alter = 0;   // Do not change strategy during approximation
  int threshold = 0;  // Do not limit the collection of new knot line segments
  // from the ones defined by the choosen strategy
  double swap = -100.0;  // Default value when alter = 0
  approx->setRefinementStrategy(refcat, alter, threshold, swap,
				refcat, threshold);

  // Ensure that the approximating surface does not overshoot with relation
  // to the point heights. This can be a problem with high variation in the
  // points and a point domain that does not cover the entire surface domain.
  double zrange = extent[5] - extent[4];
  double zfac = std::max(eps, 0.005*zrange);
  approx->addLowerConstraint(extent[4] - zfac);
  approx->addUpperConstraint(extent[5] + zfac);
  approx->setLocalConstraint(zfac);

  // Write current approximation accuracy to standard out for each iteration
  approx->setVerbose(true);
  
  // Perform approxatimation
  // Accuracy information to be set by LRSurfApprox
  double maxdist; // Maximum distance between points and surface. 
  double avdist_high;  // Average distance between points and surface for
  // the point where the distance exceeds to the tolerance. 
  double avdist; // Average distance between points and surface. 
  int nmb_out;  // Number of points with a distance larger than the tolerance.

  shared_ptr<LRSplineSurface> surf;
  try {
    surf = approx->getApproxSurf(maxdist, avdist, avdist_high,
				 nmb_out, num_iter);
  }
  catch (...)
    {
      std::cout << "ERROR: Surface approximation failed" << std::endl;
      return 1;
    }

  // Approximation accuracy
  std::cout << std::endl;
  std::cout << "Total number of points: " << nmb_pts << std::endl;
  std::cout << "Number of elements: " << surf->numElements() << std::endl;
  std::cout << "Number of coefficients: " << surf->numBasisFunctions() << std::endl;
  std::cout << "Maximum distance: " << maxdist << std::endl;
  std::cout << "Average distance: " << avdist << std::endl;
  std::cout << "Average distance for points outside of the tolerance: " << avdist_high << std::endl;
  std::cout << "Number of points outside the tolerance: " << nmb_out << std::endl;


  // Update parameter domain of surface to fit with the x- and y-coordinates
  // of the given point cloud
  // First enquire the parameter domain of the surface
  double umin = surf->paramMin(XFIXED);
  double umax = surf->paramMax(XFIXED);
  double vmin = surf->paramMin(YFIXED);
  double vmax = surf->paramMax(YFIXED);

  // Set updated parameter domain
  surf->setParameterDomain(umin + mid[0], umax + mid[0],
			   vmin + mid[1], vmax + mid[1]);

  // Write surface to file
  // The file can be visualized by goview_vol_and_lr
  // Note that the surface quality is reduced in the steep area. A function
  // typically has problems handing this steepness. The function is highly
  // refine in the area, but the diagonal configuartion of this feature
  // contradicts the properties of the spline function.
  surf->writeStandardHeader(of1);
  surf->write(of1);

  // Write mesh to file
  writePostscriptMesh(*surf, of2);
  
  return 0;
}
