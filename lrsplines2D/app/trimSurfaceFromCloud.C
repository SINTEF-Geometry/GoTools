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
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/lrsplines2D/TrimSurface.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

using namespace Go;

using std::cout;
using std::endl;
using std::vector;
using std::string;

//#define DEBUG

void print_help_text()
{
  std::cout << "TrimSurfWithPoints trims an LR B-spline surface with respect to a point cloud. \n";
  std::cout << "The surface must be a hight function parameterized on xy and the domain of the point cloud must relate to the domain of the surface. \n";
  std::cout << "A tightness parameter governs how close to the point cloud the trimming curve will pass. \n";
  std::cout << "The density and uniformity of the cloud should be reflected in this parameter. \n";
  std::cout << "The parameter should be in the range [1:8] and higher for denser point clouds. \n";
  std::cout << "The smaller number should be used only when the point cloud is very ragged towards the boundary. \n";
  std::cout << "Usage: Input surface (.g2), input cloud (.txt, .g2), tightness parameter, output surface (.g2) \n";
  std::cout << "Optional: only outer trimming (0/1). Default 0 \n";
  std::cout << "-h or --help for help text" << std::endl;
}


int main(int argc, char* argv[])
{
  for (int kii=1; kii<argc; ++kii)
    {
      string arg(argv[kii]);
      if (arg == "-h" || arg == "--help")
	{
	  print_help_text();
	  exit(0);
	}
    }

  if (argc != 5 && argc != 6)
    {
      std::cout << "ERROR: Number of parameters is not correct" << std::endl;
      print_help_text();
      return 1;
    }

    
  std::ifstream filein_sf(argv[1]);
  char *pointfile(argv[2]);
  int tightness = atoi(argv[3]);
  std::ofstream fileout_bd_sf(argv[4]);
  int only_outer = false;
  if (argc == 6)
    only_outer = atoi(argv[5]);
  
  GoTools go_tools;
  go_tools.init();
  Registrator<LRSplineSurface> r293;

  shared_ptr<GeomObject> geom_obj;
  ObjectHeader header;
  try {
    header.read(filein_sf);
  }
  catch (...)
    {
      std::cout << "ERROR: Input object not recognized. Exiting" << std::endl;
      return 1;
    }
  shared_ptr<ParamSurface> sf;
  try {
    geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
    geom_obj->read(filein_sf);
  }
  catch (...)
    {
      std::cout << "ERROR: Input surface could not be read. Exiting" << std::endl;
      return 1;
    }
  sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cout << "ERROR: Input file contains no surface" << std::endl;
      return 1;
    }

  // Read point cloud
  vector<double> data;
  vector<double> extent(6);   // Limits for points in all coordinates
  // Possible types of input files
  char keys[6][8] = {"g2", "txt", "TXT", "xyz", "XYZ", "dat"};
  int ptstype = FileUtils::fileType(pointfile, keys, 6);
  if (ptstype < 0)
    {
      std::cout << "ERROR: File type not recognized" << std::endl;
      return 1;
    }

  int nmb_pts = 0;
  std::ifstream pointsin(pointfile);
  if (ptstype == 0)
    {
     // g2 file
      ObjectHeader header;
      PointCloud3D points;
      try {
	header.read(pointsin);
	points.read(pointsin);
      }
      catch (...)
	{
	  std::cout << "ERROR: Not a valid point file" << std::endl;
	  return -1;
	}
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
      std::cout << "Domain: [" << extent[0] << ", " << extent[1] << "] x [" << extent[2];
      std::cout << ", " << extent[3] << "]" << std::endl;
      std::cout << "Range: " << extent[4] << " - " << extent[5] << std::endl;
    }
  else
    FileUtils::readTxtPointFile(pointsin, 3, data, nmb_pts, extent);


  // Make bounded surface
  shared_ptr<BoundedSurface> bd_sf;
  vector<shared_ptr<BoundedSurface> > bd_surfs;
  bool isotrim[4];
  isotrim[0] = isotrim[1] = isotrim[2] = isotrim[3] = false;
  try {
    TrimSurface::makeBoundedSurface(sf, isotrim, data, tightness,
				    bd_sf, only_outer);
  }
  catch (...)
    {
      std::cout << "ERROR: Trimming of surface failed" << std::endl;
      return 1;
    }

  if (bd_sf->dimension() > 1)
    {
      std::ofstream ofl("sf_loop.g2");
      SplineDebugUtils::writeBoundary(*bd_sf, ofl);
    }
  
  // Write output surface
  bd_sf->writeStandardHeader(fileout_bd_sf);
  bd_sf->write(fileout_bd_sf);

   return 0;
}
