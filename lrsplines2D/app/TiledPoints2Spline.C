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

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRApproxApp.h"
#include "GoTools/lrsplines2D/LRSurfApproxUtils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;
using std::string;

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: input meta file, output meta file, file root, tolerance, number of iterations" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  const char* file_root(argv[3]);
  double eps = atof(argv[4]);
  int max_iter = atoi(argv[5]);

  // Read meta data of tiles
  int nmb_u, nmb_v;
  double u_overlap, v_overlap;
  double total_domain[4];
  vector<double> domain;
  vector<int> nmb_points;
  vector<string> tile_files;

  LRSurfApproxUtils::readTileMeta(filein, total_domain, nmb_u, nmb_v,
				  u_overlap, v_overlap, tile_files, 
				  nmb_points, domain);

  int nmb_tiles = (int)tile_files.size();
  if (nmb_u*nmb_v != nmb_tiles)
    {
      std::cout << "Inconsistences in meta data" << std::endl;
      return -1;
    }

  // Fetch filenames for output surfaces (g2 format)
  vector<string> sf_file;
  LRSurfApproxUtils::fetchFileNames(file_root, 1, nmb_u*nmb_v, sf_file);
  
  vector<double> max_dists;
  vector<double> av_dists;
  vector<int> nmb_outside;
  double tol = 1.0e-6;   // Fuzzy domain parameter
  int ki, kj, kr;
  for (kj=0, ki=0; kj<nmb_v; ++kj)
    for (kr=0; kr<nmb_u; ++kr, ++ki)
      {
	// Read point set
	// g2-format is currently assumed
	std::ifstream infile(tile_files[ki].c_str());
	ObjectHeader header;
	header.read(infile);
	PointCloud3D points;
	points.read(infile);
	vector<double> data(points.rawData(), points.rawData()+3*nmb_points[ki]);
	
	// Surface approximation
	shared_ptr<LRSplineSurface> surf;
	double maxdist, avdist, avdist_out;
	int nmb_out;
	
	// Make restricted surface domain corresponding to non-overlapping tiles
	double dom[4];
	dom[0] = (kr == 0) ? domain[4*ki] : domain[4*ki] + u_overlap;
	dom[1] = (kr == nmb_u-1) ? domain[4*ki+1] : domain[4*ki+1] - u_overlap;
	dom[2] = (kj == 0) ? domain[4*ki+2] : domain[4*ki+2] + v_overlap;
	dom[3] = (kj == nmb_v-1) ? domain[4*ki+3] : domain[4*ki+3] - v_overlap;

	std::cout << "Approximating tile nr " << ki+1 << std::endl;
	int dim = 1;
	LRApproxApp::pointCloud2Spline(data, dim, &domain[4*ki], dom, eps, 
				       max_iter, surf, maxdist, avdist, 
				       avdist_out, nmb_out);
	max_dists.push_back(maxdist);
	av_dists.push_back(avdist);
	nmb_outside.push_back(nmb_out);
	
	if (false)
	  {
	    // If there are no points in a part of the domain, the surface may be
	    // smaller than the domain. This must be changed, but for the time
	    // being ...
	    dom[0] = std::max(dom[0], surf->startparam_u());
	    dom[1] = std::min(dom[1], surf->endparam_u());
	    dom[2] = std::max(dom[2], surf->startparam_v());
	    dom[3] = std::min(dom[3], surf->endparam_v());
	  }

	std::cout << "Extracting sub surface" << std::endl;
	shared_ptr<LRSplineSurface> sub_sf(surf->subSurface(dom[0], dom[2],
							    dom[1], dom[3], tol));

	// Surface to file
	std::ofstream of(sf_file[ki].c_str());
	(void)of.precision(15);
	sub_sf->writeStandardHeader(of);
	sub_sf->write(of);
    }

  // Write metadata
  LRSurfApproxUtils::writeSurfMeta(fileout, total_domain, nmb_u, nmb_v, 
				   eps, max_iter, nmb_points, max_dists, 
				   av_dists, nmb_outside, sf_file);
      
}
