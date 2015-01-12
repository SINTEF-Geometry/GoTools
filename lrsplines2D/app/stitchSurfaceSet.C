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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfStitch.h"
#include "GoTools/lrsplines2D/LRSurfApproxUtils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;
using std::string;
using std::pair;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: input surface metafile(.json), output meta file, file root" << std::endl;
    return -1;
  }

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);
  const char* file_root(argv[3]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  //Registrator<BoundedSurface> r210;

  // Read meta information about input surfaces
  double total_domain[4];
  int nmb_u, nmb_v;
  double eps;
  int max_iter;
  vector<int> nmb_points;
  vector<int> nmb_points_outside;
  vector<double> max_dists;
  vector<double> av_dists;
  vector<string> surf_files;
  
  LRSurfApproxUtils::readSurfMeta(input, total_domain, nmb_u, nmb_v,
				  eps, max_iter,
				  nmb_points, max_dists, av_dists,
				  nmb_points_outside, surf_files);

  // Read all surfaces. In the future intermediate file storage
  // will probably by necessary, but for the time being we assume
  // that all surfaces can be held in memory
  // Note that some surfaces may be trimmed. The underlying surfaces
  // are assumed to be LR B-spline surfaces.
  vector<shared_ptr<ParamSurface> > sfs(surf_files.size());
  for (size_t kh=0; kh<surf_files.size(); ++kh)
    {
      std::ifstream infile(surf_files[kh].c_str());
      ObjectHeader header;
      header.read(infile);
      shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
      geom_obj->read(infile);
  
      sfs[kh] = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
      if (!sfs[kh].get())
	{
	  std::cout << "Expected surface not found" << std::endl;
	  return -1;
	}
    }

  // Fetch filenames for output surfaces (g2 format)
  vector<string> sf_file;
  LRSurfApproxUtils::fetchFileNames(file_root, 1, nmb_u*nmb_v, sf_file);

  // Stitch surfaces along common edges.
  // Assosiated corners will be handled first
  double tol = eps;
  int kj, kr;
  int nmb_modified;
  bool matched = false;
  for (kj=1; kj<nmb_v; ++kj)
    {
      for (kr=0; kr<=nmb_u; ++kr)
      {
	if (kj == 1 && kr > 0 && kr < nmb_u)
	  {
	    // Match corner along the lower boundary
	    // Corners are numbered: 0=lower left, 1=lower right,
	    // 2=upper left, 3=upper right
	    vector<pair<shared_ptr<ParamSurface>, int> > corner_match1(2);
	    corner_match1[0] = make_pair(sfs[kr-1], 1);
	    corner_match1[1] = make_pair(sfs[kr], 0);
	    nmb_modified = LRSurfStitch::averageCorner(corner_match1, tol);
	    std::cout << "Nmb surfaces = " << corner_match1.size();
	    std::cout << ", nmb modified = " << nmb_modified << std::endl;
	  }
	if (kj == nmb_v-1 && kr > 0 && kr < nmb_u)
	  {
	    // Match corner along the upper boundary
	    // Corners are numbered: 0=lower left, 1=lower right,
	    // 2=upper left, 3=upper right
	    vector<pair<shared_ptr<ParamSurface>, int> > corner_match2(2);
	    corner_match2[0] = make_pair(sfs[kj*nmb_u+kr-1], 3);
	    corner_match2[1] = make_pair(sfs[kj*nmb_u+kr], 2);
	    nmb_modified = LRSurfStitch::averageCorner(corner_match2, tol);
	    std::cout << "Nmb surfaces = " << corner_match2.size();
	    std::cout << ", nmb modified = " << nmb_modified << std::endl;
	  }
	
	// Match inner corner
	vector<pair<shared_ptr<ParamSurface>, int> > corner_match3;
	if (kr > 0)
	  corner_match3.push_back(make_pair(sfs[(kj-1)*nmb_u+kr-1], 3));
	if (kr < nmb_u)
	  corner_match3.push_back(make_pair(sfs[(kj-1)*nmb_u+kr], 2));
	if (kr > 0)
	  corner_match3.push_back(make_pair(sfs[kj*nmb_u+kr-1], 1));
	if (kr < nmb_u)
	  corner_match3.push_back(make_pair(sfs[kj*nmb_u+kr], 0));
	nmb_modified = LRSurfStitch::averageCorner(corner_match3, tol);
	std::cout << "Nmb surfaces = " << corner_match3.size();
	std::cout << ", nmb modified = " << nmb_modified << std::endl;

	// Match vertical edges
	// Edges are numbered: 0=left, 1=right, 2=lower, 3=upper
	if (kr > 0 && kr < nmb_u)
	  matched = LRSurfStitch::averageEdge(sfs[(kj-1)*nmb_u+kr-1], 1,
					      sfs[(kj-1)*nmb_u+kr], 0, tol);
	if (kj == nmb_v-1 && kr > 0 && kr < nmb_u)
	  matched = LRSurfStitch::averageEdge(sfs[kj*nmb_u+kr-1], 1,
					      sfs[kj*nmb_u+kr], 0, tol);

	// Match horizontal edge
	if (kr < nmb_u)
	  matched = LRSurfStitch::averageEdge(sfs[(kj-1)*nmb_u+kr], 3,
					      sfs[kj*nmb_u+kr], 2, tol);
      }
    }

  for (size_t kh=0; kh<sfs.size(); ++kh)
    {
      // Surface to file
      std::ofstream of(sf_file[kh].c_str());
      (void)of.precision(15);
      sfs[kh]->writeStandardHeader(of);
      sfs[kh]->write(of);
    }
  
  // Write metadata
  LRSurfApproxUtils::writeSurfMeta(output, total_domain, nmb_u, nmb_v, 
				   eps, max_iter, nmb_points, max_dists, 
				   av_dists, nmb_points_outside, sf_file);
}
