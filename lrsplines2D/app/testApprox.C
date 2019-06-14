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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 12 && argc != 13) {
    std::cout << "Usage: point cloud (.g2), lrspline_out.g2, tol, maxiter, grid (0/1), smoothing weight, MBA(0/1), toMBA(n), initMBA(0/1), set minsize(0/1), to3D(-1/n), optional:output distance field (x,y,z,d)" << std::endl;
    return -1;
  }
  int ki;

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  double AEPSGE = atof(argv[3]);
  int max_iter = atoi(argv[4]);
  int grid = atoi(argv[5]);
  double smoothwg = atof(argv[6]);//1.0e-10; //atof(argv[11]);
  int mba = atoi(argv[7]);
  int tomba = atoi(argv[8]);
  int initmba = atoi(argv[9]);
  int setmin = atoi(argv[10]);
  int to3D = atoi(argv[11]);
  char *field_out = 0;
  if (argc == 13)
    field_out = argv[12];

  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  // if (mba)
  //   to3D = -1;

  double limit[2];
  double cell_del[2];
  if (grid==1)
    {
      std::cout << "Give domain start (umin, umax): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> limit[ki];
      std::cout << "Cell size (u, v): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> cell_del[ki];

      to3D = -1;
    }

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid = 0.5*(low + high);
  Vector3D vec(-mid[0], -mid[1], 0.0);
  points.translate(vec);
  if (grid == 1)
    {
      limit[0] += vec[0];
      limit[1] += vec[1];
    }
  else if (grid==2)
    {
      bool rotate = true;
      if (rotate)
  	{
  	  double tmp[3];
  	  std::cout << "from vec:" << std::endl;
  	  for (ki=0; ki<3; ++ki)
  	    std::cin >> tmp[ki];
  	  Vector3D vec1(tmp[0], tmp[1], tmp[2]);
  	  std::cout << "to vec: " << std::endl;
  	  for (ki=0; ki<3; ++ki)
  	    std::cin >> tmp[ki];
  	  Vector3D vec2(tmp[0], tmp[1], tmp[2]);
	  vec1.normalize();
	  vec2.normalize();
  	  points.rotate(vec1, vec2);
  	}
      grid = 0;
    }

  std::cout<< "x min-max: " << low[0] << " " << high[0] << std::endl;
  std::cout<< "y min-max: " << low[1] << " " << high[1] << std::endl;
  std::cout<< "Elevation min-max: " << low[2] << " " << high[2] << std::endl;

#ifdef DEBUG
  std::ofstream of("translated_cloud.g2");
  points.writeStandardHeader(of);
  points.write(of);
#endif

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  //int dim = 1;
  int nmb_coef = 6;
  int order = 3; //4;
  double mba_coef = 0.0;
  if (initmba)
    mba_coef = 0.5*(low[2]+high[2]);
  LRSurfApprox approx(nmb_coef, order, nmb_coef, order, data, 1, AEPSGE, 
		      initmba ? true : false, mba_coef, true, true);
  //LRSurfApprox approx(4, 4, 4, 4, data, 1, AEPSGE, true, true /*false*/);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  approx.setTurn3D(to3D);
  if (grid)
    approx.setGridInfo(limit, cell_del);
  if (mba)
    approx.setUseMBA(true);
  else
    {
      // if (initmba)
      // 	approx.setInitMBA(initmba, 0.5*(low[2]+high[2]));
      approx.setSwitchToMBA(tomba);
      approx.setMakeGhostPoints(true);
    }
  approx.setVerbose(true);

  // TESTING
  //approx.addLowerConstraint(0.0);
  approx.addLowerConstraint(low[2] - 0.1*(high[2]-low[2]));
  approx.addUpperConstraint(high[2] + 0.1*(high[2]-low[2]));
  approx.setLocalConstraint(0.01);

  if (setmin)
    {
      double min_el_u, min_el_v;
      std::cout << "Size of domain in u direction: " << high[0]-low[0];
      std::cout << " Give minimum element size: " << std::endl;
      std::cin >> min_el_u;
      std::cout << "Size of domain in v direction: " << high[1]-low[1];
      std::cout << " Give minimum element size: " << std::endl;
      std::cin >> min_el_v;

      approx.setMinimumElementSize(min_el_u, min_el_v);
    }
      

  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = 
    approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, max_iter);

  std::cout << "No. elements: " << surf->numElements();
  std::cout << ", maxdist= " << maxdist << ", avdist= " << avdist_total;
  std::cout << ", avdist(out)= " << avdist;
  std::cout << ", nmb out= " << nmb_out_eps << std::endl;

  if (surf.get())
    {
      if (field_out)
	{
	  // Fetch data points with distance information
	  vector<double> pnts_dist;
	  pnts_dist.reserve(4*nmb_pts);
	   LRSplineSurface::ElementMap::const_iterator elem = surf->elementsBegin();
	   LRSplineSurface::ElementMap::const_iterator last = surf->elementsEnd();
	  for (; elem != last; ++elem)
	    {
	      if (!elem->second->hasDataPoints())
		continue;
	      vector<double>& points = elem->second->getDataPoints();
	      pnts_dist.insert(pnts_dist.end(), points.begin(), points.end());
	    }

	  // Translate to initial domain
	  for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	    {
	      pnts_dist[kj] += mid[0];
	      pnts_dist[kj+1] += mid[1];
	    }

	  // Write to file
	  std::ofstream field_info(field_out);
	  (void)field_info.precision(15);
	  for (size_t kj=0; kj<pnts_dist.size(); kj+=4)
	    {
	      for (ki=0; ki<4; ++ki)
		field_info << pnts_dist[kj+ki] << " ";
	      field_info << std::endl;
	    }
	}
	      
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (surf->dimension() == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  std::ofstream of3("translated_sf.g2");
	  surf->writeStandardHeader(of3);
	  surf->write(of3);

	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	  
	}
#endif
	  
      // Translate back
      if (surf->dimension() == 3)
	{
	  Point tmp_vec(vec.begin(), vec.end());
	  surf->translate(-tmp_vec);
	}
      else
	{
	  // Update parameter domain
	  double umin = surf->paramMin(XFIXED);
	  double umax = surf->paramMax(XFIXED);
	  double vmin = surf->paramMin(YFIXED);
	  double vmax = surf->paramMax(YFIXED);

	  surf->setParameterDomain(umin - vec[0], umax - vec[0],
				   vmin - vec[1], vmax - vec[1]);
	}

      surf->writeStandardHeader(fileout);
      surf->write(fileout);
#ifdef DEBUG
      if (surf->dimension() == 1)
	{
	  std::ofstream of2("surf_3D.g2");
	  surf->to3D();
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
#endif
    }
}

