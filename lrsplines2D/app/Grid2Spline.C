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
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

//#define DEBUG

using namespace Go;
using std::vector;

int compare(const char *str1, char str2[][15], int nmb)
{
  for (int ki=0; ki<nmb; ++ki)
    if (strcmp(str1, str2[ki]) == 0)
      return ki;
  return -1;
}

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: grid(.asc) surface_out(.g2) info_out(.txt) tol maxiter " << std::endl;
    return -1;
  }
  int ki, kj;

  std::ifstream gridin(argv[1]);
  std::ofstream sfout(argv[2]);
  std::ofstream infoout(argv[3]);
  double AEPSGE = atof(argv[4]);
  int max_iter = atoi(argv[5]);

  double smoothwg = 0.000000001;
  int dim = 1;  // Create function
  int del = 3;  // xyz

  // Read grid and construct points
  int nmb_pts = 0;
  vector<double> data;
  double x, y, z;
  char tmp_char[20];
  int in, im;
  double x1, y1, xmid, ymid;
  double no_data = -9999;
  double size;
  bool mid1 = false, mid2 = false, nodata = false;
  
  char keywords[16][15] = {"ncols", "NCOLS", "nrows", "NROWS", "xllcenter",
			   "XLLCENTER", "xllcorner", "XLLCORNER", "yllcenter",
			   "YLLCENTER", "yllcorner", "YLLCORNER", "cellsize",
			   "CELLSIZE", "nodata_value", "NODATA_VALUE"};

  for (ki=0; ki<6; ++ki)
    {
      gridin >> tmp_char;
      int key_idx = compare(tmp_char, keywords, 16);
      int key_idx2 = (key_idx < 0) ? -1 : key_idx/2;
      
      switch (key_idx2)
	{
	case 0:
	  gridin >> in;
	  break;
	case 1:
	  gridin >> im;
	  break;
	case 2:
	  gridin >> xmid;
	  mid1 = true;
	  break;
	case 3:
	  gridin >> x1;
	  break;
	case 4:
	  gridin >> ymid;
	  mid2 = true;
	  break;
	case 5:
	  gridin >> y1;
	  break;
	case 6:
	  gridin >> size;
	  break;
	case 7:
	  gridin >> no_data;
	  nodata = true;
	  break;
	default:
	  break;
	}
    }
  if (mid1)
    x1 = xmid - 0.5*in*size;
  if (mid2)
    y1 = ymid - 0.5*im*size;
  if (!nodata)
    z = atof(tmp_char);
      
  for (kj=0, y=y1-0.5*size; kj<im; ++kj, y+=size)
    for (ki=0, x=x1-0.5*size; ki<in; ++ki, x+=size)
      {
	if (nodata)
	  gridin >> z;
	else
	  nodata = true;
	if (z == no_data)
	  continue;
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
	nmb_pts++;
      }

  double limit[2], cell_del[2];
  limit[0] = x1;
  limit[1] = y1;
  cell_del[0] = size;
  cell_del[1] = size;

  // Compute bounding box
  Point low(data[del-3], data[del-2], data[del-1]);
  Point high(data[del-3], data[del-2], data[del-1]);
  for (ki=1; ki<nmb_pts; ++ki)
    for (kj=0; kj<3; ++kj)
      {
	double tmp = data[del*ki+kj+del-3];
	low[kj] = std::min(low[kj], tmp);
	high[kj] = std::max(high[kj], tmp);
      }
  Point mid = 0.5*(low + high);
  mid[2] = 0.0;  // Only translate the xy-values

  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      data[del*ki+kj] -= mid[kj-del+3];
  limit[0] -= mid[0];
  limit[1] -= mid[1];

#ifdef DEBUG
  // Write translated surface and points to g2 format
  vector<double> data2;
  data2.reserve(nmb_pts*3);
  for (ki=0, kj=0; ki<nmb_pts; ++ki, kj+=del)
    data2.insert(data2.end(), data.begin()+kj, data.begin()+kj+3);
  PointCloud3D cloud(data2.begin(), nmb_pts);
  std::ofstream of1("translated_sf.g2");
  std::ofstream of2("translated_points.g2");
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
#endif

  std::cout << std::endl;
  std::cout << "Input grid read and pre processed. Ready for surface creation.";
  std::cout << std::endl << std::endl;
  
   int nmb_coef = 10;
  int order = 3; 
  LRSurfApprox approx(nmb_coef, order, nmb_coef, order, data, 1, AEPSGE, true, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  approx.setGridInfo(limit, cell_del);
  approx.setVerbose(true);

  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = 
    approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, max_iter);

  std::cout << std::endl;
  std::cout << "Approximation completed. Writing to output files." << std::endl;
 
  infoout << "Total number of points: " << nmb_pts << std::endl;
  infoout << "Number of elements: " << surf->numElements() << std::endl;
  infoout << "Maximum distance: " << maxdist << std::endl;
  infoout << "Average distance: " << avdist_total << std::endl;
  infoout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
  infoout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;

  if (surf.get())
    {
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (surf->dimension() == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	  
	}
#endif
	  
      // Translate
      // Update parameter domain
      double umin = surf->paramMin(XFIXED);
      double umax = surf->paramMax(XFIXED);
      double vmin = surf->paramMin(YFIXED);
      double vmax = surf->paramMax(YFIXED);

      surf->setParameterDomain(umin + mid[0], umax + mid[0],
			       vmin + mid[1], vmax + mid[1]);

      surf->writeStandardHeader(sfout);
      surf->write(sfout);

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

