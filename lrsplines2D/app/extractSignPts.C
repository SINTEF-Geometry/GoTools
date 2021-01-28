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
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 7) {
    std::cout << "Usage: input surface(.g2), point cloud(.g2), input extremal points(.g2), output significant points(.g2), selection (1/2), min/max level" << std::endl;
    return -1;
  }

  std::ifstream insf(argv[1]);
  std::ifstream inpts(argv[2]);
  std::ifstream inext(argv[3]);
  std::ofstream output(argv[4]);
  int select = atoi(argv[5]);
  double level = atoi(argv[6]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  //Registrator<BoundedSurface> r210;

  ObjectHeader header;
  header.read(insf);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(insf);
  
  shared_ptr<LRSplineSurface> sf = 
    dynamic_pointer_cast<LRSplineSurface, GeomObject>(geom_obj);

  if (!sf.get() || sf->dimension() != 1)
    {
      std::cout << "Dimension not equal to 1" << std::endl;
      return 1;
    }

  header.read(inpts);
  PointCloud3D points;
  points.read(inpts);

  PointCloud3D extpts;
  for (int ki=0; ki<2; ++ki)
    {
      header.read(inext);
      PointCloud3D tmp;
      tmp.read(inext);
      if (ki == select-1)
	{
	  extpts = tmp;
	  break;
	}
    }

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  if (data.size() > 0)
    LRSplineUtils::distributeDataPoints(sf.get(), data, false, 
					LRSplineUtils::REGULAR_POINTS);
  else
    return 1;


  vector<Point> sign_pts;
  int num_ext = extpts.numPoints();
  for (int kj=0; kj<num_ext; ++kj)
    {
      // Identify element
      Vector3D curr = extpts.point(kj);
      Element2D *elem = sf->coveringElement(curr[0],curr[1]);
      
      // Fetch data points
      vector<double> elem_pts = elem->getDataPoints();
      double maxh = level;
      Point curr_sign;
      if (elem_pts.size() == 0)
	{
	  // Search in neighbours
	  vector<Element2D*> next_elems;
	  elem->fetchNeighbours(next_elems);
	  for (size_t kr=0; kr<next_elems.size(); ++kr)
	    {
	      vector<double> elem_pts2 = next_elems[kr]->getDataPoints();
	      for (size_t kh=0; kh<elem_pts2.size(); kh+=3)
		{
		  if (select == 1 && elem_pts2[kh+2] < maxh)
		    {
		      maxh = elem_pts2[kh+2];
		      curr_sign = Point(elem_pts2[kh],elem_pts2[kh+1],elem_pts2[kh+2]);
		    }
		  else if (select == 2 && elem_pts2[kh+2] > maxh)
		    {
		      maxh = elem_pts2[kh+2];
		      curr_sign = Point(elem_pts2[kh],elem_pts2[kh+1],elem_pts2[kh+2]);
		    }
		}
	    }
	}
      else
	{
	  for (size_t kh=0; kh<elem_pts.size(); kh+=3)
	    {
	      if (select == 1 && elem_pts[kh+2] < maxh)
		{
		  maxh = elem_pts[kh+2];
		  curr_sign = Point(elem_pts[kh],elem_pts[kh+1],elem_pts[kh+2]);
		}
	      else if (select == 2 && elem_pts[kh+2] > maxh)
		{
		  maxh = elem_pts[kh+2];
		  curr_sign = Point(elem_pts[kh],elem_pts[kh+1],elem_pts[kh+2]);
		}
	    }
	}
      if (curr_sign.dimension() != 0)
	sign_pts.push_back(curr_sign);
      int stop_break = 1;
    }

  output << "400 1 0 4 100 100 55 255" << std::endl;
  output << sign_pts.size() << std::endl;
  for (size_t kr=0; kr<sign_pts.size(); ++kr)
    output << sign_pts[kr] << std::endl;
}
