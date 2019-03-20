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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRMinMax.h"
#include <fstream>
#include <iostream>

using namespace Go;
using std::vector;
using std::pair;

int main(int argc, char* argv[] )
{

  if (argc != 5)
    {
      std::cout << "Usage: " " infile surface, dir (+/-1), outfile, tolerance" << std::endl;
      exit(1);
    }

  // Open input files
  std::ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  int dir = atoi(argv[2]);

  // Open outfile
  std::ofstream os(argv[3]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  double tol = atof(argv[4]);
  double eps = 1.0e-6;

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  vector<shared_ptr<ParamSurface> > sfs;
  while (!is.eof())
    {
      // Read surface
      ObjectHeader header;
      shared_ptr<GeomObject> geom_obj;
      try {
	header.read(is);
	geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
	geom_obj->read(is);
  }
      catch(...)
	{
	  std::cout << "WARNING: No surface found" << std::endl;
	  exit(0);
	}
  
      shared_ptr<ParamSurface> surf = 
	dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
      sfs.push_back(surf);

      Utils::eatwhite(is);
    }
  
  vector<pair<Point, Point> > extpoints;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      vector<pair<Point, Point> > extpoints2;
      
      LRMinMax::computeExtremalPoints(sfs[ki], dir, tol, eps, extpoints2);
      extpoints.insert(extpoints.end(), extpoints2.begin(), extpoints2.end());
    }

  if (extpoints.size() > 0)
    {
      os << "400 1 0 4 155 100 0 255" << std::endl;
      os << extpoints.size() << std::endl;
      for (size_t kj=0; kj<extpoints.size(); ++kj)
	{
	  os << extpoints[kj].second << " ";
	  os << extpoints[kj].first << std::endl;
	}
    }

  int stop_break = 1;
}

