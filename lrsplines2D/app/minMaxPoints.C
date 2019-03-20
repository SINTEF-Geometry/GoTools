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

  if (argc != 6)
    {
      std::cout << "Usage: " " infile curves 2D, infile curves 3D, infile surface, outfile, tolerance" << std::endl;
      exit(1);
    }

  // Open input files
  std::ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  std::ifstream is2(argv[2]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");

  std::ifstream is3(argv[3]);
  ALWAYS_ERROR_IF(is3.bad(), "Bad or no input filename");

  // Open outfile
  std::ofstream os(argv[4]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  double tol = atof(argv[5]);
  double eps = 1.0e-6;

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ObjectHeader header;
  shared_ptr<GeomObject> geom_obj;
  try {
  header.read(is3);
  geom_obj = shared_ptr<GeomObject>(Factory::createObject(header.classType()));
  geom_obj->read(is3);
  }
  catch(...)
    {
      std::cout << "WARNING: No surface found" << std::endl;
      exit(0);
    }
  
  shared_ptr<ParamSurface> surf = 
    dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  
  // Fetch curve loop
  CurveLoop loop = SurfaceTools::outerBoundarySfLoop(surf, eps);

  // Read curves
  vector<pair<shared_ptr<ParamCurve>, double> > cvs;
  while (!is.eof())
    {
      // Read curves from file
      shared_ptr<SplineCurve> cv(new SplineCurve());
      is >> header;
      cv->read(is);

      if (cv->dimension() != 2)
	{
	  std::cout << "Curve dimension different from 2" << std::endl;
	  exit(1);
	}

      shared_ptr<SplineCurve> cv2(new SplineCurve());
      is2 >> header;
      cv2->read(is2);
      Point val = cv2->ParamCurve::point(0.5*(cv2->startparam()+cv2->endparam()));
      cvs.push_back(std::make_pair(cv, val[val.dimension()-1]));
 
      Utils::eatwhite(is);
      Utils::eatwhite(is2);
    }
  
  vector<pair<Point, Point> > minpoints;
  vector<pair<Point, Point> > maxpoints;
  LRMinMax::computeMinMaxPoints(surf, cvs, tol, eps, minpoints, maxpoints);

  if (minpoints.size() > 0)
    {
      os << "400 1 0 4 0 255 0 255" << std::endl;
      os << minpoints.size() << std::endl;
      for (size_t kj=0; kj<minpoints.size(); ++kj)
	{
	  os << minpoints[kj].second << " ";
	  os << minpoints[kj].first << std::endl;
	}
    }

  if (maxpoints.size() > 0)
    {
      os << "400 1 0 4 255 0 0  255" << std::endl;
      os << maxpoints.size() << std::endl;
      for (size_t kj=0; kj<maxpoints.size(); ++kj)
	{
	  os << maxpoints[kj].second << " ";
	  os << maxpoints[kj].first << std::endl;
	}
    }
  int stop_break = 1;
}

