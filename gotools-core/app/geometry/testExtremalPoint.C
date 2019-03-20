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
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ExtremalPoint.h"
#include <fstream>
#include <iostream>

using namespace Go;
using std::vector;
using std::pair;

int main(int argc, char* argv[] )
{

  if (argc < 5 || argc > 7)
    {
      std::cout << "Usage: " " infile, outfile, tolerance, direction" << std::endl;
      exit(1);
    }

  // Open input files
  std::ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  // Open outfile
  std::ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  double tol = atof(argv[3]);

  Point dir(argc-4);
  for (int ki=0; ki<argc-4; ++ki)
    dir[ki] = atof(argv[ki+4]);

  GoTools::init();

  vector<shared_ptr<ParamSurface> > sfs;
  ObjectHeader header;
  while (!is.eof())
    {
      // Read curves from file
      is >> header;
      shared_ptr<GeomObject> geom_obj  = 
	shared_ptr<GeomObject>(Factory::createObject(header.classType()));
      geom_obj->read(is);

      shared_ptr<ParamSurface> surf = 
	dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
      if (surf.get())
	sfs.push_back(surf);
      Utils::eatwhite(is);
    }

  vector<pair<Point, Point> > extrempt;
  ExtremalPoint::computeExtremalPoints(sfs, dir, tol, extrempt);

  if (extrempt.size() > 0)
    {
      os << "400 1 0 4 155 100 0 255" << std::endl;
      os << extrempt.size() << std::endl;
      for (size_t kj=0; kj<extrempt.size(); ++kj)
	{
	  if (dir.dimension() != 3)
	    os << extrempt[kj].second << " ";
	  os << extrempt[kj].first << std::endl;
	}
    }
}

