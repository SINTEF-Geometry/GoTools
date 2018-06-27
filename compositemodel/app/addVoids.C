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

#include <fstream>

#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/ParamSurface.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

int main(int argc, char* argv[] )
{
  if (argc < 4)
    {
      cout << "Usage: " << "<infile CAD (g22)> <infile voids (g22)> <outfile (g22)> <Check voids (0/1)" << endl;
      exit(-1);
    }

  char* infile(argv[1]);
  ofstream outfile(argv[argc-2]);
  int check_orient = atoi(argv[argc-1]);
  shared_ptr<SurfaceModel> outershell;
  int material_id = -1;

  CompositeModelFileHandler fileread1;
  shared_ptr<Body> body = fileread1.readBody(infile);
  if (!body.get())
    {
      outershell = fileread1.readShell(infile);
      if (!outershell.get())
	exit(1);
    }
  else
    {
      outershell = body->getOuterShell();
      material_id = body->getMaterial();
    }

  vector<shared_ptr<SurfaceModel> > allshells;
  allshells.push_back(outershell);

  for (int ki=2; ki<argc-2; ++ki)
    {
      char* invoids(argv[ki]);
      CompositeModelFileHandler fileread2;
      vector<shared_ptr<SurfaceModel> > voids = 
	fileread2.readSurfModels(invoids);

      allshells.insert(allshells.end(), voids.begin(), voids.end());
    }

  if (check_orient)
    {
      double eps = outershell->getTolerances().gap;
      shared_ptr<Body> bd(new Body(allshells[0]));
      int nmb_bd0 = allshells[0]->nmbBoundaries();
      if (nmb_bd0 > 0)
	{
	  // Outer shell not closed
	  cout << "Model not closed. Returning" << std::endl;
	  exit(1);
	}
      for (size_t ki=1; ki<allshells.size(); ++ki)
	{
	  // Fetch an arbitrary point and surface normal in the shell
	  shared_ptr<ParamSurface> surf = allshells[ki]->getSurface(0);
	  double upar, vpar;
	  Point pos = surf->getInternalPoint(upar, vpar);
	  Point norm;
	  surf->normal(norm, upar, vpar);
	  
	  // Localize the first hit with respect to all shells
	  // First modify point slightly to make sure that it is inside
	  // the shell
	  norm.normalize();
	  Point pos2 = pos + 2.0*eps*norm;

	  // Define body to use in inside test
	  // First check if the void is closed
	  int nmb_bd = allshells[ki]->nmbBoundaries();
	  if (nmb_bd > 0)
	    {
	      cout << "Void not closed. Returning" << std::endl;
	      exit(1);
	    }
	  shared_ptr<Body> bd2(new Body(allshells[ki]));
	  
	  // Check if the void is inside the outer shell
	  double dist1;
	  bool inside_outer = allshells[0]->isInside(pos2, dist1);
	  if (!inside_outer)
	    {
	      // No intersection with the outer shell
	      cout << "Void not inside outer shell. Returning" << std::endl;
	      exit(1);
	    }

	  double dist2;
	  bool inside_void = allshells[ki]->isInside(pos2, dist2);
	  if (!inside_void)
	    {
	      // Outwards pointing normal. Turn surface orientation
	      int nmb = allshells[ki]->nmbEntities();
	      vector<shared_ptr<ParamSurface> > all_sfs;
	      for (int ka=0; ka<nmb; ++ka)
		{
		  shared_ptr<ParamSurface> surf = allshells[ki]->getSurface(ka);
		  surf->swapParameterDirection();
		  all_sfs.push_back(surf);
		}
	      
	      tpTolerances tptol = allshells[ki]->getTolerances();
	      shared_ptr<SurfaceModel> shell2(new SurfaceModel(tptol.gap, 
							       tptol.gap,
							       tptol.neighbour,
							       tptol.kink, 
							       tptol.bend,
							       all_sfs));
	      allshells[ki] = shell2;
	    }
	}
    }

  shared_ptr<Body> result(new Body(allshells, material_id));

  CompositeModelFileHandler filewrite;
  filewrite.writeStart(outfile);
  filewrite.writeHeader("Body with voids", outfile);
  filewrite.writeBody(result, outfile);
  filewrite.writeEnd(outfile);

  int stop_break = 1;
}

  
