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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines3D/LinDepUtils3D.h"
#include "GoTools/lrsplines3D/DefineRefs3D.h"
#include "GoTools/lrsplines3D/LRSplinePlotUtils3D.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read LR B-spline volume from specified file.
/// Identify B-splines involved in a linear dependence situation.
/// Remove the dependence by applying structured mesh refinement to the
/// main B-splines in the dependence relation.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  if (argc != 3)
    {
      std::cout << "usage: Input LR spline volume, output LR spline volume" << std::endl;
    }

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);

  // Read LR B-spline volume from file
  // Read header specifying the type of geometry entity
  GoTools::init();
  ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  // The following assumes that the specified file contains an LR B-spline
  // volume
  // Create empty volume
  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(input);
  if (!vol.get())
    {
      std::cerr << "The file contains no LR B-spline volume" << std::endl;
      exit(-1);
    }

  // Write LR mesh to file.
  std::ofstream of1("mesh_in.g2");
  writeElementLineCloud(*vol, of1);

  // Fetch candidates for a linear dependence relationship
  int deplin_lim = 8;  // Minimum number of B-splines involved 
  vector<LRBSpline3D*> funs = LinDepUtils3D::fetchUnpeelable(*vol, deplin_lim);
  std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
  
  vector<vector<LRBSpline3D*> > lindep;
  if (funs.size() >= deplin_lim)
    {
      // A possibility for a linear dependence relationship. Check if there
      // really is such a situation
      // A linear dependency relation requires at least seven B-splines situated
      // entirely within the domain of one large B-splines
      // Write the support of the large B-spline to file
      LinDepUtils3D::checkOverloaded(deplin_lim, funs, lindep);
    }
  if (lindep.size() > 0)
    std::cout << "Number linear dependence situations: " << lindep.size() << std::endl;

  if (lindep.size() > 0)
    {
      // Resolve the linear dependence by applying structured mesh refinement
      // to the main B-spline in the dependence relation
      
      int nmb_iter = 10;   // The refinement is intended to resolve the linear
      // dependence situation but might create new ones. Thus, more than one
      // refinement iteration might be required.
      for (int ka=0; ka<nmb_iter; ++ka)
	{
	  // // Remove internal B-splines in the relation
	  vector<LRBSpline3D*> source;
	  for (size_t kr=0; kr<lindep.size(); ++kr)
	    source.push_back(lindep[kr][0]);

	  // Apply structured mesh refinement to linear dependency sources
	  // First define new mesh lines
	  vector<LRSplineVolume::Refinement3D> refs_x, refs_y, refs_z;
	  bool adjust = true;   // Prefer elongation of existing mesh lines
	  // if possible
	  bool reduced = true;  // Avoid dense mesh lines if a sufficient
	  // number of refinements are defined otherwise (deviates from
	  // structured mesh)
	  DefineRefs3D::refineStructuredMesh(*vol, source, refs_x, refs_y,
					     refs_z, adjust, reduced);
	  std::cout << "Number of new mesh rectangles: " << refs_x.size() + refs_y.size() + refs_z.size() << std::endl;

	  // Perform refinement
	  for (size_t kr=0; kr<refs_x.size(); ++kr)
	    vol->refine(refs_x[kr], true);

	  for (size_t kr=0; kr<refs_y.size(); ++kr)
	    vol->refine(refs_y[kr], true);
	  
	  for (size_t kr=0; kr<refs_z.size(); ++kr)
	    vol->refine(refs_z[kr], true);
	  
	  // Check again for linear independence
	  vector<LRBSpline3D*> funs2 = LinDepUtils3D::fetchUnpeelable(*vol, deplin_lim);
	  std::cout << "Number of unpeelable functions after extra refinement: " << funs2.size() << std::endl;
	  if (funs2.size() < 8)
	    break;
	  vector<vector<LRBSpline3D*> > lindep2;
	  LinDepUtils3D::checkOverloaded(deplin_lim, funs2, lindep2);
	  std::cout << "Number of linear dependence relations: " << lindep2.size() << std::endl;

	  if (lindep2.size() > 0)
	    lindep = lindep2;
	  else
	    break;
	}
    }

  // Write final surface and to file
  vol->writeStandardHeader(output);
  vol->write(output);

  std::ofstream of2("mesh_out.g2");
  writeElementLineCloud(*vol, of2);

}

