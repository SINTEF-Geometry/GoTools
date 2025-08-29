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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/lrsplines2D/DefineRefs2D.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read LR B-spline surface from file.
/// Identify B-splines involved in a linear dependence situation.
/// Remove the dependence by applying structured mesh refinement to the
/// main B-splines in the dependence relation.
/// The support of the identified B-splines are written to file along
/// with the LR mesh. The output can be visualized with
/// gotools/viewlib/app/goview_vol_and_lr or gotools/viewlib/app/goview.
/// The refined surface and the corresponding mesh is written to file and can
/// be visualized with gotools/viewlib/app/goview_vol_and_lr as well as the
/// input surface.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  std::string outfile("data/lin_dep.g2");
  std::ofstream lin_out(outfile.c_str());
  std::string outsurf("data/lrsurf_resolved.g2");
  std::ofstream surf_out(outsurf.c_str());
  std::string outmesh("data/lrmesh_resolved.g2");
  std::ofstream mesh_out(outmesh.c_str());

  // Read LR B-spline surface from file
  std::string infile("../../gotools-data/lrsplines2D/examples/data/lrsurf_lindep.g2");
  std::ifstream input(infile.c_str());

    // Read header specifying the type of geometry entity
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
  // surface
  // Create empty surface
  shared_ptr<LRSplineSurface> surf(new LRSplineSurface());
  surf->read(input);
  if (!surf.get())
    {
      std::cerr << "The file contains no LR B-spline surface" << std::endl;
      exit(-1);
    }

  // Fetch candidates for a linear dependence relationship
  int deplin_lim = 8;  // Minimum number of B-splines involved 
  vector<LRBSpline2D*> funs = LinDepUtils::fetchUnpeelable(*surf, deplin_lim);
  std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
  
  vector<vector<LRBSpline2D*> > lindep;
  if (funs.size() >= deplin_lim)
    {
      // A possibility for a linear dependence relationship. Check if there
      // really is such a situation
      LinDepUtils::checkOverloaded(deplin_lim, funs, lindep);
      if (lindep.size() > 0)
	{
	  std::cout << "Number linear dependence situations: " << lindep.size() << std::endl;
	  // A linear dependency relation requires at least seven B-splines situated
	  // entirely within the domain of one large B-splines
	  // Write the support of the large B-spline to file
	  for (size_t kr=0; kr<lindep.size(); ++kr)
	    {
	      LRBSpline2D *cb = lindep[kr][0];
	      lin_out << "410 1 0 4 0 0 0 255" << std::endl;
	      lin_out << "4" << std::endl;
	      lin_out << cb->umin() << " " << cb->vmin() << " 0 ";
	      lin_out << cb->umax() << " " << cb->vmin() << " 0" << std::endl;
	      lin_out << cb->umin() << " " << cb->vmin() << " 0 ";
	      lin_out << cb->umin() << " " << cb->vmax() << " 0" << std::endl;
	      lin_out << cb->umax() << " " << cb->vmin() << " 0 ";
	      lin_out << cb->umax() << " " << cb->vmax() << " 0" << std::endl;
	      lin_out << cb->umin() << " " << cb->vmax() << " 0 ";
	      lin_out << cb->umax() << " " << cb->vmax() << " 0" << std::endl;
	    }
	}
    }

  // Write LR mesh to file.
  int colour = 1;
  writeg2Mesh(*surf, lin_out, colour);

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
	  vector<LRBSpline2D*> source;
	  for (size_t kr=0; kr<lindep.size(); ++kr)
	    source.push_back(lindep[kr][0]);

	  // Apply structured mesh refinement to linear dependency sources
	  // First define new mesh lines
	  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
	  bool adjust = true;   // Prefer elongation of existing mesh lines
	  // if possible
	  bool reduced = true;  // Avoid dense mesh lines if a sufficient
	  // number of refinements are defined otherwise (deviates from
	  // structured mesh)
	  DefineRefs2D::refineStructuredMesh(*surf, source, refs_x, refs_y,
					     adjust, reduced);
	  std::cout << "Number of new mesh lines: " << refs_x.size() + refs_y.size() << std::endl;

	  // Perform refinement
	  for (size_t kr=0; kr<refs_x.size(); ++kr)
	    surf->refine(refs_x[kr], true);

	  for (size_t kr=0; kr<refs_y.size(); ++kr)
	    surf->refine(refs_y[kr], true);
	  
	  // Check again for linear independence
	  vector<LRBSpline2D*> funs2 = LinDepUtils::fetchUnpeelable(*surf, deplin_lim);
	  std::cout << "Number of unpeelable functions after extra refinement: " << funs2.size() << std::endl;
	  if (funs2.size() < 8)
	    break;
	  vector<vector<LRBSpline2D*> > lindep2;
	  LinDepUtils::checkOverloaded(deplin_lim, funs2, lindep2);
	  std::cout << "Number of linear dependence relations: " << lindep2.size() << std::endl;

	  if (lindep2.size() > 0)
	    lindep = lindep2;
	  else
	    break;
	}
    }

  // Write final surface and to file
  surf->writeStandardHeader(surf_out);
  surf->write(surf_out);

   writeg2Mesh(*surf, mesh_out, colour);

}

