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
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/BoundingBox.h"

#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description: Read LR spline function or trimmed LR spline function from
/// file. Compute equally spaced contour curves.
///
/// Input to the example is read from
/// gotools/gotools-data/lrsplines2D/examples/data. The path is hardcoded.
/// Note that the example must be run from a build directory directly
/// placed under gotools to find the data. With another location, modify
/// the path.
/// The resulting curves in geometry and parameter space is written to
/// data/iso_cvs_3D.g2 and data/iso_cvs_par.g2, respectively. In the first
/// case, the two first coordinates are given by the parameter
/// domain curves while the height value constitutes the last coordinate.
//
//===========================================================================

int main(int argc, char *argv[])
{
  // Prepare for reading the surface file
  std::string infile("../../gotools-data/lrsplines2D/examples/data/Fjoloy_surf.g2");
  std::ifstream sfin(infile.c_str());
  
  // Prepare for output
  std::string outfile1("data/iso_cvs_3D.g2");
  std::string outfile2("data/iso_cvs_par.g2");
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());
  
   // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ObjectHeader header;
  header.read(sfin);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));

  // The surface is either an LR B-spline surface (current) of a trimmed
  // version thereof. In both case, the surface is a parametric surface.
  shared_ptr<ParamSurface> surf =
    dynamic_pointer_cast<ParamSurface,GeomObject>(geom_obj);
  if (!surf.get())
    {
      std::cout << "Object one is not a surface" << std::endl;
      exit(1);
    }
  surf->read(sfin);
  sfin.close();

  int dim = surf->dimension();
  if (dim != 1)
    {
      std::cout << "The given surface is not in 1D (a function). Exiting" << std::endl;
      exit(1);
    }

  // Define height values for the contour curves
  // Find an upper bound on the extent of the height values by computing
  // the bounding box of the surface
  BoundingBox bbox = surf->boundingBox();

  const double minval = bbox.low()[0];
  const double maxval = bbox.high()[0];
  double minval2 = 10.0*((int)minval/10);
  if (minval > 0.0)
    minval2 += 1.0;
  double maxval2 = 10.0*((int)maxval/10);

  int num_val = (maxval2 - minval2)/10 + 1;
  if (num_val <= 2)
    num_val = 3;

  vector<double> iso_vals(num_val, 0);
  for (int ki = 0; ki != num_val; ++ki)
    iso_vals[ki] = minval2 + ((maxval2 - minval2)/(num_val-1)) * ki;

  // Compute iso contours
  double tol = 0.1;   // The tolerance applies to the position in the
  // x- and y- plane. The height value is excact. As the domain of
  // the surface is large, a high tolerance is feasible
  int threshold_missing = 100; // To govern an internal split of the surface
  // into a set of simpler LR B-spline surfaces

  // CurveVec is defined as
  //std::vector<std::pair<std::shared_ptr<const SplineCurve>, std::shared_ptr<const SplineCurve>>>.
  // The first curve pointer of the pair represents a 2D curve in the parameter
  // domain of the investigated LR spline function.  The second curve pointer
  // represents the corresponding 3D curve, if requested.
  const vector<CurveVec> curves = LRTraceIsocontours(surf,
  						     iso_vals,
						     threshold_missing,
  						     tol);

  // Write curves in geometry space, file. The curves can be visualized in
  // goview or goview_vol_and_lr
  for (size_t kj=0; kj<curves.size(); ++kj)
    {
      std::cout << "Height: " << iso_vals[kj] << ". Number of identified curves: " << curves[kj].size() << std::endl;
      for (size_t kr=0; kr<curves[kj].size(); ++kr)
	{
	  if (curves[kj][kr].second.get())
	    {
	      curves[kj][kr].second->writeStandardHeader(of1);
	      curves[kj][kr].second->write(of1);
	    }
	}
    }

  // Write parameter space curves to file
  for (size_t kj=0; kj<curves.size(); ++kj)
    {
      for (size_t kr=0; kr<curves[kj].size(); ++kr)
	{
	  if (curves[kj][kr].first.get())
	    {
	      curves[kj][kr].first->writeStandardHeader(of2);
	      curves[kj][kr].first->write(of2);
	    }
	}
    }

  return 0;
}

  

  

