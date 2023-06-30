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

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/utils/BoundingBox.h";
#include "GoTools/utils/Point.h";
#include "GoTools/igeslib/IGESconverter.h"
#include <iostream>
#include <fstream>

//===========================================================================
//                                                                           
// File: intersect_with_plane.C                                                      
//                                                                           
/// Description:
///
/// This program demonstrates how to read from and write to IGES files.
/// It shows how to fetch geometry entities from the IGES converter, compute the
/// bounding box of a surface and compute the intersection between a surface and
/// a plane using functionality in the namespace BoundedUtils.
///
/// The input surface and constructed intersection curve can be visualized in the
/// application goview in viewlib.
///
//===========================================================================

using namespace Go;

int main(int argc, char** argv)
{
  std::cout << "\nRunning program " << argv[0] << std::endl;

    // Define input file and IGES converter
  std::string file_in("data/surface.igs");
  std::ifstream cfile(file_in.c_str());
  if (!cfile) {
    std::cout << "\nFile error. Could not open file: " << file_in.c_str()
	      << std::endl;
    return 1;
  }
    
    IGESconverter conv;

    // Read data into the IGES converter
    try {
    conv.readIGES(cfile);
    } catch (...) {
      std::cout << "Failed reading input IGES file, exiting." << std::endl;
      return 1;  
    }

    // Fetch the geometry entities stored in the converter.
    // The entities are stored as GeomObject, which is the root of the geometry
    // inheritance tree.
    std::vector<shared_ptr<GeomObject> > geom = conv.getGoGeom();

    // We expect one surface
    if (geom.size() != 1)
      {
	std::cout << "One entity expected, exiting."  << std::endl;
	return 1;
      }

    // Safe cast to parametric surface
    shared_ptr<ParamSurface> surf =
      dynamic_pointer_cast<ParamSurface,GeomObject>(geom[0]);

    // Compute bounding box of surface
    BoundingBox bbox = surf->boundingBox();

    // Define the plane with which to interset
    // Fetch the mid point of the bounding box
    Point mid = 0.5*(bbox.low() + bbox.high());
    Point norm(0.0, 0.0, 1.0);  // The z-axis

    // Intersect
    double tol = 1.0e-5;   // Intersection tolerance
    std::vector<shared_ptr<CurveOnSurface> > intcvs =
		BoundedUtils::intersectWithPlane(surf,
						 mid, norm,
						 tol);

    // A CurveOnSurface instance stores a curve in the parameter domain of a
    // surface and the corresponding geometry space curve. Both curve instances
    // need not to be defined.
    // Check if any geometry curve exist and write them to an IGES file
    std::ofstream ofile("data/intersection_curves.igs");
    IGESconverter conv2;
    for (size_t ki=0; ki<intcvs.size(); ++ki)
      {
	if (intcvs[ki]->hasSpaceCurve())
	  {
	    conv2.addGeom(intcvs[ki]->spaceCurve());
	  }
      }

    if (conv2.num_geom() > 0)
      conv2.writeIGES(ofile);
   
}
