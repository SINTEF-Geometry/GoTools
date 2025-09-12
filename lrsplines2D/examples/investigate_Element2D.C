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
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read LR B-spline surface from file.
/// Iterate through all elements in a surface and demonstrate available 
/// enquiries excluding those connected to scattered data approximation.
///  
/// Input to the example is the surface constructed in the example
/// program refine_lrsurf.
/// For each element, the Bezier coefs of the corresponding patch are computed
/// and stored in the file data/Bezier_coefs.g2. The corresponding patches
/// are represented as spline surfaces and written to data/Bezier_patches.g2.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Read LR B-spline surface from file
  std::string infile("data/lrsurf_fin.g2");
  std::ifstream input1(infile.c_str());
  std::ifstream input2(infile.c_str());

  // Prepare for output
  std::string outfile1("data/Bezier_coefs.g2");
  std::string outfile2("data/Bezier_patches.g2");
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());
  
  // Read header specifying the type of geometry entity
  // The function throws if the entity header is invalid
  ObjectHeader header1;
  try {
    header1.read(input1);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  // The following assumes that the specified file contains an LR spline
  // surface
  // Create empty surface
  shared_ptr<LRSplineSurface> surf1(new LRSplineSurface());
  surf1->read(input1);
  if (!surf1.get())
    {
      std::cerr << "The file contains no LR B-spline surface" << std::endl;
      exit(-1);
    }

  // Alternative reading procedure if we don't know that the file contains
  // an LR spline surface
  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read the geometry object specified in the header
 ObjectHeader header2;
  try {
    header2.read(input2);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }

  // Read specified entity
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header2.classType()));
  geom_obj->read(input2);

  // Safe cast to LRSplineSurface
  shared_ptr<LRSplineSurface> surf2 =
    dynamic_pointer_cast<LRSplineSurface, GeomObject>(geom_obj);
  if (!surf2.get())
    {
      std::cerr << "The file contains no LR spline surface" << std::endl;
      exit(-1);
    }

  // Fetch parameter domain of surface
  RectDomain dom = surf1->containingDomain();

  // Mid point of parameter domain
  double umid = 0.5*(dom.umin() + dom.umax());
  double vmid = 0.5*(dom.vmin() + dom.vmax());

  // Fetch degrees and dimension of geometry space for the
  // surface
  int deg1 = surf1->degree(XFIXED);
  int deg2 = surf1->degree(YFIXED);
  int dim = surf1->dimension();

  // Knot vector for Bezier patch
  vector<double> knots1(2*(deg1+1));
  vector<double> knots2(2*(deg2+1));
  for (int ki=0; ki<=deg1; ++ki)
    {
      knots1[ki] = 0.0;
      knots1[deg1+1+ki] = 1.0;
    }
  for (int ki=0; ki<=deg2; ++ki)
    {
      knots2[ki] = 0.0;
      knots2[deg2+1+ki] = 1.0;
    }
  
  // Iterate through all elements
  for (LRSplineSurface::ElementMap::const_iterator el = surf1->elementsBegin();
       el != surf1->elementsEnd(); ++el)
    {
      // Enquire parameter domain of the current element
      double umin = el->second->umin();
      double umax = el->second->umax();
      double vmin = el->second->vmin();
      double vmax = el->second->vmax();

      // Parameter domain area of element
      double area = el->second->area();

      // Number of B-splines having this element in its support
      int num_Bspline = el->second->nmbBasisFunctions();  // Alternative nmbSupport()

      // Fetch all B-splines having this element in its support
      // Iterators are also available
      const vector<LRBSpline2D*> Bsplines = el->second->getSupport();

      // Fetch first B-spline
      LRBSpline2D* curr_bspline = 0;
      if (num_Bspline > 0)
	{
	  curr_bspline = el->second->supportFunction(0);
	}

      // Check if the mid parameter belongs to this element
      bool is_inside = el->second->contains(umid, vmid);

      // Get all neighbouring elements to this one
      vector<Element2D*> neighbours;
      el->second->fetchNeighbours(neighbours);

      // Translate the surface patch corresponding to the current element
      // to a Bezier surface and fetch the associated coefficients
      vector<double> Bezier_coefs = el->second->unitSquareBernsteinBasis();

      // Write coefficients to file as a PointCloud.
      // NB! The current surface is a 3D surface, for a function the
      // coefficients needs to be represented in 3D before defining the
      // PointCloud.
      PointCloud3D coef_cloud(&Bezier_coefs[0], Bezier_coefs.size()/3);
      coef_cloud.writeStandardHeader(of1);
      coef_cloud.write(of1);

      // Define Bezier patch as spline surface
      shared_ptr<SplineSurface> Bezier_patch(new SplineSurface(deg1+1, deg2+1,
							       deg1+1, deg2+1,
							       &knots1[0],
							       &knots2[0],
							       &Bezier_coefs[0],
							       dim));
      Bezier_patch->writeStandardHeader(of2);
      Bezier_patch->write(of2);

    }
  
}

