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
#include "GoTools/lrsplines2D/LRBSpline2D.h"
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
using std::pair;

//===========================================================================
//                                                                           
/// Description:
/// Read LR spline surface from file.
/// Iterate trough all elements and B-splines in the surface.
/// Enquire various information about the surface and represent the
/// surface as a tensor product surface and as a collection of simpler
/// LR B-spline surfaces.
/// The purpose is to demonstrate various functionaliy.
/// Parts of this example is similar to investigate_Element2D.
/// For evaluation, see the example evaluateLRSurface.
/// More functionality can be found in LRSplineSurface.h
///
/// The input surface is computed by the example program
/// approximateParPointsWithLRSurf and expected to be found in
/// data/approx_lrsurf.g2
/// The tensor product surface is written to data/TP_surface.g2 and
/// the collection of simpler surfaces to data/sub_surfs.g2.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Read LR B-spline surface from file
  std::string infile("data/lrsurf_ref1.g2"); //approx_lrsurf.g2");
  std::ifstream input1(infile.c_str());
  std::ifstream input2(infile.c_str());

  // Prepare for output
  std::string outfile1("data/TP_surface.g2");
  std::string outfile2("data/sub_surf.g2");
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());

  // Two approaches is used to read the surface, one where the type of the
  // geometric entity in the file is known and one where this is not the case
  
  // Read header specifying the type of geometry entity
  ObjectHeader header1;
  try {
    header1.read(input1);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  
  // The following assumes that the specified file contains an LR B-spline
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
  // an LR B-spline surface
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

  // Read specified element
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

  // Fetch degrees and dimension of geometry space for the
  // surface
  int deg1 = surf1->degree(XFIXED);
  int deg2 = surf1->degree(YFIXED);
  int dim = surf1->dimension();
  std::cout << "Dimension: " << dim << ", degree in first parameter direction: ";
  std::cout << deg1 << ", degree in second parameter direction: " << deg2 << std::endl;

  // Enquire limits of parameter domain
  double umin = surf1->paramMin(XFIXED);  // First parameter direction
  double umax = surf1->paramMax(XFIXED);
  double vmin = surf1->paramMin(YFIXED);  // Second parameter direction
  double vmax = surf1->paramMax(YFIXED);
  std::cout << "Parameter domain: [" << umin << "," << umax << "] x [";
  std::cout << vmin << "," << vmax << "]" << std::endl;

  int num_elem = surf1->numElements();
  std::cout << "Number of elements: " << num_elem << std::endl;
  
  // Iterate through all elements
  for (LRSplineSurface::ElementMap::const_iterator el = surf1->elementsBegin();
       el != surf1->elementsEnd(); ++el)
    {
      // Example functionality, see investigate_Element2D for more
      // Enquire parameter domain of the current element
      double umin_el = el->second->umin();
      double umax_el = el->second->umax();
      double vmin_el = el->second->vmin();
      double vmax_el = el->second->vmax();
    }

  int num_bspl = surf1->numBasisFunctions();
  std::cout << "Number of B-splines: " << num_bspl << std::endl;
  
  // Iterate through all B-splines
  for (LRSplineSurface::BSplineMap::const_iterator bsp = surf1->basisFunctionsBegin();
       bsp != surf1->basisFunctionsEnd(); ++bsp)
    {
      // Example functionaly
      // Enquire knot vector in the two parameter directions
      // First fetch index of knots in the LR Mesh
      vector<int> kvec1 = bsp->second->kvec(XFIXED);
      vector<int> kvec2 = bsp->second->kvec(YFIXED);

      // Populate with actual knot values
      vector<double> knots1(kvec1.size());
      vector<double> knots2(kvec2.size());
      for (size_t ki=0; ki<kvec1.size(); ++ki)
	knots1[ki] = bsp->second->knotval(XFIXED, (int)ki);
      for (size_t ki=0; ki<kvec2.size(); ++ki)
	knots2[ki] = bsp->second->knotval(YFIXED, (int)ki);
    }

  // Split the surface into simpler pieces
  int threshold_missing = 100;  // Govern requested degree of simplicity
  double tol = 0.001;   // Used in the context of a trimmed surface
  vector<pair<shared_ptr<LRSplineSurface>, LRSplineSurface::PatchStatus> > sub_sfs =
    surf1->subdivideIntoSimpler(threshold_missing, tol);
  std::cout << "Number of sub surfaces: " << sub_sfs.size() << std::endl;
  if (sub_sfs.size() > 0)
    {
      for (size_t ki=0; ki<sub_sfs.size(); ++ki)
	{
	  sub_sfs[ki].first->writeStandardHeader(of2);
	  sub_sfs[ki].first->write(of2);
	}
    }
  
  // Check if the LR spline surface is a tensor-product surface
  bool TP = surf1->isFullTensorProduct();
  std::cout << "Surface is tensor product: " << TP << std::endl;

  // Create tensor product spline surface
  // First add knots such to obtain a tensor product mesh
  // NB! This functionality is not reversible
  surf1->expandToFullTensorProduct();

  // Fetch spline surface
  shared_ptr<SplineSurface> tpsurf(surf1->asSplineSurface());
  if (tpsurf.get())
    {
      tpsurf->writeStandardHeader(of1);
      tpsurf->write(of1);
    }
}


