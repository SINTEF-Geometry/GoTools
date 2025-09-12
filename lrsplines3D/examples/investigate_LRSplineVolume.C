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
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Read LR spline volume from file.
/// Iterate trough all elements and B-splines in the volume.
/// Enquire various information about the volume and represent the
/// volume as a tensor product volume. 
/// The purpose is to demonstrate various functionaliy.
/// Parts of this example is similar to investigate_Element3D.
/// For evaluation, see the example
/// More functionality can be found in LRSplineVolume.h
///
/// The input volume is computed by the example program
/// approximateParPointsWithLRSurf and expected to be found in
/// data/approx_lrsurf.g2
/// The tensor product volume is written to data/TP_volume.g2 and
/// a sub volume is written to data/sub_vol.g2. A constant parameter
/// surface is written to data/iso_surf.g2
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Prepare for input
  std::string infile("data/lrvol_fin.g2");
  std::ifstream input1(infile.c_str());
  std::ifstream input2(infile.c_str());

  // Prepare for output
  std::string outfile1("data/TP_volume.g2");
  std::string outfile2("data/sub_vol.g2");
  std::string outfile3("data/iso_surf.g2");
  std::ofstream of1(outfile1.c_str());
  std::ofstream of2(outfile2.c_str());
  std::ofstream of3(outfile3.c_str());

  //  Read LR spline volume from file
    // Two approaches is used to read the volume, one where the type of the
  // geometric entity in the file is known and one where this is not the case
  
  shared_ptr<LRSplineVolume> vol1(new LRSplineVolume());  // Empty volume
  
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
  
  try {
    vol1->read(input1);
  } catch (const std::exception& e) {
    std::cerr << "Error reading LRSplineVolume: " << e.what() << std::endl;
    return 1;
  }

  // Alternative reading procedure if we don't know that the file contains
  // an LR B-spline volume
  // Create the default factory
  GoTools::init();
  Registrator<LRSplineVolume> r793;

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

  // Safe cast to LRSplineVolume
  shared_ptr<LRSplineVolume> vol2 =
    dynamic_pointer_cast<LRSplineVolume, GeomObject>(geom_obj);
  if (!vol2.get())
    {
      std::cerr << "The file contains no LR spline volume" << std::endl;
      exit(-1);
    }

  // Fetch degrees and dimension of geometry space for the volume
  int deg1 = vol2->degree(XDIR);
  int deg2 = vol2->degree(YDIR);
  int deg3 = vol2->degree(ZDIR);
  int dim = vol2->dimension();
  std::cout << "Dimension: " << dim << ", degrees in the three parameter directions: ";
  std::cout << deg1 << ", " << deg2 << ", " << deg3 << std::endl;
  
  // Enquire limits of parameter domain
  double umin = vol2->paramMin(XDIR);  // First parameter direction
  double umax = vol2->paramMax(XDIR);
  double vmin = vol2->paramMin(YDIR);  // Second parameter direction
  double vmax = vol2->paramMax(YDIR);
  double wmin = vol2->paramMin(ZDIR);  // Third parameter direction
  double wmax = vol2->paramMax(ZDIR);
  std::cout << "Parameter domain: [" << umin << "," << umax << "] x [";
  std::cout << vmin << "," << vmax << "] x [" << wmin << "," << wmax << "]" << std::endl;

  // Alternative approaches to fetch the limits of the parameter domain
  const Array<double,6> dom = vol2->parameterSpan();  // The sequence: umin, umax, vmin, vmax, wmin, wmax
  double umin2 = vol2->startparam_u();
  double umax2 = vol2->endparam_u();  // Same for v and w
    
  
  int num_elem = vol2->numElements();
  std::cout << "Number of elements: " << num_elem << std::endl;
  
  // Iterate through all elements
  for (LRSplineVolume::ElementMap::const_iterator el = vol2->elementsBegin();
       el != vol2->elementsEnd(); ++el)
    {
      // Example functionality, see investigate_Element3D for more
      // Enquire parameter domain of the current element
      double umin_el = el->second->umin();
      double umax_el = el->second->umax();
      double vmin_el = el->second->vmin();
      double vmax_el = el->second->vmax();
      double wmin_el = el->second->wmin();
      double wmax_el = el->second->wmax();
    }

  int num_bspl = vol2->numBasisFunctions();
  std::cout << "Number of B-splines: " << num_bspl << std::endl;
  
  // Iterate through all B-splines
  for (LRSplineVolume::BSplineMap::const_iterator bsp = vol2->basisFunctionsBegin();
       bsp != vol2->basisFunctionsEnd(); ++bsp)
    {
      // Example functionaly
      // Enquire knot vector in the three parameter directions
      // First fetch index of knots in the LR Mesh
      vector<int> kvec1 = bsp->second->kvec(XDIR);
      vector<int> kvec2 = bsp->second->kvec(YDIR);
      vector<int> kvec3 = bsp->second->kvec(ZDIR);

      // Populate with actual knot values
      vector<double> knots1(kvec1.size());
      vector<double> knots2(kvec2.size());
      vector<double> knots3(kvec3.size());
      for (size_t ki=0; ki<kvec1.size(); ++ki)
	knots1[ki] = bsp->second->knotval(XDIR, (int)ki);
      for (size_t ki=0; ki<kvec2.size(); ++ki)
	knots2[ki] = bsp->second->knotval(YDIR, (int)ki);
      for (size_t ki=0; ki<kvec3.size(); ++ki)
	knots3[ki] = bsp->second->knotval(ZDIR, (int)ki);
    }


  // Evaluate volume in one point. Let the algorithm find the associated element
  // Define parameter tripple
  double upar = 0.5*(umin + umax);
  double vpar = 0.2*vmin + 0.8*vmax;
  double wpar = 0.65*wmin + 0.35*wmax;
  Point pos;
  vol2->point(pos, upar, vpar, wpar);
  std::cout << "Position in parameter (" << upar << "," << vpar << ",";
  std::cout << wpar <<") is: " << pos << std::endl;

  // Closest point
  double clo_u, clo_v, clo_w;  // Parameter value to closest point
  Point clo_pt;     // Closest point to input point in the volume
  double clo_dist;  // Distance between input point and closest point
  double eps = 1.0e-6;  // Used in computation

  // Do the computation. A seed could alternatively be provided
  vol2->closestPoint(pos, clo_u, clo_v, clo_w, clo_pt, clo_dist, eps);
  std::cout << "Closest parameter: " << clo_u << ", " << clo_v << ", " << clo_w << std::endl;
  std::cout << "Distance between input point and closest point: " << clo_dist << std::endl;
  
  // Pick a part of the volume
  // The sub volume is specified as: umin, vmin, wmin, umax, vmax, wmax
  // The last parameter indicated equality of parameter values to avoid subdividing
  // very close to a knot
  shared_ptr<LRSplineVolume> sub_vol(vol2->subVolume(umin, vpar, wmin, upar, vmax, wmax, eps));
  if (sub_vol.get())
    {
      sub_vol->writeStandardHeader(of2);
      sub_vol->write(of2);
    }

  // Fetch a constant parameter surface in the third parameter direction
  shared_ptr<LRSplineSurface> iso_surf(vol2->constParamSurface(wpar, 2));
  if (iso_surf.get())
    {
      iso_surf->writeStandardHeader(of3);
      iso_surf->write(of3);
    }
  
  // Check if the LR spline volume is a tensor-product volume
  bool TP = vol2->isFullTensorProduct();
  std::cout << "Volume is tensor product: " << TP << std::endl;

  // Create tensor product spline surface
  // First add knots such to obtain a tensor product mesh
  // NB! This functionality is not reversible
  // Alternatively the function asSplineVolume() can be called directly. Then
  // the initial volume is not altered
  vol2->expandToFullTensorProduct();

  // Fetch spline volume
  shared_ptr<SplineVolume> tpvol(vol2->asSplineVolume());
  if (tpvol.get())
    {
      tpvol->writeStandardHeader(of1);
      tpvol->write(of1);
    }
}
