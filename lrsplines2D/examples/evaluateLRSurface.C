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
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/Element2D.h"  // Not really needed as it is
// included in LRSplineSurface.h
#include "GoTools/lrsplines2D/LRBSpline2D.h" // As for Element2D
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Demonstrates the various evaluation possibilities for an LR B-spline
/// surface
///  
/// Input to the example is the surface constructed in the example
/// program refine_surface.
/// Evaluation results are written to standard output.
/// Result of grid evaluation is written to the file data/grid_on_surf.g2.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Read LR B-spline surface from file
  std::string infile("data/lrsurf_fin.g2");
  std::ifstream input(infile.c_str());
  
  // Prepare for output
  std::string outfile("data/grid_on_surf.g2");
  std::ofstream of(outfile.c_str());
  
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

  // Dimension of geometry space
  int dim = surf->dimension();
  
  // Fetch parameter domain of surface
  double umin = surf->paramMin(XFIXED);   // First parameter direction
  double umax = surf->paramMax(XFIXED);   // First parameter direction
  double vmin = surf->paramMin(YFIXED);   // Second parameter direction
  double vmax = surf->paramMax(YFIXED);   // Second parameter direction

  // Alternative enquiry
  double umin2 = surf->startparam_u();
  double umax2 = surf->endparam_u();
  double vmin2 = surf->startparam_v();
  double vmax2 = surf->endparam_v();

  // Parameter in which to evaluate
  double upar = 0.3*umin + 0.7*umax;
  double vpar = 0.5*(vmin + vmax);
  std::cout << "Parameter pair in which to evaluate: " << upar << ", " << vpar << std::endl;

  // Evaluate without information of the element in which the parameter is
  // situated
  // Only position
  Point pos1;
  surf->point(pos1, upar, vpar);
  std::cout << "Position at specified parameter: " << pos1 << std::endl;

  // Alternative interface
  Point pos2 = surf->ParamSurface::point(upar, vpar);

  // Position and first derivatives
  // Given the surface S, the sequence is: S, Su, and Sv. If more derivatives
  // are requested the sequence continues: Suu, Suv, Svv, Suuu, Suuv, Suvv,
  // Svvv, ....
  vector<Point> pts1(3);
  int no_derivs = 1;  // Number of derivatives to evaluate, = 0 means only
  // position
  bool upar_from_right = false;  // True only if the parameter is placed
  // at the upper boundary in the first parameter direction
  bool vpar_from_right = false;  // True only if the parameter is placed
  // at the upper boundary in the second parameter direction
  surf->point(pts1, upar, vpar, no_derivs, upar_from_right, vpar_from_right);
  std::cout << "Position and first derivatives: " << pts1[0] << std::endl;
  std::cout << pts1[1] << std::endl;
  std::cout << pts1[2] << std::endl;

  // If the element corresponding to the parameter pair is known
  // Identify element. Note that this call may be time consuming if the
  // surface is large. If more than one evaluation is performed, the
  // positioning of the previous parameter pair will be used to speed up
  // the search of the current
  Element2D *elem = surf->coveringElement(upar, vpar);

   // Only position
  Point pos3;
  surf->point(pos3, upar, vpar, elem);
  
  vector<Point> pts2(3);
  surf->point(pts2, upar, vpar, no_derivs, elem, upar_from_right,
	      vpar_from_right);


  // Evaluate surface normal with and without specified element
  Point norm1, norm2;
  surf->normal(norm1, upar, vpar);
  surf->normal(norm2, upar, vpar, elem);

  // Fetch all Bsplines with the specified elemnt in its support
  const vector<LRBSpline2D*> Bsplines = elem->getSupport();
     
  // Evaluate all identified B-splines in the given parameter pair
  vector<double> Bval;
  LRSplineUtils::evalAllBSplines(Bsplines, upar, vpar, upar_from_right,
				 vpar_from_right, Bval);

  // Evaluate all identified B-splines in the given parameter pair and
  // multiply with the corresponding coefficient and scaling factor
  vector<Point> Bvalpos;
  LRSplineUtils::evalAllBSplinePos(Bsplines, upar, vpar, upar_from_right, 
				   vpar_from_right, Bvalpos);

  // Sum up the contribution from all B-splines to compute position
  Point pos4(dim), pos5(dim);
  pos4.setValue(0.0);  // Initialize
  pos5.setValue(0.0);
  for (size_t ki=0; ki<Bval.size(); ++ki)
    {
      // Fetch coefficient times scaling factor from B-spline
      Point scaled_coef = Bsplines[ki]->coefTimesGamma();

      // Compute contribution to position
      pos4 += Bval[ki]*scaled_coef;
    }
  
  for (size_t ki=0; ki<Bval.size(); ++ki)
    {
      pos5 += Bvalpos[ki];
    }
  
  std::cout << "Position from Bspline values: " << pos4 << std::endl;
  std::cout << "Position from Bspline values times coefficient: " << pos5 << std::endl;
  
  // Grid evaluation
  int num_u = 50, num_v = 50;
  double nodata = -9999;  // This parameter is obsolete in this case.
  // The function is inherited from ParamSurface, and if the surface
  // is bounded, points outside the surface domain will be given this value
  vector<double> points;  // The points are stored consecutive. The geometry
  // dimension runs fastest, then the first parameter direction and finally
  // the second direction. That is: p(0,0,x), p(0,0,y), p(0,0,z), p(1,0,x),
  // ... p(num_u-1, num_v-1, z)
  // The domain of the grid can be reduced compared to the surface
  surf->evalGrid(num_u, num_v, umin, umax, vmin, vmax, points, nodata);
  
  // Write points to file as a PointCloud.
  // NB! The current surface is a 3D surface, for a function the
  // points needs to be represented in 3D before defining the
  // PointCloud.
  PointCloud3D point_cloud(&points[0], points.size()/3);
  point_cloud.writeStandardHeader(of);
  point_cloud.write(of);

}
