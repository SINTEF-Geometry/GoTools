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
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines3D/Element3D.h"  // Not really needed as it is
// included in LRSplineVolume.h
#include "GoTools/lrsplines3D/LRBSpline3D.h" // As for Element3D
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// Demonstrates the various evaluation possibilities for an LR spline volume
///  
/// Input to the example is the volume constructed in the example
/// program refine_lrvol.
/// Evaluation results are written to standard output.
/// Result of grid evaluation is written to data/grid_on_vol.g and grie
/// evaluation for one element is written to the file data/grid_on_elem.g2.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Read LR B-spline volume from file
  std::string infile("data/lrvol_fin.g2");
  std::ifstream input(infile.c_str());
  
  // Prepare for output
  std::string outfile("data/grid_on_vol.g2");
  std::ofstream of(outfile.c_str());
  std::string outfile2("data/grid_on_elem.g2");
  std::ofstream of2(outfile2.c_str());
  
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
  
  // The following assumes that the specified file contains an LR spline
  // volume
  // Create empty volume
  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(input);
  if (!vol.get())
    {
      std::cerr << "The file contains no LR B-spline volume" << std::endl;
      exit(-1);
    }

  // Dimension of geometry space
  int dim = vol->dimension();
  
  // Fetch parameter domain of volume
  double umin = vol->paramMin(XDIR);   // First parameter direction
  double umax = vol->paramMax(XDIR);   // First parameter direction
  double vmin = vol->paramMin(YDIR);   // Second parameter direction
  double vmax = vol->paramMax(YDIR);   // Second parameter direction
  double wmin = vol->paramMin(ZDIR);   // Third parameter direction
  double wmax = vol->paramMax(ZDIR);   // Third parameter direction

  // Alternative enquiry
  double umin2 = vol->startparam_u();
  double umax2 = vol->endparam_u();
  double vmin2 = vol->startparam_v();
  double vmax2 = vol->endparam_v();
  double wmin2 = vol->startparam_w();
  double wmax2 = vol->endparam_w();

  // Parameter in which to evaluate
  double upar = 0.3*umin + 0.7*umax;
  double vpar = 0.5*(vmin + vmax);
  double wpar = 0.55*wmin + 0.45*wmax;
  std::cout << "Parameter tripple in which to evaluate: " << upar << ", " << vpar;
  std::cout << ", " << wpar << std::endl;
 
  // Evaluate without information of the element in which the parameter is
  // situated
  // Only position
  Point pos1;
  vol->point(pos1, upar, vpar, wpar);
  std::cout << "Position: " << pos1 << std::endl;

  // Position and first derivatives
  // Given the volume S, the sequence is: S, Su, Sv, and Sw. If more derivatives
  // are requested the sequence continues: Suu, Suv, Suw, Svv, Svw, Sww, ...
  vector<Point> pts1(4);
  int no_derivs = 1;  // Number of derivatives to evaluate, = 0 means only
  // position
  bool upar_from_right = true;  // False only if the parameter is placed
  // at the upper boundary in the first parameter direction
  bool vpar_from_right = true;  // False only if the parameter is placed
  // at the upper boundary in the second parameter direction
  bool wpar_from_right = true;  // False only if the parameter is placed
  // at the upper boundary in the third parameter direction
  vol->point(pts1, upar, vpar, wpar, no_derivs, upar_from_right, vpar_from_right,
	     wpar_from_right);
  std::cout << "Position and first derivatives: " << pts1[0] << std::endl;
  std::cout << pts1[1] << std::endl;
  std::cout << pts1[2] << std::endl;
  std::cout << pts1[3] << std::endl;
  
  // If the element corresponding to the parameter tripple is known
  // Identify element. Note that this call may be time consuming if the
  // volume is large. If more than one evaluation is performed, the
  // positioning of the previous parameter pair will be used to speed up
  // the search of the current
  Element3D *elem = vol->coveringElement(upar, vpar, wpar);

   // Only position
  Point pos2;
  vol->point(pos2, upar, vpar, wpar, elem);
  
  vector<Point> pts2(4);
  vol->point(pts2, upar, vpar, wpar, no_derivs, elem, upar_from_right,
	     vpar_from_right, wpar_from_right);

  // Fetch all Bsplines with the specified elemnt in its support
  const vector<LRBSpline3D*> Bsplines = elem->getSupport();
     
  // Evaluate all identified B-splines in the given parameter pair
  vector<double> Bval;
  LRSpline3DUtils::evalAllBSplines(Bsplines, upar, vpar, wpar, upar_from_right,
				   vpar_from_right, wpar_from_right, Bval);

  // Evaluate all identified B-splines in the given parameter pair and
  // multiply with the corresponding coefficient and scaling factor
  vector<Point> Bvalpos;
  LRSpline3DUtils::evalAllBSplinePos(Bsplines, upar, vpar, wpar, upar_from_right, 
				     vpar_from_right, wpar_from_right, Bvalpos);

  // Sum up the contribution from all B-splines to compute position
  Point pos3(dim), pos4(dim);
  pos3.setValue(0.0);  // Initialize
  pos4.setValue(0.0);
  for (size_t ki=0; ki<Bval.size(); ++ki)
    {
      // Fetch coefficient times scaling factor from B-spline
      Point scaled_coef = Bsplines[ki]->coefTimesGamma();

      // Compute contribution to position
      pos3 += Bval[ki]*scaled_coef;
    }
  
  for (size_t ki=0; ki<Bval.size(); ++ki)
    {
      pos4 += Bvalpos[ki];
    }
  
  std::cout << "Position from Bspline values: " << pos3 << std::endl;
  std::cout << "Position from Bspline values times coefficient: " << pos4 << std::endl;
  
  // Grid evaluation
  int unum = 20, vnum = 20, wnum = 20;
  vector<double> points;  // The points are stored consecutive. The geometry
  // dimension runs fastest, then the first parameter direction, the second
  // parameter direction and and finally the third parameter direction. 
  // The domain of the grid can be reduced compared to the volume
  vol->evalGrid(unum, vnum, wnum, umin, umax, vmin, vmax, wmin, wmax, points);
  
  // Write points to file as a PointCloud.
  // NB! This is only possible if the geometry space of the volume is 3D
  PointCloud3D point_cloud(&points[0], points.size()/3);
  point_cloud.writeStandardHeader(of);
  point_cloud.write(of);

  // Grid evaluation in element
  // First define parameter values for evaluation in the three directions
  double umin_el = elem->umin();
  double umax_el = elem->umax();
  double vmin_el = elem->vmin();
  double vmax_el = elem->vmax();
  double wmin_el = elem->wmin();
  double wmax_el = elem->wmax();

  int num_u = 3, num_v = 4, num_w = 5;
  vector<double> par_u(num_u);
  vector<double> par_v(num_v);
  vector<double> par_w(num_w);
  double del_u = (umax_el - umin_el)/(double)(num_u-1);
  double del_v = (vmax_el - vmin_el)/(double)(num_v-1);
  double del_w = (wmax_el - wmin_el)/(double)(num_w-1);
  int ki;
  double par;
  for (ki=0, par=umin_el; ki<num_u; ++ki, par+=del_u)
    par_u[ki] = par;
  for (ki=0, par=vmin_el; ki<num_v; ++ki, par+=del_v)
    par_v[ki] = par;
  for (ki=0, par=wmin_el; ki<num_w; ++ki, par+=del_w)
    par_w[ki] = par;

  // Evaluate
  vector<double> points2;
  vol->elementGridEvaluate(elem, par_u, par_v, par_w, points2);
  
  // Write points to file as a PointCloud.
  PointCloud3D point_cloud2(&points2[0], points2.size()/3);
  point_cloud2.writeStandardHeader(of2);
  point_cloud2.write(of2);

  // Grid evaluation in element including first derivatives
  vector<double> points3;
  vol->elementGridEvaluate(elem, &par_u[0], par_u.size(), &par_v[0], par_v.size(),
			   &par_w[0], par_w.size(), no_derivs, points3);

}
