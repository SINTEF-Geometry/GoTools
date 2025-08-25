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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <fstream>

using namespace Go;

//===========================================================================
//                                                                           
/// Description:
/// The program defines a quadratic spline surface and represents this 
/// surface as an LR B-spline surface. This surface is refined, first in the 
/// second parameter direction, then in the first parameter direction.
///  
/// Input to the geometry construction is hardcoded
/// Current surfaces and meshes are written to g2- and eps-files as we go
/// along. LR B-spline surfaces are visualized by
/// gotools/viewlib/app/goview_vol_and_lr. Meshes are visualized by
/// ghostview (gv). The final surface is written to the file
/// data/lrsurf_fin.g2. The final mesh can be found at data/lrmesh_fin.eps. 
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Prepare for output files
  std::string outfile1("data/spline_sf_in.g2");
  std::string outfile2("data/lrsurf_in.g2");
  std::string outfile3("data/lrsurf_ref1.g2");
  std::string outfile4("data/lrsurf_fin.g2");
  std::string outfile5("data/lrmesh_1.eps");
  std::string outfile6("data/lrmesh_2.eps");
  std::string outfile7("data/lrmesh_fin.eps");

  // Define biquadratic spline surface with three inner knots in the first
  // parameter direction and two in the second direction
  int dim = 3;  // Dimension of geometry space
  int in1 = 6, in2 = 5;  // Number of coefficients in the two parameter directions
  int deg1 = 2, deg2 = 2;   // Polynomial degree
  double knots1[9] = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0};
  double knots2[8] = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0};
  double coefs[90] = {0.0, 0.0, 0.0,
		      1.0, 0.0, 0.0,
		      2.0, 0.0, 0.0,
		      3.0, 0.0, 0.0,
		      4.0, 0.0, 0.0,
		      5.0, 0.0, 0.0,
		      0.0, 1.0, 0.0,
		      1.0, 1.0, 0.0,
		      2.0, 1.0, 0.0,
		      3.0, 1.0, 0.0,
		      4.0, 1.0, 0.0,
		      5.0, 1.0, 0.0,
		      0.0, 2.0, 0.0,
		      1.0, 2.0, 0.0,
		      2.0, 2.0, 0.0,
		      3.0, 2.0, 2.0,
		      4.0, 2.0, 2.0,
		      5.0, 2.0, 0.0,
		      0.0, 3.0, 0.0,
		      1.0, 3.0, 0.0,
		      2.0, 3.0, 3.0,
		      3.0, 3.0, 3.0,
		      4.0, 3.0, 0.0,
		      5.0, 3.0, 0.0,
		      0.0, 4.0, 0.0,
		      1.0, 4.0, 0.0,
		      2.0, 4.0, 1.0,
		      3.0, 4.0, 1.0,
		      4.0, 4.0, 0.0,
		      5.0, 4.0, 0.0};

  // Create spline surface
  shared_ptr<SplineSurface> spline_init(new SplineSurface(in1, in2,
							  deg1+1, deg2+1,
							  knots1, knots2, 
							  coefs, dim));
  std::ofstream of1(outfile1.c_str());
  spline_init->writeStandardHeader(of1);
  spline_init->write(of1);

  // Represent the initial spline surface as an LR B-spline surface
  // The second parameter indicates when two consequtive knots are
  // regarded as identical
  shared_ptr<LRSplineSurface> lr_surf(new LRSplineSurface(spline_init.get(),
							  1.0e-6));
  std::ofstream of2(outfile2.c_str());
  lr_surf->writeStandardHeader(of2);
  lr_surf->write(of2);

  std::ofstream of5(outfile5.c_str());
  writePostscriptMesh(*lr_surf, of5);

  // Refine in the second parameter direction
  Direction2D dir = YFIXED;    // The constant direction of the new
  // knot line segment
  double par_y1 = 0.5;  // Position of new knot in the corresponding knot vector
  double par_y2 = 1.5;  // Position of new knot in the corresponding knot vector
  double start_y = 0.0; // The start value of the new segment
  double end_y = 3.0;   // The end value of the new segment
  int mult = 1;         // Knot multiplicity
  bool absolute = true; // true = Do not increment multiplicity during
                        // multiple knot insertions,
                        // false = increment multiplicity

  // Perform knot insertion, one knot line segment at the time
  lr_surf->refine(dir, par_y1, start_y, end_y, mult, absolute);

  // Collect knot insertion information in a structure. This call has the
  // same effect as the previous
  LRSplineSurface::Refinement2D ref_y2;
  ref_y2.setVal(par_y2, start_y, end_y, dir, mult);
  lr_surf->refine(ref_y2, absolute);
  
  std::ofstream of3(outfile3.c_str());
  lr_surf->writeStandardHeader(of3);
  lr_surf->write(of3);

  std::ofstream of6(outfile6.c_str());
  writePostscriptMesh(*lr_surf, of6);

  // Refine in the first parameter direction
  dir = XFIXED;
  double par_x1 = 1.5, par_x2 = 2.5;
  double start_x = 0.0;
  double end_x = 1.5;
  lr_surf->refine(dir, par_x1, start_x, end_x, mult, absolute);
  lr_surf->refine(dir, par_x2, start_x, end_x, mult, absolute);
  
  std::ofstream of4(outfile4.c_str());
  lr_surf->writeStandardHeader(of4);
  lr_surf->write(of4);

  std::ofstream of7(outfile7.c_str());
  writePostscriptMesh(*lr_surf, of7);
}

