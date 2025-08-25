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
#include "GoTools/lrsplines3D/LRSpline3DBezierCoefs.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Go;
using namespace std;

//===========================================================================
//                                                                           
/// Description:
/// Read LR spline volume from file.
/// Input to the example is a volume constructed in the example
/// program refine_lrvol.
/// Extract coefficients of Bezier patches.
/// Each patch is extended with knot vectors to be applicable for the g2-format
/// and written to file. The output can be visualized with
/// gotools/viewlib/app/goview_vol_and_lr.
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{
  // Output file
  std::string outfile("data/Bezier_patches.g2");
  std::ofstream of(outfile);
  
  // Read LR B-spline volume from file
  std::string infile("data/lrvol_fin.g2"); // Input file name
  std::cout << "\nAttempting to read LRSplineVolume from file: " << infile << std::endl;

  std::unique_ptr<LRSplineVolume> vol(new LRSplineVolume());

  // Specific file reading and volume creation.
  // In GoTools, you can use a factory or a specific reader method:
  std::ifstream input_stream(infile.c_str());
  if (!input_stream.is_open()) {
    std::cerr << "Error: Could not open file " << infile << std::endl;
    return 1;
  }
    
  // Read header specifying the type of geometry entity
  // The function throws if the entity header is invalid
  ObjectHeader header;
  try {
    header.read(input_stream);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  
  try {
    vol->read(input_stream);
  } catch (const std::exception& e) {
    std::cerr << "Error reading LRSplineVolume: " << e.what() << std::endl;
    return 1;
  }

  // Compute Bezier coefficients. The coefficients are stored in within bez
  LRSpline3DBezierCoefs bez(*vol);
  bez.getBezierCoefs();

  std::cout << "Finished converting volume " << std::endl;

  // Fetch information about Bezier patches
  int num = bez.numPatches();
  int dim = bez.dimension();  // Dimension of geometry space
  int deg1 = bez.degree(0);   // Polynomial degree in 1. parameter direction
  int deg2 = bez.degree(1);
  int deg3 = bez.degree(2);
  std::cout << "Number of patches: " << num << std::endl;
  std::cout << "Dimension of geometry space: " << dim << std::endl;
  std::cout << "Polynomial degrees: " << deg1 << ", " << deg2 << ", " << deg3 << std::endl;

  vector<double> coefs;   // Bezier coefficents for all patches. Each patch has
  // dim*(deg1+1)*(deg2+1)*(deg3+1) coeffients
  bez.BezierCoefficients(coefs);

  // For each Bezier patch, construct tensor-product spline volume and write to file
  // First define knot vectors
  vector<double> knots1(2*(deg1+1), 0.0);
  vector<double> knots2(2*(deg2+1), 0.0);
  vector<double> knots3(2*(deg3+1), 0.0);
  for (int ki=0; ki<=deg1; ++ki)
    knots1[deg1+1+ki] = 1.0;
  for (int ki=0; ki<=deg2; ++ki)
    knots2[deg2+1+ki] = 1.0;
  for (int ki=0; ki<=deg3; ++ki)
    knots3[deg3+1+ki] = 1.0;

  int n_coef = dim*(deg1+1)*(deg2+1)*(deg3+1);
  for (int ki=0; ki<num; ++ki)
    {
      // Define spline volume
      shared_ptr<SplineVolume> tpvol(new SplineVolume(deg1+1, deg2+1, deg3+1,
						      deg1+1, deg2+1, deg3+1,
						      &knots1[0], &knots2[0], &knots3[0],
						      &coefs[ki*n_coef], dim));

      // Write volume to file
      tpvol->writeStandardHeader(of);
      tpvol->write(of);
    }

  // One patch
  int ix = num/2;
  vector<double> coefs2;
  bez.BezierCoefficientsOfPatch(ix, coefs2);
  std::cout << "Bezier coefficients of patch number (starting count from zero): " << ix << ":" << std::endl;
  for (int ka=0; ka<(int)coefs2.size(); ka+=dim)
    {
      for (int kb=0; kb<dim; ++kb)
	std::cout << coefs2[ka+kb] << " ";
      std::cout << std::endl;
    }
 
}

