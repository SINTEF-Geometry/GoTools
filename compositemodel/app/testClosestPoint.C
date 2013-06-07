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

#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftMessage.h"
#include <istream>
#include <fstream>
#include <stdlib.h> // For atof()


using std::vector;
using std::istream;


using namespace Go;



int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if (argc < 5 || (argc % 3) != 2)
    {
      std::cout << "Input arguments : Input file on IGES format, ";
      std::cout << "p1X, p1Y, p1Z, ... , pnX, pxY, pnZ" << std::endl;

      exit(-1);
    }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");


  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;

  vector<double> coord;
  int nmb_points = (argc - 2) / 3;
  for (int i = 0; i < nmb_points * 3; ++i) coord.push_back (atof(argv[i+2]));

  vector<int> rgb;
  rgb.push_back(255);  rgb.push_back(0);  rgb.push_back(0);
  rgb.push_back(0);  rgb.push_back(255);  rgb.push_back(0);
  rgb.push_back(255);  rgb.push_back(255);  rgb.push_back(0);
  rgb.push_back(0);  rgb.push_back(255);  rgb.push_back(255);
  rgb.push_back(255);  rgb.push_back(128);  rgb.push_back(0);
  rgb.push_back(255);  rgb.push_back(64);  rgb.push_back(128);
  rgb.push_back(255);  rgb.push_back(0);  rgb.push_back(0);
  rgb.push_back(128);  rgb.push_back(255);  rgb.push_back(64);

  
  // Create the Surface model
  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromIges(file1);


  // Left out from fairingToolbox
  // Make clean degenerate surfaces and remove surfaces degenerated to a line
  // int nmb_deg=0, nmb_bd=0, nmb_removed=0;
  // status = tool->ensureCleanDegeneracy(nmb_deg, nmb_bd, nmb_removed);
  // std::cout << "Clean degeneracy. Status message : " << status.getMessage() << std::endl;
  // std::cout << "Number of degenerate triangular or banana surfaces: " << nmb_deg << std::endl;
  // std::cout << "Number of modified degenerate boundaries: " << nmb_bd << std::endl;
  // std::cout << "Number of removed surfaces: " << nmb_removed << std::endl;
  
  std::ofstream baseStr("data/basePoint.g2");

  for (int i = 0; i < nmb_points; ++i)
    {

      Point p(coord[i*3], coord[i*3+1], coord[i*3+2]);
      Point clp;
      int idx;
      double clo_p[2];
      double dist;

      model->closestPoint(p, clp, idx, clo_p, dist);
      if (nmb_points == 1)
	{
	  //ftSurface* closestSurface;
	  //Point norm;
	  //Point normCross;

	  //closestSurface = model->getSurface2(idx);
	  //norm = closestSurface -> normal(clo_p[0], clo_p[1]);
	  //normCross = norm % (p - clp);

	  std::cout << "Closest point is " << clp << std::endl;
	  std::cout << "Surface index is " << idx << std::endl;
	  std::cout << "Parameters on surface are (" << clo_p[0] << ", " << clo_p[1] << ")" << std::endl;
	  std::cout << "Distance is " << dist << std::endl;
	  //std::cout << "Norm cross length is " << normCross.length() << std::endl;
	}

      baseStr << "400 1 0 4 \n";
      int colPos = (i & 7) * 3;
      baseStr << rgb[colPos] << " " << rgb[colPos + 1] << " " << rgb[colPos + 2] << " 255\n";
      baseStr << "1\n";
      baseStr << p << "\n";
      baseStr << "400 1 0 4 \n";
      baseStr << (rgb[colPos]>>1) << " " << (rgb[colPos + 1]>>1) << " " << (rgb[colPos + 2]>>1) << " 255\n";
      baseStr << "1\n";
      baseStr << clp << "\n";
  
    }
}
