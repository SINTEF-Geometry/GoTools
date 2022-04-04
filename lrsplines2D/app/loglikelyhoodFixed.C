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

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/lrsplines2D/LogLikelyhood.h"
#include <iostream>
#include <fstream>

using std::vector;
using namespace Go;

int main(int argc, char *argv[])
{

  if (argc != 3)
    {
      std::cout << "Parameters: residuals (x,y,z,r), degrees of freedom in T-distribution" << std::endl;
      return 1;
    }
  
  std::ifstream input(argv[1]);
  double Tny = atof(argv[2]);

  // Read point cloud with residuals
  vector<double> data;
  vector<double> extent(8);   // Limits for points in all coordinates
  int nmb_pts = 0;
  FileUtils::readTxtPointFile(input, 4, data, nmb_pts, extent);

  std::cout << "Number of points: " << nmb_pts << std::endl;
  // Extract residuals
  vector<double> residuals(nmb_pts);
  for (size_t ki=0; ki<nmb_pts; ++ki)
    residuals[ki] = data[4*ki+3];
  std::cout << "Min: " << extent[6] << ", max: " << extent[7] << std::endl;

  double loglh2 = 0.0;
  double loglh = LogLikelyhood::compute(residuals, Tny, false, loglh2);
  printf("Loglikelyhood: %7.13f, %7.13f \n",loglh, loglh2);
}
