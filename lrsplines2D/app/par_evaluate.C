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
 
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/FileUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include <iostream>
#include <fstream>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Parameters : Input points (.txt), input surface (.g2), output points (.txt)"  << std::endl;
    exit(-1);
  }
  std::ifstream pointfile(argv[1]);
  std::ifstream surffile(argv[2]);
  std::ofstream pointout(argv[3]);

  // Read parameter values (and points)
  int del = 5;
  int nmb_pts = 0;
  vector<double> data;
  vector<double> extent(2*del);   // Limits for points in all coordinates
  FileUtils::readTxtPointFile(pointfile, del, data, nmb_pts, extent);

  // Read surface
  ObjectHeader header;
  header.read(surffile);
  shared_ptr<SplineSurface> surf(new SplineSurface());
  surf->read(surffile);
  
  // Evaluate
  for (int ki=0; ki<nmb_pts; ++ki)
    {
      Point pos = surf->ParamSurface::point(data[ki*del], data[ki*del+1]);
      pointout << data[ki*del] << " " << data[ki*del+1] << " " << pos << std::endl;
    }
}

