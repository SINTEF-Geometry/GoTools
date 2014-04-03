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

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, nmb, output file"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  int nmb = atoi(argv[2]);
  std::ofstream file2(argv[3]);

  shared_ptr<ParamSurface> surf(new SplineSurface());
  ObjectHeader header;
  header.read(file1);
  surf->read(file1);

  shared_ptr<ftPointSet> triang = shared_ptr<ftPointSet>(new ftPointSet());
  vector<int> local_corner;
  RectDomain dom = surf->containingDomain();
  AdaptSurface::createTriangulation(surf, dom, triang, local_corner, false, nmb);

  vector<vector<int> > tri;
  vector<Vector3D> pos;

  // Get points and triangles
  triang->getPoints(pos);
  triang->getTriangles(tri);

  file2 << pos.size() << std::endl;
  file2 << tri.size() << std::endl;
  for (size_t ki=0; ki<pos.size(); ++ki)
    {
      pos[ki].write(file2);
      file2 << std::endl;
    }
  file2 << std::endl;
  for (size_t ki=0; ki<tri.size(); ++ki)
    {
      file2 << tri[ki][0] << " ";
      file2 << tri[ki][1] << " ";
      file2 << tri[ki][2] << " ";
      file2 << std::endl;
   }

}

