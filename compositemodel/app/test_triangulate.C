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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, IGES or g2 (1/0), density"  << std::endl;
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
  int useIGES = atoi(argv[2]);
  double density = atof(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  if (useIGES)
      model = factory.createFromIges(file1);
  else
      model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);


  shared_ptr<ftPointSet> triang;
  triang = sfmodel->triangulate(density);

  vector<vector<int> > tri;
  vector<Vector3D> pos;

  // Get points and triangles
  triang->getPoints(pos);
  triang->getOrientedTriangles(tri);

  std::ofstream out("triangles.raw");
  streamsize prev = out.precision(15);
  std::ofstream of("triangles1.g2");
  prev = of.precision(15);
   out << pos.size() << std::endl;
  out << tri.size() << std::endl;
  of << "400 1 0 4 0 255 0 255 " << std::endl;
  of << pos.size() << std::endl;
  for (size_t ki=0; ki<pos.size(); ++ki)
    {
      pos[ki].write(out);
      out << std::endl;
      pos[ki].write(of);
      of << std::endl;
    }
  out << std::endl;
  of << std::endl;
  of << "410 1 0 4 0 255 0 255 " << std::endl;
  of << 3*tri.size() << std::endl;
  for (size_t ki=0; ki<tri.size(); ++ki)
    {
      // if ((ki%2) == 0)
      // 	{
	  out << tri[ki][0] << " ";
	  out << tri[ki][1] << " ";
	  out << tri[ki][2] << " ";
	  out << std::endl;

	  pos[tri[ki][0]].write(of);
	  of << " ";
	  pos[tri[ki][1]].write(of);
	  of << std::endl;
	  pos[tri[ki][2]].write(of);
	  of << " ";
	  pos[tri[ki][1]].write(of);
	  of << std::endl;
	  pos[tri[ki][0]].write(of);
	  of << " ";
	  pos[tri[ki][2]].write(of);
	  of << std::endl;

      // 	}
      // else
      // 	{
      // 	  out << tri[ki][0] << " ";
      // 	  out << tri[ki][2] << " ";
      // 	  out << tri[ki][1] << " ";
      // 	  out << std::endl;
      // 	}

      // out << "410 1 0 0" << std::endl;
      // out << 3 << std::endl;
      // Vector3D x1 = pos[tri[ki][0]];
      // Vector3D x2 = pos[tri[ki][1]];
      // Vector3D x3 = pos[tri[ki][2]];

      // x1.write(out);
      // out << " ";
      // x2.write(out);
      // out << std::endl;
      // x2.write(out);
      // out << " ";
      // x3.write(out);
      // out << std::endl;
      // x3.write(out);
      // out << " ";
      // x1.write(out);
      // out << std::endl;
   }

}

