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

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/LineCloud.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, n, density"  << std::endl;
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
  int n = atoi(argv[2]);
  double density = atof(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  model = factory.createFromIges(file1);
  
  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);


  std::vector<shared_ptr<GeneralMesh> > meshes;
  cvmodel->tesselate(&n, meshes);

  size_t ki;
  std::ofstream out1("lines1.g2");
  for (ki=0; ki<meshes.size(); ++ki)
    {
      double *nodes = meshes[ki]->vertexArray();
      int num_vx = meshes[ki]->numVertices();
      vector<double> seg;
      Point pt1(&nodes[0],&nodes[3]);
      Point pt2;
      for (int kj=1; kj<num_vx; kj++)
	{
	  pt2 = Point(&nodes[3*kj], &nodes[3*(kj+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  pt1 = pt2;
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out1);
      line_seg.write(out1);
    }
  
  std::vector<shared_ptr<GeneralMesh> > meshes2;
  cvmodel->tesselate(density, meshes2);

  std::ofstream out2("lines2.g2");
  for (ki=0; ki<meshes2.size(); ++ki)
    {
      double *nodes = meshes2[ki]->vertexArray();
      int num_vx = meshes2[ki]->numVertices();
      vector<double> seg;
      Point pt1(&nodes[0],&nodes[3]);
      Point pt2;
      for (int kj=1; kj<num_vx; kj++)
	{
	  pt2 = Point(&nodes[3*kj], &nodes[3*(kj+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  pt1 = pt2;
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out2);
      line_seg.write(out2);
     }
  
  int break_point;
  break_point = 1;
}

