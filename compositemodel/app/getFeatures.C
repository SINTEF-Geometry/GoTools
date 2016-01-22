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
#include "GoTools/tesselator/LineStrip.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  if (sfmodel.get())
    {
      int nmb_bd = sfmodel->nmbBoundaries();
      int ki;
      std::ofstream out1("boundaries.g2");
      for (ki=0; ki<nmb_bd; ++ki)
	{
	  ftCurve bd = sfmodel->getBoundary(ki);
	  std::vector<shared_ptr<LineStrip> > bd_mesh;
	  bd.tesselate(20, bd_mesh);
	  for (size_t kj=0; kj<bd_mesh.size(); ++kj)
	    {
	      LineCloud lines;
	      std::vector<double> tmp_lines;
	      int nmb_vert = bd_mesh[kj]->numVertices();
	      double *vertices = bd_mesh[kj]->vertexArray();
	      for (int kr=0; kr<nmb_vert-1; ++kr)
		tmp_lines.insert(tmp_lines.end(), vertices+kr*3, 
				 vertices+(kr+2)*3);

	      lines.setCloud(&tmp_lines[0], nmb_vert-1);
	      lines.writeStandardHeader(out1);
	      lines.write(out1);
	    }
	}

      std::ofstream out2("g1Disc.g2");
      ftCurve g1disc = sfmodel->getG1Disconts();
      std::vector<shared_ptr<LineStrip> > g1d_mesh;
      g1disc.tesselate(20, g1d_mesh);
      for (size_t kj=0; kj<g1d_mesh.size(); ++kj)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = g1d_mesh[kj]->numVertices();
	  double *vertices = g1d_mesh[kj]->vertexArray();
	  for (int kr=0; kr<nmb_vert-1; ++kr)
	    tmp_lines.insert(tmp_lines.end(), vertices+kr*3, 
			     vertices+(kr+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out2);
	  lines.write(out2);
	}
    }
}

