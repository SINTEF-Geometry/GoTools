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

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 14) {
      std::cout << "Input parameters : mid point (x,y,z), direction (x,y,z), major radius, minor radius, start par, angle (in degrees), axis (x,y,z)" << std::endl;
    exit(-1);
  }

#ifdef __BORLANDC__
  using Go::Point;
#endif

  Point midpt(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  Point dir(atof(argv[4]), atof(argv[5]), atof(argv[6]));
  double r1 = atof(argv[7]);
  double r2 = atof(argv[8]);
  double par1 = atof(argv[9]);
  double angle = atof(argv[10]);
  Point axis(atof(argv[11]), atof(argv[12]), atof(argv[13]));
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  angle = (angle*M_PI)/180.0;
  par1 = (par1*M_PI)/180.0;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createEllipticArc(midpt, dir, r1, r2, par1,
						    angle, axis);

  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  if (cvmodel)
  {
      std::ofstream out_file("ellipse.g2");
      int nmb = cvmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamCurve> crv = cvmodel->getCurve(ki);

	  crv->writeStandardHeader(out_file);
	  crv->write(out_file);
      }

      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "2" << std::endl;
      out_file << midpt << std::endl;
      out_file << midpt+dir*r1 << std::endl;


      std::ofstream out1("elliptic_mesh.g2");
      std::vector<shared_ptr<GeneralMesh> > meshes;
      int res[2];
      res[0] = res[1] = 100;
      cvmodel->tesselate(res, meshes);
      for (size_t kj=0; kj<meshes.size(); ++kj)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = meshes[kj]->numVertices();
	  double *vertices = meshes[kj]->vertexArray();
	  for (int kj=0; kj<nmb_vert-1; ++kj)
	    tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out1);
	  lines.write(out1);
	}
  }

  delete model;
}
