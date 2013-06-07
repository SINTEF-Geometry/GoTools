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


#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/IntResultsModel.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/tesselator/LineStrip.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::ofstream;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if (argc != 12)
    {
      std::cout << "Input arguments : Input file on IGES format, Input file on IGES format, ";
      std::cout << "x-value of first point on plane, y-value f first point, ... , z-value of third point" << std::endl;
      exit(-1);
    }


  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.1;
  double approx = 0.001;


  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");
  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  CompositeModel *model1 = factory.createFromIges(file1);
  CompositeModel *model2 = factory.createFromIges(file2);
  SurfaceModel* sm1 = dynamic_cast<SurfaceModel*>(model1);
  SurfaceModel* sm2 = dynamic_cast<SurfaceModel*>(model2);
  CompositeCurve* cm1 = dynamic_cast<CompositeCurve*>(model1);
  CompositeCurve* cm2 = dynamic_cast<CompositeCurve*>(model2);

  Point
    p_x (atof(argv[3]), atof(argv[4]), atof(argv[5])),
    p_y (atof(argv[6]), atof(argv[7]), atof(argv[8])),
    p_z (atof(argv[9]), atof(argv[10]), atof(argv[11]));

  ftPlane pl(p_x, p_y, p_z);

  ofstream osCrv ("dumpCurve.g2");
  ofstream osPl ("dumpPlane.g2");

  osPl << "200 1 0 0" << std::endl;
  osPl << "3 0" << std:: endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << p_x << std::endl;
  osPl << p_x << std::endl;
  osPl << p_y << std::endl;
  osPl << p_z << std::endl;

  if (sm1 && sm2)
    {
      int nmb_faces = (int)sm2 -> nmbEntities();

      for (int i = 0; i < nmb_faces; ++i)
	{
	  std::cout << "Appending " << i << " of " << nmb_faces << std::endl;
	  sm1 -> append(sm2 -> getFace(i));
	}

    }
  else if (cm1 && cm2)
    {
      int nmb_curves = (int)cm2->nmbEntities();
      for (int i = 0; i < nmb_curves; ++i)
	{
	  std::cout << "Appending " << i << " of " << nmb_curves << std::endl;
	  cm1 -> append(cm2 -> getCurve(i));
	}
      
      double start, end;
      cm1->parameterRange(start, end);
      std::cout << "Parameter range: " << start << ", " << end << std::endl;

      while (true)
	{
	  std::cout << "Get global parameter: ";
	  double par;
	  std::cin >> par;
	  if (par < start || par > end)
	    break;
	  double local;
	  int local_idx;
	  cm1->getLocalPar(par, local_idx, local);
	  std::cout << "Local parameter: " << local << " in curve ";
	  std::cout << local_idx << std::endl;

	  double global = cm1->getGlobalPar(local_idx, local);
	  std::cout << "Global par: " << global << std::endl;

	  Point pt1, pt2;
	  cm1->evaluate(local_idx, &local, pt1);
	  cm1->evaluateCurve(global, pt2);
	  std::cout << "Evaldist: " << pt1.dist(pt2) << std::endl;

	  shared_ptr<ParamCurve> curve = cm1->getCurve(local_idx);
	  int idx2 = cm1->getIndex(curve.get());
	  std::cout << "Index: " << idx2 << std::endl;
	}
    }
  double density = 1.0;
  std::vector<shared_ptr<LineStrip> > line_seg;
  PointCloud3D points;
  shared_ptr<IntResultsModel> results = model1->intersect_plane(pl);

  results->tesselate(density, line_seg, points);

  for (size_t ki=0; ki<line_seg.size(); ++ki)
    {
      LineCloud lines;
      vector<double> tmp_lines;
      int nmb_vert = line_seg[ki]->numVertices();
      double *vertices = line_seg[ki]->vertexArray();
      for (int kj=0; kj<nmb_vert-1; ++kj)
	tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

      lines.setCloud(&tmp_lines[0], nmb_vert-1);
      lines.writeStandardHeader(osCrv);
      lines.write(osCrv);
    }

  if (points.numPoints() > 0)
    {
      points.writeStandardHeader(osCrv);
      points.write(osCrv);
    }
}
