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
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPoint.h"
#include <fstream>
#include <stdlib.h> // For atof()


const double EPS_GAP = 1e-7;
const double EPS_NEIGHBOUR = 1e-4;
const double EPS_KINK = 1e-2;
const double EPS_BEND = 4e-2;

const double EPS_APPROX = 1e-4;


using namespace std;
using namespace Go;


int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  if (argc != 8) {
    std::cout << "Input parameters : Input file on iges format, p0_x, p0_y, .., p1_x" << std::endl;
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


  ifstream file(argv[1]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromIges(file);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);
  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  Point p0(atof(argv[2]), atof(argv[3]), atof(argv[4]));
  Point p1(atof(argv[5]), atof(argv[6]), atof(argv[7]));

  Point dir = p1-p0;

  std::cout << "Starting search" << std::endl;

  if (sfmodel)
    {
      ftPoint result(0.0, 0.0, 0.0);
      bool hit = sfmodel->hit(p0, dir, result);
      std::cout << "Hit found " << hit << std::endl;

      std::ofstream out_hit("hit.g2");
      out_hit << "400 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << std::endl;
      out_hit << "410 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << " " << p1 << std::endl;
      if (hit)
	{
	  out_hit << "400 1 0 4 0 255 0 255" << std::endl;
	  out_hit << 1 << std::endl;
	  out_hit << result.position() << std::endl;
	}
    }
  else if (cvmodel)
    {
      PointOnCurve result;
      bool hit = cvmodel->hit(p0, dir, result);
      std::cout << "Hit found " << hit << std::endl;

      std::ofstream out_hit("hit.g2");
      out_hit << "400 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << std::endl;
      out_hit << "410 1 0 4 255 0 0 255" << std::endl;
      out_hit << "1" << std::endl;
      out_hit << p0 << " " << p1 << std::endl;
      if (hit)
	{
	  out_hit << "400 1 0 4 0 255 0 255" << std::endl;
	  out_hit << 1 << std::endl;
	  out_hit << result.getPos() << std::endl;
	}
    }
}
