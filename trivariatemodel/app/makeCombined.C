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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 13)
      cout << "Usage: " << "<infile1> <point> <axis> <angle> <bd idx> <bd idx2> <length> <outfile>" << endl;

  ifstream infile1(argv[1]);
  ALWAYS_ERROR_IF(infile1.bad(), "Bad or no input filename");

  ObjectHeader header;
  header.read(infile1);

  shared_ptr<SplineVolume> vol1(new SplineVolume());
  vol1->read(infile1);

  Point pos(atof(argv[2]),atof(argv[3]),atof(argv[4]));
  Point axis(atof(argv[5]),atof(argv[6]),atof(argv[7]));
  double angle = atof(argv[8]);
  int bd = atoi(argv[9]);
  int bd2 = atoi(argv[10]);
  double len = atof(argv[11]);

  ofstream outfile(argv[12]);

  shared_ptr<SplineSurface> bd_sf = vol1->getBoundarySurface(bd);

  SweepVolumeCreator creator;
  shared_ptr<SplineVolume> vol2 = 
    shared_ptr<SplineVolume>(creator.rotationalSweptVolume(*bd_sf,
							   angle,
							   pos, 
							   axis));

  shared_ptr<SplineSurface> bd_sf2 = vol2->getBoundarySurface(bd2);
  
  double u1 = 0.5*(bd_sf2->startparam_u() + bd_sf2->endparam_u());
  double v1 = 0.5*(bd_sf2->startparam_v() + bd_sf2->endparam_v());
  vector<Point> res(3);
  bd_sf2->point(res, u1, v1, 1);
  Point norm = res[1].cross(res[2]);
  norm.normalize();
  Point pnt2 = res[0] + len*norm;
  shared_ptr<SplineCurve> crv
    = shared_ptr<SplineCurve>(new SplineCurve(res[0], pnt2));
  shared_ptr<SplineVolume> vol3 = 
    shared_ptr<SplineVolume>(creator.linearSweptVolume(*bd_sf2,
						       *crv,
						       res[0]));
  


  vol1->writeStandardHeader(outfile);  
  vol1->write(outfile);
  vol2->writeStandardHeader(outfile);  
  vol2->write(outfile);
  vol3->writeStandardHeader(outfile);  
  vol3->write(outfile);
}
