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



#include <vector>
#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"



using namespace Go;
using namespace std;


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc != 4) {
    cout << "Usage:  " << argv[0] << " surfaceFile pointFile use_id_reg" << endl;
    return 1;
  }

  ifstream in_surf(argv[1]);
  ObjectHeader header;
  vector<shared_ptr<GeomObject> > surfaces;

  while (!in_surf.eof())
    {
      header.read(in_surf);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(in_surf);
      surfaces.push_back(obj);
      Utils::eatwhite(in_surf);
    }
  in_surf.close();

  ifstream in_pts(argv[2]);
  vector<float> pts;

  while (!in_pts.eof())
    {
      for (int j = 0; j < 3; ++j)
	{
	  float f;
	  in_pts >> f;
	  pts.push_back(f);
	}
      Utils::eatwhite(in_pts);
    }
  in_pts.close();

  bool use_id_reg = atoi(argv[3]) == 1;
  vector<vector<double> > regRotation;
  Point regTranslation;

  if (use_id_reg)
    {
      regTranslation = Point(0.0, 0.0, 0.0);
      regRotation.resize(3);
      for(int i = 0; i < 3; i++)
	{
	  regRotation[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    regRotation[i][j] = i==j;
	}
    }
  else
    {
      vector<Point> reg_pts_pc;
      vector<Point> reg_pts_sf;

      reg_pts_pc.push_back(Point(535.215, -516.438, -561.635));
      reg_pts_sf.push_back(Point(526.994, -461.239, -540.883));

      reg_pts_pc.push_back(Point(1419.07, -1378.61, -513.867));
      reg_pts_sf.push_back(Point(1533.7, -1447.84, -523.269));

      reg_pts_pc.push_back(Point(1899.35, 762.057, 345.142));
      reg_pts_sf.push_back(Point(1961.96, 775.782, 365.377));

      RegistrationInput regParameters;
      RegistrationResult regResult = registration(reg_pts_sf, reg_pts_pc, false, regParameters);
      /*
      RegistrationResult rescaleRegResult = registration(reg_pts_sf, reg_pts_pc, true, regParameters);
      cout << "Rescaling set to " << rescaleRegResult.rescaling_ << endl;
      */
      regRotation = regResult.rotation_matrix_;
      regTranslation = regResult.translation_;
    }

  shared_ptr<boxStructuring::BoundingBoxStructure> structure = preProcessClosestVectors(surfaces, 1000.0);

  // vector<float> distances = closestVectors(pts, structure, regRotation, regTranslation, 4, 70000, 10000000, 10000000);
  //vector<float> distances = closestVectors(pts, structure, regRotation, regTranslation, 4, 0, 100, 10000000);
  //distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 4, 0, 100, 10000000);

  vector<float> distances = closestVectors(pts, structure, regRotation, regTranslation, 4, 0, 10, 10000000);
  //  distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 4, 0, 100, 10000000);

  // vector<float> distances = closestVectors(pts, structure, regRotation, regTranslation, 4, 3400, 10000000, 100000000);
  // distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 4, 3400, 10000000, 100000000);
}
