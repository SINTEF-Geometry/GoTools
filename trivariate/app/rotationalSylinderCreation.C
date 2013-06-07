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

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

  // Test input parameters

  if (argc != 3)
    {
      cout << "Usage: " << argv[0] << " <parameters input filename> <final volume output filename>" << endl;
      exit(-1);
    }

  // Open input file
  ifstream infile(argv[1]);
  if (infile.bad())
    {
      cout << "Error when opening input file " << argv[1] << endl;
      exit(-1);
    }

  // Open volume output file
  ofstream outfile(argv[2]);
  if (outfile.bad())
    {
      cout << "Error when creating output file " << argv[1] << endl;
      exit(-1);
    }

  // The base point (intersection of bottom edge and rotation axis)
  double base_pt_x, base_pt_y, base_pt_z;
  infile >> base_pt_x >> base_pt_y >> base_pt_z;
  Point base_pt (base_pt_x, base_pt_y, base_pt_z);

  // The vector of the rotation axis direction, from bottom to top face of the solid
  // It does not have to be a unit, but it must be non-zero
  double rot_axis_x, rot_axis_y, rot_axis_z;
  infile >> rot_axis_x >> rot_axis_y >> rot_axis_z;
  Point rot_axis (rot_axis_x, rot_axis_y, rot_axis_z);
  if (rot_axis.length2() == 0.0)
    {
      cout << "Error: Vector describing the rotation axis is zero" << endl;
      exit(-1);
    }

  // Inner and outer radius, must be positive, but the first need not be smaller than
  // the other
  double radius_1, radius_2;
  infile >> radius_1 >> radius_2;
  if (radius_1 < 0.0)
    {
      cout << "Error: First radius is " << radius_1 << ", should be positive" << endl;
      exit(-1);
    }
  if (radius_2 < 0.0)
    {
      cout << "Error: Second radius is " << radius_2 << ", should be positive" << endl;
      exit(-1);
    }

  //Heigth, i.e. distance between bottom and top face. Should be positive
  double heigth;
  infile >> heigth;
  if (heigth < 0.0)
    {
      cout << "Error: Second radius is " << radius_2 << ", should be positive" << endl;
      exit(-1);
    }

  // First we create a vector normal to the rotation axis. Together with the
  // base point, it defines the line containing the bottom edge of the rectangle to be
  // swept when creating the solid

  // To create the vector, we pick one of the coordinate system base vectors, and take the
  // cross product with the axis.
  Point unit_vector;   // Either 1,0,0 or 0,1,0 or 0,0,1. Choose one to make sure the cross product is non-zero
  if (fabs (rot_axis_x) < fabs (rot_axis_y) && 
      fabs (rot_axis_x) < fabs (rot_axis_z))
    unit_vector = Point (1.0, 0.0, 0.0);
  else if (fabs (rot_axis_y) < fabs (rot_axis_z))
    unit_vector = Point (0.0, 1.0, 0.0);
  else
    unit_vector = Point (0.0, 0.0, 1.0);

  //Now create the vector normal to the the roation axis, and of unit length
  Point axis_normal = rot_axis % unit_vector;
  axis_normal.normalize();

  // Create the start and end point of the rectangle edge on the bottom face
  Point start_baseLine = base_pt + axis_normal * radius_1;
  Point end_baseLine = base_pt + axis_normal * radius_2;

  // Create the bottom line segment of the rectangle
  SplineCurve* baseLine = new SplineCurve(start_baseLine,
					  end_baseLine);

  // Then create the heigth of the rectangle, along the inner side of the solid
  SplineCurve* heigthLine = new SplineCurve(start_baseLine,
					    start_baseLine + (rot_axis * heigth / rot_axis.length()));

  // Now use a linear sweep to create the rectangle
  SplineSurface* verticalSection
    = SweepSurfaceCreator::linearSweptSurface(*baseLine,
					      *heigthLine,
					      start_baseLine);

  // Then rotate the rectangle to create the solid
  SplineVolume* cylinder
    = SweepVolumeCreator::rotationalSweptVolume(*verticalSection,
						2 * M_PI,
						base_pt,
						rot_axis);

  // Fianlly store the solid on the output file
  cylinder -> writeStandardHeader(outfile);
  cylinder -> write(outfile);

  return 0;
}
