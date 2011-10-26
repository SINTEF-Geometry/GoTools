//===========================================================================
//
// File : rotationalSylinderCreation.C
//
// Created: Wed Dec 10 10:53:58 2008
//
// Author: Kjell Fredrik Pettersen
//         SINTEF IKT
//
// Revision: $Id: rotationalSylinderCreation.C,v 1.1 2008-12-11 09:10:23 kfp Exp $
//
// Description: Program for creating a solid cylindrical tube with open ends
//
//===========================================================================



// The purpose of this program is to create a volume object representing a
// solid cylindrical tube with open ends (SCTOE), i.e. the object obtained as
// the set difference between two solid cylinders with the same rotation axis
// and the same end face planes, but with different radius.
//
// The shape and size of a SCTOE are described by three parameters:
//
//    - The heigth, i.e. distance between the two end faces
//    - The inner radius
//    - The outer radius
//
// In addition, the geometrical position and orientation of the SCTOE are given
// by two geometrical objects:
//
//    - The intersection point between the plane of the bottom face and the
//      rotation axis of the SCTOE, called the base point
//    - A vector describing the direction of rotation axis, from the bottom
//      face twoards the top face
//
// Any plane containing the rotation axis will intersect the SCTOE into two
// disjoint rectangles of the same shape and orientation, where the height
// is the heigth of the SCTOE and the widht is the difference between the inner
// and outer radius of the SCTOE. We fix such a rectangle R.
//
// The program uses two sweep methods when makeing the SCTOE.
//
// 1. First we make the rectangle R. This is given as a linear sweep of two
//    line segments beeing two non-opposite sides of R. Our choice of line
//    segments will be
//      a. S1, the side of R lying on the bottom face of the SCTOE as one side
//      b. S2, the side parallel to the rotation axis and with distance from the
//         rotation axis equal to the inner radius
//
// 2. Finally we create the SCTOE as a full rotational sweep of R around the
//    rotation axis
//
// NB! Notice that we could very well do the operations the other way: First
// a rotational sweep of S1 around the rotation axis to create the bottom face
// of the SCTOE. Then a linear sweep of the bottom face along S2.
//
// Input and output data:
//
// The program takes two arguments. One input file with the parameters describing
// the SCTOE and one output file for the SCTOE on g2-format.
//
// The data in the input files has this format:
//
//      base_point_x_coordinate,     base_point_y_coordinate,     base_point_z_coordinate,
//      rotation_axis_x_coordinate,  rotation_axis_y_coordinate,  rotation_axis_z_coordinate,
//      inner_radius,
//      outer_radius,
//      heigth



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
