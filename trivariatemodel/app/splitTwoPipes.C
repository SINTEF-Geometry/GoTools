#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;


int main(int argc, char* argv[] )
{
  if (argc != 5)
      cout << "Usage: " << "<pipe data 1> <pipe data 2> <rational?> <output>" << endl;

  ifstream is1(argv[1]);
  ALWAYS_ERROR_IF(is1.bad(), "Bad or no input filename");

  ifstream is2(argv[2]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");

  int rat = atoi(argv[3]);
  double eps = 5.0e-3;

  std::ofstream of(argv[4]);

  // First pipe
  // The base point (intersection of bottom edge and rotation axis)
  double base_pt_x, base_pt_y, base_pt_z;
  is1 >> base_pt_x >> base_pt_y >> base_pt_z;
  Point base_pt (base_pt_x, base_pt_y, base_pt_z);

  // The vector of the rotation axis direction, from bottom to top face of the solid
  // It does not have to be a unit, but it must be non-zero
  double rot_axis_x, rot_axis_y, rot_axis_z;
  is1 >> rot_axis_x >> rot_axis_y >> rot_axis_z;
  Point rot_axis (rot_axis_x, rot_axis_y, rot_axis_z);
  if (rot_axis.length2() == 0.0)
    {
      cout << "Error: Vector describing the rotation axis is zero" << endl;
      exit(-1);
    }

  // Inner and outer radius, must be positive, but the first need not be smaller than
  // the other
  double radius_1, radius_2;
  is1 >> radius_1 >> radius_2;
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
  is1 >> heigth;
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

  shared_ptr<SplineVolume> cylinder1;
  if (rat)
    {
      // Then rotate the rectangle to create the solid
      cylinder1 =
	shared_ptr<SplineVolume>(SweepVolumeCreator::rotationalSweptVolume(*verticalSection,
									   2 * M_PI,
									   base_pt,
									   rot_axis));
    }
  else
    {
      int kstat = 0;
      SISLCurve *circle1 = NULL, *circle2 = NULL;
      s1303(start_baseLine.begin(), eps, 2*M_PI, base_pt.begin(), 
	    rot_axis.begin(), 3, &circle1, &kstat);
      s1303(end_baseLine.begin(), eps, 2*M_PI, base_pt.begin(), 
	    rot_axis.begin(), 3, &circle2, &kstat);
      shared_ptr<SplineCurve> rot1 = 
	shared_ptr<SplineCurve>(SISLCurve2Go(circle1));
      freeCurve(circle1);
      shared_ptr<SplineCurve> rot2 = 
	shared_ptr<SplineCurve>(SISLCurve2Go(circle2));
      freeCurve(circle2);

      vector<shared_ptr<SplineCurve> > rot_crvs(2); 
      rot_crvs[0] = 
	shared_ptr<SplineCurve>(rot1->subCurve(rot1->startparam(),
					       rot1->endparam()));
      rot_crvs[1] = 
	shared_ptr<SplineCurve>(rot2->subCurve(rot2->startparam(),
					       rot2->endparam()));
      unifyCurveSplineSpace(rot_crvs, 1.0e-3);

      vector<double> loft_par(2);
      loft_par[0] = 0.0;
      loft_par[1] = fabs(radius_2 - radius_1);

      shared_ptr<SplineSurface> seg =
	shared_ptr<SplineSurface>(CoonsPatchGen::loftSurface(rot_crvs.begin(),
							     loft_par.begin(), 2));
      Point end_pt = start_baseLine + (rot_axis * heigth / rot_axis.length());
      shared_ptr<SplineCurve> lin =
	shared_ptr<SplineCurve>(new SplineCurve(start_baseLine, 0.0, end_pt, heigth));
      cylinder1 =
	shared_ptr<SplineVolume>(SweepVolumeCreator::linearSweptVolume(*seg,
								       *lin,
								       start_baseLine));
    }

  // Second pipe
  // The base point (intersection of bottom edge and rotation axis)
  is2 >> base_pt_x >> base_pt_y >> base_pt_z;
  base_pt = Point(base_pt_x, base_pt_y, base_pt_z);

  // The vector of the rotation axis direction, from bottom to top face of the solid
  // It does not have to be a unit, but it must be non-zero
  is2 >> rot_axis_x >> rot_axis_y >> rot_axis_z;
  rot_axis = Point(rot_axis_x, rot_axis_y, rot_axis_z);
  if (rot_axis.length2() == 0.0)
    {
      cout << "Error: Vector describing the rotation axis is zero" << endl;
      exit(-1);
    }

  // Inner and outer radius, must be positive, but the first need not be smaller than
  // the other
  is2 >> radius_1 >> radius_2;
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
  is2 >> heigth;
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
  if (fabs (rot_axis_x) < fabs (rot_axis_y) && 
      fabs (rot_axis_x) < fabs (rot_axis_z))
    unit_vector = Point (1.0, 0.0, 0.0);
  else if (fabs (rot_axis_y) < fabs (rot_axis_z))
    unit_vector = Point (0.0, 1.0, 0.0);
  else
    unit_vector = Point (0.0, 0.0, 1.0);

  //Now create the vector normal to the the roation axis, and of unit length
  axis_normal = rot_axis % unit_vector;
  axis_normal.normalize();

  // Create the start and end point of the rectangle edge on the bottom face
  start_baseLine = base_pt + axis_normal * radius_1;
  end_baseLine = base_pt + axis_normal * radius_2;

  // Create the bottom line segment of the rectangle
  baseLine = new SplineCurve(start_baseLine,
			     end_baseLine);

  // Then create the heigth of the rectangle, along the inner side of the solid
  heigthLine = new SplineCurve(start_baseLine,
			       start_baseLine + (rot_axis * heigth / rot_axis.length()));

  // Now use a linear sweep to create the rectangle
  verticalSection
    = SweepSurfaceCreator::linearSweptSurface(*baseLine,
					      *heigthLine,
					      start_baseLine);

  // Then rotate the rectangle to create the solid
  shared_ptr<SplineVolume> cylinder2;
  if (rat)
    cylinder2 =
      shared_ptr<SplineVolume>(SweepVolumeCreator::rotationalSweptVolume(*verticalSection,
									 2 * M_PI,
									 base_pt,
									 rot_axis));
  else
    {
      int kstat = 0;
      SISLCurve *circle1 = NULL, *circle2 = NULL;
      s1303(start_baseLine.begin(), eps, 2*M_PI, base_pt.begin(), 
	    rot_axis.begin(), 3, &circle1, &kstat);
      s1303(end_baseLine.begin(), eps, 2*M_PI, base_pt.begin(), 
	    rot_axis.begin(), 3, &circle2, &kstat);
      shared_ptr<SplineCurve> rot1 = 
	shared_ptr<SplineCurve>(SISLCurve2Go(circle1));
      freeCurve(circle1);
      shared_ptr<SplineCurve> rot2 = 
	shared_ptr<SplineCurve>(SISLCurve2Go(circle2));
      freeCurve(circle2);

      vector<shared_ptr<SplineCurve> > rot_crvs(2); 
      rot_crvs[0] = 
	shared_ptr<SplineCurve>(rot1->subCurve(rot1->startparam(),
					       rot1->endparam()));
      rot_crvs[1] = 
	shared_ptr<SplineCurve>(rot2->subCurve(rot2->startparam(),
					       rot2->endparam()));
      unifyCurveSplineSpace(rot_crvs, 1.0e-3);

      vector<double> loft_par(2);
      loft_par[0] = 0.0;
      loft_par[1] = fabs(radius_2 - radius_1);

      shared_ptr<SplineSurface> seg =
	shared_ptr<SplineSurface>(CoonsPatchGen::loftSurface(rot_crvs.begin(),
							     loft_par.begin(), 2));
      Point end_pt = start_baseLine + (rot_axis * heigth / rot_axis.length());
      shared_ptr<SplineCurve> lin =
	shared_ptr<SplineCurve>(new SplineCurve(start_baseLine, 0.0, end_pt, heigth));
      cylinder2 =
	shared_ptr<SplineVolume>(SweepVolumeCreator::linearSweptVolume(*seg,
								       *lin,
								       start_baseLine));
    }
  double gap_eps = 1.0e-4;
  double kink_eps = 1.0e-2;

  shared_ptr<ftVolume> pipe1 = 
    shared_ptr<ftVolume>(new ftVolume(cylinder1, gap_eps, kink_eps));

  shared_ptr<ftVolume> pipe2 = 
    shared_ptr<ftVolume>(new ftVolume(cylinder2, gap_eps, kink_eps));

  pipe1->getVolume()->writeStandardHeader(of);
  pipe1->getVolume()->write(of);
  pipe2->getVolume()->writeStandardHeader(of);
  pipe2->getVolume()->write(of);

  std::ofstream of1("Pipe_1.g2");
  shared_ptr<SurfaceModel> model1 = pipe1->getOuterShell();
  int nmb1 = model1->nmbEntities();
  int ki;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamSurface> sf = model1->getSurface(ki);
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }

  std::ofstream of2("Pipe_2.g2");
  shared_ptr<SurfaceModel> model2 = pipe2->getOuterShell();
  int nmb2 = model2->nmbEntities();
  for (ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamSurface> sf = model2->getSurface(ki);
      sf->writeStandardHeader(of2);
      sf->write(of2);
    }


//   vector<shared_ptr<ftVolume> > vols = 
//     ftVolumeTools::splitVolumes(pipe1, pipe2, gap_eps);
//   std::cout << "Number of volumes: " << vols.size() << std::endl;

}
