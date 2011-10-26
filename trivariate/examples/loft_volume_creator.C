//===========================================================================
//                                                                           
// File: loft_volume_creator.C                                            
//                                                                           
//                                                                           
// Description:
//
// This program demonstrates the use of a static function 'loftVolume'
// in namespace 'LoftVolumeCreator'.
// The function use lofting to create a new 'SplineVolume' based on a set of
// surfaces.  The surfaces are not changed during the lofting process.
// The surfaces must lie in the same space.
//
// This program create a set of three SplineSurfaces as input to the function.
// The surfaces are rectangles perpendicular to the x-axis. The second rectangle
// is rotated and lifted.
// A vector is filled with shared pointers to the faces, and an iterator to the
// first surface and the number of surfaces are input arguments to the function.
// The function returns a pointer to the created lofting volume.
// 
// Output from this program is a file in Go-format. The file name is hard-coded to
// 'loft_volume_creator.g2'. The program 'goview' can't display volumes, but you
// can use the programs 'makeShield' or 'getBoundarySfs' to extract the boundary
// faces. They write a new file which can be used by 'goview'.
// Both programs have inputfilename and outputfilename as arguments, but
// 'makeShield' has a third optional parameter. If this third parameter is set
// to 0 (zero), and the file has more then one volume, the volumes will be
// displayed in different colours.
// 
//===========================================================================

#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/utils/Point.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include <memory>

using namespace Go;
using namespace std;
using std::shared_ptr;

int main(int argc, char* argv[] )
{
  cout << "\nRunning program " << argv[0]
       << ".\nNo input from user." << endl;

  Point x_axis(1.0, 0.0, 0.0);

  // Define a surface. Plane perpendicular to the x-axis.
  Point loc = Point(0.0, 0.0, 0.0);
  Point normal = x_axis;
  Plane plane(loc, normal);

  // Make a rectangle
  plane.setParameterBounds(0.0, 0.0, 2.0, 1.0); // (u1, v1, u2, v2)
  SplineSurface* surface1 = plane.geometrySurface(); // Spline representation

  // Place the rectangles along the x-axis.  
  Point trans_vec_x(5.0, 0.0, 0.0);
  SplineSurface* surface2 = surface1->clone();
  translateSplineSurf(trans_vec_x, *surface2);

  SplineSurface* surface3 = surface2->clone();
  translateSplineSurf(trans_vec_x, *surface3);

  // Rotate surface2 minus 30 degrees and lift it .
  rotateSplineSurf(x_axis, -M_PI/6.0, *surface2);
  Point trans_vec_z(0.0, 0.0, 2.0);
  translateSplineSurf(trans_vec_z, *surface2);

  // Create the lofting volume.
  vector<shared_ptr<SplineSurface> > surfs(3);
  surfs[0] = shared_ptr<SplineSurface>(surface1);
  surfs[1] = shared_ptr<SplineSurface>(surface2);
  surfs[2] = shared_ptr<SplineSurface>(surface3);
  SplineVolume* vol = Go::LoftVolumeCreator::loftVolume(surfs.begin(),
							(int)surfs.size());

  cout << "\nLoftVolume"
       << "\nBounding box = " << vol->boundingBox()
       << "\nParameter span =  " << vol->parameterSpan() << endl;
  Point startpnt(3), endpnt(3);
  vol->point(startpnt, vol->startparam(0), vol->startparam(1),
	     vol->startparam(2)); 
  vol->point(endpnt, vol->endparam(0), vol->endparam(1), vol->endparam(2)); 
  cout << "Point at parameter start = " <<startpnt << " and end = " << endpnt
       << endl; 
  cout << "Number of control points = u:" << vol->numCoefs(0) << "  v:"
       << vol->numCoefs(1) << "   w:" << vol->numCoefs(2) << endl;
  cout << "Order =  u: " << vol->order(0) << "  v: " << vol->order(0)
       <<  "  w: " << vol->order(0) << endl;
  cout << "Is rational? " << boolalpha << vol->rational() << endl;
  cout << "Is left handed? " << boolalpha << vol->isLeftHanded() << endl;

  // Open output volume file.
  ofstream os("loft_volume_creator.g2");

  // Write volume to file.
  vol->writeStandardHeader(os);
  vol->write(os);

  os.close();

  cout << "\nRun: makeShield loft_volume_creator.g2 loft_volume_creator_sf.g2\n"
       << "and open the file 'loft_volume_creator_sf.g2' in 'goview' to look at"
       << " the result.\n"
       << endl;

  // Clean up
  delete vol;

  return 0;
}

