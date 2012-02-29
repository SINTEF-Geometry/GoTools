
#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/utils/Point.h"

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: rotational_swept_volume.C
//                                                                           
/// Description:
///
/// This program demonstrates the use of the static function
/// 'rotationalSweptVolume' in the class 'SweepVolumeCreator'.
/// The function can generate a SplineVolume by rotating a surface around an axis.
/// The surface must be such that it doesn't lead to self-intersection.
///
/// This program creates a surface. The surface is part of a plane restricted in
/// x- and y-direction perpendicular to the z-axis . The rotational axis is
/// defined by an arbitrary point on the axis and the axis' direction vector.
/// Then it uses this surface and axis to create a SplineVolume .
///
/// Output is two files in Go-format. The file names is are hard-coded to
/// 'vol_rot_sweep_surf.g2' and 'rotational_swept_volume.g2'. The program 'goview'
/// can't display volumes, but you can use the programs 'makeShield' or
/// 'getBoundarySfs' to extract the boundary faces. They write a new file which
/// can be used by 'goview'. 
/// Both programs have inputfilename and outputfilename as arguments, but
/// 'makeShield' has a third optional parameter. If this third parameter is set
/// to 0 (zero) and the file has more then one volume, the volumes will be
/// displayed in different colours.
/// 
//===========================================================================

int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0]
	 << ".\nNo input from user." << endl;

    Point z_axis(0.0, 0.0, 1.0);

    // Define the surface. Plane perpendicular to the z-axis.
    Point loc = Point(0.0, 0.0, -2.0);
    Point normal = z_axis;
    Plane plane(loc, normal);
    plane.setParameterBounds(-4.0, -1.0, -2.0, 1.0); // (u1, v1, u2, v2)
    SplineSurface* surface = plane.geometrySurface(); // Spline representation
    cout << "Surface bounding box = " << surface->boundingBox() << endl;

   // Write surface to file. Colour=red.
    ofstream fout("vol_rot_sweep_surf.g2"); 
    fout << "200 1 0 4 255 0 0  255\n" << *surface << endl;
    fout.close();

    // Create a volume by rotating the surface around the axis an angle of 1.5PI.
    double angle = 1.5*M_PI;
    Point point_on_axis(-1.0, 0.0, -2.0); 
    Point axis_dir(0.0, 1.0, 0.0);
    SplineVolume* vol =
	SweepVolumeCreator::rotationalSweptVolume(*surface, angle,
						  point_on_axis, axis_dir);
    cout << "Rotational axis direction vector= " << axis_dir
	 << "\t Point on axis= " << point_on_axis << endl;

    cout << "\nRotationalSweptVolume"
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

    // Write volume to file.
    ofstream vout("rotational_swept_volume.g2");
    vol->writeStandardHeader(vout);
    vol->write(vout);

    // cout << "\nRun: makeShield rotational_swept_volume.g2 rotational_swept_"
    //      << "volume_sf.g2\nand open the files 'vol_rot_sweep_surf.g2' and "
    //      << "'rotational_swept_volume_sf.g2'\nin this order in 'goview' to look"
    //      << " at the results.\n" << endl;

    delete surface;
    delete vol;

    return 0;
}
