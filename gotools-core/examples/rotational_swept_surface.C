//===========================================================================
//                                                                           
// File: rotational_swept_surface.C
//                                                                           
// Description:
//
// This program demonstrates the use of the static function
// 'rotationalSweptSurface' in the class 'SweepSurfaceCreator'.
// The function can generate a B-spline surface by rotating a curve around 
// an axis.
// This program reads a spline curve from file, prints it's bounding box and
// start and end points. The curve must be such that it doesn't lead to a
// self-intersecting surface.
// Then it creates a SplineSurface by calling 'rotationalSweptSurface' and
// prints it's bounding box, axis and the point on the axis.
// The spline curve's file name is hardcoded to 'approj_curve.g2', and the the
// axis direction and a point on the axis are also hardcoded in the program.
//
// Output is a file in Go-format for plotting the curve and the surface.
// The file name is hard-coded to "rotational_swept_surface.g2"
//
// See also example program 'surface_of_revolution.C'
//
//===========================================================================

#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include <fstream>

using std::string;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using namespace Go;
using std::shared_ptr;

int main(int argc, char** argv)
{
    const string inp_curve_filename("approj_curve.g2");

    cout << "\nRunning program '" << argv[0]
	 << "'\nSpline curve filename= '"
	 << inp_curve_filename.c_str() << "'." << endl;

    // Read spline curve file
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: "
	     << inp_curve_filename.c_str() << endl;
	return 1;
    }
    shared_ptr<SplineCurve> curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*curve);
    cfile.close();

    // Print some curve information
    Point pnt3d(3);
    curve->point(pnt3d, curve->startparam());     
    cout << "\nSplineCurve:  Dim= " << curve->dimension()
	 << "\nStart.  Param= " << curve->startparam() << "  Point= "
	 << pnt3d << endl;
    curve->point(pnt3d, curve->endparam());    
    cout << "End.  Param= " << curve->endparam() << "  Point= "
	 << pnt3d << endl;
    cout << "Bounding box =   " << curve->boundingBox() << endl;

    // Create a surface by rotating the curve around the axis an angle of 2PI.
    double angle = 2.0*M_PI;
    Point point_on_axis(0.0, 5.0, 200.0); 
    Point axis_dir(1.0, 0.0, 0.0);
    SplineSurface* surf =
	SweepSurfaceCreator::rotationalSweptSurface(*curve, angle,
						    point_on_axis, axis_dir);
    cout << "\nSurface:  Dim= " << surf->dimension() << endl;
    cout << "Bounding box =  "  << surf->boundingBox() << endl;
    cout << "Point on axis =  " << point_on_axis << endl;
    cout << "Axis direction = " << axis_dir << endl;

    // Open output  file
    ofstream fout("rotational_swept_surface.g2");
    // Write curve to file. Colour=red.
    fout << "100 1 0 4 255 0 0  255" << endl;
    curve->write(fout);

    // Write surface to file. Default colour=blue.    
    surf->writeStandardHeader(fout);
    surf->write(fout);

    // Write axis to file. Colour=green.
    double dlength = 1.2*(surf->boundingBox().high()[0] -
			  surf->boundingBox().low()[0]);
    Point endp = point_on_axis + dlength*axis_dir;
    SplineCurve* axis =  new SplineCurve(point_on_axis, endp);
    fout << "100 1 0 4 0 255 0  255" << endl;
    axis->write(fout);

    cout << "Open the file 'rotational_swept_surface.g2' in 'goview' to look"
	 << " at the results" << endl;

    delete surf;
    delete axis;

    return 0;
}
