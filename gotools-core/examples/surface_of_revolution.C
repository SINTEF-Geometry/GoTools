//===========================================================================
//                                                                           
// File: surface_of_revolution.C
//                                                                           
// Description:
//  
// This program demonstrates the use of the class SurfaceOfRevolution.
//
// SurfaceOfRevolution is swept out by a SplineCurve that is rotated
// around an axis with a complete revolution, and is thereby a
// parametric surface. The space dimension is 3.
// The curve must be such that it doesn't lead to a self-intersecting surface.
//
// This program reads a spline curve from file, prints it's bounding box and
// start and end points. Then it creates a SurfaceOfRevolution object and
// prints it's bounding box and axis and location. 
//
// Input/Output
// The spline curve's file name is hardcoded to 'approj_curve.g2', and the
// location and the axis direction are also hardcoded in the program.
//
// Output is a file in Go-format for plot of the surface.
// The file name is hard-coded to 'surface_of_revolution.g2'
//   
//===========================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/ProjectCurve.h"
#include "GoTools/geometry/SurfaceOfRevolution.h"
#include <fstream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using namespace Go;
using std::shared_ptr;

int main(int argc, char** argv)
{
    string inp_curve_filename("approj_curve.g2"); 
    Point location(0.0, 5.0, 200.0);
    Point axis_dir(1.0, 0.0, 0.0);
    cout << "\nRunning program '" << argv[0] << "' with spline curve filename= '"
	 << inp_curve_filename.c_str() << "'." << endl;


    // Read spline curve file
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: " << inp_curve_filename.c_str()
	     << endl;
	return 1;
    }
    shared_ptr<SplineCurve> spline_curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*spline_curve);
    cfile.close();

    // Print some curve information
    Point pnt3d(3);
    spline_curve->point(pnt3d, spline_curve->startparam());     
    cout << "\nSplineCurve:  Dim= " << spline_curve->dimension()
	 << "\nStart.  Param= " << spline_curve->startparam() << "  Point= "
	 << pnt3d << endl;
    spline_curve->point(pnt3d, spline_curve->endparam());    
    cout << "End.  Param= " << spline_curve->endparam() << "  Point= "
	 << pnt3d << endl;
    cout << "Bounding box =   " << spline_curve->boundingBox() << endl;

    // Create the SurfaceOfRevolution object.
    SurfaceOfRevolution surf_of_revolution(location, axis_dir, spline_curve);
    cout << "\nSurface:  Dim= " << surf_of_revolution.dimension() << endl;
    cout << "Bounding box =  " << surf_of_revolution.boundingBox() << endl;
    cout << "Location =  " << surf_of_revolution.getLocation() << endl;
    cout << "Axis direction =  " << surf_of_revolution.getAxisDir() << endl;


    // Make a SplineSurface representation and write to file.
    SplineSurface* spline_surf = surf_of_revolution.geometrySurface();
    ofstream fout("surface_of_revolution.g2");
    spline_surf->writeStandardHeader(fout);
    spline_surf->write(fout);
    fout.close();
    cout << "\nOpen the files 'surface_of_revolution.g2' and 'approj_curve.g2'"
	 << " in 'goview' to look at the results.\n" << endl;
    delete spline_surf;

    return 0;
}
