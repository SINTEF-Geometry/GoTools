#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: ellipse.C                                                   
///                                                                           
/// Description:
///
/// This program demonstrates the use of the class Ellipse.
/// It is a subclass of ElementaryCurve.
/// The space dimension of an ellipse is either 2 or 3.
/// The default parametrization is an angle from 0 to 2*PI.
///
/// This program constructs a full ellipse in the xy-plane and computes and
/// prints the ellipse's length, bounding box, and some points on the ellipse
/// at some parameter values. These points are written to the file
/// "spline_ellipse.g2".
/// The closest point on the ellipse to a given point gives an wrong answer, and
/// the parameter search area must be bounded to get the correct answer.
/// The given point and the closest point on the ellipse are written to the file
/// "spline_ellipse.g2" and are plotted with green colour.

/// Then it creates an ellipse segment (from 120 to 210 degrees) and computes
/// and prints the ellipse segments's length, bounding box and closest point.
///
/// The function reverseParameterDirection() is under construction and does not
/// create the correct curve yet, and the call to this function is commented out.
/// 
/// Then this program makes a SplineCurve representation of the ellipse segment.
/// It's endpoints will have the same parameter values as the ellipse's segment.
/// The SplineCurve ellipse segment is appended to the file "spline_ellipse.g2"
/// for plotting.
///
///
//===========================================================================


int main(int argc, char** argv)
{
    // Construct a full ellipse in the xy-plane.
    // The constructor takes the centre, the axis direction, the normal to
    // the plane of the ellipse and the lengths of the two semi-axis.
    double twopi = 2.0 * M_PI;
    double epsilon = 1.0e-10;
    double r1 = 2.0;
    double r2 = 1.0;
    Point centre(0.0, 0.0, 0.0);
    Point x_axis(1.0, 0.0, 0.0);  // Direction vector
    Point z_axis(0.0, 0.0, 1.0);  // Direction vector
    Point close_test_point(0.5, -0.5, 0.0);  // Point for test of closestPoint()

    Point plane_normal = z_axis;
    Point direction = x_axis;
 
    cout << "\n*** Full Ellipse ***" << endl;
    Ellipse ellipse(centre, direction, plane_normal, r1, r2);
    cout << "Start parameter = " << ellipse.startparam() << ".  "
	 << "End   parameter = " << ellipse.endparam() << endl;
    cout << "Ellipse length  = " << ellipse.length(epsilon) << endl;
    cout << "Bounding box    = " << ellipse.boundingBox() << endl;

    // Evaluate the ellipse's position at some parameter values t,
    // and write the points to a file for plotting.
    const int numpt = 13;
    ofstream fout("spline_ellipse.g2");  // Open output stream to file.
    // Class_PointCloud=400 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "400 1 0 4  255 0 0 255" << endl;   // Write header.
    fout << numpt << endl;

    cout << "\nEllipse position at some parameter values.\n";
    for (int i=0; i<numpt; ++i) {
	double t = i*twopi/(numpt-1.0);  // Parameter value
	Point pnt;                       // Ellipse position at parameter t.
	ellipse.point(pnt, t);
	cout << "Parameter=" << t << "\t Position= " << pnt << endl;
	fout << pnt << endl;  // Write input point coordinates to file.
    }

    // Closest point(). Full ellipse. Parameter search from 0 to 2*PI 
    double tmin, tmax, clo_t, clo_dist;
    Point pnt, clo_pt;
    pnt = close_test_point;
    tmin = 0.0;
    tmax = twopi;
    ellipse.closestPoint(pnt, tmin, tmax, clo_t, clo_pt, clo_dist);
    cout << "\nclosestPoint () between parameters " << tmin << " and " << tmax
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nt = "    << clo_t << "   distance = "     << clo_dist
	 << "\nThis is not really the closest point. Restrict the search area!"
	 << endl;

    // This was not really the closest point. Restrict the search area!
    // Closest point(). Parameter search from 1.1*PI to 2*PI 
    pnt = close_test_point;
    tmin = 1.1*M_PI;
    tmax = twopi;
    ellipse.closestPoint(pnt, tmin, tmax, clo_t, clo_pt, clo_dist);
    cout << "\nclosestPoint () between parameters " << tmin << " and " << tmax
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nt = "    << clo_t << "   distance = "     << clo_dist
	 << "\nThis is really the closest point."
	 << endl;

    // Write the closestPoint()'s points to the plotting file. Green colour.
    fout << "\n400 1 0 4  0 255 0 255" << endl;   // Write header.
    fout << "2\n" << close_test_point << '\n' << clo_pt << "\n\n";

    // Closest point(). Parameter search from 120 to 210 deg.
    // Parameter restrictions:
    //  -2*PI <= tmin < tmax <= 2*PI, and tmin-tmax <= 2*PI.
    tmin = 4*twopi/12;    tmax = 7*twopi/12;  // 120 to 210 deg.
    ellipse.closestPoint(pnt, tmin, tmax, clo_t, clo_pt, clo_dist);
    cout << "\nclosestPoint () between parameters " << tmin << " and " << tmax
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nt = "    << clo_t << "   distance = "     << clo_dist
	 << "\nThis closest point is the start point of the search area"
	 << endl;


    // Ellipse segment
    Ellipse ellipse_segm(ellipse);  // Default copy constructor
    //Set bounds for the parametrization of the ellipse.
    ellipse_segm.setParamBounds(tmin, tmax);
    Point start_point, end_point;
    ellipse_segm.point(start_point, tmin);    
    ellipse_segm.point(end_point, tmax);
    cout << "\n*** Ellipse segment ***" << endl;
    cout << "Start parameter = "  << ellipse_segm.startparam() << ".  "
	 << "Start point = "      << start_point << ".  " << endl
 	 << "End parameter   = "  << ellipse_segm.endparam() <<  ".  "
	 << "End point = "        << end_point << ".  " << endl;
    cout << "Ellipse length  = "  << ellipse_segm.length(epsilon) << endl;
    cout << "Bounding box    = "  << ellipse_segm.boundingBox() << endl;

    /*      reverseParameterDirection() is not correctly implemented yet!
    // Reverse parameter direction
    ellipse_segm.reverseParameterDirection();
    ellipse_segm.point(start_point, tmin);    
    ellipse_segm.point(end_point, tmax);  
    cout << "\n*** Ellipse segment after reversion ***" << endl;
    cout << "Start parameter = "  << ellipse_segm.startparam() << ".  "
	 << "Start point = "      << start_point << ".  " << endl
	 << "End parameter   = "  << ellipse_segm.endparam() <<  ".  "
    	 << "End point = "        << end_point << ".  " << endl;
    cout << "Ellipse length  = "  << ellipse_segm.length(epsilon) << endl;
    cout << "Bounding box    = "  << ellipse_segm.boundingBox() << endl;
    */

    // Spline representation of an ellipse.
    SplineCurve* spline_ellipse_segm = ellipse_segm.geometryCurve();
    spline_ellipse_segm->point(start_point, tmin);    
    spline_ellipse_segm->point(end_point, tmax);    
    cout << "\n*** Spline curve ellipse segment ***" << endl;
    cout << "Start parameter = " << spline_ellipse_segm->startparam() << ".  "
	 << "Start point = "     << start_point << ".  " << endl
	 << "End parameter   = " << spline_ellipse_segm->endparam() << ".  "
	 << "End point = "       << end_point << endl;
    cout << "Curve length    = " << spline_ellipse_segm->length(epsilon) << endl;
    cout << "Bounding box    = " << spline_ellipse_segm->boundingBox() << endl;

    // Write subcurve to file
    // Class_SplineCurve=100 MAJOR_VERSION=1 MINOR_VERSION=0 auxillary_data=0
    spline_ellipse_segm->writeStandardHeader(fout); // write header.
    fout << *spline_ellipse_segm;    // write spline curve data.
    fout.close();
    // cout << "\nOpen the file 'spline_ellipse.g2' in 'goview' to look at the results."
    //      << endl;
    delete spline_ellipse_segm;
}










