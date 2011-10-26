//===========================================================================
//                                                                           
// File: torus.C                                                   
//                                                                           
// Description:
//
// This program demonstrates the use of the class Torus.
// It is a subclass of ElementarySurface.
// The space dimension of a torus is 3.
// The default parametrization is the angles u along the major circle and v
// along the minor circle from 0 to 2*PI.
//
// If the minor radius is greater than the major radius, the Torus is
// degenerate. In this case there are parameters you can use to select the 
// inner or outer (non-selfintersecting) part of the surface. See the program
// documentation for more information about this case.
//
// This program constructs a full torus and prints the values used to construct
// it.
// It computes and prints the position, first derivatives and normal vector at
// a parameter pair u and v, and the closest point on the torus to a given point.
//
// Then it modifies the torus by bounding the v parameters from PI/4 to 7/4PI,
// and the u parameters from 0 to 1.85*PI radians and prints the new bounding
// box.
//
// It makes a SplineSurface representation of this torus and prints
// number of control points, spline orders, area and an approximate bounding
// box in addition to position, derivatives and closest point.
// Since the SplineSurface is not a perfect Torus, there is a slight difference
// between these values.
// A minor and a major circle is constructed, and their parameter values and
// their bounding boxes are printed.
//
// Finally the SplineSurface torus is written to the file "spline_torus.g2".
// The minor and major circles, the computed tangent vectors and the normal
// vector, and the point and it's closest point are also written to this file.
// Open the file in 'goview' to look at the results.
//
//===========================================================================


#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0] << "\nNo input arguments." << endl;

    // Construct a full torus..
    // The constructor takes the major and minor radius, the location, the
    // direction of the (local) z-axis, the (approximate) direction of the
    // (local) x-axis. The local z-axis will be the torus' axis, and location
    // will be the torus' centre point.
    //
    // The use of the z- and x-axis in this constructor is based on the use of
    // the 'axis2_placement_3d' entity in STEP.

    double epsilon = 1.0e-10;
    double major_radius = 2.0;
    double minor_radius = 0.5;
    Point location(0.0, 0.0, 0.0);
    Point x_axis(1.0, 0.0, 0.0);  // Direction vector
    Point z_axis(0.0, 0.0, 1.0);  // Direction vector
    double inp_angle = x_axis.angle(z_axis)*180/M_PI;

    cout << "\n*** Torus ***" << endl;
    // Make a full torus.
    Torus torus(major_radius, minor_radius, location, z_axis, x_axis);
    // Get u and v parameter domain
    const RectDomain& rect_dom = torus.parameterDomain();
    double startparam_u = rect_dom.umin(); 
    double endparam_u   = rect_dom.umax();
    double startparam_v = rect_dom.vmin(); 
    double endparam_v   = rect_dom.vmax();

    // Get axes
    Point y_axis; // Has been computed in the constructor
    torus.getCoordinateAxes(x_axis, y_axis, z_axis);

    cout << "Location = " << torus.getLocation()
	 << "\nMajor radius = " << torus.getMajorRadius()
	 << "\nMinor radius = " << torus.getMinorRadius()
	 << "\nx-axis = " << x_axis << "\ny-axis= " << y_axis
	 << "\nz-axis = " << z_axis << endl;
    cout << "Angle between input x and z axes = " << inp_angle << " degrees.\n";
    cout << "Parameter u from "
	 << startparam_u << " to " << endparam_u << endl;
    cout << "Parameter v from "
	 << startparam_v << " to " << endparam_v << endl;
    cout << "Bounding box    = "; 
    cout << torus.boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    double upar = 1.5*M_PI;
    double vpar = M_PI;
    vector<Point> pts(3, Point(3));  // Three 3D points.
    Point normal;
    torus.point(pts, upar, vpar, 1);
    torus.normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;	  

    // Closest point().
    double clo_u, clo_v, clo_dist;
    Point pnt, clo_pt;
    pnt = Point(1.0, 1.0, -1.0);
    torus.closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nClosest point:  Point= " << pnt
	 << "\nClosest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v
	 << "   distance = " << clo_dist << endl;


    // Make the torus bounded.
    startparam_v = 0.25*M_PI;
    endparam_v   = 1.75*M_PI;
    startparam_u = 0.0;
    endparam_u   = 1.85*M_PI;
    torus.setParameterBounds(startparam_u, startparam_v,
			     endparam_u, endparam_v);
    cout << "\nThe torus is now made bounded by bounding v-parameters from "
	 << startparam_v << "  to " << endparam_v
	 << "\nu-parameters from "
	 << startparam_u << "  to " << endparam_u << endl; 
    cout << "Bounding box    = ";
    cout << torus.boundingBox() << endl;

    // Spline representation of a torus.
    SplineSurface* spline_torus = torus.geometrySurface();
    cout << "\n*** Torus Spline surface  ***" << endl;
    cout << "Parameter u from " << spline_torus->startparam_u() << " to "
	 << spline_torus->endparam_u()
	 << "\nParameter v from " << spline_torus->startparam_v() << " to "
	 << spline_torus->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_torus->numCoefs_u()
	 << "   v= " << spline_torus->numCoefs_v()<< endl;
    cout << "Order:  u= " << spline_torus->order_u()
	 << "   v= " << spline_torus->order_v()<< endl;
    cout << "Surface area    = " << spline_torus->area(epsilon) << endl;
    cout << "Bounding box    = " << spline_torus->boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    spline_torus->point(pts, upar, vpar, 1);
    spline_torus->normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;

    // Closest point().
    spline_torus->closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nClosest point:  Point= " << pnt
	 << "\nClosest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v
	 << "   distance = " << clo_dist << endl;

    // Get the minor and major circles for parameter PI/2.
    // If the other parameter is bounded, only a segment of a full circle is
    // returned.
    upar = M_PI/2;
    vpar = M_PI/2;
    const std::shared_ptr<Circle> minor_circle = torus.getMinorCircle(upar);
    const std::shared_ptr<Circle> major_circle = torus.getMajorCircle(vpar);

    cout << "\nMinor circle. u parameter = " << upar
	 << "  Circle parameter from " << minor_circle->startparam() << " to "
	 << minor_circle->endparam()
	 << "\nBounding box = " << minor_circle->boundingBox() << endl;

    cout << "\nMajor circle. v parameter = " << vpar
	 << "  Circle parameter from " << major_circle->startparam() << " to "
	 << major_circle->endparam()
	 << "\nBounding box = " << major_circle->boundingBox() << endl;

    // Write spline torus to file for plotting.
    ofstream fout("spline_torus.g2");
    spline_torus->writeStandardHeader(fout); // write header.
    fout << *spline_torus;    // write spline surface data.

    // Write tangent vectors and the normal vector.
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "410 1 0 4 255 0 0 255" << endl; // Header.
    fout << 4 << endl;  // Four lines
    fout << pts[0] << ' ' << pts[0] + pts[1] << endl;
    fout << pts[0] << ' ' << pts[0] + pts[2] << endl;
    // Write the normal vector
    fout << pts[0] << ' ' << pts[0] + normal << endl;    
    // Write a line from point to the closest point on the torus.
    fout << pnt << ' ' << clo_pt << endl;

    // Write points
    fout << "400 1 0 4 0 255 0 255" << endl; // Header. Point cloud.
    fout << 3 << endl;       // Three points
    fout << pts[0] << endl;  // Point on surface
    fout << pnt << endl;     // Point
    fout << clo_pt << endl;  // It's closest point on the torus

    // Write circles. Colour=green.
    SplineCurve* spline_minor_circle = minor_circle->geometryCurve();
    fout << "\n100 1 0 4 0 255 0  255" << endl;
    spline_minor_circle->write(fout);

    SplineCurve* spline_major_circle = major_circle->geometryCurve();
    fout << "100 1 0 4 0 255 0  255" << endl;
    spline_major_circle->write(fout);
    fout.close();

    cout << "\nOpen the file 'spline_torus.g2' in 'goview' to look at the results\n"
	 << endl;

    delete spline_torus;
    delete spline_minor_circle;
    delete spline_major_circle;

    return 0;
}
