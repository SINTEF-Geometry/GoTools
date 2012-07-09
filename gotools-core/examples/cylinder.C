
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: cylinder.C                                                   
//                                                                           
/// Description:
///
/// This program demonstrates the use of the class Cylinder.
/// It is a subclass of ElementarySurface.
/// The space dimension of a cylinder is 3.
/// The default parametrization is the angle u from 0 to 2*PI and the
/// distance v from minus infinity to plus infinity.
///
/// This program constructs a cylinder of infinite length and computes and
/// prints the position, first derivatives and normal vector at a parameter
/// pair u and v, and the closest point on the cylinder to a given point.
///
/// Then it modifies the cylinder to have a length of five by bounding  the
/// v parameters to -2.5 and 2.5.
/// Then it makes a SplineSurface representation of this cylinder and prints
/// number of control points, spline orders, area and bounding box.
/// Finally the SplineSurface cylinder is written to the file
/// "spline_cylinder_segm.g2".
///
//===========================================================================


int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0] << "\nNo input arguments." << endl;

    // Construct a cylinder with infinite length.
    // The constructor takes the radius, the location, the direction of the
    // (local) z-axis and the (approximate) direction of the (local) x-axis
    // as input. The z-axis will be the cylinder's axis. The use of the z- and
    // x-axis in this constructor is based on the use of the
    // 'axis2_placement_3d' entity in STEP.

    double epsilon = 1.0e-10;
    double radius = 1.0;
    Point location(0.0, 0.0, 0.0);
    const double sqrt2 = sqrt(2.0);
    Point x_axis(1.0, 1.0, -sqrt2);  // Direction vector
    Point z_axis(1.0, 1.0,  sqrt2);  // Direction vector
    double inp_angle = x_axis.angle(z_axis)*180/M_PI;

    cout << "\n*** Cylinder ***" << endl;
    // Make a cylinder. Unbounded in parameter v direction and local z-axis
    // direction.
    Cylinder cylinder(radius, location, z_axis, x_axis);
    // Get u and v parameter domain
    const RectDomain& rect_dom = cylinder.parameterDomain();
    double startparam_u = rect_dom.umin(); 
    double endparam_u   = rect_dom.umax();
    double startparam_v = rect_dom.vmin(); 
    double endparam_v   = rect_dom.vmax();

    // Get axes
    Point y_axis; // Have been computed in the constructor
    cylinder.getCoordinateAxes(x_axis, y_axis, z_axis);

    cout << "Location = " << location
	 << "\nRadius = " << radius
	 << "\nx-axis: " << x_axis << "\ny-axis= " << y_axis
	 << "\nz-axis: " << z_axis << endl;
    cout << "Angle between input x and z axes = " << inp_angle << " degrees.\n";
    cout << "Parameter u from "
	 << startparam_u << " to " << endparam_u << endl;
    cout << "Parameter v from "
	 << startparam_v << " to " << endparam_v << endl;

    // Bounding box of an infinite cylinder does not make sense.
    // It is not initialized, so trying to print it : 
    //   cout << cylinder.boundingBox() << endl;
    // will terminate the program by throwing an instance of 'std::exception'.

    // Evaluate position and first derivatives and normal vector at upar and vpar.
    double upar = 1.5*M_PI;
    double vpar = 0.0;
    vector<Point> pts(3, Point(3));  // Three 3D points.
    Point normal;
    cylinder.point(pts, upar, vpar, 1);
    cylinder.normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;	  

    // Closest point().
    double clo_u, clo_v, clo_dist;
    Point pnt, clo_pt;
    pnt = Point(0.0, 0.0, 2.0);
    cylinder.closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nclosestPoint()"
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v << "   distance = "
	 << clo_dist << endl;


    // Make the cylinder bounded in the v parameter direction.
    startparam_v = -2.5;  // Cylinder of length five.
    endparam_v   =  2.5; 
    cylinder.setParameterBounds(startparam_u, startparam_v,
				endparam_u, endparam_v);
    cout << "\nThe cylinder is now made finite by bounding v-parameters from "
	 << startparam_v << "  to " << endparam_v << endl;
    cout << "Bounding box    = ";
    cout << cylinder.boundingBox() << endl;

    // Spline representation of a cylinder.
    SplineSurface* spline_cylinder = cylinder.geometrySurface();
    
    cout << "\n*** Cylinder Spline surface  ***" << endl;
    cout << "Parameter u from " << spline_cylinder->startparam_u() << " to "
	 << spline_cylinder->endparam_u() << endl;
    cout << "Parameter v from " << spline_cylinder->startparam_v() << " to "
	 << spline_cylinder->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_cylinder->numCoefs_u()
	 << "   v= " << spline_cylinder->numCoefs_v()<< endl;
    cout << "Order:  u= " << spline_cylinder->order_u()
	 << "   v= " << spline_cylinder->order_v()<< endl;

    cout << "Surface area    = " << spline_cylinder->area(epsilon) << endl;
    cout << "Bounding box    = " << spline_cylinder->boundingBox() << endl;

    // Evaluate position and first derivatives and normal vector at upar and vpar.
    upar = 1.5*M_PI;
    vpar = 0.0;
    spline_cylinder->point(pts, upar, vpar, 1);
    spline_cylinder->normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;

    // Closest point().
    spline_cylinder->closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nclosestPoint()"
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v << "   distance = "
	 << clo_dist << endl;

    // Write spline cylinder to file for plotting.
    ofstream fout("spline_cylinder.g2");
    spline_cylinder->writeStandardHeader(fout); // write header.
    fout << *spline_cylinder;    // write spline surface data.

    // Write tangent vectors and the normal vector.
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "410 1 0 4 255 0 0 255" << endl; // Header.
    fout << 4 << endl;  // Four lines
    fout << pts[0] << ' ' << pts[0] + pts[1] << endl;
    fout << pts[0] << ' ' << pts[0] + pts[2] << endl;
    // Write the normal vector
    fout << pts[0] << ' ' << pts[0] + normal << endl;    
    // Write a line from point to the closest point on cylinder.
    fout << pnt << ' ' << clo_pt << endl;

    // Write points
    fout << "400 1 0 4 0 255 0 255" << endl; // Header. Point cloud.
    fout << 3 << endl;       // Three points
    fout << pts[0] << endl;  // Point on surface
    fout << pnt << endl;     // Point
    fout << clo_pt << endl;  // It's closest point on the cylinder
    fout.close();

    // cout << "\nOpen the file 'spline_cylinder.g2' in 'goview' to look at the results"
    //      << endl;
    delete spline_cylinder;

    return 0;
}










