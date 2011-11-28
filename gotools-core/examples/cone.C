//===========================================================================
//                                                                           
// File: cone.C                                                   
//                                                                           
// Description:
//
// This program demonstrates the use of the class Cone.
// It is a subclass of ElementarySurface.
// The space dimension of a cone is 3.
// The default parametrization is the angle u from 0 to 2*PI and the
// distance v from minus infinity to plus infinity.
//
// This program constructs a cone of infinite length and computes and
// prints the position, first derivatives and normal vector at a parameter
// pair u and v, and the closest point on the cone to a given point.
//
// Then it modifies the cone to have a length of four by bounding  the
// v parameters from -4 to 0, and the u parameters from 0 to 5 radians,
//
// Then it makes a SplineSurface representation of this cone and prints
// number of control points, spline orders, area and bounding box in
// addition to position, derivatives and closest point.
// Since the SplineSurface is not a perfect Cone, there is a slight difference
// between these values.
// Finally the SplineSurface cone is written to the file
// "spline_cone_segm.g2".
//
//===========================================================================


#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>

using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0] << "\nNo input arguments." << endl;

    // Construct a double cone with infinite length.
    // The constructor takes the radius, the location, the direction of the
    // (local) z-axis, the (approximate) direction of the (local) x-axis and
    // the cone angle in radians as input.
    // The z-axis will be the cone's axis. Location is the point where the
    // parameter in v-direction is zero, and this is where the radius of a
    // cross section perpendicular to the cone's axis is equal to the given
    // radius. The cone angle is the angle between the cone's axis and it's
    // side.
    // The use of the z- and x-axis in this constructor is based on the use of
    // the 'axis2_placement_3d' entity in STEP.

    double epsilon = 1.0e-10;
    double radius = 2.0;
    Point location(0.0, 0.0, 0.0);
    double cone_angle = 0.25*M_PI;
    Point x_axis(1.0, 0.0, 0.0);  // Direction vector
    Point z_axis(0.0, 0.0, 1.0);  // Direction vector
    double inp_angle = x_axis.angle(z_axis)*180/M_PI;

    cout << "\n*** Cone ***" << endl;
    // Make a cone. Unbounded in parameter v direction and local z-axis
    // direction.
    Cone cone(radius, location, z_axis, x_axis, cone_angle);
    // Get u and v parameter domain
    const RectDomain& rect_dom = cone.parameterDomain();
    double startparam_u = rect_dom.umin(); 
    double endparam_u   = rect_dom.umax();
    double startparam_v = rect_dom.vmin(); 
    double endparam_v   = rect_dom.vmax();

    // Get axes
    Point y_axis; // Has been computed in the constructor
    cone.getCoordinateAxes(x_axis, y_axis, z_axis);

    cout << "Location = " << cone.getLocation()
	 << "\nRadius = " << cone.getRadius()
	 << "\nx-axis = " << x_axis << "\ny-axis= " << y_axis
	 << "\nz-axis = " << z_axis
	 << "\ncone_angle = " << cone.getConeAngle() << endl;
    cout << "Angle between input x and z axes = " << inp_angle << " degrees.\n";
    cout << "Parameter u from "
	 << startparam_u << " to " << endparam_u << endl;
    cout << "Parameter v from "
	 << startparam_v << " to " << endparam_v << endl;
    cout << "Bounding box    = "; 
    cout << cone.boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    double upar = 1.5*M_PI;
    double vpar = -0.5;
    vector<Point> pts(3, Point(3));  // Three 3D points.
    Point normal;
    cone.point(pts, upar, vpar, 1);
    cone.normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;	  

    // Closest point().
    double clo_u, clo_v, clo_dist;
    Point pnt, clo_pt;
    pnt = Point(1.0, 1.0, -1.0);
    cone.closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nclosestPoint()"
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v << "   distance = "
	 << clo_dist << endl;


    // Make the cone bounded in the v parameter direction.
    startparam_v = -4.0;  // Cone of length four.
    endparam_v   =  0.0;
    // And in the u parameter direction.
    startparam_u = 0.0;
    endparam_u   = 5.0; 
    cone.setParameterBounds(startparam_u, startparam_v,
    		    endparam_u, endparam_v);
    cout << "\nThe cone is now made finite by bounding v-parameters from "
	 << startparam_v << "  to " << endparam_v << endl;
    cout << "\nu-parameters from "
	 << startparam_u << "  to " << endparam_u << endl; 
    cout << "Bounding box    = ";
    cout << cone.boundingBox() << endl;

    // Spline representation of a cone.
    SplineSurface* spline_cone = cone.geometrySurface();
    
    cout << "\n*** Cone Spline surface  ***" << endl;
    cout << "Parameter u from " << spline_cone->startparam_u() << " to "
	 << spline_cone->endparam_u() << endl;
    cout << "Parameter v from " << spline_cone->startparam_v() << " to "
	 << spline_cone->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_cone->numCoefs_u()
	 << "   v= " << spline_cone->numCoefs_v()<< endl;
    cout << "Order:  u= " << spline_cone->order_u()
	 << "   v= " << spline_cone->order_v()<< endl;

    cout << "Surface area    = " << spline_cone->area(epsilon) << endl;
    cout << "Bounding box    = " << spline_cone->boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    upar = 1.5*M_PI;
    vpar = -0.5;
    spline_cone->point(pts, upar, vpar, 1);
    spline_cone->normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;

    // Closest point().
    spline_cone->closestPoint(pnt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    cout << "\nclosestPoint()"
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nu = "    << clo_u << "   v = "  << clo_v << "   distance = "
	 << clo_dist << endl;

    // Write spline cone to file for plotting.
    ofstream fout("spline_cone.g2");
    spline_cone->writeStandardHeader(fout); // write header.
    fout << *spline_cone;    // write spline surface data.

    // Write tangent vectors and the normal vector.
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "410 1 0 4 255 0 0 255" << endl; // Header.
    fout << 4 << endl;  // Four lines
    fout << pts[0] << ' ' << pts[0] + pts[1] << endl;
    fout << pts[0] << ' ' << pts[0] + pts[2] << endl;
    // Write the normal vector
    fout << pts[0] << ' ' << pts[0] + normal << endl;    
    // Write a line from point to the closest point on cone.
    fout << pnt << ' ' << clo_pt << endl;

    // Write points
    fout << "400 1 0 4 0 255 0 255" << endl; // Header. Point cloud.
    fout << 3 << endl;       // Three points
    fout << pts[0] << endl;  // Point on surface
    fout << pnt << endl;     // Point
    fout << clo_pt << endl;  // It's closest point on the cone
    fout.close();

    // cout << "\nOpen the file 'spline_cone.g2' in 'goview' to look at the results"
    //      << endl;
    delete spline_cone;

    return 0;
}










