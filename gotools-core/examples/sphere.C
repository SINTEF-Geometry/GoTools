//===========================================================================
//                                                                           
// File: sphere.C                                                   
//                                                                           
// Description:
//
// This program demonstrates the use of the class Sphere.
// It is a subclass of ElementarySurface.
// The space dimension of a sphere is 3.
// The default parametrization is the angles u from 0 to 2*PI and v from
// -PI/2 to +PI/2.
//
// This program constructs a full sphere and prints the values used to
// construct it.
// It computes and prints the position, first derivatives and normal vector at
// a parameter pair u and v, and the closest point on the sphere to a given
// point.
//
// Then it modifies the sphere by bounding the u parameters from -PI/2 to PI/2
// and the v parameters from  -PI/3 to PI/3 radians and prints the new bounding
// box and gets the circle on the sphere at a constant longitude.
//
// It makes a SplineSurface representation of this sphere and prints
// number of control points, spline orders, area and an approximate bounding
// box in addition to position, derivatives and closest point.
// Since the SplineSurface is not a perfect Sphere, there is a slight difference
// between these values.
// Finally the SplineSurface sphere is written to the file "spline_sphere.g2".
// The longitude circle, the computed tangent vectors and the normal
// vector, and the point and it's closest point are also written to this file.
// Open the file in 'goview' to look at the results.
//
//===========================================================================

#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using namespace Go;
using std::shared_ptr;


int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0] << "\nNo input arguments." << endl;

    // Construct a full sphere.
    // The constructor takes the radius, the centre, the direction of the
    // (local) z-axis and the (approximate) direction of the (local) x-axis
    // as input. The use of the z- and x-axis in this constructor is based on
    // the use of the 'axis2_placement_3d' entity in STEP.

    double epsilon = 1.0e-10;
    double radius = 1.0;
    Point x_axis(1.0, 0.0, 0.0);   // Direction vector
    Point z_axis(0.0, 0.0, 1.0);   // Direction vector
    Point centre(2.0, 2.0, 2.0);
    double inp_angle = x_axis.angle(z_axis)*180/M_PI;

    cout << "\n*** Sphere ***" << endl;
    // Make a full sphere
    Sphere sphere(radius, centre, z_axis, x_axis);

    // Get u and v parameter domain
    const RectDomain& domain1 = sphere.parameterDomain();

    // Get axes
    Point y_axis; // Has been computed in the constructor
    sphere.getCoordinateAxes(x_axis, y_axis, z_axis);

    cout << "Centre = " << sphere.getLocation()
	 << "\nRadius = " << sphere.getRadius()
	 << "\nx-axis = " << x_axis << "\ny-axis= " << y_axis
	 << "\nz-axis = " << z_axis << endl;
    cout << "Angle between input x and z axes = " << inp_angle << " degrees.\n";
    cout << "Parameter u from " << domain1.umin() << " to "
	 << domain1.umax() << endl;
    cout << "Parameter v from " << domain1.vmin() << " to "
	 << domain1.vmax() << endl;
    cout << "Bounding box    = " << sphere.boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    double upar = M_PI/4;
    double vpar = 0.0;
    vector<Point> pts(3, Point(3));  // Three 3D points.
    Point normal;
    sphere.point(pts, upar, vpar, 1);
    sphere.normal(normal, upar, vpar);
    cout << "Point position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;


    // Closest point(). Search the whole sphere.
    double clo_u, clo_v; // Parameter values of closest point
    double clo_dist;     // Distance to closest point
    Point clo_pt1;       // Closest point on the sphere
    Point pnt(1.0, 2.0, 3.0); // Find closest point on the sphere to this point
    sphere.closestPoint(pnt, clo_u, clo_v, clo_pt1, clo_dist, epsilon);
    cout << "\nClosestPoint on whole sphere:  Point= " << pnt
	 << "\nClosest point= " << clo_pt1
	 << "\nu-parmeter = " << clo_u << "  v-parmeter = " << clo_v
	 << "   distance = "     << clo_dist << endl;


    // Closest point(). Parameter search from u = -PI/2 to PI/2, v from -PI/3
    // to PI/3
    Point clo_pt2;
    Array<double, 2> low(-M_PI/2, -M_PI/3);
    Array<double, 2> high(M_PI/2,  M_PI/3);
    RectDomain search_domain(low, high);
    sphere.closestPoint(pnt, clo_u, clo_v, clo_pt2, clo_dist, epsilon,
			 &search_domain);
    cout << "\nClosest Point. Search between u-parameters -PI/2 and PI/2"
	 << " and v-parameters -PI/3 and PI/3:\n"
	 << "Point= " << pnt
	 << "\nClosest point= " << clo_pt2
	 << "\nu-parmeter = " << clo_u << "  v-parmeter = " << clo_v
	 << "   distance = "     << clo_dist << endl;


    // Make the sphere bounded.
    sphere.setParameterBounds(-M_PI/2, -M_PI/3, M_PI/2, M_PI/3);
    const RectDomain domain = sphere.parameterDomain();
    cout << "\nThe sphere is now made bounded by bounding u-parameters from "
	 << domain.umin() << " to " << domain.umax()
	 << "\nv-parameters from " << domain.vmin() << " to " << domain.vmax()
	 << endl;
    cout << "Bounding box    = " << sphere.boundingBox() << endl;

    // Get a circle at a constant longitude(u-parameter).
    double longitude = 0;
    shared_ptr<Circle> longitudinal_circle = 
	sphere.getLongitudinalCircle(longitude);

    // Spline representation of the sphere.
    SplineSurface* spline_sphere = sphere.geometrySurface();
    cout << "\n*** Sphere spline surface ***" << endl;
    cout << "Parameter u from " << spline_sphere->startparam_u() << " to "
	 << spline_sphere->endparam_u() << endl;
    cout << "Parameter v from " << spline_sphere->startparam_v() << " to "
	 << spline_sphere->endparam_v() << endl;
    cout << "Number of control points:  u= " << spline_sphere->numCoefs_u()
	 << "   v= " << spline_sphere->numCoefs_v()<< endl;
    cout << "Order:  u= " << spline_sphere->order_u()
	 << "   v= " << spline_sphere->order_v()<< endl;
    cout << "Surface area    = " << spline_sphere->area(epsilon) << endl;
    cout << "Bounding box    = " << spline_sphere->boundingBox() << endl;

    // Evaluate position, first derivatives and normal vector at upar and vpar.
    spline_sphere->point(pts, upar, vpar, 1);
    spline_sphere->normal(normal, upar, vpar);
    cout << "\nPoint position at u=" << upar << " v=" << vpar << " : " << pts[0]
	 << "\nFirst derivatives : " << pts[1] << "  and  " << pts[2]
	 << "\nNormal vector     : " << normal << endl;

    // Closest point(). New point.
    pnt = Point(2.2, 1.8, 2.2); // Find closest point on the sphere to this point
    spline_sphere->closestPoint(pnt, clo_u, clo_v, clo_pt1, clo_dist, epsilon);
    cout << "\nClosestPoint on spline sphere segment:  Point= " << pnt
	 << "\nClosest point= " << clo_pt1
	 << "\nu-parmeter = " << clo_u << "  v-parmeter = " << clo_v
	 << "   distance = "     << clo_dist << endl;

    // Write spline surface to file
    ofstream fout("spline_sphere.g2");
    spline_sphere->writeStandardHeader(fout); // write header.
    spline_sphere->write(fout);  // write spline surface data.

    // Write tangent vectors and the normal vector.
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout << "410 1 0 4 255 0 0 255" << endl; // Header. LineCloud
    fout << 4 << endl;  // Four lines
    fout << pts[0] << ' ' << pts[0] + pts[1] << endl;
    fout << pts[0] << ' ' << pts[0] + pts[2] << endl;
    // Write the normal vector
    fout << pts[0] << ' ' << pts[0] + normal << endl; 
    // Write a line from point to the closest point on the sphere.
    fout << pnt << ' ' << clo_pt1 << endl;

    // Write longitudinal circle. Colour=red.
    SplineCurve* spline_long_circle = longitudinal_circle->geometryCurve();
    fout << "\n100 1 0 4 255 0 0  255" << endl; // Header. SplineCurve
    spline_long_circle->write(fout);

    // Write points
    fout << "400 1 0 4 0 255 0 255" << endl; // Header. Point cloud.
    fout << 3 << endl;       // Three points
    fout << pts[0] << endl;  // Point on surface
    fout << pnt << endl;     // Point
    fout << clo_pt1 << endl; // It's closest point on the sphere

    fout.close();
    // cout << "\nOpen the file 'spline_sphere.g2' in 'goview' to look at the results\n"
    //      << endl;

    delete spline_sphere;
    delete spline_long_circle;

    return 0;
}










