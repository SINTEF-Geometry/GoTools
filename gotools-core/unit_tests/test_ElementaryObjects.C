/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#define BOOST_TEST_MODULE gotools-core/testElementaryObjects
#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"


using namespace std;
using namespace Go;


// Copied from app - should really unit test each object type separately


//int main(int argc, char** argv)
BOOST_AUTO_TEST_CASE(testElementaryObjects)
{
    double radius = 1.0;
    Point centre(0.0, 0.0, 0.0);
    Point z_axis(0.0, 0.0, 1.0);
    Point x_axis(1.0, 0.0, 0.0);
    Point normal = z_axis;
    Point location = centre;
    Point pt;
    double clo_t, clo_u, clo_v, clo_dist;
    Point clo_pt;
    double epsilon = 1.0e-10;

    // Plane
    //cout << endl << "*** Plane ***" << endl;
    Plane plane;

    // Line
    //cout << endl << "*** Line ***" << endl;
    Point loc(0.0, 0.0, 0.0);
    Point dir(1.0, 1.0, 1.0);
    Line line(loc, dir);
    line.setParamBounds(0.5, 1.5);
    line.point(pt, 0.5);
    //cout << pt << endl;
    BOOST_CHECK_EQUAL(pt, Point(0.5, 0.5, 0.5));

    line.point(pt, 1.5);
    //cout << pt << endl;
    BOOST_CHECK_EQUAL(pt, Point(1.5, 1.5, 1.5));

    line.reverseParameterDirection();
    line.point(pt, 0.5);
    //cout << pt << endl;
    BOOST_CHECK_EQUAL(pt, Point(1.5, 1.5, 1.5));

    line.point(pt, 1.5);
    //cout << pt << endl;
    BOOST_CHECK_EQUAL(pt, Point(0.5, 0.5, 0.5));

    loc = Point(1.0, 0.0, 0.0);
    dir = Point(0.0, 2.0, 0.0);
    line = Line(loc, dir);
    pt = Point(-1.0, -1.0, 0.0);
    double tmin = -5.0;
    double tmax = 5.0;
    line.closestPoint(pt, tmin, tmax, clo_t, clo_pt, clo_dist);
  //  cout << "closestPoint:" << endl
	 //<< clo_t << endl
	 //<< clo_pt << endl
	 //<< clo_dist << endl;
    BOOST_CHECK_CLOSE(clo_t, -0.5, epsilon);
    BOOST_CHECK_EQUAL(clo_pt, Point(1, -1, 0));
    BOOST_CHECK_CLOSE(clo_dist, 2, epsilon);

    // Circle
    //cout << endl << "*** Circle ***" << endl;
    Circle circle(radius, centre, normal, x_axis);
  //  cout << "Circle:" << endl
	 //<< circle << endl;
    double twopi = 2.0 * M_PI;
    double t0 = (1.0 / 12.0) * twopi;
    double t1 = (1.0 / 6.0) * twopi;
    Point pt0, pt1;
    circle.point(pt0, t0);
    circle.point(pt1, t1);
  //  cout << t0 << "\t" << pt0 << endl
	 //<< t1 << "\t" << pt1 << endl;
    BOOST_CHECK_SMALL(pt0.dist(Point(0.866025403784439, 0.5, 0)), epsilon);
    BOOST_CHECK_SMALL(pt1.dist(Point(0.5, 0.866025403784439, 0)), epsilon);

    // Check that the endpoints of a circle segment has correct
    // parametrization
    Circle* c0 = circle.subCurve(t0, t1);
    SplineCurve* sc0 = c0->geometryCurve();
    sc0->point(pt0, t0);
    sc0->point(pt1, t1);
  //  cout << t0 << "\t" << pt0 << endl
	 //<< t1 << "\t" << pt1 << endl;
    BOOST_CHECK_SMALL(pt0.dist(Point(0.866025403784439, 0.5, 0)), epsilon);
    BOOST_CHECK_SMALL(pt1.dist(Point(0.5, 0.866025403784439, 0)), epsilon);

    // The spline representation of a full circle does not have
    // arc-length parametrization...
    SplineCurve* sc1 = circle.geometryCurve();
    sc1->point(pt0, t0);
    sc1->point(pt1, t1);
  //  cout << t0 << "\t" << pt0 << endl
	 //<< t1 << "\t" << pt1 << endl;
    BOOST_CHECK_SMALL(pt0.dist(Point(0.872260419102717, 0.489041676410868, 0)), epsilon);
    BOOST_CHECK_SMALL(pt1.dist(Point(0.489041676410868, 0.872260419102717, 0)), epsilon);

    // closestPoint()
    circle = Circle(radius, centre, normal, x_axis);
    tmin = -0.25 * M_PI;
    tmax = 1.75 * M_PI;
    circle.setParamBounds(tmin, tmax);
    pt = Point(1.0, -1.0, 0.0);
    circle.closestPoint(pt, tmin, tmax, clo_t, clo_pt, clo_dist);
  //  cout << "closestPoint:" << endl
	 //<< clo_t << endl
	 //<< clo_pt << endl
	 //<< clo_dist << endl;

    double radius1180 = 1180.0;
    Point centre1180(-1320.0, -480.0, 0.0);
    Circle c1000(radius1180, centre1180, normal, x_axis);
    double t = 2.02;
    c1000.point(pt, t);
    c1000.closestPoint(pt, 0.0, 2.0*M_PI, clo_t, clo_pt, clo_dist);
    BOOST_CHECK_SMALL(clo_dist, 0.001);


    // reverseParameterDirection()
    tmin = 3.5;
    tmax = 4.0;
    circle.setParamBounds(tmin, tmax);
    //cout << "Before reversion:" << endl;
    circle.point(pt0, tmin);
    circle.point(pt1, tmax);
  //  cout << pt0 << endl
	 //<< pt1 << endl;
    circle.reverseParameterDirection();
    //cout << "After reversion:" << endl;
    circle.point(pt0, tmin);
    circle.point(pt1, tmax);
  //  cout << pt0 << endl
	 //<< pt1 << endl;

    // subCurve()
    circle = Circle(radius, centre, normal, x_axis);
    sc0 = circle.geometryCurve();
    tmin = 0.0;
    tmax = 0.33 * M_PI;
    sc1 = sc0->subCurve(tmin, tmax);
  //  cout << "Before cubCurve():" << endl
	 //<< *sc0 << endl;
    sc0->point(pt, tmax);
  //  cout << pt << endl;
  //  cout << "After cubCurve():" << endl
	 //<< *sc1 << endl;
    sc1->point(pt, tmax);
    //cout << pt << endl;

    // createSplineCurve()
    circle = Circle(radius, centre, normal, x_axis);
    tmin = 1.0;
    tmax = 2.0;
    circle.setParamBounds(tmin, tmax);
    //cout << "Before createSplineCurve():" << endl
    //    << circle << endl;
    sc0 = circle.createSplineCurve();
    //cout << "After createSplineCurve():" << endl
    //        << *sc0 << endl;
    sc0->equalBdWeights(false);
    //cout << "After equalBdWeights(false):" << endl
    //        << *sc0 << endl;


    // Cylinder
    //cout << endl << "*** Cylinder ***" << endl;
    radius = 1.0;
    location = Point(0.0, 0.0, 0.0);
    z_axis = Point(0.0, 0.0, 1.0);
    x_axis = Point(1.0, 0.0, 0.0);
    Cylinder cylinder(radius, location, z_axis, x_axis);
    cylinder.boundingBox();
    cylinder.point(pt, M_PI/2.0, 1.0);
    //cout << "(pi/2, 1) -> " << pt << endl;
    cylinder.point(pt, 1.0, M_PI/2.0);
    //cout << "(1, pi/2) -> " << pt << endl;
    // Swap parameters
    cylinder.swapParameterDirection();
    cylinder.point(pt, 1.0, M_PI/2.0);
    //cout << "swap(pi/2, 1) -> " << pt << endl;

    // Sphere
    //cout << endl << "*** Sphere ***" << endl;
    Sphere sphere(radius, centre, normal, x_axis);
    SplineSurface* sph = sphere.geometrySurface();
    //ofstream sphout("sphere_spline.g2");
    //sph->writeStandardHeader(sphout);
    //sph->write(sphout);

    // Cone
    //cout << endl << "*** Cone ***" << endl;
    radius = 6.25;
    z_axis = Point(0.0, -1.0, 0.0);
    x_axis = Point(-1.0, 0.0, 0.0);
    location = Point(0.0, 0.25, 0.0);
    double cone_angle = M_PI / 4.0;
    Cone cone(radius, location, z_axis, x_axis, cone_angle);
    pt = Point(-6.5, 0.0, 0.0);
    cone.closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
  //  cout << "closestPoint:" << endl
	 //<< clo_u << "\t" << clo_v << endl
	 //<< clo_pt << endl
	 //<< clo_dist << endl;

    cone.setParameterBounds(0.0, -20.0, 2.0*M_PI, 20.0);
    SplineSurface* scone = cone.geometrySurface();
    double upar = 1.0;
    double vpar = 2.0;
    cone.point(pt, upar, vpar);
  //  cout << "Cone:" << endl
	 //<< pt << endl;
    scone->point(pt, upar, vpar);
  //  cout << "SplineSurface cone:" << endl
	 //<< pt << endl;
    cone.setParameterBounds(upar, -20.0, 2.0*upar, 20.0);
    scone = cone.geometrySurface();
    scone->point(pt, upar, vpar);
  //  cout << "SplineSurface cone with setParameterDomain():" << endl
	 //<< pt << endl;


    // Torus
    //cout << endl << "*** Torus ***" << endl;
    double major_radius = 2.0;
    double minor_radius = 0.5;
    z_axis = Point(0.0, 0.0, 1.0);
    x_axis = Point(1.0, 0.0, 0.0);
    location = Point(0.0, 0.0, 0.0);
    Torus torus(major_radius, minor_radius, location, z_axis, x_axis);
    torus.setParameterBounds(0.25*M_PI, 0.0, 0.5*M_PI, 0.5*M_PI);
    SplineSurface* sstorus = torus.geometrySurface();
    //ofstream torout("torus_spline.g2");
    //sstorus->writeStandardHeader(torout);
    //sstorus->write(torout);

    // Disc
    //cout << endl << "*** Disc ***" << endl;
    centre = Point(0.0, 0.0, 0.0);
    radius = 10.0;
    x_axis = Point(1.0, 0.0, 0.0);
    normal = Point(0.0, 0.0, 1.0);
    Disc disc(centre, radius, x_axis, normal);
    disc.setParameterBounds(1.0, 0.0, 5.0, 0.5*M_PI);
    //cout << "Disc:" << endl
    //    << disc << endl;
    SplineSurface* sdisc = disc.asSplineSurface();

    //ofstream discout("disc_spline.g2");
    //sdisc->writeStandardHeader(discout);
    //sdisc->write(discout);
    //cout << "Disc as SplineSurface:" << endl
    //    << *sdisc << endl;

}
