#define BOOST_TEST_MODULE gotools-core/testSphere
#include <boost/test/unit_test.hpp>

#include <fstream>
#include "GoTools/geometry/Sphere.h"


using namespace std;
using namespace Go;


struct Config {
public:
    Config()
    {
        datadir = "data/"; // Relative to build/gotools-core

        //GoTools::init();
    }

public:
    string datadir;
    vector<string> infiles;
};


BOOST_FIXTURE_TEST_CASE(testSphere, Config)
{
    // A sphere
    double radius = 1.0;
    Point location(0.0, 0.0, 0.0);
    Point z_axis(0.0, 0.0, 1.0);
    Point x_axis(1.0, 0.0, 0.0);
    bool isSwapped = true;
    Sphere sphere(radius, location, z_axis, x_axis, isSwapped);

    // closestPoint()
    double clo_u, clo_v, clo_dist;
    Point clo_pt;
    double epsilon = 1.0e-6;
    Point pt_top(0.0, 0.0, 1.0);
    sphere.closestPoint(pt_top, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, 0.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);

    Point pt_bottom(0.0, 0.0, -1.0);
    sphere.closestPoint(pt_bottom, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, -0.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);

    Point pt1(0.0, -1.0, 0.0);
    sphere.closestPoint(pt1, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    BOOST_CHECK_CLOSE(clo_u, 0.0, epsilon);
    BOOST_CHECK_CLOSE(clo_v, 1.5*M_PI, epsilon);
    BOOST_CHECK_SMALL(clo_dist, epsilon);



    // A space curve
    Point normal(-1.0, 0.0, 0.0);
    Point circle_x(0.0, 0.0, 1.0);
    bool isReversed = true;
    Circle circle(radius, location, normal, circle_x, isReversed);
    double pi = 3.14159265358979; // M_PI slightly truncated
    circle.setParamBounds(pi, 2.0*M_PI);

    // getElementaryParamCurve()
    double tol = 1.0e-6;
    shared_ptr<ElementaryCurve> elem_cv 
        = sphere.getElementaryParamCurve(&circle, tol);

    BOOST_CHECK_MESSAGE(elem_cv, "Didn't get elementary curve.");

}

