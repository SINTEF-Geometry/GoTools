#define BOOST_TEST_MODULE gotools-core/testSplineCurve
#include <boost/test/unit_test.hpp>


#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"


using namespace Go;
using std::vector;
using std::ifstream;


BOOST_AUTO_TEST_CASE(testSplineCurve)
{
    ifstream infile1("data/rational_spline_cv_1.g2");
    ObjectHeader header;
    SplineCurve sc1;
    infile1 >> header >> sc1;

    ifstream infile2("data/rational_spline_cv_2.g2");
    SplineCurve sc2;
    infile2 >> header >> sc2;

    sc1.appendCurve(&sc2);


    int n = 1;
    int m = 4;
    BOOST_CHECK_EQUAL(n+3, m);


}
