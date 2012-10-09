#define BOOST_TEST_MODULE testSplineSurface
#include <boost/test/unit_test.hpp>

#include "GoTools/geometry/SplineSurface.h"


using namespace Go;


BOOST_AUTO_TEST_CASE(testSplineSurface)
{
    int n = 1;
    int m = 4;
    //BOOST_CHECK_EQUAL(n+1, m);
    BOOST_CHECK_EQUAL(n+3, m);

    SplineSurface surf;
}
