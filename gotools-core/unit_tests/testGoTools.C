#define BOOST_TEST_MODULE testGoTools
#include <boost/test/unit_test.hpp>

#include "GoTools/geometry/GoTools.h"


using namespace Go;


BOOST_AUTO_TEST_CASE(testGoTools)
{
    GoTools::init();
    SplineSurface surf;
}
