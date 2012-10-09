#define BOOST_TEST_MODULE module_testExample
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(first_pass)
{
    int n = 1;
    int m = 4;
    //BOOST_CHECK_EQUAL(n+1, m);
    BOOST_CHECK_EQUAL(n+3, m);
}


BOOST_AUTO_TEST_CASE(second_fail)
{
    int n = 1;
    int m = 4;
    //BOOST_CHECK_EQUAL(n+1, m);
    BOOST_CHECK_EQUAL(n+3, m);
}