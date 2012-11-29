#define BOOST_TEST_MODULE testSplineSurface
#include <boost/test/unit_test.hpp>

#include "GoTools/geometry/SplineSurface.h"


using namespace Go;
using std::vector;


BOOST_AUTO_TEST_CASE(testSplineSurface)
{
    int n = 1;
    int m = 4;
    //BOOST_CHECK_EQUAL(n+1, m);
    BOOST_CHECK_EQUAL(n+3, m);

    // Data from looped_surface.g2
    int dim = 3;
    int ncoefsu = 5;
    int ncoefsv = 2;
    int orderu = 4;
    int orderv = 2;
    double knotsu[] = { 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0 };
    double knotsv[] = { 0.0, 0.0, 2.0, 2.0 };
    double coefs[] = { 
        -1.0, -1.0, -1.0,
        0.5, -1.0, 0.5,
        0.0, -1.0, 2.0,
        -0.5, -1.0, 0.5,
        1.0, -1.0, -1.0,
        -1.0, 1.0, -1.0,
        0.5, 1.0, 0.5,
        0.0, 1.0, 2.0,
        -0.5, 1.0, 0.5,
        1.0, 1.0, -1.0
    };
    SplineSurface surf(ncoefsu, ncoefsv, orderu, orderv, knotsu, knotsv,
        coefs, dim);

    // Check knot vector
    vector<int> multu, multv;
    surf.basis_u().knotMultiplicities(multu);
    surf.basis_v().knotMultiplicities(multv);
    vector<double> knotvalsu, knotvalsv;
    surf.basis_u().knotsSimple(knotvalsu);
    surf.basis_v().knotsSimple(knotvalsv);
    BOOST_CHECK_EQUAL(multu.size(), 3);
    BOOST_CHECK_EQUAL(multu[0], 4);
    BOOST_CHECK_EQUAL(multu[1], 1);
    BOOST_CHECK_EQUAL(multu[2], 4);
    BOOST_CHECK_EQUAL(multv.size(), 2);
    BOOST_CHECK_EQUAL(multv[0], 2);
    BOOST_CHECK_EQUAL(multv[1], 2);
    BOOST_CHECK_EQUAL(knotvalsu.size(), 3);
    BOOST_CHECK_EQUAL(knotvalsu[0], 0.0);
    BOOST_CHECK_EQUAL(knotvalsu[1], 1.0);
    BOOST_CHECK_EQUAL(knotvalsu[2], 2.0);
    BOOST_CHECK_EQUAL(knotvalsv.size(), 2);
    BOOST_CHECK_EQUAL(knotvalsv[0], 0.0);
    BOOST_CHECK_EQUAL(knotvalsv[1], 2.0);

}
