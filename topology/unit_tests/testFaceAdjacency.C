#define BOOST_TEST_MODULE testFaceAdjacency
#include <boost/test/unit_test.hpp>

#include <fstream>
#include "GoTools/topology/FaceAdjacency.h"
#include "GoTools/topology/tpEdge.h"
#include "GoTools/topology/tpFace.h"


using namespace std;
using namespace Go;


BOOST_AUTO_TEST_CASE(testFaceAdjacency)
{
    double tol_gap = 0.1;
    double tol_neighbour = 0.01;
    double tol_kink = 0.001;
    double tol_bend = 0.0001;
        
    tpTolerances tol(tol_gap, tol_neighbour, tol_kink, tol_bend);
    FaceAdjacency<tpEdge, tpFace> adjacency(tol);
        
}

