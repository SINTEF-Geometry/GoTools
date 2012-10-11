#define BOOST_TEST_MODULE ftEdgeTest
#include <boost/test/unit_test.hpp>

#include <fstream>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"


using namespace std;
using namespace Go;


BOOST_AUTO_TEST_CASE(ConstructFromCircle)
{
    ifstream in("../data/compositemodel/ftEdge.dat");
    BOOST_CHECK_MESSAGE(!in.bad(), "Input file not found or file corrupt");

    ObjectHeader header;
    shared_ptr<Circle> circle(new Circle);
    Point pt1(3);
    Point pt2(3);

    int index = 0;
    while (!in.eof()) {

	cout << "index: " << index << endl;

	// Circle
	header.read(in);
	circle->read(in);

	// Vertices
	pt1.read(in);
	shared_ptr<Vertex> v1(new Vertex(pt1));
	pt2.read(in);
	shared_ptr<Vertex> v2(new Vertex(pt2));

	// ftEdge
	ftEdge edge(circle, v1, v2);

	// Test vertices
	shared_ptr<Vertex> va, vb;
	edge.getVertices(va, vb);
	double eps = 1.0e-6;
	double dist1 = (v1->getVertexPoint() - va->getVertexPoint()).length();
	double dist2 = (v2->getVertexPoint() - vb->getVertexPoint()).length();
	bool verticesOK = (dist1 < eps) && (dist2 < eps);
        BOOST_CHECK(verticesOK);

	// Test if not reversed
	bool isNotReversed = !edge.isReversed();
	BOOST_CHECK(isNotReversed);

	// Test parameter interval
	double parlen = edge.tMax() - edge.tMin();
	bool lessthan2pi = (parlen <= 2.0 * M_PI);
	BOOST_CHECK(lessthan2pi);

	Utils::eatwhite(in);

	++index;
    }

}

