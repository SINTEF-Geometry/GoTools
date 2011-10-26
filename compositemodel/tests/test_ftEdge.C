//==========================================================================
//                                                                          
// File: test_ftEdge.C                                            
//                                                                          
// Created: Wed Feb  3 12:23:16 2010                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id:$
//                                                                          
// Description:
//                                                                          
//==========================================================================


// Experimental use of Gtest

#if 0

#include <fstream>
#include <gtest/gtest.h>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"


using namespace std;
using namespace Go;


// int main(int argc, char** argv)
TEST(ftEdgeTest, ConstructFromCircle)
{
    ifstream in("data/ftEdge.dat");

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
	ASSERT_TRUE(verticesOK);

	// Test if not reversed
	bool isNotReversed = !edge.isReversed();
	ASSERT_TRUE(isNotReversed);

	// Test parameter interval
	double parlen = edge.tMax() - edge.tMin();
	bool lessthan2pi = (parlen <= 2.0 * M_PI);
	ASSERT_TRUE(lessthan2pi);

	eatwhite(in);

	++index;
    }

}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


#endif // if 0
