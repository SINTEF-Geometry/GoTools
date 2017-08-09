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

#define BOOST_TEST_MODULE ftEdgeTest
#include <boost/test/included/unit_test.hpp>

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
    ifstream in("data/ftEdge.g2");
    BOOST_CHECK_MESSAGE(!in.bad(), "Input file not found or file corrupt");

    ObjectHeader header;
    shared_ptr<Circle> circle(new Circle);
    Point pt1(3);
    Point pt2(3);

    int index = 0;
    while (!in.eof()) {

	//cout << "index: " << index << endl;

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

