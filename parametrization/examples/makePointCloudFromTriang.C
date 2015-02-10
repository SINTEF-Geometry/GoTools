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

#include <fstream>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include "GoTools/parametrization/PrTriangulation_OP.h"

using namespace std;

int main(int varnum, char** vararg)
{
    if (varnum != 3) {
	cerr << "Usage: makePointCloudFromTriang <input file (triangulation)> <output file ( point cloud)> \n";
	cerr << "This program converts a triangulation readable for the PrTriangulation_OP class into \n";
	cerr << "a point cloud that is readable for the PrFastUnorganized_OP class." << endl;
	return -1;
    }

    ifstream is(vararg[1]);
    if (!is) {
	cerr << "Unable to open input file.  Aborting." << endl;
	return -1;
    }
    ofstream os(vararg[2]);
    if (!os) {
	cerr << "Unable to open output file.  Aborting." << endl;
	return -1;
    }

    PrTriangulation_OP tri;
    tri.scanRawData(is);

    int num_nodes = tri.getNumNodes();
    vector<PrNode> interior_nodes;
    map<int, PrNode> boundary_nodes;
    for (int i = 0; i < num_nodes; ++i) {
	if (tri.isBoundary(i)) {
	    boundary_nodes[i] = tri.getPrNode(i);
	}  else {
	    interior_nodes.push_back(tri.getPrNode(i));
	}
    }
    const int num_interior_nodes = (int)interior_nodes.size();
    const int num_boundary_nodes = (int)boundary_nodes.size();
    const int total_num_nodes = num_interior_nodes + num_boundary_nodes;
    assert(total_num_nodes == num_nodes);
    assert(num_boundary_nodes > 0);

    // sorting boundary nodes
    vector<PrNode> sorted_boundary_nodes;
    int first_bnode_ix = boundary_nodes.begin()->first;
    int cur_bnode = first_bnode_ix;
    std::vector<int> neigh;
    do {
	sorted_boundary_nodes.push_back(tri.getPrNode(cur_bnode));	
	tri.getNeighbours(cur_bnode, neigh);
	assert(neigh.size() >= 2);
	cur_bnode = neigh[0];
	
    } while (cur_bnode != first_bnode_ix);

    assert(sorted_boundary_nodes.size() == boundary_nodes.size());
    
    os << total_num_nodes << " " << num_interior_nodes << '\n';
    
    for (int i = 0; i < num_interior_nodes;++i) {
	const PrNode& tmp = interior_nodes[i];
	os << tmp.x() << " " << tmp.y() << " " << tmp.z() << endl;
    }
    os << endl;
    for (int i = 0; i < num_boundary_nodes; ++i) {
	const PrNode& tmp = sorted_boundary_nodes[i];
	os << tmp.x() << " " << tmp.y() << " " << tmp.z() << endl;
    }
};
