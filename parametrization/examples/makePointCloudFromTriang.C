//===========================================================================
//                                                                           
// File: makePointCloudFromTriang.C                                          
//                                                                           
// Created: Wed Mar 29 11:39:53 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: makePointCloudFromTriang.C,v 1.2 2006-03-29 11:03:07 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
    
    for (int i = 0; i < num_boundary_nodes; ++i) {
	const PrNode& tmp = sorted_boundary_nodes[i];
	os << tmp.x() << " " << tmp.y() << " " << tmp.z() << endl;
    }
};
