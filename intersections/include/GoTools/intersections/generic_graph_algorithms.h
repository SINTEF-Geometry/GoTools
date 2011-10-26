//===========================================================================
//                                                                           
// File: generic_graph_algorithms.h                                          
//                                                                           
// Created: Mon Nov 29 11:57:22 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: generic_graph_algorithms.h,v 1.5 2005-08-02 13:19:32 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GENERIC_GRAPH_ALGORITHMS_H
#define _GENERIC_GRAPH_ALGORITHMS_H

#include <vector>

//===========================================================================
//Paton's method for finding a fundamental set of cycles of an
//undirected graph.  Connectivity information between nodes is
//expressed by the functor 'connectedTo'.  NB: The outcome will be
//_appended_ to the 'result' sequence, no clear() operator will be
//called.
template<class FunctorConnectedTo > // FunctorConnectedTo(x, y)
                                    // returns 'true' if the node
                                    // indexed 'x' shares an edge with
                                    // the node indexed 'y'.
void get_fundamental_cycle_set(int num_nodes, // number of nodes in
					      // set, they will be
					      // indexed from 0
			       FunctorConnectedTo connectedTo, 
			       std::vector<std::vector<int> >& result); 
//===========================================================================

//===========================================================================
//Paton's method for finding a fundamental set of cycles of an
//undirected graph.  Connectivity information between nodes is
//expressed by the 'table' sequence.  Entry 'i' in this table is a
//list of indices of the nodes connected to node 'i'.  The sequence
//'result' is NOT emptied before use!
void get_fundamental_cycle_set(int num_nodes, 
			       const std::vector<std::vector<int> >& table,
			       std::vector<std::vector<int> >& result);
//===========================================================================

//===========================================================================
// Determining whether there is a path between two nodes in a graph
template<class FunctorConnectedTo> // Functor ConnectedTo(x,y) returns
                                   // 'true' if the node indexed 'x'
                                   // shares an edge with the node
                                   // indexed 'y'.
bool is_path_connected(int node1_index, 
		       int node2_index, 
		       int num_nodes, 
		       FunctorConnectedTo connected_to); 
//===========================================================================



//===========================================================================
//Will return all paths between branch-nodes, as well as all
//remaining, closed cycles.  Branch nodes are defined as those nodes
//incident with one, three or more edges.  (Nodes incident with 0
//edges are isolated, nodes incident with two edges lie on an
//"evident" path).  Closed cycles will be stored without duplication
//of the first node at the end.
template<class FunctorConnectedTo>
void get_individual_paths(int num_nodes,
			  FunctorConnectedTo connectedTo,
			  std::vector<std::vector<int> >& paths,
			  std::vector<std::vector<int> >& cycles,
			  std::vector<int>& isolated_nodes);
//===========================================================================



			       
#include "GoTools/intersections/generic_graph_algorithms_implementation.h"
#endif // _GENERIC_GRAPH_ALGORITHMS_H
