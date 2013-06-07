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
