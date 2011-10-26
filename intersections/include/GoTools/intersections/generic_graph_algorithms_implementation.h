//===========================================================================
//                                                                           
// File: generic_graph_algorithms_implementation.h                           
//                                                                           
// Created: Mon Nov 29 12:23:31 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: generic_graph_algorithms_implementation.h,v 1.10 2005-09-06 13:28:12 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GENERIC_GRAPH_ALGORITHMS_IMPLEMENTATION_H
#define _GENERIC_GRAPH_ALGORITHMS_IMPLEMENTATION_H

#include "GoTools/utils/errormacros.h"
#include <vector>
#include <set>

//===========================================================================
// anonymous namespace 
namespace {
//===========================================================================

    enum NodeStatus {OUTSIDE_TREE, IN_TREE, IN_TREE_CYCLE_MEMBER};
    
    class EdgeSet : private std::set<std::pair<int, int> >
    {
    public:
	void insert(int i, int j) 
	{
	    least_first(i, j);
	    Base::insert(std::pair<int,int>(i, j));
	}

	bool inSet(int i, int j) const 
	{
	    least_first(i, j);
	    return Base::find(std::pair<int,int>(i, j)) != Base::end();
	}
    private:
	typedef std::set<std::pair<int, int> > Base;
	static inline void least_first(int& i, int& j) 
	{
	    ASSERT(i != j);
	    if (i > j) {
		std::swap(i, j);
	    }
	}
    };

    int find_next_root_node(const std::vector<NodeStatus>& node_status);

    void examine_node_generate_cycles(int parent_node_ix,
				      int node_ix, 
				      std::vector<NodeStatus>& node_status,
				      std::vector<int> prev_node_index,
				      const std::vector<std::vector<int> >& connection_table,
				      std::vector<std::vector<int> >& result);

    bool cycle_already_registered(int node_ix, int neigh_ix,
				  const std::vector<std::vector<int> >& registered_cycles);

    void make_cycle(int neigh_ix, int node_ix, 
		    const std::vector<int>& prev_node_index, 
		    std::vector<NodeStatus>& node_status,
		    std::vector<std::vector<int> >& result);

    template<class FunctorConnectedTo>
    void build_connectivity_table(int num_nodes, 
				  FunctorConnectedTo connectedTo,
				  std::vector<std::vector<int> >& table);

    bool is_connected_recursive(int node1_index, 
				int node2_index, 
				std::vector<bool>& visited,
				const std::vector<std::vector<int> >& connected);

    std::pair<int, int> find_next_path(const std::vector<int>& branch_point_indices,
				       const EdgeSet& eset,
				       const std::vector<std::vector<int> >& table);

    template <class Iter>
    Iter remove_if_withboolvec(Iter start, Iter end, const std::vector<bool>& flags)
    {
	Iter from = start;
	Iter to = start;
	for (; from != end; ++from) {
	    if (!flags[*from]) {
		*to = *from;
		++to;
	    }
	}
	return to;
    }

};

//===========================================================================
template<class FunctorConnectedTo>
void get_fundamental_cycle_set(int num_nodes, 
			       FunctorConnectedTo connectedTo, 
			       std::vector<std::vector<int> >& result)
//===========================================================================
{
    std::vector<std::vector<int> > table;
    build_connectivity_table(num_nodes, connectedTo, table);
    get_fundamental_cycle_set(num_nodes, table, result);
};

//===========================================================================
void get_fundamental_cycle_set(int num_nodes, 
			       const std::vector<std::vector<int> >& table,
			       std::vector<std::vector<int> >& result)
//===========================================================================
{
    //result.clear();
    std::vector<NodeStatus> node_status(num_nodes, OUTSIDE_TREE);
    std::vector<int> prev_node_index(num_nodes, -1);

     for (int i = 0; i != num_nodes; i = find_next_root_node(node_status)) {
	 node_status[i] = IN_TREE;
	 examine_node_generate_cycles(-1, i, node_status, prev_node_index, table, result);
    }
}


//===========================================================================
// Determining whether there is a path between two nodes in a graph
template<class FunctorConnectedTo>
bool is_path_connected(int node1_index, 
		       int node2_index, 
		       int num_nodes, 
		       FunctorConnectedTo connected_to)
//===========================================================================
{
    ASSERT(node1_index >= 0 && 
	   node1_index < num_nodes && 
	   node2_index >= 0 && 
	   node2_index < num_nodes);

    // shortcut if the nodes are equal or connected directly to each other
    // (prevents us from having to construct the connectivity table
    if (node1_index == node2_index || connected_to(node1_index, node2_index)) {
	return true;
    } 
    std::vector<std::vector<int> > table;
    build_connectivity_table(num_nodes, connected_to, table);
    std::vector<bool> visited(num_nodes, false);

    return is_connected_recursive(node1_index, node2_index, visited, table);
}

//===========================================================================
template<class FunctorConnectedTo>
void get_individual_paths(int num_nodes,
			  FunctorConnectedTo connectedTo,
			  std::vector<std::vector<int> >& paths,
			  std::vector<std::vector<int> >& cycles,
			  std::vector<int>& isolated_nodes)
//===========================================================================
{
    paths.clear();
    cycles.clear();
    isolated_nodes.clear();
    std::vector<std::vector<int> > table;
    build_connectivity_table(num_nodes, connectedTo, table);
    
//     // debug: print connectivity table
//     std::ofstream os("ConnectivityTable.data");
//     for (int i = 0; i < table.size(); ++i) {
// 	os << i << " connected to: " ;
// 	for (int j = 0; j < table[i].size(); ++j) {
// 	    os << table[i][j] << " ";
// 	}
// 	os << std::endl;
//     }


    std::vector<int> branch_point_indices;
    std::vector<bool> node_covered(num_nodes, false);
    branch_point_indices.reserve(num_nodes);

    // detecting branch points and isolated points
    for (int i = 0; i < num_nodes; ++i) {
	if (table[i].size() == 0) {
	    isolated_nodes.push_back(i);
	} else if (table[i].size() != 2) {
	    branch_point_indices.push_back(i); // this can be considered a branch point
	}
    }

    // register paths from branchpoints
    std::pair<int, int> path_start;
    EdgeSet exhausted_connections;
    path_start = find_next_path(branch_point_indices, exhausted_connections, table);
    int& prev_node = path_start.first;    // current branchpoint
    int& cur_node = path_start.second; // current "moving" node

    // unnesting paths
    while(prev_node != -1) {
	paths.insert(paths.end(), std::vector<int>());
	std::vector<int>& cur_path = paths.back();
	cur_path.insert(cur_path.end(), prev_node);
	cur_path.insert(cur_path.end(), cur_node);
	exhausted_connections.insert(prev_node, cur_node);
	node_covered[prev_node] = node_covered[cur_node] = true;
	
	while (table[cur_node].size() == 2) {
	    if (table[cur_node][0] != prev_node) {
		prev_node = cur_node;
		cur_node = table[cur_node][0];
	    } else {
		prev_node = cur_node;
		cur_node = table[cur_node][1];
	    }
	    cur_path.insert(cur_path.end(), cur_node);
	    node_covered[cur_node] = true;
	}
	// when we got here, we've arrived at another branch point
	// marking the end of the path as covered
	exhausted_connections.insert(prev_node, cur_node);

	path_start = find_next_path(branch_point_indices, exhausted_connections, table);
    }

    // Now we will have to register the remaining, closed loops.  
    // First, we will "unconnect" all covered nodes, then we call get_fundamental_cycle_set(...)
    for (int i = 0; i < int(table.size()); ++i) {
	// @@ Line below generated a compiler error on gcc 3.4.3. Replaced
	// with remove_if_withboolvec as a workaround.
// 	std::vector<int>::iterator new_end = 
// 	    std::remove_if(table[i].begin(), table[i].end(),
// 			   std::bind1st(std::mem_fun_ref(&std::vector<bool>::operator[]), node_covered));
	std::vector<int>::iterator new_end
	    = remove_if_withboolvec(table[i].begin(), table[i].end(), node_covered);
	if (new_end == table[i].begin()) {
	    table[i].clear();
	} else {
	    table[i].erase(new_end, table[i].end());
	}
    }

    get_fundamental_cycle_set(num_nodes, table, cycles);
}

//===========================================================================
namespace { // begin implementation of anonymous namespace
//===========================================================================

//===========================================================================
std::pair<int, int> find_next_path(const std::vector<int>& branch_point_indices,
				   const EdgeSet& eset,
				   const std::vector<std::vector<int> >& table)
//===========================================================================
{
    std::pair<int, int> result(-1, -1);
    int num_branch_points = (int)branch_point_indices.size();
    
    for (int i = 0; i < num_branch_points; ++i) {
	int branch_point = branch_point_indices[i];
	const std::vector<int>& connections = table[branch_point];
	for (int j = 0; j < int(connections.size()); ++j) {
	    int dir_point = connections[j];
	    if (!eset.inSet(branch_point, dir_point)) {
		// we found a branch point with an uncovered direction
		result.first = branch_point;
		result.second = dir_point;
		return result;
	    }
	}
	
    }
    return result; // (pair<int>(-1, -1) )
}


//===========================================================================
bool is_connected_recursive(int node1_index, 
			    int node2_index, 
			    std::vector<bool>& visited,
			    const std::vector<std::vector<int> >& connected)
//===========================================================================
{
    if (node1_index == node2_index) {
	return true;
    } 
    bool found = false;
    const std::vector<int>& node1_connections = connected[node1_index];
    for (int i = 0; i < int(node1_connections.size()) && !found; ++i) {
	int neigh = node1_connections[i];
	if (!visited[neigh]) {
	    visited[neigh] = true;
	    found = is_connected_recursive(node1_connections[i], node2_index, visited, connected);
	}
    }
    return found;
}

//===========================================================================
int find_next_root_node(const std::vector<NodeStatus>& node_status)
//===========================================================================
{
    int i;
    for (i = 0; i < int(node_status.size()); ++i) {
	if (node_status[i] == OUTSIDE_TREE) {
	    break;
	}
    }
    return i;
}


//===========================================================================
template<class FunctorConnectedTo>
void build_connectivity_table(int num_nodes, 
			      FunctorConnectedTo connectedTo,
			      std::vector<std::vector<int> >& table)
//===========================================================================
{
    table.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
	for (int j = i+1; j < num_nodes; ++j) {
	    if (connectedTo(i, j)) {
		table[i].insert(table[i].end(), j);
		table[j].insert(table[j].end(), i);
	    }
	}
    }
}


//===========================================================================
void examine_node_generate_cycles(int parent_node_ix,
				  int node_ix, 
				  std::vector<NodeStatus>& node_status,
				  std::vector<int> prev_node_index,
				  const std::vector<std::vector<int> >& connection_table,
				  std::vector<std::vector<int> >& result)
//===========================================================================
{
    const std::vector<int>& neighbours = connection_table[node_ix];
    for (int n = 0; n < int(neighbours.size()); ++n) {
	int neigh_ix = neighbours[n];
	if (neigh_ix != parent_node_ix) { // we disregard the connection back to the parent node
	    switch(node_status[neigh_ix]) {
	    case OUTSIDE_TREE:
		// add neighbour node to tree and examine further 
		node_status[neigh_ix] = IN_TREE;
		prev_node_index[neigh_ix] = node_ix;
		examine_node_generate_cycles(node_ix,
					     neigh_ix, 
					     node_status, 
					     prev_node_index, 
					     connection_table, 
					     result);
		break;
	    case IN_TREE:
		// we have detected new cycle
		make_cycle(neigh_ix, node_ix, prev_node_index, node_status, result);
		break;
	    case IN_TREE_CYCLE_MEMBER:
		if (!cycle_already_registered(node_ix, neigh_ix, result)) {
		    // we have detected a new cycle
		    make_cycle(neigh_ix, node_ix, prev_node_index, node_status, result);
		}
		break;
	    }
	}
    }
}

//===========================================================================
void make_cycle(int neigh_ix, 
		int node_ix, 
		const std::vector<int>& prev_node_index, 
		std::vector<NodeStatus>& node_status,
		std::vector<std::vector<int> >& result)
//===========================================================================
{
    result.insert(result.end(), std::vector<int>());
    std::vector<int>& new_cycle = result.back();
    new_cycle.insert(new_cycle.end(), neigh_ix);
    node_status[neigh_ix] = IN_TREE_CYCLE_MEMBER;

    for (int next_ix = node_ix; next_ix != neigh_ix; next_ix = prev_node_index[next_ix]) {
	ASSERT(next_ix != -1);
	new_cycle.insert(new_cycle.end(), next_ix);
	node_status[next_ix] = IN_TREE_CYCLE_MEMBER;
    }
}

//===========================================================================
bool cycle_already_registered(int node_ix, 
			      int neigh_ix,
			      const std::vector<std::vector<int> >& registered_cycles)
//===========================================================================
{
    for (int i = 0; i < int(registered_cycles.size()); ++i) {
	ASSERT(registered_cycles[i].size() >= 3);
	if (registered_cycles[i][0] == node_ix && registered_cycles[i][1] == neigh_ix) {
	    return true;
	}
    }
    return false;
}

}; // end implementation of anonymous namespace

#endif // _GENERIC_GRAPH_ALGORITHMS_IMPLEMENTATION_H



