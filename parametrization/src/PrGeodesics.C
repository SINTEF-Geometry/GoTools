/********************************************************************
 FILENAME    : PrGeodesics.C
 AUTHOR      : Valerie PHAM-TRONG, SINTEF
 DATE        : Mai 2002
 DESCRIPTION : Geodesic path between two vertices
 CHANGE LOG  :
*********************************************************************/

#include <iomanip>
#include <string>

#include "GoTools/parametrization/PrGeodesics.h"
#include "GoTools/parametrization/PrPathTriangleSeq.h"
#include "GoTools/parametrization/PrInterface.h"

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"
#include <algorithm>

using namespace Go;
using namespace std;

//----------------------------------------------------------------------------
vector<Vector3D> GeodesicPT(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex)
//----------------------------------------------------------------------------
/*
 computes a geodesic path (ie locally shortest path)
 between the vertices indexed by sce_vertex and dest_vertex.
 iterativ method consisting in :
 1. initialisation of a path with the Dijkstra method which computes a shortest 
    path following the edges of the triangulation
 2. improves locally the current path and associated sequence
    with an update of the sequence around pivot vertices, ie
    vertices where the path is not locally optimal. 
 properties: - decreasing length at each step
             - finite number of steps
 returns the geodesic path as a polygonal line on the triangulated surface 
 */
{
    // cout << "\nComputation of a geodesic path";
    // cout << " between the vertex " << sce_vertex;
    // cout << " and the vertex " << dest_vertex << endl;

    // check that sce_vertex and dest_vertex are in the range of existing vertices
    while ((sce_vertex<0) ||(sce_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> sce_vertex;	
    }
    while ((dest_vertex<0) ||(dest_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> dest_vertex;	
    }

    vector<Vector3D> current_path;
    vector<int> tr_seq3d ;

    // lengths evolution during the iterations leading to the geodesic path 
    vector<double> lengths;
    double length;
	
    // variables for listing the pivots and updating around it in a sequence
    bool new_pivot;
		
    std::list<int> list_pivots_prec;
    std::list<int> list_pivots;
    std::list<int>::iterator iter_pivot;
    int pivot = -1;
    std::list<double> list_deviation;

    int nb_steps;
    nb_steps = 0;

    // path initialization (Dijkstra) between sce_vertex and dest_vertex    
    vector<int> vertex_path;
    vertex_path = DijkstraPT(triangulation, sce_vertex, dest_vertex);
    current_path = point3D_from_indice(triangulation, vertex_path);
    tr_seq3d = tr_crossed_by_path_vertices(triangulation, vertex_path);
 
    // Shortest path in a triangle sequence
    current_path = 
	shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					tr_seq3d, list_pivots, length);
 
    lengths.push_back(length);       
    if (list_pivots.empty())
    {
	new_pivot = 0;
    }
    else
    {
	new_pivot = 1;
	iter_pivot = list_pivots.begin();
    }
    while (new_pivot)
    {
	nb_steps ++;

	pivot = *iter_pivot;

	list_pivots_prec.clear();
	list_pivots_prec.insert(list_pivots_prec.begin(), list_pivots.begin(), list_pivots.end());
	// Update the sequence around the pivot if possible
	if (update_triangle_sequence_around_pivot(triangulation, tr_seq3d, pivot))
	{
	    // Shortest path in a triangle sequence
	    current_path = 
	    shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					    tr_seq3d, list_pivots, length);
	    lengths.push_back(length);       
	    //  pivot search
	    new_pivot = pivot_in_new_list(list_pivots, list_pivots_prec, (int) lengths.size()-1, lengths, iter_pivot, pivot);
	}
	else // no update, see next pivot
	{
	    new_pivot = next_pivot(list_pivots, iter_pivot, pivot);
	}		
    } // end while (new_pivot)
    
    // cout << "\nAfter " << lengths.size() << " iterations, the final path has ";
    // cout << list_pivots.size() << " pivots : " ;
    // print_list_int(list_pivots);    
    // cout << "its length is : " 
    // << length_polygonal_path(current_path) << endl;

    return(current_path);
}













//----------------------------------------------------------------------------
void GeodesicPT(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		std::vector<int> &vert_ind)
//----------------------------------------------------------------------------
/*
  see previous function, outputs different informations 
  to follow the procedure	       
 */
{
    while ((sce_vertex<0) ||(sce_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> sce_vertex;	
    }
    while ((dest_vertex<0) ||(dest_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> dest_vertex;	
    }

    //  vector<Vector3D> current_path;
    vector<int> tr_seq3d ;
    vert_ind.clear();

    // lengths evolution during the iterations leading to the geodesic path 
    vector<double> lengths;
    double length;
	
    // variables for listing the pivots and updating around it in a sequence
    bool new_pivot;
		
    std::list<int> list_pivots_prec;
    std::list<int> list_pivots;
    std::list<int>::iterator iter_pivot;
    int pivot = -1;
    std::list<double> list_deviation;

    std::vector<EdgeType> edge_seq;
    std::vector<double> ratios;
    int nb_steps;
    nb_steps = 0;

    // path initialization (Dijkstra) between sce_vertex and dest_vertex    
    vector<int> vertex_path;
    vertex_path = DijkstraPT(triangulation, sce_vertex, dest_vertex);
    //current_path = point3D_from_indice(triangulation, vertex_path);
    tr_seq3d = tr_crossed_by_path_vertices(triangulation, vertex_path);
 
    // Shortest path in a triangle sequence
    shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
				    tr_seq3d, list_pivots,
				    edge_seq, ratios, length);
    lengths.push_back(length);       
    if (list_pivots.empty())
    {
	new_pivot = 0;
    }
    else
    {
	new_pivot = 1;
	iter_pivot = list_pivots.begin();
    }
    while (new_pivot)
    {
	nb_steps ++;

	pivot = *iter_pivot;
	list_pivots_prec.clear();
	list_pivots_prec.insert(list_pivots_prec.begin(), list_pivots.begin(), list_pivots.end());
	// Update the sequence around the pivot if possible
	if (update_triangle_sequence_around_pivot(triangulation, tr_seq3d, 
						  pivot))
	{
	   shortest_path_triangle_sequence(triangulation, sce_vertex, 
					   dest_vertex,
					   tr_seq3d, list_pivots,
					   edge_seq, ratios, length);
	   //  pivot search
	}
	else // no update, see next pivot
	{
	    new_pivot = next_pivot(list_pivots, iter_pivot, pivot);
	}		
    } // end while (new_pivot)

    vert_ind.push_back(sce_vertex);
    for(size_t i=0; i<edge_seq.size(); i++)
    {
      int trinum=tr_seq3d[i];
      PrTriangle tri=triangulation.getPrTriangle(trinum);
      PrTriangle tri2=triangulation.getPrTriangle(tr_seq3d[i+1]);
      int nind1=edge_seq[i].vertex_[0];
      int nind2=edge_seq[i].vertex_[1];
      PrNode n1=triangulation.getPrNode(nind1);
      PrNode n2=triangulation.getPrNode(nind2);
      // ??? What on earth...?
      //      i++;
      //      i--;
    }
    vert_ind.push_back(dest_vertex);
    
    // cout << "\nAfter " << lengths.size() << " iterations, the final path has ";
    // cout << list_pivots.size() << " pivots : " ;
    // print_list_int(list_pivots);    
    // cout << "its length is : " 
    // << length_polygonal_path(current_path) << endl;

    //    return(current_path);
}





































//----------------------------------------------------------------------------
vector<Vector3D> GeodesicPT1(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex)
//----------------------------------------------------------------------------
/*
  see previous function, outputs more informations 
  to follow the procedure	       
 */
{
    cout << "\nComputation of a geodesic path";
    cout << " between the vertex " << sce_vertex;
    cout << " and the vertex " << dest_vertex << endl;

    // check that sce_vertex and dest_vertex are in the range of existing vertices
    while ((sce_vertex<0) ||(sce_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> sce_vertex;	
    }
    while ((dest_vertex<0) ||(dest_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> dest_vertex;	
    }

    vector<Vector3D> current_path;
    vector<int> tr_seq3d ;

    // lengths evolution during the iterations leading to the geodesic path 
    vector<double> lengths;
    double length;
	
    // variables for listing the pivots and updating around it in a sequence
    bool new_pivot;	
    std::list<int> list_pivots_prec;
    std::list<int> list_pivots;
    std::list<int>::iterator iter_pivot;
    int pivot = -1;
    std::list<double> list_deviation;
    int nb_steps;
    nb_steps = 0;

    // path initialization (Dijkstra) between sce_vertex and dest_vertex    
    vector<int> vertex_path;
    vertex_path = DijkstraPT(triangulation, sce_vertex, dest_vertex);
    cout << "Shortest path in the graph : " << endl;
    print(vertex_path);
    current_path = point3D_from_indice(triangulation, vertex_path);
    std::cout << "length of this initial path : " 
	 << length_polygonal_path(current_path) << endl << endl;
    tr_seq3d = tr_crossed_by_path_vertices(triangulation, vertex_path);
 
    std::ofstream os1("geom_sequence_init.off");
    printGeomviewTriangleSequence(os1, triangulation, tr_seq3d);
    std::ofstream os2("geom_path_init.vect");
    printGeomviewPath(os2, current_path);

    // Shortest path in a triangle sequence
    current_path = 
	shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					tr_seq3d, list_pivots, length);
 
    lengths.push_back(length);       
    if (list_pivots.empty())
    {
	new_pivot = 0;
    }
    else
    {
	new_pivot = 1;
	iter_pivot = list_pivots.begin();
    }
    while (new_pivot)
    {
	nb_steps ++;

	pivot = *iter_pivot;
	cout << "pivot : " << pivot;

	list_pivots_prec.clear();
	list_pivots_prec.insert(list_pivots_prec.begin(), list_pivots.begin(), list_pivots.end());
	
	// cout << "\nSequence of 3D Triangles: \n";
	// print(tr_seq3d);


	// Update the sequence around the pivot if possible
	if (update_triangle_sequence_around_pivot(triangulation, tr_seq3d, pivot))
	{
	    // Shortest path in a triangle sequence
	    current_path.clear();

	    current_path = 
	    shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					    tr_seq3d, list_pivots, length);
	    lengths.push_back(length); 
	    cout << endl << length << " ";      
//cout << "= long du + court ch. 2d calc par shortest_path_triangle_sequence : ";
//cout << endl << length_polygonal_path(current_path) ;      
//cout << " =long du ch 3d calc par shortest_path_triangle_sequence : ";

    //  pivot search
	    new_pivot = pivot_in_new_list(list_pivots, list_pivots_prec, (int) lengths.size()-1, lengths, iter_pivot, pivot);
	}
	else // no update, see next pivot
	{
	    new_pivot = next_pivot(list_pivots, iter_pivot, pivot);
	}		
    } // end while (new_pivot)
    
    std::ofstream os3("geom_sequence.off");
    printGeomviewTriangleSequence(os3, triangulation, tr_seq3d);
    std::ofstream os4("geom_path.vect");
    printGeomviewPath(os4, current_path);
 

    cout << "\nAfter " << lengths.size() << " iterations, the final path has ";
    cout << list_pivots.size() << " pivots : " ;
    print_list_int(list_pivots);    
    cout << "its length is : " 
	 << length_polygonal_path(current_path) << endl;

    //1// print_lengths(lengths);

    return(current_path);
}


//----------------------------------------------------------------------------
void GeodesicPT(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		vector<Vector3D>& current_path, vector<int>& tr_seq3d)
//----------------------------------------------------------------------------
/*
  see previous function
  returns the geodesic path as a polygonal line on the triangulated surface 
  and the triangle sequence crossed by this path
 */
{
    while ((sce_vertex<0) ||(sce_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> sce_vertex;	
    }
    while ((dest_vertex<0) ||(dest_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> dest_vertex;	
    }

    // lengths evolution during the iterations leading to the geodesic path 
    vector<double> lengths;
    double length;
    double length_last_reset;
	
    // variables for listing the pivots and updating around it in a sequence
    bool new_pivot;
		
    std::vector<int> vector_pivots_tried;
    std::list<int> list_pivots_prec;
    std::list<int> list_pivots;
    std::list<int>::iterator iter_pivot;
    int pivot = -1;
    std::list<double> list_deviation;

    int nb_steps;
    nb_steps = 0;

    // path initialization (Dijkstra) between sce_vertex and dest_vertex    
    vector<int> vertex_path;
    vertex_path = DijkstraPT(triangulation, sce_vertex, dest_vertex);
    current_path = point3D_from_indice(triangulation, vertex_path);
    tr_seq3d = tr_crossed_by_path_vertices(triangulation, vertex_path);
 
    current_path = 
	shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					tr_seq3d, list_pivots, length);
 
    length_last_reset=length;
    lengths.push_back(length);       
    if (list_pivots.empty())
    {
	new_pivot = 0;
    }
    else
    {
	new_pivot = 1;
	iter_pivot = list_pivots.begin();
    }
    while (new_pivot)
    {
	nb_steps ++;

	pivot = *iter_pivot;
	vector_pivots_tried.push_back(pivot);
	list_pivots_prec.clear();
	list_pivots_prec.insert(list_pivots_prec.begin(), list_pivots.begin(), list_pivots.end());
	
	// Update the sequence around the pivot if possible
	if (update_triangle_sequence_around_pivot(triangulation, tr_seq3d, pivot))
	{
	    // Shortest path in a triangle sequence
	    current_path = 
	    shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					    tr_seq3d, list_pivots, length);
	    lengths.push_back(length);       
	    //  pivot search
	   std::vector<int>::const_iterator tr_beg=vector_pivots_tried.begin();
	   std::vector<int>::const_iterator tr_end=vector_pivots_tried.end();
	   new_pivot = pivot_in_new_list(list_pivots, 
					 tr_beg, tr_end,
					 iter_pivot);
	   if (!new_pivot)
	   {
	      if (length<length_last_reset)
	      {
		 iter_pivot=list_pivots.begin();
		 if (iter_pivot!=list_pivots.end())
		 {
		    new_pivot=true;
		    length_last_reset=length;
		    vector_pivots_tried.clear();
		 }
	      }
	   }
	}	
	else // no update, see next pivot
	{
	    new_pivot = next_pivot(list_pivots, iter_pivot, pivot);
	}		
    } // end while (new_pivot)
    
    return;
}

//----------------------------------------------------------------------------
void GeodesicPT1(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		vector<Vector3D>& current_path, vector<int>& tr_seq3d)
//----------------------------------------------------------------------------
/*
  see previous function, outputs more informations 
  to follow the procedure	       
 */
{
    cout << "\nComputation of a geodesic path";
    cout << " between the vertex " << sce_vertex;
    cout << " and the vertex " << dest_vertex << endl;

    while ((sce_vertex<0) ||(sce_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> sce_vertex;	
    }
    while ((dest_vertex<0) ||(dest_vertex>(int)triangulation.getNumNodes()))
    {
	cout << "\nIndice of the source vertex should belong to 0 - ";
	cout << (int)triangulation.getNumNodes() << endl;
	cin >> dest_vertex;	
    }

    // lengths evolution during the iterations leading to the geodesic path 
    vector<double> lengths;
    double length;
	
    // variables for listing the pivots and updating around it in a sequence
    bool new_pivot;
		
    std::list<int> list_pivots_prec;
    std::list<int> list_pivots;
    std::list<int>::iterator iter_pivot;
    int pivot = -1;
    std::list<double> list_deviation;

    int nb_steps;
    nb_steps = 0;

    // path initialization (Dijkstra) between sce_vertex and dest_vertex    
    vector<int> vertex_path;
    vertex_path = DijkstraPT(triangulation, sce_vertex, dest_vertex);
    cout << "Shortest path in the graph : " << endl;
    print(vertex_path);
    current_path = point3D_from_indice(triangulation, vertex_path);
    cout << "length of this initial path : " ;
    cout << length_polygonal_path(current_path) << endl << endl;
    tr_seq3d = tr_crossed_by_path_vertices(triangulation, vertex_path);
 
    std::ofstream os1("geom_sequence_init.off");
    printGeomviewTriangleSequence(os1, triangulation, tr_seq3d);
    std::ofstream os2("geom_path_init.vect");
    printGeomviewPath(os2, current_path);

    // Shortest path in a triangle sequence
    current_path = 
	shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					tr_seq3d, list_pivots, length);
 
    lengths.push_back(length);       
    if (list_pivots.empty())
    {
	new_pivot = 0;
    }
    else
    {
	new_pivot = 1;
	iter_pivot = list_pivots.begin();
    }
    while (new_pivot)
    {
	nb_steps ++;

	pivot = *iter_pivot;
	//1// cout << "Treated Pivot : " << pivot;

	list_pivots_prec.clear();
	list_pivots_prec.insert(list_pivots_prec.begin(), list_pivots.begin(), list_pivots.end());
	
	// cout << "\nSequence of 3D Triangles: \n";
	// print(tr_seq3d);

	// Update the sequence around the pivot if possible
	if (update_triangle_sequence_around_pivot(triangulation, tr_seq3d, pivot))
	{
	    // Shortest path in a triangle sequence
	    current_path = 
	    shortest_path_triangle_sequence(triangulation, sce_vertex, dest_vertex,
					    tr_seq3d, list_pivots, length);
	    lengths.push_back(length);       
	    //  pivot search
	    new_pivot = pivot_in_new_list(list_pivots, list_pivots_prec, (int) lengths.size()-1, lengths, iter_pivot, pivot);
	}	
	else // no update, see next pivot
	{
	    new_pivot = next_pivot(list_pivots, iter_pivot, pivot);
	}		
    } // end while (new_pivot)
    
 
    std::ofstream os3("geom_sequence.off");
    printGeomviewTriangleSequence(os3, triangulation, tr_seq3d);
    std::ofstream os4("geom_path.vect");
    printGeomviewPath(os4, current_path);

    cout << "\nAfter " << lengths.size() << " iterations, the final path has ";
    cout << list_pivots.size() << " pivots : " ;
    print_list_int(list_pivots);    
    cout << "its length is : " 
	 << length_polygonal_path(current_path) << endl;

    //1// print_lengths(lengths);

    return;
}

//----------------------------------------------------------------------------
std::vector<int> DijkstraPT(PrTriangulation_OP& triangulation, 
			   int sce_vertex, int dest_vertex)
//----------------------------------------------------------------------------
/*
  Dijkstra computes the shortest path between sce_vertex and dest_vertex, 
  along the edges of the triangulation. It returns the polygonal solution path 
  given by its vertices, vertices of the triangulation. 
 */
{
    Dijkstra D;
    D.setGraph(&triangulation);
    D.initialize();
    D.setSource(sce_vertex);
    D.run(dest_vertex);
        
    vector<int> vertex_path;
    int curr = dest_vertex;
    vertex_path.push_back(curr);
    while (curr != sce_vertex) 
    {
	curr = D.closestNeighbour(curr);
	vertex_path.insert(vertex_path.begin(), curr);
    }
    // cout << "Shortest path in the graph : " << endl;
    // print(vertex_path);
    return(vertex_path);
}

//----------------------------------------------------------------------------
vector<Vector3D> point3D_from_indice(PrTriangulation_OP& t,
				       const vector<int> &ind_path)
//----------------------------------------------------------------------------
/*
  3D points from their indices in the triangulation
 */
{
    vector<Vector3D> path;
    std::vector<int>::const_iterator iter;

    for (iter=ind_path.begin(); iter!=ind_path.end(); iter++)
    {
	path.push_back(t.get3dNode(*iter));
    }
    // cout << "3D points from the previous indices:" << endl;
    // print(path);
    return(path);
}

//----------------------------------------------------------------------------
vector<int> tr_crossed_by_path_vertices(PrTriangulation_OP& t, 
					vector<int> &vert_path)
//----------------------------------------------------------------------------
/*
  triangle sequence crossed by a path 
  where the path is given by a list of vert. indices
 */
{
    // traiter le cas ou le sommet courant est un sommet du bord
    vector<int> tr_seq;
    std::vector<int>::const_iterator iter;
    std::vector<int>::iterator i, j;
    int curr_tr;
    int curr_vertex;
    int next_vertex;
    
    // VertPath contains at least 2 vertices
    if (vert_path.size()<2)
    {
	cout << "The initial path should contain at least 2 vertices." << endl;
	exit(0);
    }
    iter = vert_path.begin();
    curr_vertex = *iter;

    // first triangle = triangle containing the two first vertices
    curr_tr = (t.getPrNode(curr_vertex)).tr();
    iter ++;
    next_vertex = *iter;

    if (t.getPrTriangle(curr_tr).isVertex(next_vertex))
    {
      /*      if (t.isBoundary(curr_vertex))
      {
	next_vertex=t.getPrTriangle(curr_tr).getClockwiseNode(curr_vertex);
	vert_path.insert(vert_path.begin()+1, next_vertex);
	iter=vert_path.begin()+1;
	} */
      curr_tr = t.getPrTriangle(curr_tr).getLeftTriangle(curr_vertex);
    }
    // Les sommets successifs doivent etre voisins
    while (!((t.getPrTriangle(curr_tr)).isVertex(next_vertex)))
    {
	curr_tr = t.getPrTriangle(curr_tr).getLeftTriangle(curr_vertex);
    }
    tr_seq.push_back(curr_tr);
    iter ++;
    while (iter != vert_path.end())
    {
	curr_vertex = next_vertex;
	next_vertex = *iter;
	while (!((t.getPrTriangle(curr_tr)).isVertex(next_vertex)))
	{
	    tr_seq.push_back(curr_tr);
	    // subsequence turning left by default
	    if (t.getPrTriangle(curr_tr).getLeftTriangle(curr_vertex) != -1)
	    {
		curr_tr = t.getPrTriangle(curr_tr).getLeftTriangle(curr_vertex);
	    }
	    // turning right if boundary reached
	    else
	    {
		while (!((t.getPrTriangle(curr_tr)).isVertex(next_vertex)))
		{
		    curr_tr = t.getPrTriangle(curr_tr).getRightTriangle(curr_vertex);
		    tr_seq.push_back(curr_tr);
		}
	    }
	}
	tr_seq.push_back(curr_tr);
	iter ++;    
    }

    // simplify the list such that two successic triangles are not equal
    i = tr_seq.begin();
    j = tr_seq.begin()+1;
    while (j != tr_seq.end())
    {
	if ((*i)==(*j))
	{
	    j = tr_seq.erase(j);
	}
	else 
	{
	    i++; j++;
	}
    }

    // cout << "Triangle sequence (3D) crossed by this path : " << endl;
    // print(tr_seq);
    
    // print_edges_lengths(tr_seq, t);
    return(tr_seq);    
}

//----------------------------------------------------------------------------
bool update_triangle_sequence_around_pivot(
    PrTriangulation_OP& t, 
    vector<int>& tr_seq, int pivot)
//----------------------------------------------------------------------------
/*
  a triangle sequence is modified/updated around a given pivot vertex
  returns 
  - 0 if it is not possible to update 
  (if pivot belongs to the boundary or pivot does not belong to the sequence)
  - 1 otherwise
 */
{
    std::vector<int>::iterator iter_tr = tr_seq.begin();
    std::vector<int>::iterator iter_pivot_first; iter_pivot_first = tr_seq.end();

    // triangles around the pivot 
    int first_tr = -1; 
    int last_tr = -1; 
    std::list<int> tr_subseq;
    int right_left;

    // boundary vertex : no updating possible
    if (t.isBoundary(pivot))
    {
	cout << "\nNo possible update around boundary vertex " << pivot;
	cout << endl;
	return(0);
    }

    while ( (!(t.getPrTriangle(*iter_tr)).isVertex(pivot)) &&
	    (iter_tr < tr_seq.end() ) ) 
    {
	iter_tr ++;
    }
    if (iter_tr == tr_seq.end())
    {
	std::cout << "\n Updating the sequence around this vertex is impossible" 
		  << "\n because this vertex does not belong to the sequence !" 
		  << std::endl;
	return(0);
    }

    iter_pivot_first = iter_tr;
    first_tr = *iter_tr;

    // triangles containing the pivot, ordered clockwise, starting from first_tr
    tr_subseq = tr_sequence_around_vertex(t, pivot, first_tr);
    // cout << "\ntr_sub_seq around the Pivot " << pivot << endl;
    // print_list_int(tr_subseq);
    iter_tr ++;
    
    if ((*iter_tr) == *(++ tr_subseq.begin()))
    {
	right_left = 1;
    }
    else right_left = 0;

    while ( (iter_tr != tr_seq.end()) && ((t.getPrTriangle(*iter_tr)).isVertex(pivot)) ) 
    {
	// cout << "\n" << (*iter_tr);
 	iter_tr ++;
    }
    iter_tr --;
    last_tr =(*iter_tr);

    // erase the subsequence to modify
    iter_tr = ++iter_pivot_first;
    while ((*iter_tr)!=last_tr)
    {
	// cout << "\n" << (*iter_tr);
	iter_tr = tr_seq.erase(iter_tr);
    }

    // insert the new subsequence
    std::list<int>::iterator iter_tr_subseq;
    if (right_left == 0)
    {
	iter_tr_subseq = ++tr_subseq.begin();
	while ((*iter_tr_subseq)!=last_tr)
	{
	    iter_tr = tr_seq.insert(iter_tr, *iter_tr_subseq);
	    iter_tr++;
	    iter_tr_subseq++;
	}
    }
    else // right_left = 1
    {
	iter_tr_subseq = --tr_subseq.end();
	while ((*iter_tr_subseq)!=last_tr)
	{
	    iter_tr = tr_seq.insert(iter_tr, *iter_tr_subseq);
	    iter_tr ++;
	    iter_tr_subseq--;
	}
    }
    // cout << "\nNew sequence of triangles : " << endl;
    // print(tr_seq);
    // print_edges_lengths(tr_seq, t);
    return(1);   
}

//----------------------------------------------------------------------------
std::list<int> tr_sequence_around_vertex(PrTriangulation_OP& t, 
				   int vertex)
//----------------------------------------------------------------------------
/*
 triangles containing a vertex, ordered clockwise
 */  
{
    std::list<int> tr_seq;
    int tr = t.getPrNode(vertex).tr();

    do
    {
	tr_seq.push_back(tr);
	tr = t.getPrTriangle(tr).getLeftTriangle(vertex);
    }
    while (tr != (t.getPrNode(vertex).tr()));

    return(tr_seq);
}


//----------------------------------------------------------------------------
std::list<int> tr_sequence_around_vertex(PrTriangulation_OP& t, 
				   int vertex, int first_tr)
//----------------------------------------------------------------------------
/*
 triangles containing a vertex, ordered clockwise, starting from first_tr
 */  
{
    std::list<int> tr_seq;
    int tr = first_tr;

    do
    {
	tr_seq.push_back(tr);
	tr = t.getPrTriangle(tr).getLeftTriangle(vertex);
    }
    while (tr != first_tr);

    return(tr_seq);
}


//----------------------------------------------------------------------------
bool pivot_in_new_list(std::list<int> &list_pivots, 
		       std::vector<int>::const_iterator &tried_begin,
		       std::vector<int>::const_iterator &tried_end,
		       std::list<int>::iterator& iter_pivot)
//----------------------------------------------------------------------------
{
   if (list_pivots.empty())
   {
      return(0);
   }
   iter_pivot = list_pivots.begin();
   while (iter_pivot != list_pivots.end())
   {
      std::vector<int>::const_iterator ind_it;
      ind_it=find(tried_begin, tried_end, *iter_pivot);
      if (ind_it==tried_end)
      {
	 break;
      }
      iter_pivot++;
   }
   if (iter_pivot==list_pivots.end())
   {
      return (0);
   } else
      return (1);
}

//----------------------------------------------------------------------------
bool pivot_in_new_list(std::list<int> &list_pivots, 
		       const std::list<int> &list_pivots_prec, 
		       int nb_steps, 
		       const vector<double>& lengths,
		       std::list<int>::iterator& iter_pivot, int prec_pivot)
//----------------------------------------------------------------------------
/*
  decide if there are pivot vertices in the new sequence 
  set the pivot vertex where the sequence is to update in iter_pivot.
  compare the pivot lists before and after updating :
  if equal : incrementation of the pivot
  if not : first pivot in the new list
 */
{
    if (list_pivots.empty())
    {
	return(0);
    }
    iter_pivot = list_pivots.begin();
    if (list_pivots == list_pivots_prec)
    {
	while ((*iter_pivot) != prec_pivot)
	{
	    iter_pivot ++;
	}
	if (decreasing_length(nb_steps, lengths))
	{
	    iter_pivot ++;
	    if (iter_pivot == list_pivots.end())
	    {
		return(0);
	    }    
       	}
    }
    // case where the update has not changed the path locally but because
    // of numerical unstability, changed the pivot list
    else
    {
	std::list<int>::iterator i;
	i = list_pivots.begin();
	while ((i != list_pivots.end()) && (*i != prec_pivot))
	{
	    i++;
	}
	if (*i == prec_pivot)
	{
	    i++;
	    if (i != list_pivots.end())
	    {
		iter_pivot = i;
	    }
	    else return(0);
	}
    }
    return(1);
}

//----------------------------------------------------------------------------
bool next_pivot(std::list<int> &list_pivots, 
		std::list<int>::iterator& iter_pivot, int prec_pivot)
//----------------------------------------------------------------------------
/*
  the treated pivot belongs to the boundary, consider the next one in the list.
 */
{
    std::list<int>::iterator i;
    i = list_pivots.begin();
    while ((i != list_pivots.end()) && (*i != prec_pivot))
    {
	i++;
    }
    if (*i == prec_pivot)
    {
	i++;
	if (i != list_pivots.end())
	{
	    iter_pivot = i;
	    return(1);
	}
    }
    return(0);
}

//----------------------------------------------------------------------------
bool decreasing_length(const int nb_steps, const vector<double>& lengths)
//----------------------------------------------------------------------------
/*
  if updating leaves the shortest path in the sequence unchanged,
  this function checks that the path has (numerically) decreasing length.
  if not, comes back to the previous sequence
 */
{
    if (nb_steps>0)
    {
	if (lengths[nb_steps]>lengths[nb_steps-1])
	{				
	    //1// fprintf(stderr, "\n***----- numerical increasing length ! -----*** : ");
	    //1// fprintf(stderr, "%5.15e \n", lengths[nb_steps]-lengths[nb_steps-1]);
	    return(0);
       	}
    }
    return(1);
}

//----------------------------------------------------------------------------
bool decreasing_length(const vector<double>& lengths)
//----------------------------------------------------------------------------
/*
  if updating leaves the shortest path in the sequence unchanged,
  this function checks that the path has (numerically) decreasing length.
  if not, comes back to the previous sequence
 */
{
    int l = (int) lengths.size();
    if (l>1)
    {
	double a, b;
	a = lengths[l-1];
	b = lengths[l];
	if (b>a)
	{				
	    //1// fprintf(stderr, "\n***----- numerical increasing length ! -----*** : ");
	    //1// fprintf(stderr, "%5.15e \n", b-a);
	    return(0);
       	}
    }
    return(1);
}

//----------------------------------------------------------------------------
double dist3D(Vector3D& a, Vector3D& b) 
//----------------------------------------------------------------------------
/* distance between two R3 points */
{
    return(sqrt(sqr(a.x()-b.x())+sqr(a.y()-b.y())+sqr(a.z()-b.z())));
}

//----------------------------------------------------------------------------
double dist2D(Vector2D& a, Vector2D& b) 
//----------------------------------------------------------------------------
/* distance between two R2 points */
{
    return(sqrt(sqr(a.x()-b.x())+sqr(a.y()-b.y())));
}

//----------------------------------------------------------------------------
double sqr(double x) 
//----------------------------------------------------------------------------
{
    return(x*x);
}

//----------------------------------------------------------------------------
double length_polygonal_path(const vector<Vector3D> &path)
//----------------------------------------------------------------------------
{
    double l; l = 0;
    Vector3D pt_i, pt_f;
    std::vector<Vector3D>::const_iterator iter_i, iter_f;

    if (path.size()<2)
    {
	fprintf(stderr, "\n this path should "
		"at least contain the sce and dest. vertices");
	exit(0);
    } 

    iter_i = path.begin();
    iter_f = path.begin(); iter_f ++;
    while (iter_f != path.end())
    {
	pt_i = *iter_i;
	pt_f = *iter_f;
	l += dist3D(pt_i,pt_f);
	iter_i ++;
	iter_f ++;
    }
    return(l);
}
