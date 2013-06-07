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
#include <list>
#include <string>

#include "GoTools/parametrization/PrGeodesics.h"
#include "GoTools/parametrization/PrPathTriangleSeq.h"
#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrInterface.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrTriangle.h"
//#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"

using namespace Go;
using namespace std;

//----------------------------------------------------------------------------
void EdgeType::initEdgeType()
//----------------------------------------------------------------------------
{
    vertex_[0] = -1;
    vertex_[1] = -1;
    return;
}

//-----------------------------------------------------------------------------
bool EdgeType::isVertex(int node)
//-----------------------------------------------------------------------------
{
    if (vertex_[0] == node || vertex_[1] == node) return true;
    else return false;
}

//-----------------------------------------------------------------------------
void PathType::print()
//-----------------------------------------------------------------------------
{
    std::vector<int>::iterator iter;
    cout << "\nl_path_: ";
    for (iter = l_path_.begin(); iter < l_path_.end(); iter++)
    {
      cout << *iter << " - ";
    }
    cout << endl;
    cout << "r_path_: ";
    for (iter = r_path_.begin(); iter < r_path_.end(); iter++)
    {
	cout << *iter << " - ";
    }
    cout << endl;
    cout << "funnel_: ";
    for (iter = funnel_.begin(); iter < funnel_.end(); iter++)
    {
	cout << *iter << " - ";
    }
    cout << endl;
    return;
}

//-----------------------------------------------------------------------------
void UnfNodeType::printNode()
//-----------------------------------------------------------------------------
{
    cout << n_.x() << " " << n_.y() << endl;
}

//----------------------------------------------------------------------------
void printNode(vector<UnfNodeType>& v)
//----------------------------------------------------------------------------
{
    std::vector<UnfNodeType>::iterator iter;
    for (iter = v.begin(); iter < v.end(); iter++)
    {
	iter->printNode();
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void PathType::initPathType()
//----------------------------------------------------------------------------
{
    l_path_.clear();
    r_path_.clear();
    funnel_.clear();
    return;
}

// --------------------------------------------------------------------
void PathType::funnelOfPath(PathType prec)
// --------------------------------------------------------------------
/*
  computes the funnel for the funnel algorithm
 */
{
    if ((l_path_.size() == 0) || (r_path_.size() == 0))
    {
	funnel_.clear();
	return;
    }
    if (l_path_[0] != r_path_[0])
    {
	fprintf(stderr, "the two paths do not have the same source vertex: ");
	fprintf(stderr, "impossible to compute the corresponding funnel_ ");
	cout << endl;
	return;
    }
    funnel_.insert(funnel_.begin(), prec.funnel_.begin(), prec.funnel_.end());
    int i; i = 0;
    bool test; test = true;
    while (test)
    {
	if (i!=0)
	{
	    funnel_.push_back(*(l_path_.begin()));
	}
	l_path_.erase(l_path_.begin());
	r_path_.erase(r_path_.begin());
	i++;
	if ((l_path_.size() == 0) || (r_path_.size() == 0))
	{
	    test = false;
	}
	else 
	{
	    test = (*(l_path_.begin()) == *(r_path_.begin()));
	}
    }    
    l_path_.insert(l_path_.begin(), funnel_.end()-1, (funnel_.end()));
    r_path_.insert(r_path_.begin(), funnel_.end()-1, (funnel_.end()));
    return;    
}

//----------------------------------------------------------------------------
void PathType::modification_l_path_(int new_pt, const vector<UnfNodeType>& nodes_unf)
//----------------------------------------------------------------------------
/*
  modification of the left path in the funnel algorithm
 */
{
    Vector2D a, b, c;
    c = nodes_unf[new_pt].n_;

    l_path_.push_back(new_pt);
    if (l_path_.size() > 2)
    {
	a =  nodes_unf[l_path_[l_path_.size()-3]].n_;
	b =  nodes_unf[l_path_[l_path_.size()-2]].n_;
	while ( (vect_prod2D(a, b, b, c) < 0) && (l_path_.size() > 2))
	{
	    l_path_.erase(l_path_.end()-2);
	    if (l_path_.size() > 2)
	    {
		a =  nodes_unf[l_path_[l_path_.size()-3]].n_;
		b =  nodes_unf[l_path_[l_path_.size()-2]].n_;
	    }	
	}
    }
    if (l_path_.size() == 2)
    {
	a = nodes_unf[l_path_[0]].n_;
	if (r_path_.size()>1)
	{
	    size_t j = 1;
	    b = nodes_unf[r_path_[j]].n_;
	    while ( (vect_prod2D(a, b, a, c) < 0) && (j < r_path_.size()))
	    {
		l_path_.insert((l_path_.end()-1), r_path_[j]);
		j++;
		if (j < r_path_.size())
		{
		    b = nodes_unf[r_path_[j]].n_;
		}
	    }
	}
    }    
}

//----------------------------------------------------------------------------
void PathType::modification_r_path_(int new_pt, const vector<UnfNodeType>& nodes_unf)
//----------------------------------------------------------------------------
/*
  modification of the right path in the funnel algorithm
 */
{
    Vector2D a, b, c;
    c = nodes_unf[new_pt].n_;

    r_path_.push_back(new_pt);
    if (r_path_.size() > 2)
    {
	a =  nodes_unf[r_path_[r_path_.size()-3]].n_;
	b =  nodes_unf[r_path_[r_path_.size()-2]].n_;
	while ( (vect_prod2D(a, b, b, c) > 0) && (r_path_.size() > 2))
	{
	    r_path_.erase(r_path_.end()-2);
	    if (r_path_.size() > 2)
	    {
		a =  nodes_unf[r_path_[r_path_.size()-3]].n_;
		b =  nodes_unf[r_path_[r_path_.size()-2]].n_;
	    }	
	}
    }
    if (r_path_.size() == 2)
    {
	a = nodes_unf[r_path_[0]].n_;
	if (l_path_.size()>1)
	{
	    size_t j = 1;
	    b = nodes_unf[l_path_[j]].n_;
	    while ( (vect_prod2D(a, b, a, c) > 0) && (j < l_path_.size()))
	    {
		r_path_.insert((r_path_.end()-1), l_path_[j]);
		j++;
		if (j < l_path_.size())
		{
		    b = nodes_unf[l_path_[j]].n_;
		}
	    }
	}
    }  
    return;  
}

//----------------------------------------------------------------------------
vector<Vector3D> shortest_path_triangle_sequence(
    const PrTriangulation_OP& t, 
    int sce_vertex, int dest_vertex,
    const vector<int> &tr_seq3d, list<int>& list_pivots, 
    double& length)
//----------------------------------------------------------------------------
/*
  input: triangulation and triangle sequence 
  output: 
  - shortest path in that sequence = polygonal line given by its 3D vertices
  - list of pivot vertices, ie vertices of the triangulation 
        belonging to this shortest path
  - length of this shortest path
*/
{
    vector<Vector3D> sh_path;
    vector<double> ratio_on_edge;
    vector<EdgeType> edge_seq3d;
    vector<EdgeType> edge_seq_unf;
    int sce_vertex_unf, dest_vertex_unf;
    vector<PrTriangle> tr_seq_unf;
    list<int> list_pivots_unf;
    vector<UnfNodeType> nodes_unf;
    list<double> list_deviation;

    edge_seq3d = edge_sequence(tr_seq3d, t);
    // print_tr_vertices(tr_seq3d, t);

    unfolding_triangle_sequence(t, tr_seq3d, tr_seq_unf, edge_seq3d, edge_seq_unf, nodes_unf);
    unfolding_vertex(t, tr_seq3d, tr_seq_unf, sce_vertex, sce_vertex_unf);
    unfolding_vertex(t, tr_seq3d, tr_seq_unf, dest_vertex, dest_vertex_unf);
    
    // cout << "\nShortest Path in the sequence";    
    // cout << endl << nodes_unf.size() << " nodes_unf\n";
    // printNode(nodes_unf);
    // cout << endl << tr_seq_unf.size() << " tr_seq_unf\n";
    // print(tr_seq_unf);
    // cout << endl << tr_seq3d.size() << " tr_seq3d\n";
    // print_tr_vertices(tr_seq3d, t);
         
    if (tr_seq3d.size()>2)
    {
	ratio_on_edge = shortest_path_triangle_sequence_2D(
	    t, sce_vertex, dest_vertex, sce_vertex_unf, dest_vertex_unf,
	    tr_seq3d, tr_seq_unf, 
	    list_pivots_unf, nodes_unf, length);
	//1// std::ofstream geom_unf_objects("geom_unf_objects");
	//1// printUnfoldedObjects(geom_unf_objects, nodes_unf, tr_seq_unf, list_pivots_unf);


	list_deviation = deviation_from_pivots(list_pivots_unf, nodes_unf);
	//1// print_deviation(list_deviation);

	list_pivots = pivots3D_from_pivots_in_unf_seq(list_pivots_unf, nodes_unf);
	//1// cout << "\nPivots : " ;
	//1// print_list_int(list_pivots);

        sh_path.clear();	
	sh_path.push_back(t.get3dNode(sce_vertex));
	path_3Dpts_from_ratio_on_edge(t, sh_path, ratio_on_edge, edge_seq3d);
	sh_path.push_back(t.get3dNode(dest_vertex));
	// cout << "\nShortest Path given by the 3D points, vertices of the polygonal path  :" << endl;
	// print(sh_path);

     }
    else // the source and destination vertices belong to the same triangle
    {
	length = dist2D(nodes_unf[dest_vertex_unf].n_,nodes_unf[sce_vertex_unf].n_);
	
	list_deviation.clear();
	list_pivots.clear();
	//1// 
cout << "\nNo pivot, no deviation";

	sh_path.clear();
	sh_path.push_back(t.get3dNode(sce_vertex));
	sh_path.push_back(t.get3dNode(dest_vertex));
        // cout << "\nShortest Path given by the 3D points, vertices of the polygonal path  :" << endl;
	// print(sh_path);
    }

    return(sh_path);
}




//----------------------------------------------------------------------------
void shortest_path_triangle_sequence(
    const PrTriangulation_OP& t, 
    const int& sce_vertex, const int& dest_vertex,
    vector<int> &tr_seq3d, list<int>& list_pivots, 
    std::vector<EdgeType> &edge_seq3d, std::vector<double> &ratio_on_edge,
    double& length)
//----------------------------------------------------------------------------
/*
  input: triangulation and triangle sequence 
  output: 
  - shortest path in that sequence = polygonal line given by its 3D vertices
  - list of pivot vertices, ie vertices of the triangulation 
        belonging to this shortest path
  - length of this shortest path
*/
{
  //vector<Vector3D> sh_path;
    vector<EdgeType> edge_seq_unf;
    int sce_vertex_unf, dest_vertex_unf;
    vector<PrTriangle> tr_seq_unf;
    list<int> list_pivots_unf;
    vector<UnfNodeType> nodes_unf;
    list<double> list_deviation;

    edge_seq3d = edge_sequence(tr_seq3d, t);
    // print_tr_vertices(tr_seq3d, t);

    unfolding_triangle_sequence(t, tr_seq3d, tr_seq_unf, edge_seq3d, edge_seq_unf, nodes_unf);
    unfolding_vertex(t, tr_seq3d, tr_seq_unf, sce_vertex, sce_vertex_unf);
    unfolding_vertex(t, tr_seq3d, tr_seq_unf, dest_vertex, dest_vertex_unf);
    
    // cout << "\nShortest Path in the sequence";    
    // cout << endl << nodes_unf.size() << " nodes_unf\n";
    // printNode(nodes_unf);
    // cout << endl << tr_seq_unf.size() << " tr_seq_unf\n";
    // print(tr_seq_unf);
    // cout << endl << tr_seq3d.size() << " tr_seq3d\n";
    // print_tr_vertices(tr_seq3d, t);
         
    if (tr_seq3d.size()>2)
    {
	ratio_on_edge = shortest_path_triangle_sequence_2D(
	    t, sce_vertex, dest_vertex, sce_vertex_unf, dest_vertex_unf,
	    tr_seq3d, tr_seq_unf, 
	    list_pivots_unf, nodes_unf, length);
	//1// std::ofstream geom_unf_objects("geom_unf_objects");
	//1// printUnfoldedObjects(geom_unf_objects, nodes_unf, tr_seq_unf, list_pivots_unf);


	list_deviation = deviation_from_pivots(list_pivots_unf, nodes_unf);
	//1// print_deviation(list_deviation);

	list_pivots = pivots3D_from_pivots_in_unf_seq(list_pivots_unf, nodes_unf);
	//1// cout << "\nPivots : " ;
	//1// print_list_int(list_pivots);

	// cout << "\nShortest Path given by the 3D points, vertices of the polygonal path  :" << endl;
	// print(sh_path);

     }
    else // the source and destination vertices belong to the same triangle
    {
	length = dist2D(nodes_unf[dest_vertex_unf].n_,nodes_unf[sce_vertex_unf].n_);
	
	list_deviation.clear();
	list_pivots.clear();
	//1// 
cout << "\nNo pivot, no deviation";

        // cout << "\nShortest Path given by the 3D points, vertices of the polygonal path  :" << endl;
	// print(sh_path);
    }

}







//----------------------------------------------------------------------------
void path_3Dpts_from_ratio_on_edge(
    const PrTriangulation_OP& t, vector<Vector3D>& sh_path,
    const vector<double>& ratio_on_edge, const vector<EdgeType>& edge_seq3d)
//----------------------------------------------------------------------------
/*
  computes the 3D points of the path from their position (ratio) on
  the edges of the sequence
 */
{
    Vector3D pt;
    size_t i;

    for (i = 0; i < edge_seq3d.size(); i++)
    {
	pt = pt3D_from_ratio_on_edge(ratio_on_edge[i], 
				     t.get3dNode(edge_seq3d[i].vertex_[0]), 
				     t.get3dNode(edge_seq3d[i].vertex_[1]));
	sh_path.push_back(pt);
    }
    return;
}

//----------------------------------------------------------------------------
Vector3D pt3D_from_ratio_on_edge(double r, Vector3D a, Vector3D b)
//----------------------------------------------------------------------------
/*
  computes a 3D point from its position (ratio) on an edge
 */
{
    return(a + r*(b-a));
}


//----------------------------------------------------------------------------
void path_2Dpts_from_ratio_on_edge(
    const vector<UnfNodeType> &nodes_unf, vector<Vector2D>& path,
    const vector<double>& ratio_on_edge, const vector<EdgeType>& edge_seq2d)
//----------------------------------------------------------------------------
/*
  computes the 2D points of the path from their position (ratio) on
  the edges of the sequence
 */
{
    Vector2D pt;
    size_t i;

    for (i = 0; i < edge_seq2d.size(); i++)
    {
	pt = pt2D_from_ratio_on_edge(ratio_on_edge[i], 
				     nodes_unf[edge_seq2d[i].vertex_[0]].n_, 
				     nodes_unf[edge_seq2d[i].vertex_[1]].n_);
	path.push_back(pt);
    }
    return;
}

//----------------------------------------------------------------------------
Vector2D pt2D_from_ratio_on_edge(double r, Vector2D a, Vector2D b)
//----------------------------------------------------------------------------
/*
  computes a 2D point from its position (ratio) on an edge
 */
{
    return(a + r*(b-a));
}


//----------------------------------------------------------------------------
list<int> pivots3D_from_pivots_in_unf_seq(
    const list<int> &list_pivots_unf, const vector<UnfNodeType> &nodes_unf)
//----------------------------------------------------------------------------
/*
  computes the 3d pivot list corresponding to the pivot in the unfolded sequence.
  The pivots are the vertices of the shortest path.
  The pivots in the unfolded sequence contain all the vertices of the path (source and dest. too).
  The 3d pivots do not contain the source and dest. vertices.
 */
{
    list<int> list_pivots;
    std::list<int>::const_iterator iter;

    iter = list_pivots_unf.begin();
    iter++;
    int i;
    for (i=1; i<(int)list_pivots_unf.size()-1; i++)
    {
	list_pivots.push_back((nodes_unf[*iter]).origine_);
	iter ++;
    }
    return(list_pivots);
}

//----------------------------------------------------------------------------
list<double> deviation_from_pivots(
    const list<int> &list_pivots_unf, const vector<UnfNodeType>& nodes_unf)
//----------------------------------------------------------------------------
/*
  At each interior pivot, the path is deviated from the straight line.
  This function computes the deviation angles.
 */
{
    list<double> list_deviation;
    list_deviation.clear();
    if (list_pivots_unf.size() > 2)
    {
	std::list<int>::const_iterator i, j, k;
	i = list_pivots_unf.begin();
	j = list_pivots_unf.begin(); j++;
	k = list_pivots_unf.begin(); k++; k++;

	Vector2D a, b, c;
	a = nodes_unf[*i].n_;
	b = nodes_unf[*j].n_;
	c = nodes_unf[*k].n_;
	list_deviation.push_back(abs_val(angle(a, b, b, c)));
	i++; j++; k++;
	while (k != list_pivots_unf.end())
	{
	    a = nodes_unf[*i].n_;
	    b = nodes_unf[*j].n_;
	    c = nodes_unf[*k].n_;
	    list_deviation.push_back(abs_val(angle(a, b, b, c)));
	    i++; j++; k++;
	}
    }
    return(list_deviation);
}

//----------------------------------------------------------------------------
double abs_val(double a)
//----------------------------------------------------------------------------
/* absolute value */
{
    if (a>0)
	return(a);
    else return(-a);   
}

// --------------------------------------------------------------------
vector<double> shortest_path_triangle_sequence_2D(
    const PrTriangulation_OP &t, 
    const int sce_vertex, const int dest_vertex,
    const int sce_vertex_unf, int dest_vertex_unf,
    const vector<int>& tr_seq3d, vector<PrTriangle>& tr_seq_unf, 
    list<int>& list_pivots_unf, 
    vector<UnfNodeType> nodes_unf, double& length)
// --------------------------------------------------------------------
/*
  computes the shortes path in an ordered 2d seuqence of triangles, 
  corresponding to a 3d seuqence that was unfolded
 */
{
    vector<double> ratio_on_edge;
    vector<PathType> path_in_edge_seq;
    PathType p;
    vector<EdgeType> edge_seq_unf;

    edge_seq_unf = edge_sequence(tr_seq_unf);

    // check that the source vertex belongs to the first triangle
    if (! (tr_seq_unf[0].isVertex(sce_vertex_unf)))
    {
	fprintf(stderr, "the source vertex is not in the first triangle of the sequence");
	cout << endl;
	return(ratio_on_edge);
    }
    
    size_t i = 0;
    // look for the first edge in the seq, not containing the source vertex
    while ((edge_seq_unf[i].isVertex(sce_vertex_unf)) && (i < edge_seq_unf.size()))
    {
	p.initPathType();
	path_in_edge_seq.push_back(p);
	path_in_edge_seq[i].l_path_.push_back(sce_vertex_unf);
	path_in_edge_seq[i].r_path_.push_back(sce_vertex_unf);
	path_in_edge_seq[i].funnel_.push_back(sce_vertex_unf);

	// cout << "\nedge" << i << ": ";
	// path_in_edge_seq[i].print();
	i ++;
    } 

    p.initPathType();
    path_in_edge_seq.push_back(p);
    path_in_edge_seq[i].l_path_.push_back(sce_vertex_unf);
    path_in_edge_seq[i].r_path_.push_back(sce_vertex_unf);
    Vector2D a, b, c;
    a = nodes_unf[sce_vertex_unf].n_;
    b = nodes_unf[edge_seq_unf[i].vertex_[0]].n_;
    c = nodes_unf[edge_seq_unf[i].vertex_[1]].n_;
    if (vect_prod2D(a, b, a, c) > 0)
    {
	path_in_edge_seq[i].l_path_.push_back(edge_seq_unf[i].vertex_[1]);
	path_in_edge_seq[i].r_path_.push_back(edge_seq_unf[i].vertex_[0]);
    }
    else
    {
	path_in_edge_seq[i].l_path_.push_back(edge_seq_unf[i].vertex_[0]);
	path_in_edge_seq[i].r_path_.push_back(edge_seq_unf[i].vertex_[1]);
    }
    path_in_edge_seq[i].funnel_.push_back(sce_vertex_unf);

    // cout << "origine edge unf " << i << ": ";
    // path_in_edge_seq[i].print();
   i++;
    // treat each edge until reaching the dest vertex
    while ((!(edge_seq_unf[i].isVertex(dest_vertex_unf))) && (i < edge_seq_unf.size()))
    {
	p.initPathType();
	path_in_edge_seq.push_back(p);
	path_to_edge(sce_vertex_unf, (int)i, path_in_edge_seq, nodes_unf, tr_seq_unf, edge_seq_unf);
	i++;
    }
    while (i < edge_seq_unf.size())
    {
	p.initPathType();
	path_in_edge_seq.push_back(p);
	path_in_edge_seq[i] = path_in_edge_seq[i-1];
	i++;
    }
    p = path_in_edge_seq[path_in_edge_seq.size()-1];
    list_pivots_unf.clear();

    // Dest on the last edge
    if  (dest_vertex_unf == (*(p.l_path_.end()-1))) 
    {
	p.l_path_.erase(p.l_path_.begin());
	list_pivots_unf.insert(list_pivots_unf.begin(), p.funnel_.begin(), p.funnel_.end());
	list_pivots_unf.insert(list_pivots_unf.end(), p.l_path_.begin(), p.l_path_.end());
    }
    else
    {
	if (dest_vertex_unf == (*(p.r_path_.end()-1)))
	{
	    p.r_path_.erase(p.r_path_.begin());
	    list_pivots_unf.insert(list_pivots_unf.begin(), p.funnel_.begin(), p.funnel_.end());
	    list_pivots_unf.insert(list_pivots_unf.end(), p.r_path_.begin(), p.r_path_.end());
	}
	// Dest vertex = 3rd vertex of the last triangle
	else
	{
 	    p.modification_l_path_(dest_vertex_unf, nodes_unf);
	    p.funnel_.clear();
	    p.funnelOfPath(*(path_in_edge_seq.end()-1));

	    // cout << "\ndest_vertex" << ": ";
	    // P.print();

	    p.l_path_.erase(p.l_path_.begin());
	    list_pivots_unf.insert(list_pivots_unf.begin(), p.funnel_.begin(), p.funnel_.end());
	    list_pivots_unf.insert(list_pivots_unf.end(), p.l_path_.begin(), p.l_path_.end());
	}
    }
    // cout << "\nList pivots unf : ";
    // print_list_int(list_pivots_unf);

    ratio_on_edge = 
	path_ratio_on_edge_from_pivots(list_pivots_unf, edge_seq_unf, nodes_unf);
	
    // print_edges_lengths(tr_seq_unf, nodes_unf);
    // print_triangle_sequence(tr_seq_unf, nodes_unf);

    length = length_polygonal_path(list_pivots_unf, nodes_unf);
    //1// cout << "\nLength of the path: " << length;
    return(ratio_on_edge);
}

//----------------------------------------------------------------------------
void print_deviation(list<double>& list_deviation)
//----------------------------------------------------------------------------
{
    cout << "\nDeviation" << ": ";
    std::list<double>::iterator it;
    for (it = list_deviation.begin(); it != list_deviation.end(); it++)
    {
	cout << *it << " ";
    }
    cout << endl;
}

//----------------------------------------------------------------------------
void print_list_int(list<int>& list_pivots)
//----------------------------------------------------------------------------
{
    std::list<int>::iterator it;
    for (it = list_pivots.begin(); it != list_pivots.end(); it++)
    {
	cout << *it << " ";
    }
    cout << endl;
}

//----------------------------------------------------------------------------
double length_polygonal_path(
    const list<int> &list_pivots_unf, vector<UnfNodeType>& nodes_unf)
//----------------------------------------------------------------------------
{
    double l; l = 0;
    Vector2D pt_i, pt_f;
    std::list<int>::const_iterator iter_i, iter_f;

    if (list_pivots_unf.size()<2)
    {
	fprintf(stderr, "\n this pivot list in the unfolded doamin should "
		"at least contain the sce and dest. vertices");
	exit(0);
    } 

    iter_i = list_pivots_unf.begin();
    iter_f = list_pivots_unf.begin(); iter_f ++;
    while (iter_f != list_pivots_unf.end())
    {
	pt_i = nodes_unf[*(iter_i)].n_;
	pt_f = nodes_unf[*(iter_f)].n_;
	l += dist2D(pt_i,pt_f);
	iter_i ++;
	iter_f ++;
    }
    return(l);
}

//----------------------------------------------------------------------------
vector<double> path_ratio_on_edge_from_pivots(
    list<int>& list_pivots_unf, vector<EdgeType>& edge_seq_unf, 
    const vector<UnfNodeType> &nodes_unf)
//----------------------------------------------------------------------------
/*
  computes the intersection of the path, given by the pivots, with the
  edge sequence, returns the ratio/position on each edge
 */
{
    vector<double> p;
    int pt_i, pt_f;
    std::list<int>::iterator iter_i, iter_f, iter_lim;
    std::vector<EdgeType>::iterator iter_e;
    double ratio;

    if (list_pivots_unf.size()<2)
    {
	fprintf(stderr, "\n this pivot list in the unfolded domain should "
		"at least contain the sce and dest. vertices");
	exit(0);
    } 

    iter_i = list_pivots_unf.begin();
    iter_f = list_pivots_unf.begin();
    iter_f ++;
    pt_i = *(iter_i);
    pt_f = *(iter_f);

    for (iter_e = edge_seq_unf.begin(); iter_e != edge_seq_unf.end(); iter_e ++)
    {
	// interior intersection of [pt_i, pt_f] with the edge
	if ((!((*iter_e).isVertex(pt_i))) && (!((*iter_e).isVertex(pt_f)))) 
	{
	    ratio = segment_intersection(nodes_unf[pt_i].n_, nodes_unf[pt_f].n_,
		nodes_unf[(*iter_e).vertex_[0]].n_, nodes_unf[(*iter_e).vertex_[1]].n_);
	    p.push_back(ratio);
	}
	else
	{
	    // pt_i belongs to the edge
	    if (
		((iter_e->isVertex(pt_i)) && (!(iter_e->isVertex(pt_f))))
		||
		((iter_e->isVertex(pt_i)) && ((iter_e->isVertex(pt_f))) && (p.empty()))
		)
	    {
		if (pt_i == (iter_e->vertex_[0]))
		{
		    p.push_back(0);
		}
		else 
		{
		    if (pt_i == (iter_e->vertex_[1]))
		    {
			p.push_back(1);
		    }
		}
	    }
	    iter_lim = list_pivots_unf.end();
	    iter_lim --;
	    // pt_f belongs to the edge
	    if ( (iter_e->isVertex(pt_f)) 
		 &&
		 (! ((iter_e->isVertex(pt_i)) && ((int) p.size() == 1)) )
		)
	    {
		if (pt_f == (iter_e->vertex_[0]))
		{
		    p.push_back(0);
		}
		else
		{
		    if (pt_f == (iter_e->vertex_[1]))
		    {
			p.push_back(1);
		    }
		}
		if (iter_f != iter_lim)
		{
		    iter_i ++;
		    iter_f ++;
		    pt_i = *(iter_i);
		    pt_f = *(iter_f);
		}
	    }
	}
    }
    return (p);    
}

//----------------------------------------------------------------------------
double segment_intersection(Vector2D a, Vector2D b, 
			    Vector2D c, Vector2D d)
//----------------------------------------------------------------------------
/* 
   intersection between [ab] and [cd].
   returns the position of the intersection on the oriented segment [cd] given by a ratio
 */
{
    double ratio = -1.0; // indicates error
    int found;
    Vector2D i;

    double det, det1, det2;
    det = vect_prod2D(a, b, c, d);
    if (fabs(det/dist2D(a,b)/dist2D(c,d)) < EPS)
    {
	det1 = vect_prod2D(a, d, c, b)/dist2D(a,d)/dist2D(c,b);
	det2 = vect_prod2D(a, c, b, d)/dist2D(a,c)/dist2D(b,d);

	if ((fabs(det1) < EPS)&&(fabs(det2) < EPS))
	{
	    // the lines are the same, infinity of solutions 
	    found = 2;
	    ratio=0.5;
	}
	else
	{
	    // no solution, the lines are parallel
	    found = 0;
	}
    }
    else
    {	
	ratio = vect_prod2D(a, c, a, b)/det; 
	i = c + ratio*(d-c);
	found = 1;
    }
    if ((found==0) || (ratio<0) || (ratio>1))
    {
	fprintf(stderr, "\n Error: no intersection of the path "
		"and the edge in the sequence");
	/*if (found!=1)
	exit(0);
	  else*/ if (ratio<0)
	  ratio=0;
	else
	  ratio=1;
    }

    return(ratio);  
}

//----------------------------------------------------------------------------
int line_intersection( Vector2D a, Vector2D b, 
		       Vector2D c, Vector2D d, 
		       Vector2D& i)
//----------------------------------------------------------------------------
/*
  computes i, the intersection between (ab) and (cd)
  returns 1 if there is a unique solution,
  2 if the lines are the same, 0 if they are parallel
 */
{
    int found;

    double det, den, det1, det2;
    det = vect_prod2D(a, b, c, d);
    if (abs_val(det/dist2D(a,b)/dist2D(c,d)) < EPS)
    {
	det1 = vect_prod2D(a, d, c, b)/dist2D(a,d)/dist2D(c,b);
	det2 = vect_prod2D(a, c, b, d)/dist2D(a,c)/dist2D(b,d);

	if ((fabs(det1) < EPS)&&(fabs(det2) < EPS))
	{
	    // the lines are the same, infinity of solutions 
	    found = 2;
	}
	else
	{
	    // no solution, the lines are parallel
	    found = 0;
	}
    }
    else
    {	
	den = (b - a).x();
	if (abs_val(den) > EPS)
	{	
	    i.x() = ((d.x()-c.x())*(b.x()*a.y()-a.x()*b.y()) + den*(d.y()*c.x()-d.x()*c.y()))/det;
	    i.y() = (b.y()-a.y())*(i.x()-a.x())/den+ a.y();
	}
	/* if (AB) vertical, (CD) is not (otherwise they are parallel) */
	else	
	{
	    den = d.x() - c.x();
	    i.x() = a.x();
	    i.y() = (d.y()-c.y())*(i.x()-c.x())/den+c.y();
	}
	found = 1;
    }
    return(found);
}

// --------------------------------------------------------------------
void path_to_edge( 
    const int sce_vertex_unf, int i,
    vector<PathType>& p, const vector<UnfNodeType>& nodes_unf, 
    const vector<PrTriangle>& tr_seq_unf, vector<EdgeType>& edge_seq_unf)
// --------------------------------------------------------------------
/* computes the path from a source point to an edge of the 2d sequence:
   this path is given by 
   - a l_path_, polygonal chain of vertices to the left vertex of the edge
   - a r_path_, polygonal chain of vertices to the right vertex of the edge
   - funnel_, common part of the path 
 */
{    

    p[i].initPathType();
    p[i].l_path_.insert(p[i].l_path_.begin(), p[i-1].l_path_.begin(), p[i-1].l_path_.end());
    p[i].r_path_.insert(p[i].r_path_.begin(), p[i-1].r_path_.begin(), p[i-1].r_path_.end());

    // considers the new point brought by edge i
    int old_pt, new_pt;
    if ( (!(edge_seq_unf[i-1].isVertex(edge_seq_unf[i].vertex_[0]))) 
	&& (edge_seq_unf[i-1].isVertex(edge_seq_unf[i].vertex_[1])) )
    {
	old_pt = edge_seq_unf[i].vertex_[1];
	new_pt = edge_seq_unf[i].vertex_[0];
    }
    else
    {
	if ( (!(edge_seq_unf[i-1].isVertex(edge_seq_unf[i].vertex_[1]))) 
	    && (edge_seq_unf[i-1].isVertex(edge_seq_unf[i].vertex_[0])) )
	{
	    old_pt = edge_seq_unf[i].vertex_[0];
	    new_pt = edge_seq_unf[i].vertex_[1];
	}
	else
	{
	    fprintf(stderr, "The edge sequence is not valide: "
		   "two consecutiv edges should have a common vertex");
	    cout << endl;
	    return;
	}
    }
    if (old_pt == (p[i-1].l_path_[(int)p[i-1].l_path_.size()-1]))
    {
	p[i].modification_r_path_(new_pt, nodes_unf);
    }
    else
    { 
	if (old_pt == (p[i-1].r_path_[p[i-1].r_path_.size()-1]))
	{
	    p[i].modification_l_path_(new_pt, nodes_unf);
	}
	else
	{
	    fprintf(stderr, "Error in the computation of l_path_ and r_path_ ");
	    cout << endl;
	    return;
	}
    }
    p[i].funnelOfPath(p[i-1]);

    // cout << "\nedge " << i << ": ";
    // p[i].print();

    return;    
}

//----------------------------------------------------------------------------
vector<EdgeType> edge_sequence(
    const vector<int> &tr_seq, const PrTriangulation_OP& t)
//----------------------------------------------------------------------------
/*
  given a triangle sequence, computes an edge sequence consisting of 
  the common edges to two successiv triangles in the sequence
 */
{
    vector<EdgeType> edge_seq;
    std::vector<int>::const_iterator iter;
 
    for (iter=tr_seq.begin()+1; iter<tr_seq.end(); iter++)
    {
	// cout << *(iter-1) << " and " << *(iter) << " --- ";
	edge_seq.push_back(
	    common_edge(t.getPrTriangle(*(iter-1)), t.getPrTriangle(*(iter))));
	// cout << common_edge(t.getPrTriangle(*(iter-1)), t.getPrTriangle(*(iter))).vertex_[0];
	// cout << " " ;
	// cout << common_edge(t.getPrTriangle(*(iter-1)), t.getPrTriangle(*(iter))).vertex_[1];
	// cout << endl;
    }
    return(edge_seq);
}

//----------------------------------------------------------------------------
vector<EdgeType> edge_sequence(vector<PrTriangle> tr_seq)
//----------------------------------------------------------------------------
/*
  given a triangle sequence, computes an edge sequence consisting of 
  the common edges to two successiv triangles in the sequence
 */
{
    vector<EdgeType> edge_seq;
    std::vector<PrTriangle>::const_iterator iter;
 
    for (iter=tr_seq.begin()+1; iter != tr_seq.end(); iter++)
    {
	edge_seq.push_back(common_edge(*(iter-1), *(iter)));
    }
    return(edge_seq);
}

//----------------------------------------------------------------------------
EdgeType common_edge(PrTriangle t1, PrTriangle t2)
//----------------------------------------------------------------------------
/*
  common edge between two triangles, with the nodes given in the order of t1
 */
{
    EdgeType e;
    e.initEdgeType();
    if (t2.isVertex(t1.n1()))
    {
	if (t2.isVertex(t1.n2()))
	{
	    e.vertex_[0]=t1.n1();
	    e.vertex_[1]=t1.n2();
	}
	else
	{
	    if (t2.isVertex(t1.n3()))
	    {
		e.vertex_[0]=t1.n3();
		e.vertex_[1]=t1.n1();
	    }
	    else 
	    {
		fprintf(stderr, "These two triangles do not share a common edge");
		return(e);
	    }
	}
    }
    else
    {
	if ((t2.isVertex(t1.n2())) && (t2.isVertex(t1.n3())))
	{
	    e.vertex_[0]=t1.n2();
	    e.vertex_[1]=t1.n3();
	}
	else 
	{
	    fprintf(stderr, "These two triangles do not share a common edge");
	    return(e);
	}
    } 
    return(e);
}

//----------------------------------------------------------------------------  
void unfolding_triangle_sequence(
    const PrTriangulation_OP& t, 
    const vector<int>& tr_seq3d, vector<PrTriangle>& tr_seq_unf,
    const vector<EdgeType>& edge_seq3d, vector<EdgeType>& edge_seq_unf,
    vector<UnfNodeType>& nodes_unf)             
//----------------------------------------------------------------------------
/* 
   unfolding of a sequence of triangle 3D from a triangulation.
   returns a set of unfolded nodes, a sequence of unfolded triangles and 
   a sequence of unfolded edges
 */
{
    unfolding_first_triangle( t, tr_seq3d, tr_seq_unf, nodes_unf);    
    
    size_t i; 
    Vector3D U, V, W;
    
    // unfolding other triangles in the sequence
    for (i = 1; i < tr_seq3d.size() ; i++)
    {        
	// local frame of this triangle
        local_frame(t.get3dNode((t.getPrTriangle(tr_seq3d[i])).n1()), 
		    t.get3dNode((t.getPrTriangle(tr_seq3d[i])).n2()), 
		    t.get3dNode((t.getPrTriangle(tr_seq3d[i])).n3()),
		    U, V, W);
	
	// local coordinates in this triangle
	Vector3D vertex1, vertex2, vertex3;
       	vertex1 = local_coordinates(t.get3dNode(((t.getPrTriangle(tr_seq3d[i])).n1())), U, V, W);
	vertex2 = local_coordinates(t.get3dNode(((t.getPrTriangle(tr_seq3d[i])).n2())), U, V, W);    
	vertex3 = local_coordinates(t.get3dNode(((t.getPrTriangle(tr_seq3d[i])).n3())), U, V, W);
                                                                                       
        // check that the vertices have same z-coordinate
	if ((abs_val(vertex1.z()-vertex2.z())>EPS) || 
	    (abs_val(vertex2.z()-vertex3.z())>EPS) || //cout << "\nShortest path in the sequence given by the ratio on the sequence of edges :";
//cout << endl;
//print(p);

	    (abs_val(vertex3.z()-vertex1.z())>EPS))
	{
	    fprintf(stderr, "\n problem in the computation of the local coordinates \n");
	    exit(0);
        }

        // local coordinates in the plane 
	Vector2D v1, v2, v3;

	v1.x() = vertex1.x(); 	v1.y() = vertex1.y(); 
	v2.x() = vertex2.x(); 	v2.y() = vertex2.y(); 
	v3.x() = vertex3.x(); 	v3.y() = vertex3.y(); 

        // vertices of the common edge between triangles (i-1) and i  
	EdgeType e;
        if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i-1]).n1()) )
	{
                e.vertex_[0] = (tr_seq_unf[i-1].n1());
		e.vertex_[1] = (tr_seq_unf[i-1].n2());
	}
	else
	{ 
	    if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i-1]).n2()))
            {
                e.vertex_[0] = (tr_seq_unf[i-1].n2());
                e.vertex_[1] = (tr_seq_unf[i-1].n3());
            }
	    else
	    { 
		if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i-1]).n3()))
		{
		    e.vertex_[0] = tr_seq_unf[i-1].n3(); 
		    e.vertex_[1] = tr_seq_unf[i-1].n1();
		}
	    }
	}
	edge_seq_unf.push_back(e);

        // vertices of the current triangle corresponding to the common edge
	int ve[2];
	ve[0] = -1; ve[1] = -1;
        Vector2D vect_translation;
	double theta;
        Vector2D v;
	UnfNodeType new_node;
	PrTriangle new_tr;

        if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i]).n1()))
        {
	    ve[0] = 1;
	    ve[1] = 3;
	    new_tr.n1() = edge_seq_unf[i-1].vertex_[0];
	    new_tr.n3() = edge_seq_unf[i-1].vertex_[1];
	    // translation v1 = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_
	    vect_translation = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_ - v1;
	    v3 += vect_translation;
	    v2 += vect_translation;
	    v1 = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_;
	    // rotation v3 = edge_seq_unf[i-1].vertex_[1]
	    theta = 
                angle(v1, v3, v1, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_);
	    rotation2D(v2, theta, v1, v2);
	    v = nodes_unf[third_vertex(tr_seq_unf[i-1], edge_seq_unf[i-1])].n_;
	    if (same_side(v, v2, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_))
	    {
		axial_symmetry(v2, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, 
			       nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_, v2);
	    }
	    new_node.n_ = v2;
	    new_node.origine_ = t.getPrTriangle(tr_seq3d[i]).n2();
	    new_tr.n2() = (int)i+2;
	}
	else
	{
	    if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i]).n2()))
            { 
                ve[0] = 2;
                ve[1] = 1;
                new_tr.n2() = edge_seq_unf[i-1].vertex_[0];
                new_tr.n1() = edge_seq_unf[i-1].vertex_[1];
                // translation v2 = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_
	        vect_translation = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_ - v2;
	        v1 += vect_translation;
	        v3 += vect_translation;
		v2 = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_;
                // rotation v1 =edge_seq_unf[i-1].vertex_[1]
                theta = 
                angle(v2, v1, v2, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_);
	        rotation2D(v3, theta, v2, v3);
                v = nodes_unf[third_vertex(tr_seq_unf[i-1], edge_seq_unf[i-1])].n_;
                if (same_side(v, v3, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_))
                { 
                    axial_symmetry(v3, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, 
				   nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_, v3);
                }
                new_node.n_ = v3;
                new_node.origine_ = t.getPrTriangle(tr_seq3d[i]).n3();
                new_tr.n3() = (int)i + 2;
            }
	    else
	    {
		if (edge_seq3d[i-1].vertex_[0] == (t.getPrTriangle(tr_seq3d[i]).n3()))
		{ 
		    ve[0] = 3;
		    ve[1] = 2;
		    new_tr.n3() = edge_seq_unf[i-1].vertex_[0];
		    new_tr.n2() = edge_seq_unf[i-1].vertex_[1];
		    // translation v3 = edge_seq_unf[i-1].vertex_[0].n_
		    vect_translation = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_ - v3;
		    v2 += vect_translation;
		    v1 += vect_translation;
		    v3 = nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_;
		    // rotation v2 =edge_seq_unf[i-1].vertex_[1]
		    theta = 
			angle(v3, v2, v3, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_);
		    rotation2D(v1, theta, v3, v1);
		    v = nodes_unf[third_vertex(tr_seq_unf[i-1], edge_seq_unf[i-1])].n_;
		    if (same_side(v, v1, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_))
		    {
			axial_symmetry(v1, nodes_unf[edge_seq_unf[i-1].vertex_[0]].n_, 
				       nodes_unf[edge_seq_unf[i-1].vertex_[1]].n_, v1);
		    }
		    new_node.n_ = v1;
	            new_node.origine_ = t.getPrTriangle(tr_seq3d[i]).n1();
		    new_tr.n1() = (int)i + 2;
		}
		else
		{ 
		    fprintf(stderr, "\n vertex identification with the common edge impossible");
		    exit(0);
		}
	    }
        }
	nodes_unf.push_back(new_node);
	tr_seq_unf.push_back(new_tr);
    }
    return;
}


//----------------------------------------------------------------------------
void unfolding_first_triangle( const PrTriangulation_OP& t, 
			       const vector<int>& tr_seq, 
			       vector<PrTriangle>& tr_seq_unf,
			       vector<UnfNodeType>& nodes_unf)
//----------------------------------------------------------------------------
/*
  first triangle unfolded with 
  1st vertex = (0,0), 2nd on the x-axis, 3rd in the y>0 half plane 
 */
{
    if (tr_seq.size()<1)
        return;

    // local frame of this triangle
    Vector3D U, V, W;
    local_frame(t.get3dNode((t.getPrTriangle(tr_seq[0])).n1()), 
	        t.get3dNode((t.getPrTriangle(tr_seq[0])).n2()), 
		t.get3dNode((t.getPrTriangle(tr_seq[0])).n3()),
	        U, V, W);

    // local coordinates in this triangle
    Vector3D vertex0, vertex1, vertex2;
    vertex0 = local_coordinates(t.get3dNode((t.getPrTriangle(tr_seq[0])).n1()), U, V, W);
    vertex1 = local_coordinates(t.get3dNode((t.getPrTriangle(tr_seq[0])).n2()), U, V, W);
    vertex2 = local_coordinates(t.get3dNode((t.getPrTriangle(tr_seq[0])).n3()), U, V, W);
    
    // check that the z-coordinate is zero
    if ((abs_val(vertex1.z()-vertex2.z())>EPS) || 
	    (abs_val(vertex2.z()-vertex0.z())>EPS) || 
	    (abs_val(vertex0.z()-vertex1.z())>EPS))
    {
        fprintf(stderr, "\n problem in the computation of the local coordinates \n");
        exit(0);
    }

    // local coordinates in the plane 
    Vector2D v0, v1, v2;
 
    v0.x() = vertex0.x(); 	v0.y() = vertex0.y(); 
    v1.x() = vertex1.x(); 	v1.y() = vertex1.y(); 
    v2.x() = vertex2.x(); 	v2.y() = vertex2.y(); 
	
    // triangle translation (of vector -v0) such that v0 is in (0,0) 
    v1 -= v0;
    v2 -= v0;
    v0 = Vector2D(0, 0);

    // rotation such that v1 is on the x-axis
    Vector2D O = Vector2D(0, 0);
    Vector2D U1 = Vector2D(1,0);
    double theta;
    theta = angle(O, v1, O, U1);
    rotation2D(v1, theta, O, v1);
    rotation2D(v2, theta, O, v2);

    if (abs_val(v1.y())>EPS)
    {
        fprintf(stderr,"\n vertex v1 is not on the x-axis");
	exit(0);
    }

    // record in the tables
    UnfNodeType new_node;
    new_node.n_ = v0;
    new_node.origine_ = t.getPrTriangle(tr_seq[0]).n1();
    nodes_unf.push_back(new_node);
    new_node.n_ = v1;
    new_node.origine_ = t.getPrTriangle(tr_seq[0]).n2();
    nodes_unf.push_back(new_node);
    new_node.n_ = v2;
    new_node.origine_ = t.getPrTriangle(tr_seq[0]).n3();
    nodes_unf.push_back(new_node);
    
    tr_seq_unf.push_back(PrTriangle(0,1,2, -1, -1, -1));

    return;
}

// --------------------------------------------------------------------
void unfolding_vertex(
    const PrTriangulation_OP& t, 
    const vector<int>& tr_seq3d, vector<PrTriangle>& tr_seq_unf,
    const int vertex, int& vertex_unf)
// --------------------------------------------------------------------
/*
  associates the unfolded node to a vertex of the sequence
  (for example for the source and destination vertices)
 */
{
    if (tr_seq3d.size()<1)
        return;

    int i, i_max;
    i = -1;
    i_max = (int)(tr_seq3d.size()-1);
    PrTriangle tr;
    
    // look for the triangle containing the vertex
    if ((t.getPrTriangle(tr_seq3d[i_max])).isVertex(vertex))
    {
	i = i_max;
	tr = t.getPrTriangle(tr_seq3d[i_max]);
    }
    else
    {
	i = 0;
	while (( !(t.getPrTriangle(tr_seq3d[i]).isVertex(vertex)) )
	       && (i < i_max ))
	{
	    i++;
	}
	if (t.getPrTriangle(tr_seq3d[i]).isVertex(vertex))
	{
	    tr = t.getPrTriangle(tr_seq3d[i]);
	}
	else
	{
	    cout << "This vertex does not belong to the triangle sequence and ";
	    cout <<     "can therefore not be unfolded"; 
	    cout << endl;
	return;
	}    
    }
    
    if (vertex == (tr.n1())) 
    {
	vertex_unf = tr_seq_unf[i].n1();
    }
    else
    {
	if (vertex == (tr.n2()))
	{
	    vertex_unf = tr_seq_unf[i].n2();
	}
	else
	{
	    if (vertex == (tr.n3()))
	    {
		vertex_unf = tr_seq_unf[i].n3();
	    }
	    else 
	    {
		fprintf(stderr, "\n problem : this vertex does not belong to the triangle");
		cout << endl;
	    }
	}
    }
    return;
}

// --------------------------------------------------------------------
void local_frame( 
    Vector3D a, Vector3D b, Vector3D c,  
    Vector3D& U, Vector3D& V, Vector3D& W)
// --------------------------------------------------------------------
/* local frame of triangle ABC */
{
    W = vector_W(a, b, c);
    U = vector_U(W);
    V = vector_V(W, U);	
    return;
}

// --------------------------------------------------------------------
Vector3D vector_W(
    Vector3D a, Vector3D b, Vector3D c)
// --------------------------------------------------------------------
/* axis and normal vector in the current triangle plane for the local frame */
{	
    /* normal vector to triangle */
    Vector3D vect;
    Vector3D o = Vector3D(0,0,0);

    vect = (b-a).cross(c-a);
	
    double norm_vect;
    norm_vect = dist3D(o,vect);

    if (norm_vect < EPS)
    {
	vect = (c-b).cross(a-b);
	norm_vect = dist3D(o,vect);
	if (norm_vect < EPS)
	{
	    vect = (a-c).cross(b-c);
	    norm_vect = dist3D(o,vect);
	    if (norm_vect < EPS)
	    {
		fprintf(stderr, "\nthe normal vector W in the local frame"
			"has a length under the epsilon tolerance\n");
		exit(0);
	    }
			
	}
    }
    return(vect/norm_vect);
}

// --------------------------------------------------------------------
Vector3D vector_U(const Vector3D W)
// --------------------------------------------------------------------
/* first axis/vector in the current triangle plane for the local frame */
{
    Vector3D v;
    Vector3D o = Vector3D(0,0,0);
    double norm_v;
    
    if (abs_val(W.x()) <= abs_val(W.y()) )
    {
	if (abs_val(W.x()) <= abs_val(W.z()))
	{ // x-coordinate is minimal
	    v = Vector3D(0, - W.z(), W.y());
	}
	else
	{ // z-coordinate is minimal
	    v = Vector3D(- W.y(), W.x(), 0);
	}
    }
    else
    {
	if (abs_val(W.y()) <= abs_val(W.z()) )
	{ // y-coordinate is minimal
	    v = Vector3D(- W.z(), 0, W.x());
	}
	else
	{ // z-coordinate is minimal 
	    v = Vector3D(-W.y(), W.x(), 0);
	}
    }
    
    norm_v = dist3D(o,v);
    if ( norm_v < EPS )
    {
	fprintf(stderr, "\nthe first vector U in the plane "
		" of the triangle for the local frame has a length"
		" under the epsilon allowed tolerance \n");
	exit(0);
    }
    
    return(v/norm_v);
}

// --------------------------------------------------------------------
Vector3D vector_V(const Vector3D W, const Vector3D U)
// --------------------------------------------------------------------
/* second axis/vector in the current triangle plane for the local frame */
{
    return(W.cross(U));
}

// --------------------------------------------------------------------
Vector3D local_coordinates(const Vector3D a, 
	 const Vector3D& U, const Vector3D& V, const Vector3D& W)
// --------------------------------------------------------------------
/* local coordinates of a vector A in an axis frame given by (U, V, W) */
{
    return(Vector3D(a*U, a*V, a*W));
}

// --------------------------------------------------------------------
double angle(const Vector2D& a, const Vector2D& b, 
	     const Vector2D& c, const Vector2D& d)
// --------------------------------------------------------------------
/* angle between two vectors AB and CD in the plane */
{
    double theta1, theta2;
	
    theta1 = atan2(b.y() - a.y(), b.x() - a.x());
    theta2 = atan2(d.y() - c.y(), d.x() - c.x()); 
    return(theta2 - theta1);
}

// --------------------------------------------------------------------
void rotation2D(const Vector2D a, const double theta, 
		const Vector2D c, Vector2D& vect)
// --------------------------------------------------------------------
{
    vect.x() = c.x() + cos(theta)*(a.x() - c.x()) - sin(theta)*(a.y() - c.y()); 
    vect.y() = c.y() + sin(theta)*(a.x() - c.x()) + cos(theta)*(a.y() - c.y()); 
    return;
}

// --------------------------------------------------------------------
void axial_symmetry(const Vector2D c, const Vector2D a, 
		    const Vector2D b, Vector2D& v)
// --------------------------------------------------------------------
{
    double theta;
    theta = 2*angle(a, c, a, b);
    rotation2D(c, theta, a, v);
    return;
}

// --------------------------------------------------------------------
int third_vertex(PrTriangle t, EdgeType e)
// --------------------------------------------------------------------
/* third vertex pf triangle t not belonging to edge e */
{
    if ((t.n1() != e.vertex_[0]) && (t.n1() != e.vertex_[1]))
	return (t.n1());
    if ((t.n2() != e.vertex_[0]) && (t.n2() != e.vertex_[1]))
	return (t.n2());
    if ((t.n3() != e.vertex_[0]) && (t.n3() != e.vertex_[1]))
	return (t.n3());
    fprintf(stderr, "\n this triangle is degenerated, two of its vertices are the same");
    cout << endl;
    return (0);
}

// --------------------------------------------------------------------
int same_side(const Vector2D pt1, const Vector2D pt2, 
	      const Vector2D ptA, const Vector2D ptB)
// --------------------------------------------------------------------
/*  check if pt1 and pt2 are both on the same side of [ptA ptB]
 */
{
	double prod1, prod2;
	prod1 = vect_prod2D(ptA, ptB, ptA, pt1);
	prod2 = vect_prod2D(ptA, ptB, ptA, pt2);
			
	if (prod1*prod2 < 0) 	
		return(0);
	else 	return(1);
}

// --------------------------------------------------------------------
double vect_prod2D(Vector2D a, Vector2D b, Vector2D c, Vector2D d)
// --------------------------------------------------------------------
{
    return((b.x() - a.x())*(d.y() - c.y()) - (b.y() - a.y())*(d.x() -c.x()));
}

//----------------------------------------------------------------------------
void printUnfoldedSequence(std::ofstream& os, 
			   vector<UnfNodeType>& nodes_unf, vector<PrTriangle>& tr_seq_unf)
//-----------------------------------------------------------------------------
{  
    std::cout << "Writing file for the unfolded triangle sequence visualisation with gnuplot ...";
    std::cout << std::endl;

    os.precision(5);
    std::vector<PrTriangle>::iterator i;
    for (i = tr_seq_unf.begin(); i != tr_seq_unf.end(); i++)
    {
	if (abs_val((nodes_unf[i->n1()].n_).x())<EPS)
	    (nodes_unf[i->n1()].n_).x() = 0;
	if (abs_val((nodes_unf[i->n2()].n_).x())<EPS)
	    (nodes_unf[i->n2()].n_).x() = 0;
	if (abs_val((nodes_unf[i->n3()].n_).x())<EPS)
	    (nodes_unf[i->n3()].n_).x() = 0;
	if (abs_val((nodes_unf[i->n1()].n_).y())<EPS)
	    (nodes_unf[i->n1()].n_).y() = 0;
	if (abs_val((nodes_unf[i->n2()].n_).y())<EPS)
	    (nodes_unf[i->n2()].n_).y() = 0;
	if (abs_val((nodes_unf[i->n3()].n_).y())<EPS)
	    (nodes_unf[i->n3()].n_).y() = 0;

	os << (nodes_unf[i->n1()].n_).x() << " " << (nodes_unf[i->n1()].n_).y() << endl;
	os << (nodes_unf[i->n2()].n_).x() << " " << (nodes_unf[i->n2()].n_).y() << endl;
	os << endl;
	os << (nodes_unf[i->n2()].n_).x() << " " << (nodes_unf[i->n2()].n_).y() << endl;
	os << (nodes_unf[i->n3()].n_).x() << " " << (nodes_unf[i->n3()].n_).y() << endl;
	os << endl;
 	os << (nodes_unf[i->n3()].n_).x() << " " << (nodes_unf[i->n3()].n_).y() << endl;
	os << (nodes_unf[i->n1()].n_).x() << " " << (nodes_unf[i->n1()].n_).y() << endl;
	os << endl;
   }
}

//----------------------------------------------------------------------------
void printUnfoldedPath(std::ofstream& os, vector<Vector2D> path)
//----------------------------------------------------------------------------
{  
    std::cout <<
      "Writing file for the unfolded path visualisation with gnuplot ...";
    std::cout << std::endl;

    os.precision(5);
    std::vector<Vector2D>::iterator i, j;
    i = path.begin();
    j = path.begin(); j++;
  
    if (abs_val(i->x())<EPS)
	i->x() = 0;
    if (abs_val(i->y())<EPS)
	i->y() = 0;

    while (j != path.end())
    {
	if (abs_val(j->x())<EPS)
	    j->x() = 0;
	if (abs_val(j->y())<EPS)
	    j->y() = 0;

	os << i->x() << " " << i->y() << endl;
	os << j->x() << " " << j->y() << endl;
	os << endl;
	i++; j++;
   }
}

//----------------------------------------------------------------------------
void printUnfoldedPath(std::ofstream& os, 
		       vector<UnfNodeType>& nodes_unf,     
		       list<int>& list_pivots_unf)
//----------------------------------------------------------------------------
/*
 */
{  
    std::cout <<
      "Writing file for the unfolded path visualisation with gnuplot ...";
    std::cout << std::endl;

    os.precision(5);
    std::list<int>::iterator i;
    i = list_pivots_unf.begin();
    std::list<int>::iterator j;
    j = list_pivots_unf.begin(); j++;
   
    if (abs_val((nodes_unf[*i].n_).x())<EPS)
	(nodes_unf[*i].n_).x() = 0;
    if (abs_val((nodes_unf[*i].n_).y())<EPS)
	(nodes_unf[*i].n_).y() = 0;

    while (j != list_pivots_unf.end())
    {
	if (abs_val((nodes_unf[*j].n_).x())<EPS)
	    (nodes_unf[*j].n_).x() = 0;
	if (abs_val((nodes_unf[*j].n_).y())<EPS)
	    (nodes_unf[*j].n_).y() = 0;

	os << (nodes_unf[*i].n_).x() << " " << (nodes_unf[*i].n_).y() << endl;
	os << (nodes_unf[*j].n_).x() << " " << (nodes_unf[*j].n_).y() << endl;
	os << endl;
	i++; j++;
   }
}

//----------------------------------------------------------------------------
void printUnfoldedObjects(std::ofstream& os, 
		       vector<UnfNodeType>& nodes_unf, vector<PrTriangle>& tr_seq_unf,     
		       list<int>& list_pivots_unf)
//----------------------------------------------------------------------------
{  
    // cout << "\nWriting file for the unfolded objects visualisation with gnuplot ...";

    os.precision(5);
    std::vector<PrTriangle>::iterator it;
    for (it = tr_seq_unf.begin(); it != tr_seq_unf.end(); it++)
    {
	if (abs_val((nodes_unf[it->n1()].n_).x())<EPS)
	    (nodes_unf[it->n1()].n_).x() = 0;
	if (abs_val((nodes_unf[it->n2()].n_).x())<EPS)
	    (nodes_unf[it->n2()].n_).x() = 0;
	if (abs_val((nodes_unf[it->n3()].n_).x())<EPS)
	    (nodes_unf[it->n3()].n_).x() = 0;
	if (abs_val((nodes_unf[it->n1()].n_).y())<EPS)
	    (nodes_unf[it->n1()].n_).y() = 0;
	if (abs_val((nodes_unf[it->n2()].n_).y())<EPS)
	    (nodes_unf[it->n2()].n_).y() = 0;
	if (abs_val((nodes_unf[it->n3()].n_).y())<EPS)
	    (nodes_unf[it->n3()].n_).y() = 0;

	os << (nodes_unf[it->n1()].n_).x() << " " << (nodes_unf[it->n1()].n_).y() << endl;
	os << (nodes_unf[it->n2()].n_).x() << " " << (nodes_unf[it->n2()].n_).y() << endl;
	os << endl;
	os << (nodes_unf[it->n2()].n_).x() << " " << (nodes_unf[it->n2()].n_).y() << endl;
	os << (nodes_unf[it->n3()].n_).x() << " " << (nodes_unf[it->n3()].n_).y() << endl;
	os << endl;
 	os << (nodes_unf[it->n3()].n_).x() << " " << (nodes_unf[it->n3()].n_).y() << endl;
	os << (nodes_unf[it->n1()].n_).x() << " " << (nodes_unf[it->n1()].n_).y() << endl;
	os << endl;
   }

    std::list<int>::iterator i;
    i = list_pivots_unf.begin();
    std::list<int>::iterator j;
    j = list_pivots_unf.begin(); j++;
   
    if (abs_val((nodes_unf[*i].n_).x())<EPS)
	(nodes_unf[*i].n_).x() = 0;
    if (abs_val((nodes_unf[*i].n_).y())<EPS)
	(nodes_unf[*i].n_).y() = 0;

    while (j != list_pivots_unf.end())
    {
	if (abs_val((nodes_unf[*j].n_).x())<EPS)
	    (nodes_unf[*j].n_).x() = 0;
	if (abs_val((nodes_unf[*j].n_).y())<EPS)
	    (nodes_unf[*j].n_).y() = 0;

	os << (nodes_unf[*i].n_).x() << " " << (nodes_unf[*i].n_).y() << endl;
	os << (nodes_unf[*j].n_).x() << " " << (nodes_unf[*j].n_).y() << endl;
	os << endl;
	i++; j++;
   }
}


//----------------------------------------------------------------------------
void print_edges_lengths(vector<PrTriangle> tr, vector<UnfNodeType> nodes_unf)
//----------------------------------------------------------------------------
/*
  print the lengths of the edges of the triangles in the unfolded sequence
 */
{
    std::vector<PrTriangle>::iterator iter;
    int i = 0;
    Vector2D a, b, c;

    cout << "\nedges lengths in the unfolded sequence : " << endl;
    for (iter = tr.begin(); iter != tr.end(); iter++)
    {
	a = nodes_unf[iter->n1()].n_;
	b = nodes_unf[iter->n2()].n_;
	c = nodes_unf[iter->n3()].n_;

	i++;
	cout << "triangle " << i << " : ";
	cout << dist2D(a, b) << " - ";
	cout << dist2D(b, c) << " - ";
	cout << dist2D(c, a) << endl;
    }
    cout << endl;
    
}

//----------------------------------------------------------------------------
void print_triangle_sequence(vector<PrTriangle> tr, vector<UnfNodeType> nodes_unf)
//----------------------------------------------------------------------------
{
    std::vector<PrTriangle>::iterator iter;
    int i = 0;
    int a, b, c;

    cout << "\ntriangles in the unfolded sequence and corresponding 3d vertices : ";
    cout << endl;
    for (iter = tr.begin(); iter != tr.end(); iter++)
    {
	a = nodes_unf[iter->n1()].origine_;
	b = nodes_unf[iter->n2()].origine_;
	c = nodes_unf[iter->n3()].origine_;

	i++;
	cout << "triangle " << i << " : ";
	cout << a << " - " << b << " - " << c << endl;
    }
    cout << endl;
    
}

//----------------------------------------------------------------------------
void print_tr_vertices(vector<int>& tr, PrTriangulation_OP& t)
//----------------------------------------------------------------------------
{
    std::vector<int>::iterator iter;
    int i = 0;
    int a, b, c;

    std::cout << "\ntriangles in the 3d sequence and  3d vertices : ";
    std::cout << std::endl;
    for (iter = tr.begin(); iter != tr.end(); iter++)
    {
	a = t.getPrTriangle(*iter).n1();
	b = t.getPrTriangle(*iter).n2();
	c = t.getPrTriangle(*iter).n3();

	i++;
	std::cout << "triangle " << i << " : ";
	std::cout << a << " - " << t.get3dNode(a);
	std::cout << b << " - " << t.get3dNode(b);
	std::cout << c << " - " << t.get3dNode(c) << std::endl;
    }
    std::cout << std::endl;
    
}
