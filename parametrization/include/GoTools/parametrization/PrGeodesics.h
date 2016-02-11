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

#ifndef _PRGEODESICS_H
#define _PRGEODESICS_H

#include <stdlib.h>
#include <stdio.h>
#ifndef __APPLE_CC__
#include <malloc.h>
#endif
#include <math.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <list>

#include "GoTools/parametrization/PrPathTriangleSeq.h"
#include "GoTools/parametrization/PrInterface.h"

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"

using Go::Vector3D;
using Go::Vector2D;

#define EPS 1e-8


/*
void GeodesicPT(PrTriangulation_OP& triangulation);

void GeodesicPT(PrTriangulation_OP& triangulation, 
			   int sce_vertex, int dest_vertex);

vector<GoVector3D> GeodesicPT(PrTriangulation_OP& triangulation, 
			   int sce_vertex, int dest_vertex);
vector<GoVector3D> GeodesicPT1(PrTriangulation_OP& triangulation, 
			   int sce_vertex, int dest_vertex);
*/

/** Computes a geodesic path (ie locally shortest path)
 * between the vertices indexed by sce_vertex and dest_vertex.
 * iterativ method consisting in :
 * 1. Initialisation of a path with the Dijkstra method which computes a
 * shortest path following the edges of the triangulation.
 * 2. Improves locally the current path and associated sequence
 *    with an update of the sequence around pivot vertices, ie
 *    vertices where the path is not locally optimal. 
 * properties: 
 * - decreasing length at each step.
 * - finite number of steps.
 *
 * Returns the geodesic path as a polygonal line on the triangulated surface
 * and the triangle sequence crossed by this path. 
 */
void GeodesicPT(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		vector<Vector3D>& current_path, vector<int>& tr_seq3d);
/*
void GeodesicPT(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		vector<Vector3D>& current_path, vector<int>& tr_seq3d,
		std::vector<int> &vert_ind);
void GeodesicPT1(PrTriangulation_OP& triangulation, 
		int sce_vertex, int dest_vertex,
		vector<Vector3D>& current_path, vector<int>& tr_seq3d);
*/

/** Dijkstra computes the shortest path between sce_vertex and dest_vertex, 
 * along the edges of the triangulation. It returns the polygonal solution path
 * given by its vertices, vertices of the triangulation. 
 */
std::vector<int> DijkstraPT(PrTriangulation_OP& triangulation, 
			   int sce_vertex, int dest_vertex);

/// 3D points from their indices in the triangulation
vector<Vector3D> point3D_from_indice(PrTriangulation_OP& t,
				     const vector<int> &ind_path);

/// triangle sequence crossed by a path where the path is given by a list of
/// vert. indices
vector<int> tr_crossed_by_path_vertices(PrTriangulation_OP& t, 
					vector<int> &vert_path);

/** a triangle sequence is modified/updated around a given pivot vertex.
 * returns 
 * - 0 if it is not possible to update 
 * (if pivot belongs to the boundary or pivot does not belong to the sequence)
 * - 1 otherwise
 */
bool update_triangle_sequence_around_pivot(
    PrTriangulation_OP& t, 
    vector<int>& tr_seq, int pivot);

/// triangles containing a vertex, ordered clockwise.
std::list<int> tr_sequence_around_vertex(PrTriangulation_OP& t, 
				    int vertex);

/// triangles containing a vertex, ordered clockwise, starting from first_tr.
std::list<int> tr_sequence_around_vertex(PrTriangulation_OP& t, 
				    int vertex, int first_tr);

bool pivot_in_new_list(std::list<int> &list_pivots, 
		       std::vector<int>::const_iterator &tried_begin,
		       std::vector<int>::const_iterator &tried_end,
		       std::list<int>::iterator& iter_pivot);

/**  decide if there are pivot vertices in the new sequence 
 * set the pivot vertex where the sequence is to update in iter_pivot.
 * compare the pivot lists before and after updating :
 * if equal : incrementation of the pivot.
 * if not : first pivot in the new list.
 */
bool pivot_in_new_list(std::list<int> &list_pivots, 
		       const std::list<int> &list_pivots_prec, 
		       int nb_steps, 
		       const vector<double>& lengths, 
		       std::list<int>::iterator& iter_pivot, int pivot);

/// the treated pivot belongs to the boundary, consider the next one in the list.
bool next_pivot(std::list<int> &list_pivots, 
		std::list<int>::iterator& iter_pivot, int prec_pivot);

/* if updating leaves the shortest path in the sequence unchanged,
 * this function checks that the path has (numerically) decreasing length.
 * if not, comes back to the previous sequence
 */
bool decreasing_length(const int nb_steps, const vector<double>& lengths);

/// distance between two R3 points
double dist3D(Vector3D& a, Vector3D& b);

/// distance between two R2 points 
double dist2D(Vector2D& a, Vector2D& b);
double sqr(double x); 

double length_polygonal_path(const vector<Vector3D> &path);

#endif







