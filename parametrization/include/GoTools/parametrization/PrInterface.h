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

#ifndef _PRINTERFACE_H
#define _PRINTERFACE_H


#include <iomanip>
#include <stdlib.h>
//#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <math.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <list>

#include "GoTools/parametrization/PrPathTriangleSeq.h"
#include "GoTools/parametrization/PrGeodesics.h"

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"

using Go::Vector3D;
using Go::Vector2D;
using std::list;

void print(vector<int>& v);
void print(vector<double>& v);
void print(vector<Vector3D>& v);
void print(vector<Vector2D>& v);
void print(vector<PrTriangle>& t);

/// Writing file for the triangulation visualisation with geomview
void printGeomviewTriangulation(PrTriangulation_OP& t);
/// Writing file for the triangulation visualisation with geomview.
/// Problem in giving the name function as parameter.
void printGeomviewTriangulation(std::ofstream& os, PrTriangulation_OP& t);

/// Writing file for the triangle sequence visualisation with geomview
void printGeomviewTriangleSequence(PrTriangulation_OP& t, vector<int> tr_seq);
/// Writing file for the triangle sequence visualisation with geomview.
/// Problem in giving the name function as parameter.
void printGeomviewTriangleSequence(std::ofstream& os, PrTriangulation_OP& t, vector<int> tr_seq);

/// Writing file for the shortest path visualisation with geomview
void printGeomviewPath(vector<Vector3D>& path);
/// Writing file for the shortest path visualisation with geomview.
/// Problem in giving the name function as parameter.
void printGeomviewPath(std::ofstream& os, vector<Vector3D>& path);

void print_edges_lengths(vector<int> tr_seq, PrTriangulation_OP& t);
void print_lengths(const vector<double>& lengths);

#endif










