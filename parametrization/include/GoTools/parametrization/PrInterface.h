/********************************************************************
 FILENAME    : PrGeomview.h
 AUTHOR      : Valerie PHAM-TRONG, SINTEF
 DATE        : Mai 2002
 DESCRIPTION : 
 CHANGE LOG  :
*********************************************************************/

#ifndef _PRINTERFACE_H
#define _PRINTERFACE_H


#include <iomanip>
#include <stdlib.h>
//#include <stdio.h>
#include <malloc.h>
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










