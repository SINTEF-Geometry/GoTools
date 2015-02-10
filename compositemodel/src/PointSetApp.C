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

#include "GoTools/compositemodel/PointSetApp.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrPlanarGraph_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrOrganizedPoints.h"
#include <fstream>

#define DEBUG

using namespace Go;
using std::vector;

//===========================================================================
void PointSetApp::parameterizeTriang(const double *xyz_points, int nmbp, 
				     const int *triangles, int nmbt,
				     vector<double>& uv_pars)
//===========================================================================
{

  // First create ftPointSet.
  // Insert points
  shared_ptr<ftPointSet> triang(new ftPointSet());
  int ki;
  for (ki=0; ki<nmbp; ++ki)
    {
      Vector3D pos(xyz_points[3*ki], xyz_points[3*ki+1],xyz_points[3*ki+2]);
      shared_ptr<ftSamplePoint> curr_pt(new ftSamplePoint(pos, 0));
      triang->addEntry(curr_pt);
    }

  // Add triangle information
  for (ki=0; ki<nmbt; ++ki)
    {
      ftSamplePoint* pt1 = (*triang)[triangles[3*ki]];
      ftSamplePoint* pt2 = (*triang)[triangles[3*ki+1]];
      ftSamplePoint* pt3 = (*triang)[triangles[3*ki+2]];
      pt1->addNeighbour(pt2);
      pt2->addNeighbour(pt1);
      pt2->addNeighbour(pt3);
      pt3->addNeighbour(pt2);
      pt3->addNeighbour(pt1);
      pt1->addNeighbour(pt3);
    }

#ifdef DEBUG

  std::ofstream of0("triang1.g2");
  of0 << "400 1 0 4 0 0 255 255"<< std::endl;
  of0 << std::endl;
  triang->printXYZEdges(of0);

#endif 

  // Recognize boundary and 4 corner nodes
  vector<int> corner_ix;
  bool found = recognizeBoundary(triang, corner_ix);
  if (!found)
    return;

  // Parameterize points
  // We must make sure that the ftPointSet has the neighbour
  // structure the way the parametrization code expects it.
  triang->checkAndUpdateTriangCorners();
  triang->orderNeighbours();

  // Parameterize
  PrPrmUniform par;
  PrParametrizeBdy bdy;
  shared_ptr<PrOrganizedPoints> op = shared_ptr<PrOrganizedPoints>(triang);
  double umin = 0.0;
  double umax = 1.0;
  double vmin = 0.0;
  double vmax = 1.0;

  try {
    bdy.attach(op);
    bdy.parametrize(corner_ix[0], corner_ix[1], corner_ix[2], corner_ix[3], 
    		    umin, umax, vmin, vmax);
    // bdy.parametrize(corner_ix[3], corner_ix[2], corner_ix[1], corner_ix[0], 
    // 		    umin, umax, vmin, vmax);

    par.attach(op);
    par.setBiCGTolerance(0.00001);
    par.parametrize();
  } catch(...) {
    return;
  }
  
#ifdef DEBUG
  std::ofstream of("param.g2");
  triang->write2D(of);

  std::ofstream of2("triang2.g2");
  of2 << "400 1 0 4 0 0 255 255"<< std::endl;
  of2 << std::endl;
  triang->printXYZEdges(of2);

  vector<Vector3D> bd_nodes;
  vector<Vector3D> inner_nodes;
  int k2;
  for (k2=0; k2<(int)triang->size(); ++k2)
    {
      if ((*triang)[k2]->isOnBoundary())
	bd_nodes.push_back((*triang)[k2]->getPoint());
      else
	inner_nodes.push_back((*triang)[k2]->getPoint());
    }
  
	of2 << "400 1 0 4 255 0 0 255" << std::endl;
	of2 << bd_nodes.size() << std::endl;
	for (k2=0; k2<(int)bd_nodes.size(); ++k2)
	    of2 << bd_nodes[k2][0] << " " << bd_nodes[k2][1] << " " << bd_nodes[k2][2] << std::endl;
	of2 << "400 1 0 4 0 255 0 255" << std::endl;
	of2 << inner_nodes.size() << std::endl;
	for (k2=0; k2<(int)inner_nodes.size(); ++k2)
	    of2 << inner_nodes[k2][0] << " " << inner_nodes[k2][1] << " " << inner_nodes[k2][2] << std::endl;
#endif 

  // Extract parameter values
  int num = triang->size();
  uv_pars.reserve(2*num);
  for (ki=0; ki<num; ++ki)
    {
      Vector2D par = (*triang)[ki]->getPar();
      uv_pars.insert(uv_pars.end(), par.begin(), par.end());
    }
}


//===========================================================================
bool PointSetApp::recognizeBoundary(shared_ptr<ftPointSet>& triang,
				    vector<int>& corner_ix)
//===========================================================================
{
  // Recognize boundary nodes
  // First extract oriented triangle information
  int ki;
  vector<vector<int> > tri;
  triang->getOrientedTriangles(tri);
  vector<int> tri2;
  tri2.reserve(3*tri.size());
  for (ki=0; ki<(int)tri.size(); ++ki)
    tri2.insert(tri2.end(), tri[ki].begin(), tri[ki].end());

  // Use PrTriangulation_OP to identify boundary nodes
  int nb = triang->getNumNodes();
  vector<double> xyz_points;
  xyz_points.reserve(3*nb);
  for (ki=0; ki<nb; ++ki)
    {
      Vector3D pos = (*triang)[ki]->getPoint();
      xyz_points.insert(xyz_points.end(), pos.begin(), pos.end());
    }
  PrTriangulation_OP prtri(&xyz_points[0], nb, &tri2[0], (int)tri.size());
  int num_nodes = prtri.getNumNodes();
  int first_bd = -1;
  for (ki=0; ki<num_nodes; ++ki)
    if (prtri.isBoundary(ki))
      {
	if (first_bd < 0)
	  first_bd = ki;
	(*triang)[ki]->setBoundary(1);
      }
      
  if (first_bd == -1)
    return false;   // Nothing more to do

  // Sort boundary nodes
  vector<int> bd_nodes;
  int cur_bnode = first_bd;
  std::vector<int> neigh;
  do {
    bd_nodes.push_back(cur_bnode);
    prtri.getNeighbours(cur_bnode, neigh);
    cur_bnode = neigh[0];
    
  } while (cur_bnode != first_bd);
  
  // Select 4 corner nodes
  bool found = recognizeCornerNodes(triang, bd_nodes, corner_ix);
  triang->setFirst((*triang)[bd_nodes[corner_ix[0]]]);
  triang->setSecond((*triang)[bd_nodes[corner_ix[0]+1]]);

  // The corner indices are changed to reflect the triangulation nodes
  for (ki=0; ki<4; ++ki)
    corner_ix[ki] = bd_nodes[corner_ix[ki]];
  //std::sort(corner_ix.begin(), corner_ix.end());
   
  return found;
}

//===========================================================================
bool PointSetApp::recognizeCornerNodes(const shared_ptr<ftPointSet>& triang,
				       const vector<int>& bd_nodes,
				       vector<int>& corner_ix)
//===========================================================================
{
  if (bd_nodes.size() < 4)
    return false;  // Do not make a suggestion for a degenerate surface

  // Compute angles between consequtive triangle edges at the boundar
  int ki, kj;
  vector<double> bd_ang;
  vector<double> bd_ang0;
  Vector3D prev, curr, next;
  prev = (*triang)[bd_nodes[bd_nodes.size()-1]]->getPoint();
  curr = (*triang)[bd_nodes[0]]->getPoint();
  for (ki=0; ki<(int)bd_nodes.size(); ++ki)
    {
      kj = (ki+1)%((int)bd_nodes.size());
      next = (*triang)[bd_nodes[kj]]->getPoint();

      Vector3D vec1 = curr - prev;
      Vector3D vec2 = next - curr;
      double angle = vec1.angle(vec2);
      bd_ang.push_back(angle);

      vec1[2] = 0.0;
      vec2[2] = 0.0;
      angle = vec1.angle(vec2);
      bd_ang0.push_back(angle);
      prev = curr;
      curr = next;
    }

  // @@@ VSK, 0214. Sort corner angles after the sum of the angles and the
  // angles projected onto the xy-plane and uses the nodes with the 4 largest
  // angles as corners. This is probably a too simple solution that needs to
  // be revised after getting some experience with the functionality
  vector<int> perm(bd_nodes.size());
  for (ki=0; ki<(int)bd_nodes.size(); ++ki)
    perm[ki] = ki;

  for (ki=0; ki<(int)bd_nodes.size(); ++ki)
    {
      double ang2 = bd_ang[perm[ki]] + bd_ang0[perm[ki]];
      for (kj=ki+1; kj<(int)bd_nodes.size(); ++kj)
	{
	  double ang3 = bd_ang[perm[kj]] + bd_ang0[perm[kj]];
	  if (ang3 > ang2)
	    {
	      std::swap(perm[ki], perm[kj]);
	      ang2 = ang3;
	    }
	}
    }
      
  // The corner indices reflects the boundary node array
  corner_ix.resize(4);
  for (ki=0; ki<4; ++ki)
    corner_ix[ki] = perm[ki];

  std::sort(corner_ix.begin(), corner_ix.end());
  return true;
}
