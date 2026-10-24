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

#include "GoTools/compositemodel/ParameterizeUtils.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrPrmExperimental.h"
#include "GoTools/parametrization/PrPrmLeastSquare.h"
#include "GoTools/parametrization/PrPrmMeanValue.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/utils/Point.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>  
//#define DEBUG

using std::vector;
using std::make_shared;
using namespace Go;

//===========================================================================
void ParameterizeUtils::parameterizeGraph(shared_ptr<ftPointSet> graph,
					  vector<double>& uv_pars)
//===========================================================================
{
  int num = graph->size();
  ftSamplePoint *first = NULL;
  ftSamplePoint *second = NULL;
  for (int kj=0; kj<num; ++kj)
    {
      ftSamplePoint *curr = (*graph)[kj];
      if (curr->isOnBoundary())
	{
	  first = curr;
	  vector<ftSamplePoint*> adj = curr->getNeighbours();
	  for (size_t kr=0; kr<adj.size(); ++kr)
	    {
	      if (adj[kr]->isOnBoundary())
		{
		  second = adj[kr];
		  break;
		}
	    }
	  if (second)
	    break;
	}
    }

  // first and second may have the wrong order
  if (first)
    graph->setSecond(first);
  if (second)
    graph->setFirst(second);

  // We must make sure that the ftPointSet has the neighbour
  // structure the way the parametrization code expects it.
  graph->orderNeighbours();
  
  // Parameterize
  PrPrmUniform par;
  PrParametrizeBdy bdy;
  shared_ptr<PrOrganizedPoints> op = shared_ptr<PrOrganizedPoints>(graph);

  doParameterize(op, uv_pars);
}

//===========================================================================
void ParameterizeUtils::parameterizeTriang(const double *xyz_points, int nmbp, 
					   const int *triangles, int nmbt,
					   vector<double>& uv_pars)
//===========================================================================
{
  shared_ptr<PrTriangulation_OP> prt = 
    make_shared<PrTriangulation_OP>(xyz_points, nmbp, triangles, nmbt);
  PrParametrizeBdy bdy;
  shared_ptr<PrOrganizedPoints> op = shared_ptr<PrOrganizedPoints>(prt);

  doParameterize(op, uv_pars);
}

//===========================================================================
void ParameterizeUtils::doParameterize(shared_ptr<PrOrganizedPoints>& op,
				       vector<double>& uv_pars)
//===========================================================================
{
  PrParametrizeBdy bdy;
  bdy.attach(op);
  bdy.parametrize();

  vector<int> bd_nodes;
  vector<Vector3D> bd_pnts;
  int bnode = bdy.findBdyNode();
  bd_nodes.push_back(bnode);
  bd_pnts.push_back(op->get3dNode(bnode));
  while (true)
    {
      int bnode2 = bdy.getNextBdyNode(bnode);
      if (bnode2 == bd_nodes[0])
	break;
      bnode = bnode2;
      bd_nodes.push_back(bnode);
      bd_pnts.push_back(op->get3dNode(bnode));
    }
#ifdef DEBUG
  std::ofstream ofb("bd.g2");
  for (size_t kj=0; kj<bd_pnts.size(); ++kj)
    {
      ofb << "400 1 0 4 0 0 255 255" << std::endl;
      ofb << "1" << std::endl;
      ofb << bd_pnts[kj] << std::endl;
    }
#endif
  // Recognize corner nodes
  vector<int> cc(4);
  double bd_len = bdy.boundaryLength(bd_nodes[0], bd_nodes[0]);
  bool found = recognizeCornerNodes(bd_nodes, bd_pnts, bd_len, cc);
  if (!found)
    {
      MESSAGE("WARNING: Corners of parameter domain not found");
      //return;
   

      cc[0] = bd_nodes[0];
      cc[1] = bd_nodes[bd_nodes.size()/4];
      cc[2] = bd_nodes[bd_nodes.size()/2];
      cc[3] = bd_nodes[3*bd_nodes.size()/4];
    }
    
  bdy.parametrize(cc[0],cc[1],cc[2],cc[3]);

#ifdef DEBUG
  vector<Vector3D> bd;
  vector<Vector3D> inner;
  int num = op->getNumNodes();
  for (int kj=0; kj<num; ++kj)
    {
      if (op->isBoundary(kj))
	bd.push_back(op->get3dNode(kj));
      else
	inner.push_back(op->get3dNode(kj));
    }
   std::ofstream of2("triangvx.g2");
   (void)of2.precision(15);
   int k2;
   of2 << "400 1 0 4 0 255 0 255"<< std::endl;
   of2 << inner.size() << std::endl;
   for (k2=0; k2<(int)inner.size(); ++k2)
     of2 << inner[k2][0] << " " << inner[k2][1] << " " << inner[k2][2] << std::endl;
   of2 << std::endl;
   of2 << "400 1 0 4  255 0 0 255"<< std::endl;
   of2 << bd.size() << std::endl;
   for (k2=0; k2<(int)bd.size(); ++k2)
     of2 << bd[k2][0] << " " << bd[k2][1] << " " << bd[k2][2] << std::endl;
   of2 << "400 1 0 4  255 0 255 0 "<< std::endl;
   of2 << "4" << std::endl;
   for (k2=0; k2<4; ++k2)
     of2 << op->get3dNode(cc[k2]) << std::endl;
#endif

  //PrPrmUniform intr;
  PrPrmMeanValue intr;
  //PrPrmLeastSquare intr;
  //PrPrmShpPres intr;
  //PrPrmExperimental intr;
  intr.attach(op);
  intr.parametrize();

  int nmbp = op->getNumNodes();
  uv_pars.reserve(2*nmbp);
  for (int ki=0; ki<nmbp; ++ki)
    {
      uv_pars.push_back(op->getU(ki));
      uv_pars.push_back(op->getV(ki));
    }
}

//===========================================================================
bool ParameterizeUtils::recognizeCornerNodes(vector<int>& bd_nodes,
					     vector<Vector3D>& bd_pnts,
					     double bd_len,
					     vector<int>& cc)
//===========================================================================
{
  if (bd_nodes.size() < 4)
    return false;  // Do not make a suggestion for a degenerate surface
#ifdef DEBUG
  std::ofstream ofbd("boundary_nodes.g2");
  ofbd << "400 1 0 4 255 0 0 255" << std::endl;
  ofbd << bd_pnts.size() << std::endl;
  for (size_t kh=0; kh<bd_pnts.size(); ++kh)
    ofbd << bd_pnts[kh] << std::endl;

  ofbd << "410 1 0 4 50 0 200 255" << std::endl;
   ofbd << bd_pnts.size() << std::endl;
  for (size_t kh=0; kh<bd_pnts.size(); ++kh)
    ofbd << bd_pnts[kh] << " " << bd_pnts[(kh+1)%bd_pnts.size()] << std::endl;
 
#endif
  // Compute angles between consequtive triangle edges at the boundary
  int size = (int)bd_nodes.size();
  int nmb = 1; //(bd_nodes.size() > 40) ? 2 : 1;
  int ki, kj, kr;
  vector<double> bd_ang;
  //vector<double> bd_ang0;
  //for (ki=nmb; ki<size; ++ki)
   for (ki=0; ki<size; ++ki)
   {
      kj = (ki+1)%((int)bd_pnts.size());
      kr = (ki == 0) ? (int)bd_pnts.size()-1 : ki-1;
      Vector3D vec1 = bd_pnts[kj] - bd_pnts[ki];
      Vector3D vec2 = bd_pnts[ki] - bd_pnts[kr];
      // for (kr=0; kr<nmb; ++kr)
      // 	{
      // 	  vec1 += (bd_pnts[ki-kr] - bd_pnts[ki-kr-1]);
      // 	  vec2 += (bd_pnts[(ki+kr+1)%size] - bd_pnts[(ki+kr)%size]);
      // 	}
      // vec1 /= (double)nmb;
      // vec2 /= (double)nmb;

      // // TEST
      // vec1[2] = 0.0;
      // vec2[2] = 0.0;
      // // END TEST

      double angle = vec1.angle(vec2);
      angle = std::min(angle, M_PI-angle);
      bd_ang.push_back(angle);

      // vec1[2] = 0.0;
      // vec2[2] = 0.0;
      // angle = vec1.angle(vec2);
      // bd_ang0.push_back(angle);
    }
  // for (kr=0; kr<nmb; ++kr)
  //   {
  //     bd_ang.insert(bd_ang.begin(), bd_ang[bd_ang.size()-1]);
  //     bd_ang0.insert(bd_ang0.begin(), bd_ang0[bd_ang0.size()-1]);
  //     bd_ang.pop_back();
  //     bd_ang0.pop_back();
  //   }

  // @@@ VSK, 0214. Sort corner angles after the sum of the angles and the
  // angles projected onto the xy-plane and uses the nodes with the 4 largest
  // angles as corners. This is probably a too simple solution that needs to
  // be revised after getting some experience with the functionality
  vector<int> perm(bd_nodes.size());
  for (ki=0; ki<(int)bd_nodes.size(); ++ki)
    perm[ki] = ki;

  for (ki=0; ki<(int)bd_nodes.size(); ++ki)
    {
      double ang2 = bd_ang[perm[ki]]; // + bd_ang0[perm[ki]];
      for (kj=ki+1; kj<(int)bd_nodes.size(); ++kj)
	{
	  double ang3 = bd_ang[perm[kj]]; // + bd_ang0[perm[kj]];
	  if (ang3 > ang2)
	    {
	      std::swap(perm[ki], perm[kj]);
	      ang2 = ang3;
	    }
	}
    }
      
  // Dismiss very close corners
  double threshold = 0.05*bd_len;
  int nmb_perm = (int)perm.size();
  for (ki=0; ki<4; ki++)
    {
      for (kj=ki+1; kj<nmb_perm; )
	{
	  double len = bd_pnts[perm[ki]].dist(bd_pnts[perm[kj]]);
	  if (len < threshold)
	    {
	      perm.push_back(perm[kj]);
	      perm.erase(perm.begin()+kj);
	      nmb_perm--;
	    }
	  else
	    ++kj;
	}
    }
  
  // The corner indices reflects the boundary node array
  cc.resize(4);
  for (ki=0; ki<4; ++ki)
    cc[ki] = perm[ki];

  std::sort(cc.begin(), cc.end());
  for (ki=0; ki<4; ++ki)
    cc[ki] = bd_nodes[cc[ki]];

  return true;
}

//===========================================================================
void ParameterizeUtils::getPointsAndNormal(const double *xyz_points, int nmbp, 
					   const int *triangles, int nmbt,
					   vector<double>& point_and_normal)
//===========================================================================
{
  shared_ptr<PrTriangulation_OP> prt = 
    make_shared<PrTriangulation_OP>(xyz_points, nmbp, triangles, nmbt);

  int num = prt->getNumNodes();
  for (int ki=0; ki<num; ++ki)
    {
      Vector3D node = prt->get3dNode(ki);
      Point pos(node[0], node[1], node[2]);
      vector<int> neighbours;
      prt->getNeighbours(ki, neighbours);
      Vector3D node2 = prt->get3dNode(neighbours[neighbours.size()-1]);
      Point pos2(node2[0], node2[1], node2[2]);
      Point vec1 = pos2 - pos;
      Point norm;
      for (size_t kj=0; kj<neighbours.size(); ++kj)
	{
	  Vector3D node2 = prt->get3dNode(neighbours[kj]);
	  Point pos2(node2[0], node2[1], node2[2]);
	  Point vec2 = pos2 - pos;
	  Point norm2 = vec1.cross(vec2);
	  norm2.normalize();
	  if (kj == 0)
	    norm = norm2;
	  else
	    norm += norm2;
	  vec1 = vec2;
	}
      norm.normalize();
      point_and_normal.insert(point_and_normal.end(), pos.begin(), pos.end());
      point_and_normal.insert(point_and_normal.end(), norm.begin(), norm.end());
    }
}
