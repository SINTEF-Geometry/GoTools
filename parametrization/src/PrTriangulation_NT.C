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

#include "GoTools/parametrization/PrTriangulation_NT.h"
using std::cerr;
using std::cout;
using std::endl;

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrTriangulation_NT::PrTriangulation_NT(PrTriangulation_OP& t)
//-----------------------------------------------------------------------------
//   Construct the PrTriangulation_NT from a PrTriangulation_OP,
//   using just one level in the hierarchy.
//   This can later be refined by refine().
{
  int npts = t.getNumNodes();
  node_.resize(npts);

  int j;
  for(j=0; j< npts; j++)
  {
    PrNode& node = t.getPrNode(j);
    node_[j].init(node.x(),node.y(),node.z(),node.u(),node.v(),0);
    node_[j].addTrianglePtr(node.tr());
  }

  vector<PrTriangle>& nt = t.getTriangleArray();
  triang_.push_back(new PrLevelTriangulation_OP(&node_,nt,0));
}

/*
//-----------------------------------------------------------------------------
PrTriangulation_NT::PrTriangulation_NT(PrParamTriangulation& t, int level)
//-----------------------------------------------------------------------------
//   Construct the PrTriangulation_NT from a PrTriangulation_OP,
//   using just one level in the hierarchy.
//   This can later be refined by refine().
{
  PrTriangulation_NT(t.basemesh_);
  for(int k=1; k<= level; k++)
  {
    refine();
  }
}
*/

//-----------------------------------------------------------------------------
void PrTriangulation_NT::refine()
//-----------------------------------------------------------------------------
//   Refine the nested triangulation s by one level
//   (new nodes ar placed at mid points of their parents
//    and triangles are split into four congruent subtriangles).
{
    int level = (int)triang_.size() - 1;
  PrLevelTriangulation_OP* lt = triang_[level];
  int numTriangles = lt->findNumFaces();
  int numNodes = (int)node_.size();

  vector<PrTriangle> nt(4*numTriangles);
  int i;
  for(i=0; i< numTriangles; i++)
  {
    PrTriangle& t0 = lt->getPrTriangle(i);

    int tn[3],vn[3];
    tn[0] = t0.t1();
    tn[1] = t0.t2();
    tn[2] = t0.t3();
    vn[0] = t0.n1();
    vn[1] = t0.n2();
    vn[2] = t0.n3();
    int v_index[3];

    int j;
    for(j=0; j<3; j++)
    {
      if(tn[j] > i || tn[j] == -1)
      {
        int jj = (j+1) % 3;
        int jjj = (j+2) % 3;
        Vector3D vec = 0.5 * (lt->get3dNode(vn[jj])+lt->get3dNode(vn[jjj]));
        node_.push_back(PrNestedNode(vec.x(),vec.y(),vec.z(),0.0,0.0,level+1));
        node_[numNodes].addTrianglePtr(4*i+jjj+1);
        v_index[j] = numNodes;
        numNodes++;
      }
      else
      {
        if(lt->getPrTriangle(tn[j]).t1() == i)
          v_index[j] = nt[4*tn[j]].n1();
        else if(lt->getPrTriangle(tn[j]).t2() == i)
          v_index[j] = nt[4*tn[j]].n2();
        else v_index[j] = nt[4*tn[j]].n3();
      }
    }

    // Set up the four new triangles
    // Set up the small triangle in the middle

    nt[4*i].init(v_index[0],v_index[1],v_index[2],4*i+1,4*i+2,4*i+3);

    int st1[3],st2[3];
    
    for(j=0; j<3; j++)
    {
      if(tn[j] == -1)
      {
        st1[j] = -1;
        st2[j] = -1;
      }
      else if(i == lt->getPrTriangle(tn[j]).t1())
      {
        st1[j] = 4*tn[j]+3;
        st2[j] = 4*tn[j]+2;
      }
      else if(i == lt->getPrTriangle(tn[j]).t2())
      {
        st1[j] = 4*tn[j]+1;
        st2[j] = 4*tn[j]+3;
      }
      else
      {
        st1[j] = 4*tn[j]+2;
        st2[j] = 4*tn[j]+1;
      }
    }

    nt[4*i+1].init(t0.n1(),v_index[2],v_index[1],4*i,st2[1],st1[2]);
    nt[4*i+2].init(v_index[2],t0.n2(),v_index[0],st1[0],4*i,st2[2]);
    nt[4*i+3].init(v_index[1],v_index[0],t0.n3(),st2[0],st1[1],4*i);

    // add a triangle ptr on the new level to each vertex of the
    // triangle if that is its triangle at the current level

    for(j=0; j<3; j++)
      if( i == lt->getPrNestedNode(vn[j]).tr(level))
        (lt->getPrNestedNode(vn[j])).addTrianglePtr(4*i+j+1);
  }

  triang_.push_back(new PrLevelTriangulation_OP(&node_, nt, level+1));
}

//-----------------------------------------------------------------------------
void PrTriangulation_NT::lift(PrParamTriangulation* pt)
//-----------------------------------------------------------------------------
//   lift the vertices of the nested triangulation to some
//   surface, using the parameterization defined by pt
{
  Vector3D liftedNode;
  Vector3D coarseBC;
  int currentT;
  int bc[3];
  int level;
  int p;

  for (int i=0; i<int(node_.size()); i++) {
    level = node_[i].level();
    currentT = node_[i].tr(level);
    if (currentT<0)
       continue;
    p = 1;
    bc[0] = 0; bc[1] = 0; bc[2] = 0;
    PrTriangle& tri = triang_[level]->getPrTriangle(currentT);

    if (tri.n1() == i)
      bc[0] = 1;
    else if (tri.n2() == i)
      bc[1] = 1;
    else if (tri.n3() == i)
      bc[2] = 1;
    else
      cerr << "LIFT: catastrophic error!" << endl;

    for (int j=1; j<=level; j++) {
      switch (currentT%4) {
      case 0: bc[0] = p-bc[0]; 
              bc[1] = p-bc[1]; 
	      bc[2] = p-bc[2]; 
	      break;
      case 1: bc[0] += p; break;
      case 2: bc[1] += p; break;
      case 3: bc[2] += p; break;
      }
      p *= 2;
      currentT /= 4;
    }

    coarseBC.x() = (double)bc[0]/(double)p; 
    coarseBC.y() = (double)bc[1]/(double)p; 
    coarseBC.z() = (double)bc[2]/(double)p; 

    liftedNode = pt->getSurfPoint(currentT, coarseBC);

    node_[i].x() = liftedNode.x();
    node_[i].y() = liftedNode.y();
    node_[i].z() = liftedNode.z();
  }
}

