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

#include "GoTools/parametrization/PrParamTriangulation.h"
#include "GoTools/utils/errormacros.h"
#include <set>
/*  using std::set; */
/*  using std::endl; */
/*  using std::cout; */
/*  using std::cerr; */

using namespace std;


//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
open (istream& is)
//-----------------------------------------------------------------------------
{
  int i;
  int nfp, nft, nct;
  is >> nfp >> nft >> nct;

  double* fineXYZ = new double[3*nfp];
  double* fineUV = new double[2*nfp];
  int* fineTri = new int[3*nft];
  int* coarseTri = new int[3*nft];
  baseTriangle_.resize(nfp);
  patch_triangles_.resize(nct);

  for (i=0; i<nfp; i++) {
    is >> fineXYZ[3*i] >> fineXYZ[3*i+1] >> fineXYZ[3*i+2];
    is >> fineUV[2*i] >> fineUV[2*i+1];
    is >> baseTriangle_[i];
  }

  for (i=0; i<nft; i++)
    is >> fineTri[3*i] >> fineTri[3*i+1] >> fineTri[3*i+2];

  for (i=0; i<nct; i++) {
    is >> coarseTri[3*i] >> coarseTri[3*i+1] >> coarseTri[3*i+2];
    int n;
    is >> n;
    patch_triangles_[i].resize(n);
    for (int j=0; j<n; j++)
      is >> patch_triangles_[i][j];
  }

  mesh_ = new PrTriangulation_OP(fineXYZ, fineUV, nfp, fineTri, nft);

  map<int,int> old2new;
  map<int,int>::const_iterator p;
  int ncp = 0;
  for (i=0; i<nct; i++) 
    for (int j=0; j<3; j++) {
      p = old2new.find(coarseTri[3*i+j]);
      if (p != old2new.end())
	coarseTri[3*i+j] = p->second;
      else {
	old2new[coarseTri[3*i+j]] = ncp;
	coarseTri[3*i+j] = ncp;
	ncp++;
      }
    }

  double* coarseXYZ = new double[3*ncp];
  for (p=old2new.begin(); p!=old2new.end(); ++p) {
    int oldIdx = p->first;
    int newIdx = p->second;
    coarseXYZ[3*newIdx]   = fineXYZ[3*oldIdx];
    coarseXYZ[3*newIdx+1] = fineXYZ[3*oldIdx+1];
    coarseXYZ[3*newIdx+2] = fineXYZ[3*oldIdx+2];
  }
  
  basemesh_ = new PrTriangulation_OP(coarseXYZ, ncp, coarseTri, nct);
}

//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
save (ostream& os)
//-----------------------------------------------------------------------------
{
  int i;
  int nfp = mesh_->getNumNodes();
  int nft = mesh_->findNumFaces();
  int nct = basemesh_->findNumFaces();

  os << nfp << " " << nft << " " << nct << endl;
  os << endl;

  for (i=0; i<nfp; i++) {
    const PrNode& node = mesh_->getPrNode(i);
    os << node.x() << " " << node.y() << " " << node.z() << " ";
    os << node.u() << " " << node.v() << " " << baseTriangle_[i] << endl;
  }
  os << endl;

  for (i=0; i<nft; i++) {
    const PrTriangle& tri = mesh_->getPrTriangle(i);
    os << tri.n1() << " " << tri.n2() << " " << tri.n3() << endl;
  }
  os << endl;

  for (i=0; i<nct; i++) {
    const PrTriangle& tri = basemesh_->getPrTriangle(i);
    os << tri.n1() << " " << tri.n2() << " " << tri.n3() << " ";
    int n = (int)patch_triangles_[i].size();
    os << n << " ";
    for (int j=0; j<n; j++)
      os << patch_triangles_[i][j] << " ";
    os << endl;
  }
}


//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
attach(PrTriangulation_OP* mesh, PrTriangulation_OP* basemesh) 
//-----------------------------------------------------------------------------
{
  mesh_ = mesh; 
  basemesh_ = basemesh;
  
  baseTriangle_.assign(mesh_->getNumNodes(), -1);
  patch_triangles_.resize(basemesh_->findNumFaces());
}  


//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
makeCorrespondences (int idx, PrSubTriangulation& sub_tri)
//-----------------------------------------------------------------------------
{
//    cerr << " * PrParamTriangulation: setting up correspondences for ";
//    cerr << "coarse triangle " << idx << endl;

  for (int i=0; i<sub_tri.getNumNodes(); i++)
    baseTriangle_[sub_tri.getGlobalIndex(i)] = idx;

  sub_tri.getTriangleIndices(patch_triangles_[idx]);

//    cerr << "     following triangles are inside this one: " << endl;
//    for (i=0; i<patch_triangles_[idx].size(); i++)
//      cerr << patch_triangles_[idx][i] << ", ";
//    cerr << endl;
}

//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
printBaseTriangles ()
//-----------------------------------------------------------------------------
{
  int i;
  for (i=0; i<mesh_->getNumNodes(); i++)
    if (baseTriangle_[i] == -1)
      cerr << i << "->" << baseTriangle_[i] << ", ";

  //return;

  int num_pat = (int)patch_triangles_.size();
  for (i=0; i<num_pat; i++) {
      cout << "Coarse Triangle #" << i << " contains " << (int)patch_triangles_[i].size() << " triangles: ";
      int num_tri = (int)patch_triangles_[i].size();
    for (int j=0; j<num_tri; j++) {
      cout << patch_triangles_[i][j] << ": ";
      const PrTriangle& tri = mesh_->getPrTriangle(patch_triangles_[i][j]);
      double u,v;
      getUV(tri.n1(), i,u,v);
      cout << "(" << 1-u-v << "," << u << "," << v << "), ";
      getUV(tri.n2(), i,u,v);
      cout << "(" << 1-u-v << "," << u << "," << v << "), ";
      getUV(tri.n3(), i,u,v);
      cout << "(" << 1-u-v << "," << u << "," << v << ")" << endl;
    }
    cout << endl;
  }

  vector< set<int> > cTperNode(mesh_->getNumNodes());
  for (i=0; i<num_pat; i++) {
      int num_tri = (int)patch_triangles_[i].size();
    for (int j=0; j<num_tri; j++) {
      const PrTriangle& tri = mesh_->getPrTriangle(patch_triangles_[i][j]);
      cTperNode[tri.n1()].insert(i);
      cTperNode[tri.n2()].insert(i);
      cTperNode[tri.n3()].insert(i);
    }
  }
  for (i=0; i<int(cTperNode.size()); i++) {
    cout << "Node " << i << " is connected to : " << endl;
    set<int>::const_iterator ci;
    for (ci=cTperNode[i].begin(); ci!=cTperNode[i].end(); ++ci) {
      cout << "  " << *ci << " : ";
      double u,v;
      getUV(i, *ci, u,v);
      cout << "(" << 1-u-v << "," << u << "," << v << ")" << endl;	
    }
    cout << endl;
  }
  
  cerr << " * PrParamTriangulation: original vertices correspond to ";
  cerr << " the following base triangles:\n     ";
  for (i=0; i<mesh_->getNumNodes(); i++)
    cerr << i << "->" << baseTriangle_[i] << ", ";
  cerr << endl;
}

//-----------------------------------------------------------------------------
Vector3D 
PrParamTriangulation::
getParamPoint(int i) 
//-----------------------------------------------------------------------------
{
  double u = mesh_->getU(i);
  double v = mesh_->getV(i);
  Vector3D bc (1.0-u-v, u, v);
  
  return evaluator(*basemesh_, baseTriangle_[i], bc);
}


//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
findTriangleContainingPoint(int coarseT, const Vector3D& coarseBC,
			    int& fineT, Vector3D& fineBC)
//-----------------------------------------------------------------------------
{
  const vector<int>& triList = patch_triangles_[coarseT];
  for (size_t i=0; i<triList.size(); i++) {
    if (triangleContainsPoint(triList[i], coarseT, coarseBC, fineBC)) {
      fineT = triList[i];
      return;
    }
  }
  fineT = -1;
}

//-----------------------------------------------------------------------------
void 
PrParamTriangulation::
getUV (int node, int tri, double&u, double&v)
//-----------------------------------------------------------------------------
// returns the barycentric coordinates of node w.r.t. tri
// normally, these are the u/v values obtained by the
// parameterization, but for vertices on edges ...
{
  u = mesh_->getU(node);
  v = mesh_->getV(node);
  int nodeTri = baseTriangle_[node];
  if (tri == nodeTri)
    return;
  
  bool rZero = (1-u-v < 1e-8);
  bool sZero = (u < 1e-8);
  bool tZero = (v < 1e-8);
  
  if ((rZero && (sZero || tZero)) || (sZero && tZero)) {

    // handling vertices of the base mesh...

    const PrTriangle& T1 = basemesh_->getPrTriangle(nodeTri);
    int coarseNode;
    if (!rZero) {
      coarseNode = T1.n1();
    }
    else if (!sZero) {
      coarseNode = T1.n2();
    }
    else if (!tZero) {
      coarseNode = T1.n3();
    }
    else {
      cout << "OOPS! all bc's are Zero" << endl;
      return;
    }
    
    const PrTriangle& T2 = basemesh_->getPrTriangle(tri);
    if (!T2.isVertex(coarseNode)) {
      cout << "OOPS! the common vetex was NOT a common vertex" << endl;
      return;
    }
    
    if (T2.n1() == coarseNode) {
      u = 0.0;
      v = 0.0;
    }
    else if (T2.n2() == coarseNode) {
      u = 1.0;
      v = 0.0;
    }
    else {
      u = 0.0;
      v = 1.0;
    }
    return;
  }
  
  // handling nodes on edges...
  
  const PrTriangle& T1 = basemesh_->getPrTriangle(nodeTri);
  
  int commonNode=-1, oppTri=-1;
  
  double q1,q2;
  if (rZero) {
    commonNode = T1.n2();
    oppTri = T1.t1();
    q1 = u; 
    q2 = v;
  }
  else if (sZero) {
    commonNode = T1.n3();
    oppTri = T1.t2();
    q1 = v; 
    q2 = 1.0-u-v;
  }
  else if (tZero) {
    commonNode = T1.n1();
    oppTri = T1.t3();
    q1 = 1.0-u-v; 
    q2 = u;
  }
  else {
    cout << "OOPS! There was no bc equal to Zero! for node " << node << " / Triangle " << tri << endl;
    return;
  }
  
  if (oppTri != tri) {
    cout << "OOPS! Did not find the coarse Triangle! for node " << node << " / Triangle " << tri << endl;
    return;
  }
  
  const PrTriangle& T2 = basemesh_->getPrTriangle(oppTri);
  
  if (T2.n1() == commonNode) {
    u = 0.0; 
    v = q2;
  }
  else if (T2.n2() == commonNode) {
    u = q1;
    v = 0.0;
  }
  else if (T2.n3() == commonNode) {
    u = q2; 
    v = q1;
  }
  else {
    cout << "OOPS! Triangle found but node lost" << endl;
    return;
  }
  
}

//-----------------------------------------------------------------------------
bool 
PrParamTriangulation::
triangleContainsPoint(int fineTri, int coarseTri,
		      const Vector3D& coarseBC, 
		      Vector3D& fineBC)
//-----------------------------------------------------------------------------
{
  #define DET(a1,a2,b1,b2,c1,c2) (b1*c2+c1*a2+a1*b2-a2*b1-b2*c1-c2*a1)

  const PrTriangle tri = mesh_->getPrTriangle(fineTri);
  
  double r,s,t;

  double u1,v1;
  getUV (tri.n1(), coarseTri, u1, v1);
  
  double u2,v2;
  getUV (tri.n2(), coarseTri, u2, v2);
  
  t = DET (u1,v1, u2,v2, coarseBC[1], coarseBC[2]);
  if (t < -1e-10)
    return false;

  double u3,v3;
  getUV (tri.n3(), coarseTri, u3, v3);
  
  s = DET (u1,v1, coarseBC[1],coarseBC[2], u3,v3);
  if (s < -1e-10)
    return false;

  r = DET (coarseBC[1],coarseBC[2], u2,v2, u3,v3);
  if (r < -1e-10)
    return false;

  double scale = r+s+t;
  fineBC.x() = r/scale;
  fineBC.y() = s/scale;
  fineBC.z() = t/scale;
  
  return true;
}

//-----------------------------------------------------------------------------
Vector3D 
PrParamTriangulation::
getSurfPoint(int coarseTri, const Vector3D& coarseBC) 
//-----------------------------------------------------------------------------
{
  Vector3D fineBC;
  int fineTri;
  findTriangleContainingPoint(coarseTri, coarseBC, fineTri, fineBC);
  ALWAYS_ERROR_IF(fineTri == -1, "couldn't find surrounding triangle!");

  Vector3D surf_point = evaluator(*mesh_, fineTri, fineBC);
  return surf_point;
}

//-----------------------------------------------------------------------------
Vector3D 
PrParamTriangulation::
evaluator(PrTriangulation_OP& mesh, int idx, Vector3D& bc)
//-----------------------------------------------------------------------------
{
  Vector3D result (0.0, 0.0, 0.0);
  
  const PrTriangle& tri = mesh.getPrTriangle(idx);
  
  result += bc[0]*mesh.get3dNode(tri.n1());
  result += bc[1]*mesh.get3dNode(tri.n2());
  result += bc[2]*mesh.get3dNode(tri.n3());
  
  return result;
}

