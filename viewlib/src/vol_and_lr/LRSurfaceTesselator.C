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

#include "GoTools/viewlib/vol_and_lr/LRSurfaceTesselator.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/LineCloud.h"
#include <fstream>

//#define DEBUG

using std::vector;
using std::pair;

using namespace Go;



//===========================================================================
LRSurfaceTesselator::~LRSurfaceTesselator()
//===========================================================================
{
}


//===========================================================================
void LRSurfaceTesselator::changeRes(int n, int m)
//===========================================================================
{
    if ((m != m_) || (n != n_)) {
	m_ = m;
	n_ = n;
	tesselate();
    }
}


//===========================================================================
void LRSurfaceTesselator::tesselate()
//===========================================================================
{
  double eps = 1.0e-9;
  const LRSplineSurface *lrsf = dynamic_cast<const LRSplineSurface*>(&surf_);
  if (!lrsf)
    return;   // Not the right surface type

#ifdef DEBUG
  std::cout << "Tesselate LR surface" << std::endl;
#endif
  RectDomain dom = lrsf->containingDomain();
  double umin = lrsf->startparam_u();
  double umax = lrsf->endparam_u();
  double vmin = lrsf->startparam_v();
  double vmax = lrsf->endparam_v();
  
  // Construct mesh of element pointers
  vector<Element2D*> elements;
  lrsf->constructElementMesh(elements);

  Mesh2D lrmesh = lrsf->mesh();
  int dim = lrsf->dimension();
    
  // Get all knot values in the u-direction
  const double* const uknots = lrmesh.knotsBegin(XFIXED);
  const double* const uknots_end = lrmesh.knotsEnd(XFIXED);
  int nmb_knots_u = lrmesh.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = lrmesh.knotsBegin(YFIXED);
  const double* const vknots_end = lrmesh.knotsEnd(YFIXED);
  int nmb_knots_v = lrmesh.numDistinctKnots(YFIXED);
  const double* knotv;

  double udel = (umax - umin)/(double)(n_-1);
  double vdel = (vmax - vmin)/(double)(m_-1);
  double upar = umin;
  double vpar = vmin;
  double tolu = std::max(1.0e-8, 1.0e-8*udel);
  double tolv = std::max(1.0e-8, 1.0e-8*vdel);
  double pdel = 0.2*std::min(udel, vdel);

  // Evaluate regular mesh
  int ki, kj, kr, kh;
  Point pos, normal, par;
  vector<Point> vertex, vertex_norm, vertex_par, vertex_tex;
  vertex.reserve(n_*m_);
  vertex_par.reserve(n_*m_);
  vertex_norm.reserve(n_*m_);
  vertex_tex.reserve(n_*m_);

  double tdiv_u = (double)(n_-1);
  double tdiv_v = (double)(m_-1);
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<m_ && vpar <= (*knotv)+lastv*tolv; ++kr, vpar+=vdel)
	{
	  if (lastv)
	    vpar = std::min(vpar, *knotv);
	  upar = umin;
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<n_ && upar <= (*knotu)+lastu*tolu; ++kh, upar+=udel)
		{
		  if (lastu)
		    upar = std::min(upar, *knotu);
		  lrsf->point(pos, upar, vpar, elem);
		  vertex.push_back(pos);
		  vertex_par.push_back(Point(upar, vpar));
		  if (mesh_->useNormals())
		    {
		      lrsf->normal(normal, upar, vpar, elem);
		      vertex_norm.push_back(normal);
		    }
		  if (mesh_->useTexCoords())
		    vertex_tex.push_back(Point((double)kh/tdiv_u, (double)kr/tdiv_v));
		}
	    }
	}
    }

  // Define preliminary triangles
  int num_tri = 2*(n_-1)*(m_-1);
  vector<unsigned int> tri;
  tri.reserve(3*num_tri);

  for (kr=0; kr<m_-1; ++kr)
    for (kh=0; kh<n_-1; ++kh)
      {
	tri.push_back(kr*n_ + kh);
	tri.push_back(kr*n_ + kh + 1);
	tri.push_back((kr+1)*n_ + kh + 1);
	
	tri.push_back(kr*n_ + kh);
	tri.push_back((kr+1)*n_ + kh + 1);
	tri.push_back((kr+1)*n_ + kh);
      }

  // Insert gaps due to multiple knot lines
  // First identify knot line segments of high multiplicity (deg + 1)
  int deg1 = lrsf->degree(XFIXED);
  int deg2 = lrsf->degree(YFIXED);
  int nknots1 = lrmesh.numDistinctKnots(XFIXED);
  vector<pair<Point, Point> > gaps;
  for (int ka=1; ka<nknots1-1; ++ka)
    {
      double xval = lrmesh.kval(XFIXED, ka);
      vector<pair<int,int> > x_segs = lrmesh.segments(XFIXED, ka, deg1+1);
      for (size_t ki=0; ki<x_segs.size(); ++ki)
	{
	  double yval1 = lrmesh.kval(YFIXED, x_segs[ki].first);
	  double yval2 = lrmesh.kval(YFIXED, x_segs[ki].second);
	  gaps.push_back(std::make_pair(Point(xval,yval1), Point(xval,yval2)));
	}
    }

  int nknots2 = lrmesh.numDistinctKnots(YFIXED);
  for (int ka=1; ka<nknots2-1; ++ka)
    {
      double yval = lrmesh.kval(YFIXED, ka);
      vector<pair<int,int> > y_segs = lrmesh.segments(YFIXED, ka, deg2+1);
      for (size_t ki=0; ki<y_segs.size(); ++ki)
	{
	  double xval1 = lrmesh.kval(XFIXED, y_segs[ki].first);
	  double xval2 = lrmesh.kval(XFIXED, y_segs[ki].second);
	  gaps.push_back(std::make_pair(Point(xval1,yval), Point(xval2,yval)));
	}
    }

#ifdef DEBUG
  std::ofstream of1("par_tri1.g2");
  std::cout << "High multiplicity: " << std::endl;
  for (size_t ki=0; ki<gaps.size(); ++ki)
    {
      std::cout << gaps[ki].first << " " << gaps[ki].second << std::endl;

      of1 << "410 1 0 4 255 0 0 255" << std::endl;
      of1 << "1" << std::endl;
      of1 << gaps[ki].first << " 0.0 " << gaps[ki].second << " 0.0" << std::endl;
    }
  of1 << "410 1 0 0" << std::endl;
  of1 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      of1 << vertex_par[tri[kh]] << " 0.0 " << vertex_par[tri[kh+1]] << " 0.0" << std::endl;
      of1 << vertex_par[tri[kh+1]] << " 0.0 " << vertex_par[tri[kh+2]] << " 0.0" << std::endl;
      of1 << vertex_par[tri[kh+2]] << " 0.0 " << vertex_par[tri[kh]] << " 0.0" << std::endl;
    }
  
  std::ofstream ofm1("init_mesh.g2");
  ofm1 << "410 1 0 0" << std::endl;
  ofm1 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      ofm1 << vertex[tri[kh]] << "  " << vertex[tri[kh+1]] << " " << std::endl;
      ofm1 << vertex[tri[kh+1]] << "  " << vertex[tri[kh+2]] << " " << std::endl;
      ofm1 << vertex[tri[kh+2]] << "  " << vertex[tri[kh]] << " " << std::endl;
    }
#endif

  vector<Point> joints;
  if (gaps.size() > 1)
    {
      // Split gaps at joints
      for (size_t ki=0; ki<gaps.size(); ++ki)
	{
	  int i1 = (gaps[ki].second[0] -  gaps[ki].first[0] >
		    gaps[ki].second[1] -  gaps[ki].first[1]) ? 0 : 1;
	  double a1 = std::min(gaps[ki].first[i1],gaps[ki].second[i1]);
	  double a2 = std::max(gaps[ki].first[i1],gaps[ki].second[i1]);
	  double b = gaps[ki].first[1-i1];
	  for (size_t kj=ki+1; kj<gaps.size(); ++kj)
	    {
	      int j1 = (gaps[kj].second[0] -  gaps[kj].first[0] >
			gaps[kj].second[1] -  gaps[kj].first[1]) ? 0 : 1;
	      if (i1 == j1)
		continue;
	      double c1 = std::min(gaps[kj].first[j1],gaps[kj].second[j1]);
	      double c2 = std::max(gaps[kj].first[j1],gaps[kj].second[j1]);
	      double d = gaps[kj].first[1-j1];
	      if (d > a1+eps && d < a2-eps && b > c1+eps && b < c2+eps)
		{
		  double mid[2];
		  mid[j1] = b;
		  mid[i1] = d;
		  Point mid_gap(mid[0], mid[1]);
		  gaps.push_back(std::make_pair(mid_gap, gaps[ki].second));
		  gaps.push_back(std::make_pair(mid_gap, gaps[kj].second));
		  gaps[ki].second = mid_gap;
		  gaps[kj].second = mid_gap;
		  joints.push_back(mid_gap);
		}
	    }
	}
    }

  // Split triangulation at joints
#ifdef DEBUG
  std::ofstream oftg("gap_tri.g2");
  std::ofstream oftg0("pre_joint.g2");
#endif
  vector<int> pre_joint_vx;
   for (size_t ki=0; ki<joints.size(); ++ki)
    {
      vector<size_t> tri_near;
      vector<BoundingBox> bb;
      identifyJointTriangles(joints[ki], vertex_par, tri, pdel, tri_near, bb);

      // Check configuration
      int ix1 = -1, ix2 = -1;
      int edge0[2];
      Point pp;
      int config = parConfiguration(joints[ki], vertex_par, tri, pdel, tri_near,
				    ix1, ix2, edge0, pp);

      //size_t nvx = vertex_par.size();
      int ix = (config == 1) ? tri[ix1+ix2] : -1;
      Element2D* elem[4];
      elem[0] = lrsf->coveringElement(joints[ki][0]-pdel, joints[ki][1]-pdel);
      elem[1] = lrsf->coveringElement(joints[ki][0]+pdel, joints[ki][1]-pdel);
      elem[2] = lrsf->coveringElement(joints[ki][0]-pdel, joints[ki][1]+pdel);
      elem[3] = lrsf->coveringElement(joints[ki][0]+pdel, joints[ki][1]+pdel);
      vector<vector<int> > joint_vxs;
      joint_vxs.resize(4);
      for (int kb=0; kb<4; ++kb)
	{
	  if (ix >= 0)
	    {
#ifdef DEBUG
	      writeGapTri(oftg0, vertex, tri[ix1], tri[ix1+1], tri[ix1+2]);
#endif
	      updateVertex(joints[ki], lrsf, dom, vertex_par[ix], vertex[ix], 
			   vertex_norm[ix], vertex_tex[ix], elem[kb]);
	      joint_vxs[kb].push_back(ix);
	      pre_joint_vx.push_back(ix);
	      ix = -1;
#ifdef DEBUG
	      writeGapTri(oftg0, vertex, tri[ix1], tri[ix1+1], tri[ix1+2]);
#endif
	    }
	  else
	    {
	      Point vx, vx_par, vx_norm, vx_tex;
	      updateVertex(joints[ki], lrsf, dom, vx_par, vx, vx_norm, vx_tex, elem[kb]);
	      vertex_par.push_back(vx_par);
	      vertex.push_back(vx);
	      if (mesh_->useNormals())
		vertex_norm.push_back(vx_norm);
	      if (mesh_->useTexCoords())
		vertex_tex.push_back(vx_tex);
	      size_t vx_ix = vertex_par.size() - 1;
	      joint_vxs[kb].push_back((int)vx_ix);
	    }
	}

      vector<int> edge;
      vector<int> ix_tri;
      vector<int> gap_vx;
      if (config == 1)
	{
	  // For all associated triangles, connect to appropriate vertex
	  ix = (int)tri[ix1 + ix2];
	  vector<Point> pvx(2);
	  vector<int> vx_at(2);
	  for (size_t kj=0; kj<tri_near.size(); ++kj)
	    {
	      int kb;
	      for (kb=0; kb<3; ++kb)
		if (tri[tri_near[kj]+kb] == ix)
		  break;

	      if (kb < 3)
		{
		  ix_tri.push_back(tri_near[kj]);
		  edge.push_back(-1);
		  gap_vx.push_back(ix);
		}
	    }
	}
      else
	{
	  if (ix1 >= 0)
	    {
	      ix_tri.push_back(ix1);
	      edge.push_back(edge0[0]);
	      gap_vx.push_back(-1);
	    }
	  if (ix2 >= 0)
	    {
	      ix_tri.push_back(ix2);
	      edge.push_back(edge0[1]);
	      gap_vx.push_back(-1);
	    }
	}

      vector<Point> pvx(2);
      vector<int> vx_at(2);
      vector<Point> par;
      vector<pair<int, int> > quart;
      vector<int> gap_vx_ix;
      for (size_t kj=0; kj<ix_tri.size(); ++kj)
	{
	  int idx = ix_tri[kj];
	  
	  for (int kc=0; kc<3; ++kc)
	    {
	      if (edge[kj] == kc+1  || (int)tri[idx+kc] == gap_vx[kj] ||
				      (int)tri[idx+((kc+1)%3)] == gap_vx[kj])
		continue;

	      size_t pfirst = par.size();
	      pvx[0] = vertex_par[tri[idx+kc]];
	      pvx[1] = vertex_par[tri[idx+((kc+1)%3)]];
	      vector<int> par_at;
	      jointSplitAtEdge(joints[ki], pvx, vx_at, pdel, par, par_at, quart);
	      for (size_t kh=0; kh<par_at.size(); ++kh)
		gap_vx_ix.push_back((par_at[kh] < 0) ? -1 : tri[idx+(kc+par_at[kh])%3]);

	      if (vx_at[0] < 0)
		{
		  // Add existing vertex to joint element pool
		  int kx = (pvx[0][0] < joints[ki][0]) ? 0 : 1;
		  if (pvx[0][1] > joints[ki][1])
		    kx += 2;
		  auto it = std::find(joint_vxs[kx].begin(), joint_vxs[kx].end(),
				      (int)tri[idx+kc]);
		  if (it != joint_vxs[kh].end())
		    joint_vxs[kx].push_back((int)tri[idx+kc]);
		}
	    }
	}

      for (size_t kr=0; kr<gap_vx_ix.size(); ++kr)
	if (gap_vx_ix[kr] >= 0)
	  {
	    size_t kh;
	    for (kh=0; kh<kr; ++kh)
	      if (gap_vx_ix[kh] == gap_vx_ix[kr])
		break;
	    if (kh < kr)
	      gap_vx_ix[kr] = -1;
	  }
      
      for (size_t kr=0; kr<par.size(); ++kr)
	{
	  int kd, kx;
	  for (kd=0, kx=quart[kr].first; kd<2; ++kd, kx=quart[kr].second)
	    {
	      if (kx < 0)
		continue;

	      Element2D *elem2 = elem[kx];
	      if (!elem2->contains(par[kr][0], par[kr][1]))
		{
		  double u = par[kr][0];
		  double v = par[kr][1];
		  if ((kx==0 || kx == 2) && fabs(elem2->umax()-u) < eps)
		    u -= pdel;
		  else if ((kx==1 || kx == 3) && fabs(u-elem2->umin()) < eps)
		    u += pdel;
		  if ((kx==0 || kx == 1) && fabs(elem2->vmax()-v) < eps)
		    v -= pdel;
		  else if ((kx==2 || kx == 3) && fabs(v-elem2->vmin()) < eps)
		    v += pdel;
		  elem2 = lrsf->coveringElement(u, v);
		}
		
	      int ix = gap_vx_ix[kr];
	      if (kd == 0 && quart[kr].second < 0 && ix >= 0 &&
		  (kx == 0 || kx == 3))
		{
		  // Update existing vertex
		  updateVertex(par[kr], lrsf, dom, vertex_par[ix], vertex[ix],
			       vertex_norm[ix], vertex_tex[ix], elem2);
		  joint_vxs[kx].push_back(ix);
		  pre_joint_vx.push_back(ix);
		}
	      else
		{
		  // New vertex at gap line
		  Point vx, vx_par, vx_norm, vx_tex;
		  updateVertex(par[kr], lrsf, dom, vx_par, vx, vx_norm, vx_tex, elem2);
		  vertex_par.push_back(vx_par);
		  vertex.push_back(vx);
		  if (mesh_->useNormals())
		    vertex_norm.push_back(vx_norm);
		  if (mesh_->useTexCoords())
		    vertex_tex.push_back(vx_tex);
		  size_t vx_ix = vertex_par.size() - 1;
		  joint_vxs[kx].push_back((int)vx_ix);
		}
	    }
	}
	

      size_t kr = 0;
      for (size_t kj=0; kj<joint_vxs.size(); ++kj)
	{
	  if (joint_vxs[kj].size() >= 3)
	    {
	      vector<unsigned int*> triang(joint_vxs[kj].size()-2);
	      size_t kh=0;
	      for (; kr<ix_tri.size() && kh<triang.size(); ++kr, ++kh)
		triang[kh] = &tri[ix_tri[kr]];

	      if (kh < triang.size())
		{
		  vector<unsigned int> tmp(3*(triang.size()-kh));
		  tri.insert(tri.end(), tmp.begin(), tmp.end());
		}
	      for (; kh<triang.size(); ++kh)
		triang[kh] = &tri[tri.size()-3*(triang.size()-kh)];
		  
	      sortAndConnect(vertex_par, joint_vxs[kj], triang);

#ifdef DEBUG
	      std::ofstream os("triang_at_joint.g2");
	      os << "400 1 0 4 255 0 0 255" << std::endl;
	      os << joint_vxs[kj].size() << std::endl;
	      for (size_t kr=0; kr<joint_vxs[kj].size(); ++kr)
		os << vertex[joint_vxs[kj][kr]] << std::endl;
	      for (size_t kh=0; kh<triang.size(); ++kh)
		{
		  os << "410 1 0 4 55 200 0 255" << std::endl;
		  os << "3" << std::endl;
		  os << vertex[triang[kh][0]] << " " << vertex[triang[kh][1]] << std::endl;
		  os << vertex[triang[kh][1]] << " " << vertex[triang[kh][2]] << std::endl;
		  os << vertex[triang[kh][2]] << " " << vertex[triang[kh][0]] << std::endl;
		}
	      for (size_t kh=0; kh<triang.size(); ++kh)
		writeGapTri(oftg, vertex, triang[kh][0], triang[kh][1], triang[kh][2]);
#endif
	    }
	}
    }
  
#ifdef DEBUG
  std::ofstream ofm2("joints_mesh.g2");
  ofm2 << "410 1 0 0" << std::endl;
  ofm2 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      ofm2 << vertex[tri[kh]] << "  " << vertex[tri[kh+1]] << " " << std::endl;
      ofm2 << vertex[tri[kh+1]] << "  " << vertex[tri[kh+2]] << " " << std::endl;
      ofm2 << vertex[tri[kh+2]] << "  " << vertex[tri[kh]] << " " << std::endl;
    }
  
  std::ofstream of0("par_tri_joint.g2");
  of0 << "410 1 0 0" << std::endl;
  of0 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      of0 << vertex_par[tri[kh]] << " 0.0 " << vertex_par[tri[kh+1]] << " 0.0" << std::endl;
      of0 << vertex_par[tri[kh+1]] << " 0.0 " << vertex_par[tri[kh+2]] << " 0.0" << std::endl;
      of0 << vertex_par[tri[kh+2]] << " 0.0 " << vertex_par[tri[kh]] << " 0.0" << std::endl;
    }
#endif

   for (size_t ki=0; ki<gaps.size(); ++ki)
    {
      int dir = (gaps[ki].second[0] -  gaps[ki].first[0] >
		 gaps[ki].second[1] -  gaps[ki].first[1]) ? 1 : 0;
      vector<size_t> tri_near;
      vector<BoundingBox> bb;
      identifyNearTriangles(gaps[ki], vertex_par, tri, pdel, tri_near, bb);

 #ifdef DEBUG
     std::cout << "Near triangles identified" << std::endl;
#endif
      // Treat end points
      Point curr;
      int ka;
      int end_ix[2];
      end_ix[0] = end_ix[1] = -1;
      for (ka=0, curr=gaps[ki].first; ka<2; ++ka, curr=gaps[ki].second)
	{
	  size_t kr;
	  for (kr=0; kr<joints.size(); ++kr)
	    if (fabs(curr[0]-joints[kr][0]) < eps && fabs(curr[1]-joints[kr][1]) < eps)
	      break;
	  if (kr < joints.size())
	    continue;
	  if (dir == 1 && (fabs(curr[0]-umin) < eps || fabs(umax-curr[0]) < eps))
	    continue;
	  if (dir == 0 && (fabs(curr[1]-vmin) < eps || fabs(vmax-curr[1]) < eps))
	    continue;
	  
	  vector<size_t> near2;
	  size_t kj;
	  for (kj=0; kj<bb.size(); ++kj)
	    if (bb[kj].containsPoint(curr, pdel))
	      near2.push_back(tri_near[kj]);

#ifdef DEBUG
	  std::cout << "Near2 set, ka = " << ka << std::endl;
#endif

	  int ix1 = -1, ix2 = -1;
	  int edge[2];
	  Point pp;
	  
	  int config = parConfiguration(curr, vertex_par, tri, pdel, near2, ix1, ix2,
					edge, pp);

	  if (config == 1)
	    {
	      // Endpoint close to existing vertex. Update vertex information
	      //unsigned int ix = tri[near2[kj]+kb];
	      unsigned int ix = tri[ix1 + ix2];
	      end_ix[ka] = (int)ix;
	      updateVertex(curr, lrsf, dom, vertex_par[ix], vertex[ix], 
			   vertex_norm[ix], vertex_tex[ix]);
	      std::cout << "Existing vertx updated" << std::endl;
	    }
	  else if (config > 1)
	    {
#ifdef DEBUG
	      std::cout << "Info to create new vertex" << std::endl;
#endif
	      // Make new vertex
	      Point vx, vx_par, vx_norm, vx_tex;
	      updateVertex(curr, lrsf, dom, vx_par, vx, vx_norm, vx_tex);
	      vertex_par.push_back(vx_par);
	      vertex.push_back(vx);
	      if (mesh_->useNormals())
		vertex_norm.push_back(vx_norm);
	      if (mesh_->useTexCoords())
		vertex_tex.push_back(vx_tex);

	      size_t vx_ix = vertex_par.size() - 1;
	      end_ix[ka] = (int)vx_ix;

#ifdef DEBUG
	      std::cout << "New vertex defined" << std::endl;
#endif
	      if (config == 2)
		{
		  // Update triangle(s)
		  for (int kb=0, ix=ix1; kb<2; ++kb, ix=ix2)
		    {
		      if (ix < 0)
			break;
#ifdef DEBUG
		      std::cout << "To update, new vertex at edge, ix = " << ix << std::endl;
#endif
		      tri_near.push_back(tri.size());
		      if (edge[kb] == 1)
			{
			  tri.push_back(tri[ix]);
			  tri.push_back(vx_ix);
			  tri.push_back(tri[ix+2]);
			  tri[ix] = vx_ix;
			}
		      else if (edge[kb] == 2)
			{
			  tri.push_back(tri[ix]);
			  tri.push_back(tri[ix+1]);
			  tri.push_back(vx_ix);
			  tri[ix+1] = vx_ix;
			}
		      else
			{
			  tri.push_back(tri[ix]);
			  tri.push_back(tri[ix+1]);
			  tri.push_back(vx_ix);
			  tri[ix] = vx_ix;
			}
#ifdef DEBUG
		      writeGapTri(oftg, vertex, tri[ix], tri[ix+1], tri[ix+2]);
		      size_t kr=tri.size()-3;
		      writeGapTri(oftg, vertex, tri[kr], tri[kr+1], tri[kr+2]);
#endif
		      int stop_break = 1;
    		    }
		}
	      else if (config == 3) //(ix3 >= 0)
		{
#ifdef DEBUG
		  std::cout << "To update, new vertex inside triangle" << std::endl;
#endif
		  // Split triangle containing the new vertex
		  tri_near.push_back(tri.size());
		  tri.push_back(tri[ix1]);
		  tri.push_back(tri[ix1+1]);
		  tri.push_back(vx_ix);
		  tri_near.push_back(tri.size());
		  tri.push_back(tri[ix1+1]);
		  tri.push_back(tri[ix1+2]);
		  tri.push_back(vx_ix);
		  tri[ix1+1] = vx_ix;
		  
#ifdef DEBUG
		  writeGapTri(oftg, vertex, tri[ix1], tri[ix1+1], tri[ix1+2]);
		  size_t kr=tri.size()-6;
		  writeGapTri(oftg, vertex, tri[kr], tri[kr+1], tri[kr+2]);
		  kr += 3;
		  writeGapTri(oftg, vertex, tri[kr], tri[kr+1], tri[kr+2]);
#endif
		}
#ifdef DEBUG
	      std::cout << "Triangulation updated" << std::endl;
#endif
	    }
	}

 #ifdef DEBUG
     std::ofstream of("par_tri.g2");
      of << "410 1 0 0" << std::endl;
      of << tri.size() << std::endl;
      for (size_t kh=0; kh<tri.size(); kh+=3)
	{
	  of << vertex_par[tri[kh]] << " 0.0 " << vertex_par[tri[kh+1]] << " 0.0" << std::endl;
	  of << vertex_par[tri[kh+1]] << " 0.0 " << vertex_par[tri[kh+2]] << " 0.0" << std::endl;
	  of << vertex_par[tri[kh+2]] << " 0.0 " << vertex_par[tri[kh]] << " 0.0" << std::endl;
	}
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << gaps[ki].first << " 0.0 " << gaps[ki].second << " 0.0" << std::endl;
      
      std::cout << "Endpoints inserted" << std::endl;
#endif

      // Intersect triangulation with gap
      // Find first triangle(s) crossed by the gap

bool do_swap = ((dir == 0 && end_ix[0] >= 0) || (dir == 1 && end_ix[0] < 0));
      if (end_ix[0] <  0)
	std::swap(end_ix[0], end_ix[1]);
      int xt1 = end_ix[0], xt2 = end_ix[0];
      int tri_ix1 = -1, tri_ix2 = -1;
      Point par1 = vertex_par[end_ix[0]];
      Point par2 = (par1.dist(gaps[ki].first) > par1.dist(gaps[ki].second)) ? 
		    gaps[ki].first : gaps[ki].second;
      Point pv(0.0, 0.0);
      pv[dir] = pdel;
      Point pp;
      int vx_ix = -1;
      int next1, next2;
      getNextTriangle(end_ix[0], end_ix[0], par1, par2, tri_near, tri, 
		      vertex_par, pdel, tri_ix1, tri_ix2, pp, vx_ix,
		      next1, next2);
      while (tri_ix1 >= 0 && pp.dimension() == 2)
	{
	  // Compute double vertices
	  Point v1 = Point(pp[1]-par1[1], par1[0]-pp[0]);
	  Point p3 = pp - pv;
	  Point p4 = pp + pv;
	  Element2D *elem1 = lrsf->coveringElement(p3[0], p3[1]);
	  Element2D *elem2 = lrsf->coveringElement(p4[0], p4[1]);
	  int xt3, xt4;

	  bool update = true;
	  if (end_ix[1] >= 0 && end_ix[1] == vx_ix)
	    update = false;
	  else
	    {
	      for (size_t kh=0; kh<pre_joint_vx.size(); ++kh)
		if (vx_ix == pre_joint_vx[kh])
		  update = false;
	    }
      
 	  if (vx_ix >= 0)
	    {
	      if (update)
		updateVertex(pp, lrsf, dom, vertex_par[vx_ix], vertex[vx_ix],
			     vertex_norm[vx_ix], vertex_tex[vx_ix], elem1);
	      xt3 = vx_ix;
	    }
	  else
	    {
	    
	      Point vx1, vx_par1, vx_norm1, vx_tex1;
	      updateVertex(pp, lrsf, dom, vx_par1, vx1, vx_norm1, vx_tex1, elem1);
	      vertex_par.push_back(vx_par1);
	      vertex.push_back(vx1);
	      if (mesh_->useNormals())
		vertex_norm.push_back(vx_norm1);
	      if (mesh_->useTexCoords())
		vertex_tex.push_back(vx_tex1);
	      xt3 = (int)vertex_par.size() - 1;
	    }
	  
	  Point vx2, vx_par2, vx_norm2, vx_tex2;
	  updateVertex(pp, lrsf, dom, vx_par2, vx2, vx_norm2, vx_tex2, elem2);
	  vertex_par.push_back(vx_par2);
	  vertex.push_back(vx2);
	  if (mesh_->useNormals())
	    vertex_norm.push_back(vx_norm2);
	  if (mesh_->useTexCoords())
	    vertex_tex.push_back(vx_tex2);
	  xt4 = (int)vertex_par.size() - 1;

	  // Identify extra triangles
	  vector<int> ntri;
	  if (tri_ix2 > 0 && vx_ix >= 0)
	    {
	      for (size_t kh=0; kh<tri_near.size(); ++kh)
		{
		  for (int kc=0; kc<3; ++kc)
		    if (tri[tri_near[kh]+kc] == vx_ix && (int)tri_near[kh] != tri_ix1 &&
			(int)tri_near[kh] != tri_ix2)
		      {
			ntri.push_back((int)tri_near[kh]);
			break;
		      }
		}
	    }

	  // Update triangulation
	  int tri_ix;
	  int kb;
	  for (kb=0, tri_ix=tri_ix1; kb<2; ++kb, tri_ix=tri_ix2)
	    {
	      if (tri_ix < 0)
		break;
	      vector<int> pre_ix, post_ix, zero_ix, at_ix;
	      vector<Point> ppar(3);
	      for (int kc=0; kc<3; ++kc)
		{
		  if (tri[tri_ix+kc] == xt1 || tri[tri_ix+kc] == xt2 ||
		      tri[tri_ix+kc] == end_ix[1])
		    zero_ix.push_back(kc);
		  else
		    {
		      ppar[kc] = vertex_par[tri[tri_ix+kc]];
		      if (fabs(ppar[kc][dir] - pp[dir]) < eps)
			at_ix.push_back(kc);
		      else if (ppar[kc][dir] < pp[dir])
			pre_ix.push_back(kc);
		      else
			post_ix.push_back(kc);
		    }
		}

	      if (pre_ix.size() >= 1)
		{
		  vector<int> vx_ixs;
		  for (size_t kh=0; kh<pre_ix.size(); ++kh)
		    vx_ixs.push_back(tri[tri_ix+pre_ix[kh]]);
		  vx_ixs.push_back(xt1);
		  vx_ixs.push_back(xt3);
		  vector<unsigned int*> triang(pre_ix.size());
		  size_t first = 0;
		  if (at_ix.size() > 0)
		    triang[first++] = &tri[tri_ix];
		  if (first < triang.size())
		    {
		      vector<unsigned int> tmp(3*(triang.size()-first));
		      tri.insert(tri.end(), tmp.begin(), tmp.end());
		      for (size_t kh=first; kh<triang.size(); ++kh)
			triang[kh] = &tri[tri.size()-3*(triang.size()-kh)];
		    }
		  sortAndConnect(vertex_par, vx_ixs, triang);
#ifdef DEBUG
		  for (size_t kh=0; kh<triang.size(); ++kh)
		    writeGapTri(oftg, vertex, triang[kh][0], triang[kh][1], triang[kh][2]);
#endif
		}

	      if (post_ix.size() >= 1)
		{
		  vector<int> vx_ixs;
		  for (size_t kh=0; kh<post_ix.size(); ++kh)
		    vx_ixs.push_back(tri[tri_ix+post_ix[kh]]);
		  vx_ixs.push_back(xt2);
		  vx_ixs.push_back(xt4);
		  vector<unsigned int*> triang(post_ix.size());
		  size_t first = 0;
		  triang[first++] = &tri[tri_ix];
		  if (first < triang.size())
		    {
		      vector<unsigned int> tmp(3*(triang.size()-first));
		      tri.insert(tri.end(), tmp.begin(), tmp.end());
		      for (size_t kh=first; kh<triang.size(); ++kh)
			triang[kh] = &tri[tri.size()-3*(triang.size()-kh)];
		    }
		  sortAndConnect(vertex_par, vx_ixs, triang);
#ifdef DEBUG
		  for (size_t kh=0; kh<triang.size(); ++kh)
		    writeGapTri(oftg, vertex, triang[kh][0], triang[kh][1], triang[kh][2]);
#endif
		}
	    }

	  for (size_t kh=0; kh<ntri.size(); ++kh)
	    {
	      vector<int> at_ix;
	      int npost = 0;
	      for (int kc=0; kc<3; ++kc)
		{
		  Point ppar = vertex_par[tri[ntri[kh]+kc]];
		  if (fabs(ppar[dir] - pp[dir]) < eps)
		    at_ix.push_back(kc);
		  else if (ppar[dir] > pp[dir])
		    npost++;
		}

	      if (at_ix.size() == 1 && npost == 2)
		tri[ntri[kh]+at_ix[0]] = xt4;
	    }

	  par1 = vertex_par[xt3];
	  int tri_ix3 = -1, tri_ix4 = -1, vx_ix2 = -1;
	  Point pq;
	  int next3, next4;
	  getNextTriangle(next1, next2, par1, par2, tri_near, tri,
			  vertex_par, pdel,
			  tri_ix3, tri_ix4, pq, vx_ix2, next3, next4);
	  
	  xt1 = xt3;
	  xt2 = xt4;
	  tri_ix1 = tri_ix3;
	  tri_ix2 = tri_ix4;
	  if (end_ix[1] >= 0 && end_ix[1] == vx_ix)
	    tri_ix1 = tri_ix2 = -1;
	  else
	    {
	      for (size_t kh=0; kh<pre_joint_vx.size(); ++kh)
		if (vx_ix == pre_joint_vx[kh])
		  tri_ix1 = tri_ix2 = -1;
	    }
	  
	  vx_ix = vx_ix2;
	  pp = pq;
	  next1 = next3;
	  next2 = next4;
	}
    }
      
#ifdef DEBUG
  std::ofstream ofm3("fin_mesh.g2");
  ofm3 << "410 1 0 0" << std::endl;
  ofm3 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      ofm3 << vertex[tri[kh]] << "  " << vertex[tri[kh+1]] << " " << std::endl;
      ofm3 << vertex[tri[kh+1]] << "  " << vertex[tri[kh+2]] << " " << std::endl;
      ofm3 << vertex[tri[kh+2]] << "  " << vertex[tri[kh]] << " " << std::endl;
    }
  
  std::ofstream of2("par_tri2.g2");
  of2 << "410 1 0 0" << std::endl;
  of2 << tri.size() << std::endl;
  for (size_t kh=0; kh<tri.size(); kh+=3)
    {
      of2 << vertex_par[tri[kh]] << " 0.0 " << vertex_par[tri[kh+1]] << " 0.0" << std::endl;
      of2 << vertex_par[tri[kh+1]] << " 0.0 " << vertex_par[tri[kh+2]] << " 0.0" << std::endl;
      of2 << vertex_par[tri[kh+2]] << " 0.0 " << vertex_par[tri[kh]] << " 0.0" << std::endl;
    }

  for (size_t kh=0; kh<gaps.size(); ++kh)
    {
      of2 << "410 1 0 4 255 0 0 255" << std::endl;
      of2 << "1" << std::endl;
      of2 << gaps[kh].first << " 0.0 " << gaps[kh].second << " 0.0" << std::endl;
    }

  std::cout << "Ready to transfer mesh" << std::endl;
 #endif
   
  // Transfer information to mesh_
  mesh_->resize((int)vertex.size(), (int)tri.size()/3);
  for (size_t ki=0; ki<vertex.size(); ++ki)
    {
      int ka;
      for (ka=0; ka<dim; ++ka)
	mesh_->vertexArray()[ki*3+ka] = vertex[ki][ka];
      for (; ka<3; ++ka)
	mesh_->vertexArray()[ki*3+ka] = 0.0;
      mesh_->paramArray()[ki*2] = vertex_par[ki][0];
      mesh_->paramArray()[ki*2+1] = vertex_par[ki][1];
      if (mesh_->useNormals())
	{
	  for (ka=0; ka<3; ++ka)
	    mesh_->normalArray()[ki*3+ka] = vertex_norm[ki][ka];
	}
      if (mesh_->useTexCoords())
	{
	  mesh_->texcoordArray()[ki*2] = vertex_tex[ki][0];
	  mesh_->texcoordArray()[ki*2+1] = vertex_tex[ki][1];
	}
    }

  for (size_t ki=0; ki<tri.size(); ++ki)
    mesh_->triangleIndexArray()[ki] = tri[ki];
}

  
  
//===========================================================================
void LRSurfaceTesselator::identifyJointTriangles(const Point& joint,
						vector<Point>& vx_par,
						vector<unsigned int>& tri,
						double del,
						vector<size_t>& tri_near,
						vector<BoundingBox>& bb_near)
//===========================================================================
{
  BoundingBox joint_bb(2);
  joint_bb.addUnionWith(joint);
  for (size_t ki=0; ki<tri.size(); ki+=3)
    {
      Point pa = vx_par[tri[ki]];
      Point pb = vx_par[tri[ki+1]];
      Point pc = vx_par[tri[ki+2]];
      BoundingBox bb(2);
      bb.addUnionWith(pa);
      bb.addUnionWith(pb);
      bb.addUnionWith(pc);
      if (bb.overlaps(joint_bb, del))
	{
	  tri_near.push_back(ki);
	  bb_near.push_back(bb);
	}
	  
    }
  
}

//===========================================================================
void LRSurfaceTesselator::identifyNearTriangles(pair<Point,Point> gap_par,
						vector<Point>& vx_par,
						vector<unsigned int>& tri,
						double del,
						vector<size_t>& tri_near,
						vector<BoundingBox>& bb_near)
//===========================================================================
{
  Point p1(std::min(gap_par.first[0],gap_par.second[0]),
	   std::min(gap_par.first[1],gap_par.second[1]));
  Point p2(std::max(gap_par.first[0],gap_par.second[0]),
	   std::max(gap_par.first[1],gap_par.second[1]));
  BoundingBox gap_bb(p1, p2);
  for (size_t ki=0; ki<tri.size(); ki+=3)
    {
      Point pa = vx_par[tri[ki]];
      Point pb = vx_par[tri[ki+1]];
      Point pc = vx_par[tri[ki+2]];
      BoundingBox bb(2);
      bb.addUnionWith(pa);
      bb.addUnionWith(pb);
      bb.addUnionWith(pc);
      if (bb.overlaps(gap_bb, del))
	{
	  tri_near.push_back(ki);
	  bb_near.push_back(bb);
	}
	  
    }
  
}

//===========================================================================
void  LRSurfaceTesselator::updateVertex(Point& par, const LRSplineSurface* lrsf,
					RectDomain& dom, Point& vertex_par,
					Point& vertex, Point& vertex_norm,
					Point& vertex_tex, Element2D *elem)
//===========================================================================
{
  lrsf->point(vertex, par[0], par[1], elem);
  if (mesh_->useNormals())
    lrsf->normal(vertex_norm, par[0], par[1], elem);
  if (mesh_->useTexCoords())
    {
      double t1 = (par[0]-dom.umin())*(double)(n_-1)/(dom.umax()-dom.umin());
      double t2 = (par[1]-dom.vmin())*(double)(m_-1)/(dom.vmax()-dom.vmin());
      vertex_par = Point(t1, t2);
    }
  vertex_par = par;
}

//===========================================================================
double  LRSurfaceTesselator::edgeDist(Point curr, Point pa, Point pb, Point& pp)
//===========================================================================
{
  Point vec = pb - pa;
  vec.normalize_checked();
  pp = pa + ((curr - pa)*vec)*vec;
  return curr.dist(pp);
}

//===========================================================================
bool LRSurfaceTesselator::isInsideTri(Point par, Point pa, Point pb, Point pc)
//===========================================================================
{
  Point v1 = pc - pa;
  Point v2 = pb - pa;
  Point vp = par - pa;
  double v11 = v1*v1;
  double v12 = v1*v2;
  double v22 = v2*v2;
  double v1p = v1*vp;
  double v2p = v2*vp;
  double det = v11*v22 - v12*v12;
  double beta = (v1p*v22 - v12*v2p)/det;
  double gamma = (v11*v2p - v1p*v12)/det;
  return (beta > 0.0 && gamma > 0.0 && beta+gamma < 1.0);
}

//===========================================================================
void LRSurfaceTesselator::getNextTriangle(int start_ix1, int start_ix2,
					  Point par1, Point par2,
					  vector<size_t>& tri_near,
					  vector<unsigned int>& tri,
					  vector<Point>& vx_par, double pdel,
					  int& tri_ix1, int& tri_ix2, Point& pp,
					  int& vx_ix, int& next_ix1,
					  int& next_ix2)
//===========================================================================
 {
   double eps = 1.0e-9;
   vx_ix = -1;
   
   for (size_t kj=0; kj<tri_near.size(); ++kj)
     {
       bool found1 = false, found2 = false;
       if ((int)tri[tri_near[kj]] == start_ix1 || (int)tri[tri_near[kj]+1] == start_ix1 ||
	   (int)tri[tri_near[kj]+2] == start_ix1)
	 found1 = true;
       if ((int)tri[tri_near[kj]] == start_ix2 || (int)tri[tri_near[kj]+1] == start_ix2 ||
	   (int)tri[tri_near[kj]+2] == start_ix2)
	 found2 = true;

       if (!(found1 && found2))
	 continue;

       vector<unsigned int> idvx;
       for (size_t kr=0; kr<3; ++kr)
	 if (((int)tri[tri_near[kj]+kr] != start_ix1 &&
	      (int)tri[tri_near[kj]+kr] != start_ix2) || start_ix1 != start_ix2)
	   idvx.push_back(tri[tri_near[kj]+kr]);
       vector<Point> pvx(idvx.size());
       for (size_t kr=0; kr<idvx.size(); ++kr)
	 pvx[kr] = vx_par[idvx[kr]];

       double dd = std::numeric_limits<double>::max();
       int tmp_ix = -1;
       Point pp1;
       for (size_t kr=0; kr<idvx.size(); ++kr)
	 {
	   Point pp0;
	   double dd0 = edgeDist(pvx[kr], par1, par2, pp0);
	   if (pp0.dimension() != par1.dimension())
	     continue;
	   if (dd0 < dd && (par2-par1)*(pp0-par1) > 0.0)
	     {
	       dd = dd0;
	       pp1 = pp0;
	       tmp_ix = (int)kr;
	     }
	 }

       if (dd < pdel)
	 {
	   pp = pp1;
	   vx_ix = idvx[tmp_ix];
	   next_ix1 = next_ix2 = (int)idvx[tmp_ix];
	   if (tri_ix1 < 0)
	     tri_ix1 = (int)tri_near[kj];
	   else
	     tri_ix2 = (int)tri_near[kj];
	 }
       else
	 {
	   // Next gap point is not at an existing vertex
	   Point v1 = par2 - par1;
	   for (size_t kr=0; kr<pvx.size(); ++kr)
	     {
	       size_t kh = (kr+1)%pvx.size();
	       Point v2 = pvx[kh] - pvx[kr];
	       Point v3 = pvx[kr] - par1;
	       double det = v1[0]*v2[1] - v1[1]*v2[0];
	       if (fabs(det) < eps)
		 {
		   if (fabs(v1[0]*v3[1] - v1[1]*v3[0]) < eps)
		     {
		       if (tri_ix1 < 0)
			 tri_ix1 = (int)tri_near[kj];
		       else
			 tri_ix2 = (int)tri_near[kj];
		       next_ix1 = next_ix2 = (int)idvx[kr];
		       std::swap(tri_near[kj], tri_near[tri_near.size()-1]);
		       tri_near.pop_back();
		     }
		 }
	       else
		 {
		   double t1 = (v3[0]*v2[1]-v3[1]*v2[0])/det;
		   double t2 = (v3[0]*v1[1]-v3[1]*v1[0])/det;
		   double ppdist = par1.dist( pvx[kr] + t2*v2);
		   if (t1 >= 0.0 && t1 <= 1.0 && t2 >= 0.0 && t2 <= 1.0 && ppdist > pdel)
		     {
		       if (tri_ix1 < 0)
			 tri_ix1 = (int)tri_near[kj];
		       else
			 tri_ix2 = (int)tri_near[kj];
		       std::swap(tri_near[kj], tri_near[tri_near.size()-1]);
		       tri_near.pop_back();
		       
		       pp = pvx[kr] + t2*v2;
		       next_ix1 = (int)idvx[kr];
		       next_ix2 = (int)idvx[kh];
		       break;
		     }
		 }
	     }
	 }
       if (tri_ix1 >=0 && tri_ix2 >= 0)
	 break;
     }
 }


 
//===========================================================================
int LRSurfaceTesselator::parConfiguration(const Point& param,
					  vector<Point>& vx_par,
					  vector<unsigned int>& tri,
					  double del,
					  vector<size_t>& tri_near,
					  int& ix1, int& ix2, 
					  int edge[], Point& pp)
//===========================================================================
 {
   // Check distance to vertex
   size_t kj;
   for (kj=0; kj<tri_near.size(); ++kj)
     {
       int kb;
       for (kb=0; kb<3; ++kb)
	 {
	   double dd = param.dist(vx_par[tri[tri_near[kj]+kb]]);
	   if (dd < del)
	     break;
	 }
       if (kb < 3)
	 {
	   ix1 = (int)tri_near[kj];
	   ix2 = kb;
	   break;
	 }
     }
   
   if (kj < tri_near.size())
     return 1;

   // Look for a close edge
   ix1 = ix2 = -1;
   int ix3 = -1;
   double dist[2];
   for (size_t kr=0; kr<tri_near.size(); ++kr)
     {
       Point pa = vx_par[tri[tri_near[kr]]];
       Point pb = vx_par[tri[tri_near[kr]+1]];
       Point pc = vx_par[tri[tri_near[kr]+2]];
       Point pp1, pp2, pp3;
       double dd1 = edgeDist(param, pa, pb, pp1);
       double dd2 = edgeDist(param, pa, pc, pp2);
       double dd3 = edgeDist(param, pb, pc, pp3);
       if (std::min(dd1, std::min(dd2, dd3)) < del)
	 {
	   if (ix1 >= 0 && ix2 >= 0)
	     {
	       double dmin = std::min(dd1, std::min(dd2, dd3));
	       if (dmin < dist[0] && dist[0] > dist[1])
		 ix1 = -1;
	       else if (dmin < dist[1])
		 ix2 = -1;
	     }
	   if (ix1 < 0 || ix2 < 0)
	     pp = (dd1 < std::min(dd2, dd3)) ? pp1 :
	       ((dd2 < dd3) ? pp2 : pp3);
	   if (ix1 < 0)
	     {
	       ix1 = (int)tri_near[kr];
	       edge[0] = (dd1 < std::min(dd2, dd3)) ? 1 :
		 ((dd2 < dd3) ? 3 : 2);
	       dist[0] = std::min(dd1, std::min(dd2, dd3));
	     }
	   else if (ix2 < 0)
	     {
	       ix2 = (int)tri_near[kr];
	       edge[1] =  (dd1 < std::min(dd2, dd3)) ? 1 :
		 ((dd2 < dd3) ? 3 : 2);
	       dist[1] = std::min(dd1, std::min(dd2, dd3));
	     }
	 }

       else if (isInsideTri(param, pa, pb, pc))
	 ix3  = (int)tri_near[kr];
     }

   if (ix1 > 0)
     return 2;
   else if (ix3 > 0)
     {
       ix1 = ix3;
       return 3;
     }

   return 0;
 }

//===========================================================================
void LRSurfaceTesselator::sortAndConnect(vector<Point>& vx_par, vector<int>& vx_ix,
					 vector<unsigned int*> triangle)
//===========================================================================
{
  if (vx_ix.size() < 3 || vx_ix.size() > 4)
    return;   // Wrong use
  if ((vx_ix.size() == 3 && triangle.size() != 1) ||
      (vx_ix.size() == 4 && triangle.size() != 2))
    return;

  Point mid(0.0, 0.0);
  double fac = 1.0/(double)vx_ix.size();
  for (size_t ki=0; ki<vx_ix.size(); ++ki)
    mid += fac*vx_par[vx_ix[ki]];

  vector<double> angle(vx_ix.size());
  for (size_t ki=0; ki<vx_ix.size(); ++ki)
    angle[ki] = atan2(vx_par[vx_ix[ki]][1]-mid[1], vx_par[vx_ix[ki]][0]-mid[0]);

  for (size_t ki=0; ki<vx_ix.size(); ++ki)
    for (size_t kj=ki+1; kj<vx_ix.size(); ++kj)
      if (angle[kj] < angle[ki])
	{
	  std::swap(vx_ix[ki], vx_ix[kj]);
	  std::swap(angle[ki], angle[kj]);
	}

  if (vx_ix.size() == 3)
    {
      for (size_t ki=0; ki<vx_ix.size(); ++ki)
	triangle[0][ki] = vx_ix[ki];
    }
  else // if (vx_ix.size() == 4
    {
      double d1 = vx_par[vx_ix[0]].dist(vx_par[vx_ix[2]]);
      double d2 = vx_par[vx_ix[1]].dist(vx_par[vx_ix[3]]);
      if (d1 < d2)
	{
	  triangle[0][0] = vx_ix[0];
	  triangle[0][1] = vx_ix[1];
	  triangle[0][2] = vx_ix[2];
	  triangle[1][0] = vx_ix[0];
	  triangle[1][1] = vx_ix[2];
	  triangle[1][2] = vx_ix[3];
	}
      else
	{
	  triangle[0][0] = vx_ix[0];
	  triangle[0][1] = vx_ix[1];
	  triangle[0][2] = vx_ix[3];
	  triangle[1][0] = vx_ix[1];
	  triangle[1][1] = vx_ix[2];
	  triangle[1][2] = vx_ix[3];
	}
    }
  
  
}

//===========================================================================
void LRSurfaceTesselator::jointSplitAtEdge(const Point& joint,
					   vector<Point>& pvx,
					   vector<int>& at,
					   double del, vector<Point>& par,
					   vector<int>& par_at,
					   vector<pair<int, int> >& parquart)
//===========================================================================
{
  if (pvx.size() != 2 || at.size() != 2)
    return;

  int quart[2];
  at[0] = at[1] = -1;
  quart[0] = quart[1] = 0;
  for (int kb=0; kb<2; kb++)
    {
      for (int kc=0; kc<2; ++kc)
	{
	  if (fabs(pvx[kb][kc]-joint[kc]) < del)
	    at[kb] = kc;
	  if (pvx[kb][kc] > joint[kc])
	    quart[kb] += ((kc == 1) ? 2 : 1);
	}
    }

  if (at[0] >= 0 && at[1] == at[0]);
  else
    {
      double low[2], high[2];
      for (int kc=0; kc<2; ++kc)
	{
	  low[kc]= std::min(pvx[0][kc],pvx[1][kc]);
	  high[kc] = std::max(pvx[0][kc],pvx[1][kc]);
	}
      for (int kc=0; kc<2; ++kc)
	{
	  if (at[kc] >= 0)
	    {
	      int ix = at[kc];
	      int jx = 1 - ix;

	      Point ppar(2);
	      ppar[ix] = joint[ix];
	      ppar[1-ix] = (at[kc] >= 0) ? pvx[kc][1-ix] : pvx[1-kc][1-ix];
	      par.push_back(ppar);
	      double pval1 = (at[0] >= 0) ? std::min(pvx[0][0], pvx[1][0]) :
		ppar[0];
	      int qval = (low[kc] <= joint[kc]-del) ? 0 : ((kc == 0) ? 1 : 2);
	      double pval2 = (fabs(ppar[1-kc]-joint[1-kc]) < del) ? high[1-kc] : ppar[1-kc];
	      if (pval2 >= joint[1-kc]+del)
		qval += (kc == 0) ? 2 : 1;

	      parquart.push_back(std::make_pair(qval, -1));
	      par_at.push_back(kc);
	    }
	  else
	    {
	      if (low[kc] <= joint[kc]-del && high[kc] >= joint[kc]+del)
		{
		  Point ppar(2);
		  ppar[kc] = joint[kc];
		  double z = (joint[kc] - low[kc])/(high[kc] - joint[kc]);
		  ppar[1-kc] = low[1-kc] + z*(high[1-kc]-low[1-kc])/(1.0+z);
		  par.push_back(ppar);
		  int qval1 = (pvx[0][kc] < pvx[1][kc]) ? 0 : ((kc == 0) ? 1 : 2);
		  int qval2 = (qval1 > 0) ? 0 : ((kc == 0) ? 1 : 2);
		  if (at[1-kc] >= 0 && 0.5*(low[1-kc] + high[1-kc]) > joint[1-kc])
		    {
		      qval1 += (kc == 0) ? 2 : 1;
		      qval2 += (kc == 0) ? 2 : 1;
		    }
		  else if (at[1-kc] < 0)
		    {
		      if (pvx[0][1-kc] > joint[1-kc])
			qval1 += (kc == 0) ? 2 : 1;
		      if (pvx[1][1-kc] > joint[1-kc])
			qval2 += (kc == 0) ? 2 : 1;
		    }		    
		  parquart.push_back(std::make_pair(qval1, qval2));
		  par_at.push_back(-1);
		}
	    }
	}
				    
    }
}

//===========================================================================
void LRSurfaceTesselator::writeGapTri(std::ofstream& oftg, vector<Point>& vertex,
				      unsigned int ix1, unsigned int ix2,
				      unsigned int ix3)
//===========================================================================
{
  oftg << "410 1 0 4 155 50 50 255" << std::endl;
  oftg << "3" << std::endl;
  oftg << vertex[ix1] << " " << vertex[ix2] << std::endl;
  oftg << vertex[ix2] << " " << vertex[ix3] << std::endl;	
  oftg << vertex[ix3] << " " << vertex[ix1] << std::endl;	
}
