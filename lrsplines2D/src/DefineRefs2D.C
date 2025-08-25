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

#include "GoTools/lrsplines2D/DefineRefs2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"

using namespace Go;
using std::vector;

//==============================================================================
// Prepare input for structured mesh refinement corresponding to a selection
// of B-splines
//
void DefineRefs2D::refineStructuredMesh(const LRSplineSurface& surf,
					vector<LRBSpline2D*>& source,
					vector<LRSplineSurface::Refinement2D>& refs_x,
					vector<LRSplineSurface::Refinement2D>& refs_y,
					bool adjust, bool reduced)
//==============================================================================
{
  double tol = surf.getKnotTol();
  if (source.size() == 0)
    return;
  const Mesh2D& mesh = surf.mesh();
  for (size_t ki=0; ki<source.size(); ++ki)
    {
      // Refine all knot spans in both parameter directions
      int size1 = source[ki]->degree(XFIXED)+1;
      int size2 = source[ki]->degree(YFIXED)+1;
      const vector<int>& kvec1 = source[ki]->kvec(XFIXED);
      const vector<int>& kvec2 = source[ki]->kvec(YFIXED);
      double umin = source[ki]->umin();
      double umax = source[ki]->umax();
      int uxmin = source[ki]->suppMin(XFIXED);
      int uxmax = source[ki]->suppMax(XFIXED);
      double vmin = source[ki]->vmin();
      double vmax = source[ki]->vmax();
      int vxmin = source[ki]->suppMin(YFIXED);
      int vxmax = source[ki]->suppMax(YFIXED);
      for (int kj=0; kj<size1; ++kj)
	{
	  double u1 = mesh.kval(XFIXED, kvec1[kj]);
	  double u2 = mesh.kval(XFIXED, kvec1[kj+1]);
	  if (u2-u1 <= tol)
	    continue;
	  double upar;
	  if (kvec1[kj+1]-kvec1[kj] > 1 && adjust)
	    {
	      int ix = kj;
	      double maxdist = 0.0;
	      double mid = (mesh.kval(XFIXED, kvec1[kj+1]) +
			    mesh.kval(XFIXED, kvec1[kj]))/2;
	      for (int ka=kvec1[kj]+1; ka<kvec1[kj+1]; ++ka)
		{
		  const vector<GPos> mr = mesh.mrects(XFIXED, ka);
		  int nu1, nu2;
		  for (nu1=0; nu1<mr.size()-1; ++nu1)
		    if (mr[nu1+1].ix > vxmin)
		      break;
		  for (nu2=0; nu2<mr.size(); ++nu2)
		    if (mr[nu2].ix > vxmax)
		      break;
		  double dist = 0.0;
		  for (int kb=nu1; kb<nu2; ++kb)
		    if (mr[kb].mult > 0)
		      {
			int l1 = std::max(vxmin, mr[kb].ix);
			int l2 = (kb < mr.size()-1) ? std::min(vxmax, mr[kb+1].ix) : vxmax;
			dist += fabs(mesh.kval(YFIXED, l2) - mesh.kval(YFIXED, l1));
		      }
		  if (fabs(mesh.kval(XFIXED,ka)-mid) < fabs(mesh.kval(XFIXED,ix)-mid) &&
		      dist > 0.0)
		    {
		      ix = ka;
		      if (dist > maxdist)
			maxdist = dist;
		    }
		  else if (dist > 2.0*maxdist)
		    {
		      ix = ka;
		      maxdist = dist;
		    }
		}
	      if (maxdist > 0.0)
		upar =  mesh.kval(XFIXED, ix);
	      else if (reduced)
		continue;
	      else
		upar = 0.5*(u1+u2);
	    }
	  else if (reduced)
	    continue;
	  else
	    upar = 0.5*(u1+u2);
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(upar, source[ki]->vmin(),
			  source[ki]->vmax(), XFIXED, 1);
	  appendRef(refs_x, curr_ref, tol);
	}
      
      for (int kj=0; kj<size2; ++kj)
	{
	  double v1 = mesh.kval(YFIXED, kvec2[kj]);
	  double v2 = mesh.kval(YFIXED, kvec2[kj+1]);
	  if (v2-v1 <= tol)
	    continue;
	  double vpar;
	  if (kvec2[kj+1]-kvec2[kj] > 1 && adjust)
	    {
	      int ix = kj;
	      int maxdist = 0.0;
	      double mid = (mesh.kval(YFIXED, kvec2[kj+1]) +
			    mesh.kval(YFIXED, kvec2[kj]))/2;
	      for (int ka=kvec2[kj]+1; ka<kvec2[kj+1]; ++ka)
		{
		  const vector<GPos> mr = mesh.mrects(YFIXED, ka);
		  int nu1, nu2;
		  for (nu1=0; nu1<mr.size()-1; ++nu1)
		    if (mr[nu1+1].ix > uxmin)
		      break;
		  for (nu2=0; nu2<mr.size(); ++nu2)
		    if (mr[nu2].ix > uxmax)
		      break;
		  double dist = 0.0;
		  for (int kb=nu1; kb<nu2; ++kb)
		    if (mr[kb].mult > 0)
		      {
			int l1 = std::max(uxmin, mr[kb].ix);
			int l2 = (kb < mr.size()-1) ? std::min(uxmax, mr[kb+1].ix) : uxmax;
			dist += fabs(mesh.kval(XFIXED, l2) - mesh.kval(XFIXED, l1));
		      }
		  if (fabs(mesh.kval(YFIXED,ka)-mid) < fabs(mesh.kval(YFIXED,ix)-mid) &&
		      dist > 0.0)
		    {
		      ix = ka;
		      if (dist > maxdist)
			maxdist = dist;
		    }
		  else if (dist > 2.0*maxdist)
		    {
		      ix = ka;
		      maxdist = dist;
		    }
		}
	      if (maxdist > 0.0)
		vpar =  mesh.kval(YFIXED, ix);
	      else if (reduced)
		continue;
	      else
		vpar = 0.5*(v1+v2);
	    }
	  else if (reduced)
	    continue;
	  else
	    vpar = 0.5*(v1+v2);
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(vpar, source[ki]->umin(),
			  source[ki]->umax(), YFIXED, 1);
	  appendRef(refs_y, curr_ref, tol);
	}
    }

}

//==============================================================================
void DefineRefs2D::appendRef(vector<LRSplineSurface::Refinement2D>& refs,
			     LRSplineSurface::Refinement2D& curr_ref, double tol)
//==============================================================================
{
  // Check if the current refinement can be combined with an existing one
  size_t ki;
  for (ki=0; ki<refs.size(); ++ki)
    {
      // Check direction and knot value
      if (/*refs[ki].d == curr_ref.d &&*/ 
	  fabs(refs[ki].kval-curr_ref.kval) < tol)
	{
	  // Check extent of refinement
	  if (!(refs[ki].start > curr_ref.end+tol ||
		curr_ref.start > refs[ki].end+tol))
	    {
	      // Merge new knots
	      refs[ki].start = std::min(refs[ki].start, curr_ref.start);
	      refs[ki].end = std::max(refs[ki].end, curr_ref.end);
	      break;
	    }
	}
    }

  if (ki == refs.size())
    refs.push_back(curr_ref);
}

