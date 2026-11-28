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
#include "GoTools/lrsplines2D/Mesh2DUtils.h"

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
void DefineRefs2D::defineRefsAtEdge(const Mesh2D& mesh, Direction2D dir,
				    int edge, int ix, const vector<double>& knot_vals, 
				    int element_width,
				    vector<LRSplineSurface::Refinement2D>& refs)
//==============================================================================
{
  Direction2D curr_dir = flip(dir);
  for (size_t kj=0; kj<knot_vals.size(); ++kj)
    {
      double curr_knot = knot_vals[kj];
      int other_ix = 
	Mesh2DUtils::last_nonlarger_knotvalue_ix(mesh, dir, curr_knot);

      double p1 = mesh.kval(curr_dir, ix);
      double p2;
      if (edge == 1 || edge == 3)
      {
	  int c_ix = ix;
	  for (int ki=0; ki<element_width; ++ki)
	  {
	      int p_ix = c_ix;
	      c_ix = // We search for the next line which contains curr_knot.
		  Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh, curr_dir,
									 c_ix-1, other_ix);
	      // // If segment already exists we must decrease p1.
	      // int mult = mesh.nu(dir, other_ix, c_ix, p_ix);
	      // if (mult > 0)
	      // {
	      // 	  MESSAGE("Do something!");
	      // 	  p1 = mesh.kval(curr_dir, c_ix);
	      // }
	  }
	  p2 = mesh.kval(curr_dir, c_ix);
	}
      else
	{
	  int c_ix = ix;
	  for (int ki=0; ki<element_width; ++ki)
	  {
	      int p_ix = c_ix;
	      c_ix =
		  Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh, curr_dir,
								       c_ix+1, other_ix);
	      // // If segment already exists we must increase p1.
	      // int mult = mesh.nu(dir, other_ix, p_ix, c_ix);
	      // if (mult > 0)
	      // {
	      // 	  MESSAGE("Do something!");
	      // 	  p1 = mesh.kval(curr_dir, c_ix);
	      // }
	  }
	  p2 = mesh.kval(curr_dir, c_ix);
	}

      if (p1 > p2)
	std::swap(p1, p2);

      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(curr_knot, p1, p2, dir, 1);
      refs.push_back(curr_ref);
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

//==============================================================================
//  Structured mesh refinement given bspline
//
void DefineRefs2D::refineStructuredMesh(const LRSplineSurface& surf,
					LRBSpline2D* bspline, Direction2D fixdir,
					int mult,
					vector<LRSplineSurface::Refinement2D>& refs)
//==============================================================================
{
  double tol = surf.getKnotTol();
  Direction2D d2 = flip(fixdir);
  int size = bspline->degree(fixdir)+1;
  const Mesh2D* mesh = bspline->getMesh();
  double pmin = bspline->pmin(d2);
  double pmax = bspline->pmax(d2);
  const vector<int>& kvec = bspline->kvec(fixdir);
  for (int kj=0; kj<size; ++kj)
    {
      double p1 = mesh->kval(fixdir, kvec[kj]);
      double p2 = mesh->kval(fixdir, kvec[kj+1]);
      if (p2-p1 <= tol)
	continue;
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(0.5*(p1+p2), pmin, pmax, fixdir, mult);
      appendRef(refs, curr_ref, tol);
    }
}

//==============================================================================
//  Full span refinement given parameter at which to refine
//
void DefineRefs2D::refineFullSpan(const LRSplineSurface& surf,
				  double upar, double vpar, Direction2D fixdir,
				  int mult,
				  LRSplineSurface::Refinement2D& refs)
//==============================================================================
{
  Direction2D d2 = flip(fixdir);

  Element2D *elem = surf.coveringElement(upar, vpar);
  if (!elem)
    return;
  double pmin = elem->pmin(d2);
  double pmax = elem->pmax(d2);

  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();
  
  // All overlapping B-splines
  for (size_t ki=0; ki<nmb; ++ki)
    {
      double bmin = bsplines[ki]->pmin(d2);
      double bmax = bsplines[ki]->pmax(d2);
      pmin = std::min(pmin, bmin);
      pmax = std::max(pmax, bmax);
    }

  refs.setVal((fixdir == XFIXED) ? upar : vpar, pmin, pmax, fixdir, mult);
}
 
//==============================================================================
//  Element based refinement 
//  strategies: fullSpan, minSpanSize (minimum span based on B-spline support),
//  minSpanEror (minimum span based on error in associated data points),
//  minSpanCombined (minimum span based on combination fof size and error)
void DefineRefs2D::refineFromElement(const LRSplineSurface& surf,
				     Element2D *elem, Direction2D fixdir,
				     RefineStrategy strategy, int mult,
				     LRSplineSurface::Refinement2D& refs)
//==============================================================================
{
  Direction2D d2 = flip(fixdir);
  
  // Fetch B-splines
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  double ppar = (fixdir == XFIXED) ? 0.5*(elem->umin() + elem->umax()) :
    0.5*(elem->vmin() + elem->vmax());
  double pmin = elem->pmin(d2);
  double pmax = elem->pmax(d2);
  double par2 = 0.5*(pmin + pmax); 
  if (strategy == fullSpan)
    {
      // All overlapping B-splines
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  double bmin = bsplines[ki]->pmin(d2);
	  double bmax = bsplines[ki]->pmax(d2);
	  pmin = std::min(pmin, bmin);
	  pmax = std::max(pmax, bmax);
	}
    }
  else if (strategy == minSpanSize)
    {
      // "Largest" overlapping B-spline. Weighted with mid split
      double tol = surf.getKnotTol();
      double max_size = 0.0;
      int ix = -1;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Compute size of B-spline
	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin = bsplines[ki]->pmin(d2);
	  double bmax = bsplines[ki]->pmax(d2);

	  double bsize = 2.0*(bmax0 - bmin0) + bmax - bmin;
	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin, bmax-par2)/std::max(par2-bmin, bmax-par2);
	  double frac = 2.0*frac0 + frac1;

	  if (frac*bsize > max_size)
	    {
	      max_size = frac*bsize;
	      ix = (int)ki;
	    }
	}
      pmin = bsplines[ix]->pmin(d2);
      pmax = bsplines[ix]->pmax(d2);
    }
  else if (strategy == minSpanError)
    {
      // "Best" overlapping B-spline
      double tol = 0.1;
      double max_wgt = 0.0;
      double max_frac = 0.0;
      int ix = -1;
     for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Count the number of elements with large error affected
	  double curr_wgt = 0.0;
	  const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = curr_el[kj]->pmin(fixdir);
	      double emax = curr_el[kj]->pmax(fixdir);
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside, nmb_out_sign;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
					   nmb_out_sign);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin = bsplines[ki]->pmin(d2);
	  double bmax = bsplines[ki]->pmax(d2);

	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin, bmax-par2)/std::max(par2-bmin, bmax-par2);
	  double frac = 2.0*frac0 + frac1;

	  if ((fabs(max_wgt-curr_wgt) < tol && frac > max_frac) ||
	      curr_wgt > max_wgt)
	    {
	      max_wgt = curr_wgt;
	      max_frac = frac;
	      ix = (int)ki;
	    }
	}
      pmin = bsplines[ix]->pmin(d2);
      pmax = bsplines[ix]->pmax(d2);
    }
  else  if (strategy == minSpanCombined)
    {
      // Combination of size and error
      double tol = 0.1;
      double max_wgt = 0.0;
      double max_size = 0.0;
      double max_frac = 0.0;
      double maxfrac_combined = 0.0;
      int ix = -1;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Count the number of elements with large error affected
	  double curr_wgt = 0.0;
	  const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = curr_el[kj]->pmin(fixdir);
	      double emax = curr_el[kj]->pmax(fixdir);
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside, nmb_out_sign;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
					   nmb_out_sign);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  // Compute size of B-spline
	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin = bsplines[ki]->pmin(d2);
	  double bmax = bsplines[ki]->pmax(d2);

	  double bsize = 2.0*(bmax0 - bmin0) + bmax - bmin;
	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin, bmax-par2)/std::max(par2-bmin, bmax-par2);
	  double frac = 2.0*frac0 + frac1;

	  double frac_combined = (ix < 0) ? 1.0 : curr_wgt/max_wgt + bsize/max_size;
	  if ((fabs(maxfrac_combined-frac_combined) < tol && frac > max_frac) ||
	      frac_combined > maxfrac_combined)
	    {
	      maxfrac_combined = frac_combined;
	      max_wgt = curr_wgt;
	      max_size = bsize;
	      max_frac = frac;
	      ix = (int)ki;
	    }
	}
      
      pmin = bsplines[ix]->pmin(d2);
      pmax = bsplines[ix]->pmax(d2);
   }
  
  refs.setVal(ppar, pmin, pmax, fixdir, mult);
}
