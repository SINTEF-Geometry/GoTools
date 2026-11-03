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

#include "GoTools/lrsplines3D/DefineRefs3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/Mesh3DUtils.h"

using namespace Go;
using std::vector;

//==============================================================================
// Prepare input for structured mesh refinement corresponding to a selection
// of B-splines
//
void DefineRefs3D::refineStructuredMesh(const LRSplineVolume& vol,
					vector<LRBSpline3D*>& source,
					vector<LRSplineVolume::Refinement3D>& refs_x,
					vector<LRSplineVolume::Refinement3D>& refs_y,
					vector<LRSplineVolume::Refinement3D>& refs_z,
					bool adjust, bool reduced)
//==============================================================================
{
  double tol = vol.getKnotTol();
  if (source.size() == 0)
    return;
  const Mesh3D& mesh = vol.mesh();
  Direction3D dir[3] = {XDIR, YDIR, ZDIR};
  for (size_t ki=0; ki<source.size(); ++ki)
    {
      // Refine all knot spans in all parameter directions
      for (int ka=0; ka<3; ++ka)
	{
	  int ka2 = (ka+1)%3;
	  int ka3 = (ka>0) ? ka-1 : 2;
	  int size = source[ki]->degree(dir[ka]) + 1;
	  const vector<int>& kvec = source[ki]->kvec(dir[ka]);
	  double pmin = source[ki]->pmin(dir[ka]);
	  double pmax = source[ki]->pmax(dir[ka]);
	  int ixmin1 = source[ki]->suppMin(dir[ka2]);
	  int ixmax1 = source[ki]->suppMax(dir[ka2]);
	  int ixmin2 = source[ki]->suppMin(dir[ka3]);
	  int ixmax2 = source[ki]->suppMax(dir[ka3]);

	  for (int kj=0; kj<size; ++kj)
	    {
	      double p1 = mesh.kval(dir[ka], kvec[kj]);
	      double p2 = mesh.kval(dir[ka], kvec[kj+1]);
	      if (p2-p1 <= tol)
		continue;
	      double ppar;
	      if (kvec[kj+1]-kvec[kj] > 1 && adjust)
		{
		  int ix = kj;
		  double maxdist = 0.0;
		  double mid = (mesh.kval(dir[ka], kvec[kj+1]) +
				mesh.kval(dir[ka], kvec[kj]))/2;
		  for (int kb=kvec[kj]+1; kb<kvec[kj+1]; ++kb)
		    {
		      double dist = 0.0;
		      const vector<GPos2D> mr = mesh.mrects(dir[ka], kb);
		      size_t kr;
		      for (kr=0; kr<mr.size(); ++kr)
			{
			  if (mr[kr].mult == 0)
			    continue;
			  
			  if (mr[kr].ll[0] < ixmax1 && mr[kr].ur[0] > ixmin1 &&
			      mr[kr].ll[1] < ixmax2 && mr[kr].ur[1] > ixmin2)
			    {
			      int l1 = std::max(ixmin1, mr[kr].ll[0]);
			      int l2 = std::min(ixmax1, mr[kr].ur[0]);
			      int l3 = std::max(ixmin2, mr[kr].ll[1]);
			      int l4 = std::min(ixmax2, mr[kr].ur[1]);

			      dist += fabs(mesh.kval(dir[ka2],l2) - mesh.kval(dir[ka2],l1));
			      dist += fabs(mesh.kval(dir[ka3],l4) - mesh.kval(dir[ka3],l3));
			    }
			}

		      if (fabs(mesh.kval(dir[ka],kb)-mid) < fabs(mesh.kval(dir[ka],ix)-mid) &&
			  dist < 0.0)
			{
			  ix = kb;
			  if (dist > maxdist)
			    maxdist = dist;
			}
		      else if (dist > 2.0*maxdist)
			{
			  ix = kb;
			   maxdist = dist;
			}
		    }
		  if (maxdist > 0.0)
		    ppar =  mesh.kval(dir[ka], ix);
		  else if (reduced)
		    continue;
		  else
		    ppar = 0.5*(p1+p2);
		}
	      else if (reduced)
		continue;
	      else
		ppar = 0.5*(p1+p2);
	      LRSplineVolume::Refinement3D curr_ref;
	      curr_ref.setVal(ppar, source[ki]->pmin(dir[ka2]),
			      source[ki]->pmax(dir[ka2]),
			      source[ki]->pmin(dir[ka3]),
			      source[ki]->pmax(dir[ka3]), dir[ka], 1);
	      if (ka == 0)
		appendRef(refs_x, curr_ref, tol);
	      else if (ka == 1)
		appendRef(refs_y, curr_ref, tol);
	      else
		appendRef(refs_z, curr_ref, tol);
	    }
	}
      
    }
}


//==============================================================================
void DefineRefs3D::appendRef(vector<LRSplineVolume::Refinement3D>& refs,
			     LRSplineVolume::Refinement3D& curr_ref, double tol)
//==============================================================================
{
  // Check if the current refinement can be combined with an existing one
  size_t ki;
  for (ki=0; ki<refs.size(); ++ki)
    {
      // Check direction and knot value
      if (fabs(refs[ki].kval-curr_ref.kval) < tol)
	{
	  // Check extent of refinement
	  if (fabs(refs[ki].start1-curr_ref.start1) < tol &&
	      fabs(refs[ki].end1-curr_ref.end1) < tol &&
	      (!(refs[ki].start2 > curr_ref.end2+tol ||
		 curr_ref.start2 > refs[ki].end2+tol)))
	    {
	      // Merge new knots
	      refs[ki].start2 = std::min(refs[ki].start2, curr_ref.start2);
	      refs[ki].end2 = std::max(refs[ki].end2, curr_ref.end2);
	      break;
	    }
	  else if (fabs(refs[ki].start2-curr_ref.start2) < tol &&
		   fabs(refs[ki].end2-curr_ref.end2) < tol &&
		   (!(refs[ki].start1 > curr_ref.end1+tol ||
		      curr_ref.start1 > refs[ki].end1+tol)))
	    {
	      // Merge new knots
	      refs[ki].start1 = std::min(refs[ki].start1, curr_ref.start1);
	      refs[ki].end1 = std::max(refs[ki].end1, curr_ref.end1);
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
void DefineRefs3D::refineStructuredMesh(const LRSplineVolume& vol,
					LRBSpline3D* bspline, Direction3D fixdir,
					int mult,
					vector<LRSplineVolume::Refinement3D>& refs)
//==============================================================================
{
  double tol = vol.getKnotTol();
  Direction3D d2 = next(fixdir);
  Direction3D d3 = prev(fixdir);
  int size = bspline->degree(fixdir)+1;
  const Mesh3D* mesh = bspline->getMesh();
  double pmin1 = bspline->pmin(d2);
  double pmax1 = bspline->pmax(d2);
  double pmin2 = bspline->pmin(d3);
  double pmax2 = bspline->pmax(d3);
  const vector<int>& kvec = bspline->kvec(fixdir);
  for (int kj=0; kj<size; ++kj)
    {
      double p1 = mesh->kval(fixdir, kvec[kj]);
      double p2 = mesh->kval(fixdir, kvec[kj+1]);
      if (p2-p1 <= tol)
	continue;
      LRSplineVolume::Refinement3D curr_ref;
      curr_ref.setVal(0.5*(p1+p2), pmin1, pmax1, pmin2, pmax2, fixdir, mult);
      appendRef(refs, curr_ref, tol);
    }
}

//==============================================================================
//  Full span refinement given parameter at which to refine
//
void DefineRefs3D::refineFullSpan(const LRSplineVolume& vol,
				  double upar, double vpar, double wpar,
				  Direction3D fixdir, int mult,
				  LRSplineVolume::Refinement3D& refs)
//==============================================================================
{
  Direction3D d2 = next(fixdir);
  Direction3D d3 = prev(fixdir);
  Element3D *elem = vol.coveringElement(upar, vpar, wpar);
  if (!elem)
    return;

  double ppar = (fixdir == XDIR) ? upar : ((fixdir == YDIR) ? vpar : wpar);
  double pmin1 = elem->pmin(d2);
  double pmax1 = elem->pmax(d2);
  double pmin2 = elem->pmin(d3);
  double pmax2 = elem->pmax(d3);

  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();
  
  // All overlapping B-splines
  for (size_t ki=0; ki<nmb; ++ki)
    {
      double bmin1 = bsplines[ki]->pmin(d2);
      double bmax1 = bsplines[ki]->pmax(d2);
      double bmin2 = bsplines[ki]->pmin(d3);
      double bmax2 = bsplines[ki]->pmax(d3);

      pmin1 = std::min(pmin1, bmin1);
      pmax1 = std::max(pmax1, bmax1);
      pmin2 = std::min(pmin2, bmin2);
      pmax2 = std::max(pmax2, bmax2);
    }

  refs.setVal(ppar, pmin1, pmax1, pmin2, pmax2, fixdir, mult);
}
 
//==============================================================================
//  Element based refinement 
//  strategies: fullSpan, minSpanSize (minimum span based on B-spline support),
//  minSpanEror (minimum span based on error in associated data points),
//  minSpanCombined (minimum span based on combination fof size and error)
void DefineRefs3D::refineFromElement(const LRSplineVolume& vol,
				     Element3D *elem, Direction3D fixdir,
				     RefineStrategy3D strategy, int mult,
				     LRSplineVolume::Refinement3D& refs)
//==============================================================================
{
  Direction3D d2 = next(fixdir);
  Direction3D d3 = prev(fixdir);
  
  // Fetch B-splines
  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  double ppar = 0.5*(elem->pmin(fixdir) + elem->pmax(fixdir));
  double pmin1 = elem->pmin(d2);
  double pmax1 = elem->pmax(d2);
  double pmin2 = elem->pmin(d3);
  double pmax2 = elem->pmax(d3);

  double par2 = 0.5*(pmin1 + pmax1); 
  double par3 = 0.5*(pmin2 + pmax2); 
  if (strategy == fullSpan)
    {
      // All overlapping B-splines
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  double bmin1 = bsplines[ki]->pmin(d2);
	  double bmax1 = bsplines[ki]->pmax(d2);
	  double bmin2 = bsplines[ki]->pmin(d3);
	  double bmax2 = bsplines[ki]->pmax(d3);

	  pmin1 = std::min(pmin1, bmin1);
	  pmax1 = std::max(pmax1, bmax1);
	  pmin2 = std::min(pmin2, bmin2);
	  pmax2 = std::max(pmax2, bmax2);
 	}
    }
  else if (strategy == minSpanSize)
    {
      // "Largest" overlapping B-spline. Weighted with mid split
      double tol = vol.getKnotTol();
      double max_size = 0.0;
      int ix = -1;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Compute size of B-spline
	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin1 = bsplines[ki]->pmin(d2);
	  double bmax1 = bsplines[ki]->pmax(d2);
	  double bmin2 = bsplines[ki]->pmin(d3);
	  double bmax2 = bsplines[ki]->pmax(d3);

	  double bsize = 2.0*(bmax0 - bmin0) + (bmax1 - bmin1) + (bmax2 - bmin2);
	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin1, bmax1-par2)/std::max(par2-bmin1, bmax1-par2);
	  double frac2 = std::min(par3-bmin2, bmax2-par3)/std::max(par3-bmin2, bmax2-par3);
	  double frac = 2.0*frac0 + frac1 + frac2;

	  if (frac*bsize > max_size)
	    {
	      max_size = frac*bsize;
	      ix = (int)ki;
	    }
	}
      pmin1 = bsplines[ix]->pmin(d2);
      pmax1 = bsplines[ix]->pmax(d2);
      pmin2 = bsplines[ix]->pmin(d3);
      pmax2 = bsplines[ix]->pmax(d3);
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
	  const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = curr_el[kj]->pmin(fixdir);
	      double emax = curr_el[kj]->pmax(fixdir);
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin1 = bsplines[ki]->pmin(d2);
	  double bmax1 = bsplines[ki]->pmax(d2);
	  double bmin2 = bsplines[ki]->pmin(d3);
	  double bmax2 = bsplines[ki]->pmax(d3);

	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin1, bmax1-par2)/std::max(par2-bmin1, bmax1-par2);
	  double frac2 = std::min(par3-bmin2, bmax2-par3)/std::max(par3-bmin2, bmax2-par3);
	  double frac = 2.0*frac0 + frac1 + frac2;

	  if ((fabs(max_wgt-curr_wgt) < tol && frac > max_frac) ||
	      curr_wgt > max_wgt)
	    {
	      max_wgt = curr_wgt;
	      max_frac = frac;
	      ix = (int)ki;
	    }
	}
      pmin1 = bsplines[ix]->pmin(d2);
      pmax1 = bsplines[ix]->pmax(d2);
      pmin2 = bsplines[ix]->pmin(d3);
      pmax2 = bsplines[ix]->pmax(d3);
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
	  const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = curr_el[kj]->pmin(fixdir);
	      double emax = curr_el[kj]->pmax(fixdir);
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  double bmin0 = bsplines[ki]->pmin(fixdir);
	  double bmax0 = bsplines[ki]->pmax(fixdir);
	  double bmin1 = bsplines[ki]->pmin(d2);
	  double bmax1 = bsplines[ki]->pmax(d2);
	  double bmin2 = bsplines[ki]->pmin(d3);
	  double bmax2 = bsplines[ki]->pmax(d3);
	  double bsize = 2.0*(bmax0 - bmin0) + (bmax1 - bmin1) + (bmax2 - bmin2);

	  double frac0 = std::min(ppar-bmin0, bmax0-ppar)/std::max(ppar-bmin0, bmax0-ppar);
	  double frac1 = std::min(par2-bmin1, bmax1-par2)/std::max(par2-bmin1, bmax1-par2);
	  double frac2 = std::min(par3-bmin2, bmax2-par3)/std::max(par3-bmin2, bmax2-par3);
	  double frac = 2.0*frac0 + frac1 + frac2;
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
      
      pmin1 = bsplines[ix]->pmin(d2);
      pmax1 = bsplines[ix]->pmax(d2);
      pmin2 = bsplines[ix]->pmin(d3);
      pmax2 = bsplines[ix]->pmax(d3);
   }
  
  refs.setVal(ppar, pmin1, pmax1, pmin2, pmax2, fixdir, mult);
}
