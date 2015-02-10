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


#include "GoTools/lrsplines2D/LRSurfStitch.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <iostream> // @@ debug
#include <fstream> // @@ debug

#define DEBUG

using namespace Go;
using std::vector;
using std::pair;

//==============================================================================
int LRSurfStitch::averageCorner(vector<pair<shared_ptr<ParamSurface>,int> >& sfs,
				double tol)
//==============================================================================
{
  int nmb_mod = (int)sfs.size(); // Initially all surfaces are expected to be handled
  int dim = sfs[0].first->dimension();

  // Pass through the surfaces and extract LR B-spline surfaces
  size_t kr, kh;
  vector<shared_ptr<LRSplineSurface> > lrsfs(sfs.size());
  for (kr=0; kr<sfs.size(); ++kr)
    {
      if (sfs[kr].first->dimension() != dim)
	{
	  return 0;
	}

      lrsfs[kr] = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sfs[kr].first);
      if (!lrsfs[kr].get())
	{
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[kr].first);
	  if (bd_sf.get())
	    {
	      shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	      lrsfs[kr] = 
		dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	    }
	}
    }

  // Remove the surfaces that are too far from a common corner to be modified
  // For a 1D surface, the distance is computed in the parameter domain of the surface
  // otherwise in geometry space
  vector<int> nmb_match(sfs.size(), 0);
  vector<Point> corner(sfs.size());
  vector<Point> corner_val(sfs.size());
  vector<Point> param(sfs.size());
  for (kr=0; kr<lrsfs.size(); ++kr)
    {
      if (!lrsfs[kr].get())
	{
	  nmb_mod--;
	  continue;
	}

      int corner_pos = sfs[kr].second;
      double u = (corner_pos == 0 || corner_pos == 2) ? lrsfs[kr]->startparam_u() :
	lrsfs[kr]->endparam_u();
      double v = (corner_pos == 0 || corner_pos == 1) ? lrsfs[kr]->startparam_v() :
	lrsfs[kr]->endparam_v();
      Point pos = lrsfs[kr]->ParamSurface::point(u,v);
      if (dim == 1)
	corner[kr] = Point(u,v);
      else
	corner[kr] = pos;
      corner_val[kr] = pos;
      param[kr] = Point(u,v);
    }

  for (kr=0; kr<lrsfs.size()-1; ++kr)
    {
      if (!lrsfs[kr].get())
	continue;
      for (kh=kr+1; kh<lrsfs.size(); ++kh)
	{
	  if (!lrsfs[kh].get())
	    continue;

	  // Check distance
	  double dist = corner[kr].dist(corner[kh]);
	  if (dist <= tol)
	    {
	      nmb_match[kr]++;
	      nmb_match[kh]++;
	    }
	}
    }

  // Compute average corner
  for (kr=0; kr<nmb_match.size(); ++kr)
    {
      if (nmb_match[kr] == 0)
	nmb_mod--;
    }
   for (kr=0; kr<nmb_match.size(); ++kr)
    {
      if (nmb_match[kr] > 0 && nmb_match[kr] < nmb_mod-1)
	break;
    }
   if (kr < nmb_match.size())
     {
       // Several disjunct corners. Special treatment
       // Currently, no modification
       return 0;
     }

   if (nmb_mod == 0)
     return 0;  // Nothing to do

   Point av_corner(dim);
   av_corner.setValue(0.0);
   for (kr=0; kr<lrsfs.size(); ++kr)
     {
       if (lrsfs[kr].get() && nmb_match[kr] > 0)
	 av_corner += corner_val[kr];
     }
   av_corner /= (double)nmb_mod;

  // Replace corner coefficient
   for (kr=0; kr<lrsfs.size(); ++kr)
     {
       if (lrsfs[kr].get() && nmb_match[kr] > 0)
	 {
	   // Identify B-spline
	   vector<LRBSpline2D*> cand_bsplines = 
	     lrsfs[kr]->basisFunctionsWithSupportAt(param[kr][0], param[kr][1]);
	   for (kh=0; kh<cand_bsplines.size(); ++kh)
	     {
	       int deg1 = cand_bsplines[kh]->degree(XFIXED);
	       int deg2 = cand_bsplines[kh]->degree(YFIXED);
	       int mult1 = cand_bsplines[kh]->endmult_u(true);
	       int mult2 = cand_bsplines[kh]->endmult_u(false);
	       int mult3 = cand_bsplines[kh]->endmult_v(true);
	       int mult4 = cand_bsplines[kh]->endmult_v(false);
	       if ((mult1 == deg1+1 || mult2 == deg1+1) &&
		   (mult3 == deg2+1 || mult4 == deg2+1))
		 break;
	     }

	   if (kh == cand_bsplines.size())
	     {
	       nmb_mod--;
	       continue;
	     }

	   // Replace coefficient (multiple knots at boundaries are assumed)
	   double gamma = cand_bsplines[kh]->gamma();
	   cand_bsplines[kh]->setCoefAndGamma(av_corner, gamma);
	 }
     }
   return nmb_mod;
}

//==============================================================================
bool LRSurfStitch::averageEdge(shared_ptr<ParamSurface> surf1, int edge1,
			       shared_ptr<ParamSurface> surf2, int edge2, 
			       double tol)
//==============================================================================
{
  int dim = surf1->dimension();
  if (surf2->dimension() != dim)
    return false;

 // Extract LR B-spline surfaces
  shared_ptr<LRSplineSurface> lrsf1 = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(surf1);
  if (!lrsf1.get())
    {
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf1);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	  lrsf1 = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	}
    }
  if (!lrsf1.get())
    return false;

  shared_ptr<LRSplineSurface> lrsf2 = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(surf2);
  if (!lrsf2.get())
    {
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf2);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	  lrsf2 = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	}
    }
  if (!lrsf2.get())
    return false;

  return averageEdge(lrsf1, edge1, lrsf2, edge2, tol);
}

//==============================================================================

bool bspline_sort_in_u(LRBSpline2D* bspline1, LRBSpline2D* bspline2)
{
  //std::cout << 1 << " " << bspline1 << " " << bspline2 << std::endl;
  if (bspline1->umin() == bspline2->umin())
    return (bspline1->umax() < bspline2->umax());
  return (bspline1->umin() < bspline2->umin());
}

bool bspline_sort_in_v(LRBSpline2D* bspline1, LRBSpline2D* bspline2)
{
  //std::cout << 2 << " " << bspline1 << " " << bspline2 << std::endl;
  if (bspline1->vmin() == bspline2->vmin())
    return (bspline1->vmax() < bspline2->vmax());
  return (bspline1->vmin() < bspline2->vmin());
}

//==============================================================================
bool LRSurfStitch::averageEdge(shared_ptr<LRSplineSurface> surf1, int edge1,
			       shared_ptr<LRSplineSurface> surf2, int edge2, 
			       double tol)
//==============================================================================
{
  int dim = surf1->dimension();
  if (surf2->dimension() != dim)
    return false;

  // Evaluate surface corners and check consistency
  double upar[4], vpar[4];
  Point corner[4];
  
  fetchEdgeCorners(surf1, edge1, upar[0], vpar[0], upar[1], vpar[1]);
  fetchEdgeCorners(surf2, edge2, upar[2], vpar[2], upar[3], vpar[3]);

  for (int ki=0; ki<4; ++ki)
    {
      if (dim == 1)
	corner[ki] = Point(upar[ki], vpar[ki]);
      else
	corner[ki] = (ki<2) ? surf1->ParamSurface::point(upar[ki], vpar[ki]) :
	  surf2->ParamSurface::point(upar[ki], vpar[ki]);
    }

  if (corner[0].dist(corner[2]) > tol && corner[0].dist(corner[3]) > tol)
    return false;
  if (corner[1].dist(corner[2]) > tol && corner[1].dist(corner[3]) > tol)
    return false;

  // Fetch knot vectors along the edge
  // Edges are numbered: 0=left, 1=right, 2=lower, 3=upper
  const Mesh2D& m1 = surf1->mesh();
  const Mesh2D& m2 = surf2->mesh();
  Direction2D dir1 = (edge1 == 0 || edge1 == 1) ? YFIXED : XFIXED;
  Direction2D dir2 = (edge2 == 0 || edge2 == 1) ? YFIXED : XFIXED;

  // Raise surface degree if necessary
  // This functionality is not implemented for LRSplineSurface. Must be done
  // if this function is to be used in a more general setting
  if (surf1->degree(dir1) != surf2->degree(dir2))
    return false;

  int ix1 = (edge1 == 0 || edge1 == 2) ? 0 : m1.numDistinctKnots(flip(dir1)) - 1;
  int ix2 = (edge2 == 0 || edge2 == 2) ? 0 : m2.numDistinctKnots(flip(dir2)) - 1;

  // The knot vectors are included multiplicities (no value at zero multiplicity)
  vector<vector<double> > knots(2);
  knots[0] = m1.getKnots(dir1, ix1, (ix1==0));
  knots[1] = m2.getKnots(dir2, ix2, (ix2==0));

  if (corner[0].dist(corner[2]) > corner[0].dist(corner[3]))
    {
      // The direction of the boundary curves are turned. Turn knot vector
    }

  double knot_tol = 1.0e-4;   // Should relate to existing knot intervals
  if (fabs(knots[0][0]-knots[1][0]) > knot_tol || 
      fabs(knots[0][knots[0].size()-1]-knots[1][knots[1].size()-1]) > knot_tol)
    {
      // Knot intervals are not matching. Change interval of second knot vector
    }

  // Compute knots to insert
  // First compute union knot vector
  vector<double> union_knots;
  GeometryTools::makeUnionKnots(knots, knot_tol, union_knots);

  // Extract the knots existing in the union vector, but not in the two knot vectors
  int order = surf1->degree(dir1) + 1;
  vector<double> new_knots1;
  extractMissingKnots(union_knots, knots[0], knot_tol, order, new_knots1);
  
  vector<double> new_knots2;
  extractMissingKnots(union_knots, knots[1], knot_tol, order, new_knots2);
  
  // If knot vector 2 has been altered, represent the corresponding missing
  // knots in the original interval

  // Define end parameters of knot intervals to insert. Set up refinement info
  int element_width = 2;
  vector<LRSplineSurface::Refinement2D> refs1;
  defineRefinements(m1, dir1, edge1, ix1, new_knots1, element_width, refs1);

  vector<LRSplineSurface::Refinement2D> refs2;
  defineRefinements(m2, dir2, edge2, ix2, new_knots2, element_width, refs2);

  // To ensure equally sized corresponding B-spline domains along the boundary
  // curve after refinement, let the knot line closest to the boundary traverse
  // the entire domain
  double par1 = m1.kval(flip(dir1), (ix1==0) ? ix1+1 : ix1-1);
  vector<pair<int, int> > zero_segm1 =
    m1.zeroSegments(flip(dir1), (ix1==0) ? ix1+1 : ix1-1);
  for (size_t kj=0; kj<zero_segm1.size(); ++kj)
    {
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(par1, m1.kval(dir1, zero_segm1[kj].first),
		      m1.kval(dir1, zero_segm1[kj].second), 
		      flip(dir1), 1);
      refs1.push_back(curr_ref);
    }

  double par2 = m2.kval(flip(dir2), (ix2==0) ? ix2+1 : ix2-1);
  vector<pair<int, int> > zero_segm2 =
    m2.zeroSegments(flip(dir2), (ix2==0) ? ix2+1 : ix2-1);
  for (size_t kj=0; kj<zero_segm2.size(); ++kj)
    {
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(par2, m2.kval(dir2, zero_segm2[kj].first),
		      m2.kval(dir2, zero_segm2[kj].second), 
		      flip(dir2), 1);
      refs2.push_back(curr_ref);
    }

  // Refine surfaces
 #ifdef DEBUG
  std::ofstream ofmesh1("mesh1.eps");
  writePostscriptMesh(*surf1, ofmesh1);
  std::ofstream ofmesh2("mesh2.eps");
  writePostscriptMesh(*surf2, ofmesh2);
#endif
 
  for (size_t kj=0; kj<refs1.size(); ++kj)
    surf1->refine(refs1[kj]);

  for (size_t kj=0; kj<refs2.size(); ++kj)
    surf2->refine(refs2[kj]);

 #ifdef DEBUG
  std::ofstream ofmesh3("mesh1_2.eps");
  writePostscriptMesh(*surf1, ofmesh3);
  std::ofstream ofmesh4("mesh2_2.eps");
  writePostscriptMesh(*surf2, ofmesh4);
#endif

  // Traverse the boundary B-splines of the two surfaces and compute the
  // average coefficients
  // Assumes now that the surfaces are oriented equally
  
  // Fetch B-splines along edges of surfaces
  vector<LRBSpline2D*> bsplines1;
  extractBoundaryBsplines(surf1, edge1, bsplines1);

  vector<LRBSpline2D*> bsplines2;
  extractBoundaryBsplines(surf2, edge2, bsplines2);

  // Traverse B-splines and compute average coefficients
  vector<double> knots1 = 
    m1.getKnots(dir1, (ix1==0) ? 0 : m1.numDistinctKnots(flip(dir1)) - 1, (ix1==0));
  vector<double> knots2 = 
    m2.getKnots(dir2, (ix2==0) ? 0 : m2.numDistinctKnots(flip(dir2)) - 1, (ix2==0));
  if (bsplines1.size() != bsplines2.size())
    return false;
  
  // Sort B-splines to ensure corresponding sequences
  std::sort(bsplines1.begin(), bsplines1.end(), (dir1==XFIXED) ? bspline_sort_in_u :
	    bspline_sort_in_v);
  std::sort(bsplines2.begin(), bsplines2.end(), (dir2==XFIXED) ? bspline_sort_in_u :
  	    bspline_sort_in_v);

#ifdef DEBUG
  std::ofstream of("bspline.txt");
#endif
  for (size_t kj=0; kj<bsplines1.size(); ++kj)
    {
      Point coef1 = bsplines1[kj]->Coef();
      double gamma1 = bsplines1[kj]->gamma();
      Point coef2 = bsplines2[kj]->Coef();
      double gamma2 = bsplines2[kj]->gamma();
      Point coef = 0.5*(coef1 + coef2);
#ifdef DEBUG
      of << kj << std::endl;
      of << " bspline1: " << bsplines1[kj]->umin() << "  " << bsplines1[kj]->umax();
      of << " " << bsplines1[kj]->vmin() << "  " << bsplines1[kj]->vmax();
      of << " " << coef1[0] << " " << gamma1 << std::endl;
      of << " bspline2: "  << bsplines2[kj]->umin() << "  " << bsplines2[kj]->umax();
      of << " " << bsplines2[kj]->vmin() << "  " << bsplines2[kj]->vmax();
      of << " " << coef2[0] << " " << gamma2 << std::endl << std::endl;
#endif
      bsplines1[kj]->setCoefAndGamma(coef, gamma1);
      bsplines2[kj]->setCoefAndGamma(coef, gamma2);
    }

  return true;
}


//==============================================================================
void LRSurfStitch::fetchEdgeCorners(shared_ptr<LRSplineSurface> surf, int edge,
				    double& u1, double& v1, double& u2, double& v2)
//==============================================================================
{
  // Edges are numbered: 0=left, 1=right, 2=lower, 3=upper
  if (edge <= 1)
    u1 = u2 = (edge == 0) ? surf->startparam_u() : surf->endparam_u();
  else
    {
      u1 = surf->startparam_u();
      u2 = surf->endparam_u();
    }
  if (edge >= 2)
    v1 = v2 = (edge == 2) ? surf->startparam_v() : surf->endparam_v();
  else
    {
      v1 = surf->startparam_v();
      v2 = surf->endparam_v();
    }
}

//==============================================================================
void LRSurfStitch::extractMissingKnots(vector<double>& union_vec, 
				       vector<double>& vec,
				       double tol, int order,
				       vector<double>& resvec)
//==============================================================================
{
  int ki, kj;
  int size1 = (int)vec.size() - order;
  int size2 = (int)union_vec.size() - order;
  for (ki=order, kj=order; ki<size1 || kj<size2; )
    {
      if (fabs(vec[ki]-union_vec[kj]) < tol)
	{
	  ki++;
	  kj++;
	}
      else if (union_vec[kj] < vec[ki])
	{
	  resvec.push_back(union_vec[kj]);
	  kj++;
	}
      else
	ki++;
    }
}

//==============================================================================
void LRSurfStitch::defineRefinements(const Mesh2D& mesh, Direction2D dir,
				     int edge, int ix, vector<double>& knot_vals, 
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
	    c_ix =
	      Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh, curr_dir,
								     c_ix-1, other_ix);
	  p2 = mesh.kval(curr_dir, c_ix);
	}
      else
	{
	  int c_ix = ix;
	  for (int ki=0; ki<element_width; ++ki)
	    c_ix =
	      Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh, curr_dir,
								   c_ix+1, other_ix);
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
void LRSurfStitch::extractBoundaryBsplines(shared_ptr<LRSplineSurface> surf,
					   int edge,
					   vector<LRBSpline2D*>& bsplines)
//==============================================================================
{
  // edge is given as: 0=left, 1=right, 2=lower, 3=upper

  // Traverse B-splines and extract those that have k-tupple multiplicity
  // at the specified edge. Since the B-splines are sorted according to the
  // lower left corner primarily and the upper right corner secondarly, the
  // B-splines will be order from left to right or lower to upper

  bsplines.clear();
  bool at_start = (edge == 0 || edge == 2);
  for (LRSplineSurface::BSplineMap::iterator it=surf->basisFunctionsBeginNonconst(); 
       it != surf->basisFunctionsEndNonconst(); ++it)
	{
	  int deg = it->second->degree((edge == 0 || edge == 1) ? YFIXED : XFIXED);
	  int mult = (edge <= 1) ? it->second->endmult_u(at_start) :
	    it->second->endmult_v(at_start);

	  if (mult == deg+1)
	    bsplines.push_back(it->second.get());
	}
}
  
