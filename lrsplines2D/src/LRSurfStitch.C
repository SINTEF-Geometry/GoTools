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
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <iostream> // @@ debug
#include <fstream> // @@ debug

//#define DEBUG

using namespace Go;
using std::vector;
using std::set;
using std::pair;

//==============================================================================
void LRSurfStitch::stitchRegSfs(vector<shared_ptr<ParamSurface> >& sfs,
				int nmb_u, int nmb_v, double eps,
				int cont)
//==============================================================================
{
  // Represent as LR B-spline surfaces. Trimmed surfaces are replaces
  // by their undrlying surface
  vector<shared_ptr<LRSplineSurface> > lr_sfs(sfs.size());
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      if (!sfs[ki].get())
	continue;
      shared_ptr<LRSplineSurface> sflr =
	dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sfs[ki]);
      if (!sflr.get())
	{
	  shared_ptr<BoundedSurface> sfbd = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[ki]);
	  if (sfbd.get())
	    {
	      shared_ptr<ParamSurface> tmp_sf = sfbd->underlyingSurface();
	      sflr = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	    }
	}
      
      if (sflr.get())
	lr_sfs[ki] = sflr;
    }

  stitchRegSfs(lr_sfs, nmb_u, nmb_v, eps, cont);
}

//==============================================================================
void LRSurfStitch::stitchRegSfs(vector<shared_ptr<LRSplineSurface> >& sfs,
				int nmb_u, int nmb_v, double eps,
				int cont)
//==============================================================================
{
#ifdef DEBUG
  std::ofstream ofmesh1("mesh_1.eps");
  writePostscriptMesh(*sfs[0], ofmesh1);
  std::ofstream ofmesh2("mesh_2.eps");
  writePostscriptMesh(*sfs[1], ofmesh2);
  if (sfs.size() > 2)
    {
      std::ofstream ofmesh3("mesh_3.eps");
      writePostscriptMesh(*sfs[2], ofmesh3);
    }
  if (sfs.size() > 3)
    {
      std::ofstream ofmesh4("mesh_4.eps");
      writePostscriptMesh(*sfs[3], ofmesh4);
    }
  if (sfs.size() > 4)
    {
      std::ofstream ofmesh5("mesh_5.eps");
      writePostscriptMesh(*sfs[4], ofmesh5);
    }
  if (sfs.size() > 5)
    {
      std::ofstream ofmesh6("mesh_6.eps");
      writePostscriptMesh(*sfs[5], ofmesh6);
    }
#endif
  // Ensure corresponding spline spaces along common boundaries
  // We also make sure that the surfaces are full tensor product surfaces
  // along adjacent edges (the first couple of element rows).
  // Note that the surfaces (i.e. the coefs) are not altered.
  consistentSplineSpaces(sfs, nmb_u, nmb_v, eps, cont);
  consistentSplineSpaces(sfs, nmb_u, nmb_v, eps, cont);
  

  // Stitch surfaces along common edges (by altering the coefs).
  // Assosiated corners will be handled first
  int kj, kr;
  int nmb_modified;
  bool matched = false;
  for (kj=0; kj<nmb_v; ++kj)
    {
      for (kr=0; kr<=nmb_u; ++kr)
      {
	if (kj == 0 && kr > 0 && kr < nmb_u &&
	    sfs[kr-1].get() && sfs[kr].get())
	  {
	    // Match corner along the lower boundary
	    // Corners are numbered: 0=lower left, 1=lower right,
	    // 2=upper left, 3=upper right
	    vector<pair<shared_ptr<LRSplineSurface>, int> > corner_match1(2);
	    corner_match1[0] = make_pair(sfs[kr-1], 1);
	    corner_match1[1] = make_pair(sfs[kr], 0);
	    if (cont == 0)
	      nmb_modified = averageCorner(corner_match1, eps);
	    else
	      nmb_modified = makeCornerC1(corner_match1, eps);
	  }
	if (kj == nmb_v-1 && kr > 0 && kr < nmb_u &&
	    sfs[kj*nmb_u+kr-1].get() && sfs[kj*nmb_u+kr].get())
	  {
	    // Match corner along the upper boundary
	    // Corners are numbered: 0=lower left, 1=lower right,
	    // 2=upper left, 3=upper right
	    vector<pair<shared_ptr<LRSplineSurface>, int> > corner_match2(2);
	    corner_match2[0] = make_pair(sfs[kj*nmb_u+kr-1], 3);
	    corner_match2[1] = make_pair(sfs[kj*nmb_u+kr], 2);
	    if (cont == 0)
	      nmb_modified = averageCorner(corner_match2, eps);
	    else
	      nmb_modified = makeCornerC1(corner_match2, eps);
	  }
	
	// Match inner corner
	if (kj > 0)
	  {
	    vector<pair<shared_ptr<LRSplineSurface>, int> > corner_match3;
	    if (kr > 0 && sfs[(kj-1)*nmb_u+kr-1].get())
	      corner_match3.push_back(make_pair(sfs[(kj-1)*nmb_u+kr-1], 3));
	    if (kr < nmb_u && sfs[(kj-1)*nmb_u+kr].get())
	      corner_match3.push_back(make_pair(sfs[(kj-1)*nmb_u+kr], 2));
	    if (kr > 0 && sfs[kj*nmb_u+kr-1].get())
	      corner_match3.push_back(make_pair(sfs[kj*nmb_u+kr-1], 1));
	    if (kr < nmb_u && sfs[kj*nmb_u+kr].get())
	      corner_match3.push_back(make_pair(sfs[kj*nmb_u+kr], 0));
	    if (corner_match3.size() > 1)
	      {
		if (cont == 0)
		  nmb_modified = averageCorner(corner_match3, eps);
		else
		  nmb_modified = makeCornerC1(corner_match3, eps);
	      
	      }

	    // Match vertical edges
	    // Edges are numbered: 0=left, 1=right, 2=lower, 3=upper
	    if (kr > 0 && kr < nmb_u &&
		sfs[(kj-1)*nmb_u+kr-1].get() && sfs[(kj-1)*nmb_u+kr].get())
	      {
		matched = averageEdge(sfs[(kj-1)*nmb_u+kr-1], 1,
				      sfs[(kj-1)*nmb_u+kr], 0, 
				      cont, eps);
		if (!matched)
		  {
#ifdef DEBUG
		    std::cout << "Failed vertical match! ki = " << kr-1 << ", kj = " << kj - 1 << std::endl;
		    //	      MESSAGE("Failed vertical match!");
#endif
		  }
	      }
	  }
	if (kj == nmb_v-1 && kr > 0 && kr < nmb_u &&
	    sfs[kj*nmb_u+kr-1].get() && sfs[kj*nmb_u+kr].get())
	{
	  matched = averageEdge(sfs[kj*nmb_u+kr-1], 1,
				sfs[kj*nmb_u+kr], 0, 
				cont, eps);
	  if (!matched)
	  {
#ifdef DEBUG
	      std::cout << "Failed vertical match! ki = " << kr-1 << ", kj = " << kj << std::endl;
//	      MESSAGE("Failed vertical match!");
#endif
	  }
	}

	// Match horizontal edge
	if (kj > 0)
	  {
	    if (kr < nmb_u &&
		sfs[(kj-1)*nmb_u+kr].get() && sfs[kj*nmb_u+kr].get())
	      {
		matched = averageEdge(sfs[(kj-1)*nmb_u+kr], 3,
				      sfs[kj*nmb_u+kr], 2, 
				      cont, eps);
		if (!matched)
		  {
#ifdef DEBUG
		    std::cout << "Failed horizontal match! ki = " << kr << ", kj = " << kj-1 << std::endl;
		    //	      MESSAGE("Failed horizontal match!");
#endif
		  }
	      }
	  }
      }
    }
#ifdef DEBUG
  std::ofstream ofmesh_1("mesh2_1.eps");
  writePostscriptMesh(*sfs[0], ofmesh_1);
  std::ofstream ofmesh_2("mesh2_2.eps");
  writePostscriptMesh(*sfs[1], ofmesh_2);
  if (sfs.size() > 2)
    {
      std::ofstream ofmesh_3("mesh2_3.eps");
      writePostscriptMesh(*sfs[2], ofmesh_3);
    }
  if (sfs.size() > 3)
    {
      std::ofstream ofmesh_4("mesh2_4.eps");
      writePostscriptMesh(*sfs[3], ofmesh_4);
    }
  if (sfs.size() > 4)
    {
      std::ofstream ofmesh_5("mesh2_5.eps");
      writePostscriptMesh(*sfs[4], ofmesh_5);
    }
  if (sfs.size() > 5)
    {
      std::ofstream ofmesh_6("mesh2_6.eps");
      writePostscriptMesh(*sfs[5], ofmesh_6);
    }
#endif
}

//==============================================================================
vector<double> LRSurfStitch::analyzeContinuity(vector<shared_ptr<ParamSurface> >& sfs,
					       int nmb_u, int nmb_v, int derivs,
					       int nmb_edge_samples)
//==============================================================================
{
  assert(derivs < 2);

  vector<shared_ptr<LRSplineSurface> > lr_sfs(sfs.size());
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      if (!sfs[ki].get())
	continue;
      shared_ptr<LRSplineSurface> sflr =
	dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sfs[ki]);
      if (!sflr.get())
	{
	  shared_ptr<BoundedSurface> sfbd = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[ki]);
	  if (sfbd.get())
	    {
	      shared_ptr<ParamSurface> tmp_sf = sfbd->underlyingSurface();
	      sflr = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	    }
	}
      
      if (sflr.get())
	lr_sfs[ki] = sflr;
    }

    vector<double> global_max_cont_dist(derivs + 1, -1.0);

    // We first analyze the cont in the nmb_u direction.
    for (int kj = 0; kj < nmb_v; ++kj)
    {
	for (int ki = 0; ki < nmb_u - 1; ++ki)
	{
	    vector<double> max_cont_dist(derivs + 1, -1.0);
	    vector<int> sample_ind(derivs + 1, -1);

	    LRSplineSurface* sf_left = lr_sfs[kj*nmb_u + ki].get();
	    LRSplineSurface* sf_right = lr_sfs[kj*nmb_u + ki + 1].get();

	    if (sf_left == NULL || sf_right == NULL)
	      continue;

	    // We may assume the the lr_sfs are oriented such that sf_left->umax matches sf_right->umin.
	    double upar_left = sf_left->endparam_u();
	    double upar_right = sf_right->startparam_u();

	    // We may also assume that they share domain along the shared edge.
	    double vmin = sf_left->startparam_v();
	    double vmax = sf_left->endparam_v();
	    double vstep = (vmax - vmin)/(nmb_edge_samples - 1);
	    vector<Point> pt_left((derivs + 1)*(derivs + 2));
	    vector<Point> pt_right((derivs + 1)*(derivs + 2));
	    for (int kk = 0; kk < nmb_edge_samples; ++kk)
	    {
		double vpar = vmin + kk*vstep;
		sf_left->point(pt_left, upar_left, vpar, derivs);
		sf_right->point(pt_right, upar_right, vpar, derivs);
		for (int kr = 0; kr < derivs + 1; ++kr)
		{
		    int ind = (kr == 1) ? 1 : 0;
		    double dist = pt_left[ind].dist(pt_right[ind]);
		    if (dist > max_cont_dist[kr])
		    {
			max_cont_dist[kr] = dist;
			sample_ind[kr] = kk;
		    }
		}
	    }

	    for (size_t kk = 0; kk < max_cont_dist.size(); ++kk)
	    {
		// std::cout << "ki: " << ki << ", kj: " << kj << ", u-dir, sample #: " << sample_ind[kk] << 
		//     ", max_cont_dist[" << kk << "]: " << max_cont_dist[kk] << std::endl;
		if (max_cont_dist[kk] > global_max_cont_dist[kk])
		{
		    global_max_cont_dist[kk] = max_cont_dist[kk];
		}
	    }
	}
    }

    // We then analyze the cont in the nmb_v direction.
    for (int kj = 0; kj < nmb_u; ++kj)
    {
	for (int ki = 0; ki < nmb_v - 1; ++ki)
	{
	    vector<double> max_cont_dist(derivs + 1, -1.0);
	    vector<int> sample_ind(derivs + 1, -1);

	    LRSplineSurface* sf_bottom = lr_sfs[ki*nmb_u + kj].get();
	    LRSplineSurface* sf_top = lr_sfs[(ki + 1)*nmb_u + kj].get();

	    if (sf_bottom == NULL || sf_top == NULL)
	      continue;

	    // We may assume the the lr_sfs are oriented such that sf_bottom->vmax matches sf_top->vmin.
	    double vpar_bottom = sf_bottom->endparam_v();
	    double vpar_top = sf_top->startparam_v();

	    // We may also assume that they share domain along the shared edge.
	    double umin = sf_bottom->startparam_u();
	    double umax = sf_bottom->endparam_u();
	    double ustep = (umax - umin)/(nmb_edge_samples - 1);
	    vector<Point> pt_bottom((derivs + 1)*(derivs + 2));
	    vector<Point> pt_top((derivs + 1)*(derivs + 2));
	    for (int kk = 0; kk < nmb_edge_samples; ++kk)
	    {
		double upar = umin + kk*ustep;
		sf_bottom->point(pt_bottom, upar, vpar_bottom, derivs);
		sf_top->point(pt_top, upar, vpar_top, derivs);
		for (int kr = 0; kr < derivs + 1; ++kr)
		{
		    int ind = (kr == 1) ? 2 : 0;
		    double dist = pt_bottom[ind].dist(pt_top[ind]);
		    if (dist > max_cont_dist[kr])
		    {
			max_cont_dist[kr] = dist;
			sample_ind[kr] = kk;
		    }
		}
	    }

	    for (size_t kk = 0; kk < max_cont_dist.size(); ++kk)
	    {
		// std::cout << "ki: " << ki << ", kj: " << kj << ", v-dir, sample #: " << sample_ind[kk] << 
		//     ", max_cont_dist[" << kk << "]: " << max_cont_dist[kk] << std::endl;
		if (max_cont_dist[kk] > global_max_cont_dist[kk])
		{
		    global_max_cont_dist[kk] = max_cont_dist[kk];
		}
	    }
	}
    }

    return global_max_cont_dist;
}

//==============================================================================
void LRSurfStitch::consistentSplineSpaces(vector<shared_ptr<LRSplineSurface> >& sfs,
					  int nmb_u, int nmb_v, double eps,
					  int cont)
//==============================================================================
{
  int kj, kr;

  int max_nmb = std::max(nmb_u, nmb_v);
  // The current method refines globally around surface corners. In theory this means that 'max_nmb - 1'
  // corner adjustments must be performed for the method to propagate along the nmb_u x nmb_v-grid.  Once
  // there is no change in the number of basis functions we exit.
  int sum_basis_functions = 0;
  for (int kk = 0; kk < sfs.size(); ++kk)
  {
      if (sfs[kk].get() != NULL)
      {
	  sum_basis_functions += sfs[kk]-> numBasisFunctions();
      }
  }
  for (int ki = 0; ki < max_nmb - 1; ++ki)
  {
      // First ensure tensor-product structure close to the boundaries
      bool edges[4];
      for (kj=0; kj<nmb_v; ++kj)
      {
	  edges[2] = (kj != 0);
	  edges[3] = (kj != nmb_v-1);
	  for (kr=0; kr<nmb_u; ++kr)
	  {
	      if (!sfs[kj*nmb_u+kr].get())
		  continue;
	      edges[0] = (kr != 0);
	      edges[1] = (kr != nmb_u-1);
	      // Other parts of the code requires 2 inner rows for c0.
	      int num_inner_rows = std::max(cont + 1, 2); // I.e. the number of elements that are affected.
	      tensorStructure(sfs[kj*nmb_u+kr], num_inner_rows, edges);
	  }
      }

      // Then make corresponding spline spaces across boundaries (i.e. insert knots along the common edge).
      for (kj=0; kj<nmb_v; ++kj)
      {
	  for (kr=0; kr<=nmb_u; ++kr)
	  {
	      // Vertical edge
	      // Edges are numbered: 0=left, 1=right, 2=lower, 3=upper
	      const int element_width = cont + 2;//std::max(2, cont + 1);//cont + 2; // Number of rows with inner knots.
	      if (kj > 0 && kr > 0 && kr < nmb_u &&
		  sfs[(kj-1)*nmb_u+kr-1].get() && sfs[(kj-1)*nmb_u+kr].get())
		  matchSplineSpace(sfs[(kj-1)*nmb_u+kr-1], 1,
				   sfs[(kj-1)*nmb_u+kr], 0, element_width, eps);
	      if (kj == nmb_v-1 && kr > 0 && kr < nmb_u &&
		  sfs[kj*nmb_u+kr-1].get() && sfs[kj*nmb_u+kr].get())
		  matchSplineSpace(sfs[kj*nmb_u+kr-1], 1,
				   sfs[kj*nmb_u+kr], 0, element_width, eps);

	      // Horizontal edge
	      if (kj > 0 && kr < nmb_u &&
		  sfs[(kj-1)*nmb_u+kr].get() && sfs[kj*nmb_u+kr].get())
		  matchSplineSpace(sfs[(kj-1)*nmb_u+kr], 3,
				   sfs[kj*nmb_u+kr], 2, element_width, eps);
	  }
      }

      int new_sum_basis_functions = 0;
      for (int kk = 0; kk < sfs.size(); ++kk)
      {
	  if (sfs[kk].get() != NULL)
	  {
	      new_sum_basis_functions += sfs[kk]-> numBasisFunctions();
	  }
      }
      if (new_sum_basis_functions == sum_basis_functions)
      {
	  break;
      }
      else
      {
	  sum_basis_functions = new_sum_basis_functions;
      }
  }
}

//==============================================================================
int LRSurfStitch::averageCorner(vector<pair<shared_ptr<ParamSurface>,int> >& sfs,
				double tol)
//==============================================================================
{
  int dim = sfs[0].first->dimension();

  // Pass through the surfaces and extract LR B-spline surfaces
  size_t kr;
  vector<pair<shared_ptr<LRSplineSurface>,int> > lrsfs(sfs.size());
  for (kr=0; kr<sfs.size(); ++kr)
    {
      if (sfs[kr].first->dimension() != dim)
	{
	  return 0;
	}

      shared_ptr<LRSplineSurface> sflr =
	dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sfs[kr].first);
      if (!sflr.get())
	{
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[kr].first);
	  if (bd_sf.get())
	    {
	      shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	      sflr = 
		dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
	    }
	}
      lrsfs[kr] = make_pair(sflr, sfs[kr].second);
    }

  return averageCorner(lrsfs, tol);
}

//==============================================================================
void LRSurfStitch::checkCornerMatch(vector<pair<shared_ptr<LRSplineSurface>,int> >& sfs,
				   double tol, vector<int>& nmb_match,
				   vector<Point>& corner_val, 
				   vector<Point>& param)
//==============================================================================
{
  // Remove the surfaces that are too far from a common corner to be modified
  // For a 1D surface, the distance is computed in the parameter domain of the surface
  // otherwise in geometry space
  size_t kr, kh;
  vector<Point> corner(sfs.size());
  for (kr=0; kr<sfs.size(); ++kr)
    {
      if (!sfs[kr].first.get())
	{
	  continue;
	}

      int corner_pos = sfs[kr].second;
      double u = (corner_pos == 0 || corner_pos == 2) ? sfs[kr].first->startparam_u() :
	sfs[kr].first->endparam_u();
      double v = (corner_pos == 0 || corner_pos == 1) ? sfs[kr].first->startparam_v() :
	sfs[kr].first->endparam_v();
      Point pos = sfs[kr].first->ParamSurface::point(u,v);
      if (sfs[kr].first->dimension() == 1)
	corner[kr] = Point(u,v);
      else
	corner[kr] = pos;
      corner_val[kr] = pos;
      param[kr] = Point(u,v);
    }

  for (kr=0; kr<sfs.size()-1; ++kr)
    {
      if (!sfs[kr].first.get())
	continue;
      for (kh=kr+1; kh<sfs.size(); ++kh)
	{
	  if (!sfs[kh].first.get())
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
}

//==============================================================================
int LRSurfStitch::averageCorner(vector<pair<shared_ptr<LRSplineSurface>,int> >& sfs,
				double tol)
//==============================================================================
{
  // Check if the given corners match
  vector<int> nmb_match(sfs.size(), 0);
  vector<Point> corner_val(sfs.size());
  vector<Point> param(sfs.size());
  checkCornerMatch(sfs, tol, nmb_match, corner_val, param);

  // Compute average corner
  size_t kr;
  int dim;
  int nmb_mod = (int)nmb_match.size();
  for (kr=0; kr<nmb_match.size(); ++kr)
    {
      if (nmb_match[kr] == 0)
	nmb_mod--;
      else
	dim = sfs[kr].first->dimension();
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
   for (kr=0; kr<sfs.size(); ++kr)
     {
       if (sfs[kr].first.get() && nmb_match[kr] > 0)
	 av_corner += corner_val[kr];
     }
   av_corner /= (double)nmb_mod;

  // Replace corner coefficient
   for (size_t kr=0; kr<sfs.size(); ++kr)
     {
       if (sfs[kr].first.get() && nmb_match[kr] > 0)
	 {
	   // Identify B-spline
	   vector<LRBSpline2D*> cand_bsplines = 
	     sfs[kr].first->basisFunctionsWithSupportAt(param[kr][0], param[kr][1]);
	   size_t kh;
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
int LRSurfStitch::makeCornerC1(vector<pair<shared_ptr<LRSplineSurface>,int> >& sfs,
			       double tol)
//==============================================================================
{
  // Check if the given corners match
  vector<int> nmb_match(sfs.size(), 0);
  vector<Point> corner_val(sfs.size());
  vector<Point> param(sfs.size());
  checkCornerMatch(sfs, tol, nmb_match, corner_val, param);

  // Check configuration
  int nmb_mod = (int)nmb_match.size();
  size_t kr;
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

  // Extract corner B-splines for all matching surfaces
  vector<LRBSpline2D*> bsplines;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      if (sfs[ki].first.get() && nmb_match[ki]>0)
	{
	  bool at_start1 = (sfs[ki].second == 0 || sfs[ki].second == 2) ? true :
	    false;
	  bool at_start2 = (sfs[ki].second <= 1 ) ? true : false;
	  int deg1 = sfs[ki].first->degree(XFIXED);
	  int deg2 = sfs[ki].first->degree(YFIXED);

	   vector<LRBSpline2D*> cand_bsplines = 
	     sfs[ki].first->basisFunctionsWithSupportAt(param[ki][0], param[ki][1]);
	   for (kr=0; kr<cand_bsplines.size(); ++kr)
	     {
	       int mult1 = cand_bsplines[kr]->endmult_u(at_start1);
	       int mult2 = cand_bsplines[kr]->endmult_v(at_start2);
	       if ((mult1 == deg1+1 || mult1 == deg1) &&
		   (mult2 == deg2+1 || mult2 == deg2))
		 bsplines.push_back(cand_bsplines[kr]);
	     }
	}
    }

  // Find least squares plane for the corner coefficients of all
  // surfaces (1D only)
  // Set up equation system
  double ls1=0, ls24=0, ls37=0, ls5=0, ls68=0;
  double hs1=0, hs2=0, hs3=0;
  double ls9 = 1.0;

  int kj;
  vector<Point> par(bsplines.size());
  vector<Point> coef(bsplines.size());
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      par[kj] = bsplines[kj]->getGrevilleParameter();
      double x = par[kj][0];
      double y = par[kj][1];
      coef[kj] = bsplines[kj]->Coef();
      double z = coef[kj][0];
      ls1 += (x*x);
      ls24 += (x*y);
      ls37 += x;
      ls5 += (y*y);
      ls68 += y;
      hs1 += (z*x);
      hs2 += (z*y);
      hs3 += z;
    }

  MatrixXD<double,3> M;
  M(0,0) = ls1;
  M(0,1) = ls24;
  M(0,2) = ls37;
  M(1,0) = ls24;
  M(1,1) = ls5;
  M(1,2) = ls68;
  M(2,0) = ls37;
  M(2,1) = ls68;
  M(2,2) = ls9;

  // Solve equation system by Cramers rule
  double det = M.det();
  double detk[3];
  for (kj=0; kj<3; ++kj)
    {
      MatrixXD<double,3> M1 = M;
      M1(kj,0) = hs1;
      M1(kj,1) = hs2;
      M1(kj,2) = hs3;
      detk[kj] = M1.det();
    }
  double a = detk[0]/det;
  double b = detk[1]/det;
  double c = detk[2]/det;

  // Project the surfaces onto the found plane
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      coef[kj][0] = a*par[kj][0] + b*par[kj][1] + c;
      double gamma = bsplines[kj]->gamma();
      bsplines[kj]->setCoefAndGamma(coef[kj], gamma);
    }

#if 0//def DEBUG
  {
      // We check if the corner continuity is satisfied.
      for (kj = 0; kj < (int)sfs.size() - 1; ++kj)
      {
	  const LRSplineSurface* sf1 = sfs[kj].first.get();
	  const int edge1 = sfs[kj].second;
	  const LRSplineSurface* sf2 = sfs[kj+1].first.get();
	  const int edge2 = sfs[kj+1].second;
	  const bool dir_is_u = (edge1 > 1);
	  double sf_dist = -1.0;
	  double tang1_ang = -1.0;
	  double tang2_ang = -1.0;
//	  double normal_ang = -1.0;
	  Point par1 = param[kj];
	  Point par2 = param[kj+1];
	  try
	  {
	      surfaceDifference(sf1, par1, sf2, par2,// dir_is_u,
				sf_dist, tang1_ang, tang2_ang);
	  }
	  catch (...)
	  {
	      MESSAGE("Fail ...");
	  }
	  std::cout << "sf_dist = " << sf_dist << ", tang1_ang = " << tang1_ang << ", tang2_ang = " << tang2_ang << std::endl;
      }

  }
#endif

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
void LRSurfStitch::tensorStructure(shared_ptr<LRSplineSurface> surf, 
				   int element_width, bool edges[4])
//==============================================================================
{
  const Mesh2D& mesh = surf->mesh();
  vector<LRSplineSurface::Refinement2D> refs;
  for (int ki=0; ki<4; ++ki)
    {
      if (!edges[ki])
	continue;    // Do not imply tensor product structure along this edge

      Direction2D dir = (ki <= 1) ? XFIXED : YFIXED;
      bool at_start = (ki ==0 || ki == 2);
      int nmb = mesh.numDistinctKnots(dir);
      for (int kj=1; kj<=element_width; ++kj)
	{
	  int ix = at_start ? kj : nmb-kj-1;
	  double par = mesh.kval(dir, ix);
	  vector<pair<int, int> > zero_segm = mesh.zeroSegments(dir, ix);
	  for (size_t kr=0; kr<zero_segm.size(); ++kr)
	    {
	      LRSplineSurface::Refinement2D curr_ref;
	      curr_ref.setVal(par, mesh.kval(flip(dir), zero_segm[kr].first),
			      mesh.kval(flip(dir), zero_segm[kr].second), 
			      dir, 1);
	      refs.push_back(curr_ref);
	    }
	}
    }

  // Refine surfaces
 #ifdef DEBUG
  std::ofstream ofmesh("mesh.eps");
  writePostscriptMesh(*surf, ofmesh);
#endif
 
  for (size_t kr=0; kr<refs.size(); ++kr)
    surf->refine(refs[kr]);

// #ifdef DEBUG
//   std::ofstream ofmesh_2("mesh_2.eps");
//   writePostscriptMesh(*surf, ofmesh_2);
// #endif

  return;
}

//==============================================================================
bool LRSurfStitch::matchSplineSpace(shared_ptr<LRSplineSurface> surf1, 
				    int edge1,
				    shared_ptr<LRSplineSurface> surf2,
				    int edge2, 
				    int element_width, double tol)
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
  int num_knot_rows = element_width;// + 1;
  vector<vector<double> > knots(2*(num_knot_rows));
  int sgn1 = (ix1 == 0) ? 1 : -1;
  int sgn2 = (ix2 == 0) ? 1 : -1;
  for (int ki=0; ki<num_knot_rows; ++ki)
    {
      knots[ki] = m1.getKnots(dir1, ix1+sgn1*ki, (ix1==0));
      knots[num_knot_rows+ki] = m2.getKnots(dir2, ix2+sgn2*ki, (ix2==0));
    }

  if (corner[0].dist(corner[2]) > corner[0].dist(corner[3]))
    {
      // The direction of the boundary curves are turned. Turn knot vector
	MESSAGE("We should turn the knot vector! Not implemented yet.");
    }

  double knot_tol = 1.0e-4;   // Should relate to existing knot intervals
  if (fabs(knots[0][0]-knots[num_knot_rows][0]) > knot_tol || 
      fabs(knots[0][knots[0].size()-1]-knots[num_knot_rows][knots[num_knot_rows].size()-1]) > knot_tol)
    {
      // Knot intervals are not matching. Change interval of second knot vector
    }

  // Compute knots to insert
  // First compute union knot vector
  vector<double> union_knots;
  GeometryTools::makeUnionKnots(knots, knot_tol, union_knots);

  // Extract the knots existing in the union vector, but not in the two knot vectors
  int order = surf1->degree(dir1) + 1;
  set<double> all_new_knots1;
  set<double> all_new_knots2;
  for (int ki=0; ki<num_knot_rows; ++ki)
    {
      vector<double> tmp1;
      extractMissingKnots(union_knots, knots[ki], knot_tol, order, tmp1);
      all_new_knots1.insert(tmp1.begin(), tmp1.end());
      vector<double> tmp2;
      extractMissingKnots(union_knots, knots[num_knot_rows+ki], knot_tol, 
			  order, tmp2);
      all_new_knots2.insert(tmp2.begin(), tmp2.end());
    }
  
  vector<double> new_knots1(all_new_knots1.begin(), all_new_knots1.end());
  vector<double> new_knots2(all_new_knots2.begin(), all_new_knots2.end());
  
  // If knot vector 2 has been altered, represent the corresponding missing
  // knots in the original interval

  // Define end parameters of knot intervals to insert. Set up refinement info
  vector<LRSplineSurface::Refinement2D> refs1;
  defineRefinements(m1, dir1, edge1, ix1, new_knots1, element_width, refs1);

  vector<LRSplineSurface::Refinement2D> refs2;
  defineRefinements(m2, dir2, edge2, ix2, new_knots2, element_width, refs2);

  // // To ensure equally sized corresponding B-spline domains along the boundary
  // // curve after refinement, let the element_width knot lines closest to 
  // // the boundary traverse the entire domain
  // for (int ki=1; ki<=element_width; ++ki)
  //   {
  //     double par1 = m1.kval(flip(dir1), (ix1==0) ? ix1+ki : ix1-ki);
  //     vector<pair<int, int> > zero_segm1 =
  // 	m1.zeroSegments(flip(dir1), (ix1==0) ? ix1+ki : ix1-ki);
  //     for (size_t kj=0; kj<zero_segm1.size(); ++kj)
  // 	{
  // 	  LRSplineSurface::Refinement2D curr_ref;
  // 	  curr_ref.setVal(par1, m1.kval(dir1, zero_segm1[kj].first),
  // 			  m1.kval(dir1, zero_segm1[kj].second), 
  // 			  flip(dir1), 1);
  // 	  refs1.push_back(curr_ref);
  // 	}

  //     double par2 = m2.kval(flip(dir2), (ix2==0) ? ix2+ki : ix2-ki);
  //     vector<pair<int, int> > zero_segm2 =
  // 	m2.zeroSegments(flip(dir2), (ix2==0) ? ix2+ki : ix2-ki);
  //     for (size_t kj=0; kj<zero_segm2.size(); ++kj)
  // 	{
  // 	  LRSplineSurface::Refinement2D curr_ref;
  // 	  curr_ref.setVal(par2, m2.kval(dir2, zero_segm2[kj].first),
  // 			  m2.kval(dir2, zero_segm2[kj].second), 
  // 			  flip(dir2), 1);
  // 	  refs2.push_back(curr_ref);
  // 	}
  //   }
  // Refine surfaces
 #ifdef DEBUG
  std::ofstream ofmesh1("mesh1.eps");
  writePostscriptMesh(*surf1, ofmesh1);
  std::ofstream ofmesh2("mesh2.eps");
  writePostscriptMesh(*surf2, ofmesh2);
#endif
 
  // For both refinements we make sure that we end up with at most 1 lines (not multiplicities).
  for (size_t kj=0; kj<refs1.size(); ++kj)
      surf1->refine(refs1[kj], true);

  for (size_t kj=0; kj<refs2.size(); ++kj)
      surf2->refine(refs2[kj], true);

 #ifdef DEBUG
  std::ofstream ofmesh3("mesh1_2.eps");
  writePostscriptMesh(*surf1, ofmesh3);
  std::ofstream ofmesh4("mesh2_2.eps");
  writePostscriptMesh(*surf2, ofmesh4);
#endif

  return true;
}

//==============================================================================
bool LRSurfStitch::averageEdge(shared_ptr<LRSplineSurface> surf1, int edge1,
			       shared_ptr<LRSplineSurface> surf2, int edge2, 
			       int cont, double tol)
//==============================================================================
{
  int dim = surf1->dimension();
  // @@sbr201504 In matchSplineSpace() 'cont + 2' is used as input ...
  // @@sbr 201504 It may make sense to require tensor-product-structure along common boundaries to avoid
  // problems with over-determined basis functions (gamma != 1.0). But no need to use c1-methods for c0 in this function.
  int element_width = cont+1;//std::max(cont+1, 2); // Yes, even for c0 we require 2 elements along border.
  if (surf2->dimension() != dim)
    return false;
  if (dim != 1)
    element_width = 1;   // C1 stitching only implemented for 1D surfaces

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

  // Traverse the boundary B-splines of the two surfaces and compute the
  // average coefficients
  // Assumes now that the surfaces are oriented equally
  
  // Fetch B-splines along edges of surfaces
  vector<vector<LRBSpline2D*> > bsplines1(element_width);
  extractBoundaryBsplines(surf1, edge1, bsplines1);

  vector<vector<LRBSpline2D*> > bsplines2(element_width);
  extractBoundaryBsplines(surf2, edge2, bsplines2);

  // Traverse B-splines and compute average coefficients
  
  // Sort B-splines to ensure corresponding sequences
  Direction2D dir1 = (edge1 == 0 || edge1 == 1) ? YFIXED : XFIXED;
  Direction2D dir2 = (edge2 == 0 || edge2 == 1) ? YFIXED : XFIXED;
  for (int ki=0; ki<element_width; ++ki)
    {
      if (bsplines1[ki].size() != bsplines2[ki].size())
      {
#ifndef NDEBUG
	  {
	      std::ofstream ofmesh("mesh1.eps");
	      writePostscriptMesh(*surf1, ofmesh);
	      std::ofstream ofmesh2("mesh2.eps");
	      writePostscriptMesh(*surf2, ofmesh2);
	  }
#endif
	return false;
      }
      std::sort(bsplines1[ki].begin(), bsplines1[ki].end(), 
		(dir1==XFIXED) ? bspline_sort_in_u : bspline_sort_in_v);
      std::sort(bsplines2[ki].begin(), bsplines2[ki].end(), 
		(dir2==XFIXED) ? bspline_sort_in_u : bspline_sort_in_v);
    }

#ifdef DEBUG
  std::ofstream of("bspline.txt");
#endif
  for (size_t kj=element_width; kj<bsplines1[0].size()-element_width; ++kj)
    { // The corner coefs should be handled separately.
      Point coef1 = bsplines1[0][kj]->Coef();
      double gamma1 = bsplines1[0][kj]->gamma();
      Point coef2 = bsplines2[0][kj]->Coef();
      double gamma2 = bsplines2[0][kj]->gamma();
#ifdef DEBUG
      of << kj << std::endl;
      of << " bspline1: " << bsplines1[0][kj]->umin() << "  " << bsplines1[0][kj]->umax();
      of << " " << bsplines1[0][kj]->vmin() << "  " << bsplines1[0][kj]->vmax();
      of << " " << coef1[0] << " " << gamma1 << std::endl;
      of << " bspline2: "  << bsplines2[0][kj]->umin() << "  " << bsplines2[0][kj]->umax();
      of << " " << bsplines2[0][kj]->vmin() << "  " << bsplines2[0][kj]->vmax();
      of << " " << coef2[0] << " " << gamma2 << std::endl << std::endl;
#endif
      if (element_width == 1)
	{
	  Point coef = 0.5*(coef1 + coef2);
	  bsplines1[0][kj]->setCoefAndGamma(coef, gamma1);
	  bsplines2[0][kj]->setCoefAndGamma(coef, gamma2);
	}
      else // element_width == 2
	{
	  // @@@ VSK, 150313. Must be modified if dim > 1
	  // C1 continuity is satisfied if all 4 coefficients lies
	  // along a line. 
	  // Compute least squares linear polynomial: f(t) = at + b
	  LRBSpline2D* bsp[4];
	  bsp[0] = bsplines1[1][kj];
	  bsp[1] = bsplines1[0][kj];
	  bsp[2] = bsplines2[0][kj];
	  bsp[3] = bsplines2[1][kj];
	  makeLineC1(bsp, flip(dir1));
	}
    }

  return true;
}


//==============================================================================
void LRSurfStitch::makeLineC1(LRBSpline2D* bsp[4], Direction2D dir)
//==============================================================================
{
  double par1 = bsp[0]->getGrevilleParameter(dir);
  double par2 = bsp[3]->getGrevilleParameter(dir);

  // C1 continuity is satisfied if all 4 coefficients lies
  // along a line. 
  // Compute least squares linear polynomial: f(t) = at + b
  double ls1=0, ls23=0, hs1=0, hs2=0;
  double ls4 = 1.0;
  int ki;
  double par[4];
  Point coef[4];
  for (ki=0; ki<4; ++ki)
    {
      double t = bsp[ki]->getGrevilleParameter(dir);
      //t = (t - par1)/(par2 - par1);  // Reparameterize to [0,1]
      par[ki] = t;// - par1;
      coef[ki] = bsp[ki]->Coef();
    }

  // par[1] = par[2] = 0.5*(par[1]+par[2]);
  // coef[1][0] = coef[2][0] = 0.5*(coef[1][0]+coef[2][0]);
  for (ki=0; ki<4; ++ki)
    {
      double t = par[ki];
      double z = coef[ki][0];
      ls1 += (t*t);
      ls23 += t;
      hs1 += (t*z);
      hs2 += z;
    }
  double det = ls1*ls4 - ls23*ls23;
  double a = (ls4*hs1 - ls23*hs2)/det;
  double b = (ls1*hs2 - ls23*hs1)/det;
  
  // Replace coefficients
  // for (ki=0; ki<4; ++ki)
  //   {
  //     coef[ki][0] = a*par[ki] + b;
  //     double gamma = bsp[ki]->gamma();
  //     bsp[ki]->setCoefAndGamma(coef[ki], gamma);
  //   }
  //Point coefn = ((par[1]-par[0])*coef[3] + (par[3]-par[2])*coef[0])/(par[3]-par[0]);
  par[0] = (dir == XFIXED) ? bsp[0]->umin() : bsp[0]->vmin();
  par[1] = (dir == XFIXED) ? bsp[1]->umax() : bsp[1]->vmax();
  par[2] = (dir == XFIXED) ? bsp[2]->umin() : bsp[2]->vmin();
  par[3] = (dir == XFIXED) ? bsp[3]->umax() : bsp[3]->vmax();
  Point coefn = ((par[3]-par[2])*coef[0] + (par[1]-par[0])*coef[3])/(par[3]-par[0]);
  //Point coefn = 0.5*(coef[3] + coef[0]);
#ifdef DEBUG
  if (fabs(bsp[0]->gamma()-1.0)>1.0e-10)
    {
      std::cout << "Bspline0: (" << bsp[0]->umin() <<"," << bsp[0]->umax();
      std::cout << "," << bsp[0]->vmin() << "," << bsp[0]->vmax();
      std::cout << "), gamma = " << bsp[0]->gamma() << std::endl;
    }
  if (fabs(bsp[1]->gamma()-1.0)>1.0e-10)
    {
      std::cout << "Bspline1: (" << bsp[1]->umin() <<"," << bsp[1]->umax();
      std::cout << "," << bsp[1]->vmin() << "," << bsp[1]->vmax();
      std::cout << "), gamma = " << bsp[1]->gamma() << std::endl;
    }
  if (fabs(bsp[2]->gamma()-1.0)>1.0e-10)
    {
      std::cout << "Bspline2: (" << bsp[2]->umin() <<"," << bsp[2]->umax();
      std::cout << "," << bsp[2]->vmin() << "," << bsp[2]->vmax();
      std::cout << "), gamma = " << bsp[2]->gamma() << std::endl;
    }
  if (fabs(bsp[3]->gamma()-1.0)>1.0e-10)
    {
      std::cout << "Bspline3: (" << bsp[3]->umin() <<"," << bsp[3]->umax();
      std::cout << "," << bsp[3]->vmin() << "," << bsp[3]->vmax();
      std::cout << "), gamma = " << bsp[3]->gamma() << std::endl;
    }
#endif
  bsp[1]->setCoefAndGamma(coefn, 1.0);
  bsp[2]->setCoefAndGamma(coefn, 1.0);
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
void LRSurfStitch::extractBoundaryBsplines(shared_ptr<LRSplineSurface> surf,
					   int edge, 
					   vector<vector<LRBSpline2D*> >& bsplines)
//==============================================================================
{
  // edge is given as: 0=left, 1=right, 2=lower, 3=upper

  // Traverse B-splines and extract those that have sufficient multiplicity
  // at the specified edge. 
  int element_width = (int)bsplines.size();
  // Extract the number of B-splines in the opposite direction specified
  // by element width
  // First define one-sided key values
  int ki;
  double eps = 1.0e-12;
  bool at_start = (edge == 0 || edge == 2);
  Direction2D dir = (edge == 0 || edge == 1) ? XFIXED : YFIXED;
  int deg = surf->degree(dir);
  const Mesh2D& mesh = surf->mesh();
  int ix = (at_start) ? 0 : mesh.numDistinctKnots(dir) - 1;
  vector<double> lim(2*element_width);
  vector<int> mm(element_width);
  for (ki=0; ki<element_width; ++ki)
    {
      mm[ki] = deg + 1 - ki;
      lim[2*ki] = mesh.kval(dir, at_start ? ix : ix-ki-1);
      lim[2*ki+1] = mesh.kval(dir, at_start ? ix+ki+1 : ix);
    }

  for (ki=0; ki<element_width; ++ki)
    bsplines[ki].clear();
  for (LRSplineSurface::BSplineMap::iterator it=surf->basisFunctionsBeginNonconst(); 
       it != surf->basisFunctionsEndNonconst(); ++it)
	{
	  int mult = (edge <= 1) ? it->second->endmult_u(at_start) :
	    it->second->endmult_v(at_start);
	  double t1 = mesh.kval(dir, it->second->suppMin(dir));
	  double t2 = mesh.kval(dir, it->second->suppMax(dir));

	  for (ki=0; ki<element_width; ++ki)
	    {
	      if (mult == mm[ki] && fabs(t1-lim[2*ki]) < eps &&
		  fabs(t2-lim[2*ki+1]) < eps)
		bsplines[ki].push_back(it->second.get());
	    }
	}
}
  

void LRSurfStitch::surfaceDifference(const LRSplineSurface* sf1, const Point& param1,
				     const LRSplineSurface* sf2, const Point& param2,
//				     bool along_u_dir,
				     double& sf_dist, double& tang1_ang, double& tang2_ang)
{
    const int derivs = 1;
    const int totpts = (derivs + 1)*(derivs + 2)/2;
    vector<Point> sf_pt1(totpts), sf_pt2(totpts);
    sf1->point(sf_pt1, param1[0], param1[1], derivs);
    sf2->point(sf_pt2, param2[0], param2[1], derivs);

    sf_dist = sf_pt1[0].dist(sf_pt2[0]);
    tang1_ang = sf_pt1[1].angle(sf_pt2[1]);
    tang2_ang = sf_pt1[2].angle(sf_pt2[2]);

    // tangent_ang = (along_u_dir) ? sf_pt1[1].angle(sf_pt2[1]) : sf_pt1[2].angle(sf_pt2[2]);

    // Point normal1, normal2;
    // try
    // {
    // 	sf1->normal(normal1, param1[0], param1[1]);
    // 	sf2->normal(normal2, param2[0], param2[1]);
    // }
    // catch (...)
    // {
    // 	std::cout << "tangent dim: " << sf_pt1[1].dimension() << std::endl;
    // }
    // normal_ang = normal1.angle(normal2);
}
