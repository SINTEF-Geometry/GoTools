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


#include "GoTools/lrsplines3D/LRVolStitch.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Mesh3DUtils.h"
#include "GoTools/lrsplines3D/LRBSpline3DUtils.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/Utils.h"
//#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <iostream> // @@ debug
#include <fstream> // @@ debug

//#define DEBUG
#define DEBUG2

using namespace Go;
using std::vector;
using std::set;
using std::pair;



//==============================================================================
void LRVolStitch::stitchRegVols(vector<shared_ptr<LRSplineVolume> >& vols,
				int nmb_u, int nmb_v, int nmb_w, double eps,
				int cont)
//==============================================================================
{
#ifdef DEBUG2
  vector<double> bd_dist1 = analyzeContinuity(vols, nmb_u, nmb_v, nmb_w, cont);
  std::cout << "Continuity figures initial: ";
  for (size_t kb=0; kb<bd_dist1.size(); ++kb)
    std::cout << bd_dist1[kb] << " ";
  std::cout << std::endl;
#endif
  
  // Ensure corresponding spline spaces along common boundaries
  // We also make sure that the surfaces are full tensor product volumes
  // along adjacent edges (the first couple of element rows).
  // Note that the volumes (i.e. the coefs) are not altered.
  consistentSplineSpaces(vols, nmb_u, nmb_v, nmb_w, eps, cont);
  consistentSplineSpaces(vols, nmb_u, nmb_v, nmb_w, eps, cont);

  // Stitch volumes along common edges (by altering the coefs).
  // Assosiated corners will be handled first
  // A corner is associated three indices, one in each parameter direction: 0=min, 1=max
  // An edge is described by the running direction (Direction3D) and indices for the remaining
  // parameter directions: 0=min, 1=max. The sequence of directions is always equal so if the
  // running direction of the edge is YDIR, then the first remaining index corresponds to Z and
  // the second to X.

  int kj, kr, kh;
  int nmb_modified;
  bool matched = false;
  for (kh=0; kh<=nmb_w; ++kh)
    {
      for (kj=0; kj<=nmb_v; ++kj)
	{
	  for (kr=0; kr<=nmb_u; ++kr)
	    {
	      vector<VolCorner> corner_match;
	      if (kh > 0)
		{
		  if (kj > 0)
		    {
		      if (kr > 0)
			corner_match.push_back(VolCorner(vols[((kh-1)*nmb_v+(kj-1))*nmb_u+kr-1], 1, 1, 1));
		      if (kr < nmb_u)
			corner_match.push_back(VolCorner(vols[((kh-1)*nmb_v+(kj-1))*nmb_u+kr], 0, 1, 1));
		    }
		  if (kj < nmb_v)
		    {
		      if (kr > 0)
			corner_match.push_back(VolCorner(vols[((kh-1)*nmb_v+kj)*nmb_u+kr-1], 1, 0, 1));
		      if (kr < nmb_u)
			corner_match.push_back(VolCorner(vols[((kh-1)*nmb_v+kj)*nmb_u+kr], 0, 0, 1));
		    }
		}
	      if (kh < nmb_w)
		{
		  if (kj > 0)
		    {
		      if (kr > 0)
			corner_match.push_back(VolCorner(vols[(kh*nmb_v+(kj-1))*nmb_u+kr-1], 1, 1, 0));
		      if (kr < nmb_u)
			corner_match.push_back(VolCorner(vols[(kh*nmb_v+(kj-1))*nmb_u+kr], 0, 1, 0));
		    }
		  if (kj < nmb_v)
		    {
		      if (kr > 0)
			corner_match.push_back(VolCorner(vols[(kh*nmb_v+kj)*nmb_u+kr-1], 1, 0, 0));
		      if (kr < nmb_u)
			corner_match.push_back(VolCorner(vols[(kh*nmb_v+kj)*nmb_u+kr], 0, 0, 0));
		    }
		}
	      if (corner_match.size() > 1)
		{
		  if (cont == 0)
		    nmb_modified = averageCorner(corner_match, eps);
		  else
		    nmb_modified = makeCornerC1(corner_match, eps);
		}
	    }
	}
    }

  for (kh=0; kh<=nmb_w; ++kh)
    {
      for (kj=0; kj<=nmb_v; ++kj)
	{
	  for (kr=0; kr<nmb_u; ++kr)
	    {
	      // Running direction is X
	      vector<VolEdge> edge_match;
	      if (kh > 0)
		{
		  if (kj > 0)
		    edge_match.push_back(VolEdge(vols[((kh-1)*nmb_v+(kj-1))*nmb_u+kr], XDIR, 1, 1));
		  if (kj < nmb_v)
		    edge_match.push_back(VolEdge(vols[((kh-1)*nmb_v+kj)*nmb_u+kr], XDIR, 0, 1));
		}

	      if (kh < nmb_w)
		{
		  if (kj > 0)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+(kj-1))*nmb_u+kr], XDIR, 1, 0));
		  if (kj < nmb_v)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+kj)*nmb_u+kr], XDIR, 0, 0));
		}
	      if (edge_match.size() > 1)
		matched = averageEdge(edge_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed XDIR match! " << std::endl;
#endif
		}
	    }
	}
      
      for (kr=0; kr<=nmb_u; ++kr)
	{
	  for (kj=0; kj<nmb_v; ++kj)
	    {
	      // Running direction is Y
	      vector<VolEdge> edge_match;
	      if (kh > 0)
		{
		  if (kr > 0)
		    edge_match.push_back(VolEdge(vols[((kh-1)*nmb_v+kj)*nmb_u+kr-1], YDIR, 1, 1));
		  if (kr < nmb_u)
		    edge_match.push_back(VolEdge(vols[((kh-1)*nmb_v+kj)*nmb_u+kr], YDIR, 1, 0));
		}

	      if (kh < nmb_w)
		{
		  if (kr > 0)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+kj)*nmb_u+kr-1], YDIR, 0, 1));
		  if (kr < nmb_u)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+kj)*nmb_u+kr], YDIR, 0, 0));
		}
	      if (edge_match.size() > 1)
		matched = averageEdge(edge_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed YDIR match! " << std::endl;
#endif
		}
	    }

	}
    }
  
  for (kj=0; kj<=nmb_v; ++kj)
    {
      for (kr=0; kr<=nmb_u; ++kr)
	{
	  for (kh=0; kh<nmb_w; ++kh)
	    {
	      // Running direction is Z
	      vector<VolEdge> edge_match;
	      if (kj > 0)
		{
		  if (kr > 0)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+(kj-1))*nmb_u+kr-1], ZDIR, 1, 1));
		  if (kr < nmb_u)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+(kj-1))*nmb_u+kr], ZDIR, 0, 1));
		}

	      if (kj < nmb_v)
		{
		  if (kr > 0)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+kj)*nmb_u+kr-1], ZDIR, 1, 0));
		  if (kr < nmb_u)
		    edge_match.push_back(VolEdge(vols[(kh*nmb_v+kj)*nmb_u+kr], ZDIR, 0, 0));
		}
	      if (edge_match.size() > 1)
		matched = averageEdge(edge_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed ZDIR match! " << std::endl;
#endif
		}
	    }
	}
    }
  
  // Match boundary surfaces
  for (kh=0; kh<nmb_w; ++kh)
    {
      for (kj=0; kj<nmb_v; ++kj)
	{
	  for (kr=1; kr<nmb_u; ++kr)
	    {
	      // Running direction is X
	      vector<VolBdSf> bd_match;
	      bd_match.push_back(VolBdSf(vols[(kh*nmb_v+kj)*nmb_u+kr-1], XDIR, 1)); 
	      bd_match.push_back(VolBdSf(vols[(kh*nmb_v+kj)*nmb_u+kr], XDIR, 0)); 
	      matched = averageBdSf(bd_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed XDIR match surface! " << std::endl;
#endif
		}
	    }
	}
      for (kr=0; kr<nmb_u; ++kr)
	{
	  for (kj=1; kj<nmb_v; ++kj)
	    {
	      // Running direction is Y
	      vector<VolBdSf> bd_match;
	      bd_match.push_back(VolBdSf(vols[(kh*nmb_v+kj-1)*nmb_u+kr], YDIR, 1)); 
	      bd_match.push_back(VolBdSf(vols[(kh*nmb_v+kj)*nmb_u+kr], YDIR, 0)); 
	      matched = averageBdSf(bd_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed YDIR match surface! " << std::endl;
#endif
		}
	    }
	}
    }
  for (kj=0; kj<nmb_v; ++kj)
    {
      for (kr=0; kr<nmb_u; ++kr)
	{
	  for (kh=1; kh<nmb_w; ++kh)
	    {
	      // Running direction is Z
	      vector<VolBdSf> bd_match;
	      bd_match.push_back(VolBdSf(vols[((kh-1)*nmb_v+kj)*nmb_u+kr], ZDIR, 1)); 
	      bd_match.push_back(VolBdSf(vols[(kh*nmb_v+kj)*nmb_u+kr], ZDIR, 0)); 
	      matched = averageBdSf(bd_match, cont, eps);
	      if (!matched)
		{
#ifdef DEBUG
		  std::cout << "Failed ZDIR match surface! " << std::endl;
#endif
		}
	    }
	}
    }
#ifdef DEBUG2
  vector<double> bd_dist2 = analyzeContinuity(vols, nmb_u, nmb_v, nmb_w, cont);
  std::cout << "Continuity figures final: ";
  for (size_t kb=0; kb<bd_dist2.size(); ++kb)
    std::cout << bd_dist2[kb] << " ";
  std::cout << std::endl;
#endif
  int stop_break = 1;

}

//==============================================================================
vector<double> LRVolStitch::analyzeContinuity(vector<shared_ptr<LRSplineVolume> >& vols,
					      int nmb_u, int nmb_v, int nmb_w,
					      int derivs, int nmb_samples)
//==============================================================================
{
  assert(derivs < 2);

  int totder = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
  vector<double> global_max_cont_dist(derivs+1, -1.0);

  // We first analyze the cont in the nmb_u direction.
  for (int kr=0; kr<nmb_w; ++kr)
    {
      for (int kj = 0; kj < nmb_v; ++kj)
	{
	  for (int ki = 0; ki < nmb_u - 1; ++ki)
	    {
	      vector<double> max_cont_dist(derivs+1, -1.0);
	      vector<pair<int,int> > sample_ind(derivs+1, std::make_pair(-1,-1));
	      
	      LRSplineVolume* vol_left = vols[(kr*nmb_v+kj)*nmb_u + ki].get();
	      LRSplineVolume* vol_right = vols[(kr*nmb_v+kj)*nmb_u + ki + 1].get();
	      
	      if (vol_left == NULL || vol_right == NULL)
		continue;
	      
	      // We may assume the the vols are constently oriented
	      double upar_left = vol_left->endparam_u();
	      double upar_right = vol_right->startparam_u();
	      
	      // We may also assume that they share domain along the shared edge.
	      double vmin = vol_left->startparam_v();
	      double vmax = vol_left->endparam_v();
	      double vstep = (vmax - vmin)/(nmb_samples - 1);
	      double wmin = vol_left->startparam_w();
	      double wmax = vol_left->endparam_w();
	      double wstep = (wmax - wmin)/(nmb_samples - 1);
		
	      vector<Point> pt_left(totder);
	      vector<Point> pt_right(totder);
	      for (int kh=0; kh<nmb_samples; ++kh)
		{
		  double wpar = wmin + kh*wstep;
		  for (int kk = 0; kk < nmb_samples; ++kk)
		    {
		      double vpar = vmin + kk*vstep;
		      vol_left->point(pt_left, upar_left, vpar, wpar, derivs);
		      vol_right->point(pt_right, upar_right, vpar, wpar, derivs);
		      for (int ka=0; ka<=derivs; ++ka)
			{
			  int ind = (ka == 0) ? 0 : 1;
			  double dist = pt_left[ind].dist(pt_right[ind]);
			  if (dist > max_cont_dist[ind])
			    {
			      max_cont_dist[ind] = dist;
			      sample_ind[ind] = std::make_pair(kk, kh);
			    }
			}
		    }
		}
				     

	      for (int kk = 0; kk < (int)max_cont_dist.size(); ++kk)
		{
		  std::cout << "ki: " << ki << ", kj: " << kr << ", kr: " << kj << ", u-dir, sample #: ";
		  std::cout << sample_ind[kk].first << ", " << sample_ind[kk].second <<
		       ", max_cont_dist[" << kk << "]: " << max_cont_dist[kk] << std::endl;
		  std::cout << "par1: " << upar_left << ", " << vmin+sample_ind[kk].first*vstep;
		  std::cout << ", " << wmin+sample_ind[kk].second*wstep << std::endl;
		  std::cout << "par2: "  << upar_left << ", " << vmin+sample_ind[kk].first*vstep;
		  std::cout << ", " << wmin+sample_ind[kk].second*wstep << std::endl;
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
	      vector<double> max_cont_dist(derivs+1, -1.0);
	      vector<pair<int,int> > sample_ind(derivs+1, std::make_pair(-1,-1));

	      LRSplineVolume* vol_bottom = vols[(kr*nmb_v+ki)*nmb_u + kj].get();
	      LRSplineVolume* vol_top = vols[(kr*nmb_v+ ki + 1)*nmb_u + kj].get();

	      if (vol_bottom == NULL || vol_top == NULL)
		continue;
	      
	      // We may assume the the vols are oriented such that vol_bottom->vmax matches vol_top->vmin.
	      double vpar_bottom = vol_bottom->endparam_v();
	      double vpar_top = vol_top->startparam_v();

	      // We may also assume that they share domain along the shared edge.
	      double umin = vol_bottom->startparam_u();
	      double umax = vol_bottom->endparam_u();
	      double ustep = (umax - umin)/(nmb_samples - 1);
	      double wmin = vol_bottom->startparam_w();
	      double wmax = vol_bottom->endparam_w();
	      double wstep = (wmax - wmin)/(nmb_samples - 1);

	      vector<Point> pt_bottom(totder);
	      vector<Point> pt_top(totder);
	      for (int kh=0; kh<nmb_samples; ++kh)
		{
		  double wpar = wmin + kh*wstep;
		  for (int kk = 0; kk < nmb_samples; ++kk)
		    {
		      double upar = umin + kk*ustep;
		      vol_bottom->point(pt_bottom, upar, vpar_bottom, wpar, derivs);
		      vol_top->point(pt_top, upar, vpar_top, wpar, derivs);
		      for (int ka=0; ka<=derivs; ++ka)
			{
			  int ind = (ka == 0) ? 0 : 2;
			  double dist = pt_top[ind].dist(pt_bottom[ind]);
			  if (dist > max_cont_dist[ka])
			    {
			      max_cont_dist[ka] = dist;
			      sample_ind[ka] = std::make_pair(kk, kh);
			    }
			}
		    }
		}

	      for (int kk = 0; kk < (int)max_cont_dist.size(); ++kk)
		{
		  std::cout << "ki: " << ki << ", kj: " << kr << ", kr: " << kj << ", v-dir, sample #: ";
		  std::cout << sample_ind[kk].first << ", " << sample_ind[kk].second <<
		       ", max_cont_dist[" << kk << "]: " << max_cont_dist[kk] << std::endl;
		  std::cout << "par1: " << umin+sample_ind[kk].first*ustep << ", " << vpar_bottom;
		  std::cout << ", " << wmin+sample_ind[kk].second*wstep << std::endl;
		  std::cout << "par2: " << umin+sample_ind[kk].first*ustep << ", "  << vpar_top;
		  std::cout << ", " << wmin+sample_ind[kk].second*wstep << std::endl;
		  if (max_cont_dist[kk] > global_max_cont_dist[kk])
		    {
		      global_max_cont_dist[kk] = max_cont_dist[kk];
		    }
		}
	    }
	}
    }
  
  // We then analyze the cont in the nmb_w direction.
  for (int kj = 0; kj < nmb_v; ++kj)
    {
      for (int ki = 0; ki < nmb_u; ++ki)
	{
	  for (int kr=0; kr<nmb_w-1; ++kr)
	    {
	      vector<double> max_cont_dist(derivs+1, -1.0);
	      vector<pair<int,int> > sample_ind(derivs+1, std::make_pair(-1,-1));
	      
	      LRSplineVolume* vol_front = vols[(kr*nmb_v+kj)*nmb_u + ki].get();
	      LRSplineVolume* vol_back = vols[((kr+1)*nmb_v+kj)*nmb_u + ki].get();
	      
	      if (vol_front == NULL || vol_back == NULL)
		continue;
	      
	      // We may assume the the volumes are constently oriented
	      double wpar_front = vol_front->endparam_w();
	      double wpar_back = vol_back->startparam_w();
	      
	      // We may also assume that they share domain along the shared edge.
	      double vmin = vol_front->startparam_v();
	      double vmax = vol_front->endparam_v();
	      double vstep = (vmax - vmin)/(nmb_samples - 1);
	      double umin = vol_front->startparam_u();
	      double umax = vol_front->endparam_u();
	      double ustep = (umax - umin)/(nmb_samples - 1);
		
	      vector<Point> pt_front(totder);
	      vector<Point> pt_back(totder);
	      for (int kh=0; kh<nmb_samples; ++kh)
		{
		  double vpar = vmin + kh*vstep;
		  for (int kk = 0; kk < nmb_samples; ++kk)
		    {
		      double upar = umin + kk*ustep;
		      vol_front->point(pt_front, upar, vpar, wpar_front, derivs);
		      vol_back->point(pt_back, upar, vpar, wpar_back, derivs);
		      for (int ka=0; ka<=derivs; ++ka)
			{
			  int ind = (ka == 0) ? 0 : 3;
			  double dist = pt_front[ind].dist(pt_back[ind]);
			  if (dist > max_cont_dist[ka])
			    {
			      max_cont_dist[ka] = dist;
			      sample_ind[ka] = std::make_pair(kk, kh);
			    }
			}
		    }
		}
				     

	      for (int kk = 0; kk < (int)max_cont_dist.size(); ++kk)
		{
		  std::cout << "ki: " << ki << ", kj: " << kr << ", kr: " << kj << ", w-dir, sample #: ";
		  std::cout << sample_ind[kk].first << ", " << sample_ind[kk].second <<
		       ", max_cont_dist[" << kk << "]: " << max_cont_dist[kk] << std::endl;
		  std::cout << "par1: " << umin+sample_ind[kk].first*ustep << ", ";
		  std::cout << vmin+sample_ind[kk].second*vstep<< ", " << wpar_front << std::endl;
		  std::cout << "par2: " << umin+sample_ind[kk].first*ustep << ", " ;
		  std::cout << vmin+sample_ind[kk].second*vstep << ", " << wpar_back<< std::endl;
		  if (max_cont_dist[kk] > global_max_cont_dist[kk])
		    {
		      global_max_cont_dist[kk] = max_cont_dist[kk];
		    }
		}
	    }
	}
    }
   
   return global_max_cont_dist;
}

//==============================================================================
void LRVolStitch::consistentSplineSpaces(vector<shared_ptr<LRSplineVolume> >& vols,
					 int nmb_u, int nmb_v, int nmb_w, double eps,
					 int cont)
//==============================================================================
{
  int kj, kr, kh;

  int max_nmb = std::max(nmb_u, std::max(nmb_v, nmb_w)) + 1;
  // The current method refines globally around surface corners. In theory this means that 'max_nmb - 1'
  // corner adjustments must be performed for the method to propagate along the nmb_u x nmb_v-grid.  Once
  // there is no change in the number of basis functions we exit.
  int sum_basis_functions = 0;
  for (int kk = 0; kk < vols.size(); ++kk)
  {
      if (vols[kk].get() != NULL)
      {
	  sum_basis_functions += vols[kk]-> numBasisFunctions();
      }
  }

  for (int ki = 0; ki < max_nmb - 1; ++ki)
    {
      // First ensure tensor-product structure close to the boundaries
      bool bdsfs[6];
      for (kh=0; kh<nmb_w; ++kh)
	{
	  bdsfs[4] = (kh != 0);
	  bdsfs[5] = (kh != nmb_w-1);
	  for (kj=0; kj<nmb_v; ++kj)
	    {
	      bdsfs[2] = (kj != 0);
	      bdsfs[3] = (kj != nmb_v-1);
	      for (kr=0; kr<nmb_u; ++kr)
		{
		  if (!vols[(kh*nmb_v+kj)*nmb_u+kr].get())
		    continue;
		  bdsfs[0] = (kr != 0);
		  bdsfs[1] = (kr != nmb_u-1);
		  // Other parts of the code requires 2 inner rows for c0.
		  int num_inner_rows = std::max(cont + 1, 2); // I.e. the number of elements that are affected.
		  tensorStructure(vols[(kh*nmb_v+kj)*nmb_u+kr], num_inner_rows, bdsfs);
		}
	    }
	}

      // Then make corresponding spline spaces across boundaries (i.e. insert knots along the common edge).
      const int element_width = cont + 2;//std::max(2, cont + 1);//cont + 2; // Number of rows with inner knots.
      for (kh=0; kh<nmb_w; ++kh)
      {
	for (kj=0; kj<nmb_v; ++kj)
	  {
	    for (kr=0; kr<nmb_u; ++kr)
	      {
		// Vertical boundary surface, left to right
		// Boundary surfaces are numbered: 0=left, 1=right, 2=lower, 3=upper, 4=front, 5=back
		if (kr > 0 && vols[(kh*nmb_v+kj)*nmb_u+kr-1].get() && 
		    vols[(kh*nmb_v+kj)*nmb_u+kr].get())
		  matchSplineSpace(vols[(kh*nmb_v+kj)*nmb_u+kr-1], 1,
				   vols[(kh*nmb_v+kj)*nmb_u+kr], 0, element_width, eps);

		// Horizontal boundary surface
		if (kj > 0 && vols[(kh*nmb_v+kj-1)*nmb_u+kr].get() &&
		    vols[(kh*nmb_v+kj)*nmb_u+kr].get())
		  matchSplineSpace(vols[(kh*nmb_v+kj-1)*nmb_u+kr], 3,
				   vols[(kh*nmb_v+kj)*nmb_u+kr], 2, element_width, eps);
		
		// Vertical boundary surface, front to back
		if (kh > 0 && vols[((kh-1)*nmb_v+kj)*nmb_u+kr].get() &&
		    vols[(kh*nmb_v+kj)*nmb_u+kr].get())
		  matchSplineSpace(vols[((kh-1)*nmb_v+kj)*nmb_u+kr], 5,
				   vols[(kh*nmb_v+kj)*nmb_u+kr], 4, element_width, eps);
	      }
	  }
      }

      int new_sum_basis_functions = 0;
      for (int kk = 0; kk < vols.size(); ++kk)
      {
	  if (vols[kk].get() != NULL)
	  {
	      new_sum_basis_functions += vols[kk]-> numBasisFunctions();
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
void LRVolStitch::checkCornerMatch(std::vector<VolCorner>& volc, 
				   double tol, vector<int>& nmb_match,
				   vector<Point>& corner_val, 
				   vector<Point>& param)
//==============================================================================
{
  // Remove the surfaces that are too far from a common corner to be modified
  // For a 1D surface, the distance is computed in the parameter domain of the surface
  // otherwise in geometry space
  size_t kr, kh;
  vector<Point> corner(volc.size());
  for (kr=0; kr<volc.size(); ++kr)
    {
      if (!volc[kr].vol_.get())
	{
	  continue;
	}

      double u = (volc[kr].uix_ == 0) ? volc[kr].vol_->startparam_u() :	volc[kr].vol_->endparam_u();
      double v = (volc[kr].vix_ == 0) ? volc[kr].vol_->startparam_v() :	volc[kr].vol_->endparam_v();
      double w = (volc[kr].wix_ == 0) ? volc[kr].vol_->startparam_w() :	volc[kr].vol_->endparam_w();
      Point pos;
      volc[kr].vol_->point(pos, u, v, w);
      if (volc[kr].vol_->dimension() == 1)
	corner[kr] = Point(u,v,w);
      else
	corner[kr] = pos;
      corner_val[kr] = pos;
      param[kr] = Point(u,v,w);
    }

  for (kr=0; kr<volc.size()-1; ++kr)
    {
      if (!volc[kr].vol_.get())
	continue;
      for (kh=kr+1; kh<volc.size(); ++kh)
	{
	  if (!volc[kh].vol_.get())
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
int LRVolStitch::averageCorner(std::vector<VolCorner>& volc, double tol)
//==============================================================================
{
  // Check if the given corners match
  vector<int> nmb_match(volc.size(), 0);
  vector<Point> corner_val(volc.size());
  vector<Point> param(volc.size());
  checkCornerMatch(volc, tol, nmb_match, corner_val, param);

  // Compute average corner
  size_t kr;
  int dim;
  int nmb_mod = (int)nmb_match.size();
  for (kr=0; kr<nmb_match.size(); ++kr)
    {
      if (nmb_match[kr] == 0)
	nmb_mod--;
      else
	dim = volc[kr].vol_->dimension();
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
   for (kr=0; kr<volc.size(); ++kr)
     {
       if (volc[kr].vol_.get() && nmb_match[kr] > 0)
	 av_corner += corner_val[kr];
     }
   av_corner /= (double)nmb_mod;

  // Replace corner coefficient
   for (size_t kr=0; kr<volc.size(); ++kr)
     {
       if (volc[kr].vol_.get() && nmb_match[kr] > 0)
	 {
	   // Identify B-spline
	   vector<LRBSpline3D*> cand_bsplines = 
	     volc[kr].vol_->basisFunctionsWithSupportAt(param[kr][0], param[kr][1],
							param[kr][2]);
	   size_t kh;
	   for (kh=0; kh<cand_bsplines.size(); ++kh)
	     {
	       int deg1 = cand_bsplines[kh]->degree(XDIR);
	       int deg2 = cand_bsplines[kh]->degree(YDIR);
	       int deg3 = cand_bsplines[kh]->degree(ZDIR);
	       int mult1 = cand_bsplines[kh]->endmult_u(true);
	       int mult2 = cand_bsplines[kh]->endmult_u(false);
	       int mult3 = cand_bsplines[kh]->endmult_v(true);
	       int mult4 = cand_bsplines[kh]->endmult_v(false);
	       int mult5 = cand_bsplines[kh]->endmult_w(true);
	       int mult6 = cand_bsplines[kh]->endmult_w(false);
	       if ((mult1 == deg1+1 || mult2 == deg1+1) &&
		   (mult3 == deg2+1 || mult4 == deg2+1) && 
		   (mult5 == deg3+1 || mult6 == deg3+1))
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
int LRVolStitch::makeCornerC1(vector<VolCorner>& volc, double tol)
//==============================================================================
{
  // Check if the given corners match
  vector<int> nmb_match(volc.size(), 0);
  vector<Point> corner_val(volc.size());
  vector<Point> param(volc.size());
  checkCornerMatch(volc, tol, nmb_match, corner_val, param);

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

  // Extract corner B-splines for all matching volumes
  vector<LRBSpline3D*> bsplines;
  for (size_t ki=0; ki<volc.size(); ++ki)
    {
      if (volc[ki].vol_.get() && nmb_match[ki]>0)
	{
	  int deg1 = volc[ki].vol_->degree(XDIR);
	  int deg2 = volc[ki].vol_->degree(YDIR);
	  int deg3 = volc[ki].vol_->degree(ZDIR);

	   vector<LRBSpline3D*> cand_bsplines = 
	     volc[ki].vol_->basisFunctionsWithSupportAt(param[ki][0], param[ki][1],
							 param[ki][2]);
	   for (kr=0; kr<cand_bsplines.size(); ++kr)
	     {
	       int mult1 = cand_bsplines[kr]->endmult_u(true);
	       int mult2 = cand_bsplines[kr]->endmult_u(false);
	       int mult3 = cand_bsplines[kr]->endmult_v(true);
	       int mult4 = cand_bsplines[kr]->endmult_v(false);
	       int mult5 = cand_bsplines[kr]->endmult_w(true);
	       int mult6 = cand_bsplines[kr]->endmult_w(false);
	       if ((mult1 == deg1+1 || mult1 == deg1 ||
		    mult2 == deg1+1 || mult2 == deg1) &&
		   (mult3 == deg2+1 || mult3 == deg2 ||
		    mult4 == deg2+1 || mult4 == deg2) && 
		   (mult5 == deg3+1 || mult5 == deg3 ||
		    mult6 == deg3+1|| mult6 == deg3))
		 bsplines.push_back(cand_bsplines[kr]);
	     }
	}
    }

  // Find least squares plane for the corner coefficients of all
  // volumes (1D only)
  // Set up equation system
  double ls1=0, ls2_5=0, ls3_10=0, ls4_13=0, ls6=0;
  double ls7_9=0, ls8_14=0, ls11=0, ls12_15=0;
  double hs1=0, hs2=0, hs3=0, hs4=0;
  double ls16 = 1.0;

  int kj;
  vector<Point> par(bsplines.size());
  vector<Point> coef(bsplines.size());
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      par[kj] = bsplines[kj]->getGrevilleParameter();
      double x = par[kj][0];
      double y = par[kj][1];
      double z = par[kj][2];
      coef[kj] = bsplines[kj]->Coef();
      double w = coef[kj][0];
      ls1 += (x*x);
      ls2_5 += (x*y);
      ls3_10 += (x*z);
      ls4_13 += x;
      ls6 += (y*y);
      ls7_9 += (y*z);
      ls8_14 += y;
      ls11 += (z*z);
      ls12_15 += z;
      hs1 += (w*x);
      hs2 += (w*y);
      hs3 += (w*z);
      hs4 += w;
    }

  MatrixXD<double,4> M;
  M(0,0) = ls1;
  M(0,1) = ls2_5;
  M(0,2) = ls3_10;
  M(0,3) = ls4_13;
  M(1,0) = ls2_5;
  M(1,1) = ls6;
  M(1,2) = ls7_9;
  M(1,3) = ls8_14;
  M(2,0) = ls3_10;
  M(2,1) = ls7_9;
  M(2,2) = ls11;
  M(2,3) = ls12_15;
  M(3,0) = ls4_13;
  M(3,1) = ls8_14;
  M(3,2) = ls12_15;
  M(3,3) = ls16;

  // Solve equation system by Cramers rule
  double det = M.det();
  double detk[4];
  for (kj=0; kj<4; ++kj)
    {
      MatrixXD<double,4> M1 = M;
      M1(kj,0) = hs1;
      M1(kj,1) = hs2;
      M1(kj,2) = hs3;
      M1(kj,3) = hs4;
      detk[kj] = M1.det();
    }
  double a = detk[0]/det;
  double b = detk[1]/det;
  double c = detk[2]/det;
  double d = detk[3]/det;

  // Project the volumes onto the found plane
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      coef[kj][0] = a*par[kj][0] + b*par[kj][1] + c*par[kj][2]+d;
      double gamma = bsplines[kj]->gamma();
      bsplines[kj]->setCoefAndGamma(coef[kj], gamma);
    }

#ifdef DEBUG
  {
      // We check if the corner continuity is satisfied.
      for (kj = 0; kj < (int)volc.size() - 1; ++kj)
      {
	  const LRSplineVolume* vol1 = volc[kj].vol_.get();
	  const LRSplineVolume* vol2 = volc[kj+1].vol_.get();
	  double vol_dist = -1.0;
	  double tang1_ang = -1.0;
	  double tang2_ang = -1.0;
	  double tang3_ang = -1.0;
	  Point par1 = param[kj];
	  Point par2 = param[kj+1];
	  try
	  {
	    volumeDifference(vol1, par1, vol2, par2, 
			     vol_dist, tang1_ang, tang2_ang, tang3_ang);
	  }
	  catch (...)
	  {
	      MESSAGE("Fail ...");
	  }
	  std::cout << "vol_dist = " << vol_dist << ", tang1_ang = " << tang1_ang;
	  std::cout << ", tang2_ang = " << tang2_ang << ", tang3_ang = " << tang3_ang << std::endl;
      }

  }
#endif

  return nmb_mod;
}

//==============================================================================
void LRVolStitch::makeCrossC1(vector<LRBSpline3D*>& bsplines, Direction3D dir)
//==============================================================================
{
  // Find least squares plane for the current edge coefficients of all
  // volumes (1D only)
  // Set up equation system
  double ls1=0, ls24=0, ls37=0, ls5=0, ls68=0;
  double hs1=0, hs2=0, hs3=0;
  double ls9 = 1.0;

  int kj;
  int ix1 = (dir == XDIR) ? 1 : ((dir == YDIR) ? 2 : 0);
  int ix2 = (ix1 + 1)%3;
  vector<Point> par(bsplines.size());
  vector<Point> coef(bsplines.size());
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      par[kj] = bsplines[kj]->getGrevilleParameter();
      double x = par[kj][ix1];
      double y = par[kj][ix2];
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

  // Project the volume coefficientss onto the found plane
  for (kj=0; kj<(int)bsplines.size(); ++kj)
    {
      coef[kj][0] = a*par[kj][ix1] + b*par[kj][ix2] + c;
      double gamma = bsplines[kj]->gamma();
      bsplines[kj]->setCoefAndGamma(coef[kj], gamma);
    }
}

//==============================================================================

bool bspline_sort_in_u(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  //std::cout << 1 << " " << bspline1 << " " << bspline2 << std::endl;
  if (bspline1->umin() == bspline2->umin())
    return (bspline1->umax() < bspline2->umax());
  return (bspline1->umin() < bspline2->umin());
}

bool bspline_sort_in_v(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  //std::cout << 2 << " " << bspline1 << " " << bspline2 << std::endl;
  if (bspline1->vmin() == bspline2->vmin())
    return (bspline1->vmax() < bspline2->vmax());
  return (bspline1->vmin() < bspline2->vmin());
}

bool bspline_sort_in_w(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  //std::cout << 3 << " " << bspline1 << " " << bspline2 << std::endl;
  if (bspline1->wmin() == bspline2->wmin())
    return (bspline1->wmax() < bspline2->wmax());
  return (bspline1->wmin() < bspline2->wmin());
}


bool bspline_sort_in_uv(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  if (bspline1->vmin() == bspline2->vmin())
    {
      if (bspline1->umin() == bspline2->umin())
	{
	  if (bspline1->vmax() == bspline2->vmax())
	    return (bspline1->umax() < bspline2->umax());
	  return (bspline1->vmax() < bspline2->vmax());
	}
      return (bspline1->umin() < bspline2->umin());
    }
  return (bspline1->vmin() < bspline2->vmin());
}

bool bspline_sort_in_uw(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  if (bspline1->wmin() == bspline2->wmin())
    {
      if (bspline1->umin() == bspline2->umin())
	{
	  if (bspline1->wmax() == bspline2->wmax())
	    return (bspline1->umax() < bspline2->umax());
	  return (bspline1->wmax() < bspline2->wmax());
	}
      return (bspline1->umin() < bspline2->umin());
    }
  return (bspline1->wmin() < bspline2->wmin());
}

bool bspline_sort_in_vw(LRBSpline3D* bspline1, LRBSpline3D* bspline2)
{
  if (bspline1->wmin() == bspline2->wmin())
    {
      if (bspline1->vmin() == bspline2->vmin())
	{
	  if (bspline1->wmax() == bspline2->wmax())
	    return (bspline1->vmax() < bspline2->vmax());
	  return (bspline1->wmax() < bspline2->wmax());
	}
      return (bspline1->vmin() < bspline2->vmin());
    }
  return (bspline1->wmin() < bspline2->wmin());
}


//==============================================================================
void LRVolStitch::tensorStructure(shared_ptr<LRSplineVolume> vol, 
				   int element_width, bool bdsfs[6])
//==============================================================================
{
  const Mesh3D& mesh = vol->mesh();
  vector<LRSplineVolume::Refinement3D> refs1, refs2;
  for (int ki=0; ki<6; ++ki)
    {
      if (!bdsfs[ki])
	continue;    // Do not imply tensor product structure along this boundary

      Direction3D dir = (ki <= 1) ? XDIR : ((ki <= 3) ? YDIR : ZDIR);
      bool at_start = (ki%2 ==0);
      int nmb = mesh.numDistinctKnots(dir);
      for (int kr0=ki+2; kr0<ki+8; ++kr0)
      	{
      	  int kr = kr0 % 6;
      	  Direction3D dir2 = (kr <= 1) ? XDIR : ((kr <= 3) ? YDIR : ZDIR);
      	  if (dir2 == dir)
      	    continue;

      	  if (bdsfs[kr])
      	    continue; // Tensor structure at common boundary handled by adjacent
      	  // boundary surface

      	  // Insert mesh rectangles close to the volume boundary
      	  bool at_start2 = (kr%2 == 0);
      	  int nmb2 = mesh.numDistinctKnots(dir2);
      	  for (int kj=1; kj<=element_width; ++kj)
      	    {
      	      int ix2 = at_start2 ? kj : nmb2-kj-1;
      	      double par2 = mesh.kval(dir2, ix2);

      	      // Define mesh rectangle
      	      Direction3D dir3 = next(dir2);
      	      Direction3D dir4 = prev(dir2);
      	      double min3, max3, min4, max4;
      	      if (dir3 == dir)
      		{
      		  min3 = mesh.kval(dir, at_start ? 0 : nmb-element_width-2);
      		  max3 = mesh.kval(dir, at_start ? element_width+1 : nmb-1);
      		  min4 = mesh.minParam(dir4);
      		  max4 = mesh.maxParam(dir4);
      		}
      	      else
      		{
      		  min3 = mesh.minParam(dir3);
      		  max3 = mesh.maxParam(dir3);
      		  min4 = mesh.kval(dir, at_start ? 0 : nmb-element_width-2);
      		  max4 = mesh.kval(dir, at_start ? element_width+1 : nmb-1);
      		}
      	      refs1.push_back(LRSplineVolume::Refinement3D(par2, min3, max3,
      							  min4, max4, dir2, 1));
      	    }
      	}
    }
  
  // Refine volume
  vol->refine(refs1, true);

  for (int ki=0; ki<6; ++ki)
    {
      if (!bdsfs[ki])
	continue;    // Do not imply tensor product structure along this boundary

      Direction3D dir = (ki <= 1) ? XDIR : ((ki <= 3) ? YDIR : ZDIR);
      bool at_start = (ki%2 ==0);
      int nmb = mesh.numDistinctKnots(dir);
      for (int kj=1; kj<=element_width; ++kj)
	{
	  int ix = at_start ? kj : nmb-kj-1;
	  double par = mesh.kval(dir, ix);

	  // Define a mesh rectangle covering the entire domain at the specified
	  // parameter. The refinement algorithm will combine this rectangle with
	  // the previously existing ones.
	  Direction3D dir1 = next(dir);
	  Direction3D dir2 = prev(dir);
	  double min1 = mesh.minParam(dir1);
	  double max1 = mesh.maxParam(dir1);
	  double min2 = mesh.minParam(dir2);
	  double max2 = mesh.maxParam(dir2);
	  refs2.push_back(LRSplineVolume::Refinement3D(par, min1, max1, min2, max2, dir, 1));
	}

     }

  // Refine volume
  vol->refine(refs2, true);

  return;
}

//==============================================================================
bool LRVolStitch::matchSplineSpace(shared_ptr<LRSplineVolume> vol1, 
				    int bdsf1,
				    shared_ptr<LRSplineVolume> vol2,
				    int bdsf2, 
				    int element_width, double tol)
//==============================================================================
{
  int dim = vol1->dimension();
  if (vol2->dimension() != dim)
    return false;
  
#ifdef DEBUG
  std::ofstream of1("elemmesh1.g2");
  std::ofstream of2("elemmesh2.g2");
  writeElems(of1, vol1);
  writeElems(of2, vol2);
#endif
  
  // Fetch surface corners and check consistency
  Point corner1[4], corner2[4];
  
  fetchSfCorners(vol1, bdsf1, corner1);
  fetchSfCorners(vol2, bdsf2, corner2);

  for (int ki=0; ki<4; ++ki)
    {
      double dist = corner1[ki].dist(corner2[0]);
      for (int kj=1; kj<4; ++kj)
	{
	  double dist2 = corner1[ki].dist(corner2[kj]);
	  dist = std::min(dist,dist2);
	}
      if (dist > tol)
	return false;
    }

  // Fetch mesh rectangles along the bdsf
  // Edges are numbered: 0=left, 1=right, 2=lower, 3=upper, 4=front, 5=back
  const Mesh3D& m1 = vol1->mesh();
  const Mesh3D& m2 = vol2->mesh();
  Direction3D dir0 = (bdsf1 == 0 || bdsf1 == 1) ? XDIR :
    ((bdsf1 == 2 || bdsf1 == 3) ? YDIR : ZDIR);
  Direction3D dir2_0 = (bdsf2 == 0 || bdsf2 == 1) ? XDIR :
    ((bdsf2 == 2 || bdsf2 == 3) ? YDIR : ZDIR);
  if (dir0 != dir2_0)
    return false;
    
  Direction3D dir1 = next(dir0);
  Direction3D dir2 = prev(dir0);
  int start1 = (bdsf1%2 == 0) ? 0 : m1.numDistinctKnots(dir0)-1;
  int end1 = (start1 == 0) ? element_width : start1 - element_width;
  if (start1 > end1)
    std::swap(start1, end1);
  
  int start2 = (bdsf2%2 == 0) ? 0 : m2.numDistinctKnots(dir0)-1;
  int end2 = (start2 == 0) ? element_width : start2 - element_width;
  if (start2 > end2)
    std::swap(start2, end2);

  // Raise surface degree if necessary
  // This functionality is not implemented for LRSplineVolume. Must be done
  // if this function is to be used in a more general setting
  if (vol1->degree(dir0) != vol2->degree(dir0))
    return false;
  
  vector<pair<vector<GPos2D>,int> > rects1_1 =
    m1.overlapRects(dir1, 0, m1.numDistinctKnots(dir2),
		    start1, end1);
  vector<pair<vector<GPos2D>,int> > rects1_2 =
    m1.overlapRects(dir2, start1, end1,
		    0, m1.numDistinctKnots(dir1));
  vector<pair<vector<GPos2D>,int> > rects2_1 =
    m2.overlapRects(dir1, 0, m2.numDistinctKnots(dir2),
		    start2, end2);
  vector<pair<vector<GPos2D>,int> > rects2_2 =
    m2.overlapRects(dir2, start2, end2,
		    0, m2.numDistinctKnots(dir1));
#ifdef DEBUG
  std::ofstream of3("rectmesh1_1.g2");
  writeMrects(of3, m1, rects1_1, dir1, 0);
  std::ofstream of4("rectmesh1_2.g2");
  writeMrects(of4, m1, rects1_2, dir2, 1);
  std::ofstream of5("rectmesh2_1.g2");
  writeMrects(of5, m2, rects2_1, dir1, 2);
  std::ofstream of6("rectmesh2_2.g2");
  writeMrects(of6, m2, rects2_2, dir2, 3);
#endif
  
  // Adapt mesh rectangles to boundary zones and define difference rectangles
  vector<pair<vector<GPos2D>,int> > rects3_1, rects3_2, rects4_1, rects4_2;
  int num1 = m1.numDistinctKnots(dir1);
  int num2 = m1.numDistinctKnots(dir2);
  int num3 = m2.numDistinctKnots(dir1);
  int num4 = m2.numDistinctKnots(dir2);
  adaptMeshRectangles(rects1_1, 1, start1, end1, 0, num1-1, element_width-1, rects3_1);
  adaptMeshRectangles(rects1_2, 0, start1, end1, 0, num2-1, element_width-1, rects3_2);
  adaptMeshRectangles(rects2_1, 1, start2, end2, 0, num3-1, element_width-1, rects4_1);
  adaptMeshRectangles(rects2_2, 0, start2, end2, 0, num4-1, element_width-1, rects4_2);
 
#ifdef DEBUG
  std::ofstream of7("rectmesh1_1n.g2");
  writeMrects(of7, m1, rects1_1, dir1, 4);
  std::ofstream of8("rectmesh1_2n.g2");
  writeMrects(of8, m1, rects1_2, dir2, 5);
  std::ofstream of9("rectmesh2_1n.g2");
  writeMrects(of9, m2, rects2_1, dir1, 6);
  std::ofstream of10("rectmesh2_2n.g2");
  writeMrects(of10, m2, rects2_2, dir2, 7);
#endif
  
  // Define refinements
  vector<LRSplineVolume::Refinement3D> refs1_1, refs1_2;
  vector<LRSplineVolume::Refinement3D> refs2_1, refs2_2;
  vector<LRSplineVolume::Refinement3D> refs3_1, refs3_2;
  vector<LRSplineVolume::Refinement3D> refs4_1, refs4_2;
  defineRefinements(m1, m2, dir1, rects1_1,
		    start2==0, element_width, 1, refs1_1);
  defineRefinements(m1, m2, dir2, rects1_2,
		    start2==0, element_width, 0, refs1_2);
  defineRefinements(m2, m1, dir1, rects2_1,
		    start1==0, element_width, 1, refs2_1);
  defineRefinements(m2, m1, dir2, rects2_2,
		    start1==0, element_width, 0, refs2_2);
  
  defineRefinements(m1, m1, dir1, rects3_1,
		    start1==0, element_width, 1, refs3_1);
  defineRefinements(m1, m1, dir2, rects3_2,
		    start1==0, element_width, 0, refs3_2);
  defineRefinements(m2, m2, dir1, rects4_1,
		    start2==0, element_width, 1, refs4_1);
  defineRefinements(m2, m2, dir2, rects4_2,
		    start2==0, element_width, 0, refs4_2);

  // For both refinements we make sure that we end up with at most 1 lines (not multiplicities).
  vol2->refine(refs1_1, true);
  vol2->refine(refs1_2, true);
  vol1->refine(refs2_1, true);
  vol1->refine(refs2_2, true);

  vol1->refine(refs3_1, true);
  vol1->refine(refs3_2, true);
  vol2->refine(refs4_1, true);
  vol2->refine(refs4_2, true);

#ifdef DEBUG
  vector<pair<vector<GPos2D>,int> > r1_1 =
    m1.overlapRects(dir1, 0, m1.numDistinctKnots(dir2),
		    start1, end1);

  std::ofstream off11("rmesh1_1.g2");
  writeMrects(off11, m1, r1_1, dir1, 0);
  
   vector<pair<vector<GPos2D>,int> > r1_2 =
    m1.overlapRects(dir2, start1, end1,
		    0, m1.numDistinctKnots(dir1));
  std::ofstream off12("rmesh1_2.g2");
  writeMrects(off12, m1, r1_2, dir2, 1);
  
  vector<pair<vector<GPos2D>,int> > r2_1 =
    m2.overlapRects(dir1, 0, m2.numDistinctKnots(dir2),
		    start2, end2);
  std::ofstream off13("rmesh2_1.g2");
  writeMrects(off13, m2, r2_1, dir1, 2);
  
  vector<pair<vector<GPos2D>,int> > r2_2 =
    m2.overlapRects(dir2, start2, end2,
		    0, m2.numDistinctKnots(dir1));
  std::ofstream off14("rmesh2_2.g2");
  writeMrects(off14, m2, r2_2, dir2, 3);

  vector<pair<vector<GPos2D>,int> > r1_3 =
    m1.overlapRects(dir0, 0, m1.numDistinctKnots(dir1),
		    0, m1.numDistinctKnots(dir2));
  std::ofstream off15("rmesh1_3.g2");
  writeMrects(off15, m1, r1_3, dir0, 4);
  
  vector<pair<vector<GPos2D>,int> > r2_3 =
    m2.overlapRects(dir0, 0, m2.numDistinctKnots(dir1),
		    0, m2.numDistinctKnots(dir2));
  std::ofstream off16("rmesh2_3.g2");
  writeMrects(off16, m2, r2_3, dir0, 5);

  std::ofstream of11("elemmesh1_2.g2");
  std::ofstream of12("elemmesh2_2.g2");
  writeElems(of11, vol1);
  writeElems(of12, vol2);
#endif

  return true;
}

//==============================================================================
bool LRVolStitch::averageEdge(vector<VolEdge>& voledg, int cont, double tol)
//==============================================================================
{
  int dim = voledg[0].vol_->dimension();
  cont = std::min(1, cont);  // A most C1 continuity 
  if (dim != 1)
    cont = 0;   // C1 stitching only implemented for 1D volumes

  // Remove empty edge instances
  for (size_t ki=0; ki<voledg.size(); )
    {
      if (!voledg[ki].vol_.get())
	voledg.erase(voledg.begin()+ki);
      else
	++ki;
    }

  if (voledg.size() <= 1)
    return true; // Nothing to do
  
  // Check corner consistency
  vector<Point> corner(2*voledg.size());

  for (size_t ki=0; ki<voledg.size(); ++ki)
    fetchEdgeCorners(voledg[ki], corner[2*ki], corner[2*ki+1]);

  vector<bool> samedir(voledg.size(), true);
  for (size_t ki=1; ki<voledg.size(); ++ki)
    {
      if (corner[0].dist(corner[2*ki]) > tol && corner[0].dist(corner[2*ki+1]) > tol)
	return false;
      if (corner[1].dist(corner[2*ki]) > tol && corner[1].dist(corner[2*ki+1]) > tol)
	return false;
      if (corner[0].dist(corner[2*ki]) > corner[0].dist(corner[2*ki+1]))
	samedir[ki] = false;
    }

  // Traverse the boundary B-splines of the two surfaces and compute the
  // average coefficients
  // Assumes now that the surfaces are oriented equally
  
  // Fetch B-splines along edges of volumes
  vector<vector<vector<LRBSpline3D*> > > bsplines(voledg.size());
  for (size_t ki=0; ki<voledg.size(); ++ki)
    {
      bsplines[ki].resize(1+3*cont);
      extractBoundaryBsplines(voledg[ki], cont, bsplines[ki]);
    }



  // Traverse B-splines and compute average coefficients
  
  // Sort B-splines to ensure corresponding sequences
  if (bsplines[0].size() == 0)
    return false;
  int nmb_b = bsplines[0][0].size();
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[0].size() != bsplines[ki].size())
	return false;
      for (size_t kj=0; kj<bsplines[0].size(); ++kj)
	{
	  if (nmb_b != bsplines[ki][kj].size())
	    {
	      std::cout << "B-spline mismatch" << std::endl;
	      std::ofstream mis1("mis1.g2");
	      writeBSplines(mis1, bsplines[0][0]);
	      std::ofstream mis2("mis2.g2");
	      writeBSplines(mis2, bsplines[ki][kj]);
	      return false;
	    }

	  Direction3D dir = voledg[ki].dir_;
	  std::sort(bsplines[ki][kj].begin(), bsplines[ki][kj].end(), 
		    (dir==XDIR) ? bspline_sort_in_u :
		    ((dir==YDIR) ? bspline_sort_in_v : bspline_sort_in_w));
	}
    }

// #ifdef DEBUG
//   std::ofstream of("bspline.txt");
// #endif
  for (int kr=cont+1; kr<nmb_b-cont-1; ++kr)
  //for (int kr=1; kr<nmb_b-1; ++kr)
    {
      // The corner coefs should be handled separately.
      if (cont == 0)
	{
	  Point av_coef(dim);
	  av_coef.setValue(0.0);
	  for (size_t ki=0; ki<bsplines.size(); ++ki)
	    {
	      Point coef = bsplines[ki][0][kr]->Coef();
	      av_coef += coef;
	    }
	  av_coef /= (double)bsplines.size();

      
	  for (size_t ki=0; ki<bsplines.size(); ++ki)
	    {
	      double gamma = bsplines[ki][0][kr]->gamma();

	      if (cont == 0)
		{
		  bsplines[ki][0][kr]->setCoefAndGamma(av_coef, gamma);
		}
	    }
	}
      else // cont == 1
	{
	  // // @@@ VSK, 150313. Must be modified if dim > 1
	  // // C1 continuity is satisfied if all coefficients lie
	  // // in a plane
	  vector<LRBSpline3D*> bsp;
	  for (size_t ki=0; ki<bsplines.size(); ++ki)
	    for (size_t kj=0; kj<bsplines[ki].size(); ++kj)
	      bsp.push_back(bsplines[ki][kj][kr]);

	  makeCrossC1(bsp, voledg[0].dir_);

	}
    }
  
#ifdef DEBUG
  if (cont == 1)
    {
      double eps2 = 1.0e-10;
      for (int kr=1; kr<nmb_b-1; ++kr)
	{
	  for (size_t kc=0; kc<voledg.size()-1; ++kc)
	    {
	      const LRSplineVolume* vol1 = voledg[kc].vol_.get();
	      const LRSplineVolume* vol2 = voledg[kc+1].vol_.get();
	      double vol_dist = -1.0;
	      double tang1_ang = -1.0;
	      double tang2_ang = -1.0;
	      double tang3_ang = -1.0;
	      Point par1 = bsplines[kc][0][kr]->getGrevilleParameter();
	      Point par2 = bsplines[kc+1][0][kr]->getGrevilleParameter();
	      try
		{
		  volumeDifference(vol1, par1, vol2, par2, 
				   vol_dist, tang1_ang, tang2_ang, tang3_ang);
		}
	      catch (...)
		{
		  MESSAGE("Fail ...");
		}
	      if (vol_dist > eps2 || tang1_ang > eps2 || tang2_ang > eps2 || tang3_ang > eps2)
	      	{
		  std::cout << "kr" << kr << ", vol_dist = " << vol_dist << ", tang1_ang = " << tang1_ang;
		  std::cout << ", tang2_ang = " << tang2_ang << ", tang3_ang = " << tang3_ang << std::endl;
		}
	    }
	}
    }
#endif
  return true;
}

//==============================================================================
bool LRVolStitch::averageBdSf(vector<VolBdSf>& volsf, int cont, double tol)
//==============================================================================
{
  int dim = volsf[0].vol_->dimension();
  cont = std::min(1, cont);  // A most C1 continuity 
  if (dim != 1)
    cont = 0;   // C1 stitching only implemented for 1D volumes
  int element_width = cont+1;
  
  // Remove empty edge instances
  for (size_t ki=0; ki<volsf.size(); )
    {
      if (!volsf[ki].vol_.get())
	volsf.erase(volsf.begin()+ki);
      else
	++ki;
    }

  if (volsf.size() != 2)
    return true; // Expects two matching surfaces at interfaces

  Direction3D dir = volsf[0].dir_;
  if (volsf[1].dir_ != dir)
    return false;
  
  // Check corner consistency
  Point corner1[4], corner2[4];
  int tmp1 = (volsf[0].dir_ == XDIR) ? 0 : ((volsf[0].dir_ == YDIR) ? 2 : 4);
  int bd1 = tmp1 + volsf[0].ix_;
  int tmp2 = (volsf[1].dir_ == XDIR) ? 0 : ((volsf[1].dir_ == YDIR) ? 2 : 4);
  int bd2 = tmp2 + volsf[1].ix_;
  fetchSfCorners(volsf[0].vol_, bd1, corner1);
  fetchSfCorners(volsf[1].vol_, bd2, corner2);
  for (int ki=0; ki<4; ++ki)
    {
      double dist = corner1[ki].dist(corner2[0]);
      for (int kj=1; kj<4; ++kj)
	{
	  double dist2 = corner1[ki].dist(corner2[kj]);
	  dist = std::min(dist,dist2);
	}
      if (dist > tol)
	{
	std::cout << "Bd sf, corner mismatch" << std::endl;
	return false;
	}
    }
  
  // Extract boundary B-splines of the two volumes and compute the
  // average coefficients
  // Assumes now that the volumes are oriented equally
  
  // Fetch B-splines along volume boundary surfaces
  vector<vector<LRBSpline3D*> > bsplines1(element_width);
  extractBdSfBsplines(volsf[0], cont, bsplines1);

  vector<vector<LRBSpline3D*> > bsplines2(element_width);
  extractBdSfBsplines(volsf[1], cont, bsplines2);

  // Traverse B-splines and compute average coefficients
   // Sort B-splines to ensure corresponding sequences
  for (int ki=0; ki<=cont; ++ki)
    {
      if (bsplines1[ki].size() != bsplines2[ki].size())
      {
	std::cout << "Bd sf, bspline mismatch" << std::endl;
	std::ofstream of1("bdsfmis1.g2");
	writeElems(of1, volsf[0].vol_);
	std::ofstream of2("bdsfmis2.g2");
	writeElems(of2, volsf[1].vol_);
	return false;
      }
      std::sort(bsplines1[ki].begin(), bsplines1[ki].end(), 
		(dir==XDIR) ? bspline_sort_in_vw :
		((dir==YDIR) ? bspline_sort_in_uw :  bspline_sort_in_uv));
      std::sort(bsplines2[ki].begin(), bsplines2[ki].end(), 
		(dir==XDIR) ? bspline_sort_in_vw :
		((dir==YDIR) ? bspline_sort_in_uw :  bspline_sort_in_uv));
    }
  for (size_t kj=0; kj<bsplines1[0].size(); ++kj)
    {
      // The edge and corner coefs are handled separately and not included in the
      // arrays
      Point coef1 = bsplines1[0][kj]->Coef();
      double gamma1 = bsplines1[0][kj]->gamma();
      Point coef2 = bsplines2[0][kj]->Coef();
      double gamma2 = bsplines2[0][kj]->gamma();

      if (cont == 0)
	{
	  Point coef = 0.5*(coef1 + coef2);
	  bsplines1[0][kj]->setCoefAndGamma(coef, gamma1);
	  bsplines2[0][kj]->setCoefAndGamma(coef, gamma2);
	}
      else // cont == 1
	{
	  // @@@ VSK, 150313. Must be modified if dim > 1
	  // C1 continuity is satisfied if all 4 coefficients lies
	  // along a line. 
	  // Compute least squares linear polynomial: f(t) = at + b
	  LRBSpline3D* bsp[4];
	  bsp[0] = bsplines1[1][kj];
	  bsp[1] = bsplines1[0][kj];
	  bsp[2] = bsplines2[0][kj];
	  bsp[3] = bsplines2[1][kj];
	  makeLineC1(bsp, dir);
	}
    }
  
#ifdef DEBUG
  if (cont == 1)
    {
      double eps2 = 1.0e-10;
      for (size_t kj=0; kj<bsplines1[0].size(); ++kj)
	{
	  const LRSplineVolume* vol1 = volsf[0].vol_.get();
	  const LRSplineVolume* vol2 = volsf[1].vol_.get();
	  double vol_dist = -1.0;
	  double tang1_ang = -1.0;
	  double tang2_ang = -1.0;
	  double tang3_ang = -1.0;
	  Point par1 = bsplines1[0][kj]->getGrevilleParameter();
	  Point par2 = bsplines2[0][kj]->getGrevilleParameter();
	  try
	    {
	      volumeDifference(vol1, par1, vol2, par2, 
			       vol_dist, tang1_ang, tang2_ang, tang3_ang);
	    }
	  catch (...)
	    {
	      MESSAGE("Fail ...");
	    }
	  if (vol_dist > eps2 || tang1_ang > eps2 || tang2_ang > eps2 || tang3_ang > eps2)
	    {
	      std::cout << "kj" << kj << ", vol_dist = " << vol_dist << ", tang1_ang = " << tang1_ang;
	      std::cout << ", tang2_ang = " << tang2_ang << ", tang3_ang = " << tang3_ang << std::endl;
	    }
	}
    }
#endif
  
  return true;
}

//==============================================================================
void LRVolStitch::makeLineC1(LRBSpline3D* bsp[4], Direction3D dir)
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
      par[ki] = t;
      coef[ki] = bsp[ki]->Coef();
      int mult1 = bsp[ki]->endmult_u(true);
      int mult2 = bsp[ki]->endmult_u(false);
      int mult3 = bsp[ki]->endmult_v(true);
      int mult4 = bsp[ki]->endmult_v(false);
      int mult5 = bsp[ki]->endmult_w(true);
      int mult6 = bsp[ki]->endmult_w(false);
      int totmult = mult1 + mult2 + mult3 + mult4 + mult5 + mult6;
      if (totmult > 8)
	std::cout << "Check boundary B-spline" << std::endl;
    }

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
  
#ifdef DEBUG
  // Replace coefficients
  double parn[4];
  parn[0] = (dir == XDIR) ? bsp[0]->umin() : ((dir == YDIR) ? bsp[0]->vmin() : bsp[0]->wmin());
  parn[1] = (dir == XDIR) ? bsp[1]->umax() : ((dir == YDIR) ? bsp[1]->vmax() : bsp[1]->wmax());
  parn[2] = (dir == XDIR) ? bsp[2]->umin() : ((dir == YDIR) ? bsp[2]->vmin() : bsp[2]->wmin());
  parn[3] = (dir == XDIR) ? bsp[3]->umax() : ((dir == YDIR) ? bsp[3]->vmax() : bsp[3]->wmax());
  Point coefn = ((parn[3]-parn[2])*coef[0] + (parn[1]-parn[0])*coef[3])/(parn[3]-parn[0]);
  
  double gamma0 = bsp[0]->gamma();
  double gamma1 = bsp[1]->gamma();
  double gamma2 = bsp[2]->gamma();
  double gamma3 = bsp[3]->gamma();
  if (fabs(gamma1-gamma0) > 1.0e-10 || fabs(gamma2-gamma0) > 1.0e-10 ||
      fabs(gamma3-gamma0) > 1.0e-10)
    {
      std::cout << "Bspline0: (" << bsp[0]->umin() <<"," << bsp[0]->umax();
      std::cout << "," << bsp[0]->vmin() << "," << bsp[0]->vmax();
      std::cout << "), gamma = " << bsp[0]->gamma() << std::endl;
      std::cout << "Bspline1: (" << bsp[1]->umin() <<"," << bsp[1]->umax();
      std::cout << "," << bsp[1]->vmin() << "," << bsp[1]->vmax();
      std::cout << "), gamma = " << bsp[1]->gamma() << std::endl;
      std::cout << "Bspline2: (" << bsp[2]->umin() <<"," << bsp[2]->umax();
      std::cout << "," << bsp[2]->vmin() << "," << bsp[2]->vmax();
      std::cout << "), gamma = " << bsp[2]->gamma() << std::endl;
      std::cout << "Bspline3: (" << bsp[3]->umin() <<"," << bsp[3]->umax();
      std::cout << "," << bsp[3]->vmin() << "," << bsp[3]->vmax();
      std::cout << "), gamma = " << bsp[3]->gamma() << std::endl;
    }
#endif
#if 0
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
  // double gamma = bsp[1]->gamma();
  // bsp[1]->setCoefAndGamma(coefn, gamma);
  // bsp[2]->setCoefAndGamma(coefn, gamma);
  for (ki=0; ki<4; ++ki)
    {
      Point coefn2(1);
      coefn2.setValue(a*par[ki] + b);
      double gamma = bsp[ki]->gamma();
      bsp[ki]->setCoefAndGamma(coefn2, gamma);
    }
}

//==============================================================================
void LRVolStitch::fetchEdgeCorners(const VolEdge& voledg, Point& c1, Point& c2)
//==============================================================================
{
  double u1, u2, v1, v2, w1, w2;
  if (voledg.dir_ == XDIR)
    {
      u1 = voledg.vol_->startparam_u();
      u2 = voledg.vol_->endparam_u();
      v1 = v2 = (voledg.ix1_ == 0) ? voledg.vol_->startparam_v() : voledg.vol_->endparam_v();
      w1 = w2 = (voledg.ix2_ == 0) ? voledg.vol_->startparam_w() : voledg.vol_->endparam_w();
    }
  else if (voledg.dir_ == YDIR)
    {
      u1 = u2 = (voledg.ix2_ == 0) ? voledg.vol_->startparam_u() : voledg.vol_->endparam_u();
      v1 = voledg.vol_->startparam_v();
      v2 = voledg.vol_->endparam_v();
      w1 = w2 = (voledg.ix1_ == 0) ? voledg.vol_->startparam_w() : voledg.vol_->endparam_w();
    }
  else
    {
      u1 = u2 = (voledg.ix1_ == 0) ? voledg.vol_->startparam_u() : voledg.vol_->endparam_u();
      v1 = v2 = (voledg.ix2_ == 0) ? voledg.vol_->startparam_v() : voledg.vol_->endparam_v();
      w1 = voledg.vol_->startparam_w();
      w2 = voledg.vol_->endparam_w();
    }

  c1 = Point(u1, v1, w1);
  c2 = Point(u2, v2, w2);
}

//==============================================================================
void LRVolStitch::fetchSfCorners(shared_ptr<LRSplineVolume> vol, int bdsf, Point corner[])
//==============================================================================
{
  double u1, u2, v1, v2, w1, w2;
  if (bdsf <= 1)
    {
      u1 = (bdsf == 0) ? vol->startparam_u() : vol->endparam_u();
      v1 = vol->startparam_v();
      v2 = vol->endparam_v();
      w1 = vol->startparam_w();
      w2 = vol->endparam_w();
      corner[0] = Point(u1, v1, w1);
      corner[1] = Point(u1, v2, w1);
      corner[2] = Point(u1, v1, w2);
      corner[3] = Point(u1, v2, w2);
    }
  else if (bdsf <= 3)
    {
      u1 = vol->startparam_u();
      u2 = vol->endparam_u();
      v1 = (bdsf == 2) ? vol->startparam_v() : vol->endparam_v();
      w1 = vol->startparam_w();
      w2 = vol->endparam_w();
      corner[0] = Point(u1, v1, w1);
      corner[1] = Point(u2, v1, w1);
      corner[2] = Point(u1, v1, w2);
      corner[3] = Point(u2, v1, w2);
    }
  else 
    {
      u1 = vol->startparam_u();
      u2 = vol->endparam_u();
      v1 = vol->startparam_v();
      v2 = vol->endparam_v();
      w1 = (bdsf == 4) ? vol->startparam_w() : vol->endparam_w();
      corner[0] = Point(u1, v1, w1);
      corner[1] = Point(u2, v1, w1);
      corner[2] = Point(u1, v2, w1);
      corner[3] = Point(u2, v2, w1);
    }

}

//==============================================================================
void LRVolStitch::defineRefinements(const Mesh3D& mesh1, const Mesh3D& mesh2,
				    Direction3D dir,
				    vector<pair<vector<GPos2D>,int> >& mrects,
				    bool at_start, int width, int ix,
				    vector<LRSplineVolume::Refinement3D>& refs)
//==============================================================================
{
  double start2, end2;
  Direction3D dir2 = (ix == 0) ? next(dir) : prev(dir);
  Direction3D dir3 = (ix == 1) ? next(dir) : prev(dir);
  int num = mesh2.numDistinctKnots(dir2);
  if (at_start)
    {
      start2 = mesh2.kval(dir2, 0);
      end2 = mesh2.kval(dir2, width);
    }
  else
    {
      start2 = mesh2.kval(dir2, num-width-1);
      end2 = mesh2.kval(dir2, num-1);
    }
  
  for (size_t ki=0; ki<mrects.size(); ++ki)
    {
      double kval = mesh1.kval(dir, mrects[ki].second);
      for (size_t kj=0; kj<mrects[ki].first.size(); ++kj)
	{
	  double start1 = mesh1.kval(dir3, mrects[ki].first[kj].ll[1-ix]);
	  double end1 = mesh1.kval(dir3, mrects[ki].first[kj].ur[1-ix]);
	  int mult = mrects[ki].first[kj].mult;

	  LRSplineVolume::Refinement3D curr;
	  if (ix == 0)
	    curr = LRSplineVolume::Refinement3D(kval, start2, end2,
						start1, end1, dir, mult);
	  else
	    curr = LRSplineVolume::Refinement3D(kval, start1, end1,
						start2, end2, dir, mult);
	  refs.push_back(curr);
	}
    }
}

//==============================================================================
void LRVolStitch::adaptMeshRectangles(vector<pair<vector<GPos2D>,int> >& rects,
				      int ix, int start, int end, int t1, int t2,
				      int del, vector<pair<vector<GPos2D>,int> >& rects2)
//==============================================================================
{
  int ix2 = 1 - ix;
  for (size_t kr=0; kr<rects.size(); ++kr)
    {
      vector<GPos2D> mrects;
      for (size_t kh=0; kh<rects[kr].first.size(); ++kh)
	{
	  if (rects[kr].first[kh].ll[ix2] == t1 && rects[kr].first[kh].ur[ix2] < t1+del)
	    rects[kr].first[kh].ur[ix2] = t1 + del;
	  if (rects[kr].first[kh].ur[ix2] == t2 && rects[kr].first[kh].ll[ix2] > t2-del)
	    rects[kr].first[kh].ll[ix2] = t2 - del;
	  if (start < rects[kr].first[kh].ll[ix])
	    {
	      GPos2D tmp = rects[kr].first[kh];
	      tmp.ll[ix] = start;
	      tmp.ur[ix] = rects[kr].first[kh].ll[ix];
	      mrects.push_back(tmp);
	    }
	  if (end > rects[kr].first[kh].ur[ix])
	    {
	      GPos2D tmp = rects[kr].first[kh];
	      tmp.ll[ix] = rects[kr].first[kh].ur[ix];
	      tmp.ur[ix] = end;
	      mrects.push_back(tmp);
	    }
	  rects[kr].first[kh].ll[ix] = start;
	  rects[kr].first[kh].ur[ix] = end;
	}
      if (mrects.size() > 0)
	rects2.push_back(make_pair(mrects, rects[kr].second));
    }
}


//==============================================================================
void LRVolStitch::extractBoundaryBsplines(const VolEdge& voledg, int cont,
					  vector<vector<LRBSpline3D*> >& bsplines)
//==============================================================================
{
  // Traverse B-splines and extract those that have sufficient multiplicity
  // at the specified edge. 
  int nmbel = 1 + 3*cont;
  // Extract the number of B-splines in the opposite direction specified
  // by element width
  // First define one-sided key values
  int ki, kj, kr;
  double eps = 1.0e-12;
  bool at_start1 = (voledg.ix1_ == 0);
  bool at_start2 = (voledg.ix2_ == 0);
  Direction3D dir = voledg.dir_;
  Direction3D dir1 = next(dir);
  Direction3D dir2 = prev(dir);
  int deg1 = voledg.vol_->degree(dir1);
  int deg2 = voledg.vol_->degree(dir2);
  const Mesh3D& mesh = voledg.vol_->mesh();
  int ix1 = (at_start1) ? 0 : mesh.numDistinctKnots(dir1) - 1;
  int ix2 = (at_start2) ? 0 : mesh.numDistinctKnots(dir2) - 1;
  vector<double> lim(2*(cont+1));
  vector<int> mm(2*(cont+1));
  for (ki=0; ki<=cont; ++ki)
    {
      mm[2*ki] = deg1 + 1 - ki;
      mm[2*ki+1] = deg2 + 1 - ki;
      lim[2*ki] = mesh.kval(dir1, ix1);
      lim[2*ki+1] = mesh.kval(dir2, ix2);
      // lim[4*ki] = mesh.kval(dir1, at_start1 ? ix1 : ix1-ki-1);
      // lim[4*ki+1] = mesh.kval(dir1, at_start1 ? ix1+ki+1 : ix1);
      // lim[4*ki+2] = mesh.kval(dir2, at_start2 ? ix2 : ix2-ki-1);
      // lim[4*ki+3] = mesh.kval(dir2, at_start2 ? ix2+ki+1 : ix2);
    }

  for (ki=0; ki<nmbel; ++ki)
    bsplines[ki].clear();
  for (LRSplineVolume::BSplineMap::iterator it=voledg.vol_->basisFunctionsBeginNonconst(); 
       it != voledg.vol_->basisFunctionsEndNonconst(); ++it)
	{
	  int mult1 = it->second->endmult(dir1, at_start1);
	  int mult2 = it->second->endmult(dir2, at_start2);

	  double t1 = mesh.kval(dir1, at_start1 ? it->second->suppMin(dir1) :
				it->second->suppMax(dir1));
	  double t2 = mesh.kval(dir2, at_start2 ? it->second->suppMin(dir2) :
				it->second->suppMax(dir2));
	  // double t1 = mesh.kval(dir1, it->second->suppMin(dir1));
	  // double t2 = mesh.kval(dir1, it->second->suppMax(dir1));
	  // double t3 = mesh.kval(dir2, it->second->suppMin(dir2));
	  // double t4 = mesh.kval(dir2, it->second->suppMax(dir2));

	  for (ki=0, kj=0; ki<1+3*cont; ki+=2, ++kj)
	    {
	      if (mult1 == mm[2*kj] && mult2 == mm[2*kj+1] &&
		  fabs(t1-lim[2*kj]) < eps && fabs(t2-lim[2*kj+1]) < eps)
		  // fabs(t1-lim[4*kj]) < eps && fabs(t2-lim[2*kj+1]) < eps &&
		  // fabs(t3-lim[4*kj+2]) < eps && fabs(t4-lim[4*kj+3]) < eps)
		bsplines[ki].push_back(it->second.get());
	    }
	  if (cont == 1)
	    {
	      for (ki=1, kj=0, kr=1; kr<=3; ki=0, kj=1, kr+=2)
		{
		  if (mult1 == mm[2*ki] && mult2 == mm[2*kj+1] &&
		      fabs(t1-lim[2*ki]) < eps && fabs(t2-lim[2*kj+1]) < eps)
		      // fabs(t1-lim[4*ki]) < eps && fabs(t2-lim[4*ki+1]) < eps &&
		      // fabs(t3-lim[4*kj+2]) < eps && fabs(t4-lim[4*kj+3]) < eps)
		    bsplines[kr].push_back(it->second.get());
		}
	    }
	}
  int stop_break = 1;
}

//==============================================================================
void LRVolStitch::extractBdSfBsplines(const VolBdSf& volsf, int cont,
				      vector<vector<LRBSpline3D*> >& bsplines)
//==============================================================================
{
  // Traverse B-splines and extract those that have sufficient multiplicity
  // at the specified edge. 
  // Extract the number of B-splines in the opposite direction specified
  // by element width
  // First define one-sided key values
  int ki, kj, kr;
  double eps = 1.0e-12;
  bool at_start = (volsf.ix_ == 0);
  Direction3D dir = volsf.dir_;
  Direction3D dir1 = next(dir);
  Direction3D dir2 = prev(dir);
  int deg = volsf.vol_->degree(dir);
  int deg1 = volsf.vol_->degree(dir1);
  int deg2 = volsf.vol_->degree(dir2);
  const Mesh3D& mesh = volsf.vol_->mesh();
  int ix = (at_start) ? 0 : mesh.numDistinctKnots(dir) - 1;
  int num1 = mesh.numDistinctKnots(dir1);
  int num2 = mesh.numDistinctKnots(dir2);
  vector<double> lim(6*(cont+1));
  vector<int> mm(3*(cont+1));
  for (ki=0; ki<=cont; ++ki)
    {
      mm[3*ki] = deg + 1 - ki;
      mm[3*ki+1] = deg1 + 1 - ki;
      mm[3*ki+2] = deg2 + 1 - ki;
      lim[6*ki] = mesh.kval(dir, at_start ? ix : ix-ki-1);
      lim[6*ki+1] = mesh.kval(dir, at_start ? ix+ki+1 : ix);
      lim[6*ki+2] = mesh.kval(dir1, 0);
      //lim[10*ki+3] = mesh.kval(dir1, ki+1);
      //lim[10*ki+4] = mesh.kval(dir1, num1-ki-2);
      lim[6*ki+3] = mesh.kval(dir1, num1-1);
      lim[6*ki+4] = mesh.kval(dir2, 0);
      //lim[10*ki+7] = mesh.kval(dir2, ki+1);
      //lim[10*ki+8] = mesh.kval(dir2, num2-ki-2);
      lim[6*ki+5] = mesh.kval(dir2, num2-1);
    }

  for (ki=0; ki<=cont; ++ki)
    bsplines[ki].clear();
  for (LRSplineVolume::BSplineMap::iterator it=volsf.vol_->basisFunctionsBeginNonconst(); 
       it != volsf.vol_->basisFunctionsEndNonconst(); ++it)
	{
	  int mult = it->second->endmult(dir, at_start);
	  double t1 = mesh.kval(dir, it->second->suppMin(dir));
	  double t2 = mesh.kval(dir, it->second->suppMax(dir));
	  for (ki=0; ki<=cont; ++ki)
	    {
	      if (mult == mm[3*ki] && fabs(t1-lim[6*ki]) < eps &&
		  fabs(t2-lim[6*ki+1]) < eps)
		{
		  // Check if the bspline is situated at an edge or corner
		  int mult1 = it->second->endmult(dir1, true);
		  int mult2 = it->second->endmult(dir1, false);
		  int mult3 = it->second->endmult(dir2, true);
		  int mult4 = it->second->endmult(dir2, false);
		  double t3 = mesh.kval(dir1, it->second->suppMin(dir1));
		  double t4 = mesh.kval(dir1, it->second->suppMax(dir1));
		  double t5 = mesh.kval(dir2, it->second->suppMin(dir2));
		  double t6 = mesh.kval(dir2, it->second->suppMax(dir2));

		  bool accept = true;
		  if (cont == 0 &&
		      ((mult1 == mm[3*cont+1] &&
			fabs(t3-lim[6*cont+2]) < eps /*&& fabs(t4-lim[6*cont+3]) < eps*/) ||
		       (mult2 == mm[3*cont+1] &&
			fabs(t4-lim[6*cont+3]) < eps /*&& fabs(t4-lim[6*cont+5]) < eps)*/)) ||
		      ((mult3 == mm[3*cont+2] &&
			fabs(t5-lim[6*cont+4]) < eps /*&& fabs(t6-lim[6*cont+7]) < eps*/) ||
		       (mult4 == mm[3*cont+2] &&
			fabs(t6-lim[6*cont+5]) < eps /*&& fabs(t6-lim[6*cont+9]) < eps*/)))
		    {
		      accept = false;  
		    }
		  else if (cont == 1 &&
		      ((mult1 >= mm[3*cont+1] &&
			fabs(t3-lim[6*cont+2]) < eps /*&& fabs(t4-lim[6*cont+3]) < eps*/) ||
		       //t3 >= lim[6*cont+2]-eps && t4 <= lim[6*cont+3]+eps) ||
		       (mult2 >= mm[3*cont+1] &&
			fabs(t4-lim[6*cont+3]) < eps /*&& fabs(t4-lim[6*cont+5]) < eps)*/)) ||
		       //t3 >= lim[6*cont+4]-eps && t4 <= lim[6*cont+5]+eps)) &&
		      ((mult3 >= mm[3*cont+2] &&
			fabs(t5-lim[6*cont+4]) < eps /*&& fabs(t6-lim[6*cont+7]) < eps*/) ||
		       //t5 >= lim[6*cont+6]-eps && t6 <= lim[6*cont+7]+eps) ||
		       (mult4 >= mm[3*cont+2] &&
			fabs(t6-lim[6*cont+5]) < eps /*&& fabs(t6-lim[6*cont+9]) < eps*/)))
			   //t5 >= lim[6*cont+8]-eps && t6 <= lim[6*cont+9]+eps)))
		    {
		      accept = false;
		    }

		  //accept = true;
		  if (accept)
		    {
		      int totmult = mult1+mult2+mult3+mult4;
		      int multlimit = (cont == 0) ? 8 : 4;
		      if (totmult > multlimit || mult1 > 2-cont || mult2 > 2-cont ||
			  mult3 > 2-cont || mult4 > 2-cont)
			std::cout << "Check sf bd B-spline" << std::endl;
		      bsplines[ki].push_back(it->second.get());
		    }
		  // else
		  //   std::cout << "Boundary B-spline" << std::endl;
		}
	    }
	}

  int stop_break = 1;
}

void LRVolStitch::volumeDifference(const LRSplineVolume* vol1, const Point& param1,
				   const LRSplineVolume* vol2, const Point& param2,
				   double& vol_dist, double& tang1_ang,
				   double& tang2_ang, double& tang3_ang)
{
    const int derivs = 1;
    const int totpts = 4;
    vector<Point> vol_pt1(totpts), vol_pt2(totpts);
    vol1->point(vol_pt1, param1[0], param1[1], param1[2], derivs);
    vol2->point(vol_pt2, param2[0], param2[1], param2[2], derivs);

    vol_dist = vol_pt1[0].dist(vol_pt2[0]);
    tang1_ang = vol_pt1[1].angle(vol_pt2[1]);
    tang2_ang = vol_pt1[2].angle(vol_pt2[2]);
    tang3_ang = vol_pt1[3].angle(vol_pt2[3]);

}

void LRVolStitch::writeMrects(std::ostream& out, const Mesh3D& mesh,
			      vector<pair<vector<GPos2D>,int > >& mrects,
			      Direction3D dir, int ix)
{
  vector<double> bd_rect;
  for (size_t ki=0; ki<mrects.size(); ++ki)
    {
      int k_ix = mrects[ki].second;
      for (size_t kj=0; kj<mrects[ki].first.size(); ++kj)
	{
	  Point ll(3), ur(3);
	  int ix0, ix1, ix2;
	  Direction3D dir1 = next(dir);
	  Direction3D dir2 = prev(dir);
	  ix0 = (dir == XDIR) ? 0 : ((dir == YDIR) ? 1 : 2);
	  ix1 = (ix0 + 1)%3;
	  ix2 = (ix0 + 2)%3;
	  ll[ix0] = ur[ix0] = mesh.kval(dir, k_ix);
	  ll[ix1] = mesh.kval(dir1, mrects[ki].first[kj].ll[0]);
	  ll[ix2] = mesh.kval(dir2, mrects[ki].first[kj].ll[1]);
	  ur[ix1] = mesh.kval(dir1, mrects[ki].first[kj].ur[0]);
	  ur[ix2] = mesh.kval(dir2, mrects[ki].first[kj].ur[1]);

	  Point pos1 = ll;
	  Point pos2 = pos1;
	  pos2[ix1] = ur[ix1];
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos1[ka]);
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos2[ka]);

	  pos1 = pos2;
	  pos2 = ur;
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos1[ka]);
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos2[ka]);

	  pos1 = ll;
	  pos1[ix2] = ur[ix2];
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos1[ka]);
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos2[ka]);

	  pos2 = ll;
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos1[ka]);
	  for (int ka=0; ka<3; ++ka)
	    bd_rect.push_back(pos2[ka]);
	}
    }
  LineCloud rects_lines(bd_rect.begin(), bd_rect.size()/6);
  if (ix == 0)
    out << "410 1 0 4 255 0 0 255" << std::endl;
  else if (ix == 1)
    out << "410 1 0 4 20 235 0 255" << std::endl;
  else if (ix == 2)
    out << "410 1 0 4 0 0 255 255" << std::endl;
  else if (ix == 3)
    out << "410 1 0 4 100 155 0 255" << std::endl;
  else if (ix == 4)
    out << "410 1 0 4 0 155 100 255" << std::endl;
  else if (ix == 5)
    out << "410 1 0 4 81 81 81 255" << std::endl;
  else if (ix == 6)
    out << "410 1 0 4 100 0 155 255" << std::endl;
  else if (ix == 7)
    out << "410 1 0 4 50 155 50 255" << std::endl;
  rects_lines.write(out);
}

void LRVolStitch::writeElems(std::ostream& out, shared_ptr<LRSplineVolume>& vol)
{
  int num_el = vol->numElements();
  vector<double> bd_el(72*num_el);
  int ki=0, kj=0;
  Point pos, pos2, dir1, dir2, dir3;
  for (auto it=vol->elementsBegin(); it != vol->elementsEnd(); ++it)
    {
      const Element3D* elem = it->second.get();
      double umin = elem->umin();
      double umax = elem->umax();
      double vmin = elem->vmin();
      double vmax = elem->vmax();
      double wmin = elem->wmin();
      double wmax = elem->wmax();

      pos = Point(umin, vmin, wmin);
      dir1 = Point(umax-umin, 0.0, 0.0);
      dir2 = Point(0.0, vmax-vmin, 0.0);
      dir3 = Point(0.0, 0.0, wmax-wmin);

      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir1;
      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir2;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir3;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir1;
      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir2;
      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];
     }


  LineCloud elem_lines(bd_el.begin(), bd_el.size()/6);
  elem_lines.writeStandardHeader(out);
  elem_lines.write(out);
}

void LRVolStitch::writeBSplines(std::ostream& out, vector<LRBSpline3D*>& bsplines)
{
  vector<double> bd_b(72*bsplines.size());
  int ki=0, kj=0;
  Point pos, pos2, dir1, dir2, dir3;
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      double umin = bsplines[ki]->umin();
      double umax = bsplines[ki]->umax();
      double vmin = bsplines[ki]->vmin();
      double vmax = bsplines[ki]->vmax();
      double wmin = bsplines[ki]->wmin();
      double wmax = bsplines[ki]->wmax();

      pos = Point(umin, vmin, wmin);
      dir1 = Point(umax-umin, 0.0, 0.0);
      dir2 = Point(0.0, vmax-vmin, 0.0);
      dir3 = Point(0.0, 0.0, wmax-wmin);

      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos += dir1;
      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos += dir2;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos += dir3;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos -= dir1;
      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos2 = pos-dir3;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];

      pos -= dir2;
      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_b[kj++] = pos2[ka];
     }


  LineCloud bsplinelines(bd_b.begin(), bd_b.size()/6);
  bsplinelines.writeStandardHeader(out);
  bsplinelines.write(out);
}
