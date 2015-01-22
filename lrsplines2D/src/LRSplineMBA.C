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

#include "GoTools/lrsplines2D/LRSplineMBA.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Element2D.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;
using std::map;
using std::array;
using namespace Go;

//==============================================================================
void LRSplineMBA::MBAUpdate(LRSplineSurface *srf)
//==============================================================================
{
 
//   // We start the timer.
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
//   // double time_loop = 0.0;
// #endif

  double tol = 1.0e-12;  // Numeric tolerance

  double umax = srf->endparam_u();
  double vmax = srf->endparam_v();

  // Make a copy of the surface
  shared_ptr<LRSplineSurface> cpsrf(new LRSplineSurface(*srf));

  // Set all coefficients to zero (keep the scaling factors)
  int dim = srf->dimension();
  Point coef(dim);
  coef.setValue(0.0);
  for (LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
       it1 != cpsrf->basisFunctionsEnd(); ++it1)
    cpsrf->setCoef(coef, it1->second.get());
    
  // Map to accumulate numerator and denominator to compute final coefficient value
  // for each BSplineFunction
  map<const LRBSpline2D*, Array<double,4> > nom_denom; 
  //map<const LRBSpline2D*, Point> nom_denom; 

  // Temporary vector to store weights associated with a given data point
  vector<double> tmp(dim);

  vector<double> tmp_weights;

  // Traverse all elements. The two surfaces will have corresponding elements,
  // but only the source surface elements will contain point information so 
  // the elements in both surfaces must be traversed
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  LRSplineSurface::ElementMap::const_iterator el1 = srf->elementsBegin();
  LRSplineSurface::ElementMap::const_iterator el2 = cpsrf->elementsBegin();
  for (; el1!=srf->elementsEnd(); ++el1, ++el2)
    {
      if (!el1->second->hasDataPoints())
	continue;  // No points to use in surface update

      // Fetch points from the source surface
      int nmb_pts = el1->second->nmbDataPoints();
      const vector<double>& points = el1->second->getDataPoints();
      int nmb_ghost = 0; //el1->second->nmbGhostPoints();
      //vector<double>& ghost_points = el1->second->getGhostPoints();
      vector<double> ghost_points;

      // Fetch associated B-splines belonging to the difference surface
      const vector<LRBSpline2D*>& bsplines = el2->second->getSupport();

      const int bsplines_size = bsplines.size();

      // Check if the element needs to be updated
      size_t nb;
      for (nb=0; nb<bsplines.size(); ++nb)
	if (!bsplines[nb]->coefFixed())
	  break;

      if (nb == bsplines.size())
	continue;   // Element satisfies accuracy requirements

      
      // Compute contribution from all points
      int ki, kk;
      size_t kj;
      const double *curr;
      bool u_at_end, v_at_end;
      double total_squared_inv, val, wgt, wc, phi_c, gamma;
      tmp_weights.resize(bsplines.size());
      // std::cout << "tmp_weight.size(): " << tmp_weights.size() << std::endl;
      // std::cout << "nmb_pts: " << nmb_pts << std::endl;
      // std::cout << "points.size(): " << points.size() << std::endl;
      // std::cout << "dim: " << dim << std::endl;
      // std::cout << "del: " << del << std::endl;
      int threadId = 0;

//      printf("tmp_weights.size(): %i\n", tmp_weights.size());
      // @@sbr Not thread safe!
#ifndef _OPENMP
      {
      for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
      {
#else
	  assert(dim == 1);
	  // omp_set_num_threads(1);
	  // pthread_attr_t attr;
	  // size_t stacksize;
	  // int status = pthread_attr_getstacksize(&attr, &stacksize);
	  // std::cout << "status: " << status << std::endl;
	  // status = pthread_attr_setstacksize(&attr, 3);//&attr, &stacksize);
	  // std::cout << "status: " << status << std::endl;
	  // status = pthread_attr_getstacksize(&attr, &stacksize);
	  // std::cout << "status: " << status << std::endl;
	  // std::cout << "stacksize (in MB): " << (double)stacksize/(1024.0*1024.0) << std::endl;
#pragma omp parallel default(none) private(ki, kj, kk, u_at_end, v_at_end, total_squared_inv, val, wgt, wc, phi_c, curr, nom_denom, gamma, threadId) firstprivate(tmp, tmp_weights) shared(nmb_pts, points, umax, vmax, tol, dim, del, bsplines)
	  {
	      // printf("bsplines_size: %i\n", bsplines_size);
	      // threadId = omp_get_thread_num();
	      // snprintf("ki: ", ki, "omp-%02d", threadId);
//	      printf("threadId: %i\n", threadId);
//	    vector<double> tmp_weights(bsplines_size);

#pragma omp for schedule(runtime) //dynamic, 4)//static, 4)//runtime)//guided)//auto)
	      for (ki=0; ki<nmb_pts; ++ki)
	      {
		  curr = &points[ki*del];
		  threadId = omp_get_thread_num();
#endif
//		  snprintf("ki: ", ki, "omp-%02d", threadId);
		  // printf("ki: %i\n", ki);
		  // printf("del: %i\n", del);
		  // printf("points.size(): %i\n", points.size());
		   // printf("curr[0]: %d\n", curr[0]);
		  // printf("curr[1]: %f\n", curr[1]);
//#endif//_OPENMP
		  // Computing weights for this data point
		  u_at_end = (curr[0] > umax-tol) ? true : false;
		  v_at_end = (curr[1] > vmax-tol) ? true : false;
		  total_squared_inv = 0.0;
		  for (kj=0; kj<bsplines.size(); ++kj) 
		  {
		      // printf("umin: %f\n", bsplines[kj]->umin());
		      // printf("vmin: %f\n", bsplines[kj]->vmin());
		      val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
							    u_at_end, v_at_end);
		      gamma = bsplines[kj]->gamma();
		      wgt = val*gamma;//bsplines[kj]->gamma();
		      // printf("kj: %i\n", kj);
		      // printf("tmp_weights.size(): %i\n", tmp_weights.size());
		      tmp_weights[kj] = wgt;
		      total_squared_inv += wgt*wgt;
		      // printf("total_squared_inv: %f\n", total_squared_inv);
		  }
		  // printf("Done with for loop.\n");
		  total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;
		  // printf("total_squared_inv: %f\n", total_squared_inv);

		  // Compute contribution
		  for (kj=0; kj<bsplines.size(); ++kj)
		  {
		      // printf("kj: %i\n", kj);
		      wc = tmp_weights[kj]; 
		      for (kk=0; kk<dim; ++kk)
		      {
			  phi_c = wc * curr[del-dim+kk] * total_squared_inv;
			  tmp[kk] = wc * wc * phi_c;
		      } // @@sbr201412 But is this threadsafe? Do not think so.
 // This construction makes all assignments in the for loop run in serial ...
#pragma omp critical
// #pragma omp barrier // Takes forever ...
//#pragma omp atomic // Single operation (not function) only.
			  add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
			  		   wc * wc);
		  }
		  // printf("Done with for loop.\n");
	      }
	  }
	  // printf("Done with OpenMP.\n");

	  // Compute contribution from ghost points
	  for (ki=0, curr=&ghost_points[0]; ki<nmb_ghost; ++ki, curr+=del)
	  {
	      // Computing weights for this data point
	      bool u_at_end = (curr[0] > umax-tol) ? true : false;
	      bool v_at_end = (curr[1] > vmax-tol) ? true : false;
	      double total_squared_inv = 0;
	      for (kj=0; kj<bsplines.size(); ++kj) 
	      {
		  double val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
							       u_at_end, v_at_end);
		  const double wgt = val*bsplines[kj]->gamma();
		  tmp_weights[kj] = wgt;
		  total_squared_inv += wgt*wgt;
	      }
	      total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;

	      // Compute contribution
	      for (kj=0; kj<bsplines.size(); ++kj)
	      {
		  const double wc = tmp_weights[kj]; 
		  for (int ki=0; ki<dim; ++ki)
		  {
		      const double phi_c = wc * curr[del-dim+ki] * total_squared_inv;
		      tmp[ki] = wc * wc * phi_c;
		  }
		  add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
				   wc * wc);
	      }
	  }
      }

  // Compute coefficients of difference surface
  for (LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
       it1 != cpsrf->basisFunctionsEnd(); ++it1) 
    {
      auto nd_it = nom_denom.find(it1->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ki=0; ki<dim; ++ki)
	coef[ki] = (fabs(entry[dim]<tol)) ? 0 : entry[ki] / entry[dim];
      cpsrf->setCoef(coef, it1->second.get());
    }
 
  // Update initial surface
  double fac = 1.0; //1.01;
  srf->addSurface(*cpsrf, fac);

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in MBAUpdate(): " << time_spent << std::endl;
//   // std::cout << "time_spent in for loop: " << time_loop << std::endl;
// #endif
 }

//------------------------------------------------------------------------------
void LRSplineMBA::add_contribution(int dim, 
				   //map<const LRBSpline2D*, Point>& target, 
				   map<const LRBSpline2D*, Array<double,4> >& target, 
				   const LRBSpline2D* bspline, double nom[], 
				   double denom)
//------------------------------------------------------------------------------
{
   auto it = target.find(bspline);
   if (it != target.end()) 
     {
       // already in map
       for (int ki=0; ki<dim; ++ki)
	 it->second[ki] += nom[ki];
       it->second[dim] += denom;
     } 
   else 
     {
     // not already in map.  Insert it
       vector<double> tmp(dim+1, 0.0);
       for (int ki=0; ki<dim; ++ki)
	 tmp[ki] = nom[ki];
       tmp[dim] = denom;
       //target.insert({bspline, Point(tmp.begin(), tmp.end())});
       //target.insert({bspline, Array<double,2>{tmp[0], denom}});
       target.insert({bspline, Array<double,4>(tmp.begin())});
     }
 }

