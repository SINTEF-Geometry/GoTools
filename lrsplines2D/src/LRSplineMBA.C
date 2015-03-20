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
#include "GoTools/geometry/Utils.h"

#include <iostream>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;
using std::set;
using std::map;
using std::array;
using std::cout;
using std::endl;
using namespace Go;

//==============================================================================
void LRSplineMBA::MBADistAndUpdate(LRSplineSurface *srf)
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
  int order2 = (srf->degree(XFIXED)+1)*(srf->degree(YFIXED)+1);

  // Make a copy of the surface
  shared_ptr<LRSplineSurface> cpsrf(new LRSplineSurface(*srf));

  // Set all coefficients to zero (keep the scaling factors)
  int dim = srf->dimension();
  vector<double> ptval(dim);
  Point coef(dim);
  coef.setValue(0.0);
  for (LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
       it1 != cpsrf->basisFunctionsEnd(); ++it1)
    cpsrf->setCoef(coef, it1->second.get());
    
  // Map to accumulate numerator and denominator to compute final coefficient value
  // for each BSplineFunction
  map<const LRBSpline2D*, Array<double,2> > nom_denom; 

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

      // Fetch associated B-splines belonging to the difference surface
      const vector<LRBSpline2D*>& bsplines = el1->second->getSupport();

      // Check if the element needs to be updated
      size_t nb;
      for (nb=0; nb<bsplines.size(); ++nb)
	if (!bsplines[nb]->coefFixed())
	  break;

      if (nb == bsplines.size())
	continue;   // Element satisfies accuracy requirements

     // Fetch points from the source surface
      int nmb_pts = el1->second->nmbDataPoints();
      vector<double>& points = el1->second->getDataPoints();
      int nmb_ghost = 0; //el1->second->nmbGhostPoints();
      //vector<double>& ghost_points = el1->second->getGhostPoints();
      //vector<double> ghost_points;

      tmp_weights.resize(bsplines.size());
      
      // Compute contribution from all points
      // First compute distance in the data sets and store 
      // basis function values
      int ki, kr;
      size_t kj;
      double *curr;
      vector<double> Bval;
      Bval.reserve(1.5*nmb_pts*order2);  // This vector is probably too large
      for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
	{
	  // Computing weights for this data point
	  bool u_at_end = (curr[0] > umax-tol) ? true : false;
	  bool v_at_end = (curr[1] > vmax-tol) ? true : false;
	  double total_squared_inv = 0;
	  std::fill(ptval.begin(), ptval.end(), 0.0);
	  for (kj=0; kj<bsplines.size(); ++kj) 
	    {
	      double val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
							   u_at_end, v_at_end);
	      Bval.push_back(val);
	      const Point& tmp = bsplines[kj]->coefTimesGamma();
	      for (int ka=0; ka<dim; ++ka)
		ptval[ka] += val*tmp[ka];
	    }
	  double dist;
	  if (dim == 1)
	    dist = curr[2] - ptval[0];
	  else
	    {
	      dist = Utils::distance_squared(ptval.begin(), ptval.end(),
					     points.begin()+ki*del+2); 
	      //ptval.dist(Point(curr+2, curr+del));
	      dist = sqrt(dist);
	    }
	  curr[del-1] = dist;
	}

      for (ki=0, kr=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
	{
	  // Computing weights for this data point
	  double total_squared_inv = 0;
	  for (kj=0; kj<bsplines.size(); ++kj, ++kr) 
	    {
	      double val = Bval[kr];
	      const double wgt = val*bsplines[kj]->gamma();
	      tmp_weights[kj] = wgt;
	      total_squared_inv += wgt*wgt;
	    }
	  total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;

	  // Compute contribution
	  for (kj=0; kj<bsplines.size(); ++kj)
	    {
	      const double wc = tmp_weights[kj]; 
	      for (int ka=0; ka<dim; ++ka)
		{
		  const double phi_c = wc * curr[del-dim+ka] * total_squared_inv;
		  tmp[ka] = wc * wc * phi_c;
		}
	      add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
			       wc * wc);
	    }
	}

     //  // Compute contribution from ghost points
     //  for (ki=0, curr=&ghost_points[0]; ki<nmb_ghost; ++ki, curr+=del)
     // 	{
     // 	  // Computing weights for this data point
     // 	  bool u_at_end = (curr[0] > umax-tol) ? true : false;
     // 	  bool v_at_end = (curr[1] > vmax-tol) ? true : false;
     // 	  double total_squared_inv = 0;
     // 	  for (kj=0; kj<bsplines.size(); ++kj) 
     // 	    {
     // 	      double val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
     // 							   u_at_end, v_at_end);
     // 	      const double wgt = val*bsplines[kj]->gamma();
     // 	      tmp_weights[kj] = wgt;
     // 	      total_squared_inv += wgt*wgt;
     // 	    }
     // 	  total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;

     // 	  // Compute contribution
     // 	  for (kj=0; kj<bsplines.size(); ++kj)
     // 	    {
     // 	      const double wc = tmp_weights[kj]; 
     // 	      for (int ka=0; ka<dim; ++ka)
     // 		{
     // 		  const double phi_c = wc * curr[del-dim+ka] * total_squared_inv;
     // 		  tmp[ka] = wc * wc * phi_c;
     // 		}
     // 	      add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
     // 			       wc * wc);
     // 	    }
     // 	}
    }

  // Compute coefficients of difference surface
  LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it2 = srf->basisFunctionsBegin();
  for (; it1 != cpsrf->basisFunctionsEnd(); ++it1, ++it2) 
    {
      auto nd_it = nom_denom.find(it2->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ka=0; ka<dim; ++ka)
	coef[ka] = (fabs(entry[dim]<tol)) ? 0 : entry[ka] / entry[dim];
      cpsrf->setCoef(coef, it1->second.get());
    }
 
  // Update initial surface
  double fac = 1.0; //1.01;
  srf->addSurface(*cpsrf, fac);
}


//==============================================================================
void LRSplineMBA::MBADistAndUpdate_omp(LRSplineSurface *srf)
//==============================================================================
{

  double tol = 1.0e-12;  // Numeric tolerance

  double umax = srf->endparam_u();
  double vmax = srf->endparam_v();
  int order2 = (srf->degree(XFIXED)+1)*(srf->degree(YFIXED)+1);

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
  map<const LRBSpline2D*, Array<double,2> > nom_denom; 


  vector<LRSplineSurface::ElementMap::const_iterator> el1_vec;
  const int num_elem = srf->numElements();
  el1_vec.reserve(num_elem);
  int max_num_bsplines = 0;
  for (LRSplineSurface::ElementMap::const_iterator iter = srf->elementsBegin(); iter != srf->elementsEnd(); ++iter)
  {
      el1_vec.push_back(iter);
      if (iter->second->nmbBasisFunctions() > max_num_bsplines)
      {
	  max_num_bsplines = iter->second->nmbBasisFunctions();
      }
  }
//  std::cout << "max_num_bsplines: " << max_num_bsplines << std::endl;
  vector<LRSplineSurface::ElementMap::const_iterator> el2_vec;
  el2_vec.reserve(cpsrf->numElements());
  for (LRSplineSurface::ElementMap::const_iterator iter = cpsrf->elementsBegin(); iter != cpsrf->elementsEnd(); ++iter)
  {
      el2_vec.push_back(iter);
  }

  int kdim = dim + 1;
  vector<double> elem_bspline_contributions(num_elem*max_num_bsplines*kdim, 0.0);


  // Traverse all elements. The two surfaces will have corresponding elements,
  // but only the source surface elements will contain point information so 
  // the elements in both surfaces must be traversed
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  LRSplineSurface::ElementMap::const_iterator el1;// = srf->elementsBegin();
  LRSplineSurface::ElementMap::const_iterator el2;// = cpsrf->elementsBegin();
  int kl, kk;
  // const int num_threads = 1;
  // omp_set_num_threads(num_threads);
#pragma omp parallel default(none) private(kl, kk, el1, el2) shared(nom_denom, tol, dim, el1_vec, el2_vec, umax, vmax, del, max_num_bsplines, elem_bspline_contributions, kdim, order2)
  {
      size_t nb;
      // Temporary vector to store weights associated with a given data point
      vector<double> tmp(dim);
      vector<double> tmp_weights;
      vector<double> ptval(dim);
      int nmb_pts;
      int ki, kr, ka;
      size_t kj;
      double *curr;
      vector<double> Bval;
      bool u_at_end, v_at_end;
      double val, dist, wc, wgt, phi_c, total_squared_inv;
#pragma omp for schedule(auto)//guided)//static,8)//runtime)//dynamic,4)
      for (kl = 0; kl < num_elem; ++kl)
      {
	  el1 = el1_vec[kl];
	  if (!el1->second->hasDataPoints())
	      continue;  // No points to use in surface update
	  el2 = el2_vec[kl];

	  // Fetch associated B-splines belonging to the difference surface
	  const vector<LRBSpline2D*>& bsplines = el1->second->getSupport();

	  // Check if the element needs to be updated
	  for (nb=0; nb<bsplines.size(); ++nb)
	      if (!bsplines[nb]->coefFixed())
		  break;

	  if (nb == bsplines.size())
	      continue;   // Element satisfies accuracy requirements

	  // Fetch points from the source surface
	  nmb_pts = el1->second->nmbDataPoints();
	  vector<double>& points = el1->second->getDataPoints();
//      int nmb_ghost = 0; //el1->second->nmbGhostPoints();
	  //vector<double>& ghost_points = el1->second->getGhostPoints();
	  //vector<double> ghost_points;

	  tmp_weights.resize(bsplines.size());
      
	  // Compute contribution from all points
	  // First compute distance in the data sets and store 
	  // basis function values
	  Bval.clear();
	  Bval.reserve(1.5*nmb_pts*order2);  // This vector is probably too large
	  for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
	  {
	      // Computing weights for this data point
	      u_at_end = (curr[0] > umax-tol) ? true : false;
	      v_at_end = (curr[1] > vmax-tol) ? true : false;
//	  double total_squared_inv = 0;
	      std::fill(ptval.begin(), ptval.end(), 0.0);
	      for (kj=0; kj<bsplines.size(); ++kj) 
	      {
		  val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
							u_at_end, v_at_end);
		  Bval.push_back(val);
		  const Point& tmp_pt = bsplines[kj]->coefTimesGamma();
		  for (ka=0; ka<dim; ++ka)
		      ptval[ka] += val*tmp_pt[ka];
	      }
	      dist;
	      if (dim == 1)
		  dist = curr[2] - ptval[0];
	      else
	      {
		  dist = Utils::distance_squared(ptval.begin(), ptval.end(),
						 points.begin()+ki*del+2); 
		  //ptval.dist(Point(curr+2, curr+del));
		  dist = sqrt(dist);
	      }
	      curr[del-1] = dist;
	  }

	  for (ki=0, kr=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
	  {
	      // Computing weights for this data point
	      total_squared_inv = 0;
	      for (kj=0; kj<bsplines.size(); ++kj, ++kr) 
	      {
		  val = Bval[kr];
		  wgt = val*bsplines[kj]->gamma();
		  tmp_weights[kj] = wgt;
		  total_squared_inv += wgt*wgt;
	      }
	      total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;

	      // Compute contribution
	      for (kj=0; kj<bsplines.size(); ++kj)
	      {
		  wc = tmp_weights[kj]; 
		  for (ka=0; ka<dim; ++ka)
		  {
		      phi_c = wc * curr[del-dim+ka] * total_squared_inv;
		      tmp[ka] = wc * wc * phi_c;
		  }
#if 0
#pragma omp critical
		  add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
				   wc * wc);
#else
		  for (kk = 0; kk < dim; ++kk)
		  {
		      elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + kk] += tmp[kk];
		  }
		  elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + dim] += wc*wc;
#endif
	      }
	  }
      }
  }

#if 1
  // We add the contributions sequentially.
  for (kl = 0; kl < num_elem; ++kl)
  {
      el1 = el1_vec[kl];
      const vector<LRBSpline2D*>& bsplines = el1->second->getSupport();
      int num_basis_funcs = bsplines.size();
      for (int ki = 0; ki < num_basis_funcs; ++ki)
      {
	  const LRBSpline2D* bspline = bsplines[ki];
	  auto it = nom_denom.find(bspline);
	  if (it != nom_denom.end()) 
	  {
	      // already in map
	      for (int kj=0; kj<kdim; ++kj)
		  it->second[kj] += elem_bspline_contributions[kl*max_num_bsplines*kdim + ki*kdim + kj];
	  } 
	  else 
	  {
	      double* tmp = &elem_bspline_contributions[kl*max_num_bsplines*kdim + ki*kdim];
	      nom_denom.insert({bspline, Array<double,2>{tmp[0], tmp[1]}});
	  }
      }
  }
#endif

  // Compute coefficients of difference surface
  LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it2 = srf->basisFunctionsBegin();
  for (; it1 != cpsrf->basisFunctionsEnd(); ++it1, ++it2) 
    {
      auto nd_it = nom_denom.find(it2->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ka=0; ka<dim; ++ka)
	coef[ka] = (fabs(entry[dim]<tol)) ? 0 : entry[ka] / entry[dim];
      cpsrf->setCoef(coef, it1->second.get());
    }

#if 0
  {
      std::cout << "Writing to file the srfs." << std::endl;
      std::ofstream fileout1("tmp/srf1.g2");
      srf->writeStandardHeader(fileout1);
      srf->write(fileout1);
      std::ofstream fileout2("tmp/srf2.g2");
      cpsrf->writeStandardHeader(fileout2);
      cpsrf->write(fileout2);
      double debug_val= 0.0;
  }
#endif
 
  // Update initial surface
  double fac = 1.0; //1.01;
  srf->addSurface(*cpsrf, fac);
}


//==============================================================================
void LRSplineMBA::MBAUpdate(LRSplineSurface *srf)
//==============================================================================
{
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
  map<const LRBSpline2D*, Array<double,2> > nom_denom; 

  // Temporary vector to store weights associated with a given data point
  vector<double> tmp_weights;  
  vector<double> tmp(dim);

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

      // Fetch points from the source surface
      int nmb_pts = el1->second->nmbDataPoints();
      vector<double>& points = el1->second->getDataPoints();
      int nmb_ghost = 0; //el1->second->nmbGhostPoints();
      //vector<double>& ghost_points = el1->second->getGhostPoints();
      vector<double> ghost_points;

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

      for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
      {
	  // printf("ki: %i\n", ki);
	  // printf("del: %i\n", del);
	  // printf("points.size(): %i\n", points.size());
	  // printf("curr[0]: %d\n", curr[0]);
	  // printf("curr[1]: %f\n", curr[1]);
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
	      }
	      add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
			       wc * wc);
	  }
	  // printf("Done with for loop.\n");
      }

      // Compute contribution from ghost points
      if (ghost_points.size() > 0)
      {
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
		  for (int ka=0; ka<dim; ++ka)
		  {
		      const double phi_c = wc * curr[del-dim+ka] * total_squared_inv;
		      tmp[ka] = wc * wc * phi_c;
		  }
		  add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
				   wc * wc);
	      }
	  }
      }
    }

  // Compute coefficients of difference surface
  LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it2 = srf->basisFunctionsBegin();
  for (; it1 != cpsrf->basisFunctionsEnd(); ++it1, ++it2) 
    {
      auto nd_it = nom_denom.find(it2->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ka=0; ka<dim; ++ka)
	coef[ka] = (fabs(entry[dim]<tol)) ? 0 : entry[ka] / entry[dim];
      cpsrf->setCoef(coef, it1->second.get());
    }
 
  // Update initial surface
  double fac = 1.0; //1.01;
  srf->addSurface(*cpsrf, fac);

 }


//==============================================================================
void LRSplineMBA::MBAUpdate_omp(LRSplineSurface *srf)
//==============================================================================
{
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
  map<const LRBSpline2D*, Array<double,2> > nom_denom; 

  vector<LRSplineSurface::ElementMap::const_iterator> el1_vec;
  const int num_elem = srf->numElements();
  el1_vec.reserve(num_elem);
  int max_num_bsplines = 0;
  for (LRSplineSurface::ElementMap::const_iterator iter = srf->elementsBegin(); iter != srf->elementsEnd(); ++iter)
  {
      el1_vec.push_back(iter);
      if (iter->second->nmbBasisFunctions() > max_num_bsplines)
      {
	  max_num_bsplines = iter->second->nmbBasisFunctions();
      }
  }
//  std::cout << "max_num_bsplines: " << max_num_bsplines << std::endl;
  vector<LRSplineSurface::ElementMap::const_iterator> el2_vec;
  el2_vec.reserve(cpsrf->numElements());
  for (LRSplineSurface::ElementMap::const_iterator iter = cpsrf->elementsBegin(); iter != cpsrf->elementsEnd(); ++iter)
  {
      el2_vec.push_back(iter);
  }

  int kdim = dim + 1;
  vector<double> elem_bspline_contributions(num_elem*max_num_bsplines*kdim, 0.0);

  // Traverse all elements. The two surfaces will have corresponding elements,
  // but only the source surface elements will contain point information so 
  // the elements in both surfaces must be traversed
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  LRSplineSurface::ElementMap::const_iterator el1;
  LRSplineSurface::ElementMap::const_iterator el2;
  int kl;
#pragma omp parallel default(none) private(kl, el1, el2) shared(nom_denom, tol, dim, el1_vec, el2_vec, umax, vmax, del, max_num_bsplines, elem_bspline_contributions, kdim)
  {
      vector<double> tmp(dim);
      // Temporary vector to store weights associated with a given data point
      vector<double> tmp_weights;  
      int nmb_pts, bsplines_size;
      int nmb_ghost;
      size_t nb;
      vector<double> ghost_points;
      bool u_at_end, v_at_end;
      int ki, kk, ka;
      size_t kj;
      const double *curr;
      double total_squared_inv, val, wgt, wc, phi_c, gamma;

#pragma omp for schedule(auto)//guided)//static,8)//runtime)//dynamic,4)
      for (kl = 0; kl < num_elem; ++kl)
      {
	  el1 = el1_vec[kl];
	  if (!el1->second->hasDataPoints())
	      continue;  // No points to use in surface update
	  el2 = el2_vec[kl];

	  // Fetch associated B-splines belonging to the difference surface
	  const vector<LRBSpline2D*>& bsplines = el2->second->getSupport();

	  bsplines_size = bsplines.size();

	  // Check if the element needs to be updated
	  for (nb=0; nb<bsplines.size(); ++nb)
	      if (!bsplines[nb]->coefFixed())
		  break;

	  if (nb == bsplines.size())
	      continue;   // Element satisfies accuracy requirements

	  // Fetch points from the source surface
	  nmb_pts = el1->second->nmbDataPoints();
	  vector<double>& points = el1->second->getDataPoints();
	  nmb_ghost = 0; //el1->second->nmbGhostPoints();
	  ghost_points.clear();
	  //vector<double>& ghost_points = el1->second->getGhostPoints();
	  // Compute contribution from all points
	  tmp_weights.resize(bsplines.size());
	  // std::cout << "tmp_weight.size(): " << tmp_weights.size() << std::endl;
	  // std::cout << "nmb_pts: " << nmb_pts << std::endl;
	  // std::cout << "points.size(): " << points.size() << std::endl;
	  // std::cout << "dim: " << dim << std::endl;
	  // std::cout << "del: " << del << std::endl;

	  for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
	  {
	      // printf("ki: %i\n", ki);
	      // printf("del: %i\n", del);
	      // printf("points.size(): %i\n", points.size());
	      // printf("curr[0]: %d\n", curr[0]);
	      // printf("curr[1]: %f\n", curr[1]);
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
		  }
#if 0
#pragma omp critical // Needed to make for loop thread safe. Slows down performance, should be circumvented.
		  add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
				   wc * wc);
#else
		  for (kk = 0; kk < dim; ++kk)
		  {
		      elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + kk] += tmp[kk];
		  }
		  elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + dim] += wc*wc;
#endif
	      }
	      // printf("Done with for loop.\n");
	  }

	  // Compute contribution from ghost points
	  if (ghost_points.size() > 0)
	  {
	      for (ki=0, curr=&ghost_points[0]; ki<nmb_ghost; ++ki, curr+=del)
	      {
		  // Computing weights for this data point
		  u_at_end = (curr[0] > umax-tol) ? true : false;
		  v_at_end = (curr[1] > vmax-tol) ? true : false;
		  total_squared_inv = 0;
		  for (kj=0; kj<bsplines.size(); ++kj) 
		  {
		      val = bsplines[kj]->evalBasisFunction(curr[0], curr[1], 0, 0,
							    u_at_end, v_at_end);
		      wgt = val*bsplines[kj]->gamma();
		      tmp_weights[kj] = wgt;
		      total_squared_inv += wgt*wgt;
		  }
		  total_squared_inv = (total_squared_inv < tol) ? 0.0 : 1.0/total_squared_inv;

		  // Compute contribution
		  for (kj=0; kj<bsplines.size(); ++kj)
		  {
		      wc = tmp_weights[kj]; 
		      for (ka=0; ka<dim; ++ka)
		      {
			  phi_c = wc * curr[del-dim+ka] * total_squared_inv;
			  tmp[ka] = wc * wc * phi_c;
		      }
#if 0
#pragma omp critical // Needed to make for loop thread safe. Slows down performance, should be circumvented.
		      add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
				       wc * wc);
#else
		      for (kk = 0; kk < dim; ++kk)
		      {
			  elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + kk] += tmp[kk];
		      }
		      elem_bspline_contributions[kl*max_num_bsplines*kdim + kj*kdim + dim] += wc*wc;
#endif
		  }
	      }
	  }
      }
  }

  // We add the contributions sequentially.
  for (kl = 0; kl < num_elem; ++kl)
  {
      el2 = el2_vec[kl];
      const vector<LRBSpline2D*>& bsplines = el2->second->getSupport();
      int num_basis_funcs = bsplines.size();
      for (int ki = 0; ki < num_basis_funcs; ++ki)
      {
	  const LRBSpline2D* bspline = bsplines[ki];
	  auto it = nom_denom.find(bspline);
	  if (it != nom_denom.end()) 
	  {
	      // already in map
	      for (int kj=0; kj<kdim; ++kj)
		  it->second[kj] += elem_bspline_contributions[kl*max_num_bsplines*kdim + ki*kdim + kj];
	  } 
	  else 
	  {
	      double* tmp = &elem_bspline_contributions[kl*max_num_bsplines*kdim + ki*kdim];
	      nom_denom.insert({bspline, Array<double,2>{tmp[0], tmp[1]}});
	  }
      }
  }
  

  // Compute coefficients of difference surface
  LRSplineSurface::BSplineMap::const_iterator it1 = cpsrf->basisFunctionsBegin();
  LRSplineSurface::BSplineMap::const_iterator it2 = srf->basisFunctionsBegin();
  for (; it1 != cpsrf->basisFunctionsEnd(); ++it1, ++it2) 
    {
      auto nd_it = nom_denom.find(it2->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ka=0; ka<dim; ++ka)
	coef[ka] = (fabs(entry[dim]<tol)) ? 0 : entry[ka] / entry[dim];
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


//==============================================================================
  void LRSplineMBA::MBAUpdate(LRSplineSurface *srf,
			      vector<Element2D*>& elems,
			      vector<Element2D*>& elems2)
//==============================================================================
{
  double tol = 1.0e-12;  // Numeric tolerance

  int dim = srf->dimension();
  double umax = srf->endparam_u();
  double vmax = srf->endparam_v();

  // Map to accumulate numerator and denominator to compute final coefficient value
  // for each BSplineFunction
  map<const LRBSpline2D*, Array<double,2> > nom_denom; 

  // Temporary vector to store weights associated with a given data point
  vector<double> tmp_weights;  
  vector<double> tmp(dim);

  // Collect all influenced element
  set<Element2D*> all_elems;
  elems2.clear();
  for (size_t ix1=0; ix1<elems.size(); ++ix1)
    {
      const vector<LRBSpline2D*>& bsplines = elems[ix1]->getSupport();
      for (size_t ix2=0; ix2<bsplines.size(); ++ix2)
	{
	  vector<Element2D*> el2 = bsplines[ix2]->supportedElements();
	  all_elems.insert(el2.begin(), el2.end());
	}
    }
  elems2.insert(elems2.end(), all_elems.begin(), all_elems.end());
  all_elems.clear();

  // Traverse all elements. The two surfaces will have corresponding elements,
  // but only the source surface elements will contain point information so 
  // the elements in both surfaces must be traversed
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  for (size_t ix_el=0; ix_el<elems2.size(); ++ix_el)
    {
      if (!elems2[ix_el]->hasDataPoints())
	continue;  // No points to use in surface update

      // Fetch associated B-splines belonging to the difference surface
      const vector<LRBSpline2D*>& bsplines = elems2[ix_el]->getSupport();

      // Fetch points from the source surface
      int nmb_pts = elems2[ix_el]->nmbDataPoints();
      vector<double>& points = elems2[ix_el]->getDataPoints();

       tmp_weights.resize(bsplines.size());
      
      // Compute contribution from all points
      int ki;
      size_t kj;
      double *curr;
      for (ki=0, curr=&points[0]; ki<nmb_pts; ++ki, curr+=del)
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
	      for (int ka=0; ka<dim; ++ka)
		{
		  const double phi_c = wc * curr[del-dim+ka] * total_squared_inv;
		  tmp[ka] = wc * wc * phi_c;
		}
	      add_contribution(dim, nom_denom, bsplines[kj], &tmp[0], 
			       wc * wc);
	    }
	}
     }

  // Compute coefficients of difference surface and update surface
  for (LRSplineSurface::BSplineMap::const_iterator it1 = srf->basisFunctionsBegin();
       it1 != srf->basisFunctionsEnd(); ++it1) 
    {
      auto nd_it = nom_denom.find(it1->second.get());
      const auto& entry = nd_it->second;
      Point coef(dim);
      for (int ka=0; ka<dim; ++ka)
	coef[ka] = (fabs(entry[dim]<tol)) ? 0 : entry[ka] / entry[dim];
      Point curr_coef = it1->second->Coef();
      srf->setCoef(curr_coef+coef, it1->second.get());
    }
 
}

//------------------------------------------------------------------------------
void LRSplineMBA::add_contribution(int dim, 
				   map<const LRBSpline2D*, Array<double,2> >& target, 
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
       target.insert({bspline, Array<double,2>{tmp[0], denom}});
     }
 }
