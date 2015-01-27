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

using std::vector;
using std::set;
using std::map;
using std::array;
using namespace Go;

//==============================================================================
void LRSplineMBA::MBADistAndUpdate(LRSplineSurface *srf)
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

