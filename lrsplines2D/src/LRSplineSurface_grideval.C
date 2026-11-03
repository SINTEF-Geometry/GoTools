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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
//#include "GoTools/lrsplines2D/LRBSpline2D.h"
//#include "GoTools/lrsplines2D/LRBasisEvalResult.h"

using namespace Go;
using std::vector;
using std::pair;

  //============================================================================
  // Constructs a map from a B-spline to a simple Index.
  // ============================================================================
  LRSplineSurface::Bspline2DIndexMap construct_index_map( const LRSplineSurface& lrs )
  {
    LRSplineSurface::Bspline2DIndexMap result;
    size_t i = 0;
    for (LRSplineSurface::BSplineMap::const_iterator it_bs=lrs.basisFunctionsBegin(); 
        it_bs!=lrs.basisFunctionsEnd(); ++it_bs)
      {
	result[it_bs->second.get()] = i++;
      }
    return result;
  }

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       Dmatrix& basisValues,
				       Bspline2DIndexMap& bmap) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  basisValues.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    basisValues[ki].assign(num_coefs, 0.0);

  if ((int)bmap.size() != num_coefs)
    bmap = construct_index_map(*this);
  
  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  vector<pair<LRBSpline2D*, double> > bvals;
		  getBasisValues(param_u[kh], param_v[kr], elem, bvals);
		  for (size_t kv=0; kv<bvals.size(); ++kv)
		    basisValues[kr*num_u+kh][bmap.at(bvals[kv].first)] =
		      bvals[kv].second;
		}
	    }
	}
    }

}

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       vector<LRBasisPtsSf>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh, kv;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  kv = 0;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  computeBasis(param_u[kh], param_v[kr], elem, result[kv++]);
		}
	    }
	}
    }

}

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       Dmatrix& basisValues,
				       Dmatrix& basisDerivs_u,
				       Dmatrix& basisDerivs_v,
				       Bspline2DIndexMap& bmap) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  basisValues.resize(num_par);
  basisDerivs_u.resize(num_par);
  basisDerivs_v.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    {
      basisValues[ki].assign(num_coefs, 0.0);
      basisDerivs_u[ki].assign(num_coefs, 0.0);
      basisDerivs_v[ki].assign(num_coefs, 0.0);
    }
  
  if ((int)bmap.size() != num_coefs)
    bmap = construct_index_map(*this);
  
  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  vector<pair<LRBSpline2D*, vector<double> > > bvals;
		  getBasisDerivs(param_u[kh], param_v[kr], elem, bvals);
		  for (size_t kv=0; kv<bvals.size(); ++kv)
		    {
		      int ix = bmap.at(bvals[kv].first);
		      basisValues[kr*num_u+kh][ix] = bvals[kv].second[0];
		      basisDerivs_u[kr*num_u+kh][ix] = bvals[kv].second[1];
		      basisDerivs_v[kr*num_u+kh][ix] = bvals[kv].second[2];
		    }
		}
	    }
	}
    }
}

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       vector<LRBasisDerivsSf>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh, kv;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  kv = 0;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  computeBasis(param_u[kh], param_v[kr], elem, result[kv++]);
		}
	    }
	}
    }

}

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       vector<LRBasisDerivsSf2>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh, kv;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  kv = 0;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  computeBasis(param_u[kh], param_v[kr], elem, result[kv++]);
		}
	    }
	}
    }

}

//===========================================================================
void LRSplineSurface::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       vector<LRBasisDerivsSf3>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int uorder = degree(XFIXED) + 1;
  int vorder = degree(YFIXED) + 1;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  const double* const uknots_end = mesh_.knotsEnd(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  const double* const vknots_end = mesh_.knotsEnd(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
  const double* knotv;

  // Initiate to zero
  int ki, kj, kr, kh, kv;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  kv = 0;
  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
	{
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
		{
		  computeBasis(param_u[kh], param_v[kr], elem, result[kv++]);
		}
	    }
	}
    }

}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, 
				LRBasisPtsSf& result) const
//===========================================================================
{
  Element2D *elem = coveringElement(upar, vpar);
  computeBasis(upar, vpar, elem, result);
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, Element2D* elem,
				LRBasisPtsSf& result) const
//===========================================================================
{
  const vector<LRBSpline2D*>& bspl = elem->getSupport();
  result.preparePts(upar, vpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<double> bvals(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      bvals[ki] = bspl[ki]->evalBasisFunc(upar, vpar);
    }
  
  double wgt_sum = 1.0;
  if (rational())
    {
      wgt_sum = 0.0;
      for (size_t ki=0; ki<bspl.size(); ++ki)
	wgt_sum += wgt[ki]*bvals[ki];
    }

  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      double val = (rational_) ? wgt[ki]*gamma[ki]*bvals[ki]/wgt_sum :
	gamma[ki]*bvals[ki];
      result.basisValues[ki] = val;
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, 
				LRBasisDerivsSf& result) const
//===========================================================================
{
  Element2D *elem = coveringElement(upar, vpar);
  computeBasis(upar, vpar, elem, result);
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, Element2D* elem,
				LRBasisDerivsSf& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  bool u_at_end = false, v_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  
  const vector<LRBSpline2D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals(bspl.size());  // Position and first derivatives
  double wgt_sum0 = 1.0, wgt_sum1 = 1.0, wgt_sum2 = 1.0;
  if (rational())
    wgt_sum0 = wgt_sum1 = wgt_sum2 = 0.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[2], der2[2];
      bspl[ki]->getUnivariate(XFIXED)->evalBasisFunctions(upar, 1, der1, u_at_end);
      bspl[ki]->getUnivariate(YFIXED)->evalBasisFunctions(vpar, 1, der2, v_at_end);
      bvals[ki].resize(3);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1];
  
      if (rational())
	{
	  wgt_sum0 += bvals[ki][0];
	  wgt_sum1 += bvals[ki][1];
	  wgt_sum2 += bvals[ki][2];
	}
      
       for (int ka=0; ka<3; ++ka)
	bvals[ki][ka] *= gamma[ki];
   }

  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      if (rational_)
	{
	  double sum2 = wgt_sum0*wgt_sum0;
	  result.basisValues[ki] = bvals[ki][0]/wgt_sum0;
	  result.basisDerivs_u[ki] =
	    (wgt_sum0*bvals[ki][1] - wgt_sum1*bvals[ki][0])/sum2;
	  result.basisDerivs_v[ki] =
	    (wgt_sum0*bvals[ki][2] - wgt_sum2*bvals[ki][0])/sum2;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, 
				LRBasisDerivsSf2& result) const
//===========================================================================
{
  Element2D *elem = coveringElement(upar, vpar);
  computeBasis(upar, vpar, elem, result);
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, Element2D* elem,
				LRBasisDerivsSf2& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  bool u_at_end = false, v_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  
  const vector<LRBSpline2D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals(bspl.size());  // Position and first derivatives
  double wgt_sum[6];
  for (int ka=0; ka<6; ++ka)
    wgt_sum[ka] = (rational_) ? 0.0 : 1.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[3], der2[3];
      bspl[ki]->getUnivariate(XFIXED)->evalBasisFunctions(upar, 2, der1, u_at_end);
      bspl[ki]->getUnivariate(YFIXED)->evalBasisFunctions(vpar, 2, der2, v_at_end);
      bvals[ki].resize(6);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1];
      bvals[ki][3] = wgt[ki]*der1[2]*der2[0];
      bvals[ki][4] = wgt[ki]*der1[1]*der2[1];
      bvals[ki][5] = wgt[ki]*der1[0]*der2[2];
 
      if (rational())
	{
	  for (int ka=0; ka<6; ++ka)
	    wgt_sum[ka] += bvals[ki][ka];
	}
      
      for (int ka=0; ka<6; ++ka)
	bvals[ki][ka] *= gamma[ki];
    }

  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      if (rational_)
	{
	  double sum2 = wgt_sum[0]*wgt_sum[0];
	  double sum4 = sum2*sum2;
	  result.basisValues[ki] = bvals[ki][0]/wgt_sum[0];
	  double du = wgt_sum[0]*bvals[ki][1] - wgt_sum[1]*bvals[ki][0];
	  double dv = wgt_sum[0]*bvals[ki][2] - wgt_sum[2]*bvals[ki][0];
	  result.basisDerivs_u[ki] = du/sum2;
	  result.basisDerivs_v[ki] = dv/sum2;
	  double dwu = 2.0*wgt_sum[0]*wgt_sum[1];
	  double dwv = 2.0*wgt_sum[0]*wgt_sum[2];
	  double duu = wgt_sum[0]*bvals[ki][3] - wgt_sum[3]*bvals[ki][0];
	  double duv = wgt_sum[0]*bvals[ki][4] + wgt_sum[2]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][2] - wgt_sum[4]*bvals[ki][0];
	  double dvv = wgt_sum[0]*bvals[ki][5] - wgt_sum[5]*bvals[ki][0];
	  result.basisDerivs_uu[ki] = duu/sum2 - du*dwu/sum4;
	  result.basisDerivs_uv[ki] = duv/sum2 - du*dwv/sum4;
	  result.basisDerivs_vv[ki] = dvv/sum2 - dv*dwv/sum4;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	  result.basisDerivs_uu[ki] = bvals[ki][3];
	  result.basisDerivs_uv[ki] = bvals[ki][4];
	  result.basisDerivs_vv[ki] = bvals[ki][5];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, 
				LRBasisDerivsSf3& result) const
//===========================================================================
{
  Element2D *elem = coveringElement(upar, vpar);
  computeBasis(upar, vpar, elem, result);
}

//===========================================================================
void
  LRSplineSurface::computeBasis(double upar, double vpar, Element2D* elem,
				LRBasisDerivsSf3& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  bool u_at_end = false, v_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  
  const vector<LRBSpline2D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals(bspl.size());  // Position and first derivatives
  double wgt_sum[10];
  for (int ka=0; ka<10; ++ka)
    wgt_sum[ka] = (rational_) ? 0.0 : 1.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[4], der2[4];
      bspl[ki]->getUnivariate(XFIXED)->evalBasisFunctions(upar, 3, der1, u_at_end);
      bspl[ki]->getUnivariate(YFIXED)->evalBasisFunctions(vpar, 3, der2, v_at_end);
      bvals[ki].resize(10);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1];
      bvals[ki][3] = wgt[ki]*der1[2]*der2[0];
      bvals[ki][4] = wgt[ki]*der1[1]*der2[1];
      bvals[ki][5] = wgt[ki]*der1[0]*der2[2];
      bvals[ki][6] = wgt[ki]*der1[3]*der2[0];
      bvals[ki][7] = wgt[ki]*der1[2]*der2[1];
      bvals[ki][8] = wgt[ki]*der1[1]*der2[2];
      bvals[ki][9] = wgt[ki]*der1[0]*der2[3];
 
      if (rational())
	{
	  for (int ka=0; ka<10; ++ka)
	    wgt_sum[ka] += bvals[ki][ka];
	}
      
      for (int ka=0; ka<10; ++ka)
	bvals[ki][ka] *= gamma[ki];
    }

  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      if (rational_)
	{
	  double sum2 = wgt_sum[0]*wgt_sum[0];
	  double sum4 = sum2*sum2;
	  double sum3 = wgt_sum[0]*sum2;
	  result.basisValues[ki] = bvals[ki][0]/wgt_sum[0];
	  
	  double du = wgt_sum[0]*bvals[ki][1] - wgt_sum[1]*bvals[ki][0];
	  double dv = wgt_sum[0]*bvals[ki][2] - wgt_sum[2]*bvals[ki][0];
	  result.basisDerivs_u[ki] = du/sum2;
	  result.basisDerivs_v[ki] = dv/sum2;
	  
	  double dwu = 2.0*wgt_sum[0]*wgt_sum[1];
	  double dwv = 2.0*wgt_sum[0]*wgt_sum[2];
	  double duu = wgt_sum[0]*bvals[ki][3] - wgt_sum[3]*bvals[ki][0];
	  double duv = wgt_sum[0]*bvals[ki][4] + wgt_sum[2]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][2] - wgt_sum[4]*bvals[ki][0];
	  double dvu = wgt_sum[0]*bvals[ki][4] + wgt_sum[1]*bvals[ki][2] -
	    wgt_sum[2]*bvals[ki][1] - wgt_sum[4]*bvals[ki][0];
	  double dvv = wgt_sum[0]*bvals[ki][5] - wgt_sum[5]*bvals[ki][0];
	  result.basisDerivs_uu[ki] = duu/sum2 - du*dwu/sum4;
	  result.basisDerivs_uv[ki] = duv/sum2 - du*dwv/sum4;
	  result.basisDerivs_vv[ki] = dvv/sum2 - dv*dwv/sum4;

	  double duuu = wgt_sum[0]*bvals[ki][6] + wgt_sum[1]*bvals[ki][3] -
	    wgt_sum[3]*bvals[ki][1] - wgt_sum[6]*bvals[ki][0];
	  double duuv = wgt_sum[0]*bvals[ki][7] + wgt_sum[2]*bvals[ki][3] -
	    wgt_sum[3]*bvals[ki][2] - wgt_sum[7]*bvals[ki][0];
	  double duvv = wgt_sum[0]*bvals[ki][7] + wgt_sum[2]*bvals[ki][3] -
	    wgt_sum[3]*bvals[ki][2] - wgt_sum[7]*bvals[ki][0];
	  double dvvu = wgt_sum[0]*bvals[ki][8] + wgt_sum[1]*bvals[ki][5] -
	    wgt_sum[5]*bvals[ki][1] - wgt_sum[8]*bvals[ki][0];
	  double dvvv = wgt_sum[0]*bvals[ki][9] + wgt_sum[2]*bvals[ki][5] -
	    wgt_sum[5]*bvals[ki][2] - wgt_sum[9]*bvals[ki][0];
	  result.basisDerivs_uuu[ki] = duuu/sum2 - 2.0*duu*dwu/sum4 -
	    2.0*du*wgt_sum[3]/sum3 + 6.0*du*wgt_sum[1]*wgt_sum[1]/sum4;
	  result.basisDerivs_vvv[ki] = dvvv/sum2 - 2.0*dvv*dwv/sum4 -
	    2.0*dv*wgt_sum[5]/sum3 + 6.0*dv*wgt_sum[2]*wgt_sum[2]/sum4;
	  result.basisDerivs_uuv[ki] = duuv/sum2 - duu*dwv/sum4 -
	    duv*dwu/sum4 - 2.0*du*wgt_sum[4]/sum3 + 6.0*du*wgt_sum[1]*wgt_sum[2]/sum4;
	  result.basisDerivs_uvv[ki] = dvvu/sum2 - dvv*dwu/sum4 -
	    dvu*dwv/sum4 - 2.0*dv*wgt_sum[4]/sum3 + 6.0*dv*wgt_sum[1]*wgt_sum[2]/sum4;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	  result.basisDerivs_uu[ki] = bvals[ki][3];
	  result.basisDerivs_uv[ki] = bvals[ki][4];
	  result.basisDerivs_vv[ki] = bvals[ki][5];
	  result.basisDerivs_uuu[ki] = bvals[ki][6];
	  result.basisDerivs_uuv[ki] = bvals[ki][7];
	  result.basisDerivs_uvv[ki] = bvals[ki][8];
	  result.basisDerivs_vvv[ki] = bvals[ki][9];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineSurface::getBasisValues(double upar, double vpar, Element2D* elem,
				  vector<pair<LRBSpline2D*, double> >& bvals) const
//===========================================================================
{
  const vector<LRBSpline2D*>& bspl = elem->getSupport();

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<double> bvals2(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      bvals2[ki] = bspl[ki]->evalBasisFunc(upar, vpar);
    }
  
  double wgt_sum = 1.0;
  if (rational())
    {
      wgt_sum = 0.0;
      for (size_t ki=0; ki<bspl.size(); ++ki)
	wgt_sum += (wgt[ki]*bvals2[ki]);
    }

  bvals.resize(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      double val = (rational_) ? wgt[ki]*gamma[ki]*bvals2[ki]/wgt_sum :
	gamma[ki]*bvals2[ki];
      bvals[ki] = std::make_pair(bspl[ki], val);
    }
}

//===========================================================================
void
LRSplineSurface::getBasisDerivs(double upar, double vpar, Element2D* elem,
				vector<pair<LRBSpline2D*, vector<double> > >& bvals) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  bool u_at_end = false, v_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  
  const vector<LRBSpline2D*>& bspl = elem->getSupport();

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals2(bspl.size());  // Position and first derivatives
  double wgt_sum0 = 1.0, wgt_sum1 = 1.0, wgt_sum2 = 1.0;
  if (rational())
    wgt_sum0 = wgt_sum1 = wgt_sum2 = 0.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[2], der2[2];
      bspl[ki]->getUnivariate(XFIXED)->evalBasisFunctions(upar, 1, der1, u_at_end);
      bspl[ki]->getUnivariate(YFIXED)->evalBasisFunctions(vpar, 1, der2, v_at_end);
      bvals2[ki].resize(3);
      bvals2[ki][0] = der1[0]*der2[0];
      bvals2[ki][1] = der1[1]*der2[0];
      bvals2[ki][2] = der1[0]*der2[1];
  
      if (rational())
	{
	  wgt_sum0 += (wgt[ki]*bvals2[ki][0]);
	  wgt_sum1 += (wgt[ki]*bvals2[ki][1]);
	  wgt_sum2 += (wgt[ki]*bvals2[ki][2]);
	}
    }

  bvals.resize(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      vector<double> bder(3);
      if (rational_)
	{
	  double sum2 = wgt_sum0*wgt_sum0;
	  bder[0] = wgt[ki]*gamma[ki]*bvals2[ki][0]/wgt_sum0;
	  bder[1] = (wgt_sum0*wgt[ki]*gamma[ki]*bvals2[ki][1] -
		     wgt_sum1*wgt[ki]*gamma[ki]*bvals2[ki][0])/sum2;
	  bder[2] = (wgt_sum0*wgt[ki]*gamma[ki]*bvals2[ki][2] -
		     wgt_sum2*wgt[ki]*gamma[ki]*bvals2[ki][0])/sum2;
	}
      else
	{
	  bder[0] = gamma[ki]*bvals2[ki][0];
	  bder[1] = gamma[ki]*bvals2[ki][1];
	  bder[2] = gamma[ki]*bvals2[ki][2];
	}
      bvals[ki] = std::make_pair(bspl[ki], bder);
    }
}
