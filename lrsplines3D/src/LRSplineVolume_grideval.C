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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
//#include "GoTools/lrsplines3D/LRBSpline3D.h"
//#include "GoTools/lrsplines3D/LRBasisEvalResult.h"

using namespace Go;
using std::vector;
using std::pair;

  //============================================================================
  // Constructs a map from a B-spline to a simple Index.
  // ============================================================================
  LRSplineVolume::Bspline3DIndexMap construct_index_map( const LRSplineVolume& lrs )
  {
    LRSplineVolume::Bspline3DIndexMap result;
    size_t i = 0;
    for (LRSplineVolume::BSplineMap::const_iterator it_bs=lrs.basisFunctionsBegin(); 
        it_bs!=lrs.basisFunctionsEnd(); ++it_bs)
      {
	result[it_bs->second.get()] = i++;
      }
    return result;
  }

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       Dmatrix& basisValues,
				       Bspline3DIndexMap& bmap) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  basisValues.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    basisValues[ki].assign(num_coefs, 0.0);

  if ((int)bmap.size() != num_coefs)
    bmap = construct_index_map(*this);
  
  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  vector<pair<LRBSpline3D*, double> > bvals;
			  getBasisValues(param_u[kh], param_v[kr], param_w[km],
					 elem, bvals);
			  for (size_t kv=0; kv<bvals.size(); ++kv)
			    basisValues[(km*num_w+kr)*num_u+kh][bmap.at(bvals[kv].first)] =
			      bvals[kv].second;
			}
		    }
		}
	    }
	}
    }

}

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       vector<LRBasisPtsVol>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kv, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  kv = 0;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  computeBasis(param_u[kh], param_v[kr], param_w[km],
				       elem, result[kv++]);
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       Dmatrix& basisValues,
				       Dmatrix& basisDerivs_u,
				       Dmatrix& basisDerivs_v,
				       Dmatrix& basisDerivs_w,
				       Bspline3DIndexMap& bmap) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
 int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  basisValues.resize(num_par);
  basisDerivs_u.resize(num_par);
  basisDerivs_v.resize(num_par);
  basisDerivs_w.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    {
      basisValues[ki].assign(num_coefs, 0.0);
      basisDerivs_u[ki].assign(num_coefs, 0.0);
      basisDerivs_v[ki].assign(num_coefs, 0.0);
      basisDerivs_w[ki].assign(num_coefs, 0.0);
    }
  
  if ((int)bmap.size() != num_coefs)
    bmap = construct_index_map(*this);
  
  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  vector<pair<LRBSpline3D*, vector<double> > > bvals;
			  getBasisDerivs(param_u[kh], param_v[kr], param_w[km],
					 elem, bvals);
			  for (size_t kv=0; kv<bvals.size(); ++kv)
			    {
			      int ix = bmap.at(bvals[kv].first);
			      basisValues[(km*num_w+kr)*num_u+kh][ix] =
				bvals[kv].second[0];
			      basisDerivs_u[(km*num_w+kr)*num_u+kh][ix] =
				bvals[kv].second[1];
			      basisDerivs_v[(km*num_w+kr)*num_u+kh][ix] =
				bvals[kv].second[2];
			      basisDerivs_w[(km*num_w+kr)*num_u+kh][ix] =
				bvals[kv].second[3];
			    }
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       vector<LRBasisDerivsVol>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kv, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  kv = 0;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  computeBasis(param_u[kh], param_v[kr], param_w[km],
				       elem, result[kv++]);
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       vector<LRBasisDerivsVol2>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kv, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  kv = 0;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  computeBasis(param_u[kh], param_v[kr], param_w[km],
				       elem, result[kv++]);
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void LRSplineVolume::computeBasisGrid(const Dvector& param_u,
				       const Dvector& param_v,
				       const Dvector& param_w,
				       vector<LRBasisDerivsVol3>& result) const 
//===========================================================================
{
  int num_u = (int)param_u.size();
  int num_v = (int)param_v.size();
  int num_w = (int)param_w.size();
  int uorder = degree(XDIR) + 1;
  int vorder = degree(YDIR) + 1;
  int worder = degree(ZDIR) + 1;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  constructElementMesh(elements);
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  const double* const uknots_end = mesh_.knotsEnd(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);
  const double* knotu;
    
  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  const double* const vknots_end = mesh_.knotsEnd(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  const double* const wknots_end = mesh_.knotsEnd(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);
  const double* knotw;

  // Initiate to zero
  int ki, kj, kr, kh, kv, kk, km;
  int num_coefs = numBasisFunctions();
  int num_par = num_u*num_v*num_w;
  result.resize(num_par);

  // For all parameters
  double tolu = 1.0e-9;
  double tolv = 1.0e-9;
  double tolw = 1.0e-9;
  kv = 0;
  for (kk=0, km=0, knotw=wknots, ++knotw; knotw!=wknots_end; ++knotw, ++kk)
    {
      int lastw = (knotw+1 == wknots_end);
      for (; km<num_w && param_w[km] <= (*knotw)+lastw*tolw; ++km)
	{
	  for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
	    {
	      int lastv = (knotv+1 == vknots_end);
	      for (; kr<num_v && param_v[kr] <= (*knotv)+lastv*tolv; ++kr)
		{
		  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
		       ++knotu, ++ki)
		    {
		      int lastu = (knotu+1 == uknots_end);
		      Element3D *elem =
			elements[(kk*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];
		      for (; kh<num_u && param_u[kh] <= (*knotu)+lastu*tolu; ++kh)
			{
			  computeBasis(param_u[kh], param_v[kr], param_w[km],
				       elem, result[kv++]);
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, 
				LRBasisPtsVol& result) const
//===========================================================================
{
  Element3D *elem = coveringElement(upar, vpar, wpar);
  computeBasis(upar, vpar, wpar, elem, result);
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, Element3D* elem,
				LRBasisPtsVol& result) const
//===========================================================================
{
  const vector<LRBSpline3D*>& bspl = elem->getSupport();
  result.preparePts(upar, vpar, wpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<double> bvals(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      bvals[ki] = bspl[ki]->evalBasisFunc(upar, vpar, wpar);
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
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, 
				LRBasisDerivsVol& result) const
//===========================================================================
{
  Element3D *elem = coveringElement(upar, vpar, wpar);
  computeBasis(upar, vpar, wpar, elem, result);
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar,
			       Element3D* elem, LRBasisDerivsVol& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  double wmax = endparam_w();
  bool u_at_end = false, v_at_end = false, w_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  wpar = std::min(wpar, wmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  if (wmax - wpar < eps)
    w_at_end = true;
 
  const vector<LRBSpline3D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, wpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals(bspl.size());  // Position and first derivatives
  double wgt_sum0 = 1.0, wgt_sum1 = 1.0, wgt_sum2 = 1.0, wgt_sum3 = 1.0;
  if (rational())
    wgt_sum0 = wgt_sum1 = wgt_sum2 = wgt_sum3 = 0.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[2], der2[2], der3[2];
      bspl[ki]->getUnivariate(XDIR)->evalBasisFunctions(upar, 1, der1, u_at_end);
      bspl[ki]->getUnivariate(YDIR)->evalBasisFunctions(vpar, 1, der2, v_at_end);
      bspl[ki]->getUnivariate(ZDIR)->evalBasisFunctions(wpar, 1, der3, w_at_end);
      bvals[ki].resize(4);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0]*der3[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0]*der3[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1]*der3[0];
      bvals[ki][3] = wgt[ki]*der1[0]*der2[0]*der3[1];

      if (rational())
	{
	  wgt_sum0 += bvals[ki][0];
	  wgt_sum1 += bvals[ki][1];
	  wgt_sum2 += bvals[ki][2];
	  wgt_sum3 += bvals[ki][3];
	}

      for (int ka=0; ka<4; ++ka)
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
	  result.basisDerivs_w[ki] =
	    (wgt_sum0*bvals[ki][3] - wgt_sum3*bvals[ki][0])/sum2;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	  result.basisDerivs_w[ki] = bvals[ki][3];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, 
				LRBasisDerivsVol2& result) const
//===========================================================================
{
  Element3D *elem = coveringElement(upar, vpar, wpar);
  computeBasis(upar, vpar, wpar, elem, result);
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, 
			       Element3D* elem,	LRBasisDerivsVol2& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  double wmax = endparam_w();
  bool u_at_end = false, v_at_end = false, w_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  wpar = std::min(wpar, wmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  if (wmax - wpar < eps)
    w_at_end = true;
   
  const vector<LRBSpline3D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, wpar, elem, bspl.size());

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
      
      double der1[3], der2[3], der3[3];
      bspl[ki]->getUnivariate(XDIR)->evalBasisFunctions(upar, 2, der1, u_at_end);
      bspl[ki]->getUnivariate(YDIR)->evalBasisFunctions(vpar, 2, der2, v_at_end);
      bspl[ki]->getUnivariate(ZDIR)->evalBasisFunctions(wpar, 2, der3, w_at_end);
      bvals[ki].resize(10);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0]*der3[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0]*der3[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1]*der3[0];
      bvals[ki][3] = wgt[ki]*der1[0]*der2[0]*der3[1];
      bvals[ki][4] = wgt[ki]*der1[2]*der2[0]*der3[0];
      bvals[ki][5] = wgt[ki]*der1[1]*der2[1]*der3[0];
      bvals[ki][6] = wgt[ki]*der1[1]*der2[0]*der3[1];
      bvals[ki][7] = wgt[ki]*der1[0]*der2[2]*der3[0];
      bvals[ki][8] = wgt[ki]*der1[0]*der2[1]*der3[1];
      bvals[ki][9] = wgt[ki]*der1[0]*der2[0]*der3[2];
 
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
	  result.basisValues[ki] = bvals[ki][0]/wgt_sum[0];
	  double du = wgt_sum[0]*bvals[ki][1] - wgt_sum[1]*bvals[ki][0];
	  double dv = wgt_sum[0]*bvals[ki][2] - wgt_sum[2]*bvals[ki][0];
	  double dw = wgt_sum[0]*bvals[ki][3] - wgt_sum[3]*bvals[ki][0];
	  result.basisDerivs_u[ki] = du/sum2;
	  result.basisDerivs_v[ki] = dv/sum2;
	  result.basisDerivs_w[ki] = dw/sum2;
	  
	  double dhu = 2.0*wgt_sum[0]*wgt_sum[1];
	  double dhv = 2.0*wgt_sum[0]*wgt_sum[2];
	  double dhw = 2.0*wgt_sum[0]*wgt_sum[3];
	  double duu = wgt_sum[0]*bvals[ki][4] - wgt_sum[4]*bvals[ki][0];
	  double duv = wgt_sum[0]*bvals[ki][5] + wgt_sum[2]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][2] - wgt_sum[5]*bvals[ki][0];
	  double duw = wgt_sum[0]*bvals[ki][6] + wgt_sum[3]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][3] - wgt_sum[6]*bvals[ki][0];
	  double dvv = wgt_sum[0]*bvals[ki][7] - wgt_sum[7]*bvals[ki][0];
	  double dvw = wgt_sum[0]*bvals[ki][8] + wgt_sum[3]*bvals[ki][2] -
	    wgt_sum[2]*bvals[ki][3] - wgt_sum[8]*bvals[ki][0];
	  double dww = wgt_sum[0]*bvals[ki][9] - wgt_sum[9]*bvals[ki][0];
	  result.basisDerivs_uu[ki] = duu/sum2 - du*dhu/sum4;
	  result.basisDerivs_uv[ki] = duv/sum2 - du*dhv/sum4;
	  result.basisDerivs_uw[ki] = duw/sum2 - du*dhw/sum4;
	  result.basisDerivs_vv[ki] = dvv/sum2 - dv*dhv/sum4;
	  result.basisDerivs_vw[ki] = dvw/sum2 - dv*dhw/sum4;
	  result.basisDerivs_ww[ki] = dww/sum2 - dw*dhw/sum4;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	  result.basisDerivs_w[ki] = bvals[ki][3];
	  result.basisDerivs_uu[ki] = bvals[ki][4];
	  result.basisDerivs_uv[ki] = bvals[ki][5];
	  result.basisDerivs_uw[ki] = bvals[ki][6];
	  result.basisDerivs_vv[ki] = bvals[ki][7];
	  result.basisDerivs_vw[ki] = bvals[ki][8];
	  result.basisDerivs_ww[ki] = bvals[ki][9];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, 
				LRBasisDerivsVol3& result) const
//===========================================================================
{
  Element3D *elem = coveringElement(upar, vpar, wpar);
  computeBasis(upar, vpar, wpar, elem, result);
}

//===========================================================================
void
  LRSplineVolume::computeBasis(double upar, double vpar, double wpar, Element3D* elem,
				LRBasisDerivsVol3& result) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  double wmax = endparam_w();
  bool u_at_end = false, v_at_end = false, w_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  wpar = std::min(wpar, wmax);
  if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
   if (wmax - wpar < eps)
    w_at_end = true;
 
  const vector<LRBSpline3D*>& bspl = elem->getSupport();
  result.prepareDerivs(upar, vpar, wpar, elem, bspl.size());

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals(bspl.size());  // Position and first derivatives
  double wgt_sum[20];
  for (int ka=0; ka<20; ++ka)
    wgt_sum[ka] = (rational_) ? 0.0 : 1.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[4], der2[4], der3[4];
      bspl[ki]->getUnivariate(XDIR)->evalBasisFunctions(upar, 3, der1, u_at_end);
      bspl[ki]->getUnivariate(YDIR)->evalBasisFunctions(vpar, 3, der2, v_at_end);
      bspl[ki]->getUnivariate(ZDIR)->evalBasisFunctions(wpar, 3, der3, w_at_end);
      bvals[ki].resize(20);
      bvals[ki][0] = wgt[ki]*der1[0]*der2[0]*der3[0];
      bvals[ki][1] = wgt[ki]*der1[1]*der2[0]*der3[0];
      bvals[ki][2] = wgt[ki]*der1[0]*der2[1]*der3[0];
      bvals[ki][3] = wgt[ki]*der1[0]*der2[0]*der3[1];
      bvals[ki][4] = wgt[ki]*der1[2]*der2[0]*der3[0];
      bvals[ki][5] = wgt[ki]*der1[1]*der2[1]*der3[0];
      bvals[ki][6] = wgt[ki]*der1[1]*der2[0]*der3[1];
      bvals[ki][7] = wgt[ki]*der1[0]*der2[2]*der3[0];
      bvals[ki][8] = wgt[ki]*der1[0]*der2[1]*der3[1];
      bvals[ki][9] = wgt[ki]*der1[0]*der2[0]*der3[2];
      bvals[ki][10] = wgt[ki]*der1[3]*der2[0]*der3[0];
      bvals[ki][11] = wgt[ki]*der1[2]*der2[1]*der3[0];
      bvals[ki][12] = wgt[ki]*der1[2]*der2[0]*der3[1];
      bvals[ki][13] = wgt[ki]*der1[1]*der2[2]*der3[0];
      bvals[ki][14] = wgt[ki]*der1[1]*der2[1]*der3[1];
      bvals[ki][15] = wgt[ki]*der1[1]*der2[0]*der3[2];
      bvals[ki][16] = wgt[ki]*der1[0]*der2[3]*der3[0];
      bvals[ki][17] = wgt[ki]*der1[0]*der2[2]*der3[1];
      bvals[ki][18] = wgt[ki]*der1[0]*der2[1]*der3[2];
      bvals[ki][19] = wgt[ki]*der1[0]*der2[0]*der3[3];
 
      if (rational())
	{
	  for (int ka=0; ka<20; ++ka)
	    wgt_sum[ka] += bvals[ki][ka];
	}
      
      for (int ka=0; ka<20; ++ka)
	bvals[ki][ka] *= gamma[ki];
    }

  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      if (rational_)
	{
	  double sum2 = wgt_sum[0]*wgt_sum[0];
	  double sum3 = sum2*wgt_sum[0];
	  double sum4 = sum2*sum2;
	  result.basisValues[ki] = bvals[ki][0]/wgt_sum[0];
	  double du = wgt_sum[0]*bvals[ki][1] - wgt_sum[1]*bvals[ki][0];
	  double dv = wgt_sum[0]*bvals[ki][2] - wgt_sum[2]*bvals[ki][0];
	  double dw = wgt_sum[0]*bvals[ki][3] - wgt_sum[3]*bvals[ki][0];
	  result.basisDerivs_u[ki] = du/sum2;
	  result.basisDerivs_v[ki] = dv/sum2;
	  result.basisDerivs_w[ki] = dw/sum2;
	  
	  double dhu = 2.0*wgt_sum[0]*wgt_sum[1];
	  double dhv = 2.0*wgt_sum[0]*wgt_sum[2];
	  double dhw = 2.0*wgt_sum[0]*wgt_sum[3];
	  double duu = wgt_sum[0]*bvals[ki][4] - wgt_sum[4]*bvals[ki][0];
	  double duv = wgt_sum[0]*bvals[ki][5] + wgt_sum[2]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][2] - wgt_sum[5]*bvals[ki][0];
	  double dvu = wgt_sum[0]*bvals[ki][5] + wgt_sum[1]*bvals[ki][2] -
	    wgt_sum[2]*bvals[ki][1] - wgt_sum[5]*bvals[ki][0];
	  double duw = wgt_sum[0]*bvals[ki][6] + wgt_sum[3]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][3] - wgt_sum[6]*bvals[ki][0];
	  double dwu = wgt_sum[0]*bvals[ki][6] + wgt_sum[1]*bvals[ki][3] -
	    wgt_sum[3]*bvals[ki][1] - wgt_sum[6]*bvals[ki][0];
	  double dvv = wgt_sum[0]*bvals[ki][7] - wgt_sum[7]*bvals[ki][0];
	  double dvw = wgt_sum[0]*bvals[ki][8] + wgt_sum[3]*bvals[ki][2] -
	    wgt_sum[2]*bvals[ki][3] - wgt_sum[8]*bvals[ki][0];
	  double dwv = wgt_sum[0]*bvals[ki][8] + wgt_sum[2]*bvals[ki][3] -
	    wgt_sum[3]*bvals[ki][2] - wgt_sum[8]*bvals[ki][0];
	  double dww = wgt_sum[0]*bvals[ki][9] - wgt_sum[9]*bvals[ki][0];
	  result.basisDerivs_uu[ki] = duu/sum2 - du*dhu/sum4;
	  result.basisDerivs_uv[ki] = duv/sum2 - du*dhv/sum4;
	  result.basisDerivs_uw[ki] = duw/sum2 - du*dhw/sum4;
	  result.basisDerivs_vv[ki] = dvv/sum2 - dv*dhv/sum4;
	  result.basisDerivs_vw[ki] = dvw/sum2 - dv*dhw/sum4;
	  result.basisDerivs_ww[ki] = dww/sum2 - dw*dhw/sum4;

	  double duuu = wgt_sum[0]*bvals[ki][10] + wgt_sum[1]*bvals[ki][4] -
	    wgt_sum[4]*bvals[ki][1] - wgt_sum[10]*bvals[ki][0];
	  double duuv = wgt_sum[0]*bvals[ki][11] + wgt_sum[2]*bvals[ki][4] -
	    wgt_sum[4]*bvals[ki][2] - wgt_sum[11]*bvals[ki][0];
	  double duuw = wgt_sum[0]*bvals[ki][12] + wgt_sum[3]*bvals[ki][4] -
	    wgt_sum[4]*bvals[ki][3] - wgt_sum[12]*bvals[ki][0];
	  double duvv = wgt_sum[0]*bvals[ki][7] + wgt_sum[2]*bvals[ki][4] -
	    wgt_sum[4]*bvals[ki][2] - wgt_sum[7]*bvals[ki][0];
	  double duww = wgt_sum[0]*bvals[ki][10] + wgt_sum[3]*bvals[ki][4] -
	    wgt_sum[4]*bvals[ki][3] - wgt_sum[10]*bvals[ki][0];
	  double dvvu = wgt_sum[0]*bvals[ki][13] + wgt_sum[1]*bvals[ki][7] -
	    wgt_sum[7]*bvals[ki][1] - wgt_sum[13]*bvals[ki][0];
	  double dwwu = wgt_sum[0]*bvals[ki][15] + wgt_sum[1]*bvals[ki][9] -
	    wgt_sum[9]*bvals[ki][1] - wgt_sum[15]*bvals[ki][0];
	  double duvw = wgt_sum[0]*bvals[ki][14] + wgt_sum[3]*bvals[ki][5] +
	    wgt_sum[2]*bvals[ki][6] + wgt_sum[8]*bvals[ki][1] -
	    wgt_sum[1]*bvals[ki][8] - wgt_sum[6]*bvals[ki][2] -
	    wgt_sum[5]*bvals[ki][3] - wgt_sum[14]*bvals[ki][0];
	  double dvvv = wgt_sum[0]*bvals[ki][16] + wgt_sum[2]*bvals[ki][7] -
	    wgt_sum[7]*bvals[ki][2] - wgt_sum[16]*bvals[ki][0];
	  double dvvw = wgt_sum[0]*bvals[ki][17] + wgt_sum[3]*bvals[ki][7] -
	    wgt_sum[7]*bvals[ki][3] - wgt_sum[17]*bvals[ki][0];
	  double dvww = wgt_sum[0]*bvals[ki][18] + wgt_sum[3]*bvals[ki][8] -
	    wgt_sum[8]*bvals[ki][3] - wgt_sum[18]*bvals[ki][0];
	  double dwwv = wgt_sum[0]*bvals[ki][18] + wgt_sum[2]*bvals[ki][9] -
	    wgt_sum[9]*bvals[ki][2] - wgt_sum[18]*bvals[ki][0];
	  double dwww = wgt_sum[0]*bvals[ki][19] + wgt_sum[3]*bvals[ki][9] -
	    wgt_sum[9]*bvals[ki][3] - wgt_sum[19]*bvals[ki][0];
	  
	  result.basisDerivs_uuu[ki] = duuu/sum2 - 2.0*duu*dhu/sum4 -
	    2.0*du*wgt_sum[4]/sum3 + 6.0*du*wgt_sum[1]*wgt_sum[1]/sum4;
	  result.basisDerivs_vvv[ki] = dvvv/sum2 - 2.0*dvv*dhv/sum4 -
	    2.0*dv*wgt_sum[7]/sum3 + 6.0*dv*wgt_sum[2]*wgt_sum[2]/sum4;
	  result.basisDerivs_www[ki] = dwww/sum2 - 2.0*dww*dhw/sum4 -
	    2.0*dw*wgt_sum[9]/sum3 + 6.0*dw*wgt_sum[3]*wgt_sum[3]/sum4;
	  result.basisDerivs_uuv[ki] = duuv/sum2 - duu*dhv/sum4 -
	    duv*dhu/sum4 - 2.0*du*wgt_sum[5]/sum3 + 6.0*du*wgt_sum[1]*wgt_sum[2]/sum4;
	  result.basisDerivs_uvv[ki] = dvvu/sum2 - dvv*dhu/sum4 -
	    dvu*dhv/sum4 - 2.0*dv*wgt_sum[5]/sum3 + 6.0*dv*wgt_sum[1]*wgt_sum[2]/sum4;
	  result.basisDerivs_uuw[ki] = duuw/sum2 - duu*dhw/sum4 -
	    duw*dhu/sum4 - 2.0*du*wgt_sum[6]/sum3 + 6.0*du*wgt_sum[1]*wgt_sum[3]/sum4;
	  result.basisDerivs_uww[ki] = dwwu/sum2 - dww*dhu/sum4 -
	    dwu*dhw/sum4 - 2.0*dw*wgt_sum[6]/sum3 + 6.0*dw*wgt_sum[1]*wgt_sum[3]/sum4;
	  result.basisDerivs_vvw[ki] = dvvw/sum2 - dvv*dhw/sum4 -
	    dvw*dhv/sum4 - 2.0*dv*wgt_sum[8]/sum3 + 6.0*dv*wgt_sum[2]*wgt_sum[3]/sum4;
	  result.basisDerivs_vww[ki] = dwwv/sum2 - dww*dhv/sum4 -
	    dwv*dhw/sum4 - 2.0*dw*wgt_sum[8]/sum3 + 6.0*dw*wgt_sum[2]*wgt_sum[3]/sum4;
	   result.basisDerivs_uvw[ki] = duvw/sum2 - duv*dhw/sum4 - duw*dhv/sum4 -
	     2.0*du*wgt_sum[8]/sum3 + 6.0*du*wgt_sum[2]*wgt_sum[3]/sum4;
	}
      else
	{
	  result.basisValues[ki] = bvals[ki][0];
	  result.basisDerivs_u[ki] = bvals[ki][1];
	  result.basisDerivs_v[ki] = bvals[ki][2];
	  result.basisDerivs_w[ki] = bvals[ki][3];
	  result.basisDerivs_uu[ki] = bvals[ki][4];
	  result.basisDerivs_uv[ki] = bvals[ki][5];
	  result.basisDerivs_uw[ki] = bvals[ki][6];
	  result.basisDerivs_vv[ki] = bvals[ki][7];
	  result.basisDerivs_vw[ki] = bvals[ki][8];
	  result.basisDerivs_ww[ki] = bvals[ki][9];
	  result.basisDerivs_uuu[ki] = bvals[ki][10];
	  result.basisDerivs_uuv[ki] = bvals[ki][11];
	  result.basisDerivs_uuw[ki] = bvals[ki][12];
	  result.basisDerivs_uvv[ki] = bvals[ki][13];
	  result.basisDerivs_uvw[ki] = bvals[ki][14];
	  result.basisDerivs_uww[ki] = bvals[ki][15];
	  result.basisDerivs_vvv[ki] = bvals[ki][16];
	  result.basisDerivs_vvw[ki] = bvals[ki][17];
	  result.basisDerivs_vww[ki] = bvals[ki][18];
	  result.basisDerivs_www[ki] = bvals[ki][19];
	}
      result.bsplines[ki] = bspl[ki];
    }
}

//===========================================================================
void
  LRSplineVolume::getBasisValues(double upar, double vpar, double wpar, Element3D* elem,
				  vector<pair<LRBSpline3D*, double> >& bvals) const
//===========================================================================
{
  const vector<LRBSpline3D*>& bspl = elem->getSupport();

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<double> bvals2(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      bvals2[ki] = bspl[ki]->evalBasisFunc(upar, vpar, wpar);
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
LRSplineVolume::getBasisDerivs(double upar, double vpar, double wpar, Element3D* elem,
				vector<pair<LRBSpline3D*, vector<double> > >& bvals) const
//===========================================================================
{
  double eps = 1.0e-12;
  double umax = endparam_u();
  double vmax = endparam_v();
  double wmax = endparam_w();
  bool u_at_end = false, v_at_end = false, w_at_end = false;
  upar = std::min(upar, umax);
  vpar = std::min(vpar, vmax);
  wpar = std::min(wpar, wmax);
 if (umax - upar < eps)
    u_at_end = true;
  if (vmax - vpar < eps)
    v_at_end = true;
  if (wmax - wpar < eps)
    w_at_end = true;

  const vector<LRBSpline3D*>& bspl = elem->getSupport();

  vector<double> gamma(bspl.size(), 1.0);
  vector<double> wgt(bspl.size(), 1.0);
  vector<vector<double> > bvals2(bspl.size());  // Position and first derivatives
  double wgt_sum0 = 1.0, wgt_sum1 = 1.0, wgt_sum2 = 1.0, wgt_sum3 = 1.0;
  if (rational())
    wgt_sum0 = wgt_sum1 = wgt_sum2 = wgt_sum3 = 0.0;
  
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      wgt[ki] = bspl[ki]->weight();
      gamma[ki] = bspl[ki]->gamma();
      
      double der1[2], der2[2], der3[2];
      bspl[ki]->getUnivariate(XDIR)->evalBasisFunctions(upar, 1, der1, u_at_end);
      bspl[ki]->getUnivariate(YDIR)->evalBasisFunctions(vpar, 1, der2, v_at_end);
      bspl[ki]->getUnivariate(ZDIR)->evalBasisFunctions(wpar, 1, der3, w_at_end);
      bvals2[ki].resize(4);
      bvals2[ki][0] = wgt[ki]*der1[0]*der2[0]*der3[0];
      bvals2[ki][1] = wgt[ki]*der1[1]*der2[0]*der3[0];
      bvals2[ki][2] = wgt[ki]*der1[0]*der2[1]*der3[0];
      bvals2[ki][3] = wgt[ki]*der1[0]*der2[0]*der3[1];
 
      if (rational())
	{
	  wgt_sum0 += bvals2[ki][0];
	  wgt_sum1 += bvals2[ki][1];
	  wgt_sum2 += bvals2[ki][2];
	  wgt_sum3 += bvals2[ki][3];
	}
      
      for (int ka=0; ka<4; ++ka)
	bvals2[ki][ka] *= gamma[ki];
   }

  bvals.resize(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    {
      vector<double> bder(4);
      if (rational_)
	{
	  double sum2 = wgt_sum0*wgt_sum0;
	  bder[0] = bvals2[ki][0]/wgt_sum0;
	  bder[1] = (wgt_sum0*bvals2[ki][1] - wgt_sum1*bvals2[ki][0])/sum2;
	  bder[2] = (wgt_sum0*bvals2[ki][2] - wgt_sum2*bvals2[ki][0])/sum2;
	  bder[3] = (wgt_sum0*bvals2[ki][3] - wgt_sum3*bvals2[ki][0])/sum2;
	}
      else
	{
	  bder[0] = bvals2[ki][0];
	  bder[1] = bvals2[ki][1];
	  bder[2] = bvals2[ki][2];
	  bder[3] = bvals2[ki][3];
	}
      bvals[ki] = std::make_pair(bspl[ki], bder);
    }
}
