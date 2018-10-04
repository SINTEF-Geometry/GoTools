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

#include "GoTools/creators/ModifySurf.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothSurfSet.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;
using namespace std;

//===========================================================================
void ModifySurf::smoothSurface(shared_ptr<SplineSurface>& surf,
			       int nmb_bd_keep)
//===========================================================================
{
  nmb_bd_keep = std::min(nmb_bd_keep, 2);
  nmb_bd_keep = std::max(nmb_bd_keep, 1);

  // Define fixed coefficients
  int kn1 = surf->numCoefs_u();
  int kn2 = surf->numCoefs_v();
  vector<int> coef_known(kn1*kn2, 0);  // All coefficients are free

  // Fix edges
  int kr, kj;
  for (kr = 0; kr < nmb_bd_keep; ++kr)
    for (kj = 0; kj < kn1; ++kj)
      coef_known[kj+kn1*kr] = 1;

  for (kr = 0; kr < nmb_bd_keep; ++kr)
    for (kj = 0; kj < kn2; ++kj)
      coef_known[(kj+1)*kn1-1-kr] = 1;

  for (kr = 0; kr < nmb_bd_keep; ++kr)
    for (kj = 0; kj < kn1; ++kj)
      coef_known[(kn2-1-kr)*kn1+kj] = 1;

  for (kr = 0; kr < nmb_bd_keep; ++kr)
    for (kj = 0; kj < kn2; ++kj)
      coef_known[kj*kn1+kr] = 1;
  
  // Perform smoothing 
  // First set weights
  double w1 = 0.0;
  double w2 = 0.5;
  double w3 = 0.5;

  SmoothSurf smooth;
  int seem[2];
  seem[0] = seem[1] = 0;
  smooth.attach(surf, seem, &coef_known[0]);

  smooth.setOptimize(w1, w2, w3);

  int stat;
  try {
    stat = smooth.equationSolve(surf);
  }
  catch (...)
    {
      THROW("Surface smoothing failed");
    }
}

//===========================================================================
void ModifySurf::replaceBoundary(shared_ptr<SplineSurface> surf,
				 shared_ptr<SplineCurve> curve,
				 int bd_idx, double tol)
//===========================================================================
{
  // Sample points in the initial surface
  RectDomain dom = surf->parameterDomain();
  int dim = surf->dimension();

  vector<double> sample_pts0;
  vector<double> sample_pts;
  vector<double> sample_pars;
  vector<double> par_u, par_v;
  int nmb_sample = 20;
  surf->gridEvaluator(nmb_sample, nmb_sample, sample_pts0, par_u, par_v);

  // Replace boundary curve
  surf->replaceBoundaryCurve(bd_idx, curve);

  // Collect approxation points
  sample_pars.reserve(2*nmb_sample*nmb_sample);
  sample_pts.reserve(sample_pts0.size());
  int ki, kj, kr;
  for (kj=2, kr=kj*nmb_sample; kj<nmb_sample-2; ++kj, kr=kj*nmb_sample)
    for (ki=2, kr+=2; ki<nmb_sample-2; ++ki, ++kr)
      {
	Point curr(sample_pts0.begin()+kr*dim, sample_pts0.begin()+(kr+1)*dim);
	double u1, v1, d1;
	Point pnt2;
	double seed[2];
	seed[0] = par_u[ki];
	seed[1] = par_v[kj];
	surf->closestPoint(curr, u1, v1, pnt2, d1, tol, NULL, seed);

	sample_pars.push_back(u1);
	sample_pars.push_back(v1);
	sample_pts.insert(sample_pts.end(), sample_pts0.begin()+kr*dim,
			  sample_pts0.begin()+(kr+1)*dim);
      }

  static bool mod_inner = true;
  if (mod_inner)
    {
  // Let the modified surface approximate the sample points. Parameter 
  // iteration is performed
  ApproxSurf approx(surf, sample_pts, sample_pars, surf->dimension(), tol);

  // Keep the boundary
  approx.setFixBoundary(true);

  // Do not refine surfaces
  approx.setDoRefine(false);

  double maxdist, avdist;
  int nmb_out;
  shared_ptr<SplineSurface> mod_srf = approx.getApproxSurf(maxdist, avdist, nmb_out);
				
  if (mod_srf.get())
    surf->swap(*mod_srf);
    }
}


//===========================================================================
bool 
ModifySurf::enforceCoefCoLinearity(shared_ptr<SplineSurface> sf1, int bd1, 
				   shared_ptr<SplineSurface> sf2, int bd2, 
				   double tol, 
				   std::vector<std::vector<int> >& enumeration)
//===========================================================================
{
  // Collect surfaces
  int dim = sf1->dimension();
  vector<double>::iterator c1 = sf1->coefs_begin();
  vector<double>::iterator c2 = sf2->coefs_begin();
  vector<shared_ptr<SplineSurface> > sfs(2);
  sfs[0] = sf1;
  sfs[1] = sf2;

  // Define continuity constraints
  vector<sideConstraintSet> constraints(2*enumeration.size());
  for (size_t ki=0; ki<enumeration.size(); ++ki)
    {
      // Equality of coefficients at the boundary
      sideConstraintSet curr(dim);
      curr.factor_.push_back(make_pair(make_pair(0, enumeration[ki][1]), 1.0));
      curr.factor_.push_back(make_pair(make_pair(1, enumeration[ki][2]), -1.0));
      constraints[2*ki] = curr;

      // Co-linearity of coefficients at the boundary and the adjacent
      // rows
      // Compute distance between coeffients on either side of the
      // boundary
      double d1 = sqrt(Utils::distance_squared(&c1[enumeration[ki][0]*dim],
					&c1[(enumeration[ki][0]+1)*dim],
					&c1[enumeration[ki][1]*dim]));
      double d2 = sqrt(Utils::distance_squared(&c2[enumeration[ki][2]*dim],
					&c2[(enumeration[ki][2]+1)*dim],
					&c2[enumeration[ki][3]*dim]));
      sideConstraintSet curr2(dim);
      curr2.factor_.push_back(make_pair(make_pair(0, enumeration[ki][0]), -d2));
      curr2.factor_.push_back(make_pair(make_pair(0, enumeration[ki][1]), d2));
      curr2.factor_.push_back(make_pair(make_pair(1, enumeration[ki][2]), d1));
      curr2.factor_.push_back(make_pair(make_pair(1, enumeration[ki][3]), -d1));
      constraints[2*ki+1] = curr2;
    }

  // Define which coefficients that may change
  vector<vector<int> > coef_known(2);
  int kn1 = sf1->numCoefs_u();
  int kn2 = sf1->numCoefs_v();
  int kn3 = sf2->numCoefs_u();
  int kn4 = sf2->numCoefs_v();

  // 1. surface
  vector<int> ck1(kn1*kn2, 1); // Initially all coefficients are fixed

  // Release coefficients close to the boundary
  int min_nmb = 2;
  int max_nmb = (bd1 == 0 || bd1 == 1) ? kn1 - 2 : kn2 - 2;
  int nmb = 4;
  nmb = std::max(min_nmb, std::min(max_nmb, nmb));
  if (nmb > max_nmb)
    return false;  // Not enough degrees of freedom

  int ki, kj;
  int start_u = (bd1 == 1) ? kn1-nmb-1 : 0;
  int start_v = (bd1 == 3) ? kn2-nmb-1 : 0;
  int end_u = (bd1 == 0) ? nmb : kn1;
  int end_v = (bd1 == 2) ? nmb : kn2;
  if (bd1 == 0 || bd1 == 1)
    {
      start_v += 1;
      end_v -= 1;
    }
  else
    {
      start_u += 1;
      end_u -= 1;
    }
  for (kj=start_v; kj<end_v; ++kj)
    for (ki=start_u; ki<end_u; ++ki)
      ck1[kj*kn1+ki] = 0;

  coef_known[0] = ck1;

  // 2. surface
  vector<int> ck2(kn3*kn4, 1); // Initially all coefficients are fixed

  // Release coefficients close to the boundary
  max_nmb = (bd2 == 0 || bd2 == 1)  ? kn3 - 2 : kn4 - 2;
  nmb = 4;
  nmb = std::max(min_nmb, std::min(max_nmb, nmb));
  if (nmb > max_nmb)
    return false;  // Not enough degrees of freedom

  start_u = (bd2 == 1) ? kn3-nmb-1 : 0;
  start_v = (bd2 == 3) ? kn4-nmb-1 : 0;
  end_u = (bd2 == 0) ? nmb : kn3;
  end_v = (bd2 == 2) ? nmb : kn4;
  if (bd2 == 0 || bd2 == 1)
    {
      start_v += 1;
      end_v -= 1;
    }
  else
    {
      start_u += 1;
      end_u -= 1;
    }
  for (kj=start_v; kj<end_v; ++kj)
    for (ki=start_u; ki<end_u; ++ki)
      ck2[kj*kn3+ki] = 0;

  coef_known[1] = ck2;

  // Perform smoothing across the boundary
  // First set weights
  double w1 = 0.0;
  double w2 = 0.0025;
  double w3 = 0.0025;
  double w_orig = 0.05;

  SmoothSurfSet smooth(true);
  smooth.attach(sfs, coef_known, (int)constraints.size());

  smooth.setOptimize(w1, w2, w3);

  smooth.approxOrig(w_orig);

  smooth.setSideConstraints(constraints);

  vector<shared_ptr<SplineSurface> > sfs2;
  int stat = smooth.equationSolve(sfs2);
  if (stat < 0)
    return false;

  // Modify input surfaces
  sf1->swap(*sfs2[0]);
  sf2->swap(*sfs2[1]);
  return true;
}


//===========================================================================
bool 
ModifySurf::enforceVxCoefCoLinearity(vector<shared_ptr<SplineSurface> >& sfs, 
				     vector<int>& vx_enum, 
				     vector<pair<vector<int>, pair<int,int> > >& coef_cond,
				     double tol)
//===========================================================================
{
  if (sfs.size() == 0 || sfs.size() != vx_enum.size())
    return false;  // Inconsistent input

  int dim = sfs[0]->dimension();

  // Define continuity constraints
  // It is likely (depending on the input) that some constrains are defined
  // where there are no free coefficients involved. Those will be removed
  // during smoothing
  vector<sideConstraintSet> constraints(2*coef_cond.size());
  for (size_t ki=0; ki<coef_cond.size(); ++ki)
    {
      // Fetch surface coefficients
      int sf_idx1 = coef_cond[ki].second.first;
      int sf_idx2 = coef_cond[ki].second.second;
      vector<double>::iterator c1 = sfs[sf_idx1]->coefs_begin();
      vector<double>::iterator c2 = sfs[sf_idx2]->coefs_begin();
      // Equality of coefficients at the boundary
      sideConstraintSet curr(dim);
      curr.factor_.push_back(make_pair(make_pair(sf_idx1, coef_cond[ki].first[1]), 1.0));
      curr.factor_.push_back(make_pair(make_pair(sf_idx2, coef_cond[ki].first[2]), -1.0));
      constraints[2*ki] = curr;

      // Co-linearity of coefficients at the boundary and the adjacent
      // rows
      // Compute distance between coeffients on either side of the
      // boundary
      double d1 = sqrt(Utils::distance_squared(&c1[coef_cond[ki].first[0]*dim],
					&c1[(coef_cond[ki].first[0]+1)*dim],
					&c1[coef_cond[ki].first[1]*dim]));
      double d2 = sqrt(Utils::distance_squared(&c2[coef_cond[ki].first[2]*dim],
					&c2[(coef_cond[ki].first[2]+1)*dim],
					&c2[coef_cond[ki].first[3]*dim]));
      sideConstraintSet curr2(dim);
      curr2.factor_.push_back(make_pair(make_pair(sf_idx1, coef_cond[ki].first[0]), -d2));
      curr2.factor_.push_back(make_pair(make_pair(sf_idx1, coef_cond[ki].first[1]), d2));
      curr2.factor_.push_back(make_pair(make_pair(sf_idx2, coef_cond[ki].first[2]), d1));
      curr2.factor_.push_back(make_pair(make_pair(sf_idx2, coef_cond[ki].first[3]), -d1));
      constraints[2*ki+1] = curr2;
    }

  // For all surfaces define which coefficients that may change
  vector<vector<int> > coef_known(sfs.size());

  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      int kn1 = sfs[ki]->numCoefs_u();
      int kn2 = sfs[ki]->numCoefs_v();

      vector<int> ck1(kn1*kn2, 1); // Initially all coefficients are fixed

      // Release coefficients close to the vertex
      int min_nmb = 2;
      int nmb1 = 4, nmb2 = 4;
      nmb1 = std::min(kn1-2, nmb1);
      nmb2 = std::min(kn2-2, nmb2);
      if (nmb1 < min_nmb || nmb2 < min_nmb)
	return false;  // Not enough degrees of freedom

      int coefnmb = vx_enum[ki];
      int kr, kj;
      int start_u, start_v, end_u, end_v;
      if (coefnmb == 0 || coefnmb == (kn2-1)*kn1)
	{
	  start_u = 0;
	  end_u = nmb1;
	}
      else
	{
	  start_u = kn1 - nmb1;
	  end_u = kn1;
	}
      if (coefnmb == 0 || coefnmb == kn1-1)
	{
	  start_v = 0;
	  end_v = nmb2;
	}
      else
	{
	  start_v = kn2 - nmb2;
	  end_v = kn2;
	}

      for (kj=start_v; kj<end_v; ++kj)
	for (kr=start_u; kr<end_u; ++kr)
	  ck1[kj*kn1+kr] = 0;

      coef_known[ki] = ck1;
    }

  // Perform smoothing around the vertex
  // First set weights
  double w1 = 0.0;
  double w2 = 0.0025;
  double w3 = 0.0025;
  double w_orig = 0.05;

  SmoothSurfSet smooth(true);
  smooth.attach(sfs, coef_known, (int)constraints.size());

  smooth.setOptimize(w1, w2, w3);

  smooth.approxOrig(w_orig);

  smooth.setSideConstraints(constraints);

  vector<shared_ptr<SplineSurface> > sfs2;
  int stat = smooth.equationSolve(sfs2);
  if (stat < 0)
    return false;

  // Modify input surfaces
  for (size_t ki=0; ki<sfs.size(); ++ki)
    sfs[ki]->swap(*sfs2[ki]);

  return true;
}


