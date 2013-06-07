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

#include "GoTools/trivariate/SmoothVolume.h"
#include "GoTools/creators/Integrate.h"
#include "GoTools/creators/SolveCG.h"

#include <algorithm>


using std::vector;
using std::max;
using std::min;

namespace Go
{


  //===========================================================================
  SmoothVolume::SmoothVolume() :
    copy_coefs_(true)
  //===========================================================================
  {
  }


  //===========================================================================
  SmoothVolume::SmoothVolume(bool copy_coefs) :
    copy_coefs_(copy_coefs)
  //===========================================================================
  {
  }


  //===========================================================================
  SmoothVolume::~SmoothVolume()
  //===========================================================================
  {
  }


  //===========================================================================
  void SmoothVolume::attach(shared_ptr<SplineVolume>& in_vol,
			    vector<CoefStatus> status)
  //===========================================================================
  {
    // Currently reset everything. Later, refactor code to allow reuse of
    // previous result (some precalculated arrays, like bsplineintegral_u_ etc,
    // can then be reused)

    input_volume_ = in_vol;
    int ncoefs = in_vol->numCoefs(0) * in_vol->numCoefs(1) * in_vol->numCoefs(2);

    if (in_vol->rational())
      {
	if (copy_coefs_)
	  {
	    coef_array_.resize(ncoefs * (in_vol->dimension()+1));
	    it_coefs_ = coef_array_.begin();
	    copy(in_vol->rcoefs_begin(), in_vol->rcoefs_end(), it_coefs_);
	  }
	else
	  it_coefs_ = in_vol->rcoefs_begin();
      }
    else
      {
	if (copy_coefs_)
	  {
	    coef_array_.resize(ncoefs * in_vol->dimension());
	    it_coefs_ = coef_array_.begin();
	    copy(in_vol->coefs_begin(), in_vol->coefs_end(), it_coefs_);
	  }
	else
	  it_coefs_ = in_vol->coefs_begin();
      }

    coef_status_.resize(ncoefs);
    coef_other_.resize(ncoefs);
    for (int i = 0; i < ncoefs; ++i)
      coef_status_[i] = (status[i] == CoefOther) ? CoefFree : status[i];

    bsplineintegral_u_.resize(0);
    bsplineintegral_v_.resize(0);
    bsplineintegral_w_.resize(0);
    bsplineintegral_skew_u_.resize(0);
    bsplineintegral_skew_v_.resize(0);
    bsplineintegral_skew_w_.resize(0);

    reset();
  }


  //===========================================================================
  void SmoothVolume::reset()
  //===========================================================================
  {
    weight_der_1_ = 0.0;
    weight_der_2_ = 0.0;
    weight_der_3_ = 0.0;

    weight_least_sq_ = 0.0;
    least_sq_pts_.resize(0);
    least_sq_params_.resize(0);
    least_sq_wgt_.resize(0);

    for (int i = 0; i < 3; ++i)
      seem_cont_[i] = -1;
    for (int i = 0; i < 6; ++i)
      seem_weight_[i] = 0.0;
  }


  //===========================================================================
  void SmoothVolume::setOptimize(const double weight1,
				 const double weight2,
				 const double weight3)
  //===========================================================================
  {
    weight_der_1_ = weight1;
    weight_der_2_ = weight2;
    weight_der_3_ = weight3;
  }


  //===========================================================================
  void SmoothVolume::setLeastSquares(std::vector<double>& pnts,
				     std::vector<double>& param_pnts,
				     std::vector<double>& pnt_weights,
				     const double weight)
  //===========================================================================
  {
      int nmb_points = (int)pnts.size()/geoDim();   // Number of data points. 

    least_sq_pts_.resize(nmb_points * geoDim());
    copy (pnts.begin(), pnts.end(), least_sq_pts_.begin());

    least_sq_params_.resize(nmb_points * 3);
    copy (param_pnts.begin(), param_pnts.end(), least_sq_params_.begin());

    least_sq_wgt_.resize(nmb_points);
    copy (pnt_weights.begin(), pnt_weights.end(), least_sq_wgt_.begin());

    weight_least_sq_ = weight;
  }


  //===========================================================================
  void SmoothVolume::setPeriodicity(int pardir,
				    int cont,
				    double weight1,
				    double weight2)
  //===========================================================================
  {
    if (pardir < 0 || pardir > 2)
      return;

    seem_cont_[pardir] = cont;
    seem_weight_[pardir*2] = weight1;
    seem_weight_[pardir*2 + 1] = weight2;
  }


  //===========================================================================
  int SmoothVolume::equationSolve(shared_ptr<SplineVolume>& vol)
  //===========================================================================
  {
    // Prepare equation system
    setPeriodicityConstraints();
    resetPivotAndMatrices();
    if (!rational())
      buildIntegrals();
    else
      buildBsplineVolume();

    // Build equation system
    if (weight_least_sq_ != 0.0)
      addLeastSquares();
    if (rational())
      addOptimizeRational();
    else
      addOptimizeNonrational();
    for (int i = 0; i < 3; ++i)
      if (seem_cont_[i] > 0)
	{
	  if (rational())
	    addRationalContinuityAtSeem(i);
	  else
	    addNonrationalContinuityAtSeem(i);
	}

    // Solve equation system
    int kstat = solve();
    if (kstat < 0)
      return kstat;
    if (kstat == 1)
      THROW("Failed solving system (within tolerance)!");

    // Create volume
    std::vector<double>::iterator it = it_coefs_;
    int n_u = numCoefs(0), n_v = numCoefs(1), n_w = numCoefs(2);
    int g_dim = geoDim();
    int h_dim = homogDim();
    int coef_pos = 0;

    for (int k = 0; k < n_w; ++k)
      for (int j = 0; j < n_v; ++j)
	for (int i = 0; i < n_u; ++i, ++coef_pos)
	  {
	    if (coef_status_[coef_pos] == CoefKnown || coef_status_[coef_pos] == CoefAvoid)
	      {
		it += h_dim;
		continue;
	      }

	    int pos = pivot_[coef_pos];
	    for (int d = 0; d < g_dim; ++d, ++it, pos += nmb_free_)
	      (*it) = gright_[pos];
	    if (rational())
	      ++it;
	  }

    vol = shared_ptr<SplineVolume>(new SplineVolume(n_u, n_v, n_w,
						    order(0), order(1), order(2),
						    basis(0).begin(), basis(1).begin(), basis(2).begin(),
						    it_coefs_, g_dim, rational()));

    return 0;

  }




  //===========================================================================
  void SmoothVolume::setPeriodicityConstraints()
  //===========================================================================
  {
    for (int i = 0; i < (int)coef_status_.size(); ++i)
      if (coef_status_[i] == CoefOther)
	coef_status_[i] = CoefFree;

    int num_u = numCoefs(0);
    int num_v = numCoefs(1);
    int num_w = numCoefs(2);

    bool per_u = seem_cont_[0] >= 0;
    bool per_v = seem_cont_[1] >= 0;
    bool per_w = seem_cont_[2] >= 0;

    vector<int> corners(8);
    corners[0] = 0;
    corners[1] = num_u-1;
    corners[2] = num_u*(num_v-1);
    corners[3] = corners[1] + corners[2];
    corners[4] = num_u*num_v*(num_w-1);
    corners[5] = corners[1] + corners[4];
    corners[6] = corners[2] + corners[4];
    corners[7] = corners[3] + corners[4];

    // Tripple closed volume - could this ever happen? (In R^3, volume is then self intersecting)
    if (per_u && per_v && per_w)
      setPeriodicityLocally(corners, 1, 0);

    // Double closed volume in u and v par.dir.
    if (per_u && per_v)
      setEdgePeriodicity(corners[0], corners[1], corners[2], corners[3],
			 num_w, num_u*num_v, per_w);

    // Double closed volume in u and w par.dir.
    if (per_u && per_w)
      setEdgePeriodicity(corners[0], corners[1], corners[4], corners[5],
			 num_v, num_u, per_v);

    // Double closed volume in v and w par.dir.
    if (per_v && per_w)
      setEdgePeriodicity(corners[0], corners[2], corners[4], corners[6],
			 num_u, 1, per_u);

    // Closed volume in u par.dir.
    if (per_u)
      setFacePeriodicity(corners[0], corners[1], num_v, num_w,
			 num_u, num_u*num_v, per_v, per_w);

    // Closed volume in v par.dir.
    if (per_v)
      setFacePeriodicity(corners[0], corners[2], num_u, num_w,
			 1, num_u*num_v, per_u, per_w);

    // Closed volume in w par.dir.
    if (per_w)
      setFacePeriodicity(corners[0], corners[4], num_u, num_v,
			 1, num_u, per_u, per_v);
  }


  //===========================================================================
  void SmoothVolume::setEdgePeriodicity(int corner0, int corner1, int corner2, int corner3,
					int nmb, int step, bool per_end)
  //===========================================================================
  {
    vector<int> corners(4);
    corners[0] = corner0;
    corners[1] = corner1;
    corners[2] = corner2;
    corners[3] = corner3;

    int nmb_it = nmb;
    if (per_end)
      {
	nmb_it -= 2;
	for (int i = 0; i < 4; ++i)
	  corners[i] += step;
      }

    setPeriodicityLocally(corners, nmb_it, step);
  }


  //===========================================================================
  void SmoothVolume::setFacePeriodicity(int bottom, int top, int nmb_0, int nmb_1,
					int step_0, int step_1, bool per_0, bool per_1)
  //===========================================================================
  {
    vector<int> pos(2);
    pos[0] = bottom;
    pos[1] = top;

    int nmb_it_0 = nmb_0;
    if (per_0)
      {
	nmb_it_0 -= 2;
	pos[0] += step_0;
	pos[1] += step_0;
      }

    int nmb_it_1 = nmb_1;
    if (per_1)
      {
	nmb_it_1 -= 2;
	pos[0] += step_1;
	pos[1] += step_1;
      }

    for (int i = 0; i < nmb_it_1; ++i)
      {
	setPeriodicityLocally(pos, nmb_it_0, step_0);
	pos[0] += step_1;
	pos[1] += step_1;
      }
  }


  //===========================================================================
  void SmoothVolume::setPeriodicityLocally(const vector<int>& coefs_in, int nmb, int step)
  //===========================================================================
  {
      int nmb_coefs = (int)coefs_in.size();
    int g_dim = geoDim();
    int h_dim = homogDim();
    vector<int> coefs_copy(nmb_coefs);
    copy(coefs_in.begin(), coefs_in.end(), coefs_copy.begin());

    vector<double> pnt(g_dim, 0.0);

    for (int i = 0; i < nmb; ++i)   // For each set of control points that should be equal
      {

	// First count the fixed points, and find their centre
	int fixed = 0;
	for (int j = 0; j < nmb_coefs; ++j)
	  if (coef_status_[coefs_copy[j]] == CoefKnown)
	    {
	      if (fixed==0)
		for (int k = 0; k < g_dim; ++k)
		  pnt[k] = it_coefs_[coefs_copy[j]*h_dim + k];
	      else
		for (int k = 0; k < g_dim; ++k)
		  pnt[k] += it_coefs_[coefs_copy[j]*h_dim + k];
	      ++fixed;
	    }

	if (fixed == 0)   // No fixed point. Set up constraint to let all of them be the same as the first
	  for (int j = 1; j < nmb_coefs; ++j)
	    {
	      coef_status_[coefs_copy[j]] = CoefOther;
	      coef_other_[coefs_copy[j]] = coefs_copy[0];
	    }

	else   // Some fixed point. Let the other become the centre of the fixed points, og mark them as fixed as well
	  {
	    for (int k = 0; k < g_dim; ++k)
	      pnt[k] /= (double)fixed;
	    for (int j = 0; j < nmb_coefs; ++j)
	      if (coef_status_[coefs_copy[j]] == CoefFree)
		{
		  coef_status_[coefs_copy[j]] = CoefKnown;
		  for (int k = 0; k < g_dim; ++k)
		    it_coefs_[coefs_copy[j]*h_dim + k] = pnt[k];
		}
	  }

	// Prepare for next set of points
	for (int j = 0; j < nmb_coefs; ++j)
	  coefs_copy[j] += step;
      }
  }


  //===========================================================================
  void SmoothVolume::resetPivotAndMatrices()
  //===========================================================================
  {
    int n_coefs = numCoefs();

    // Build pivot table
    nmb_free_ = 0;
    pivot_.resize(n_coefs, 0);
    for (int i = 0; i < n_coefs; ++i)
      if (coef_status_[i] == CoefFree)
	{
	  pivot_[i] = nmb_free_;
	  ++nmb_free_;
	}

    // Build pivot table for coefficients to coincide with another coefficient
    for (int i = 0; i < n_coefs; ++i)
      if (coef_status_[i] == CoefOther)
	pivot_[i] = pivot_[coef_other_[i]];

    // Resize equation system matrices
    gmat_.resize(nmb_free_ * nmb_free_, 0.0);   // Matrix at left side of equation system.
    gright_.resize(geoDim() * nmb_free_, 0.0);            // Matrix at right side of equation system.
  }


  //===========================================================================
  void SmoothVolume::buildIntegrals()
  //===========================================================================
  {
    int d_u, d_v, d_w;   // Derivation depth of itegrals to calculate
    d_u = d_v = d_w = -1;

    if (seem_cont_[0] > 0)
      d_v = d_w = 0;
    if (seem_cont_[1] > 0)
      d_u = d_w = 0;
    if (seem_cont_[2] > 0)
      d_u = d_v = 0;

    int optim_depth = -1;
    if (weight_der_1_ != 0.0)
      optim_depth = 1;
    if (weight_der_2_ != 0.0)
      optim_depth = 2;
    if (weight_der_3_ != 0.0)
      optim_depth = 3;


	d_u = max(d_u, optim_depth);
    d_v = max(d_v, optim_depth);
    d_w = max(d_w, optim_depth);

    if (d_u >= 0)
      {
	int nmb_pairs = numCoefs(0) * (2*order(0) - 1);
	int old_der_depth = (int)bsplineintegral_u_.size() / nmb_pairs;
	if (d_u+1 > old_der_depth)
	  {
	    // We need to extend bsplineintegral_u_. Notice that resize() does not change the
	    // former contents in bsplineintegral_u_, so we only need to tell GaussQuadInnerFlat()
	    // to calculate the new values
	    bsplineintegral_u_.resize((d_u+1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(0);
	    GaussQuadInnerFlat(bas, d_u, old_der_depth, 0, bas.startparam(), bas.endparam(), bsplineintegral_u_);
	  }
	int old_skew_depth = (int)bsplineintegral_skew_u_.size() / nmb_pairs;
	if (optim_depth-1 > old_skew_depth)
	  {
	    bsplineintegral_skew_u_.resize((optim_depth-1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(0);
	    GaussQuadInnerFlat(bas, optim_depth-2, old_skew_depth, 2, bas.startparam(), bas.endparam(), bsplineintegral_skew_u_);
	  }
      }

    if (d_v >= 0)
      {
	int nmb_pairs = numCoefs(1) * (2*order(1) - 1);
	int old_der_depth = (int)bsplineintegral_v_.size() / nmb_pairs;
	if (d_v+1 > old_der_depth)
	  {
	    // We need to extend bsplineintegral_v_. Notice that resize() does not change the
	    // former contents in bsplineintegral_v_, so we only need to tell GaussQuadInnerFlat()
	    // to calculate the new values
	    bsplineintegral_v_.resize((d_v+1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(1);
	    GaussQuadInnerFlat(bas, d_v, old_der_depth, 0, bas.startparam(), bas.endparam(), bsplineintegral_v_);
	  }
	int old_skew_depth = (int)bsplineintegral_skew_v_.size() / nmb_pairs;
	if (optim_depth-1 > old_skew_depth)
	  {
	    bsplineintegral_skew_v_.resize((optim_depth-1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(1);
	    GaussQuadInnerFlat(bas, optim_depth-2, old_skew_depth, 2, bas.startparam(), bas.endparam(), bsplineintegral_skew_v_);
	  }
      }

    if (d_w >= 0)
      {
	int nmb_pairs = numCoefs(2) * (2*order(2) - 1);
	int old_der_depth = (int)bsplineintegral_w_.size() / nmb_pairs;
	if (d_w+1 > old_der_depth)
	  {
	    // We need to extend bsplineintegral_w_. Notice that resize() does not change the
	    // former contents in bsplineintegral_w_, so we only need to tell GaussQuadInnerFlat()
	    // to calculate the new values
	    bsplineintegral_w_.resize((d_w+1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(2);
	    GaussQuadInnerFlat(bas, d_w, old_der_depth, 0, bas.startparam(), bas.endparam(), bsplineintegral_w_);
	  }
	int old_skew_depth = (int)bsplineintegral_skew_w_.size() / nmb_pairs;
	if (optim_depth-1 > old_skew_depth)
	  {
	    bsplineintegral_skew_w_.resize((optim_depth-1) * nmb_pairs, 0.0);
	    BsplineBasis bas = basis(2);
	    GaussQuadInnerFlat(bas, optim_depth-2, old_skew_depth, 2, bas.startparam(), bas.endparam(), bsplineintegral_skew_w_);
	  }
      }
  }


  //===========================================================================
  void SmoothVolume::buildBsplineVolume()
  //===========================================================================
  {
    int n_coefs = numCoefs();
    int h_dim = homogDim();
    vector<double> bspl_coefs(n_coefs * 2, 0.0);

    for (int pos_in = geoDim(), pos_out = 1;
	 pos_in < h_dim * n_coefs;
	 pos_in += h_dim, pos_out +=2)
      bspl_coefs[pos_out] = it_coefs_[pos_in];

    bspline_volume_ = shared_ptr<SplineVolume>
      (new SplineVolume(numCoefs(0), numCoefs(1), numCoefs(2),
			order(0), order(1), order(2),
			basis(0).begin(), basis(1).begin(), basis(2).begin(),
			bspl_coefs.begin(), 1, true));
  }


  //===========================================================================
  void SmoothVolume::addLeastSquares()
  //===========================================================================
  {
    int nmb_pts = (int)least_sq_wgt_.size();
    int ncoefs0 = numCoefs(0);
    int ncoefs1 = numCoefs(1);
    int order0 = order(0);
    int order1 = order(1);
    int order2 = order(2);
    int g_dim = geoDim();
    int h_dim = homogDim();
    vector<double> bas0(order0), bas1(order1), bas2(order2);
    vector<double> tp_basis(order0 * order1 * order2);  // Tensor product of basis functions

    vector<double>::const_iterator pnt_it = least_sq_pts_.begin();
    vector<double>::const_iterator param_it = least_sq_params_.begin();
    for (int pt_cnt = 0; pt_cnt < nmb_pts; ++pt_cnt, param_it += 3, pnt_it += g_dim)  // For every point to be approximated
      {

	int left0, left1, left2;

	// Fetch B-spline tensor products different from zero
	if (rational())
	  {
	    Point p(1);
	    vector<double>::iterator bspl_it = bspline_volume_->rcoefs_begin();
	    left0 = basis(0).knotInterval(param_it[0]);
	    left1 = basis(1).knotInterval(param_it[1]);
	    left2 = basis(2).knotInterval(param_it[2]);

	    bspl_it += 2 * (left0 - order0 + 1
			    + ncoefs0 * (left1 - order1 + 1
					 + ncoefs1 * (left2 - order2 + 1)));
	    for (int k = 0; k < order2; ++k)
	      for (int j = 0; j < order1; ++j)
		for (int i = 0; i < order0; ++i)
		  {
		    int pos = 2*(i + ncoefs0*(j + ncoefs1*k));
		    bspl_it[pos] = 1.0;
		    bspline_volume_->point(p, param_it[0], param_it[1], param_it[2]);
		    tp_basis[i + order0*(j + order1*k)] = p[0];
		    bspl_it[pos] = 0.0;
		  }
	  }
	else
	  {
	    basis(0).computeBasisValues(param_it[0], &bas0[0], 0);
	    basis(1).computeBasisValues(param_it[1], &bas1[0], 0);
	    basis(2).computeBasisValues(param_it[2], &bas2[0], 0);

	    left0 = basis(0).lastKnotInterval();
	    left1 = basis(1).lastKnotInterval();
	    left2 = basis(2).lastKnotInterval();

	    // Compute the tensor product of basis functions.
	    int pos = 0;
	    for (int k = 0; k < order2; ++k)
	      for (int j = 0; j < order1; ++j)
		for (int i = 0; i < order0; ++i, ++pos)
		  tp_basis[pos] = bas0[i] * bas1[j] * bas2[k];
	  }

	// Run through all pairs of coefficients where the B-spline
	// tensor product has support in the point
	for (int r = left2 - order2 + 1, b_pos_pqr = 0; r <= left2; ++r)    // For every w-dir B-spline, first coeff
	  for (int q = left1 - order1 + 1; q <= left1; ++q)    // For every v-dir B-spline, first coeff
	    for (int p = left0 - order0 + 1; p <= left0; ++p, ++b_pos_pqr)    // For every u-dir B-spline, first coeff
	      {
		int pos_pqr = p + ncoefs0 * (q + ncoefs1 * r);
		if (coef_status_[pos_pqr] == CoefKnown || coef_status_[pos_pqr] == CoefAvoid)
		  continue;

		int piv0 = pivot_[pos_pqr];
		double term_pqr = weight_least_sq_ * least_sq_wgt_[pt_cnt] * tp_basis[b_pos_pqr];

		// Add contribution to right hand side
		for (int d = 0; d < g_dim; ++d)
		  gright_[d*nmb_free_ + piv0] += term_pqr * pnt_it[d];

		for (int k = left2 - order2 + 1, b_pos_ijk = 0; k <= left2; ++k)    // For every w-dir B-spline, second coeff
		  for (int j = left1 - order1 + 1; j <= left1; ++j)    // For every v-dir B-spline, second coeff
		    for (int i = left0 - order0 + 1; i <= left0; ++i, ++b_pos_ijk)    // For every u-dir B-spline, second coeff
		      {
			int pos_ijk = i + ncoefs0 * (j + ncoefs1 * k);
			if (coef_status_[pos_ijk] == CoefAvoid)
			  continue;

			double term = term_pqr * tp_basis[b_pos_ijk];

			if (coef_status_[pos_ijk] == CoefKnown)
			  {
			    // Add contribution to right hand side
			    vector<double>::const_iterator coef_it = it_coefs_ + h_dim * pos_ijk;
			    for (int d = 0; d < g_dim; ++d, ++coef_it)
			      gright_[d*nmb_free_ + piv0] += term * (*coef_it);
			  }
			else
			  {
			    // Add contribution to left hand side
			    int piv1 = pivot_[pos_ijk];
			    if (piv1>piv0)
			      continue;
			    gmat_[piv0 * nmb_free_ + piv1] += term;
			    if (piv1<piv0)
			      gmat_[piv1 * nmb_free_ + piv0] += term;
			  }

		      }   // End -- For every B-spline tensor product, second coeff
	      }  // End -- For every B-spline tensor product, first coeff

      }     // End -- For every point to be approximated
  }


  //===========================================================================
  void SmoothVolume::addOptimizeNonrational()
  //===========================================================================
  {
    int derivs = 3;
    if (weight_der_3_ == 0.0)
      {
	derivs = 2;
	if (weight_der_2_ == 0.0)
	  {
	    derivs = 1;
	    if (weight_der_1_ == 0.0)
	      return;
	  }
      }

    int ncoefs0 = numCoefs(0);
    int ncoefs1 = numCoefs(1);
    int ncoefs2 = numCoefs(2);
    int order0 = order(0);
    int order1 = order(1);
    int order2 = order(2);
    int g_dim = geoDim();
    int h_dim = homogDim();

    int int_step_u = ncoefs0 * (2 * order0 - 1);
    int int_step_v = ncoefs1 * (2 * order1 - 1);
    int int_step_w = ncoefs2 * (2 * order2 - 1);

    double factor_1 = weight_der_1_;
    double factor_2 = weight_der_2_ / 5.0;
    double factor_3 = weight_der_3_ / 35.0;

    // Travers all B-splines and set up matrices of equation system.
    for (int k1 = 0, pos_1 = 0; k1 < ncoefs2; ++k1)  // For each B-spline in third direction, first B-spline tripple
      for (int j1 = 0; j1 < ncoefs1; ++j1)  // For each B-spline in second direction, first B-spline tripple
	for (int i1 = 0; i1 < ncoefs0; ++i1, ++pos_1)  // For each B-spline in first direction, first B-spline tripple
	  {
	    if (coef_status_[pos_1] == CoefKnown || coef_status_[pos_1] == CoefAvoid)
	      continue;

	    int piv_1 = pivot_[pos_1];

	    for (int k2 = 0, pos_2 = 0; k2 < ncoefs2; ++k2)  // For each B-spline in third direction, second B-spline tripple
	      {
		if (k2 <= k1-order2 || k2 >= k1+order2)
		  {
		    pos_2 += ncoefs0*ncoefs1;
		    continue;     // B-splines do not overlap
		  }

		vector<double>::const_iterator int_it_w = bsplineintegral_w_.begin() + (k2-k1+order2-1 + k1*(2*order2-1));
		vector<double>::const_iterator int_sk_it_w = bsplineintegral_skew_w_.begin() + (k2-k1+order2-1 + k1*(2*order2-1));
		vector<double>::const_iterator int_sk_op_it_w = bsplineintegral_skew_w_.begin() + (k1-k2+order2-1 + k2*(2*order2-1));

		for (int j2 = 0; j2 < ncoefs1; ++j2)  // For each B-spline in second direction, second B-spline tripple
		  {
		    if (j2 <= j1-order1 || j2 >= j1+order1)
		      {
			pos_2 += ncoefs0;
			continue;     // B-splines do not overlap
		      }

		    vector<double>::const_iterator int_it_v = bsplineintegral_v_.begin() + (j2-j1+order1-1 + j1*(2*order1-1));
		    vector<double>::const_iterator int_sk_it_v = bsplineintegral_skew_v_.begin() + (j2-j1+order1-1 + j1*(2*order1-1));
		    vector<double>::const_iterator int_sk_op_it_v = bsplineintegral_skew_v_.begin() + (j1-j2+order1-1 + j2*(2*order1-1));

		    for (int i2 = 0; i2 < ncoefs0; ++i2, ++pos_2)  // For each B-spline in first direction, second B-spline tripple
		      {
			if (i2 <= i1-order0 || i2 >= i1+order0)
			  continue;     // B-splines do not overlap

			if (coef_status_[pos_2] == CoefAvoid)
			  continue;
			int piv_2 = pivot_[pos_2];
			if (piv_2 < piv_1 && (coef_status_[pos_2] == CoefFree || coef_status_[pos_2] == CoefOther))
			  continue;

			vector<double>::const_iterator int_it_u = bsplineintegral_u_.begin() + (i2-i1+order0-1 + i1*(2*order0-1));
			vector<double>::const_iterator int_sk_it_u = bsplineintegral_skew_u_.begin() + (i2-i1+order0-1 + i1*(2*order0-1));
			vector<double>::const_iterator int_sk_op_it_u = bsplineintegral_skew_u_.begin() + (i1-i2+order0-1 + i2*(2*order0-1));

			// We have a contribution, compute it
			double eval_1 = 0.0, eval_2 = 0.0, eval_3 = 0.0;

			// First order derivatives
			if (derivs > 0)
			  eval_1 =
			    // d_u^2, d_v^2 and d_w^2
			    int_it_u[int_step_u] * int_it_v[0] * int_it_w[0] +
			    int_it_v[int_step_v] * int_it_u[0] * int_it_w[0] +
			    int_it_w[int_step_w] * int_it_u[0] * int_it_v[0];

			// Second order derivatives
			if (derivs > 1)
			  eval_2 =
			    // d_uu^2 + d_vv^2 + d_ww^2
			    3.0 * (int_it_u[2*int_step_u] * int_it_v[0] * int_it_w[0] +
				   int_it_v[2*int_step_v] * int_it_u[0] * int_it_w[0] +
				   int_it_w[2*int_step_w] * int_it_u[0] * int_it_v[0]) +
			    // d_uv^2 + d_uw^2 + d_vw^2
			    4.0 * (int_it_u[int_step_u] * int_it_v[int_step_v] * int_it_w[0] +
				   int_it_u[int_step_u] * int_it_w[int_step_w] * int_it_v[0] +
				   int_it_v[int_step_v] * int_it_w[int_step_w] * int_it_u[0]) +
			    // d_uu*d_vv + d_uu*d_ww + d_vv*d_ww
			    1.0 * (int_sk_it_u[0] * int_sk_op_it_v[0] * int_it_w[0] +
				   int_sk_it_u[0] * int_sk_op_it_w[0] * int_it_v[0] +
				   int_sk_it_v[0] * int_sk_op_it_u[0] * int_it_w[0] +
				   int_sk_it_v[0] * int_sk_op_it_w[0] * int_it_u[0] +
				   int_sk_it_w[0] * int_sk_op_it_u[0] * int_it_v[0] +
				   int_sk_it_w[0] * int_sk_op_it_v[0] * int_it_u[0]);

			// Third order derivatives
			if (derivs > 2)
			  eval_3 =
			    // d_uuu^2 + d_vvv^2 + d_www^2
			    15.0 * (int_it_u[3*int_step_u] * int_it_v[0] * int_it_w[0] +
				    int_it_v[3*int_step_v] * int_it_u[0] * int_it_w[0] +
				    int_it_w[3*int_step_w] * int_it_u[0] * int_it_v[0]) +
			    // d_uuv^2 + d_uvv^2 + d_uuw^2 + d_uww^2 + d_vvw^2 + d_vww^2
			    27.0 * (int_it_u[2*int_step_u] * int_it_v[int_step_v] * int_it_w[0] +
				    int_it_v[2*int_step_v] * int_it_u[int_step_u] * int_it_w[0] +
				    int_it_u[2*int_step_u] * int_it_w[int_step_w] * int_it_v[0] +
				    int_it_w[2*int_step_w] * int_it_u[int_step_u] * int_it_v[0] +
				    int_it_v[2*int_step_v] * int_it_w[int_step_w] * int_it_u[0] +
				    int_it_w[2*int_step_w] * int_it_v[int_step_v] * int_it_u[0]) +
			    // d_uvw^2
			    36.0 * int_it_u[int_step_u] * int_it_v[int_step_v] * int_it_w[int_step_w] +
			    // d_uuu*duvv + d_uuu*duww + d_vvv*duuv + d_vvv*d_vww + d_www*d_uuw + d_www*dvvw
			    9.0 * (int_sk_it_u[int_step_u] * int_sk_op_it_v[0] * int_it_w[0] +
				   int_sk_op_it_u[int_step_u] * int_sk_it_v[0] * int_it_w[0] +
				   int_sk_it_u[int_step_u] * int_sk_op_it_w[0] * int_it_v[0] +
				   int_sk_op_it_u[int_step_u] * int_sk_it_w[0] * int_it_v[0] +
				   int_sk_it_v[int_step_v] * int_sk_op_it_u[0] * int_it_w[0] +
				   int_sk_op_it_v[int_step_v] * int_sk_it_u[0] * int_it_w[0] +
				   int_sk_it_v[int_step_v] * int_sk_op_it_w[0] * int_it_u[0] +
				   int_sk_op_it_v[int_step_v] * int_sk_it_w[0] * int_it_u[0] +
				   int_sk_it_w[int_step_w] * int_sk_op_it_u[0] * int_it_v[0] +
				   int_sk_op_it_w[int_step_w] * int_sk_it_u[0] * int_it_v[0] +
				   int_sk_it_w[int_step_w] * int_sk_op_it_v[0] * int_it_u[0] +
				   int_sk_op_it_w[int_step_w] * int_sk_it_v[0] * int_it_u[0]) +
			    // d_uuw*d_vvw + d_uuv*d_vww + d_uww*d_vww
			    9.0 * (int_sk_it_u[0] * int_sk_op_it_v[0] * int_it_w[int_step_w] +
				   int_sk_it_v[0] * int_sk_op_it_u[0] * int_it_w[int_step_w] +
				   int_sk_it_u[0] * int_sk_op_it_w[0] * int_it_v[int_step_v] +
				   int_sk_it_w[0] * int_sk_op_it_u[0] * int_it_v[int_step_v] +
				   int_sk_it_v[0] * int_sk_op_it_w[0] * int_it_u[int_step_u] +
				   int_sk_it_w[0] * int_sk_op_it_v[0] * int_it_u[int_step_u]);

			// Add contribution to equation system
			double term =
			  (factor_1*eval_1 + factor_2*eval_2 + factor_3*eval_3);

			if (coef_status_[pos_2] == CoefKnown)
			  {
			    // The contribution of this term is added to the right
			    // side of the equation system.

			    vector<double>::const_iterator coef_it = it_coefs_ + pos_2 * h_dim;
			    for (int d = 0; d < g_dim; ++d, ++coef_it)
			      gright_[d*nmb_free_ + piv_1] -= term * (*coef_it);
			  }
			else
			  {
			    // The contribution of this term is added to the left
			    //  side of the equation system.

			    gmat_[piv_1 * nmb_free_ + piv_2] += term;
			    if (piv_1 != piv_2)
			      gmat_[piv_2 * nmb_free_ + piv_1] += term;
			  }
		      }  // End -- For each B-spline in first direction, second B-spline tripple
		  }  // End -- For each B-spline in second direction, second B-spline tripple
	      }  // End -- For each B-spline in third direction, second B-spline tripple
	  }  // End -- For each B-spline in 1., 2. and 3. direction, first B-spline tripple
  }


  //===========================================================================
  void SmoothVolume::addOptimizeRational()
  //===========================================================================
  {
    // For the rational case, we can not split the double integral into a product
    // of two separate integrals for each parameter value (like the non-rational case),
    // because of the denominator function depending on both variables.
    // Instead we use the grid evaluator to evaluate each spline product in
    // points for numerical integration by Gaussian quadrature.

    int derivs = 3;
    if (weight_der_3_ == 0.0)
      {
	derivs = 2;
	if (weight_der_2_ == 0.0)
	  {
	    derivs = 1;
	    if (weight_der_1_ == 0.0)
	      return;
	  }
      }

    int ncoefs0 = numCoefs(0);
    int ncoefs1 = numCoefs(1);
    int ncoefs2 = numCoefs(2);
    int order0 = order(0);
    int order1 = order(1);
    int order2 = order(2);
    int g_dim = geoDim();
    int h_dim = homogDim();

    BsplineBasis bas0 = basis(0);
    BsplineBasis bas1 = basis(1);
    BsplineBasis bas2 = basis(2);
    vector<double>::const_iterator knots0 = bas0.begin();
    vector<double>::const_iterator knots1 = bas1.begin();
    vector<double>::const_iterator knots2 = bas2.begin();

    // Get points and weights for integration
    vector<double> gq_points_u, gq_points_v, gq_points_w;
    vector<double> gq_weights_u, gq_weights_v, gq_weights_w;
    GaussQuadValues(bas0, gq_points_u, gq_weights_u);
    GaussQuadValues(bas1, gq_points_v, gq_weights_v);
    GaussQuadValues(bas2, gq_points_w, gq_weights_w);
    int seg_samples_u = (int)gq_weights_u.size();   // Number of samples inside each Bezier segment in first direction
    int seg_samples_v = (int)gq_weights_v.size();   // Number of samples inside each Bezier segment in second direction
    int seg_samples_w = (int)gq_weights_w.size();   // Number of samples inside each Bezier segment in third direction
    int bez_segs_u = ncoefs0 - order0 + 1;   // Number of Bezier segments in first direction
    int bez_segs_v = ncoefs1 - order1 + 1;   // Number of Bezier segments in second direction
    int bez_segs_w = ncoefs2 - order2 + 1;   // Number of Bezier segments in third direction
    int samples_u = bez_segs_u * seg_samples_u;  // Number of samples in total ( == gq_points_u.size() ) in first direction
    int samples_v = bez_segs_v * seg_samples_v;  // Number of samples in total ( == gq_points_v.size() ) in second direction
    int samples_w = bez_segs_w * seg_samples_w;  // Number of samples in total ( == gq_points_w.size() ) in third direction

    // Get B-spline evaluations in each sample point
    vector<double> basisvals_u(samples_u * order0 * (derivs + 1));   // The B-spline values and derivatives for integration points in first direction
    vector<int> left_u(samples_u);
    bas0.computeBasisValues(&gq_points_u[0], &gq_points_u[samples_u], &basisvals_u[0], &left_u[0], derivs);
    vector<double> basisvals_v(samples_v * order1 * (derivs + 1));   // The B-spline values and derivatives for integration points in second direction
    vector<int> left_v(samples_v);
    bas1.computeBasisValues(&gq_points_v[0], &gq_points_v[samples_v], &basisvals_v[0], &left_v[0], derivs);
    vector<double> basisvals_w(samples_w * order2 * (derivs + 1));   // The B-spline values and derivatives for integration points in third direction
    vector<int> left_w(samples_w);
    bas2.computeBasisValues(&gq_points_w[0], &gq_points_w[samples_w], &basisvals_w[0], &left_w[0], derivs);

    // For each triple (Bi,Cj,Dk) of B-splines in 1., 2. and 3. direction, get all needed combinations of
    // partial derivatives of the function B_i*C_j*d_k/G in all integration sample points where B_i, C_j and D_k
    // have support, where G is the denominator function of the volume

    // The results are organized in the vector sample_evaluations. It has one element for each tripple (Bi,Cj,Dk)
    // starting with (i,j,k)=(0,0,0), then (1,0,0), (2,0,0) etc, then (0,1,0), (1,1,0) etc, then (0,0,1), (1,0,1), etc
    // For each B-spline tripple, the entry is a vector of evaluations, organized in four levels. The first
    // (outermost) level is for each sample value in third paramter directions. The second level is for each
    // sample value in second direction. The third level is for each sample value in first direction.
    // The fourth (innermost) level is for all derivatives in the given sample value tripple, starting with the
    // function value, then d_u, d_v, d_w, d_uu, d_uv, d_uw, d_vv, d_vw, d_ww, d_uuu etc
    vector<vector<double> > sample_evaluations;
    sample_evaluations.resize(ncoefs0*ncoefs1*ncoefs2);

    // Some variables to describe the Bezier segments where each B-spline has support
    vector<int> left_bas_u(ncoefs0);  // The first Bezier segment where each B-spline has support
    vector<int> left_bas_v(ncoefs1);
    vector<int> left_bas_w(ncoefs2);
    vector<int> segs_bas_u(ncoefs0);  // The number of Bezier segments where each B-spline has support
    vector<int> segs_bas_v(ncoefs1);
    vector<int> segs_bas_w(ncoefs2);
    for (int i = 0; i < ncoefs0; ++i)
      {
	left_bas_u[i] = std::max(0, i-order0+1);
	int right_bas = std::min(i, ncoefs0-order0);
	segs_bas_u[i] = right_bas - left_bas_u[i] + 1;
      }
    for (int i = 0; i < ncoefs1; ++i)
      {
	left_bas_v[i] = std::max(0, i-order1+1);
	int right_bas = std::min(i, ncoefs1-order1);
	segs_bas_v[i] = right_bas - left_bas_v[i] + 1;
      }
    for (int i = 0; i < ncoefs2; ++i)
      {
	left_bas_w[i] = std::max(0, i-order2+1);
	int right_bas = std::min(i, ncoefs2-order2);
	segs_bas_w[i] = right_bas - left_bas_w[i] + 1;
      }


    // Fill sample_evaluations
    int part_derivs = ((derivs+1)*(derivs+2)*(derivs+3))/6;
    vector<double>::iterator bspl_it = bspline_volume_->rcoefs_begin();
    for (int k = 0, samp_ev_pos = 0; k < ncoefs2; ++k)
      {
	int curr_samp_w = seg_samples_w*segs_bas_w[k];  // Number of sample points where the spline has support
	vector<double>::iterator start_basisvals_w = basisvals_w.begin() + left_bas_w[k]*order2*(derivs+1)*seg_samples_w;
	vector<int>::iterator start_left_w = left_w.begin() + left_bas_w[k]*seg_samples_w;
	for (int j = 0; j < ncoefs1; ++j)
	  {
	    int curr_samp_v = seg_samples_v*segs_bas_v[j];  // Number of sample points where the spline has support
	    vector<double>::iterator start_basisvals_v = basisvals_v.begin() + left_bas_v[j]*order1*(derivs+1)*seg_samples_v;
	    vector<int>::iterator start_left_v = left_v.begin() + left_bas_v[j]*seg_samples_v;
	    for (int i = 0; i < ncoefs0; ++i, ++samp_ev_pos, bspl_it+=2)
	      {
		int curr_samp_u = seg_samples_u*segs_bas_u[i];  // Number of sample points where the spline has support
		vector<double>::iterator start_basisvals_u = basisvals_u.begin() + left_bas_u[i]*order0*(derivs+1)*seg_samples_u;
		vector<int>::iterator start_left_u = left_u.begin() + left_bas_u[i]*seg_samples_u;
		sample_evaluations[samp_ev_pos].resize(curr_samp_u*curr_samp_v*curr_samp_w*part_derivs);
		*bspl_it = 1.0;
		bspline_volume_->pointsGrid(curr_samp_u, curr_samp_v, curr_samp_w,
					    start_basisvals_u, start_basisvals_v, start_basisvals_w,
					    start_left_u, start_left_v, start_left_w,
					    derivs, sample_evaluations[samp_ev_pos]);
		*bspl_it = 0.0;
	      }
	  }
      }

    double factor_1 = weight_der_1_;
    double factor_2 = weight_der_2_ / 5.0;
    double factor_3 = weight_der_3_ / 35.0;

    // Travers all B-splines and set up matrices of equation system.
    for (int k1 = 0, pos_1 = 0; k1 < ncoefs2; ++k1)  // For each B-spline in third direction, first B-spline tripple
      for (int j1 = 0; j1 < ncoefs1; ++j1)  // For each B-spline in second direction, first B-spline tripple
	for (int i1 = 0; i1 < ncoefs0; ++i1, ++pos_1)  // For each B-spline in first direction, first B-spline tripple
	  {
	    if (coef_status_[pos_1] == CoefKnown || coef_status_[pos_1] == CoefAvoid)
	      continue;

	    int piv_1 = pivot_[pos_1];

	    for (int k2 = 0, pos_2 = 0; k2 < ncoefs2; ++k2)  // For each B-spline in third direction, second B-spline tripple
	      {
		if (left_bas_w[k1] + segs_bas_w[k1] <= left_bas_w[k2] ||
		    left_bas_w[k2] + segs_bas_w[k2] <= left_bas_w[k1])
		  {
		    pos_2 += ncoefs0*ncoefs1;
		    continue;     // B-splines do not overlap
		  }

		int first_seg_w = std::max(left_bas_w[k1], left_bas_w[k2]);
		int end_seg_w = std::min(left_bas_w[k1]+segs_bas_w[k1], left_bas_w[k2]+segs_bas_w[k2]);
		for (int j2 = 0; j2 < ncoefs1; ++j2)  // For each B-spline in second direction, second B-spline tripple
		  {
		    if (left_bas_v[j1] + segs_bas_v[j1] <= left_bas_v[j2] ||
			left_bas_v[j2] + segs_bas_v[j2] <= left_bas_v[j1])
		      {
			pos_2 += ncoefs0;
			continue;     // B-splines do not overlap
		      }

		    int first_seg_v = std::max(left_bas_v[j1], left_bas_v[j2]);
		    int end_seg_v = std::min(left_bas_v[j1]+segs_bas_v[j1], left_bas_v[j2]+segs_bas_v[j2]);
		    for (int i2 = 0; i2 < ncoefs0; ++i2, ++pos_2)  // For each B-spline in first direction, second B-spline tripple
		      {
			if (left_bas_u[i1] + segs_bas_u[i1] <= left_bas_u[i2] ||
			    left_bas_u[i2] + segs_bas_u[i2] <= left_bas_u[i1])
			  continue;     // B-splines do not overlap

			if (coef_status_[pos_2] == CoefAvoid)
			  continue;
			int piv_2 = pivot_[pos_2];
			if (piv_2 < piv_1 && (coef_status_[pos_2] == CoefFree || coef_status_[pos_2] == CoefOther))
			  continue;

			// Now we know we have a contribution, and need to evaluate all integral terms

			int first_seg_u = std::max(left_bas_u[i1], left_bas_u[i2]);
			int end_seg_u = std::min(left_bas_u[i1]+segs_bas_u[i1], left_bas_u[i2]+segs_bas_u[i2]);

			double term = 0.0;    // The final contribution for these two B-spline tripples to the equation system

			for (int seg_k = first_seg_w; seg_k < end_seg_w; ++seg_k)   // For each common Bezier segment in third direction
			  for (int seg_j = first_seg_v; seg_j < end_seg_v; ++seg_j)   // For each common Bezier segment in second direction
			    for (int seg_i = first_seg_u; seg_i < end_seg_u; ++seg_i)   // For each common Bezier segment in first direction
			      {
				double vol_size = (knots0[seg_i + order0] - knots0[seg_i + order0 - 1])
				  * (knots1[seg_j + order1] - knots1[seg_j + order1 - 1])
				  * (knots2[seg_k + order2] - knots2[seg_k + order2 - 1]);
				vector<double>::const_iterator samp_eval_1
				  = sample_evaluations[pos_1].begin()
				  + part_derivs * seg_samples_u
				  * (seg_i-left_bas_u[i1]
				     + seg_samples_v * segs_bas_u[i1]
				     * (seg_j-left_bas_v[j1]
					+ seg_samples_w * segs_bas_v[j1] * (seg_k-left_bas_w[k1])));
				vector<double>::const_iterator samp_eval_2
				  = sample_evaluations[pos_2].begin()
				  + part_derivs * seg_samples_u
				  * (seg_i-left_bas_u[i2]
				     + seg_samples_v * segs_bas_u[i2]
				     * (seg_j-left_bas_v[j2]
					+ seg_samples_w * segs_bas_v[j2] * (seg_k-left_bas_w[k2])));

				for (int samp_k = 0; samp_k < seg_samples_w; ++samp_k)   // For each sample in third direction
				  {
				    for (int samp_j = 0; samp_j < seg_samples_v; ++samp_j)   // For each sample in second direction
				      {
					for (int samp_i = 0;          // For each sample in first direction
					     samp_i < seg_samples_u;
					     ++samp_i, samp_eval_1 += part_derivs, samp_eval_2 += part_derivs)
					  {
					    // Compute contribution in given sample tripple
					    double eval_1 = 0.0, eval_2 = 0.0, eval_3 = 0.0;

					    // First order derivatives
					    if (derivs > 0)
					      eval_1 =
						// d_u^2, d_v^2 and d_w^2
						samp_eval_1[1] * samp_eval_2[1] +
						samp_eval_1[2] * samp_eval_2[2] +
						samp_eval_1[3] * samp_eval_2[3];

					    // Second order derivatives
					    if (derivs > 1)
					      eval_2 =
						// d_uu^2 + d_vv^2 + d_ww^2
						3.0 * (samp_eval_1[4] * samp_eval_2[4] +
						       samp_eval_1[7] * samp_eval_2[7] +
						       samp_eval_1[9] * samp_eval_2[9]) +
						// d_uv^2 + d_uw^2 + d_vw^2
						4.0 * (samp_eval_1[5] * samp_eval_2[5] +
						       samp_eval_1[6] * samp_eval_2[6] +
						       samp_eval_1[8] * samp_eval_2[8]) +
						// d_uu*d_vv + d_uu*d_ww + d_vv*d_ww
						1.0 * (samp_eval_1[4] * samp_eval_2[7] +
						       samp_eval_1[4] * samp_eval_2[9] +
						       samp_eval_1[7] * samp_eval_2[4] +
						       samp_eval_1[7] * samp_eval_2[9] +
						       samp_eval_1[9] * samp_eval_2[4] +
						       samp_eval_1[9] * samp_eval_2[7]);
					    // Third order derivatives
					    if (derivs > 2)
					      eval_3 =
						// d_uuu^2 + d_vvv^2 + d_www^2
						15.0 * (samp_eval_1[10] * samp_eval_2[10] +
							samp_eval_1[16] * samp_eval_2[16] +
							samp_eval_1[19] * samp_eval_2[19]) +
						// d_uuv^2 + d_uvv^2 + d_uuw^2 + d_uww^2 + d_vvw^2 + d_vww^2
						27.0 * (samp_eval_1[11] * samp_eval_2[11] +
							samp_eval_1[13] * samp_eval_2[13] +
							samp_eval_1[12] * samp_eval_2[12] +
							samp_eval_1[15] * samp_eval_2[15] +
							samp_eval_1[17] * samp_eval_2[17] +
							samp_eval_1[18] * samp_eval_2[18]) +
						// d_uvw^2
						36.0 * samp_eval_1[14] * samp_eval_2[14] +
						// d_uuu*duvv + d_uuu*duww + d_vvv*duuv + d_vvv*d_vww + d_www*d_uuw + d_www*dvvw
						9.0 * (samp_eval_1[10] * samp_eval_2[13] +
						       samp_eval_1[13] * samp_eval_2[10] +
						       samp_eval_1[10] * samp_eval_2[15] +
						       samp_eval_1[15] * samp_eval_2[10] +
						       samp_eval_1[16] * samp_eval_2[11] +
						       samp_eval_1[11] * samp_eval_2[16] +
						       samp_eval_1[16] * samp_eval_2[18] +
						       samp_eval_1[18] * samp_eval_2[16] +
						       samp_eval_1[19] * samp_eval_2[12] +
						       samp_eval_1[12] * samp_eval_2[19] +
						       samp_eval_1[19] * samp_eval_2[17] +
						       samp_eval_1[17] * samp_eval_2[19]) +
						// d_uuw*d_vvw + d_uuv*d_vww + d_uvv*d_uww
						9.0 * (samp_eval_1[12] * samp_eval_2[17] +
						       samp_eval_1[17] * samp_eval_2[12] +
						       samp_eval_1[11] * samp_eval_2[18] +
						       samp_eval_1[18] * samp_eval_2[11] +
						       samp_eval_1[13] * samp_eval_2[15] +
						       samp_eval_1[15] * samp_eval_2[13]);

					    // Sum together, and add to term
					    term +=
					      (factor_1*eval_1 + factor_2*eval_2 + factor_3*eval_3)
					      * gq_weights_u[samp_i]* gq_weights_v[samp_j] * gq_weights_w[samp_k]
					      * vol_size;
					  }  // End -- For each sample in first direction

					if (samp_j < seg_samples_v-1 || samp_k < seg_samples_w-1)
					  {
					    samp_eval_1 += part_derivs * seg_samples_u * (segs_bas_u[i1]-1);
					    samp_eval_2 += part_derivs * seg_samples_u * (segs_bas_u[i2]-1);
					  }
				      }  // End -- For each sample in second direction

				    if (samp_k < seg_samples_w-1)
				      {
					samp_eval_1 +=
					  part_derivs * seg_samples_u * seg_samples_v
					  * segs_bas_u[i1] * (segs_bas_v[j1]-1);
					samp_eval_2 +=
					  part_derivs * seg_samples_u * seg_samples_v
					  * segs_bas_u[i2] * (segs_bas_v[j2]-1);
				      }
				  }  // End -- For each sample in third direction
			      } // End -- For each common bezier segment in 1. 2. and 3. direction

			// Add contribution to equation system
			if (coef_status_[pos_2] == CoefKnown)
			  {
			    // The contribution of this term is added to the right
			    // side of the equation system.

			    vector<double>::const_iterator coef_it = it_coefs_ + pos_2 * h_dim;
			    for (int d = 0; d < g_dim; ++d, ++coef_it)
			      gright_[d*nmb_free_ + piv_1] -= term * (*coef_it);
			  }
			else
			  {
			    // The contribution of this term is added to the left
			    //  side of the equation system.

			    gmat_[piv_1 * nmb_free_ + piv_2] += term;
			    if (piv_1 != piv_2)
			      gmat_[piv_2 * nmb_free_ + piv_1] += term;
			  }
		      }  // End -- For each B-spline in first direction, second B-spline tripple
		  }  // End -- For each B-spline in second direction, second B-spline tripple
	      }  // End -- For each B-spline in third direction, second B-spline tripple
	  }  // End -- For each B-spline in 1., 2. and 3. direction, first B-spline tripple

  }


  //===========================================================================
  void SmoothVolume::addNonrationalContinuityAtSeem(int pardir)
  //===========================================================================
  {
    // For pardir=0, we look at continuity at the seem in u-direction.
    // The constraint is given by minimizing the double integral over all v,w
    // of (V_u(u_min,v,w)-V_u(u_max,v,w))^2, and (only in case cn=2) the integral
    // of (V_uu(u_min,v,w)-V_uu(u_max,v,w))^2 where V_u and V_uu are
    // first and second order differentiation of the volume w.r.t. u.
    // For pardir=1 and 2, similar formulas are used for continuity in v- and w-direction
    //
    // Variables ending with _cont are for continuity parameter direction,
    // _int1 and _int2 are for first and second integral parameter direction respectively
    int int_dir1 = (pardir == 0) ? 1 : 0;
    int int_dir2 = (pardir == 2) ? 1 : 2;

    int order_cont = order(pardir);
    int order_int1 = order(int_dir1);
    int order_int2 = order(int_dir2);

    double weight1 = seem_weight_[2*pardir];
    double weight2 = seem_weight_[2*pardir+1];

    // cn==1 if only C1-continuity should be calculated. cn==2 includes C2-continuity
    int cn = min(order_cont-1, nmb_free_);
    cn = min(cn, seem_cont_[pardir]);
    if (cn > 2)
      cn = 2;
    if (cn == 2 && weight2 == 0.0)
      cn = 1;
    if (cn < 1 || (cn == 1 && weight1 == 0.0))
      return;

    int ncoefs_cont = numCoefs(pardir);
    int ncoefs_int1 = numCoefs(int_dir1);
    int ncoefs_int2 = numCoefs(int_dir2);

    vector<double> integral_int1 = (pardir == 0) ? bsplineintegral_v_ : bsplineintegral_u_;
    vector<double> integral_int2 = (pardir == 2) ? bsplineintegral_v_ : bsplineintegral_w_;

    int g_dim = geoDim();
    int h_dim = homogDim();

    int step_cont, step_int1, step_int2;
    if (pardir == 0)
      {
	step_cont = 1;
	step_int1 = numCoefs(0);
	step_int2 = numCoefs(1)*step_int1;
      }
    else if (pardir == 1)
      {
	step_int1 = 1;
	step_cont = numCoefs(0);
	step_int2 = numCoefs(1)*step_cont;
      }
    else
      {
	step_int1 = 1;
	step_int2 = numCoefs(0);
	step_cont = numCoefs(1)*step_int2;
      }

    BsplineBasis bas_cont = basis(pardir);
    BsplineBasis bas_int1 = basis(int_dir1);
    BsplineBasis bas_int2 = basis(int_dir2);
    vector<double>::const_iterator knot_cont = bas_cont.begin();

    double d_left_0 = 1.0/(knot_cont[order_cont] - bas_cont.startparam());
    double d_right_0 = 1.0/(bas_cont.endparam() - knot_cont[ncoefs_cont-1]);
    double d_left_1 = 0.0, d_right_1 = 0.0;
    if (cn == 2)
      {
	d_left_1 = 1.0/(knot_cont[order_cont+1] - bas_cont.startparam());
	d_right_1 = 1.0/(bas_cont.endparam() - knot_cont[ncoefs_cont-2]);
      }

    // We run through the coefficients in continuity parameter direction in this order:
    // First the leftmost, then the rightmost, then next to leftmost and next to rightmost,
    // And for C2-continuity also next to next to leftmost and next to next to rightmost

    // prod_cont holds the B-spline product in continuity direction.
    double prod_cont[6][6];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
	prod_cont[i][j] = 0.0;

    // First C1-part
    vector<double> b_spl_cont(6);
    b_spl_cont[0] = -(order_cont-1) * d_left_0;
    b_spl_cont[1] = -(order_cont-1) * d_right_0;
    b_spl_cont[2] = -b_spl_cont[0];
    b_spl_cont[3] = -b_spl_cont[1];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
	prod_cont[i][j] = weight1 * b_spl_cont[i] * b_spl_cont[j];

    // Then C2-part
    if (cn == 2)
      {
	b_spl_cont[0] = (order_cont-1) * (order_cont-2) * d_left_0 * d_left_0;
	b_spl_cont[1] = -(order_cont-1) * (order_cont-2) * d_right_0 * d_right_0;
	b_spl_cont[4] = (order_cont-1) * (order_cont-2) * d_left_0 * d_left_1;
	b_spl_cont[5] = -(order_cont-1) * (order_cont-2) * d_right_0 * d_right_1;
	b_spl_cont[2] = -(b_spl_cont[0] + b_spl_cont[4]);
	b_spl_cont[3] = -(b_spl_cont[1] + b_spl_cont[5]);
	for (int i = 0; i < 6; ++i)
	  for (int j = 0; j < 6; ++j)
	    prod_cont[i][j] += weight2 * b_spl_cont[i] * b_spl_cont[j];
      }

    // Now run through all coefficient products to get contributions to equation system
    for (int r = 0; r < ncoefs_int2; ++r)   // For every position in second integration dir, first coefficient
      for (int q = 0; q < ncoefs_int1; ++q)   // For every pos in first integr dir, first coeff
	for (int p_idx = 0; p_idx < 2*cn+2; ++p_idx)   // For every pos in continuity dir, first coeff
	  {
	    int p = ((p_idx & 1) == 0) ? (p_idx>>1) : (ncoefs_cont - 1 - (p_idx>>1));
	    int pos_pqr = p * step_cont + q * step_int1 + r * step_int2;

	    if (coef_status_[pos_pqr] == CoefKnown || coef_status_[pos_pqr] == CoefAvoid)
	      continue;

	    int piv0 = pivot_[pos_pqr];

	    for (int k = 0; k < ncoefs_int2; ++k)   // For every pos in second integr dir, second coeff
	      if (k < r + order_int2 && r < k + order_int2)  // Check if B-splines have common support
		for (int j = 0; j < ncoefs_int1; ++j)   // For every pos in first integr dir, second coeff
		  if (j < q + order_int1 && q < j + order_int1)  // Check if B-splines have common support
		    {
		      double integral_prod =
			integral_int1[j-q+order_int1-1 + q*(2*order_int1-1)] *
			integral_int2[k-r+order_int2-1 + r*(2*order_int2-1)];

		      for (int i_idx = 0; i_idx < 2*cn+2; ++i_idx)   // For every pos in continuity dir, second coeff
			{
			  int i = ((i_idx & 1) == 0) ? (i_idx>>1) : (ncoefs_cont - 1 - (i_idx>>1));
			  int pos_ijk = i * step_cont + j * step_int1 + k * step_int2;
			  if (coef_status_[pos_ijk] == CoefAvoid)
			    continue;

			  double term = prod_cont[p_idx][i_idx] * integral_prod;
			  if (coef_status_[pos_ijk] == CoefKnown)
			    {
			      // Add contribution to right hand side
			      vector<double>::const_iterator coef_it = it_coefs_ + h_dim * pos_ijk;
			      for (int d = 0; d < g_dim; ++d, ++coef_it)
				gright_[d*nmb_free_ + piv0] += term * (*coef_it);
			    }
			  else
			    {
			      // Add contribution to left hand side
			      int piv1 = pivot_[pos_ijk];
			      if (piv1>piv0)
				continue;
			      gmat_[piv0 * nmb_free_ + piv1] += term;
			      if (piv1<piv0)
				gmat_[piv1 * nmb_free_ + piv0] += term;
			    }
			}   // End -- For every pos in continuity dir, second coeff
		    }  // End -- For every choice in integral directions, second coefficient
	  }  // End -- For every choice of first coefficient
  }


  //===========================================================================
  void SmoothVolume::addRationalContinuityAtSeem(int pardir)
  //===========================================================================
  {
    // To set up continuity for the rational case, the main idea is the same
    // as for non-rational, except that we can not use the integrals of
    // products of splines, as the denominator functions prevents us from
    // Splitting the tripple integral into products of univariable integrals.
    // Instead, the integration is estimated using Gaussian quadratures
    int int_dir1 = (pardir == 0) ? 1 : 0;
    int int_dir2 = (pardir == 2) ? 1 : 2;

    int order_cont = order(pardir);
    int order_int1 = order(int_dir1);
    int order_int2 = order(int_dir2);

    double weight1 = seem_weight_[2*pardir];
    double weight2 = seem_weight_[2*pardir+1];

    // cn==1 if only C1-continuity should be calculated. cn==2 includes C2-continuity
    int cn = min(order_cont-1, nmb_free_);
    cn = min(cn, seem_cont_[pardir]);
    if (cn > 2)
      cn = 2;
    if (cn == 2 && weight2 == 0.0)
      cn = 1;
    if (cn < 1 || (cn == 1 && weight1 == 0.0))
      return;

    int ncoefs_cont = numCoefs(pardir);
    int ncoefs_int1 = numCoefs(int_dir1);
    int ncoefs_int2 = numCoefs(int_dir2);

    int g_dim = geoDim();
    int h_dim = homogDim();

    int step_cont, step_int1, step_int2;
    if (pardir == 0)
      {
	step_cont = 1;
	step_int1 = numCoefs(0);
	step_int2 = numCoefs(1)*step_int1;
      }
    else if (pardir == 1)
      {
	step_int1 = 1;
	step_cont = numCoefs(0);
	step_int2 = numCoefs(1)*step_cont;
      }
    else
      {
	step_int1 = 1;
	step_int2 = numCoefs(0);
	step_cont = numCoefs(1)*step_int2;
      }

    BsplineBasis bas_cont = basis(pardir);
    BsplineBasis bas_int1 = basis(int_dir1);
    BsplineBasis bas_int2 = basis(int_dir2);
    vector<double>::const_iterator knot_cont = bas_cont.begin();
    vector<double>::const_iterator knot_int1 = bas_int1.begin();
    vector<double>::const_iterator knot_int2 = bas_int2.begin();

    double d_left_0 = 1.0/(knot_cont[order_cont] - bas_cont.startparam());
    double d_right_0 = 1.0/(bas_cont.endparam() - knot_cont[ncoefs_cont-1]);
    double d_left_1 = 0.0, d_right_1 = 0.0;

    // Define t_left_ij as D^j(B_i)(u_min) and t_right_ij as D^j(B_(ncoefs-i-1)(u_max)) for some i, j
    double t_left_01 = -(order_cont-1) * d_left_0;
    double t_left_11 = -t_left_01;
    double t_right_01 = (order_cont-1) * d_right_0;
    double t_right_11 = -t_right_01;
    double
      t_left_02 = 0.0, t_left_12 = 0.0, t_left_22 = 0.0,
      t_right_02 = 0.0, t_right_12 = 0.0, t_right_22 = 0.0;
    if (cn == 2)
      {
	d_left_1 = 1.0/(knot_cont[order_cont+1] - bas_cont.startparam());
	d_right_1 = 1.0/(bas_cont.endparam() - knot_cont[ncoefs_cont-2]);
	t_left_02 = (order_cont-1) * (order_cont-2) * d_left_0 * d_left_0;
	t_left_22 = (order_cont-1) * (order_cont-2) * d_left_0 * d_left_1;
	t_left_12 = -(t_left_02+t_left_22);
	t_right_02 = (order_cont-1) * (order_cont-2) * d_right_0 * d_right_0;
	t_right_22 = (order_cont-1) * (order_cont-2) * d_right_0 * d_right_1;
	t_right_12 = -(t_right_02+t_right_22);
      }

    // Denom_surf is a spline surface in 3 (cn==1) or 5 (cn==2) dimensions, and is used to
    // find the values of the denominator.
    // If u is the continuity direction, and D is the denominateor function,
    // the coordinates of denom_surf(v,w) are
    //
    // - D(u_min,v,w)   (the same as D(u_max,v,w))
    // - D_u(u_min,v,w)
    // - D_u(u_max,v,w)
    // - D_uu(u_min,v,w)   (only for cn=2)
    // - D_uu(u_max,v,w)   (only for cn=2)
    vector<double> denomcoefs(ncoefs_int1 * ncoefs_int2 * (1 + 2*cn));
    for (int j = 0, pos = 0; j < ncoefs_int2; ++j)
      for (int i = 0; i < ncoefs_int1; ++i)
	{
	  vector<double>::const_iterator it = it_coefs_ + (j * step_int2 + i * step_int1) * h_dim + g_dim;
	  double l_0 = it[0];
	  double l_1 = it[h_dim * step_cont];
	  double r_1 = it[h_dim * step_cont * (ncoefs_cont-2)];
	  denomcoefs[pos++] = l_0;
	  denomcoefs[pos++] = l_0*t_left_01 + l_1*t_left_11;
	  denomcoefs[pos++] = l_0*t_right_01 + r_1*t_right_11;
	  if (cn == 2)
	    {
	      double l_2 = it[2 * h_dim * step_cont];
	      double r_2 = it[h_dim * step_cont * (ncoefs_cont-3)];
	      denomcoefs[pos++] = l_0*t_left_02 + l_1*t_left_12 + l_2*t_left_22;
	      denomcoefs[pos++] = l_0*t_right_02 + r_1*t_right_12 + r_2*t_right_22;
	    }
	}
    shared_ptr<SplineSurface> denom_surf(new SplineSurface(ncoefs_int1, ncoefs_int2,
							   order_int1, order_int2,
							   bas_int1.begin(), bas_int2.begin(),
							   denomcoefs.begin(), 1 + 2*cn, false));

    // Get parameter sample values to evaluate at
    vector<double> gq_points_1, gq_points_2;   // All sample parameter values along the two integration paths
    vector<double> gq_weights_1, gq_weights_2;   // The weights used for numerical integration
    GaussQuadValues(bas_int1, gq_points_1, gq_weights_1);
    GaussQuadValues(bas_int2, gq_points_2, gq_weights_2);

    int bez_segs_1 = ncoefs_int1 - order_int1 + 1;   // Number of Bezier segments in first integration direction
    int bez_segs_2 = ncoefs_int2 - order_int2 + 1;   // Number of Bezier segments in second integration direction
    int seg_samples_1 = (int)gq_weights_1.size();   // Number of samples inside each Bezier segment in first integration direction
    int seg_samples_2 = (int)gq_weights_2.size();   // Number of samples inside each Bezier segment in second integration direction
    int samples_1 = bez_segs_1 * seg_samples_1;  // Number of samples in total ( == gq_points_1.size() ) in first int. dir.
    int samples_2 = bez_segs_2 * seg_samples_2;  // Number of samples in total ( == gq_points_2.size() ) in second int. dir.

    // Get B-spline evaluations in sample points
    vector<double> basisvals_1(samples_1 * order_int1);   // The B-spline values for sample points in first direction
    vector<int> left_1(samples_1);
    bas_int1.computeBasisValues(&gq_points_1[0], &gq_points_1[samples_1], &basisvals_1[0], &left_1[0]);
    vector<double> basisvals_2(samples_2 * order_int2);   // The B-spline values for sample points in second direction
    vector<int> left_2(samples_2);
    bas_int2.computeBasisValues(&gq_points_2[0], &gq_points_2[samples_2], &basisvals_2[0], &left_2[0]);

    // Get denominater and its derivatives in continuity parameter direction at the two isosurfaces to be glued together
    vector<double> denom_samples(samples_1 * samples_2 * (1 + 2*cn));
    denom_surf->pointsGrid(samples_1, samples_2, 0,
			   &basisvals_1[0], &basisvals_2[0],
			   &left_1[0], &left_2[0],
			   &denom_samples[0]);

    // Get som evaluations of derivatives of B_i/D in sample points at the edge surfaces to be glued together.
    // When pardir=0, bspl_denom_1[k*samples_1*4 + j*4 + i] will for the sample values v_j and w_k hold
    //
    // - B_0/D derivated w.r.t. u and evaluated in (u_min,v_j,w_k) when i = 0
    // - -B_(n-1)/D derivated w.r.t. u and evaluated in (u_max,v_j,w_k) when i = 1
    // - B_1/D derivated w.r.t. u and evaluated in (u_min,v_j,w_k) when i = 2
    // - -B_(n-2)/D derivated w.r.t. u and evaluated in (u_max,v_j,w_k) when i = 3
    // 
    // where n = number of coefficients in first parameter direction. When cn==2, then also
    // bspl_denom_2[k*samples_1*6 + j*6 + i] will hold B_0/D, -B_(n-1)/D, B_1/D, -B_(n-2)/D, B_2/D, -B_(n-3)/D
    // twice derivated, and evaluated in (u_min,v_j,w_k) when i is even, or (u_max,v_j,w_k) when i is odd
    // for i = 0,1,2,3,4,5.
    //
    // bspl_denom_1 and bspl_denom_2 are used in a similar way for pardir=1 or 2.
    vector<double> bspl_denom_1, bspl_denom_2;
    bspl_denom_1.resize(samples_1*samples_2*4);
    if (cn==2)
      bspl_denom_2.resize(samples_1*samples_2*6);
    for (int k = 0, pos1 = 0, pos2 = 0, pos_denom = 0;
	 k < samples_2; ++k)
      for (int j = 0; j < samples_1; ++j, pos_denom += 1 + 2*cn)
	{
	  double d1 = denom_samples[pos_denom];
	  double d2 = d1*d1;
	  bspl_denom_1[pos1++] = (t_left_01*d1 - denom_samples[pos_denom+1]) / d2;
	  bspl_denom_1[pos1++] = -(t_right_01*d1 - denom_samples[pos_denom+2]) / d2;
	  bspl_denom_1[pos1++] = t_left_11 / d1;
	  bspl_denom_1[pos1++] = -t_right_11 / d1;

	  if (cn==2)
	    {
	      double d3 = d2*d1;
	      bspl_denom_2[pos2++] = (t_left_02*d2 - 2*t_left_01*d1*denom_samples[pos_denom+1]
				      - d1*denom_samples[pos_denom+3] + 2*denom_samples[pos_denom+1]*denom_samples[pos_denom+1]) / d3;
	      bspl_denom_2[pos2++] = -(t_right_02*d2 - 2*t_right_01*d1*denom_samples[pos_denom+2]
				       - d1*denom_samples[pos_denom+4] + 2*denom_samples[pos_denom+2]*denom_samples[pos_denom+2]) / d3;
	      bspl_denom_2[pos2++] = (t_left_12*d1 - 2*t_left_11*denom_samples[pos_denom+1]) / d2;
	      bspl_denom_2[pos2++] = -(t_right_12*d1 - 2*t_right_11*denom_samples[pos_denom+2]) / d2;
	      bspl_denom_2[pos2++] = t_left_22/d1;
	      bspl_denom_2[pos2++] = -t_right_22/d1;
	    }
	}

    // Some variables to describe the Bezier segments where each B-spline has support
    vector<int> left_bas_1(ncoefs_int1);  // The first Bezier segment where each B-spline has support
    vector<int> left_bas_2(ncoefs_int2);
    vector<int> segs_bas_1(ncoefs_int1);  // The number of Bezier segments where each B-spline has support
    vector<int> segs_bas_2(ncoefs_int2);
    for (int i = 0; i < ncoefs_int1; ++i)
      {
	left_bas_1[i] = std::max(0, i-order_int1+1);
	int right_bas = std::min(i, ncoefs_int1-order_int1);
	segs_bas_1[i] = right_bas - left_bas_1[i] + 1;
      }
    for (int i = 0; i < ncoefs_int2; ++i)
      {
	left_bas_2[i] = std::max(0, i-order_int2+1);
	int right_bas = std::min(i, ncoefs_int2-order_int2);
	segs_bas_2[i] = right_bas - left_bas_2[i] + 1;
      }

    // Now run through all coefficient products to get contributions to equation system
    // For variable p and i, we run through the coefficients in continuity parameter direction
    // in this order:
    // First the leftmost, then the rightmost, then next to leftmost and next to rightmost,
    // and for C2-continuity also next to next to leftmost and next to next to rightmost

    // Now run through all coefficient products to get contributions to equation system
    for (int r = 0; r < ncoefs_int2; ++r)   // For every position in second integration dir, first coefficient
      for (int q = 0; q < ncoefs_int1; ++q)   // For every pos in first integr dir, first coeff
	for (int p_idx = 0; p_idx < 2*cn+2; ++p_idx)   // For every pos in continuity dir, first coeff
	  {
	    int p = ((p_idx & 1) == 0) ? (p_idx>>1) : (ncoefs_cont - 1 - (p_idx>>1));
	    int pos_pqr = p * step_cont + q * step_int1 + r * step_int2;

	    if (coef_status_[pos_pqr] == CoefKnown || coef_status_[pos_pqr] == CoefAvoid)
	      continue;

	    int piv0 = pivot_[pos_pqr];

	    for (int k = 0; k < ncoefs_int2; ++k)   // For every pos in second integr dir, second coeff
	      if (k < r + order_int2 && r < k + order_int2)  // Check if B-splines have common support
		{
		  int first_seg_2 = std::max(left_bas_2[k], left_bas_2[r]);
		  int end_seg_2 = std::min(left_bas_2[k]+segs_bas_2[k], left_bas_2[r]+segs_bas_2[r]);

		  for (int j = 0; j < ncoefs_int1; ++j)   // For every pos in first integr dir, second coeff
		    if (j < q + order_int1 && q < j + order_int1)  // Check if B-splines have common support
		      {
			int first_seg_1 = std::max(left_bas_1[j], left_bas_1[q]);
			int end_seg_1 = std::min(left_bas_1[j]+segs_bas_1[j], left_bas_1[q]+segs_bas_1[q]);

			for (int i_idx = 0; i_idx < 2*cn+2; ++i_idx)   // For every pos in continuity dir, second coeff
			  {

			    int i = ((i_idx & 1) == 0) ? (i_idx>>1) : (ncoefs_cont - 1 - (i_idx>>1));
			    int pos_ijk = i * step_cont + j * step_int1 + k * step_int2;
			    int piv1 = pivot_[pos_ijk];
			    if (piv1 < piv0 && (coef_status_[pos_ijk] == CoefFree || coef_status_[pos_ijk] == CoefOther))
			      continue;

			    // Now we know we have a contribution, and need to evaluate all integral terms

			    double term = 0.0;    // The final contribution for these two B-spline tripples to the equation system

			    for (int seg_2 = first_seg_2; seg_2 < end_seg_2; ++seg_2)   // For each common Bezier segment in second int. dir.
			      for (int seg_1 = first_seg_1; seg_1 < end_seg_1; ++seg_1)   // For each common Bez. seg. in first int. dir.
				{
				  double area
				    = (knot_int1[seg_1 + order_int1] - knot_int1[seg_1 + order_int1 - 1])
				    * (knot_int2[seg_2 + order_int2] - knot_int2[seg_2 + order_int2 - 1]);

				  // Positions in bspl_denom_1 and bspl_denom_2
				  int pbd1 = 4 * (seg_samples_1 * seg_1 + samples_1 * seg_samples_2 * seg_2);
				  int pbd2 = 6 * (seg_samples_1 * seg_1 + samples_1 * seg_samples_2 * seg_2);

				  for (int samp_2 = 0, samp_pos_2 = seg_2 * seg_samples_2;
				       samp_2 < seg_samples_2; ++samp_2, ++samp_pos_2)   // For each sample in second direction
				    {
				      int rel_pos_k = k + order_int2 - 1 -left_2[samp_pos_2];
				      int rel_pos_r = r + order_int2 - 1 -left_2[samp_pos_2];
				      if (rel_pos_k < order_int2 && rel_pos_k >= 0
					  && rel_pos_r < order_int2 && rel_pos_r >= 0)
					{
					  double term_2 = area
					    * basisvals_2[order_int2 * samp_pos_2 + rel_pos_k]
					    * basisvals_2[order_int2 * samp_pos_2 + rel_pos_r]
					    * gq_weights_2[samp_2];
					  for (int samp_1 = 0, samp_pos_1 = seg_1 * seg_samples_1;
					       samp_1 < seg_samples_1;
					       ++samp_1, pbd1 += 4, pbd2 += 6, ++samp_pos_1)  // For each sample in first direction
					    {
					      int rel_pos_j = j + order_int1 - 1 -left_1[samp_pos_1];
					      int rel_pos_q = q + order_int1 - 1 -left_1[samp_pos_1];
					      if (rel_pos_j < order_int1 && rel_pos_j >= 0
						  && rel_pos_q < order_int1 && rel_pos_q >= 0)
						{
						  double term_1 = 0.0;
						  if (i_idx < 4 && p_idx < 4)
						    term_1 += weight1 * bspl_denom_1[pbd1+i_idx] * bspl_denom_1[pbd1+p_idx];
						  if (cn == 2)
						    term_1 += weight2 * bspl_denom_2[pbd2+i_idx] * bspl_denom_2[pbd2+p_idx];
						  term += term_2 * term_1
						    * basisvals_1[order_int1 * samp_pos_1 + rel_pos_j]
						    * basisvals_1[order_int1 * samp_pos_1 + rel_pos_q]
						    * gq_weights_1[samp_1];
						}
					    }  // End -- For each sample in first direction
					}
				      if (samp_2 < seg_samples_2-1)
					{
					  pbd1 += 4 * (samples_1 - seg_samples_1);
					  pbd2 += 6 * (samples_1 - seg_samples_1);
					}
				    }   // End -- For each sample in second direction
				}   // End -- For each common Bezier segment in both integration directions

			    // Add contribution to equation system
			    if (coef_status_[pos_ijk] == CoefKnown)
			      {
				// The contribution of this term is added to the right
				// side of the equation system.

				vector<double>::const_iterator coef_it = it_coefs_ + pos_ijk * h_dim;
				for (int d = 0; d < g_dim; ++d, ++coef_it)
				  gright_[d*nmb_free_ + piv0] -= term * (*coef_it);
			      }
			    else
			      {
				// The contribution of this term is added to the left
				//  side of the equation system.

				gmat_[piv0 * nmb_free_ + piv1] += term;
				if (piv0 != piv1)
				  gmat_[piv1 * nmb_free_ + piv0] += term;
			      }

			  }   // End -- For every pos in continuity dir, second coeff
		      }  // End -- For every pos in first integr dir, second coeff
		}  // End -- For every pos in second integr dir, second coeff
	  }  // End -- For every choice of first coefficients (values of p, q, r)
  }


  //===========================================================================
  int SmoothVolume::solve()
  //===========================================================================
  {
    int g_dim = geoDim();
    int h_dim = homogDim();
    int n_coefs = numCoefs();

    // Solve the equation system by Conjugate Gradient Method.
    vector<double> eb(g_dim * nmb_free_);
    copy(gright_.begin(), gright_.end(), eb.begin());

    // Copy coefficients to array of unknowns
    for (int i = 0; i < n_coefs; ++i)
      if (coef_status_[i] != CoefKnown && coef_status_[i] != CoefAvoid)
	for (int j = 0, pos = pivot_[i]; j < g_dim; ++j, pos += nmb_free_)
	  gright_[pos] = it_coefs_[i*h_dim + j];


    // Set up CG-object
    SolveCG solveCg;

    // Create sparse matrix.
    ASSERT(gmat_.size() > 0);
    solveCg.attachMatrix(&gmat_[0], nmb_free_);

    // Attach parameters.
    solveCg.setTolerance(0.00000001);
    //solveCg.setMaxIterations(std::min(2*nmb_free_, 1000));
    solveCg.setMaxIterations(std::min(20*nmb_free_, 100000));

    // Precondition
    int precond =  1;
    if (precond) {
      double omega = 0.1;
      solveCg.precondRILU(omega);
    }

    // Solve equation systems.
    for (int i = 0; i < g_dim; ++i)
      {
	int kstat = solveCg.solve(&gright_[i*nmb_free_], &eb[i*nmb_free_], nmb_free_);
	if (kstat < 0 || kstat == 1)
	  return kstat;
      }

    // Copy result to output array. 
    for (int i = 0; i < n_coefs; ++i)
      if (coef_status_[i] != CoefKnown && coef_status_[i] != CoefAvoid)
	for (int j = 0, pos = pivot_[i]; j < g_dim; ++j, pos += nmb_free_)
	  it_coefs_[i*h_dim + j] = gright_[pos];

    return 0;
  }


  //===========================================================================
  int SmoothVolume::geoDim() const
  //===========================================================================
  {
    return input_volume_->dimension();
  }


  //===========================================================================
  int SmoothVolume::homogDim() const
  //===========================================================================
  {
    return input_volume_->dimension()
      + ((input_volume_->rational()) ? 1 : 0);
  }


  //===========================================================================
  bool SmoothVolume::rational() const
  //===========================================================================
  {
    return input_volume_->rational();
  }


  //===========================================================================
  int SmoothVolume::numCoefs() const
  //===========================================================================
  {
    return
      input_volume_->numCoefs(0)
      * input_volume_->numCoefs(1)
      * input_volume_->numCoefs(2);
  }


  //===========================================================================
  int SmoothVolume::numCoefs(int pardir) const
  //===========================================================================
  {
    return input_volume_->numCoefs(pardir);
  }


  //===========================================================================
  int SmoothVolume::order(int pardir) const
  //===========================================================================
  {
    return input_volume_->order(pardir);
  }


  //===========================================================================
  BsplineBasis SmoothVolume::basis(int pardir) const
  //===========================================================================
  {
    return input_volume_->basis(pardir);
  }



} // namespace Go
