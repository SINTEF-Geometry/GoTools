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

#include <string>
#include <algorithm>
#include "newmat.h"
#include "GoTools/compositemodel/ftSmoothSurf.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/SplineSurface.h"
#include <math.h> // Really needed?
#include <cstdio> // for debugging
#include <fstream>

//#define DEBUG

using std::vector;
using std::min;
using std::max;
using namespace NEWMAT;

namespace Go
{

//===========================================================================
ftSmoothSurf::ftSmoothSurf(shared_ptr<SplineSurface> surf, double approxtol,
			   double approx_orig_tol,
			   vector<int> ccw_edge_derivs, int maxiter,
			   bool lock_corner_points)
  : surf_(surf), approxtol_(approxtol), approx_orig_tol_(approx_orig_tol),
    init_approx_weight_(0.8), ccw_edge_derivs_(ccw_edge_derivs), max_error_(-1.0),
    mean_error_(-1.0), maxiter_(maxiter), 
    lock_corner_points_(lock_corner_points)
  //===========================================================================
{
  orig_surf_ = shared_ptr<SplineSurface>(new SplineSurface(*surf_));
  seem_[0]=0;
  seem_[1]=0;
}

//===========================================================================
ftSmoothSurf::~ftSmoothSurf()
  //===========================================================================
{
  // Nothing to do
}

//===========================================================================
void ftSmoothSurf::setSmoothU(int k)
{
  if (k>=0 && k<=1)
    seem_[1]=std::max(1,std::min(k,0));
}

//===========================================================================
void ftSmoothSurf::setSmoothV(int k)
{
  if (k>=0 && k<=1)
    seem_[0]=std::max(1,std::min(k,0));
}

//===========================================================================
void ftSmoothSurf::refineSurf(ftPointSet& points, bool reparam)
  //---------------------------------------------------------------------------
  //
  // Purpose: Express the surface to be modified on a refined knot vector
  //          depending on the distances between the given point and the
  //          given surface.
  //
  //===========================================================================
{
  double dist = points.getMaxDist();
  if (dist < 0)
    {
      // Compute the distance between the points and
      // the master surface
      if (reparam)
	points.computeDistAndRepar(surf_);
      else
	points.computeDist(surf_);	  
      dist =  points.getMaxDist();
    }
  if  (dist < approxtol_)
    return;   // All the points are within the tolerance, no
  // point in refining the surface.

#ifdef DEBUG
  std::ofstream of("error_points.g2");
  int nmb = points.size();
  vector<Point> out_pts;
  for (int ki=0; ki<nmb; ++ki)
    {
      ftSamplePoint *pt = points[ki];
      if (pt->getDist() > approxtol_)
	{
	  Vector3D pos = pt->getPoint();
	  out_pts.push_back(Point(pos.x(), pos.y(), pos.z()));
	}
    }
  if (out_pts.size() > 0)
    {
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << out_pts.size() << std::endl;
      for (int ki=0; ki<(int)out_pts.size(); ++ki)
	of << out_pts[ki] << std::endl;
    }
#endif
      

  // Following the notation in
  // Weiss, Andor, Renner, Varady: 'Advanced Surface Fitting Techniques' (2000),
  // on which the algorithm is based, we choose appropriate knots to insert
  // into our existing knot vectors. See paper for details.

  // The actual number of new knots to insert should be chosen analytically.
  // Increasingly as the function is called on same surface surf_?
  // @@ Hardcoded value!
  // We do not want to insert more knots then there are intervals.
  int nmb_new_knots = min(4, surf_->numCoefs_u() + 1 - surf_->order_u() +
			  surf_->numCoefs_v() + 1 - surf_->order_v());
  // To help deciding in which basis to insert a knot, we compute a unit
  // vector which indicates deviation relative to parameters u_r & v_r
  // (where the parameters correspond to the points point p_r).
  vector<double> m[2];

  // Fill m[1] and m[2] for each point p_r
  for (int kh = 0; kh < points.size(); ++kh) {
      Vector2D uv = points[kh]->getPar();
      vector<Point> derivs(6); // To store information about partial derivatives
      surf_->point(derivs, uv[0], uv[1], 2);
      double delta = points[kh]->getDist(); // Distance from point p_r to surf_

      // We define the matrices which we'll need
      ColumnVector n(surf_->dimension()); // Approximated normal vector in surf_(u,v)
      ColumnVector N(surf_->dimension()); // Actual normal vector in surf_(u,v)
      Matrix DS(surf_->dimension(), 2); // Differential form
      Matrix Id_2(2, 2); // The 2*2 identity matrix
      Matrix IS(2, 2); // First fundamental form
      Matrix IIS(2, 2); // Second fundamental form
      Matrix A_delta(2, 2); // Defined later
      ColumnVector grad_delta(2); // Gradient of distance function delta
      Matrix IS_delta(2, 2); // First FF of offset surface to S with distance delta
      Matrix mu(2, 2); // Our wanted values, before normalization

      // Initialize the matrices n & N
      Point N_r;

      try {
	  surf_->normal(N_r, uv[0], uv[1]);
      } catch (...) {
	  // @@ There should be attempted to estimate the normal using
	  // another method. To be implemented later.
	  MESSAGE("Surface must be degenerate; setting normal to 0.");
	  for (int kj = 0; kj < surf_->dimension(); ++kj)
	      N.element(kj) = 0.0;
      }

      for (int kj = 0; kj < surf_->dimension(); ++kj) {
	  n.element(kj) = points[kh]->getPoint()[kj] - derivs[0][kj];
	  N.element(kj) = N_r[kj];
      }

      // Cheking values to avoid dividing by zero later on.
      if (n.SumSquare() == 0) {
	  // In case of appr. normal vector of length 0, there is no real need to
	  // set values in m, but we do need to push_back an element.
	  m[0].push_back(0.0);
	  m[1].push_back(0.0);
      } else if (((n.t() * N).AsScalar() == 0)) {
	  // @@@ This should be handled in a somewhat better manner.
	  // For now assuming both parameters contributing.`
	  // Vector given by point and its closest point in the surface may be in the
	  // tangent plane to surface in the closest point, while not being 0 
	  // (typically this vector should be normal...). This may happen when
	  // closestPoint (for the surface) is given bad start parameters.
	  // We may also get in this situation when close to a degenerate edge.
	  MESSAGE("Trying to divide by zero. Check parametrization.");
	  m[0].push_back(0.5);
	  m[1].push_back(0.5);
      } else {
	  // Normalize the orthogonal vector n (N already normal)
	  n *= 1 / sqrt(n.SumSquare());

	  // debugging
	  //    std::cout << "(approx and exact normal vectors) <n,N>:  " <<
	  //  	      ((n.t() * N).AsScalar()) << std::endl;
	  // end debugging

	  // Initialize the matrix DS
	  for (int i = 0; i < 2; ++i)
	      for (int kj = 0; kj < surf_->dimension(); ++kj)
		  DS.element(kj, i) = derivs[i + 1][kj];

	  // Initialize the matrix id_2
	  Id_2.element(0, 0) = Id_2.element(1, 1) = 1;
	  Id_2.element(1, 0) = Id_2.element(0, 1) = 0;

	  //Initialize IS
	  IS.element(0, 0) = derivs[1] * derivs[1];
	  IS.element(1, 1) = derivs[2] * derivs[2];
	  IS.element(1, 0) = IS.element(0, 1) = derivs[1] * derivs[2];

	  //Initialize IIS
	  IIS.element(0, 0) = N_r * derivs[3];
	  IIS.element(1, 1) = N_r * derivs[5];
	  IIS.element(1, 0) = IIS.element(0, 1) = N_r * derivs[4];

	  // Now on to some real computations
	  A_delta = (Id_2 -  (IS.i()*delta) * IIS).t();
	  grad_delta = - ( A_delta * DS.t() * n ) / ((n.t() * N).AsScalar());
	  IS_delta = A_delta.t() * IS * A_delta; // Confirm IS_delta symmetric!

	  // We're then ready to compute mu, which is the vector describing local
	  // change of the distance.
	  mu = (grad_delta * grad_delta.t()) * IS_delta.i();
	  // As we're only interested in diagonal elements, we set others to 0
	  mu.element(1,0) = mu.element(0,1) = 0;

	  // Finally we may initialize m
	  m[0].push_back(mu.element(0, 0) / sqrt(mu.SumSquare()));
	  m[1].push_back(mu.element(1, 1) / sqrt(mu.SumSquare()));

	  // debugging
	  //  std::cout << m[0][h] << " " << m[1][h] << std::endl;
	  // end debugging
      }
  }

  // To store index i of points lying between given knots.
  // Indexing is the same as in the respective knot vectors.
  vector< vector<int> > index_u, index_v;
  int nmb_u_intervals = surf_->basis_u().numCoefs() + 1 - surf_->basis_u().order();
  int nmb_v_intervals = surf_->basis_v().numCoefs() + 1 - surf_->basis_v().order();
  index_u.resize(nmb_u_intervals);
  index_v.resize(nmb_v_intervals);
  // Run through the points, remember which 'iso-strip' it belongs to.
  for (int ki = 0; ki < points.size(); ++ki) {
      for (size_t kj = 0; kj < index_u.size(); ++kj) {
	  if ((points[ki]->getPar()[0] >=
	       surf_->basis_u().begin()[surf_->basis_u().order() - 1 + kj]) &&
	      (points[ki]->getPar()[0] <=
	       surf_->basis_u().begin()[surf_->basis_u().order() + kj])) {
	      // || surf_->basis_v().begin()[kj + 1] == surf_->basis_v().endparam()) {
	      index_u[kj].push_back(ki);
	      break;
	  }
      }
  }
  for (int ki = 0; ki < points.size(); ++ki) {
      for (size_t kj = 0; kj < index_v.size(); ++kj)
	  if ((points[ki]->getPar()[1] >=
	       surf_->basis_v().begin()[surf_->basis_v().order() - 1 + kj]) &&
	      (points[ki]->getPar()[1] <=
	       surf_->basis_v().begin()[surf_->basis_v().order() + kj])) {
	      // || surf_->basis_v().begin()[kj + 1] == surf_->basis_v().endparam()) {
	      index_v[kj].push_back(ki);
	      break;
	  }
  }

  // We next weight each knot interval (highest weight = insert knot),
  // based on total distance (to surf_) of all points in that strip.
  vector<double> weights_u, weights_v;
  for (size_t i = 0; i < index_u.size(); ++i) {
      weights_u.push_back(0.0); // Initializing
      for (size_t kj = 0; kj < index_u[i].size(); ++kj)
	  weights_u[i] += m[0][ index_u[i][kj] ] * m[0][ index_u[i][kj] ] *
	      points[ index_u[i][kj] ]->getDist() / approxtol_;
  }
  for (size_t i = 0; i < index_v.size(); ++i) {
      weights_v.push_back(0.0); // Initializing
      for (size_t kj = 0; kj < index_v[i].size(); ++kj)
	  weights_v[i] += m[1][ index_v[i][kj] ] * m[1][ index_v[i][kj] ] *
	      points[ index_v[i][kj] ]->getDist() / approxtol_;
  }

  vector<double> new_knots_u, new_knots_v;
  // We locate the nmb_new_knots largest weights. As we need the weight
  // vectors for this operation only, they'll be invalidated as soon as used.
  for (int i = 0; i < nmb_new_knots; ++i) {
      vector<double>::iterator iter_u_max =
	  std::max_element(weights_u.begin(), weights_u.end());
      vector<double>::iterator iter_v_max =
	  std::max_element(weights_v.begin(), weights_v.end());
      bool u_greater = (*iter_u_max >= *iter_v_max);
      int order = u_greater ? surf_->order_u() : surf_->order_v();
      vector<double>::iterator weight_begin = u_greater ?
	  weights_u.begin() : weights_v.begin();
      vector<double>::iterator max_weight_iter = u_greater ? iter_u_max : iter_v_max;
      int left_ki = order - 1 + (int)(max_weight_iter - weight_begin);
      vector<double>::const_iterator knots_begin = u_greater ?
	  surf_->basis_u().begin() : surf_->basis_v().begin();

      if (*max_weight_iter == 0) {
	// //Seemes we want to insert more knots than there are intervals.
      	//   MESSAGE("Number of knots to insert exceeds number of intervals, "
      	// 	     "should not happen.");
      	  break;
      }

      double numerator = 0;
      double denominator = 0;
      // We're now ready to compute where in the interval to insert the knot.
      int m = left_ki - order + 1;
      vector<vector<int> >::iterator index_iter =
	  (u_greater ? index_u.begin() : index_v.begin());
      for (size_t kj = 0; kj < index_iter[m].size(); ++kj) {
	  numerator += points[index_iter[m][kj]]->getPar()[!u_greater] *
	      (points[index_iter[m][kj]]->getDist() / approxtol_) *
	      (points[index_iter[m][kj]]->getDist() / approxtol_);
	  denominator += (points[index_iter[m][kj]]->getDist() / approxtol_) *
	      (points[index_iter[m][kj]]->getDist() / approxtol_);
      }

      double lambda_bar = ((numerator / denominator) - knots_begin[left_ki]) /
	  (knots_begin[left_ki+1] - knots_begin[left_ki]);
      double lambda; // Weight defining convex combination of interval end values.
      const double alpha = 0.2; // To avoid inserting too close to existing points.
      double num_tol = 1e-10; // We tolerate numerical errors.
      if ((-num_tol <= lambda_bar) && (lambda_bar < alpha))
	  lambda = alpha;
      else if ((alpha <= lambda_bar) && (lambda_bar <= 1 - alpha))
	  lambda = lambda_bar;
      else if ((1 - alpha < lambda_bar) && (lambda_bar <= 1 + num_tol))
	  lambda = 1 - alpha;
      else {
	  // We should debug if this message is displayed.
	  //       GO_ERROR("This should never happen. Most likely bug in refineSurf().",
	  // 		 UnknownError());
	  MESSAGE("This value should have been inside unit interval: " << lambda_bar);
	  if (lambda_bar < 0)
	      lambda = alpha;
	  else
	      lambda = 1.0 - alpha;
      }

      //     // debugging
      //     std::cout << numerator << ", " << denominator << ", " << knots_begin[left_ki]
      // 	      << ", " << knots_begin[left_ki+1] << ", " << std::endl;
      //     std::cout << lambda << ", " << lambda_bar << ", " << alpha << std::endl;
      //     // end of debugging

      // Finally we are ready to compute our new knot
      // VSK, 0410. Testing
      lambda = 0.5; 

      if (u_greater)
	  new_knots_u.push_back((1 - lambda) * knots_begin[left_ki]
				+ lambda * knots_begin[left_ki+1]);
      else
	  new_knots_v.push_back((1 - lambda) * knots_begin[left_ki]
				+ lambda * knots_begin[left_ki+1]);

      *max_weight_iter = 0.0; // We no longer need this specific weight.      
  }

  // Knots may be unordered. As of today this is handled by insertKnot().
  // In case of alterations we sort the vectors.
  if (new_knots_u.size() != 0) {
    std::sort(new_knots_u.begin(), new_knots_u.end());
    surf_->insertKnot_u(new_knots_u);
  }
  if (new_knots_v.size() != 0) {
    std::sort(new_knots_v.begin(), new_knots_v.end());
    surf_->insertKnot_v(new_knots_v);
  }

  // We must update coef_known_.

}

//===========================================================================
bool ftSmoothSurf::update(ftPointSet& points, double gapeps, bool reparam)
//---------------------------------------------------------------------------
//
// Purpose: Modify the surface according to the given point set or original surface..
//
{

  // Prepare for smoothing
  int smoothstatus;
  SmoothSurf smooth(0);

  // Check if there is any points to approximate
  int nmbpoints = points.size();

  // Set initial weights. THIS IS NOT A FINAL SOLUTION.
  double weight[4];
  //  double minsmooth = (approx_orig_tol_ > 0.0) ? 1e-05 : 0.001;
  double minsmooth = (approx_orig_tol_ > 0.0) ? 1e-08 : 0.001;
  //weight[3] = (nmbpoints > 0) ? init_approx_weight_ : 0.0;
  weight[3] = init_approx_weight_;
  weight[0] = 0.0;
  weight[1]  = 0.5*(1.0 - weight[3]);
  weight[2] = 0.5*(1.0 - weight[3]);

  // Fetch point data
  vector<double> pts;
  vector<double> params;
  vector<double> ptweight(nmbpoints, 1.0);
  pts.reserve(3*nmbpoints);
  params.reserve(2*nmbpoints);

  //    std::ofstream dump("pntdump.disp");
  //    dump << "list ptl" << std::endl;
  //    dump << "pnt" << std::endl;
  for (int ki=0; ki<nmbpoints; ki++)
    {
      Vector3D spacept = points[ki]->getPoint();
      for (int kj=0; kj<3; kj++)
	pts.push_back(spacept[kj]);
      //      for (kj=0; kj<3; kj++)
      //        dump << spacept[kj] << "   ";
      //      dump << std::endl;
    }
  //    dump << "pnttype 2" << std::endl;
  //    dump << "pntsize 4" << std::endl;
  //    dump << "diffuse 1.0 0.0 0.0" << std::endl;
  //    dump << "end" << std::endl;

  // Based on input, we mark those coefs which are to remain unaltered.
  vector<int> coef_known = getCoefKnown();

  // We may experience that no coef is free. No reason to continue.
  vector<int>::const_iterator min_iter = min_element(coef_known.begin(), coef_known.end());
  if (*min_iter == 1) {
     MESSAGE("No coefs were free to alter.");
     return false; // @@sbr What is the natural return value?
  }
  
//   std::cout << "Single surface (1), approximation weight : ";
//   std::cout << init_approx_weight_ << std::endl;

  int nmb_u_samples = 3*orig_surf_->numCoefs_u();
  int nmb_v_samples = 3*orig_surf_->numCoefs_v();
  double from_u = orig_surf_->startparam_u();
  double to_u = orig_surf_->endparam_u();
  double from_v = orig_surf_->startparam_v();
  double to_v = orig_surf_->endparam_v();
  double step_u = (to_u - from_u)/(nmb_u_samples-1);
  double step_v = (to_v - from_v)/(nmb_v_samples-1);
  double max_dist = 0.0, av_dist = 0.0;
  double dist, u2, v2, dist2;
  Point orig_pt, tmp_pt, clo_pt;
  double clo_max, clo_av; 
/* #if ((_MSC_VER > 0) && (_MSC_VER < 1300)) */
/*   const Domain domain = orig_surf_->parameterDomain(); */
/* #else */
/*   RectDomain domain = orig_surf_->parameterDomain(); */
/* #endif */
  double guess[2];

  // Iterate to an approximating surface
  double maxerr  = 10000000.0, meanerr  = 10000000.0;
  double prevmax = 10000000.0, prevmean = 10000000.0;
  bool isOK = false;
  int iter = 0;
  while (!isOK && iter < maxiter_) {
      // Get parameter values
      for (int ki=0; ki<nmbpoints; ki++) {
	  Vector2D parpt = points[ki]->getPar();
	  params.push_back(parpt[0]);
	  params.push_back(parpt[1]);
      }
    
#ifdef DEBUG
      std::ofstream of1("smoothsurf_init.g2");
      surf_->writeStandardHeader(of1);
      surf_->write(of1);
#endif
  
      // Make new surface
      smooth.attach(surf_, seem_, &coef_known[0]);
      smooth.setOptimize(weight[0], weight[1], weight[2]);

      // If approx_orig_tol_ is set, we iterate on surf, using spline coefs instead of points.
      if (approx_orig_tol_ > 0)
	  smooth.approxOrig(weight[3]);
      else
	  smooth.setLeastSquares(pts, params, ptweight, weight[3]);
      // smooth.approxOrig(0.000000001);

      shared_ptr<SplineSurface> tmp_surf;
      try {
	  smoothstatus = smooth.equationSolve(tmp_surf);
      } catch (...) {
	  MESSAGE("Failed solving equation!");
	  tmp_surf = surf_;
	  iter = maxiter_;
      }

//       double max_dist = 0.0;
      // Check distance to approximating points/original surface.
      bool b, r, t, l; // Variables not explicitly used, but mutable degen_ in surface is set.
      tmp_surf->isDegenerate(b, r, t, l, gapeps);
      if (nmbpoints > 0) {
	  if (reparam)
	      points.computeDistAndRepar(tmp_surf);
	  else
	      points.computeDist(tmp_surf);
	  maxerr = points.getMaxDist();
	  meanerr = points.getMeanDist();
#ifdef FANTASTIC_DEBUG
 	  std::cout << "iter: " << iter << ", max: " << maxerr;
 	  std::cout << ", mean: " << meanerr << std::endl;
 	  std::cout << "Smoothing weight: " << weight[0] << ", " << weight[1] << ", " << weight[2] << std::endl;;
 	  std::cout << ". Approx weight: " << weight[3] << std::endl;
#endif // FANTASTIC_DEBUG
	  isOK = (maxerr < approxtol_);

#ifdef FANTASTIC_DEBUG
	  std::ofstream pointsout("data/pointsdump.dat");
	  std::ofstream edgessout("data/triangedges.dat");
	  points.printXYZNodes(pointsout, true);
	  points.printXYZEdges(edgessout);
#endif // FANTASTIC_DEBUG

	  if (!(maxerr < 0.95*prevmax || meanerr < 0.95*prevmean)) {
	      if (maxerr < prevmax)
		  {
		      *surf_ = *tmp_surf;
		      max_error_ = maxerr;
		      mean_error_ = meanerr;
		  }
	      break;
	  }
      } else {
	  // We compute the distance from uniform sampled points.
	  max_dist = 0.0;
	  av_dist  = 0.0;
	  clo_max = 0.0;
	  clo_av = 0.0;
	  int ki, kj;
	  double u, v;
	  for (ki = 0, u = from_u; ki < nmb_u_samples; ++ki, u+=step_u)
	      for (kj = 0, v = from_v; kj < nmb_v_samples; ++kj, v+=step_v) 
		  {
		      orig_pt = orig_surf_->ParamSurface::point(u, v);
		      tmp_pt = tmp_surf->ParamSurface::point(u, v);
		      guess[0] = u;
		      guess[1] = v;
		      /* 		tmp_surf->closestPoint(orig_pt, u2, v2, clo_pt, */
		      /* 				       dist2, 0.000001, domain, */
		      /* 				       guess); */
		      tmp_surf->closestPoint(orig_pt, u2, v2, clo_pt, 
					     dist2, 0.000001, NULL,
					     guess);
		      dist = (tmp_pt - orig_pt).length();
		      max_dist = max(dist, max_dist);
		      av_dist += dist;
		      clo_max = max(clo_max, dist2);
		      clo_av += dist2;
		  }
	  av_dist /= (double)(nmb_u_samples*nmb_v_samples);
	  clo_av /= (double)(nmb_u_samples*nmb_v_samples);
	  //  	  std::cout << "Iteration : " << iter << " Max distance : " << max_dist; 
	  //  	  std::cout << " Average distance : " << av_dist << std::endl;
	  //  	  std::cout << "Closest max : " << clo_max;
	  //  	  std::cout << " Average max : " << clo_av << std::endl;

	  maxerr = max_dist;
	  meanerr = av_dist;
	  isOK = (max_dist < approx_orig_tol_);
	  //  	  if (!(maxerr < 0.95*prevmax || meanerr < 0.95*prevmean)) {
	  //  	      if (max_dist < prevmax)
	  //  		{
	  //  		  *surf_ = *tmp_surf;
	  //  		  max_error_ = maxerr;
	  //  		  mean_error_ = meanerr;
	  //  		}
	  //  	      break;
	  //      }
	  // 	  // We compute the max distance between corresponding control points.
	  // 	  int dim = surf_->dimension();
	  // 	  int num_coefs = surf_->numCoefs_u()*surf_->numCoefs_v();
	  //  	  double max_dist = 0.0;
	  // 	  for (ki = 0; ki < num_coefs; ++ki) {
	  // 	      double dist = 0.0;
	  // 	      for (kj = 0; kj < dim; ++kj)
	  // 		  dist += (orig_surf_->coefs_begin()[ki*dim+kj] - tmp_surf->coefs_begin()[ki*dim+kj])*
	  // 		      (orig_surf_->coefs_begin()[ki*dim+kj] - tmp_surf->coefs_begin()[ki*dim+kj]);
	  // 	      dist = sqrt(dist);
	  // 	      max_dist = max(dist, max_dist);
	  // 	  }
      }

      ++iter;
      if (!isOK && iter<maxiter_) {
	  // Iterate on approximation weight
	  weight[3] = 0.5*(weight[3] + 1.0);
	  weight[1] = 1.0 - weight[3];
	  //weight[0] = 0.001*(1.0 - weight[3]);
      }

#ifdef DEBUG
      std::ofstream of2("smoothsurf.g2");
      tmp_surf->writeStandardHeader(of2);
      tmp_surf->write(of2);
#endif

       if (isOK || nmbpoints > 0)
	  {
	    //*surf_ = *tmp_surf;
	    surf_->swap(*tmp_surf);
	    tmp_surf = surf_;
	      max_error_ = maxerr;
	      mean_error_ = meanerr;
#ifdef DEBUG
	      surf_->writeStandardHeader(of2);
	      surf_->write(of2);
#endif
 	  }
      else
	  *surf_ = *orig_surf_;

      if (weight[1] < minsmooth)
	  break;

      prevmax = maxerr;
      prevmean = meanerr;
  }
  //    std::cout << "Size: " << surf_->numCoefs_u() << ", ";
  //    std::cout << surf_->numCoefs_v() << std::endl;

  init_approx_weight_ = weight[3];
  //   std::cout << "Single surface (2), approximation weight : ";
  //   std::cout << init_approx_weight_ << std::endl;
  return isOK;
}


//===========================================================================
vector<int> ftSmoothSurf::getCoefKnown()
//===========================================================================
{
  int ki, kj;
  // From the information given by ccw_edge_derivs_, we set boundary conditions.
  int kn1 = surf_->numCoefs_u();
  int kn2 = surf_->numCoefs_v();
  vector<int> coef_known(kn1*kn2);
  fill(coef_known.begin(), coef_known.end(), 0);

  for (ki = 0; ki < ccw_edge_derivs_[0]; ++ki)
    for (kj = 0; kj < kn1; ++kj)
      coef_known[kn1*ki+kj] = 1;

  for (ki = 0; ki < ccw_edge_derivs_[1]; ++ki)
    for (kj = 0; kj < kn2; ++kj)
      coef_known[(kj+1)*kn1-1-ki] = 1;

  for (ki = 0; ki < ccw_edge_derivs_[2]; ++ki)
    for (kj = 0; kj < kn1; ++kj)
      coef_known[(kn2-1-ki)*kn1+kj] = 1;

  for (ki = 0; ki < ccw_edge_derivs_[3]; ++ki)
    for (kj = 0; kj < kn2; ++kj)
      coef_known[kj*kn1+ki] = 1;

  // We lock coefs denoting corner points.
  if (lock_corner_points_) {
    coef_known[0] = 1;
    coef_known[kn1-1] = 1;
    coef_known[(kn2-1)*kn1] = 1;
    coef_known[kn1*kn2-1] = 1;
  }

  return coef_known;
}

} // namespace Go
