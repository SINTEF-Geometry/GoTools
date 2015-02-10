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

#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/creators/SolveCG.h"

using namespace Go;
using std::vector;

namespace {

const int indices[] = {1, 2, 3, 4, 5, 8};

const double w_2_0 = 0.5555555556;
const double w_2_1 = 0.8888888889;
const double w_3_0 = 0.3478548451;
const double w_3_1 = 0.6521451549;
const double w_4_0 = 0.2369268851;
const double w_4_1 = 0.4786286705;
const double w_4_2 = 0.5688888889;
const double w_5_0 = 0.1012285363;
const double w_5_1 = 0.2223810345;
const double w_5_2 = 0.3137066459;
const double w_5_3 = 0.3626837834;
    
const double weight[][8] =
    { {   1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      {   1.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { w_2_0,  w_2_1,  w_2_0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { w_3_0,  w_3_1,  w_3_1,  w_3_0,    0.0,    0.0,    0.0,    0.0},
      { w_4_0,  w_4_1,  w_4_2,  w_4_1,  w_4_0,    0.0,    0.0,    0.0},
      { w_5_0,  w_5_1,  w_5_2,  w_5_3,  w_5_3,  w_5_2,  w_5_1,  w_5_0} };

const double s_0_0 = 0.5;
const double s_1_0 = -0.5773502692;
const double s_2_0 = -0.7745966692;
const double s_3_0 = -0.8611363116;
const double s_3_1 = -0.3399810436;
const double s_4_0 = -0.9061798459;
const double s_4_1 = -0.5384693101;
const double s_5_0 = -0.9602898565;
const double s_5_1 = -0.7966664774;
const double s_5_2 = -0.5255324099;
const double s_5_3 = -0.1834346425;

const double sample[][8] =
    { { s_0_0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_1_0, -s_1_0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_2_0,    0.0, -s_2_0,    0.0,    0.0,    0.0,    0.0,    0.0},
      { s_3_0,  s_3_1, -s_3_1, -s_3_0,    0.0,    0.0,    0.0,    0.0},
      { s_4_0,  s_4_1,    0.0, -s_4_1, -s_4_0,    0.0,    0.0,    0.0},
      { s_5_0,  s_5_1,  s_5_2,  s_5_3, -s_5_3, -s_5_2, -s_5_1, -s_5_0} };


}; // anonymous namespace containing the numerical variables needed for Gauss integration

  //============================================================================
  // Constructs a map from a B-spline to a simple Index.
  // ============================================================================
  BsplineIndexMap construct_approx_bsplineindex_map( const LRSplineSurface & lrs )
  {
    BsplineIndexMap result;
    size_t i = 0;
    for (LRSplineSurface::BSplineMap::const_iterator it_bs=lrs.basisFunctionsBegin(); 
        it_bs!=lrs.basisFunctionsEnd(); ++it_bs)
      {
	if (it_bs->second->coefFixed())
	  continue;
	
	result[it_bs->second.get()] = i++;
      }
    return result;
  }


//==============================================================================
LRSurfSmoothLS::LRSurfSmoothLS(shared_ptr<LRSplineSurface> surf, vector<int>& coef_known)
//==============================================================================
  : srf_(surf), coef_known_(coef_known)
{
  // Distribute information about fixed coefficients to the B-splines
  ncond_ = 0;
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      it_bs->second->setFixCoef(coef_known[ki]);
      if (coef_known[ki] == 0)
	ncond_++;
    }
  
  // Construct index map
  BSmap_ = construct_approx_bsplineindex_map(*srf_);

  // Allocate scratch for equation system
  gmat_.assign(ncond_*ncond_, 0.0);
  gright_.assign(srf_->dimension()*ncond_, 0.0);
  
}

//==============================================================================
LRSurfSmoothLS::LRSurfSmoothLS()
//==============================================================================
{
}

//==============================================================================
void LRSurfSmoothLS::setInitSf(shared_ptr<LRSplineSurface> surf, 
			       vector<int>& coef_known)
//==============================================================================
{
  srf_ = surf;
  coef_known_ = coef_known;

  // Distribute information about fixed coefficients to the B-splines
  ncond_ = 0;
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      it_bs->second->setFixCoef(coef_known[ki]);
      if (coef_known[ki] == 0)
	ncond_++;
    }
  
  // Construct index map
  BSmap_ = construct_approx_bsplineindex_map(*srf_);

  // Allocate scratch for equation system
  gmat_.assign(ncond_*ncond_, 0.0);
  gright_.assign(srf_->dimension()*ncond_, 0.0);
  
}

//==============================================================================
LRSurfSmoothLS::~LRSurfSmoothLS()
//==============================================================================
{
}

//==============================================================================
void LRSurfSmoothLS::updateLocals()
//==============================================================================
{
  ncond_ = 0;
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      if (!it_bs->second->coefFixed())
	ncond_++;
    }

  BSmap_ = construct_approx_bsplineindex_map(*srf_);

  gmat_.assign(ncond_*ncond_, 0.0);
  gright_.assign(srf_->dimension()*ncond_, 0.0);
}

//==============================================================================
bool LRSurfSmoothLS::hasDataPoints() const
//==============================================================================
{
  // For each element, check if any data points are stored
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasDataPoints())
	return true;   // Data points are found associated to the current element
    }
  return false;  // No data points are found
}

//==============================================================================
void LRSurfSmoothLS::addDataPoints(vector<double>& points, bool is_ghost_points) 
//==============================================================================
{
  LRSplineUtils::distributeDataPoints(srf_.get(), points, true, !is_ghost_points);
}

//==============================================================================
void LRSurfSmoothLS::setOptimize(const double weight1, const double weight2,
				 const double weight3)
//==============================================================================
{

  int dim = srf_->dimension();
  double eps = 1.0e-10;  // Numerical tolerance
  int der1 = (weight1 > eps) ? 1 : 0;
  int der2 = (weight2 > eps) ? 1 : 0;
  int der3 = (weight3 > eps) ? 1 : 0;

  if (der1 + der2 + der3 == 0)
    return;   // No smoothing applyed. Nothing to do.

  // Perform Bezier extraction. Not implemented yet

  // For each element
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      // For all B-splines in the support of the element
      // Compute integrals of inner products of derivatives of the B-spline
      
      // Fetch B-splines
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
      size_t nmb = bsplines.size();

      // Fetch derivative of B-splines in the Gauss points
      // Store only those entries which are used in the computations
      vector<double> basis_derivs;
      int nmbGauss;
      fetchBasisDerivs(bsplines, basis_derivs, der1, der2, der3, 
		       it->second->umin(), it->second->umax(),
		       it->second->vmin(), it->second->vmax(), nmbGauss);

      if (der1)
	{
	  // Compute contribution of integrals of d_u^2 and d_v^2
	  computeDer1Integrals(bsplines, nmbGauss, &basis_derivs[0], weight1);
	}
			       
      if (der2)
	{
	  // Compute contribution of integrals of d_uu^2, d_uv^2, d_vv^2
	  // and d_uu*d_vv
	  int idx = (der1) ? 2*bsplines.size()*nmbGauss : 0;
	  computeDer2Integrals(bsplines, nmbGauss, &basis_derivs[idx], weight2);
	}

      if (der3)
	{
	  // Compute contribution of integrals of d_uuu^2, d_uuv^2, d_uvv^2
	  // d_vvv^2 d_uuu*d_uvv and d_uuv*d_vvv
	  int idx = (der1) ? 2*bsplines.size()*nmbGauss : 0;
	  if (der2)
	    idx += 3*bsplines.size()*nmbGauss;
	  computeDer3Integrals(bsplines, nmbGauss, &basis_derivs[idx], weight3);
	}

    }
 }

//==============================================================================
void LRSurfSmoothLS::smoothBoundary(const double weight1, const double weight2,
				    const double weight3)
//==============================================================================
{

  int dim = srf_->dimension();
  double eps = 1.0e-10;  // Numerical tolerance
  int der1 = (weight1 > eps) ? 1 : 0;
  int der2 = (weight2 > eps) ? 1 : 0;
  int der3 = (weight3 > eps) ? 1 : 0;

  if (der1 + der2 + der3 == 0)
    return;   // No smoothing applyed. Nothing to do.

  // For each boundary
  Direction2D d;
  int ki;
  for (d=XFIXED, ki=0; ki<2; d=YFIXED, ++ki)
    {
      int kj;
      bool atstart;
      for (atstart=true, kj=0; kj<2; atstart=false, ++kj)
	{
	  // Fetch all non-zero B-splines along the current boundary
	  // NB! This functionality requires maximum multiplicity knots at
	  // the boundaries
	  vector<LRBSpline2D*> bsplines = srf_->getBoundaryBsplines(d, atstart);

	  // Fetch knots related to this boundary curve
	  int ix = srf_->mesh().numDistinctKnots(d) - 1;
	  vector<double> knots = srf_->mesh().getKnots(flip(d), ix);

	  for (size_t kr=1; kr<knots.size(); ++kr)
	    {
	      // Extract element boundaries and B-splines
	      double tmin = knots[kr-1];
	      double tmax = knots[kr];
	      if (tmax <= tmin)
		continue;

	      vector<LRBSpline2D*> bsplines_el = 
		bsplinesCoveringElement(bsplines, d, tmin, tmax);
	      if (bsplines_el.size() == 0)
		continue;

	      // Compute integrals of inner products of derivatives of B-splines
	      vector<double> basis_derivs;
	      int nmbGauss;
	      fetchBasisLineDerivs(bsplines_el, basis_derivs, der1, der2, der3, 
				   flip(d), tmin, tmax, nmbGauss);

	      if (der1)
		{
		  // Compute contribution of integrals of d_t^2
		  computeDer1LineIntegrals(bsplines_el, nmbGauss, 
					   &basis_derivs[0], weight1);
		}
			       
	      if (der2)
		{
		  // Compute contribution of integrals of d_tt^2
		  int idx = (der1) ? bsplines_el.size()*nmbGauss : 0;
		  computeDer2LineIntegrals(bsplines_el, nmbGauss, 
					   &basis_derivs[idx], weight2);
		}

	      if (der3)
		{
		  // Compute contribution of integrals of d_ttt^2
		  int idx = (der1) ? bsplines_el.size()*nmbGauss : 0;
		  if (der2)
		    idx += bsplines_el.size()*nmbGauss;
		  computeDer3LineIntegrals(bsplines_el, nmbGauss, 
					   &basis_derivs[idx], weight3);
		}
	    }
	  
	}
    }

}

//==============================================================================
void LRSurfSmoothLS::setLeastSquares(const double weight)
//==============================================================================
{
  int dim = srf_->dimension();

  // For each element
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      // Check if the element contains an associated least squares matrix
      bool has_LS_mat = it->second->hasLSMatrix();

      // Check if the element is changed
      bool is_modified = it->second->isModified();

      // Fetch B-splines
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
      size_t nmb = bsplines.size();

      if (!has_LS_mat || is_modified)
	{
	  // Either no pre-computed least squares matrix exists or 
	  // the element or an associated B-spline is changed.
	  // Compute the least squares matrix associated to the 
	  // element
	  // First fetch data points
	  vector<double>& elem_data = it->second->getDataPoints();

	  // Fetch ghost points (points that are included to stabilize
	  // the computation, but are not tested for accuracy
	  vector<double>& ghost_points = it->second->getGhostPoints();

	  // Compute sub matrix
	  // First get access to storage in the element
	  double *subLSmat, *subLSright;
	  int kcond;
 	  it->second->setLSMatrix();
	  it->second->getLSMatrix(subLSmat, subLSright, kcond);
 
	  localLeastSquares(elem_data, ghost_points,
			    bsplines, subLSmat, subLSright, kcond);
	  int stop_break = 1;
	}

      // Assemble stiffness matrix and right hand side based on the local least 
      // squares matrix
      // The size of the stiffness matrix is the squared number of LR B-splines
      // with a free coefficient. The size of the right hand side is equal to
      // the number of free coefficients times the dimension of the data points
      double *subLSmat, *subLSright;
      int kcond;
      it->second->getLSMatrix(subLSmat, subLSright, kcond);

      vector<size_t> in_bs(kcond);
      size_t ki, kj, kr, kh;
      for (ki=0, kj=0; ki<nmb; ++ki)
	{
	  if (bsplines[ki]->coefFixed())
	    continue;

	  // Fetch index in the stiffness matrix
	  size_t inb = BSmap_.at(bsplines[ki]);
	  in_bs[kj++] = inb;
	}

      int kk;
      for (ki=0, kr=0; ki<nmb; ++ki)
	{
	  if (bsplines[ki]->coefFixed())
	    continue;
	  size_t inb1 = in_bs[kr];
	  for (kk=0; kk<dim; ++kk)
	    gright_[kk*ncond_+inb1] += weight*subLSright[kk*kcond+kr];
	  for (kj=0, kh=0; kj<nmb; ++kj)
	    {
	      if (bsplines[kj]->coefFixed())
		continue;
	      gmat_[inb1*ncond_+in_bs[kh]] += weight*subLSmat[kr*kcond+kh];
	      kh++;
	    }
	  kr++;
	}
    }

}

//==============================================================================
void LRSurfSmoothLS::setLeastSquares(vector<double>& points, const double weight)
//==============================================================================
{
  addDataPoints(points);
  setLeastSquares(weight);
}

//==============================================================================
int
LRSurfSmoothLS::equationSolve(shared_ptr<LRSplineSurface>& surf)
//==============================================================================
{
  int kstat = 0;
  int dim = srf_->dimension();
  int ki, kk;

  // Solve the equation system by Conjugate Gradient Method.

  vector<double> eb(dim*ncond_, 0.0);
  for (ki=0; ki<dim*ncond_; ki++)
    eb[ki] = gright_[ki];

  // Copy coefficients to array of unknowns
  LRSplineSurface::BSplineMap::const_iterator it_bs;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; 
       it_bs!=srf_->basisFunctionsEnd(); ++it_bs)
    {
      if (it_bs->second->coefFixed())
	continue;
      
      Point cf = it_bs->second->coefTimesGamma();
      for (kk=0; kk<dim; kk++)
	gright_[kk*ncond_+ki] = cf[kk]; 
      ki++;
    }
       
  // Set up CG-object

  SolveCG solveCg;

  // Create sparse matrix.

  ASSERT(gmat_.size() > 0);
  solveCg.attachMatrix(&gmat_[0], ncond_);

  // Attach parameters.

  solveCg.setTolerance(0.00000001);
  // Preconditioning
  // @@sbr When using the side-constraints our system is not longer guaranteed
  //       to be symmetric and positive definite!
  //       We should then use a different solver.
  int precond = 1;
  int nmb_iter = precond ? ncond_ : 2*ncond_;
  solveCg.setMaxIterations(std::min(nmb_iter, 1000));
  if (precond) {
    double omega = 0.1;
    // 	   printf("Omega = ");
    // 	   scanf("%lf",&omega);
    solveCg.precondRILU(omega);
  }

  // Solve equation systems.
       
  for (kk=0; kk<dim; kk++)
    {
      kstat = solveCg.solve(&gright_[kk*ncond_], &eb[kk*ncond_],
			    ncond_);
      //	       printf("solveCg.solve status %d \n", kstat);
      if (kstat < 0)
	return kstat;
      if (kstat == 1)
	THROW("Failed solving system (within tolerance)!");
    }

  // Update coefficients
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; 
       it_bs!=srf_->basisFunctionsEnd(); ++it_bs)
    {
      if (it_bs->second->coefFixed())
	continue;
      
      Point cf(dim);
      for (kk=0; kk<dim; kk++)
	cf[kk] = gright_[kk*ncond_+ki];
      srf_->setCoef(cf, it_bs->second.get());
      ki++;
    }


  // Output surface
  surf = srf_;

  return 0;
}

//==============================================================================
void LRSurfSmoothLS::localLeastSquares(vector<double>& points,
				       vector<double>& ghost_points,
				       const vector<LRBSpline2D*>& bsplines,
				       double* mat, double* right, int ncond)
//==============================================================================
{
  size_t nmbb = bsplines.size();
  int dim = srf_->dimension();
  int del = dim+3;  // Parameter pair, point and distance storage
  int nmbp[2];
  nmbp[0] = (int)points.size()/del;
  nmbp[1] = (int)ghost_points.size()/del;
  double* start_pt[2];
  start_pt[0] = &points[0];
  start_pt[1] = &ghost_points[0];

  size_t ki, kj, kp, kq, kr, kk;
  double *pp;
  for (int ptype=0; ptype<2; ++ptype)
    for (kr=0, pp=start_pt[ptype]; kr<nmbp[ptype]; ++kr, pp+=del)
      {
	vector<double> sb = getBasisValues(bsplines, pp);
	for (ki=0, kj=0; ki<nmbb; ++ki)
	  {
	    if (bsplines[ki]->coefFixed())
	      continue;
	    double gamma1 = bsplines[ki]->gamma();
	    for (kk=0; kk<dim; ++kk)
	      right[kk*ncond+kj] += gamma1*pp[2+kk]*sb[ki];
	    for (kp=0, kq=0; kp<nmbb; kp++)
	      {
		int fixed = bsplines[kp]->coefFixed();
		if (fixed == 2)
		  continue;

		double gamma2 = bsplines[kp]->gamma();
		double val = gamma1*gamma2*sb[ki]*sb[kp];
		if (fixed == 1)
		  {
		    // Move contribution to the right hand side
		    const Point coef = bsplines[kp]->Coef();
		    for (kk=0; kk<dim; ++kk)
		      right[kk*ncond+kj] -= coef[kk]*val;
		  }
		else
		  {
		    mat[kq*ncond+kj] += val;
		    kq++;
		  }
	      }
	    kj++;
	  }
      }
}

//==============================================================================
vector<double> LRSurfSmoothLS::getBasisValues(const vector<LRBSpline2D*>& bsplines,
					      double *par)
//==============================================================================
{
  vector<double> bs(bsplines.size());
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      const bool u_on_end = (par[0] == bsplines[ki]->umax());
      const bool v_on_end = (par[1] == bsplines[ki]->vmax());
      bs[ki] = bsplines[ki]->evalBasisFunction(par[0], par[1], 0, 0, 
					       u_on_end, v_on_end);
    }
  return bs;
}

//==============================================================================
void LRSurfSmoothLS::fetchBasisDerivs(const vector<LRBSpline2D*>& bsplines, 
				      vector<double>& basis_derivs, 
				      int der1, int der2, int der3, 
				      double umin, double umax,
				      double vmin, double vmax, int& nmbGauss)
//==============================================================================
{
  // Note. This function will be rewritten when Bezier extraction is introduced

  // Note that rational surface are not handled.

  if (bsplines.size() == 0)
    return;  // Nothing to do
  int bsize = (int)bsplines.size();

  // Number of Gauss points
  int deg1 = bsplines[0]->degree(XFIXED);
  int deg2 = bsplines[0]->degree(YFIXED);
  int ix1 = std::min(deg1, 5);
  int ix2 = std::min(deg2, 5);
  int wgs1 = indices[ix1];
  int wgs2 = indices[ix2];
  nmbGauss = wgs1*wgs2;

  // Storage for parameters corresponding to the Gauss points
  // Define the Gauss points in the two parameter directions
  int kj;
  vector<double> gausspar1(wgs1);
  vector<double> gausspar2(wgs2);
  for (kj=0; kj<wgs1; ++kj)
    gausspar1[kj] = 0.5*(sample[ix1][kj]*(umax-umin) + umax + umin);
  for (kj=0; kj<wgs2; ++kj)
    gausspar2[kj] = 0.5*(sample[ix1][kj]*(vmax-vmin) + vmax + vmin);

  // Allocate scratch for the results of the basis evaluation. Store only those
  // entries that will be used
  int nmb = 2*der1 + 3*der2 + 4*der3;  // Number of derivatives stored
  basis_derivs.resize(nmb*bsize*nmbGauss);

  // To simplify the use of the evaluations, the derivative is the prior sequencing
  // category, the B-spline is the next and the Gauss points the last one. The
  // derivatives are stored in the following sequence: du, dv, duu, duv, dvv,
  // duuu, duuv, duvv and dvvv. Only the specified derivatives are stored. Gauss
  // points run fastest in the u-direction.

  // Number of derivatives to compute
  int nmb_der = (der3) ? 3 : ((der2) ? 2 : 1);

  // For all bsplines
  for (int ki=0; ki<bsize; ++ki)
    {
      // Compute all relevant derivatives in all Gauss points
      vector<double> derivs;  // Storage for all derivatives in
      // all points. Sequence: du for all points, then dv, duu, duv, dvv, ...
      // The position of the basis function is NOT stored.
      bsplines[ki]->evalBasisGridDer(nmb_der, gausspar1, gausspar2,
				     derivs);
      
      // Transfer result to the output array
      int curr = 0;
      if (der1)
	{
	  std::copy(derivs.begin(), derivs.begin()+nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	  std::copy(derivs.begin()+nmbGauss, derivs.begin()+2*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	}			
      if (der2)
	{
	  std::copy(derivs.begin()+2*nmbGauss, derivs.begin()+3*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	  std::copy(derivs.begin()+3*nmbGauss, derivs.begin()+4*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	  std::copy(derivs.begin()+4*nmbGauss, derivs.begin()+5*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	}			
      if (der3)
	{
	  std::copy(derivs.begin()+5*nmbGauss, derivs.begin()+6*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	  std::copy(derivs.begin()+6*nmbGauss, derivs.begin()+7*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	  std::copy(derivs.begin()+7*nmbGauss, derivs.begin()+8*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	}
     }
}

//==============================================================================
void LRSurfSmoothLS::fetchBasisLineDerivs(const vector<LRBSpline2D*>& bsplines, 
					  vector<double>& basis_derivs, 
					  int der1, int der2, int der3, 
					  Direction2D d, double tmin, 
					  double tmax, int& nmbGauss)
//==============================================================================
{
  // Note. This function will be rewritten when Bezier extraction is introduced

  // Note that rational surface are not handled.

  if (bsplines.size() == 0)
    return;  // Nothing to do
  int bsize = (int)bsplines.size();

  // Number of Gauss points
  int deg = bsplines[0]->degree(d);
  int ix = std::min(deg, 5);
  nmbGauss = indices[ix];

  // Storage for parameters corresponding to the Gauss points
  // Define the Gauss points in the two parameter directions
  int kj;
  vector<double> gausspar(nmbGauss);
  for (kj=0; kj<nmbGauss; ++kj)
    gausspar[kj] = 0.5*(sample[ix][kj]*(tmax-tmin) + tmax + tmin);

  // Allocate scratch for the results of the basis evaluation. Store only those
  // entries that will be used
  int nmb = der1 + der2 + der3;  // Number of derivatives stored
  basis_derivs.resize(nmb*bsize*nmbGauss);

  // To simplify the use of the evaluations, the derivative is the prior sequencing
  // category, the B-spline is the next and the Gauss points the last one. The
  // derivatives are stored in the following sequence: dt, dtt, dttt. Only the specified 
  // derivatives are stored. 

  // Number of derivatives to compute
  int nmb_der = (der3) ? 3 : ((der2) ? 2 : 1);

  // For all bsplines
  for (int ki=0; ki<bsize; ++ki)
    {
      // Compute all relevant derivatives in all Gauss points
      vector<double> derivs;  // Storage for all derivatives in
      // all points. Sequence: du for all points, then dv, duu, duv, dvv, ...
      // The position of the basis function is NOT stored.
      bsplines[ki]->evalBasisLineDer(nmb_der, d, gausspar, derivs);
      
      // Transfer result to the output array
      int curr = 0;
      if (der1)
	{
	  std::copy(derivs.begin(), derivs.begin()+nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	}			
      if (der2)
	{
	  std::copy(derivs.begin()+nmbGauss, derivs.begin()+2*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	  curr += bsize;
	}			
      if (der3)
	{
	  std::copy(derivs.begin()+2*nmbGauss, derivs.begin()+3*nmbGauss,
		    basis_derivs.begin()+(curr+ki)*nmbGauss);
	}
     }
}

//==============================================================================
void LRSurfSmoothLS::computeDer1Integrals(const vector<LRBSpline2D*>& bsplines, 
					  int nmbGauss, double* basis_derivs, 
					  double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
      size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double dudu = 0.0; // d_u^2
	  double dvdv = 0.0; // d_v^2
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      dudu += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	      dvdv += basis_derivs[nmbder+ki*nmbGauss+kr]*
		basis_derivs[nmbder+kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*(dudu + dvdv);
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
void LRSurfSmoothLS::computeDer1LineIntegrals(const vector<LRBSpline2D*>& bsplines, 
					      int nmbGauss, double* basis_derivs, 
					      double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
      size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double dtdt = 0.0; // d_t^2
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      dtdt += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*dtdt;
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
void LRSurfSmoothLS::computeDer2Integrals(const vector<LRBSpline2D*>& bsplines, 
					  int nmbGauss, double* basis_derivs, 
					  double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
      size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double duuduu = 0.0; // d_uu^2
	  double dvvdvv = 0.0; // d_vv^2
	  double duvduv = 0.0; // d_uv^2
	  double duudvv = 0.0; // d_uu*d_vv
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      duuduu += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	      dvvdvv += basis_derivs[2*nmbder+ki*nmbGauss+kr]*
		basis_derivs[2*nmbder+kj*nmbGauss+kr];
	      duvduv += basis_derivs[nmbder+ki*nmbGauss+kr]*
		basis_derivs[nmbder+kj*nmbGauss+kr];
	      duudvv += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[2*nmbder+kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*(3.0*(duuduu + dvvdvv) + 4.0*duvduv + 
					     2.0*duudvv);
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
void LRSurfSmoothLS::computeDer2LineIntegrals(const vector<LRBSpline2D*>& bsplines, 
					      int nmbGauss, double* basis_derivs, 
					      double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
      size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double dttdtt = 0.0; // d_tt^2
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      dttdtt += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*dttdtt;
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
void LRSurfSmoothLS::computeDer3Integrals(const vector<LRBSpline2D*>& bsplines, 
					  int nmbGauss, double* basis_derivs, 
					  double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
       size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double duuuduuu = 0.0; // d_uuu^2
	  double dvvvdvvv = 0.0; // d_vvv^2
	  double duuvduuv = 0.0; // d_uuv^2
	  double duvvduvv = 0.0; // d_uvv^2
	  double duuuduvv = 0.0; // d_uuu*d_uvv
	  double duuvdvvv = 0.0; // d_uuv*d_vvv
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      duuuduuu += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	      dvvvdvvv += basis_derivs[3*nmbder+ki*nmbGauss+kr]*
		basis_derivs[3*nmbder+kj*nmbGauss+kr];
	      duuvduuv += basis_derivs[nmbder+ki*nmbGauss+kr]*
		basis_derivs[nmbder+kj*nmbGauss+kr];
	      duvvduvv += basis_derivs[2*nmbder+ki*nmbGauss+kr]*
		basis_derivs[2*nmbder+kj*nmbGauss+kr];
	      duuuduvv += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[2*nmbder+kj*nmbGauss+kr];
	      duuvdvvv += basis_derivs[nmbder+ki*nmbGauss+kr]*
		basis_derivs[3*nmbder+kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*(5.0*(duuuduuu + dvvvdvvv) + 
					     9.0*(duuvduuv + duvvduvv) + 
					     6.0*(duuuduvv + duuvdvvv));
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
void LRSurfSmoothLS::computeDer3LineIntegrals(const vector<LRBSpline2D*>& bsplines, 
					      int nmbGauss, double* basis_derivs, 
					      double weight)
//==============================================================================
{
  int dim = srf_->dimension();
  int nmbder = (int)bsplines.size()*nmbGauss;  // Number of entries for each derivative
  size_t ki, kj;
  for (ki=0; ki<bsplines.size(); ++ki)
    {
      if (bsplines[ki]->coefFixed())
	continue;
      double gamma1 = bsplines[ki]->gamma();
       size_t ix1 = BSmap_.at(bsplines[ki]); // Index in stiffness matrix
      for (kj=ki; kj<bsplines.size(); ++kj)
	{
	  int coef_fixed = bsplines[kj]->coefFixed();
	  if (coef_fixed == 2)
	    continue;
	  double gamma2 = bsplines[kj]->gamma();
	  size_t ix2;
	  if (!coef_fixed)
	    ix2 = BSmap_.at(bsplines[kj]);

	  double dtttdttt = 0.0; // d_ttt^2
	  for (int kr=0; kr<nmbGauss; ++kr)
	    {
	      dtttdttt += basis_derivs[ki*nmbGauss+kr]*
		basis_derivs[kj*nmbGauss+kr];
	    }

	  double val = weight*gamma1*gamma2*dtttdttt;
	  if (coef_fixed)
	    {
	      // Add contribution to the right side of the equation system
	      for (int kk=0; kk<dim; ++kk)
		gright_[kk*ncond_+ix1] -= val;
	    }
	  else
	    {
	      // Add contribution to the stiffness matrix
	      gmat_[ix1*ncond_+ix2] += val;
	      if (ki != kj)
		gmat_[ix2*ncond_+ix1] += val;
	    }
	}
    }
}

//==============================================================================
vector<LRBSpline2D*>  
LRSurfSmoothLS::bsplinesCoveringElement(std::vector<LRBSpline2D*>& cand, 
					Direction2D d, double tmin, double tmax)
//==============================================================================
{
  vector<LRBSpline2D*> bsplines;
  for (size_t ki=0; ki<cand.size(); ++ki)
    {
      double t1 = (d == XFIXED) ? cand[ki]->vmin() : cand[ki]->umin();
      double t2 = (d == XFIXED) ? cand[ki]->vmax() : cand[ki]->umax();
      if (t1 >= tmax || t2 <= tmin)
	continue;
      bsplines.push_back(cand[ki]);
    }
  return bsplines;
}
