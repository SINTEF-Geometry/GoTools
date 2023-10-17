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

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineUtils.h"

using namespace std;

namespace Go
{

  namespace
  {
    /// Functor that scales the input argument.
    class ScaleBy //: public std::unary_function<double, double>
    {
      double m_scale;

    public:
      ScaleBy(const double& scale) : m_scale(scale) {}

      double operator()(const double& value) { return m_scale * value; }
    };
  } // anonymous namespace

void volume_ratder(double const eder[],int idim,int ider,double gder[]);



//===========================================================================
void  SplineVolume::point(Point& pt, double upar, double vpar, double wpar) const
//===========================================================================
{
    pt.resize(dim_);
    const int uorder = order(0);
    const int vorder = order(1);
    const int worder = order(2);
    const int unum = numCoefs(0);
    const int vnum = numCoefs(1);
    int kdim = rational_ ? dim_ + 1 : dim_;

    static ScratchVect<double, 10> Bu(uorder);
    static ScratchVect<double, 10> Bv(vorder);
    static ScratchVect<double, 10> Bw(worder);
    static ScratchVect<double, 4> tempPt(kdim);
    static ScratchVect<double, 4> tempPt2(kdim);
    static ScratchVect<double, 4> tempResult(kdim);

    Bu.resize(uorder);
    Bv.resize(vorder);
    Bw.resize(worder);
    tempPt.resize(kdim);
    tempPt2.resize(kdim);
    tempResult.resize(kdim);

    // compute tbe basis values and get some data about the spline spaces
    basis_u_.computeBasisValues(upar, Bu.begin());
    basis_v_.computeBasisValues(vpar, Bv.begin());
    basis_w_.computeBasisValues(wpar, Bw.begin());
    const int uleft = basis_u_.lastKnotInterval();
    const int vleft = basis_v_.lastKnotInterval();
    const int wleft = basis_w_.lastKnotInterval();
    
    // compute the tensor product value
    const int start_ix =  (uleft - uorder + 1 + unum * (vleft - vorder + 1 + vnum * (wleft - worder + 1))) * kdim;

    double* ptemp;
    const double* co_ptr = rational_ ? &rcoefs_[start_ix] : &coefs_[start_ix];
    fill(tempResult.begin(), tempResult.end(), double(0));

    for (double* bval_w_ptr = Bw.begin(); bval_w_ptr != Bw.end(); ++bval_w_ptr) {
      const double bval_w = *bval_w_ptr;
      fill(tempPt.begin(), tempPt.end(), 0);
      for (double* bval_v_ptr = Bv.begin(); bval_v_ptr != Bv.end(); ++bval_v_ptr) {
	const double bval_v = *bval_v_ptr;
	fill(tempPt2.begin(), tempPt2.end(), 0);
	for (double* bval_u_ptr = Bu.begin(); bval_u_ptr != Bu.end(); ++bval_u_ptr) {
	  const double bval_u = *bval_u_ptr;
	  for (ptemp = tempPt2.begin(); ptemp != tempPt2.end(); ++ptemp) {
	    *ptemp += bval_u * (*co_ptr++);
	  }
	}
	ptemp = tempPt2.begin();
	for (double* p = tempPt.begin(); p != tempPt.end(); ++p) {
	  *p += (*ptemp++) * bval_v;
	}
	co_ptr += kdim * (unum - uorder);
      }
      ptemp = tempPt.begin();
      for (double* p = tempResult.begin(); p != tempResult.end(); ++p) {
	*p += (*ptemp++) * bval_w;
      }
      co_ptr += kdim * unum * (vnum - vorder);
    }

    copy(tempResult.begin(), tempResult.begin() + dim_, pt.begin());
    if (rational_) {
	const double w_inv = double(1) / tempResult[kdim - 1];
	transform(pt.begin(), pt.end(), pt.begin(), ScaleBy(w_inv));
    }
}



//===========================================================================
void  SplineVolume::point(vector<Point>& pts, 
			  double upar, double vpar, double wpar,
			  int derivs,
			  bool u_from_right,
			  bool v_from_right,
			  bool w_from_right,
			  double resolution) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0, "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
    int rsz;
    rsz = (int)pts.size();
    DEBUG_ERROR_IF(rsz< totpts, "The vector of points must have sufficient size.");

    for (int i = 0; i < totpts; ++i) {
	if (pts[i].dimension() != dim_) {
	    pts[i].resize(dim_);
	}
    }

    if (derivs == 0) {
      point(pts[0], upar, vpar, wpar);
	return;
    }

    // Take care of the rational case
    const vector<double>& co = rational_ ? rcoefs_ : coefs_;
    int kdim = dim_ + (rational_ ? 1 : 0);

    // Make temporary storage for the basis values and a temporary
    // computation cache.
    Go::ScratchVect<double, 30> b0(basis_u_.order() * (derivs+1));
    Go::ScratchVect<double, 30> b1(basis_v_.order() * (derivs+1));
    Go::ScratchVect<double, 30> b2(basis_w_.order() * (derivs+1));
    Go::ScratchVect<double, 60> temp(kdim * totpts);
    Go::ScratchVect<double, 60> temp2(kdim * totpts);
    Go::ScratchVect<double, 60> restemp(kdim * totpts);
    fill(restemp.begin(), restemp.end(), 0.0);

    // Compute the basis values and get some data about the spline spaces
    if (u_from_right) {
	basis_u_.computeBasisValues(upar, &b0[0], derivs, resolution);
    } else {
	basis_u_.computeBasisValuesLeft(upar, &b0[0], derivs, resolution);
    }
    int uleft = basis_u_.lastKnotInterval();
    int uorder = basis_u_.order();
    int unum = basis_u_.numCoefs();
    if (v_from_right) {
	basis_v_.computeBasisValues(vpar, &b1[0], derivs, resolution);
    } else {
	basis_v_.computeBasisValuesLeft(vpar, &b1[0], derivs, resolution);
    }
    int vleft = basis_v_.lastKnotInterval();
    int vorder = basis_v_.order();
    int vnum = basis_v_.numCoefs();
    if (w_from_right) {
	basis_w_.computeBasisValues(wpar, &b2[0], derivs, resolution);
    } else {
	basis_w_.computeBasisValuesLeft(wpar, &b2[0], derivs, resolution);
    }
    int wleft = basis_w_.lastKnotInterval();
    int worder = basis_w_.order();

    // Compute the tensor product value
    int coefind = uleft-uorder+1 + unum*(vleft-vorder+1 + vnum*(wleft-worder+1));
    int derivs_plus1=derivs+1;

    for (int k = 0; k < worder; ++k) {
      int kd=k*derivs_plus1;
      fill(temp.begin(), temp.end(), 0.0);

      for (int j = 0; j < vorder; ++j) {
	int jd=j*derivs_plus1;
	fill(temp2.begin(), temp2.end(), 0.0);

	for (int i = 0; i < uorder; ++i) {
	  int id=i*derivs_plus1;
	  const double *co_p=&co[coefind*kdim];

	  for (int d = 0; d < kdim; ++d,++co_p) {
	    int temp_ind = d;

	    for (int wder = 0; wder <= derivs; ++wder) {
	      for (int vder = 0; vder <= wder; ++vder) {
		for (int uder = 0; uder <= vder; ++uder) {
		  temp2[temp_ind]
		    += b0[id + wder - vder]*(*co_p);
		  temp_ind+=kdim;
		}
	      }
	    }
	  }
	  coefind += 1;
	}

	for (int d = 0; d < kdim; ++d) {
	  int temp_ind = d;

	  for (int wder = 0; wder <= derivs; ++wder) {
	    for (int vder = 0; vder <= wder; ++vder) {
	      for (int uder = 0; uder <= vder; ++uder) {
		temp[temp_ind]
		  += temp2[temp_ind]*b1[jd + vder - uder];
		temp_ind+=kdim;
	      }
	    }
	  }
	}
	coefind += unum - uorder;
      }

      for (int d = 0; d < kdim; ++d) {
	int temp_ind = d;

	for (int wder = 0; wder <= derivs; ++wder) {
	  for (int vder = 0; vder <= wder; ++vder) {
	    for (int uder = 0; uder <= vder; ++uder) {
	      restemp[temp_ind]
		+= temp[temp_ind]*b2[kd + uder];
	      temp_ind+=kdim;
	    }
	  }
	}
      }
      coefind += unum * (vnum - vorder);
    }

    // Copy from restemp to result
    if (rational_) {
	vector<double> restemp2(totpts*dim_);
	volume_ratder(&restemp[0], dim_, derivs, &restemp2[0]);
	for (int i = 0; i < totpts; ++i) {
	    for (int d = 0; d < dim_; ++d) {
		pts[i][d] = restemp2[i*dim_ + d];
	    }
	}
    } else {
      double* restemp_it=restemp.begin();
	for (int i = 0; i < totpts; ++i) {
	    for (int d = 0; d < dim_; ++d) {
		pts[i][d] = *restemp_it;
		++restemp_it;
	    }
	}
    }
        
}



//===========================================================================
void  SplineVolume::computeBasis(double param[], 
				 vector< double > &basisValues,
				 vector< double > &basisDerivs_u,
				 vector< double > &basisDerivs_v,
				 vector< double > &basisDerivs_w,
				 bool evaluate_from_right) const
//===========================================================================
{
    int kk1 = basis_u_.order();
    int kk2 = basis_v_.order();
    int kk3 = basis_w_.order();
    int nn1 = basis_u_.numCoefs();
    int nn2 = basis_v_.numCoefs();
    vector<double> basisvals_u(2*kk1);
    vector<double> basisvals_v(2*kk2);
    vector<double> basisvals_w(2*kk3);

    // Compute basis values
    if (evaluate_from_right)
      {
	basis_u_.computeBasisValues(param[0], &basisvals_u[0], 1);
	basis_v_.computeBasisValues(param[1], &basisvals_v[0], 1);
	basis_w_.computeBasisValues(param[2], &basisvals_w[0], 1);
      }
    else 
      {
	basis_u_.computeBasisValuesLeft(param[0], &basisvals_u[0], 1);
	basis_v_.computeBasisValuesLeft(param[1], &basisvals_v[0], 1);
	basis_w_.computeBasisValuesLeft(param[2], &basisvals_w[0], 1);
       }

    // Accumulate
    int ki, kj, kh, kr;
    basisValues.resize(kk1*kk2*kk3);
    basisDerivs_u.resize(kk1*kk2*kk3);
    basisDerivs_v.resize(kk1*kk2*kk3);
    basisDerivs_w.resize(kk1*kk2*kk3);
    vector<double> weights;
    if (rational_)
    {
	int kdim = dim_ + 1;
	int uleft = basis_u_.lastKnotInterval() - kk1 + 1;
	int vleft = basis_v_.lastKnotInterval() - kk2 + 1;
	int wleft = basis_w_.lastKnotInterval() - kk3 + 1;
	weights.resize(kk1*kk2*kk3);
	for (kh=wleft, kr=0; kh<wleft+kk3; ++kh)
	    for (kj=vleft; kj<vleft+kk2; ++kj)
		for (ki=uleft; ki<uleft+kk1; ++ki)
		    weights[kr++] = rcoefs_[((kh*nn2+kj)*nn1+ki)*kdim+dim_];
    }

    accumulateBasis(&basisvals_u[0], kk1, &basisvals_v[0], kk2,
		    &basisvals_w[0], kk3, rational_ ? &weights[0] : NULL, 
		    &basisValues[0],
		    &basisDerivs_u[0], &basisDerivs_v[0],&basisDerivs_w[0]);
}



//===========================================================================
void SplineVolume::gridEvaluator(int num_u, int num_v, int num_w,
				 vector< double > &points,
				 vector< double > &der_u,
				 vector< double > &der_v,
				 vector< double > &der_w,
				 vector< double > &param_u,
				 vector< double > &param_v,
				 vector< double > &param_w,
				 bool evaluate_from_right) const
//===========================================================================
{
  const double start_u = startparam(0);
  const double start_v = startparam(1);
  const double start_w = startparam(2);
  const double du = (endparam(0) - start_u) / double(num_u-1);
  const double dv = (endparam(1) - start_v) / double(num_v-1);
  const double dw = (endparam(2) - start_w) / double(num_w-1);

  param_u.resize(num_u);
  param_v.resize(num_v);
  param_w.resize(num_w);

  for(int i = 0; i < num_u; ++i)
    param_u[i] = start_u+double(i)*du;
  for(int i = 0; i < num_v; ++i)
    param_v[i] = start_v+double(i)*dv;
  for(int i = 0; i < num_w; ++i)
    param_w[i] = start_w+double(i)*dw;

  gridEvaluator(param_u, param_v, param_w, points, der_u, der_v, der_w,
		evaluate_from_right);
}



//===========================================================================
void SplineVolume::gridEvaluator (const vector< double > &param_u,
				  const vector< double > &param_v,
				  const vector< double > &param_w,
				  vector< double > &points,
				  vector< double > &der_u,
				  vector< double > &der_v,
				  vector< double > &der_w,
				  bool evaluate_from_right) const
//===========================================================================
{
  int numb_pts = (int)(param_u.size() * param_v.size() * param_w.size());
  vector<double> pts_and_derivs(numb_pts * dim_ * 4);
  pointsGrid(param_u, param_v, param_w, 1, pts_and_derivs, evaluate_from_right);

  points.resize(numb_pts*dim_);
  der_u.resize(numb_pts*dim_);
  der_v.resize(numb_pts*dim_);
  der_w.resize(numb_pts*dim_);

  int pos_all = 0;
  int pos_result = 0;
  for (int i = 0; i < numb_pts; ++i, pos_result += dim_)
    {
      for (int j=0; j < dim_; ++j)
	points[pos_result+j] = pts_and_derivs[pos_all++];
      for (int j=0; j < dim_; ++j)
	der_u[pos_result+j] = pts_and_derivs[pos_all++];
      for (int j=0; j < dim_; ++j)
	der_v[pos_result+j] = pts_and_derivs[pos_all++];
      for (int j=0; j < dim_; ++j)
	der_w[pos_result+j] = pts_and_derivs[pos_all++];
    }
}



//===========================================================================
void SplineVolume::gridEvaluator(int num_u, int num_v, int num_w,
				 vector< double > &points,
				 vector< double > &param_u,
				 vector< double > &param_v,
				 vector< double > &param_w) const
//===========================================================================
{
  const double start_u = startparam(0);
  const double start_v = startparam(1);
  const double start_w = startparam(2);
  const double du = (endparam(0) - start_u) / double(num_u-1);
  const double dv = (endparam(1) - start_v) / double(num_v-1);
  const double dw = (endparam(2) - start_w) / double(num_w-1);

  param_u.resize(num_u);
  param_v.resize(num_v);
  param_w.resize(num_w);

  for(int i = 0; i < num_u; ++i)
    param_u[i] = start_u+double(i)*du;
  for(int i = 0; i < num_v; ++i)
    param_v[i] = start_v+double(i)*dv;
  for(int i = 0; i < num_w; ++i)
    param_w[i] = start_w+double(i)*dw;

  gridEvaluator(param_u, param_v, param_w, points);
}



//===========================================================================
void SplineVolume::gridEvaluator (const vector< double > &param_u,
				  const vector< double > &param_v,
				  const vector< double > &param_w,
				  vector< double > &points) const
//===========================================================================
{
  pointsGrid(param_u, param_v, param_w, 0, points);
}



//===========================================================================
void SplineVolume::pointsGrid(const vector< double > &param_u,
			      const vector< double > &param_v,
			      const vector< double > &param_w,
			      int derivs,
			      vector< double > &points,
			      bool evaluate_from_right) const
//===========================================================================
{
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<double> basisvals_w(numw * worder * (derivs + 1));
  vector<int>    knotinter_u(numu);
  vector<int>    knotinter_v(numv);
  vector<int>    knotinter_w(numw);

  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
				  &basisvals_u[0], &knotinter_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
				  &basisvals_v[0], &knotinter_v[0], derivs);
      basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
				  &basisvals_w[0], &knotinter_w[0], derivs);
    }
  else 
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0] + param_u.size(),
				      &basisvals_u[0], &knotinter_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0] + param_v.size(),
				      &basisvals_v[0], &knotinter_v[0], derivs);
      basis_w_.computeBasisValuesLeft(&param_w[0], &param_w[0] + param_w.size(),
				      &basisvals_w[0], &knotinter_w[0], derivs);
    }

  pointsGrid(numu, numv, numw, basisvals_u.begin(), basisvals_v.begin(), basisvals_w.begin(),
	     knotinter_u.begin(), knotinter_v.begin(), knotinter_w.begin(), derivs, points);
}


//===========================================================================
void SplineVolume::pointsGrid(int numu, int numv, int numw,
			      vector<double>::iterator basisvals_u,
			      vector<double>::iterator basisvals_v,
			      vector<double>::iterator basisvals_w,
			      vector<int>::iterator knotinter_u,
			      vector<int>::iterator knotinter_v,
			      vector<int>::iterator knotinter_w,
			      int derivs,
			      vector< double > &points) const
//===========================================================================
{
  int kdim;
  const double* scoef;

  if (rational_) {
    scoef = &rcoefs_[0];
    kdim = dim_ + 1;
  } else {
    scoef = &coefs_[0];
    kdim = dim_;
  }

  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();

  int size_dwjip = ucoefs * vcoefs * (derivs+1) * kdim;
  vector<double> temp_dwjip(size_dwjip);
  int size_dvdwip = (ucoefs * (derivs+1)*(derivs+2) * kdim) >> 1;
  vector<double> temp_dvdwip(size_dvdwip);
  int size_dudvdwp = ((derivs+1) * (derivs+2) * (derivs+3) * kdim) / 6;
  vector<double> temp_dudvdwp(size_dudvdwp);

  int points_pos = 0;
  points.resize(numu*numv*numw*dimension()*(derivs+1)*(derivs+2)*(derivs+3)/6);

  // Loop through all parameter values in third direction
  for(int idx_w = 0, basisw_pos = 0;  idx_w < numw; ++idx_w, basisw_pos += (derivs+1)*worder) {

    int basis_left = knotinter_w[idx_w];

    /* Compute the control points and derivatives of the
       w = param_w[idx_w] isosurface. Store in temp_dwjip */
    fill(temp_dwjip.begin(), temp_dwjip.end(), 0.0);

    int local_basis_pos = basisw_pos;
    for (int k = basis_left - worder + 1; k <= basis_left; ++k)
      {
	int dwjip_pos = 0;
	for (int dw = 0; dw <= derivs; ++dw)
	  {
	    double basisval = basisvals_w[local_basis_pos++];
	    int scoef_pos = kdim * ucoefs * vcoefs * k;
	    for (int k2 = 0; k2 < ucoefs*vcoefs*kdim; ++k2)
	      temp_dwjip[dwjip_pos++] += basisval * scoef[scoef_pos++];
	  }
      }

    // Loop through all parameter values in second direction
    for(int idx_v = 0, basisv_pos = 0;  idx_v < numv; ++idx_v, basisv_pos += (derivs+1)*vorder) {

      basis_left = knotinter_v[idx_v];

      /* Compute the control points and derivatives of the
	 v = param_v[idx_v], w = param_w[idx_w] isocurve.
	 Store in temp_dvdwip */
      fill(temp_dvdwip.begin(), temp_dvdwip.end(), 0.0);

      local_basis_pos = basisv_pos;
      for (int j = basis_left - vorder + 1; j <= basis_left; ++j)
	for (int dv = 0; dv <= derivs; ++dv)
	  {
	    double basisval = basisvals_v[local_basis_pos++];
	    for (int dw = 0; dw <= derivs-dv; ++dw)
	      {
		int dtot = dv+dw;
		int dvdwip_pos = ucoefs * kdim * (((dtot * (dtot+1)) >> 1) + dw);
		int dwjip_pos = ucoefs * kdim * (dw * vcoefs + j);
		for (int j2 = 0; j2 < ucoefs * kdim; ++j2)
		  temp_dvdwip[dvdwip_pos++] += basisval * temp_dwjip[dwjip_pos++];
	      }
	  }

      // Loop through all parameter values in first direction
      for(int idx_u = 0, basisu_pos = 0;  idx_u < numu; ++idx_u, basisu_pos += (derivs+1)*uorder) {

	basis_left = knotinter_u[idx_u];

	/* Compute the control points and derivatives of the point.
	   Store in temp_dudvdwp */
	fill(temp_dudvdwp.begin(), temp_dudvdwp.end(), 0.0);

	local_basis_pos = basisu_pos;
	for (int i = basis_left - uorder + 1; i <= basis_left; ++i)
	  for (int du = 0; du <= derivs; ++du)
	    {
	      double basisval = basisvals_u[local_basis_pos++];
	      for (int dv = 0; dv <= derivs-du; ++dv)
		for (int dw = 0; dw <= derivs-du-dv; ++dw)
		  {
		    int dvw = dv+dw;
		    int dvw_pos = ((dvw*(dvw+1)) >> 1) + dw;
		    int dtot = du+dvw;
		    int dudvdwp_pos = kdim * ((dtot*(dtot+1)*(dtot+2)) / 6 + dvw_pos);
		    int dvdwip_pos = kdim * (ucoefs * dvw_pos + i);
		    for (int i2 = 0; i2 < kdim; ++i2)
		      temp_dudvdwp[dudvdwp_pos++] += basisval * temp_dvdwip[dvdwip_pos++];
		  }
	    }

	// Store result in point vector. Handle rational case
	if (rational_)
	  {
	    volume_ratder(&temp_dudvdwp[0], dim_, derivs, &points[points_pos]);
	    points_pos += dim_*(derivs+1)*(derivs+2)*(derivs+3)/6;
	  }
	else
	  for (int i = 0; i < dim_*(derivs+1)*(derivs+2)*(derivs+3)/6; ++i)
	    points[points_pos++] = temp_dudvdwp[i];

      }
    }
  }
}


//===========================================================================
void SplineVolume::computeBasisGrid(const Dvector& param_u,
				    const Dvector& param_v,
				    const Dvector& param_w,
				    Dmatrix& basisValues) const 
//===========================================================================
{
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();
  int wcoefs = basis_w_.numCoefs();

  vector<double> basisvals_u(numu * uorder);
  vector<double> basisvals_v(numv * vorder);
  vector<double> basisvals_w(numw * worder);
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);
  vector<int>    left_w(numw);

  // Compute basis values
  basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
			      &basisvals_u[0], &left_u[0]);
  basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
			      &basisvals_v[0], &left_v[0]);
  basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
			      &basisvals_w[0], &left_w[0]);

  // Initiate to zero
  int ki;
  int num_coefs = ucoefs*vcoefs*wcoefs;
  int num_par = numu*numv*numw;
  basisValues.resize(num_par);
  for (ki=0; ki<num_par; ++ki)
    basisValues[ki].assign(num_coefs, 0.0);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder*worder);
      weights.resize(ucoefs*vcoefs*wcoefs);
      getWeights(weights);
  }

  // For all points
  int kj, kr, kh, kv, kw;
  int idx1, idx2, idx3, idx4, idx5;
  int numorder = uorder*vorder*worder;
  vector<double> tmpVal(numorder);
  for (kr=0, idx3=0, kh=0; kr<numw; ++kr, idx3+=worder)
  {
      for (kj=0, idx2=0; kj<numv; ++kj, idx2+=vorder)
      {
	  for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=uorder)
	  {
	      int uleft = left_u[ki] - uorder + 1;
	      int vleft = left_v[kj] - vorder + 1;
	      int wleft = left_w[kr] - worder + 1;
	      if (rational_)
	      {
		  // Collect relevant weights
		  vector<double>::iterator wgt = weights.begin() + (wleft*vcoefs + vleft)*ucoefs;
		  vector<double>::iterator currwgt = currw.begin();
		  for (kw=0; kw<worder; ++kw, wgt += (vcoefs-vorder)*ucoefs)
		  {
		      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
		      {
			  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
		      }
		  }
	      }
      
	      accumulateBasis(&basisvals_u[idx1], uorder, &basisvals_v[idx2],
			      vorder, &basisvals_w[idx3], worder, 
			      rational_ ? &currw[0] : NULL, 
			      &tmpVal[0]);

	      // Copy results into output array
	      for (kw=0, idx4=0, idx5=(wleft*vcoefs+vleft)*ucoefs; kw<worder; 
		   ++kw, idx5+=(vcoefs-vorder)*ucoefs-uleft)
		  for (kv=0, idx5+=uleft; kv<vorder; ++kv, idx4+=uorder, idx5+=ucoefs)
		  {
		      std::copy(tmpVal.begin()+idx4, tmpVal.begin()+idx4+uorder, 
				basisValues[kh].begin()+idx5);
		  }
	  }
      }
  }
  
}

//===========================================================================
void SplineVolume::computeBasisGrid(const Dvector& param_u,
				    const Dvector& param_v,
				    const Dvector& param_w,
				    Dmatrix& basisValues,
				    Dmatrix& basisDerivs_u,
				    Dmatrix& basisDerivs_v,
				    Dmatrix& basisDerivs_w,
				    bool evaluate_from_right) const 
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();
  int wcoefs = basis_w_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<double> basisvals_w(numw * worder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);
  vector<int>    left_w(numw);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
				  &basisvals_w[0], &left_w[0], derivs);
    }
  else 
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0] + param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0] + param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValuesLeft(&param_w[0], &param_w[0] + param_w.size(),
				      &basisvals_w[0], &left_w[0], derivs);
    }

  // Initiate to zero
  int ki;
  int num_coefs = ucoefs*vcoefs*wcoefs;
  int num_par = numu*numv*numw;
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

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder*worder);
      weights.resize(ucoefs*vcoefs*wcoefs);
      getWeights(weights);
  }

  // For all points
  int kj, kr, kh, kv, kw;
  int idx1, idx2, idx3, idx4, idx5;
  int numorder = uorder*vorder*worder;
  vector<double> tmpVal(numorder), tmpDer_u(numorder), tmpDer_v(numorder), tmpDer_w(numorder);
  for (kr=0, idx3=0, kh=0; kr<numw; ++kr, idx3+=2*worder)
  {
      for (kj=0, idx2=0; kj<numv; ++kj, idx2+=2*vorder)
      {
	  for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=2*uorder)
	  {
	      int uleft = left_u[ki] - uorder + 1;
	      int vleft = left_v[kj] - vorder + 1;
	      int wleft = left_w[kr] - worder + 1;
	      if (rational_)
	      {
		  // Collect relevant weights
		  vector<double>::iterator wgt = weights.begin() + (wleft*vcoefs + vleft)*ucoefs;
		  vector<double>::iterator currwgt = currw.begin();
		  for (kw=0; kw<worder; ++kw, wgt += (vcoefs-vorder)*ucoefs)
		  {
		      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
		      {
			  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
		      }
		  }
	      }
      
	      accumulateBasis(&basisvals_u[idx1], uorder, &basisvals_v[idx2],
			      vorder, &basisvals_w[idx3], worder, 
			      rational_ ? &currw[0] : NULL, 
			      &tmpVal[0], &tmpDer_u[0], &tmpDer_v[0], &tmpDer_w[0]);
	      // Copy results into output array
	      for (kw=0, idx4=0, idx5=(wleft*vcoefs+vleft)*ucoefs; kw<worder; 
		   ++kw, idx5+=(vcoefs-vorder)*ucoefs-uleft)
		  for (kv=0, idx5+=uleft; kv<vorder; ++kv, idx4+=uorder, idx5+=ucoefs)
		  {
		      std::copy(tmpVal.begin()+idx4, tmpVal.begin()+idx4+uorder, 
				basisValues[kh].begin()+idx5);
		      std::copy(tmpDer_u.begin()+idx4, tmpDer_u.begin()+idx4+uorder, 
				basisDerivs_u[kh].begin()+idx5);
		      std::copy(tmpDer_v.begin()+idx4, tmpDer_v.begin()+idx4+uorder, 
				basisDerivs_v[kh].begin()+idx5);
		      std::copy(tmpDer_w.begin()+idx4, tmpDer_w.begin()+idx4+uorder, 
				basisDerivs_w[kh].begin()+idx5);
		  }
	  }
      }
  }
  
}

//===========================================================================
void SplineVolume::computeBasis(double param_u,
				double param_v,
				double param_w,
				BasisPts& result) const
//===========================================================================
{
    int uorder = basis_u_.order();
    int vorder = basis_v_.order();
    int worder = basis_w_.order();
    int nn1 = basis_u_.numCoefs();
    int nn2 = basis_v_.numCoefs();
    vector<double> basisvals_u(uorder);
    vector<double> basisvals_v(vorder);
    vector<double> basisvals_w(worder);

    // Compute basis values
    basis_u_.computeBasisValues(param_u, &basisvals_u[0], 0);
    basis_v_.computeBasisValues(param_v, &basisvals_v[0], 0);
    basis_w_.computeBasisValues(param_w, &basisvals_w[0], 0);

    int ulast = basis_u_.lastKnotInterval();
    int vlast = basis_v_.lastKnotInterval();
    int wlast = basis_w_.lastKnotInterval();
    result.preparePts(param_u, param_v, param_w,
		      ulast, vlast, wlast,
		      uorder*vorder*worder);

    vector<double> weights;
   if (rational_)
    {
      // Collect relevant weights
      int kh, kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      int wleft = wlast - worder + 1;
      weights.resize(uorder*vorder*worder);
      for (kh=wleft, kr=0; kh<wleft+worder; ++kh)
	for (kj=vleft; kj<vleft+vorder; ++kj)
	  for (ki=uleft; ki<uleft+uorder; ++ki)
	    weights[kr++] = rcoefs_[((kh*nn2+kj)*nn1+ki)*kdim+dim_];
    }
      
  accumulateBasis(&basisvals_u[0], uorder, &basisvals_v[0],
		  vorder, &basisvals_w[0], worder, 
		  rational_ ? &weights[0] : NULL, 
		  &result.basisValues[0]);
}

//===========================================================================
void SplineVolume::computeBasis(double param_u,
				double param_v,
				double param_w,
				BasisDerivs& result,
				bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int nn1 = basis_u_.numCoefs();
  int nn2 = basis_v_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));
  vector<double> basisvals_w(worder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
      basis_w_.computeBasisValues(param_w, &basisvals_w[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
      basis_w_.computeBasisValuesLeft(param_w, &basisvals_w[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  int wlast = basis_w_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, param_w,
		       ulast, vlast, wlast,
		       uorder*vorder*worder);

  vector<double> weights;
  if (rational_)
    {
      // Collect relevant weights
      int kh, kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      int wleft = wlast - worder + 1;
      weights.resize(uorder*vorder*worder);
      for (kh=wleft, kr=0; kh<wleft+worder; ++kh)
	for (kj=vleft; kj<vleft+vorder; ++kj)
	  for (ki=uleft; ki<uleft+uorder; ++ki)
	      weights[kr++] = rcoefs_[((kh*nn2+kj)*nn1+ki)*kdim+dim_];
    }
 
  accumulateBasis(&basisvals_u[0], uorder, &basisvals_v[0],
		  vorder, &basisvals_w[0], worder, 
		  rational_ ? &weights[0] : NULL, 
		  &result.basisValues[0], &result.basisDerivs_u[0], 
		  &result.basisDerivs_v[0],&result.basisDerivs_w[0]);
}


//===========================================================================
void SplineVolume::computeBasis(double param_u,
				double param_v,
				double param_w,
				BasisDerivs2& result,
				bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 2;  // Compute position, 1. and 2. derivative
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int nn1 = basis_u_.numCoefs();
  int nn2 = basis_v_.numCoefs();
  vector<double> basisvals_u(uorder * (derivs + 1));
  vector<double> basisvals_v(vorder * (derivs + 1));
  vector<double> basisvals_w(worder * (derivs + 1));

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValues(param_v, &basisvals_v[0], derivs);
      basis_w_.computeBasisValues(param_w, &basisvals_w[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(param_u, &basisvals_u[0], derivs);
      basis_v_.computeBasisValuesLeft(param_v, &basisvals_v[0], derivs);
      basis_w_.computeBasisValuesLeft(param_w, &basisvals_w[0], derivs);
    }

  int ulast = basis_u_.lastKnotInterval();
  int vlast = basis_v_.lastKnotInterval();
  int wlast = basis_w_.lastKnotInterval();
  result.prepareDerivs(param_u, param_v, param_w,
		       ulast, vlast, wlast,
		       uorder*vorder*worder);

  vector<double> weights;
  if (rational_)
    {
      // Collect relevant weights
      int kh, kr, ki, kj;
      int kdim = dim_ + 1;
      int uleft = ulast - uorder + 1;
      int vleft = vlast - vorder + 1;
      int wleft = wlast - worder + 1;
      weights.resize(uorder*vorder*worder);
      for (kh=wleft, kr=0; kh<wleft+worder; ++kh)
	for (kj=vleft; kj<vleft+vorder; ++kj)
	  for (ki=uleft; ki<uleft+uorder; ++ki)
	      weights[kr++] = rcoefs_[((kh*nn2+kj)*nn1+ki)*kdim+dim_];
    }
 
  accumulateBasis(&basisvals_u[0], uorder, &basisvals_v[0],
		  vorder, &basisvals_w[0], worder, 
		  rational_ ? &weights[0] : NULL, 
		  &result.basisValues[0], &result.basisDerivs_u[0], 
		  &result.basisDerivs_v[0], &result.basisDerivs_w[0],
		  &result.basisDerivs_uu[0], &result.basisDerivs_uv[0],&result.basisDerivs_uw[0],
		  &result.basisDerivs_vv[0], &result.basisDerivs_vw[0],&result.basisDerivs_ww[0]);
}


//===========================================================================
void SplineVolume::computeBasis(const vector<double>::const_iterator& bas_vals_u,
				const vector<double>::const_iterator& bas_vals_v,
				const vector<double>::const_iterator& bas_vals_w,
				int left_u,
				int left_v,
				int left_w,
				vector<double>& basisValues,
				vector<double>& basisDerivs_u,
				vector<double>& basisDerivs_v,
				vector<double>& basisDerivs_w) const
//===========================================================================
{
    int kk1 = basis_u_.order();
    int kk2 = basis_v_.order();
    int kk3 = basis_w_.order();
    int nn1 = basis_u_.numCoefs();
    int nn2 = basis_v_.numCoefs();

    // Accumulate
    int ki, kj, kk, kr;
    basisValues.resize(kk1*kk2*kk3);
    basisDerivs_u.resize(kk1*kk2*kk3);
    basisDerivs_v.resize(kk1*kk2*kk3);
    basisDerivs_w.resize(kk1*kk2*kk3);
    vector<double> weights(kk1*kk2*kk3);
    if (rational_)
    {
	int kdim = dim_ + 1;
	int uleft = left_u - kk1 + 1;
	int vleft = left_v - kk2 + 1;
	int wleft = left_w - kk3 + 1;
	for (kk=wleft, kr=0; kk<wleft+kk3; ++kk)
	    for (kj=vleft; kj<vleft+kk2; ++kj)
		for (ki=uleft; ki<uleft+kk1; ++ki)
		    weights[kr++] = rcoefs_[(kk*nn1*nn2+kj*nn1+ki)*kdim+dim_];
    }

    accumulateBasis(&bas_vals_u[0], kk1,
		    &bas_vals_v[0], kk2,
		    &bas_vals_w[0], kk3,
		    &weights[0],
		    &basisValues[0],
		    &basisDerivs_u[0],
		    &basisDerivs_v[0],
		    &basisDerivs_w[0]);

}


//===========================================================================
void SplineVolume::computeBasisGrid(const Dvector& param_u,
				    const Dvector& param_v,
				    const Dvector& param_w,
				    vector<BasisPts>& result) const
//===========================================================================
{
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();
  int wcoefs = basis_w_.numCoefs();

  vector<double> basisvals_u(numu * uorder);
  vector<double> basisvals_v(numv * vorder);
  vector<double> basisvals_w(numw * worder);
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);
  vector<int>    left_w(numw);

  // Compute basis values
  basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
			      &basisvals_u[0], &left_u[0]);
  basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
			      &basisvals_v[0], &left_v[0]);
  basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
			      &basisvals_w[0], &left_w[0]);

  // Initiate output
  result.resize(numu*numv*numw);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder*worder);
      weights.resize(ucoefs*vcoefs*wcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kr, kh, kv, kw;
  int idx1, idx2, idx3;
  for (kr=0, idx3=0, kh=0; kr<numw; ++kr, idx3+=worder)
  {
      for (kj=0, idx2=0; kj<numv; ++kj, idx2+=vorder)
      {
	  for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=uorder)
	  {
	      result[kh].preparePts(param_u[ki], param_v[kj], param_w[kr],
				    left_u[ki], left_v[kj], left_w[kr], 
				    uorder*vorder*worder);


	      if (rational_)
	      {
		  // Collect relevant weights
		  int uleft = left_u[ki] - uorder + 1;
		  int vleft = left_v[kj] - vorder + 1;
		  int wleft = left_w[kr] - worder + 1;
		  vector<double>::iterator wgt = weights.begin() + (wleft*vcoefs + vleft)*ucoefs;
		  vector<double>::iterator currwgt = currw.begin();
		  for (kw=0; kw<worder; ++kw, wgt += (vcoefs-vorder)*ucoefs)
		  {
		      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
		      {
			  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
		      }
		  }
	      }
      
	      accumulateBasis(&basisvals_u[idx1], uorder, &basisvals_v[idx2],
			      vorder, &basisvals_w[idx3], worder, 
			      rational_ ? &currw[0] : NULL, 
			      &result[kh].basisValues[0]);
	  }
      }
  }
  
}

//===========================================================================
void SplineVolume::computeBasisGrid(const Dvector& param_u,
				    const Dvector& param_v,
				    const Dvector& param_w,
				    vector<BasisDerivs>& result,
				    bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 1;  // Compute position  and 1. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();
  int wcoefs = basis_w_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<double> basisvals_w(numw * worder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);
  vector<int>    left_w(numw);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
				  &basisvals_w[0], &left_w[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0] + param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0] + param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValuesLeft(&param_w[0], &param_w[0] + param_w.size(),
				      &basisvals_w[0], &left_w[0], derivs);
    }

  // Initiate output
  result.resize(numu*numv*numw);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder*worder);
      weights.resize(ucoefs*vcoefs*wcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kr, kh, kv, kw;
  int idx1, idx2, idx3;
  for (kr=0, idx3=0, kh=0; kr<numw; ++kr, idx3+=2*worder)
  {
      for (kj=0, idx2=0; kj<numv; ++kj, idx2+=2*vorder)
      {
	  for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=2*uorder)
	  {
	      result[kh].prepareDerivs(param_u[ki], param_v[kj], param_w[kr],
				      left_u[ki], left_v[kj], left_w[kr], 
				      uorder*vorder*worder);


	      if (rational_)
	      {
		  // Collect relevant weights
		  int uleft = left_u[ki] - uorder + 1;
		  int vleft = left_v[kj] - vorder + 1;
		  int wleft = left_w[kr] - worder + 1;
		  vector<double>::iterator wgt = weights.begin() + (wleft*vcoefs + vleft)*ucoefs;
		  vector<double>::iterator currwgt = currw.begin();
		  for (kw=0; kw<worder; ++kw, wgt += (vcoefs-vorder)*ucoefs)
		  {
		      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
		      {
			  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
		      }
		  }
	      }
      
	      accumulateBasis(&basisvals_u[idx1], uorder, &basisvals_v[idx2],
			      vorder, &basisvals_w[idx3], worder, 
			      rational_ ? &currw[0] : NULL, 
			      &result[kh].basisValues[0], &result[kh].basisDerivs_u[0], 
			      &result[kh].basisDerivs_v[0],&result[kh].basisDerivs_w[0]);
	  }
      }
  }
  
}



//===========================================================================
void SplineVolume::computeBasisGrid(const Dvector& param_u,
				    const Dvector& param_v,
				    const Dvector& param_w,
				    vector<BasisDerivs2>& result,
				    bool evaluate_from_right) const
//===========================================================================
{
  int derivs = 2;  // Compute position, 1. and 2. derivative
  int numu = (int)param_u.size();
  int numv = (int)param_v.size();
  int numw = (int)param_w.size();
  int uorder = basis_u_.order();
  int vorder = basis_v_.order();
  int worder = basis_w_.order();
  int ucoefs = basis_u_.numCoefs();
  int vcoefs = basis_v_.numCoefs();
  int wcoefs = basis_w_.numCoefs();

  vector<double> basisvals_u(numu * uorder * (derivs + 1));
  vector<double> basisvals_v(numv * vorder * (derivs + 1));
  vector<double> basisvals_w(numw * worder * (derivs + 1));
  vector<int>    left_u(numu);
  vector<int>    left_v(numv);
  vector<int>    left_w(numw);

  // Compute basis values
  if (evaluate_from_right)
    {
      basis_u_.computeBasisValues(&param_u[0], &param_u[0] + param_u.size(),
				  &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValues(&param_v[0], &param_v[0] + param_v.size(),
				  &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValues(&param_w[0], &param_w[0] + param_w.size(),
				  &basisvals_w[0], &left_w[0], derivs);
    }
  else
    {
      basis_u_.computeBasisValuesLeft(&param_u[0], &param_u[0] + param_u.size(),
				      &basisvals_u[0], &left_u[0], derivs);
      basis_v_.computeBasisValuesLeft(&param_v[0], &param_v[0] + param_v.size(),
				      &basisvals_v[0], &left_v[0], derivs);
      basis_w_.computeBasisValuesLeft(&param_w[0], &param_w[0] + param_w.size(),
				      &basisvals_w[0], &left_w[0], derivs);
    }

  // Initiate output
  result.resize(numu*numv*numw);

  // Fetch all weights
  vector<double> weights, currw;
  if (rational_)
  {
      currw.resize(uorder*vorder*worder);
      weights.resize(ucoefs*vcoefs*wcoefs);
      getWeights(weights);
  }

  // For all points
  int ki, kj, kr, kh, kv, kw;
  int idx1, idx2, idx3;
  for (kr=0, idx3=0, kh=0; kr<numw; ++kr, idx3+=3*worder)
  {
      for (kj=0, idx2=0; kj<numv; ++kj, idx2+=3*vorder)
      {
	  for (ki=0, idx1=0; ki<numu; ++ki, ++kh, idx1+=3*uorder)
	  {
	      result[kh].prepareDerivs(param_u[ki], param_v[kj], param_w[kr],
				      left_u[ki], left_v[kj], left_w[kr], 
				      uorder*vorder*worder);

	      if (rational_)
	      {
		  // Collect relevant weights
		  int uleft = left_u[ki] - uorder + 1;
		  int vleft = left_v[kj] - vorder + 1;
		  int wleft = left_w[kr] - worder + 1;
		  vector<double>::iterator wgt = weights.begin() + (wleft*vcoefs + vleft)*ucoefs;
		  vector<double>::iterator currwgt = currw.begin();
		  for (kw=0; kw<worder; ++kw, wgt += (vcoefs-vorder)*ucoefs)
		  {
		      for (kv=0; kv<vorder; ++kv, wgt+=ucoefs, currwgt+=uorder)
		      {
			  std::copy(wgt+uleft, wgt+uleft+uorder, currwgt);
		      }
		  }
	      }
      
	      accumulateBasis(&basisvals_u[idx1], uorder, &basisvals_v[idx2],
			      vorder, &basisvals_w[idx3], worder, 
			      rational_ ? &currw[0] : NULL, 
			      &result[kh].basisValues[0],
			      &result[kh].basisDerivs_u[0], &result[kh].basisDerivs_v[0],&result[kh].basisDerivs_w[0],
			      &result[kh].basisDerivs_uu[0], &result[kh].basisDerivs_uv[0],&result[kh].basisDerivs_uw[0],
			      &result[kh].basisDerivs_vv[0], &result[kh].basisDerivs_vw[0],&result[kh].basisDerivs_ww[0]);
	  }
      }
  }
  
}



//===========================================================================
void SplineVolume::accumulateBasis(double* basisvals_u, int uorder,
				   double* basisvals_v, int vorder,
				   double* basisvals_w, int worder,
				   double* weights, 
				   double* basisValues) const
//===========================================================================
{
    int ki, kj, kh, kr;
    if (rational_)
    {
	// Compute denominator
	double sum = 0.0;
	for (kh=0, kr=0; kh<worder; ++kh)
	  for (kj=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	      sum += basisvals_u[ki]*basisvals_v[kj]*basisvals_w[kh]*weights[kr];

	// Compute rational expression
//	double sum2 = sum*sum;
	for (kh=0, kr=0; kh<worder; ++kh)
	  for (kj=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	      basisValues[kr] = basisvals_u[ki]*basisvals_v[kj]*basisvals_w[kh]*weights[kr]/sum;
    }
    else
    {
	// Multiply basis values in the three parameter directions
        for (kh=0, kr=0; kh<worder; ++kh)
	  for (kj=0; kj<vorder; ++kj)
	    for (ki=0; ki<uorder; ++ki, ++kr)
	      basisValues[kr] = basisvals_u[ki]*basisvals_v[kj]*basisvals_w[kh];
    }

}


//===========================================================================
void SplineVolume::accumulateBasis(const double* basisvals_u, int uorder,
				   const double* basisvals_v, int vorder,
				   const double* basisvals_w, int worder,
				   const double* weights, 
				   double* basisValues,
				   double* basisDerivs_u,
				   double* basisDerivs_v,
				   double* basisDerivs_w) const
//===========================================================================
{
    int ki, kj, kh, kr;
    if (rational_)
    {
	// Compute denominator and derivatives thereof
	double sum = 0.0, dusum = 0.0, dvsum = 0.0, dwsum = 0.0;
	for (kh=0, kr=0; kh<worder; ++kh)
	    for (kj=0; kj<vorder; ++kj)
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    sum += basisvals_u[2*ki]*basisvals_v[2*kj]*basisvals_w[2*kh]*weights[kr];
		    dusum += basisvals_u[2*ki+1]*basisvals_v[2*kj]*basisvals_w[2*kh]*weights[kr];
		    dvsum += basisvals_u[2*ki]*basisvals_v[2*kj+1]*basisvals_w[2*kh]*weights[kr];
		    dwsum += basisvals_u[2*ki]*basisvals_v[2*kj]*basisvals_w[2*kh+1]*weights[kr];
		}

	// Compute rational expression
	double sum2 = sum*sum;
	for (kh=0, kr=0; kh<worder; ++kh)
	{
	    double tmp3 = (basisvals_w[2*kh+1]*sum - basisvals_w[2*kh]*dwsum)/sum2;
	    double fac = basisvals_w[2*kh]/sum2;
	    for (kj=0; kj<vorder; ++kj)
	    {
		double tmp2 = (basisvals_v[2*kj+1]*sum - basisvals_v[2*kj]*dvsum)*fac;
		double tmp = tmp3*basisvals_v[2*kj];
		double fac2 = fac*basisvals_v[2*kj];
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    basisValues[kr] = basisvals_u[2*ki]*basisvals_v[2*kj]*basisvals_w[2*kh]*weights[kr]/sum;
		    double tmp1 = (basisvals_u[2*ki+1]*sum - basisvals_u[2*ki]*dusum)*weights[kr]*fac2;
		    basisDerivs_u[kr] = tmp1;
		    basisDerivs_v[kr] = tmp2*weights[kr]*basisvals_u[2*ki];
		    basisDerivs_w[kr] = tmp*weights[kr]*basisvals_u[2*ki];
		}
	    }
	}
    }
    else
    {
	// Multiply basis values in the three parameter directions
	for (kh=0, kr=0; kh<worder; ++kh)
	    for (kj=0; kj<vorder; ++kj)
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    basisValues[kr] = basisvals_u[2*ki]*basisvals_v[2*kj]*basisvals_w[2*kh];
		    basisDerivs_u[kr] = basisvals_u[2*ki+1]*basisvals_v[2*kj]*basisvals_w[2*kh];
		    basisDerivs_v[kr] = basisvals_u[2*ki]*basisvals_v[2*kj+1]*basisvals_w[2*kh];
		    basisDerivs_w[kr] = basisvals_u[2*ki]*basisvals_v[2*kj]*basisvals_w[2*kh+1];
		}
    }

}


//===========================================================================
void SplineVolume::accumulateBasis(double* basisvals_u, int uorder,
				   double* basisvals_v, int vorder,
				   double* basisvals_w, int worder,
				   double* weights, 
				   double* basisValues,
				   double* basisDerivs_u,
				   double* basisDerivs_v,
				   double* basisDerivs_w,
				   double* basisDerivs_uu,
				   double* basisDerivs_uv,
				   double* basisDerivs_uw,
				   double* basisDerivs_vv,
				   double* basisDerivs_vw,
				   double* basisDerivs_ww) const
//===========================================================================
{
    int ki, kj, kh, kr;
    if (rational_)
    {
	// Compute denominator and derivatives thereof
	double
	  sum = 0.0,
	  dusum = 0.0, dvsum = 0.0, dwsum = 0.0,
	  duusum = 0.0, duvsum = 0.0, duwsum = 0.0,
	  dvvsum = 0.0, dvwsum = 0.0, dwwsum = 0.0;
	for (kh=0, kr=0; kh<worder; ++kh)
	    for (kj=0; kj<vorder; ++kj)
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    sum += basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh]*weights[kr];
		    dusum += basisvals_u[3*ki+1]*basisvals_v[3*kj]*basisvals_w[3*kh]*weights[kr];
		    dvsum += basisvals_u[3*ki]*basisvals_v[3*kj+1]*basisvals_w[3*kh]*weights[kr];
		    dwsum += basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh+1]*weights[kr];
		    duusum += basisvals_u[3*ki+2]*basisvals_v[3*kj]*basisvals_w[3*kh]*weights[kr];
		    duvsum += basisvals_u[3*ki+1]*basisvals_v[3*kj+1]*basisvals_w[3*kh]*weights[kr];
		    duwsum += basisvals_u[3*ki+1]*basisvals_v[3*kj]*basisvals_w[3*kh+1]*weights[kr];
		    dvvsum += basisvals_u[3*ki]*basisvals_v[3*kj+2]*basisvals_w[3*kh]*weights[kr];
		    dvwsum += basisvals_u[3*ki]*basisvals_v[3*kj+1]*basisvals_w[3*kh+1]*weights[kr];
		    dwwsum += basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh+2]*weights[kr];
		}

	// Compute rational expression
	double sum2 = sum*sum;
	double sum3 = sum*sum2;
	for (kh=0, kr=0; kh<worder; ++kh)
	{
	    double tmp_w_w = (basisvals_w[3*kh+1]*sum - basisvals_w[3*kh]*dwsum)/sum2;
	    double tmp_w_ww = (basisvals_w[3*kh+2]*sum2 + 2.0*basisvals_w[3*kh]*dwsum*dwsum
			       - basisvals_w[3*kh]*dwwsum*sum - 2.0*basisvals_w[3*kh+1]*dwsum*sum) / sum3;
	    double fac_w2 = basisvals_w[3*kh] / sum2;
	    double fac_w3 = basisvals_w[3*kh] / sum3;
	    for (kj=0; kj<vorder; ++kj)
	    {
		double tmp_vw_v = (basisvals_v[3*kj+1]*sum - basisvals_v[3*kj]*dvsum)*fac_w2;
		double tmp_vw_vv = (basisvals_v[3*kj+2]*sum2 + 2.0*basisvals_v[3*kj]*dvsum*dvsum
				    - basisvals_v[3*kj]*dvvsum*sum - 2.0*basisvals_v[3*kj+1]*dvsum*sum) * fac_w3;
		double tmp_vw_vw = (basisvals_v[3*kj+1]*basisvals_w[3*kh+1]*sum2 + 2.0*basisvals_v[3*kj]*basisvals_w[3*kh]*dvsum*dwsum
				    - basisvals_v[3*kj+1]*basisvals_w[3*kh]*dwsum*sum - basisvals_v[3*kj]*basisvals_w[3*kh+1]*dvsum*sum
				    - basisvals_v[3*kj]*basisvals_w[3*kh]*dvwsum*sum) / sum3;
		double tmp_vw_w = tmp_w_w * basisvals_v[3*kj];
		double tmp_vw_ww = tmp_w_ww * basisvals_v[3*kj];
		double fac_vw2 = fac_w2 * basisvals_v[3*kj];
		double fac_vw3 = fac_w3 * basisvals_v[3*kj];
		double fac_v3 = basisvals_v[3*kj] / sum3;
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    basisValues[kr] = basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh]*weights[kr]/sum;

		    basisDerivs_u[kr] = (basisvals_u[3*ki+1]*sum - basisvals_u[3*ki]*dusum) * weights[kr] * fac_vw2;
		    basisDerivs_v[kr] = tmp_vw_v * weights[kr] * basisvals_u[3*ki];
		    basisDerivs_w[kr] = tmp_vw_w * weights[kr] * basisvals_u[3*ki];

		    basisDerivs_uu[kr] = (basisvals_u[3*ki+2]*sum2 + 2.0*basisvals_u[3*ki]*dusum*dusum
					  - basisvals_u[3*ki]*duusum*sum - 2.0*basisvals_u[3*ki+1]*dusum*sum) * weights[kr] * fac_vw3;
		    basisDerivs_uv[kr] = (basisvals_u[3*ki+1]*basisvals_v[3*kj+1]*sum2 + 2.0*basisvals_u[3*ki]*basisvals_v[3*kj]*dusum*dvsum
					  - basisvals_u[3*ki+1]*basisvals_v[3*kj]*dvsum*sum - basisvals_u[3*ki]*basisvals_v[3*kj+1]*dusum*sum
					  - basisvals_u[3*ki]*basisvals_v[3*kj]*duvsum*sum) * weights[kr] * fac_w3;
		    basisDerivs_uw[kr] = (basisvals_u[3*ki+1]*basisvals_w[3*kh+1]*sum2 + 2.0*basisvals_u[3*ki]*basisvals_w[3*kh]*dusum*dwsum
					  - basisvals_u[3*ki+1]*basisvals_w[3*kh]*dwsum*sum - basisvals_u[3*ki]*basisvals_w[3*kh+1]*dusum*sum
					  - basisvals_u[3*ki]*basisvals_w[3*kh]*duwsum*sum) * weights[kr] * fac_v3;
		    basisDerivs_vv[kr] = tmp_vw_vv * weights[kr] * basisvals_u[3*ki];
		    basisDerivs_vw[kr] = tmp_vw_vw * weights[kr] * basisvals_u[3*ki];
		    basisDerivs_ww[kr] = tmp_vw_ww * weights[kr] * basisvals_u[3*ki];
		}
	    }
	}
    }
    else
    {
	// Multiply basis values in the three parameter directions
	for (kh=0, kr=0; kh<worder; ++kh)
	    for (kj=0; kj<vorder; ++kj)
		for (ki=0; ki<uorder; ++ki, ++kr)
		{
		    basisValues[kr] = basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh];
		    basisDerivs_u[kr] = basisvals_u[3*ki+1]*basisvals_v[3*kj]*basisvals_w[3*kh];
		    basisDerivs_v[kr] = basisvals_u[3*ki]*basisvals_v[3*kj+1]*basisvals_w[3*kh];
		    basisDerivs_w[kr] = basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh+1];
		    basisDerivs_uu[kr] = basisvals_u[3*ki+2]*basisvals_v[3*kj]*basisvals_w[3*kh];
		    basisDerivs_uv[kr] = basisvals_u[3*ki+1]*basisvals_v[3*kj+1]*basisvals_w[3*kh];
		    basisDerivs_uw[kr] = basisvals_u[3*ki+1]*basisvals_v[3*kj]*basisvals_w[3*kh+1];
		    basisDerivs_vv[kr] = basisvals_u[3*ki]*basisvals_v[3*kj+2]*basisvals_w[3*kh];
		    basisDerivs_vw[kr] = basisvals_u[3*ki]*basisvals_v[3*kj+1]*basisvals_w[3*kh+1];
		    basisDerivs_ww[kr] = basisvals_u[3*ki]*basisvals_v[3*kj]*basisvals_w[3*kh+2];
		}
    }

}



} // namespace Go
