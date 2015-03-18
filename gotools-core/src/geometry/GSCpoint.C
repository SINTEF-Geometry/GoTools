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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineUtils.h"
#include <memory>


using std::vector;


namespace Go
{

//===========================================================================
void SplineCurve::point(Point& result, double tpar) const
//===========================================================================
{
    if (result.dimension() != dim_)
	result.resize(dim_);

    // Take care of the rational case
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    int kdim = dim_ + (rational_ ? 1 : 0);

    // Make temporary storage for the basis values and a temporary
    // computation cache.
    std::vector<double> b0(basis_.order());
    std::vector<double> temp(kdim, 0.0);

    // Compute the basis values and get some data about the spline spaces
    basis_.computeBasisValues(tpar, &b0[0]);
    int left = basis_.lastKnotInterval();
    int order = basis_.order();

    // Compute the tensor product value
    int coefind = left-order+1;
    for (int ii = 0; ii < order; ++ii) {
	for (int dd = 0; dd < kdim; ++dd) {
	    temp[dd] += b0[ii]*co[coefind*kdim + dd];
	}
	coefind += 1;
    }

    // Copy from temp to result
    if (rational_) {
	for (int dd = 0; dd < dim_; ++dd) {
	    result[dd] = temp[dd]/temp[kdim-1];
	}
    } else {
	for (int dd = 0; dd < dim_; ++dd) {
	    result[dd] = temp[dd];
	}
    }
}


//===========================================================================
void
SplineCurve::point(std::vector<Point>& result, double tpar,
		   int derivs, bool from_right) const
//===========================================================================
{
    double resolution = DEFAULT_PARAMETER_EPSILON; //1.0e-12;
    DEBUG_ERROR_IF(derivs < 0, "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1);
    int rsz = (int)result.size();
    DEBUG_ERROR_IF(rsz < totpts, "The vector of points must have sufficient size.");
    int i;
    for (i = 0; i < totpts; ++i)
	if (result[i].dimension() != dim_)
	    result[i].resize(dim_);

    if (derivs == 0) {
	point(result[0], tpar);
	return;
    }

    // Take care of the rational case
    const std::vector<double>& co = rational_ ? rcoefs_ : coefs_;
    int kdim = dim_ + (rational_ ? 1 : 0);

    // Make temporary storage for the basis values and a temporary
    // computation cache.
    std::vector<double> b0(basis_.order() * (derivs+1));
    std::vector<double> temp(totpts*kdim, 0.0);

    // Compute the basis values and get some data about the spline spaces
    from_right |= (tpar - startparam() < resolution);
    if (from_right)
	basis_.computeBasisValues(tpar, &b0[0], derivs);
    else { // @@sbr By far the best solution, but a solution.
	shared_ptr<ParamCurve> temp_crv(subCurve(startparam(), tpar));
	temp_crv->point(result, tpar, derivs);
	return;
	// 	basis_.computeBasisValuesLeft(tpar, &b0[0], derivs);
    }

    int left = basis_.lastKnotInterval();
    int order = basis_.order();

    // Compute the tensor product value
    int coefind = left-order+1;
    if ((!from_right) && (basis_.begin()[left] == tpar))
	--coefind; // Returned basis values are one to the left.
    for (int ii = 0; ii < order; ++ii) {
	for (int dd = 0; dd < kdim; ++dd) {
	    for (int dercount = 0; dercount < totpts; ++dercount) {
		temp[dercount*kdim + dd]
		    += b0[dercount + ii*totpts]*co[coefind*kdim + dd];
	    }
	}
	coefind += 1;
    }

    // Copy from temp to result
    if (rational_) {
	std::vector<double> restmp(totpts*dim_);
	SplineUtils::curve_ratder(&temp[0], dim_, derivs, &restmp[0]);
	for (int i = 0; i < totpts; ++i) {
	    for (int dd = 0; dd < dim_; ++dd) {
		result[i][dd] = restmp[i*dim_ + dd];
	    }
	}
    } else {
	for (int i = 0; i < totpts; ++i) {
	    for (int dd = 0; dd < dim_; ++dd) {
		result[i][dd] = temp[i*dim_ + dd];
	    }
	}
    }
}



//===========================================================================
void SplineCurve::computeBasis(double param, 
			       std::vector<double>& basisValues,
			       std::vector<double>& basisDerivs) const
//===========================================================================
{
  int ord = basis_.order();

  basisValues.resize(ord);
  basisDerivs.resize(ord);

  std::vector<double> basisvals(2 * basis_.order());
  basis_.computeBasisValues(param, &basisvals[0], 1);

  if (rational_)
    {
      int i, pos = (dim_ + 1) * (basis_.lastKnotInterval() - ord + 1) + dim_;

      double w_func = 0.0;
      double w_der = 0.0;
      for (i = 0; i < ord; ++i, pos += dim_ + 1)
	{
	  double w = rcoefs_[pos];
	  w_func += w * basisvals[i * 2];
	  w_der += w * basisvals[i * 2 + 1];
	}
      pos = (dim_ + 1) * (basis_.lastKnotInterval() - ord + 1) + dim_;
      double w_func_2 = w_func * w_func;
      for (i = 0; i < ord; ++i, pos += dim_ + 1)
	{
          double w = rcoefs_[pos];
	  basisValues[i] = basisvals[i*2] * w / w_func;
	  basisDerivs[i] = (basisvals[i*2 + 1] * w_func - basisvals[i*2] * w_der) * w / w_func_2;
	}
    }
  else
    {
      for (int i = 0; i < ord; ++i)
	{
	  basisValues[i] = basisvals[i*2];
	  basisDerivs[i] = basisvals[i*2 + 1];
	}
    }
}





//===========================================================================
void SplineCurve::gridEvaluator(std::vector<double>& points,
				const std::vector<double>& param) const
//===========================================================================
{
    int num = (int)param.size();
    int kk = basis_.order();
    int kdim = dim_;
    vector<double> basisvals(num * kk);
    vector<int>    knotinter(num * kk);

    basis_.computeBasisValues(&param[0], &param[0]+param.size(),
				&basisvals[0], &knotinter[0], 0);

    points.resize(num * dim_);
    double *result = &points[0];


    const double *scoef;
    if (rational_) {
	scoef = &rcoefs_[0];
	kdim +=1;
    } else {
	scoef = &coefs_[0];
    }
    
    std::vector<double> temp(kdim, 0.0);
    double *b0;
    int ki;
    for (ki=0, b0=&basisvals[0]; ki<num; ++ki, result+=dim_, b0+=kk)
      {
	std::fill(temp.begin(), temp.end(), 0.0);
	int left = knotinter[ki] - kk + 1;
	for (int ii=0; ii<kk; ++ii)
	  {
	    for (int dd=0; dd<kdim; ++dd)
	      temp[dd] += b0[ii]*scoef[left*kdim + dd];
	    left++;

	    // Copy from temp to result
	    if (rational_) {
	      for (int dd = 0; dd < dim_; ++dd) {
		result[dd] = temp[dd]/temp[dim_];
	      }
	    } else {
	      for (int dd = 0; dd < dim_; ++dd) {
		result[dd] = temp[dd];
	      }
	    }
	  }
      }
}


} // namespace Go


