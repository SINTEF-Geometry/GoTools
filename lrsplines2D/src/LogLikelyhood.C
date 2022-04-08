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

#include "GoTools/lrsplines2D/LogLikelyhood.h"
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

//#define DEBUG

using std::vector;
namespace {
  const double TOLnu = 1.0e-10;
  const double TOLmu = 1.0e-7;
  const double TOLvar = 1.0e-8;
}

namespace { // anonymous namespace 

// distance function between two curves.  Used by the minimization algorithm
// initiated by ClosestPoint::closestPtCurves.
class optnuFun {
public:
    optnuFun(double mu, double var, const vector<double>& Y);
    
  inline void eval(double nu, double& f, double& df) const;
  inline double minPar() const;
  inline double maxPar() const;

private:
  double minpar_;
  double maxpar_;
  int n_;
  vector<double> delta_;
};

} // end anonymous namespace


//=============================================================================
  double LogLikelyhood::compute(const vector<double>& residual,
				double indegT, bool iterateT, double& llh2)
//=============================================================================
{
  int N = residual.size();
  double degreeT = indegT;

  double mean, variance;
  LogLikelyhood::getMeanAndVariance(residual, mean, variance);
  
  LogLikelyhood::getMeanAndVarianceIterated(residual, degreeT, mean, variance, iterateT);
#ifdef DEBUG
  std::cout << "Mean: " << mean << ", variance: " << variance << ", degreeT: " << degreeT << std::endl;
#endif

  double lgamma1 = boost::math::lgamma((degreeT+1)/2.0);
  double lgamma2 = boost::math::lgamma(degreeT/2.0);
  double stdd = sqrt(variance);

  double t5 = 0.0, t6 = 0.0;
  for (int ki=0; ki<N; ++ki)
    {
      double t5_0 = degreeT + residual[ki]*residual[ki];
      double t6_0 = degreeT + (residual[ki]-mean)*(residual[ki]-mean)/variance;
      t5 += log(t5_0);
      t6 += log(t6_0);
    }

  double llh = N*(lgamma1 - lgamma2 - 0.5*log(variance) +
		  0.5*degreeT*log(degreeT)) - 0.5*(degreeT+1)*t5;
  llh2 = N*(lgamma1 - lgamma2 - 0.5*log(variance) +
	    0.5*degreeT*log(degreeT)) - 0.5*(degreeT+1)*t6;

  return llh;
}

//=============================================================================
void LogLikelyhood::getMeanAndVariance(const vector<double>& residual,
				       double& mean, double& variance)
//=============================================================================
{
  int N = (int)residual.size();
  mean = 0.0;
  for (size_t ki=0; ki<residual.size(); ++ki)
    mean += residual[ki];
  mean /= (double)N;

  variance = 0.0;
  for (size_t ki=0; ki<residual.size(); ++ki)
    {
      double tmp = (residual[ki]-mean)*(residual[ki]-mean);
      tmp /= (double)(N-1);
      variance += tmp;
    }
}

//=============================================================================
void LogLikelyhood::getMeanAndVarianceIterated(const vector<double>& residual,
					       double& degreeT,
					       double& mean, double& variance,
					       bool iterateT)
//=============================================================================
{
  int N = (int)residual.size();

  int maxiter = 400;
  vector<double> cx(N), cxw(N), M(N), delta(N), w(N);
  for (int iter=0; iter<maxiter; ++iter)
    {
      double mean0 = mean, variance0 = variance;
      double sigma = sqrt(variance0);
      for (int ki=0; ki<N; ++ki)
	{
	  cx[ki] = residual[ki] - mean0;
	  M[ki] = cx[ki]/sigma;
	  delta[ki] = M[ki]*M[ki];
	  w[ki] = (1.0+degreeT)/(delta[ki]+degreeT);
	}

      mean = 0.0;
      double sumw = 0.0;
      for (int ki=0; ki<N; ++ki)
	{
	  mean += (residual[ki]*w[ki]);
	  sumw += w[ki];
	}
      mean /= sumw;

      variance = 0.0;
      for (int ki=0; ki<N; ++ki)
	{
	  cx[ki] = residual[ki] - mean;
	  cxw[ki] = cx[ki]*sqrt(w[ki]);
	  variance += (cxw[ki]*cxw[ki]/(double)(N-1));
	}

      if (fabs(mean-mean0) < TOLmu && fabs(variance-variance0) < TOLvar)
	break;
      
      if (iterateT /*&& iter%2 == 0*/)
	iterateDegreeT(residual, mean, variance, degreeT);

      int stop_break = 1;
    }
}

//=============================================================================
void LogLikelyhood::iterateDegreeT(const vector<double>& residual,
				   double mean, double variance,
				   double& degreeT)
//=============================================================================
{
  optnuFun opt(mean, variance, residual);
  double nu = degreeT;
  double minf = std::numeric_limits<double>::max();
  double minnu = nu;
  int maxiter = 500;
  for (int ka=0; ka<maxiter; ++ka)
    {
      double nu0= nu;
      double f, df;
      opt.eval(nu, f, df);
      if (fabs(f) < minf)
	{
	  minf = fabs(f);
	  minnu = nu;
	}
      if (fabs(df) < TOLnu)
	break;
      nu -= (f/df);
      nu = std::max(opt.minPar(), nu);
      nu = std::min(opt.maxPar(), nu);
      if (fabs(nu-nu0) < TOLnu)
	break;
    }
  
  degreeT = minnu;
}

  

namespace {

  //===========================================================================
  optnuFun::optnuFun(double mu, double var,  const vector<double>& Y)
  //===========================================================================
  {
    n_ = (int)Y.size();
    delta_.resize(n_);
    for (int ki=0; ki<n_; ++ki)
      delta_[ki] = (Y[ki]-mu)*(Y[ki]-mu)/var;
    minpar_ = 1.0;
    maxpar_ = 500.0;
  }
  
  
  //===========================================================================    
  void optnuFun::eval(double nu, double& f, double& df) const
  //===========================================================================
  {
    double nu2 = 0.5*nu;
    double pnu2 = 0.5*(1+nu);
    double pgamma = boost::math::polygamma(1,nu2);
    double ppgamma = boost::math::polygamma(1,pnu2);
    double psi = boost::math::digamma(nu2);
    double ppsi = boost::math::digamma(pnu2);
    
    vector<double> w(n_);
    for (int ki=0; ki<n_; ++ki)
      w[ki] = (1+nu)/(nu+delta_[ki]);

    double tmp = 0.0, tmpd = 0.0;
    for (int ki=0; ki<n_; ++ki)
      {
	double tmp2 = (log(w[ki]) - w[ki])/(double)n_;
	tmp += tmp2;

	double tmpd2 = (delta_[ki]-1)/(delta_[ki]+nu);
	tmpd += (tmpd2*tmpd2);
      }
   

    f = -psi + log(nu2) + tmp + 1 + ppsi - log(pnu2);
    df = -0.5*pgamma + (1.0/nu) + (tmpd/(n_*(1+nu))) + 0.5*ppgamma - (1.0/(1+nu));
    int stop_break = 1;
  }

  //===========================================================================    
  double optnuFun::minPar() const
  //===========================================================================
  {
    return minpar_;
  }

  //===========================================================================    
  double optnuFun::maxPar() const
  //===========================================================================
  {
    return maxpar_;
  }
};
