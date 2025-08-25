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
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/Utils.h"

//***************************************************************************
//
// Implementation file of the free function closestPoint defined in
// SplineVolume.h/
//
//***************************************************************************

using namespace std;
using namespace Go;

// // Anonymous namespace
// namespace {
//   const double DZERO = (double)0.0;
//   const double TOL = 1.0e-17; //1.0e-16;
//   const double REL_COMP_RES = 0.000000000000001;
//   const double ANGULAR_TOLERANCE = 0.01;
//   const double SINGULAR = 1.0e-16;
// }

// namespace { // anonymous namespace 

// // distance function between two curves.  Used by the minimization algorithm
// // initiated by ClosestPoint::closestPtCurves.
// class VolPntDistFun {
// public:
//     VolPntDistFun(const ParamVolume* vol, 
// 		  const Point& pt,
// 		  const double* const minpar = 0,
// 		  const double* const maxpar = 0);
    
//     inline double operator()(const double* arg) const;
//     inline double grad(const double* arg, double* res) const;
//     inline double minPar(int pardir) const;
//     inline double maxPar(int pardir) const;

// private:
//     double minpar_[3];
//     double maxpar_[3];
//     const ParamVolume * const vol_;
//     const Point pt_;
//     mutable Point p1_, p2_, d_;
//     mutable vector<Point> pvec_;
// };


// } // end anonymous namespace

namespace Go
{

// //===========================================================================
// void  SplineVolume::closestPoint(const Point& pt,
// 				 double&        clo_u,
// 				 double&        clo_v, 
// 				 double&        clo_w, 
// 				 Point&         clo_pt,
// 				 double&        clo_dist,
// 				 double         epsilon,
// 				 double   *seed) const
// //===========================================================================
// {
//     // Iteration 
//     double start_par[3], par[3], minpar[3], maxpar[3];
//     double dist;
//     double seed_dist = std::numeric_limits<double>::max();
//     const Array<double,6> domain = parameterSpan();
//     minpar[0] = domain[0];
//     minpar[1] = domain[2];
//     minpar[2] = domain[4];
//     maxpar[0] = domain[1];
//     maxpar[1] = domain[3];
//     maxpar[2] = domain[5];
//     if (seed)
//     {
// 	start_par[0] = seed[0];
// 	start_par[1] = seed[1];
// 	start_par[2] = seed[2];
//     }
//     else
//     {
//       seed_dist = getSeed(pt, start_par);
//       for (int ki=0; ki<3; ++ki)
// 	{
// 	  if (numCoefs(ki) <= 2 && seed_dist >= TOL)
// 	    start_par[ki] = 0.5*(minpar[ki]+maxpar[ki]);
// 	}
//     }

//     // Check if the volume is closed in any direction
//     int closed[3];
//     for (int ki=0; ki<3; ++ki)
//       closed[ki] = volumePeriodicity(ki, epsilon);

//     if (seed_dist < TOL)
//       {
// 	// Avoid closest point iteration
// 	clo_dist = seed_dist;
// 	clo_u = start_par[0];
// 	clo_v = start_par[1];
// 	clo_w = start_par[2];
// 	point(clo_pt, clo_u, clo_v, clo_w);
// 	return;
//       }

//     VolPntDistFun distfun(this, pt, minpar, maxpar);
//     FunctionMinimizer<VolPntDistFun> funmin(3, distfun, start_par, TOL);
//     try {
//       minimise_conjugated_gradient(funmin);//, 3); // number of iterations in each cycle
//     } 
//     catch (...)
//       {
// 	MESSAGE("SplineVolume::closestPoint, minimize_conjugate_gradient failed");
//       }

//     dist = sqrt(funmin.fval());
//     clo_u = par[0] = funmin.getPar(0);
//     clo_v = par[1] = funmin.getPar(1);
//     clo_w = par[2] = funmin.getPar(2);
    
//     double fac = 100.0;
//     if (dist > fac*epsilon)
//       {
// 	for (int ki=0; ki<3; ++ki)
// 	  {
// 	    if (closed[ki] >= 0)
// 	      {
// 		if (fabs(par[ki]-minpar[ki]) < fac*epsilon)
// 		  start_par[ki] = maxpar[ki];
// 		else if (fabs(maxpar[ki]-par[ki]) < fac*epsilon)
// 		  start_par[ki] = minpar[ki];
// 		else
// 		  continue;

// 		FunctionMinimizer<VolPntDistFun> funmin2(3, distfun, start_par, TOL);
// 		minimise_conjugated_gradient(funmin2);
// 		double dist2 = sqrt(funmin2.fval());
// 		if (dist2 < dist)
// 		  {
// 		    dist = dist2;
// 		    clo_u = funmin2.getPar(0);
// 		    clo_v = funmin2.getPar(1);
// 		    clo_w = funmin2.getPar(2);
// 		  }
// 	      }
// 	  }
//       }

//     point(clo_pt, clo_u, clo_v, clo_w);
//     clo_dist = pt.dist(clo_pt);
// }


//===========================================================================
    double  SplineVolume::getSeed(const Point& pt, double par[]) const
//===========================================================================
    {
    // Find the coefficient closest to this point and return the Greville
    // values of this coefficients
    int nn1 = numCoefs(0);
    int nn2 = numCoefs(1);
    int nn3 = numCoefs(2);

    vector<double>::const_iterator coefs = coefs_begin();
    int ki, kj, kr;
    int k1min, k2min, k3min;
    double dist;
    double dmin = std::numeric_limits<double>::max(); //1.0e15; // Huge
    for (kr = 0; kr < nn3; ++kr) {
        for (kj = 0; kj < nn2; ++kj) {
            for (ki = 0; ki < nn1; ++ki, coefs += dim_) {
                dist = 0.0;
                for (int d = 0; d < dim_; ++d) {
                    double tmp = coefs[d] - pt[d];
                    dist += tmp * tmp;
                }
                if (dist < dmin) {
                    dmin = dist;
                    k1min = ki;
                    k2min = kj;
                    k3min = kr;
                }
            }
        }
    }

    par[0] = (nn1 <= 2) ? 0.5*(startparam(0)+endparam(0)) :
      basis_u_.grevilleParameter(k1min);
    par[1] = (nn2 <= 2) ? 0.5*(startparam(1)+endparam(1)) :
      basis_v_.grevilleParameter(k2min);
    par[2] = (nn3 <= 2) ? 0.5*(startparam(2)+endparam(2)) :
      basis_w_.grevilleParameter(k3min);

    Point pos;
    point(pos, par[0], par[1], par[2]);
    return pt.dist(pos);
}

//===========================================================================
int  SplineVolume::closestCorner(const Point& pt,
				 double&        upar,
				 double&        vpar, 
				 double&        wpar, 
				 Point&         corner,
				 double&        dist) const
//===========================================================================
{
  const Array<double,6> domain =  parameterSpan();
  // double u1, v1, w2;
  int ki, kj, kr, kh;
  int min_idx = -1;
  dist = 1.0e8;

  // Find closest corner
  for (kr=0, kh=0; kr<2; ++kr)
    for (kj=0; kj<2; ++kj)
      for (ki=0; ki<2; ++ki, ++kh)
	{
	  Point tmp;
	  point(tmp, domain[ki], domain[kj], domain[kr]);
	  double d1 = pt.dist(tmp);
	  if (d1 < dist)
	    {
	      upar = domain[ki];
	      vpar = domain[kj];
	      wpar = domain[kr];
	      dist = d1;
	      corner = tmp;
	      min_idx = kh;
	    }
	}

  kr = min_idx/4;
  kj = (min_idx - 4*kr)/2;
  ki = min_idx%2;

  int kn[3];
  for (kh=0; kh<3; ++kh)
    kn[kh] = numCoefs(kh);

  return kr*kn[0]*kn[1]*(kn[2]-1) + kj*kn[0]*(kn[1]-1) + ki*(kn[0]-1);
}


} // namespace Go

// namespace {

// //===========================================================================
// VolPntDistFun::VolPntDistFun(const ParamVolume* vol, 
// 		       const Point& pt,
// 		       const double* const minpar,
// 		       const double* const maxpar)
// //===========================================================================
//     : vol_(vol), pt_(pt), pvec_(4)
// {
//     const Array<double,6> domain = vol_->parameterSpan();
//     if (!minpar) {
// 	minpar_[0] = domain[0];
// 	minpar_[1] = domain[2];
// 	minpar_[2] = domain[4];
//     } else {
// 	minpar_[0] = minpar[0];
// 	minpar_[1] = minpar[1];
// 	minpar_[2] = minpar[2];
//     }
//     if (!maxpar) {
// 	maxpar_[0] = domain[1];
// 	maxpar_[1] = domain[3];
// 	maxpar_[2] = domain[5];
//     } else {
// 	maxpar_[0] = maxpar[0];
// 	maxpar_[1] = maxpar[1];
// 	maxpar_[2] = maxpar[2];
//     }
// }

// //===========================================================================    
// double VolPntDistFun::operator()(const double* arg) const
// //===========================================================================
// {
//     vol_->point(p1_, arg[0], arg[1], arg[2]);
//     return p1_.dist2(pt_);
// }

// //===========================================================================
// double VolPntDistFun::grad(const double* arg, double* res) const
// //===========================================================================
// {
//     vol_->point(pvec_, arg[0], arg[1], arg[2], 1);
//     d_ = pvec_[0] - pt_;
    
//     res[0] = 2 * d_ * pvec_[1];
//     res[1] = 2 * d_ * pvec_[2];
//     res[2] = 2 * d_ * pvec_[3];
    
//     return d_.length2();
// }

// //===========================================================================
// double VolPntDistFun::minPar(int pardir) const
// //===========================================================================
// {
//   if (pardir < 0 || pardir > 2)
//     THROW("Parameter direction out of range");
//   return minpar_[pardir];
// }

// //===========================================================================
// double VolPntDistFun:: maxPar(int pardir) const
// //===========================================================================
// {
//   if (pardir < 0 || pardir > 2)
//     THROW("Parameter direction out of range");
//   return maxpar_[pardir];
// }

// };
