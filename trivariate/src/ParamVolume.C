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

#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"

using std::vector;
using namespace Go;

// Anonymous namespace
namespace {
  const double DZERO = (double)0.0;
  const double TOL = 1.0e-17; //1.0e-16;
  const double REL_COMP_RES = 0.000000000000001;
  const double ANGULAR_TOLERANCE = 0.01;
  const double SINGULAR = 1.0e-16;
}

namespace { // anonymous namespace 

// distance function between two curves.  Used by the minimization algorithm
// initiated by ClosestPoint::closestPtCurves.
class VolPntDistFun {
public:
    VolPntDistFun(const ParamVolume* vol, 
		  const Point& pt,
		  const double* const minpar = 0,
		  const double* const maxpar = 0);
    
    inline double operator()(const double* arg) const;
    inline double grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const;
    inline double maxPar(int pardir) const;

private:
    double minpar_[3];
    double maxpar_[3];
    const ParamVolume * const vol_;
    const Point pt_;
    mutable Point p1_, p2_, d_;
    mutable vector<Point> pvec_;
};


} // end anonymous namespace

namespace Go
{

//===========================================================================
ParamVolume::~ParamVolume()

//===========================================================================
{
}

//===========================================================================
void  ParamVolume::closestPoint(const Point& pt,
				double&        clo_u,
				double&        clo_v, 
				double&        clo_w, 
				Point&         clo_pt,
				double&        clo_dist,
				double         epsilon,
				double   *seed) const
//===========================================================================
{
    // Iteration 
    double start_par[3], par[3], minpar[3], maxpar[3];
    double dist;
    double seed_dist = std::numeric_limits<double>::max();
    const Array<double,6> domain = parameterSpan();
    minpar[0] = domain[0];
    minpar[1] = domain[2];
    minpar[2] = domain[4];
    maxpar[0] = domain[1];
    maxpar[1] = domain[3];
    maxpar[2] = domain[5];
    if (seed)
    {
	start_par[0] = seed[0];
	start_par[1] = seed[1];
	start_par[2] = seed[2];
    }
    else
    {
      seed_dist = getSeed(pt, start_par);
    }

    // Check if the volume is closed in any direction
    int closed[3];
    for (int ki=0; ki<3; ++ki)
      closed[ki] = volumePeriodicity(ki, epsilon);

    if (seed_dist < TOL)
      {
	// Avoid closest point iteration
	clo_dist = seed_dist;
	clo_u = start_par[0];
	clo_v = start_par[1];
	clo_w = start_par[2];
	point(clo_pt, clo_u, clo_v, clo_w);
	return;
      }

    VolPntDistFun distfun(this, pt, minpar, maxpar);
    FunctionMinimizer<VolPntDistFun> funmin(3, distfun, start_par, TOL);
    try {
      minimise_conjugated_gradient(funmin);//, 3); // number of iterations in each cycle
    } 
    catch (...)
      {
	MESSAGE("SplineVolume::closestPoint, minimize_conjugate_gradient failed");
      }

    dist = sqrt(funmin.fval());
    clo_u = par[0] = funmin.getPar(0);
    clo_v = par[1] = funmin.getPar(1);
    clo_w = par[2] = funmin.getPar(2);
    
    double fac = 100.0;
    if (dist > fac*epsilon)
      {
	for (int ki=0; ki<3; ++ki)
	  {
	    if (closed[ki] >= 0)
	      {
		if (fabs(par[ki]-minpar[ki]) < fac*epsilon)
		  start_par[ki] = maxpar[ki];
		else if (fabs(maxpar[ki]-par[ki]) < fac*epsilon)
		  start_par[ki] = minpar[ki];
		else
		  continue;

		FunctionMinimizer<VolPntDistFun> funmin2(3, distfun, start_par, TOL);
		minimise_conjugated_gradient(funmin2);
		double dist2 = sqrt(funmin2.fval());
		if (dist2 < dist)
		  {
		    dist = dist2;
		    clo_u = funmin2.getPar(0);
		    clo_v = funmin2.getPar(1);
		    clo_w = funmin2.getPar(2);
		  }
	      }
	  }
      }

    point(clo_pt, clo_u, clo_v, clo_w);
    clo_dist = pt.dist(clo_pt);
}

//===========================================================================
void ParamVolume::estimateVolSize(double& u_size, double& v_size, double& w_size,
				  int u_nmb, int v_nmb, int w_nmb)

//===========================================================================
{
  Array<double,6> dom = parameterSpan();
  double del_u = (dom[1] - dom[0])/(double)(u_nmb-1);
  double del_v = (dom[3] - dom[2])/(double)(v_nmb-1);
  double del_w = (dom[5] - dom[4])/(double)(w_nmb-1);

  int ki, kj, kr;
  double upar, vpar, wpar;
  vector<Point> pts(u_nmb*v_nmb*w_nmb);
  for (kr=0, wpar=dom[4]; kr<w_nmb; ++kr, wpar+=del_w)
    for (kj=0, vpar=dom[2]; kj<v_nmb; ++kj, vpar+=del_v)
      for (ki=0, upar=dom[0]; ki<u_nmb; ++ki, upar+=del_u)
	{
	  Point pos;
	  point(pos, upar, vpar, wpar);
	  pts[(kr*v_nmb+kj)*u_nmb+ki] = pos;
	}

  double acc_u=0.0, acc_v=0.0, acc_w=0.0;
  for (kr=0; kr<w_nmb; ++kr)
    for (kj=0; kj<v_nmb; ++kj)
      for (ki=1; ki<u_nmb; ++ki)
	acc_u += pts[(kr*w_nmb+kj)*u_nmb+ki-1].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_u /= (double)(v_nmb*w_nmb);
  
  for (kr=0; kr<w_nmb; ++kr)
    for (ki=0; ki<u_nmb; ++ki)
      for (kj=1; kj<v_nmb; ++kj)
	acc_v += pts[(kr*w_nmb+kj-1)*u_nmb+ki].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_v /= (double)(u_nmb*w_nmb);
  
  for (kj=0; kj<v_nmb; ++kj)
    for (ki=0; ki<u_nmb; ++ki)
      for (kr=1; kr<w_nmb; ++kr)
	acc_w += pts[((kr-1)*w_nmb+kj)*u_nmb+ki].dist(pts[(kr*w_nmb+kj)*u_nmb+ki]);
  acc_w /= (double)(u_nmb*v_nmb);
  
  u_size = acc_u;
  v_size = acc_v;
  w_size = acc_w;
}

//===========================================================================
  int ParamVolume::volumePeriodicity(int pardir, double epsilon) const
//===========================================================================
  {
    int nmb_sample = 10;
    if (pardir < 0 || pardir > 2)
      return -1;
    
    Array<double,6> dom = parameterSpan();
    int ix1 = (pardir+1)%3;
    int ix2 = (pardir+2)%3;
    double del1 = (dom[2*ix1+1] - dom[2*ix1])/(double)(nmb_sample-1);
    double del2 = (dom[2*ix2+1] - dom[2*ix2])/(double)(nmb_sample-1);

    double par[3];
    int ka, kb, kc;
    for (ka=0, par[ix1]=dom[2*ix1]; ka<nmb_sample; ++ka, par[ix1]+=del1)
      for (kb=0, par[ix2]=dom[2*ix2]; kb<nmb_sample; ++kb, par[ix2]+=del2)
	{
	  Point pos[2];
	  for (kc=0, par[pardir]=dom[2*pardir]; kc<2; ++kc, par[pardir]=dom[2*pardir+1])
	    point(pos[kc], par[0], par[1], par[2]);

	  double dist = pos[0].dist(pos[1]);
	  if (dist > epsilon)
	    return -1;
	}
    return 0;
  }

//===========================================================================
  double ParamVolume::getSeed(const Point& pt, double par[]) const
//===========================================================================
  {
    Array<double,6> dom = parameterSpan();
    VolumeTools::volume_seedfind(pt, *this, dom, par[0], par[1], par[2]);

    Point pt2;
    point(pt2, par[0], par[1], par[2]);
    return pt.dist(pt2);
  }

} // namespace Go


namespace {

//===========================================================================
VolPntDistFun::VolPntDistFun(const ParamVolume* vol, 
		       const Point& pt,
		       const double* const minpar,
		       const double* const maxpar)
//===========================================================================
    : vol_(vol), pt_(pt), pvec_(4)
{
    const Array<double,6> domain = vol_->parameterSpan();
    if (!minpar) {
	minpar_[0] = domain[0];
	minpar_[1] = domain[2];
	minpar_[2] = domain[4];
    } else {
	minpar_[0] = minpar[0];
	minpar_[1] = minpar[1];
	minpar_[2] = minpar[2];
    }
    if (!maxpar) {
	maxpar_[0] = domain[1];
	maxpar_[1] = domain[3];
	maxpar_[2] = domain[5];
    } else {
	maxpar_[0] = maxpar[0];
	maxpar_[1] = maxpar[1];
	maxpar_[2] = maxpar[2];
    }
}

//===========================================================================    
double VolPntDistFun::operator()(const double* arg) const
//===========================================================================
{
    vol_->point(p1_, arg[0], arg[1], arg[2]);
    return p1_.dist2(pt_);
}

//===========================================================================
double VolPntDistFun::grad(const double* arg, double* res) const
//===========================================================================
{
    vol_->point(pvec_, arg[0], arg[1], arg[2], 1);
    d_ = pvec_[0] - pt_;
    
    res[0] = 2 * d_ * pvec_[1];
    res[1] = 2 * d_ * pvec_[2];
    res[2] = 2 * d_ * pvec_[3];
    
    return d_.length2();
}

//===========================================================================
double VolPntDistFun::minPar(int pardir) const
//===========================================================================
{
  if (pardir < 0 || pardir > 2)
    THROW("Parameter direction out of range");
  return minpar_[pardir];
}

//===========================================================================
double VolPntDistFun:: maxPar(int pardir) const
//===========================================================================
{
  if (pardir < 0 || pardir > 2)
    THROW("Parameter direction out of range");
  return maxpar_[pardir];
}

};
