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

#include "GoTools/lrsplines2D/LRProjection.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/geometry/SplineSurface.h"
//#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/utils/QRFactorization.h"
#include "sislP.h"
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>

//#define DEBUG
//#define DEBUG2
//#define DEBUG3

using std::vector;
using std::pair;
using std::set;
using std::map;
using std::array;
using std::cout;
using std::endl;
using namespace Go;

//==============================================================================
double LRProjection::computeCoef(LRSplineSurface *srf, 
				 LRBSpline2D *bspl,
				 int proj_type, bool apply_smooth, double dlim,
				 Point& coef, int num_points)
//==============================================================================
{
  int del = (*bspl->supportedElementBegin())->getNmbValPrPoint();
  if (del == 0)
    del = srf->dimension() + 3;  // Parameter pair, point and distance

  bool proj_OK = true;
  double rad = -1.0;
  if (proj_type == 1)
    {
      // Collect data points
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      for (auto el=bspl->supportedElementBegin(); el!=bspl->supportedElementEnd(); ++el)
	{
	  vector<double> elem_data = (*el)->getDataPoints();
	  data.insert(data.end(), elem_data.begin(), elem_data.end());
	  max_dist = std::max(max_dist, (*el)->getMaxError());
	  av_dist += (*el)->getAccumulatedError();
	  nmb_out += (*el)->getNmbOutsideTol();
	  nmb += (*el)->nmbDataPoints();
	}
      av_dist /= (double)nmb;
  

      int deg1 = bspl->degree(XFIXED);
      int deg2 = bspl->degree(YFIXED);
      Point coef2 = bspl->Coef();
      if (nmb < deg1*deg2)
	{
	  coef = coef2;
	  proj_OK = false;
	}
      else
	{
	  // Perform projection
	  double smoothwgt = 1.0e-9; //1.0e-3; //1.0e-6;
	  try {
	  proj_OK = TPproject(bspl, data, del, smoothwgt, coef);
	  }
	  catch (...)
	    {
	      coef = coef2;
	    }
	  double dist = coef.dist(coef2);
	  if (dist > dlim)
	    {
	      proj_OK = false;
	    }
	}
    }
  else if (proj_type == 2)
    {
      // Radial basis functions with IDW
      int nmb_points = (num_points > 0) ? num_points : 15; //30;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      IDWproject(bspl, rad, data, del, coef);
    }
  else if (proj_type == 3)
    {
      // Linear approximation
      int nmb_points = (num_points > 0) ? num_points : 15; //8;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = LinProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 4)
    {
      // Bilinear approximation
      int nmb_points = (num_points > 0) ? num_points : 15; //8;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = BiLinProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 5)
    {
      // Quadratic approximation
      int nmb_points = (num_points > 0) ? num_points : 30; //15; //20; //30;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = QuadProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 6)
    {
      // Biquadratic approximation
      int nmb_points = (num_points > 0) ? num_points : 50; //25; //15; 
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = BiQuadProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 7)
    {
      // Cubic approximation
      int nmb_points = (num_points > 0) ? num_points : 50;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = CubicProject(bspl, data, del, rad, coef, apply_smooth);
    }

  else if (proj_type == 8)
    {
      // Cubic approximation
      int nmb_points = (num_points > 0) ? num_points : 60;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(srf, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = BiCubicProject(bspl, data, del, rad, coef, apply_smooth);
    }

  while (!proj_OK)
    {
      proj_type--;
      if (proj_type <= 1)
	proj_type = 2;
#ifdef DEBUG3
      std::cout << "Failure. Applying projection type " << proj_type << std::endl;
#endif
      proj_OK = computeCoef(srf, bspl, proj_type, apply_smooth, dlim, coef, num_points);
    }

  return rad;
}

//==============================================================================
void LRProjection::RDataSet(LRSplineSurface *srf, LRBSpline2D *bspl, int nmb_pts,
			    double& rad, vector<double>& data,
			    int del, int& nmb, int& nmb_out,
			    double& max_dist, double& av_dist)
//==============================================================================
{
  rad = -1.0;
  Point uv0 = bspl->getGrevilleParameter();
  vector<Element2D*> elems;
  vector<LRBSpline2D*> bsplines;
  bsplines.push_back(bspl);
  vector<double>  uv_dist;
  double umin = uv0[0], umax = uv0[0], vmin = uv0[1], vmax = uv0[1];
  double sf_umin = srf->startparam_u();
  double sf_umax = srf->endparam_u();
  double sf_vmin = srf->startparam_v();
  double sf_vmax = srf->endparam_v();
  size_t start_el = 0;
  int nmb_found = 0;
  for (size_t ki=0; ki<bsplines.size(); ++ki)
    {
      for (auto el=bsplines[ki]->supportedElementBegin();
	   el!=bsplines[ki]->supportedElementEnd(); ++el)
	{
	  auto it = std::find(elems.begin(), elems.end(), *el);
	  if (it == elems.end())
	    elems.push_back(*el);
	}

      for (size_t kj=start_el; kj<elems.size(); ++kj)
	{
	  vector<double> elem_data = elems[kj]->getDataPoints();
	  int nd = (int)elem_data.size()/del;
	  for (int ka=0; ka<nd; ++ka)
	    {
	      double dist = Utils::distance_squared(uv0.begin(), uv0.end(),
						    &elem_data[ka*del]);
	      uv_dist.push_back(sqrt(dist));
	    }
	  umin = std::min(umin, elems[kj]->umin()); 
	  umax = std::max(umax, elems[kj]->umax()); 
	  vmin = std::min(vmin, elems[kj]->vmin()); 
	  vmax = std::max(vmax, elems[kj]->vmax()); 
	}

	  std::sort(uv_dist.begin(), uv_dist.end());
	  if ((int)uv_dist.size() > nmb_pts)
	    {
	      rad = uv_dist[nmb_pts];
	      if ((umin <= uv0[0]-rad || umin == sf_umin) &&
		  (umax >= uv0[0]+rad || umax == sf_umax) &&
		  (vmin <= uv0[1]-rad || vmin == sf_vmin) &&
		  (vmax >= uv0[1]+rad || vmax == sf_vmax))
		break;
	    }


	  if (uv_dist.size() > 0)
	    rad = uv_dist[uv_dist.size()-1];
	  for (size_t kj=start_el; kj<elems.size(); ++kj)
	    {
	      for (auto b2=elems[kj]->supportBegin(); b2!=elems[kj]->supportEnd(); ++b2)
		{
		  auto itb = std::find(bsplines.begin(), bsplines.end(), *b2);
		  if (itb == bsplines.end())
		    bsplines.push_back(*b2);
		}
	    }
	  start_el = elems.size();
    }
  
  for (size_t ki=0; ki<elems.size(); ++ki)
    {
      vector<double> elem_data = elems[ki]->getDataPoints();
      data.insert(data.end(), elem_data.begin(), elem_data.end());
      max_dist = std::max(max_dist, elems[ki]->getMaxError());
      av_dist += elems[ki]->getAccumulatedError();
      nmb_out += elems[ki]->getNmbOutsideTol();
      nmb += elems[ki]->nmbDataPoints();
    }

  av_dist /= (double)nmb;
}

//==============================================================================
bool LRProjection::TPproject(LRBSpline2D *bspl,
			     vector<double>& data,
			     int del,
			     double smoothwgt,
			     Point& coef)
//==============================================================================
{
  // Fetch knot vectors of B-spline
  BSplineUniLR* uni_u = bspl->getUnivariate(XFIXED);
  BSplineUniLR* uni_v = bspl->getUnivariate(YFIXED);
  vector<double> knots1 = uni_u->getKnots();
  vector<double> knots2 = uni_v->getKnots();

  int startmult1 = uni_u->endmult(true);
  int endmult1 = uni_u->endmult(false);
  int startmult2 = uni_v->endmult(true);
  int endmult2 = uni_v->endmult(false);

  int order1 = uni_u->degree() + 1;
  int order2 = uni_v->degree() + 1;

  // Extend knot vectors to full multiplicity in endpoints
  double start1 = knots1[0];
  double end1 = knots1[knots1.size()-1];
  for (int ka=startmult1; ka<order1; ++ka)
    knots1.insert(knots1.begin(), start1);

  for (int ka=endmult1; ka<order1; ++ka)
    knots1.push_back(end1);

  double start2 = knots2[0];
  double end2 = knots2[knots2.size()-1];
  for (int ka=startmult2; ka<order2; ++ka)
    knots2.insert(knots2.begin(), start2);

  for (int ka=endmult2; ka<order2; ++ka)
    knots2.push_back(end2);

  // Extract point and parameter data
  int num_points = (int)data.size()/del;
  int dim = bspl->dimension();
  vector<double> param, points;
  param.reserve(2*num_points);
  points.reserve(dim*num_points);
  for (int ka=0; ka<(int)data.size(); ka+=del)
    {
      for (int kb=0; kb<2; ++kb)
	param.push_back(data[ka+kb]);
      for (int kb=0; kb<dim; ++kb)
	points.push_back(data[ka+kb+2]);
    }

  // Define initial spline surface
  int in1 = (int)knots1.size() - order1;
  int in2 = (int)knots2.size() - order2;
  vector<double> init_coefs(in1*in2*dim, 0.0);
  shared_ptr<SplineSurface> init_surf(new SplineSurface(in1, in2, order1, order2,
							&knots1[0], &knots2[0],
							&init_coefs[0], dim));

  // if (num_points > in1*in2)
  //   smoothwgt *= 0.01;
  
  // Approximation
  SmoothSurf approx;
  int seem[2];
  seem[0] = seem[1] = 0;
  vector<int> coef_known(in1*in2, 0);
  approx.attach(init_surf, seem, &coef_known[0]);

  double wgt1 = 0.0;
  double wgt3 = (std::min(order1,order2) < 4) ? 0.0 : 0.5*smoothwgt;
  double wgt2 = smoothwgt - wgt3 - wgt1;
  approx.setOptimize(wgt1, wgt2, wgt3);

  double approx_wgt = 1.0 - wgt1 - wgt2 - wgt3;
  vector<double> pnt_wgt(num_points, 1.0);
  approx.setLeastSquares(points, param, pnt_wgt, approx_wgt);

  shared_ptr<SplineSurface> approx_surf;
  try {
  approx.equationSolve(approx_surf);
  }
  catch (...)
    {
      coef = bspl->Coef();
      return false;
    }
  
  // Fetch coefficient corresponding to input B-spline
  int k1 = order1 - startmult1;
  int k2 = order2 - startmult2;
  int kk = k2*in1 + k1;
  coef = Point(approx_surf->coefs_begin()+kk*dim, approx_surf->coefs_begin()+(kk+1)*dim);

  return true;
}

//==============================================================================
void LRProjection::IDWproject(LRBSpline2D *bspl, double rad,
			      vector<double>& data, int del, Point& coef)
//==============================================================================
{
  double eps = 1.0e-10;
  Point greville = bspl->getGrevilleParameter();
  int dim = bspl->dimension();
  Point nom(dim);
  nom.setValue(0.0);
  double denom = 0.0;
  int pp = 2;
  int nmbd = (int)data.size()/del;
  for (int ka=0; ka<nmbd; ++ka)
    {
      Point curr(&data[ka*del+2], &data[ka*del+2+dim]);
      double tmp = sqrt(Utils::distance_squared(greville.begin(), greville.end(),
						&data[ka*del]));
      if (tmp < eps)
	{
	  coef = curr;
	  return;
	}
      double tmp2 = std::max(0.0, rad-tmp)/(rad*tmp);
      tmp2 = tmp2*tmp2;

      nom += tmp2*curr;
      denom += tmp2;
    }

  coef = nom/denom;
}



//==============================================================================
bool LRProjection::BiLinProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 1, tot_degree = 2;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, apply_smooth);
}

//==============================================================================
bool LRProjection::LinProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 1, tot_degree = 1;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, true);
}

//==============================================================================
bool LRProjection::BiQuadProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 2, tot_degree = 4;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, apply_smooth);
}

//==============================================================================
bool LRProjection::QuadProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 2, tot_degree = 2;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, apply_smooth);
}

//==============================================================================
bool LRProjection::BiCubicProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 3, tot_degree = 6;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, apply_smooth);
}

//==============================================================================
bool LRProjection::CubicProject(LRBSpline2D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 3, tot_degree = 3;
  return PolynomialProject(degree, tot_degree, bspl, data, del, rad,
			   coef, apply_smooth);
}

SISLSurf* interpolateSample(double *points, int dim, double *parvals1,
			    double *parvals2, int order1, int order2,
			    double *knots1, double *knots2)
{
  int kstat = 0;
  vector<int> type2(order2, 0);
  double *coef2 = NULL;
  int knlr = 0, knrc = 0;  // Open surface
  int kn1_2 = 0, kn2_2 = 0;
  s1891(parvals2, points, order1*dim, order2, 1, &type2[0],
	1, knots2, &coef2, &kn2_2, order2, knlr, knrc, &kstat);
  if (kstat < 0)
    return NULL;

  vector<double> coef2_2(dim*order1*order2);
  s6chpar(coef2, order1, order2, dim, &coef2_2[0]);
  vector<int> type1(order1, 0);
  double *coef1 = NULL;
  s1891(parvals1, &coef2_2[0], order2*dim, order1, 1, &type1[0],
	1, knots1, &coef1, &kn1_2, order1, knlr, knrc, &kstat);
  if (coef2 != NULL)
    free(coef2);
  if (kstat < 0)
    return NULL;
  vector<double> coef1_2(dim*order1*order2);
  s6chpar(coef1, order2, order1, dim, &coef1_2[0]);
  SISLSurf *sislsrf = newSurf(order1, order2, order1, order2, knots1,
			      knots2, &coef1_2[0], 1, dim, 1);
  if (coef1 != NULL)
    free(coef1);

#ifdef DEBUG3
  if (sislsrf != NULL)
    {
      int left1 = 0, left2 = 0;
      for (int kb=0; kb<order2; ++kb)
	{
	  for (int ka=0; ka<order1; ++ka)
	    {
	      double epar[2];
	      epar[0] = parvals1[ka];
	      epar[1] = parvals2[kb];
	      vector<double> der(dim);
	      s1424(sislsrf, 0, 0, epar, &left1, &left2, &der[0], &kstat);
	      double dd = s6dist(&der[0], &points[dim*(kb*order1+ka)], dim);
	      std::cout << "epar: " << epar[0] << " " << epar[1] << ", dd: " << dd << std::endl;
	    }
	}
    }
#endif	
  return sislsrf;
}


void polynomialTerms(double u, double v, int num1, int num2, int max_num,
		     vector<double>& terms)
{
  vector<double> u_terms(num1+1), v_terms(num2+1);
  u_terms[0] = v_terms[0] = 1.0;

  for (int ki=1; ki<=num1; ++ki)
    u_terms[ki] = u_terms[ki-1]*u;
  for (int ki=1; ki<=num2; ++ki)
     v_terms[ki] = v_terms[ki-1]*v;

  for (int kj=0, kr=0; kj<=num2; ++kj)
    for (int ki=0; ki<=num1; ++ki)
      {
	if (ki+kj <= max_num)
	  terms[kr++] = u_terms[ki]*v_terms[kj];
      }
}

//==============================================================================
bool LRProjection::PolynomialProject(int degree, int tot_degree, 
				     LRBSpline2D *bspl, vector<double>& data, 
				     int del, double rad, Point& coef,
				     bool apply_smooth)
//==============================================================================
{
  bool stat = true;
  shared_ptr<ParamSurface> dummy;
  double out_fac = 10.0;
  double out_rad = 1.25*rad;
  Point par = bspl->getGrevilleParameter();
  double u, v;
  int dim = bspl->dimension();
  int num_terms = (tot_degree == degree) ? (degree+1)*(degree+2)/2 :
    (degree+1)*(degree+1);
  vector<double> tmp(num_terms);
  int nmbd = (int)data.size()/del;
  double u1 = std::numeric_limits<double>::max();
  double u2 = std::numeric_limits<double>::lowest();
  double v1 = std::numeric_limits<double>::max();
  double v2 = std::numeric_limits<double>::lowest();
  
  // Compute eparameterization factor to let the parameter domain
  // reflect the geometry
  vector<Point> data_pts;
  int num_pt = 0;
  int num_test = 5;
  int num_id = 0;
  vector<pair<int, double> > distant_pts(num_test);
  for (int ka=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	{
	  if (num_id < num_test)
	    distant_pts[num_id++] = std::make_pair(ka, dist);
	  else
	    {
	      for (int kb=0; kb<num_test; ++kb)
		if (fabs(distant_pts[kb].second-out_rad) > (dist-out_rad))
		  {
		    distant_pts[kb] = std::make_pair(ka, dist);
		    break;
		  }
	    }
	  continue;
	}

      ++num_pt;
      
      u = data[ka*del];
      v = data[ka*del+1];
      u1 = std::min(u1, u);
      u2 = std::max(u2, u);
      v1 = std::min(v1, v);
      v2 = std::max(v2, v);

      Point pt(3);
      if (dim == 3)
	pt.setValue(&data[ka*del+2]);
      else
	pt.setValue(&data[ka*del]);
      data_pts.push_back(pt);
    }
#ifdef DEBUG2
  std::ofstream of("pol_data.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << data_pts.size() << std::endl;
  for (size_t ki=0; ki<data_pts.size(); ++ki)
    of << data_pts[ki] << std::endl;
#endif

  // Compute quadratic polynomial approximating the input points with a
  // distance to the B-spline coefficient less than the specified radius,
  // measured in parameter space
  // Define equation system using least squares
  vector<double> A(num_pt*num_terms, 0.0), b(dim*num_pt, 0.0);
  for (int ka=0, kr=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	continue;

      u = data[ka*del];
      v = data[ka*del+1];

      polynomialTerms(u, v, degree, degree, tot_degree, tmp);

      for (int kb=0; kb<num_terms; ++kb)
	A[kr*num_terms+kb] = tmp[kb];
      
      for (int kb=0; kb<dim; ++kb)
	b[kb*num_pt+kr] = data[ka*del+2+kb];
      ++kr;
    }

  vector<double> Q, R, x;
  QRFactorization::QRDecomp(A, num_terms, num_pt, Q, R);
  QRFactorization::QRSolve(Q, R, num_terms, num_pt, b, dim, x);
  
  if (apply_smooth)
    {
      // Define coefficient by evaluating the polynomial in the Greville
      // point
      coef = Point(dim);
      coef.setValue(0.0);
      polynomialTerms(par[0], par[1], degree, degree, tot_degree, tmp);
 
      for (int kb=0; kb<dim; ++kb)
	for (int kc=0; kc<num_terms; ++kc)
	  coef[kb] += x[kb*num_terms+kc]*tmp[kc];
    }
  else
    coef = Polynomial2Coef(x, degree, degree, tot_degree, u1, u2, v1, v2,
			   bspl, dummy);
      
  double maxdist = 0.0, avdist = 0.0;
  int num = 0;
  for (int ka=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      if (dist > rad)
	continue;
      u = data[ka*del];
      v = data[ka*del+1];
      
      Point pos(dim);
      pos.setValue(0.0);
      polynomialTerms(u, v, degree, degree, tot_degree, tmp);
 
      for (int kb=0; kb<dim; ++kb)
	for (int kc=0; kc<num_terms; ++kc)
	  pos[kb] += x[kb*num_terms+kc]*tmp[kc];
      double dd = sqrt(Utils::distance_squared(pos.begin(), pos.end(),
					       &data[ka*del+2]));
      maxdist = std::max(maxdist, dd);
      avdist += dd;
      num++;
    }
  avdist /= (double)num;
  Point coef2 = bspl->Coef();
  double coef_dist = coef.dist(coef2);

  if (num_id > 0)
    {
      double maxdist2 = 0.0, avdist2 = 0.0;
      double maxrad = rad;
      for (int kb=0; kb<num_id; ++kb)
	{
	  maxrad = std::max(maxrad, distant_pts[kb].second);
	  int ka = distant_pts[kb].first;
	  u = data[ka*del];
	  v = data[ka*del+1];
      
	  Point pos(dim);
	  pos.setValue(0.0);
	  polynomialTerms(u, v, degree, degree, tot_degree, tmp);
 
	  for (int kb=0; kb<dim; ++kb)
	    for (int kc=0; kc<num_terms; ++kc)
	      pos[kb] += x[kb*num_terms+kc]*tmp[kc];
	  double dd = sqrt(Utils::distance_squared(pos.begin(), pos.end(),
						   &data[ka*del+2]));
	  maxdist2 = std::max(maxdist, dd);
	  avdist2 += dd;
	}
      avdist2 /= (double)num_id;

#ifdef DEBUG3
      std::cout << "maxdist: " << maxdist << ", avdist: " << avdist << std::endl;
      std::cout << "maxdist2: " << maxdist << ", avdist2: " << avdist << std::endl;
#endif
#ifdef DEBUG3
      std::cout << "del u: " << bspl->umax()-bspl->umin() << ", del v: " << bspl->vmax()-bspl->vmin() << std::endl;
      std::cout << "nmb: " << num << ", rad: " << rad <<", maxdist: " << maxdist << ", avdist: " << avdist << ", coef_dist: " << coef_dist << std::endl << std::endl;
      if (coef_dist > 0.1)
	{
	  std::cout << "coef_dist: " << coef_dist << ", rad: " << rad;
	  if (coef_dist > 0.1)
	    std::cout << " High dist";
	    
	  std::cout << std::endl;
	}
#endif
      double facrad = maxrad/rad;
      double fac = out_fac*facrad*facrad;
      if (avdist2 > fac*avdist && maxdist2 > fac*maxdist)
	stat = false;
      //stat = true;
    }
  return stat;
}

//==============================================================================
Point LRProjection::Polynomial2Coef(vector<double>& pol, int degree1, int degree2,
				    int tot_degree,
				    double u1, double u2, double v1,
				    double v2, LRBSpline2D *bspl,
				    shared_ptr<ParamSurface> psurf)
//==============================================================================
{
  Point coef;

  int maxdeg = std::max(degree1,degree2);
  int mindeg = std::min(degree1,degree2);
  int num_terms = (tot_degree == maxdeg) ?
    maxdeg*(mindeg+1) - mindeg*(mindeg-1)/2 + 1 : (degree1+1)*(degree2+1);
  vector<double> tmp(num_terms);

  // Represent the polynomial by a Bezier surface
  // Identify inner knots
  BSplineUniLR* uni_u = bspl->getUnivariate(XFIXED);
  BSplineUniLR* uni_v = bspl->getUnivariate(YFIXED);
  double u, v;
  int dim = bspl->dimension();
  int startmult1 = bspl->endmult_u(true);
  int endmult1 = bspl->endmult_u(false);
  int startmult2 = bspl->endmult_v(true);
  int endmult2 = bspl->endmult_v(false);
  int order1 = uni_u->degree() + 1;
  int order2 = uni_v->degree() + 1;
  vector<double> kval1 = uni_u->getKnots();
  int kn1 = order1-startmult1-endmult1+1;
  vector<double> knots1(kn1);
  vector<double> kval2 = uni_v->getKnots();
  int kn2 = order2-startmult2-endmult2+1;
  vector<double> knots2(kn2);
  for (int ka=0; ka<kn1; ++ka)
    {
      knots1[ka] = kval1[ka+startmult1];
    }
  for (int ka=0; ka<kn2; ++ka)
    {
      knots2[ka] = kval2[ka+startmult2];
    }

  // Define spline space for B-spline surface interpolating the polynomial
  double u1_2 = bspl->umin();
  double u2_2 = bspl->umax();
  double v1_2 = bspl->vmin();
  double v2_2 = bspl->vmax();
  vector<double> knots1_0(2*order1);
  vector<double> knots2_0(2*order2);
  for (int ka=0; ka<order1; ++ka)
    {
      knots1_0[ka] = u1_2;
      knots1_0[order1+ka] = u2_2;
    }
  for (int ka=0; ka<order2; ++ka)
    {
      knots2_0[ka] = v1_2;
      knots2_0[order2+ka] = v2_2;
    }

  // Sample polynomial
  double u3 = std::max(u1, u1_2);
  double u4 = std::min(u2, u2_2);
  double v3 = std::max(v1, v1_2);
  double v4 = std::min(v2, v2_2);
  vector<double> points(order1*order2*dim, 0.0);
  vector<double> parvals1(order1), parvals2(order2);
  double udel = (u4 - u3)/(double)(order1-1);
  double vdel = (v4 - v3)/(double)(order2-1);
  int ka, kb, kc, kd;
  for (kb=0, v=v3; kb<order2; ++kb, v+=vdel)
    {
      parvals2[kb] = v;
      for (ka=0, u=u3; ka<order1; ++ka, u+=udel)
	{
	  parvals1[ka] = u;
	  polynomialTerms(u, v, degree1, degree2, tot_degree, tmp);

	  if (psurf.get())
	    {
	      Point pos = psurf->point(u, v);
	      for (kd=0; kd<dim; ++kd)
		points[dim*(kb*order1+ka)+kd] = pos[kd];
	    }
	  else
	    {
	      for (kd=0; kd<dim; ++kd)
		for (kc=0; kc<num_terms; ++kc)
		  points[dim*(kb*order1+ka)+kd] += pol[kd*num_terms+kc]*tmp[kc];
	    }
	}
    }

  // Interpolate
  SISLSurf *sislsrf = interpolateSample(&points[0], dim, &parvals1[0],
					&parvals2[0], order1, order2,
					&knots1_0[0], &knots2_0[0]);
  if (sislsrf == NULL)
    return coef;

  // Insert inner knots
  if (knots1.size() > 0 || knots2.size() > 0)
    {
      int kstat = 0;
      SISLSurf *sislsrf2 = NULL;
      s1025(sislsrf, &knots1[0], (int)knots1.size(),
	    &knots2[0], (int)knots2.size(), &sislsrf2, &kstat);
      if (sislsrf)
	freeSurf(sislsrf);
      if (kstat < 0)
	return coef;
      sislsrf = sislsrf2;
    }
  
  // Fetch coefficient
  int k1 = order1 - startmult1;
  int k2 = order2 - startmult2;
  int kk = k2*sislsrf->in1 + k1;
  coef = Point(sislsrf->ecoef+kk*dim, sislsrf->ecoef+(kk+1)*dim);

#ifdef DEBUG2
  std::cout << "Coef: " << coef << std::endl;
#endif
  // Free surface
  if (sislsrf != NULL)
    freeSurf(sislsrf);
  
  return coef;
}

bool next_combination(std::vector<size_t>& indices, size_t m) {
    size_t n = indices.size();
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        if (indices[i] < m - n + static_cast<size_t>(i)) {
            ++indices[i];
            for (size_t j = i + 1; j < n; ++j) {
                indices[j] = indices[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

//==============================================================================
Point LRProjection::Polynomial2Coef2(vector<double>& pol, int degree1, int degree2,
				     int tot_degree, LRBSpline2D *bspl,
				     shared_ptr<ParamSurface> psurf)
//==============================================================================
{
  int dim = bspl->dimension();
  BSplineUniLR* uni_u = bspl->getUnivariate(XFIXED);
  BSplineUniLR* uni_v = bspl->getUnivariate(YFIXED);
  vector<double> kval1 = uni_u->getKnots();
  vector<double> kval2 = uni_v->getKnots();
  int deg1 = bspl->degree(XFIXED);
  int deg2 = bspl->degree(YFIXED);

  Point L(dim);
  L.setValue(0.0);

  int ki, kj;
  vector<int> k_ix1(deg1), k_ix2(deg2);
  int db1 = 1, db2 = 1;
  for (int ki=1; ki<=deg1; ++ki)
    {
      db1 *= ki;
      k_ix1[ki-1] = ki;
    }
  for (int kj=1; kj<=deg2; ++kj)
    {
      db2 *= kj;
      k_ix2[kj-1] = kj;
    }
   
  int maxdeg = std::max(degree1,degree2);
  int mindeg = std::min(degree1,degree2);
  int num_terms = (tot_degree == maxdeg) ?
    maxdeg*(mindeg+1) - mindeg*(mindeg-1)/2 + 1 : (degree1+1)*(deg2+1);
  vector<double> tmp(num_terms);
  for (ki=1; ki<=deg1; ++ki)
    {
      vector<vector<int> > knot_ix1;
      vector<size_t> ind1(ki);
      std::iota(ind1.begin(), ind1.end(), 0);
      do {
	vector<int> curr_ix(ki);
	for (size_t kr=0; kr<ind1.size(); ++kr) 
	  curr_ix[kr] = k_ix1[ind1[kr]];
	knot_ix1.push_back(curr_ix);
      } while (next_combination(ind1, deg1));
      
      int kd1 = 1;
      for (int kc=0; kc<deg1; ++kc)
	kd1 *= ki;
      for (kj=1; kj<=deg2; ++kj)
	{
	  vector<vector<int> > knot_ix2;
	  vector<size_t> ind2(kj);
	  std::iota(ind2.begin(), ind2.end(), 0);
	  do {
	    vector<int> curr_ix(kj);
	    for (size_t kr=0; kr<ind2.size(); ++kr) 
	      curr_ix[kr] = k_ix2[ind2[kr]];
	    knot_ix2.push_back(curr_ix);
	  } while (next_combination(ind2, deg2));
	      
	  int kd2 = 1;
	  for (int kc=0; kc<deg2; ++kc)
	    kd2 *= kj;
	  int ie = deg1-ki+deg2-kj;
	  int sgn = (ie % 2 == 0) ? 1 : -1;
	  double fac2 = (double)(sgn*kd1*kd2)/(double)(db1*db2);
	  Point L2(dim);
	  L2.setValue(0.0);
	  for (size_t ka=0; ka<knot_ix1.size(); ++ka)
	    {
	      double t1 = 0.0;
	      for (size_t kc=0; kc<knot_ix1[ka].size(); ++kc)
		t1 += kval1[knot_ix1[ka][kc]];
	      t1 /= (double)knot_ix1[ka].size();
	      for (size_t kb=0; kb<knot_ix2.size(); ++kb)
		{
		  double t2 = 0.0;
		  for (size_t kc=0; kc<knot_ix2[kb].size(); ++kc)
		    t2 += kval2[knot_ix2[kb][kc]];
		  t2 /= (double)knot_ix2[kb].size();

		  Point val(dim);
		  if (psurf.get())
		    val = psurf->point(t1, t2);
		  else
		    {
		      polynomialTerms(t1, t2, degree1, degree2,
				      tot_degree, tmp);
		      val.setValue(0.0);
		      for (int kd=0; kd<dim; ++kd)
			for (int kc=0; kc<num_terms; ++kc)
			  val[kd] += pol[kd*num_terms+kc]*tmp[kc];
		    }
		    L2 += val;
		}
	    }
	  L2 *= fac2;
	  L += L2;
	}
    }
  
  return L;
}

