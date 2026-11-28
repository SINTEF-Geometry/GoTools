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

#include "GoTools/lrsplines3D/LRProjection3D.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/geometry/SplineSurface.h"
//#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/trivariate/SmoothVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "GoTools/utils/QRFactorization.h"
#include "sislP.h"
#include <cmath>
#include <iostream>
#include <fstream>

//#define DEBUG
//#define DEBUG2

using std::vector;
using std::pair;
using std::set;
using std::map;
using std::array;
using std::cout;
using std::endl;
using namespace Go;

//==============================================================================
double LRProjection3D::computeCoef(LRSplineVolume *vol,
				   LRBSpline3D *bspl,
				   int proj_type, bool apply_smooth, double dlim,
				   Point& coef, int num_points)
//==============================================================================
{
  int del = (*bspl->supportedElementBegin())->getNmbValPrPoint();
  if (del == 0)
    del = vol->dimension() + 3;  // Parameter pair, point and distance

  dlim *= 10.0;

  // if (num_points <= 0)
  //   num_points = 1000;
  
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
  

      int deg1 = bspl->degree(XDIR);
      int deg2 = bspl->degree(YDIR);
      int deg3 = bspl->degree(ZDIR);
      Point coef2 = bspl->Coef();
      if (nmb < deg1*deg2*deg3)
	{
	  coef = coef2;
	  proj_OK = false;
	}
      else
	{
	  // Perform projection
	  double smoothwgt = 1.0e-12; //1.0e-6; //1.0e-3; //1.0e-9;
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
      int nmb_points = (num_points > 0) ? num_points : 50; //15; //30;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      IDWproject(bspl, rad, data, del, coef);
    }
  else if (proj_type == 3)
    {
      // Linear approximation
      int nmb_points = (num_points > 0) ? num_points : 50; //15; //8;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = LinProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 4)
    {
      // Trilinear approximation
      int nmb_points = (num_points > 0) ? num_points : 70; //15; //8;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = TriLinProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 5)
    {
      // Quadratic approximation
      int nmb_points = (num_points > 0) ? num_points : 80; //30; //15; //20; //30;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = QuadProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 6)
    {
      // Triquadratic approximation
      int nmb_points = (num_points > 0) ? num_points : 90; //25; //15; 
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = TriQuadProject(bspl, data, del, rad, coef, apply_smooth);
    }
  else if (proj_type == 7)
    {
      // Cubic approximation
      int nmb_points = (num_points > 0) ? num_points : 100; //50;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = CubicProject(bspl, data, del, rad, coef, apply_smooth);
    }

  else if (proj_type == 8)
    {
      // Cubic approximation
      int nmb_points = (num_points > 0) ? num_points : 150; //60;
      vector<double> data;
      int nmb_out = 0, nmb = 0;
      double av_dist = 0.0, max_dist = 0.0;
      RDataSet(vol, bspl, nmb_points, rad, data, del, nmb, nmb_out, max_dist, av_dist);
      proj_OK = TriCubicProject(bspl, data, del, rad, coef, apply_smooth);
    }

  while (!proj_OK)
    {
      proj_type--;
      if (proj_type <= 1)
	proj_type = 2;
      //#ifdef DEBUG2
      std::cout << "Failure. Applying projection type " << proj_type << std::endl;
      //#endif
      proj_OK = computeCoef(vol, bspl, proj_type, apply_smooth, dlim, coef, num_points);
    }

  return rad;
}

//==============================================================================
void LRProjection3D::RDataSet(LRSplineVolume *vol, LRBSpline3D *bspl, int nmb_pts,
			    double& rad, vector<double>& data,
			    int del, int& nmb, int& nmb_out,
			    double& max_dist, double& av_dist)
//==============================================================================
{
  rad = -1.0;
  Point uvw0 = bspl->getGrevilleParameter();
  vector<Element3D*> elems;
  vector<LRBSpline3D*> bsplines;
  bsplines.push_back(bspl);
  vector<double>  uvw_dist;
  double umin = uvw0[0], umax = uvw0[0], vmin = uvw0[1], vmax = uvw0[1];
  double wmin = uvw0[2], wmax = uvw0[2];
  double sf_umin = vol->startparam_u();
  double sf_umax = vol->endparam_u();
  double sf_vmin = vol->startparam_v();
  double sf_vmax = vol->endparam_v();
  double sf_wmin = vol->startparam_w();
  double sf_wmax = vol->endparam_w();
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
	      double dist = Utils::distance_squared(uvw0.begin(), uvw0.end(),
						    &elem_data[ka*del]);
	      uvw_dist.push_back(sqrt(dist));
	    }
	  umin = std::min(umin, elems[kj]->umin()); 
	  umax = std::max(umax, elems[kj]->umax()); 
	  vmin = std::min(vmin, elems[kj]->vmin()); 
	  vmax = std::max(vmax, elems[kj]->vmax()); 
	  wmin = std::min(wmin, elems[kj]->wmin()); 
	  wmax = std::max(wmax, elems[kj]->wmax()); 
	}

	  std::sort(uvw_dist.begin(), uvw_dist.end());
	  if ((int)uvw_dist.size() > nmb_pts)
	    {
	      rad = uvw_dist[nmb_pts];
	      if ((umin <= uvw0[0]-rad || umin == sf_umin) &&
		  (umax >= uvw0[0]+rad || umax == sf_umax) &&
		  (vmin <= uvw0[1]-rad || vmin == sf_vmin) &&
		  (vmax >= uvw0[1]+rad || vmax == sf_vmax) &&
		  (wmin <= uvw0[2]-rad || wmin == sf_wmin) &&
		  (wmax >= uvw0[2]+rad || wmax == sf_wmax))
		break;
	    }


	  if (uvw_dist.size() > 0)
	    rad = uvw_dist[uvw_dist.size()-1];
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
bool LRProjection3D::TPproject(LRBSpline3D *bspl,
			     vector<double>& data,
			     int del,
			     double smoothwgt,
			     Point& coef)
//==============================================================================
{
  // Fetch knot vectors of B-spline
  BSplineUniLR* uni_u = bspl->getUnivariate(XDIR);
  BSplineUniLR* uni_v = bspl->getUnivariate(YDIR);
  BSplineUniLR* uni_w = bspl->getUnivariate(ZDIR);
  vector<double> knots1 = uni_u->getKnots();
  vector<double> knots2 = uni_v->getKnots();
  vector<double> knots3 = uni_w->getKnots();

  int startmult1 = uni_u->endmult(true);
  int endmult1 = uni_u->endmult(false);
  int startmult2 = uni_v->endmult(true);
  int endmult2 = uni_v->endmult(false);
  int startmult3 = uni_w->endmult(true);
  int endmult3 = uni_w->endmult(false);

  int order1 = uni_u->degree() + 1;
  int order2 = uni_v->degree() + 1;
  int order3 = uni_w->degree() + 1;

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

  double start3 = knots3[0];
  double end3 = knots3[knots3.size()-1];
  for (int ka=startmult3; ka<order3; ++ka)
    knots3.insert(knots3.begin(), start3);

  for (int ka=endmult3; ka<order3; ++ka)
    knots3.push_back(end3);

  // Extract point and parameter data
  int num_points = (int)data.size()/del;
  int dim = bspl->dimension();
  vector<double> param, points;
  param.reserve(3*num_points);
  points.reserve(dim*num_points);
  for (int ka=0; ka<(int)data.size(); ka+=del)
    {
      for (int kb=0; kb<3; ++kb)
	param.push_back(data[ka+kb]);
      for (int kb=0; kb<dim; ++kb)
	points.push_back(data[ka+kb+3]);
    }

  // Define initial spline surface
  int in1 = (int)knots1.size() - order1;
  int in2 = (int)knots2.size() - order2;
  int in3 = (int)knots3.size() - order3;
  vector<double> init_coefs(in1*in2*in3*dim, 0.0);
  shared_ptr<SplineVolume> init_vol(new SplineVolume(in1, in2, in3,
						      order1, order2, order2,
						      &knots1[0], &knots2[0],
						      &knots3[0],
						      &init_coefs[0], dim));

  // Approximation
  SmoothVolume approx;
  vector<CoefStatus> coef_known(in1*in2*in3, CoefFree);
  approx.attach(init_vol, coef_known);

  double wgt1 = 0.0;
  double wgt3 = (std::min(order1,std::min(order2,order3)) < 4) ? 0.0 :
    0.5*smoothwgt;
  double wgt2 = smoothwgt - wgt3 - wgt1;
  approx.setOptimize(wgt1, wgt2, wgt3);

  double approx_wgt = 1.0 - wgt1 - wgt2 - wgt3;
  vector<double> pnt_wgt(num_points, 1.0);
  approx.setLeastSquares(points, param, pnt_wgt, approx_wgt);

  shared_ptr<SplineVolume> approx_vol;
  try {
  approx.equationSolve(approx_vol);
  }
  catch (...)
    {
      coef = bspl->Coef();
      return false;
    }
  
  // Fetch coefficient corresponding to input B-spline
  int k1 = order1 - startmult1;
  int k2 = order2 - startmult2;
  int k3 = order3 - startmult3;
  int kk = (k3*in2 + k2)*in1 + k1;
  coef = Point(approx_vol->coefs_begin()+kk*dim, approx_vol->coefs_begin()+(kk+1)*dim);

  return true;
}

//==============================================================================
void LRProjection3D::IDWproject(LRBSpline3D *bspl, double rad,
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
      Point curr(&data[ka*del+3], &data[ka*del+3+dim]);
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
bool LRProjection3D::TriLinProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 1, tot_degree = 2;
  return PolynomialProject(degree, degree, degree, tot_degree, bspl, data,
			   del, rad, coef, apply_smooth);
}

//==============================================================================
bool LRProjection3D::LinProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 1, tot_degree = 1;
  return PolynomialProject(degree, degree, degree,  tot_degree, bspl, data, 
			   del, rad, coef, true);
}

//==============================================================================
bool LRProjection3D::TriQuadProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 2, tot_degree = 4;
  return PolynomialProject(degree, degree, degree, tot_degree, bspl, data, 
			   del, rad, coef, apply_smooth);
}

//==============================================================================
bool LRProjection3D::QuadProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 2, tot_degree = 2;
  return PolynomialProject(degree, degree, degree, tot_degree, bspl, data, 
			   del, rad, coef, apply_smooth);
}

//==============================================================================
bool LRProjection3D::TriCubicProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 3, tot_degree = 6;
  return PolynomialProject(degree, degree, degree, tot_degree, bspl, data, 
			   del, rad, coef, apply_smooth);
}

//==============================================================================
bool LRProjection3D::CubicProject(LRBSpline3D *bspl, vector<double>& data, int del,
			       double rad, Point& coef, bool apply_smooth)
//==============================================================================
{
  int degree = 3, tot_degree = 3;
  return PolynomialProject(degree, degree, degree, tot_degree, bspl, data, 
			   del, rad, coef, apply_smooth);
}

void interpolateSample(vector<double>& points, int dim, vector<double>& parvals1,
		       vector<double>& parvals2, vector<double>& parvals3,
		       int order1, int order2, int order3,
		       vector<double>& knots1, vector<double>& knots2,
		       vector<double>& knots3, shared_ptr<SplineVolume>& Bvol)
{
  BsplineBasis basis_u(order1, order1, &knots1[0]);
  BsplineBasis basis_v(order2, order2, &knots2[0]);
  BsplineBasis basis_w(order3, order3, &knots3[0]);
  vector<double> wgt(order1*order2*order3, 1.0);
  Bvol =
    shared_ptr<SplineVolume>(VolumeInterpolator::regularInterpolation(basis_u,
								      basis_v,
								      basis_w,
								      parvals1,
								      parvals2,
								      parvals3,
								      points,
								      dim,
								      false,
								      wgt));
}


void polynomialTerms(double u, double v, double w, int num1, int num2,
		     int num3, int max_num, vector<double>& terms)
{
  vector<double> u_terms(num1+1), v_terms(num2+1), w_terms(num3+1);
  u_terms[0] = v_terms[0] = w_terms[0] = 1.0;

  for (int ki=1; ki<=num1; ++ki)
    u_terms[ki] = u_terms[ki-1]*u;
  for (int ki=1; ki<=num2; ++ki)
     v_terms[ki] = v_terms[ki-1]*v;
  for (int ki=1; ki<=num3; ++ki)
     w_terms[ki] = w_terms[ki-1]*w;

  for (int kh=0, kr=0; kh<=num3; ++kh)
    for (int kj=0; kj<=num2; ++kj)
      for (int ki=0; ki<=num1; ++ki)
	{
	  if (ki+kj+kh <= max_num)
	    terms[kr++] = u_terms[ki]*v_terms[kj]*w_terms[kh];
	}
}

int binom(int n, int r)
{
  if (n >= r && r >= 0)
    {
      int nom = 1;
      for (int ki=2; ki<=n; ++ki)
	nom *= ki;
      int den1 = 1, den2 = 1;
      for (int ki=2; ki<=r; ++ki)
	den1 *= ki;
      for (int ki=2; ki<=n-r; ++ki)
	den2 *= ki;
      return nom/(den1*den2);
    }
  else
    return 0;
}

//==============================================================================
bool LRProjection3D::PolynomialProject(int degree1, int degree2, int degree3,
				       int tot_degree, LRBSpline3D *bspl,
				       vector<double>& data, int del,
				       double rad, Point& coef,
				       bool apply_smooth)
//==============================================================================
{
  bool stat = true;
  double out_fac = 10.0;
  double out_rad = 1.25*rad;
  Point par = bspl->getGrevilleParameter();
  double u, v, w;
  int dim = bspl->dimension();
  int maxdeg = std::max(degree1, std::max(degree2, degree3));
  int num_terms;
  if (tot_degree > maxdeg)
    num_terms = (degree1+1)*(degree2+1)*(degree3+1);
  else if (degree1 == degree2 && degree1 == degree3)
    num_terms = (degree1+1)*(degree1+2)*(degree1+3)/6;
  else
    {
      double d1 = (degree1 == maxdeg) ? maxdeg - degree2 : maxdeg - degree1;
      double d2 = (degree3 == maxdeg) ? maxdeg - degree2 : maxdeg - degree3;
      num_terms = binom(maxdeg+3,3) - binom(maxdeg-d1+2,3) -
	binom(maxdeg-d2+2,3) - binom(maxdeg-d1-d2+1,3);
    }
  vector<double> tmp(num_terms);
  int nmbd = (int)data.size()/del;
  double range[6];
  range[0] = range[2] = range[4] = std::numeric_limits<double>::max();
  range[1] = range[3] = range[5] = std::numeric_limits<double>::lowest();
  
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

      for (int kb=0; kb<3; ++kb)
	range[2*kb] = std::min(range[2*kb], data[ka*del+kb]); 
      for (int kb=0; kb<3; ++kb)
	range[2*kb+1] = std::max(range[2*kb+1], data[ka*del+kb]); 
    }

  // Compute quadratic polynomial approximating the input points with a
  // distance to the B-spline coefficient less than the specified radius,
  // measured in parameter space
  // Define equation system using least squares
  vector<double> A(num_pt*num_terms, 0.0), b(dim*num_pt, 0.0);
  for (int ka=0, kr=0; ka<nmbd; ++ka)
    {
      double dist = sqrt(Utils::distance_squared(par.begin(), par.end(),
						 &data[ka*del]));
      u = data[ka*del];
      v = data[ka*del+1];
      w = data[ka*del+2];

      if (dist > rad)
	continue;

      polynomialTerms(u, v, w, degree1, degree2, degree3, tot_degree, tmp);

      for (int kb=0; kb<num_terms; ++kb)
	A[kr*num_terms+kb] = tmp[kb];
      
      for (int kb=0; kb<dim; ++kb)
	b[kb*num_pt+kr] = data[ka*del+3+kb];
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
      polynomialTerms(par[0], par[1], par[2], degree1, degree2, degree3,
		      tot_degree, tmp);
 
      for (int kb=0; kb<dim; ++kb)
	for (int kc=0; kc<num_terms; ++kc)
	  coef[kb] += x[kb*num_terms+kc]*tmp[kc];
    }
  else
    coef = Polynomial2Coef(x, degree1, degree2, degree3, tot_degree, range, bspl);
      
  //#ifdef DEBUG2
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
      w = data[ka*del+2];
      
      Point pos(dim);
      pos.setValue(0.0);
      polynomialTerms(u, v, w, degree1, degree2, degree3, tot_degree, tmp);
 
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
	  w = data[ka*del+2];
      
	  Point pos(dim);
	  pos.setValue(0.0);
	  polynomialTerms(u, v, w, degree1, degree2, degree3, tot_degree, tmp);
 
	  for (int kb=0; kb<dim; ++kb)
	    for (int kc=0; kc<num_terms; ++kc)
	      pos[kb] += x[kb*num_terms+kc]*tmp[kc];
	  double dd = sqrt(Utils::distance_squared(pos.begin(), pos.end(),
						   &data[ka*del+2]));
	  maxdist2 = std::max(maxdist, dd);
	  avdist2 += dd;
	}
      avdist2 /= (double)num_id;
#ifdef DEBUG2
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
Point
LRProjection3D::Polynomial2Coef(vector<double>& pol, int degree1, int degree2,
				int degree3, int tot_degree, double range[],
				LRBSpline3D *bspl)
//==============================================================================
{
  Point coef;

  int maxdeg = std::max(degree1, std::max(degree2, degree3));
  int num_terms;
  if (tot_degree > maxdeg)
    num_terms = (degree1+1)*(degree2+1)*(degree3+1);
  else if (degree1 == degree2 && degree1 == degree3)
    num_terms = (degree1+1)*(degree1+2)*(degree1+3)/6;
  else
    {
      double d1 = (degree1 == maxdeg) ? maxdeg - degree2 : maxdeg - degree1;
      double d2 = (degree3 == maxdeg) ? maxdeg - degree2 : maxdeg - degree3;
      num_terms = binom(maxdeg+3,3) - binom(maxdeg-d1+2,3) -
	binom(maxdeg-d2+2,3) - binom(maxdeg-d1-d2+1,3);
    }
  vector<double> tmp(num_terms);

  // Represent the polynomial by a Bezier surface
  // Identify inner knots
  double u, v, w;
  int dim = bspl->dimension();
  Direction3D dir[3];
  dir[0] = XDIR;
  dir[1] = YDIR;
  dir[2] = ZDIR;
  int order[3];
  int startmult[3];
  double Brange[6], pdel[3];
  vector<vector<double> > knots(3);
  vector<vector<double> > knots2(3);
  for (int ka=0; ka<3; ++ka)
    {
      BSplineUniLR* uni = bspl->getUnivariate(dir[ka]);
      startmult[ka] = bspl->endmult(dir[ka], true);
      int endmult = bspl->endmult(dir[ka], false);
      order[ka] = bspl->degree(dir[ka]) + 1;

      vector<double> kval = uni->getKnots();
      int kn = order[ka] - startmult[ka] - endmult + 1;
      knots[ka].resize(kn);
      for (int kb=0; kb<kn; ++kb)
	{
	  knots[ka][kb] = kval[kb+startmult[ka]];
	}

      // Define spline space for B-spline surface interpolating the polynomial
      Brange[2*ka] = uni->min();
      Brange[2*ka+1] = uni->max();

      knots2[ka].resize(2*order[ka]);
      for (int kb=0; kb<order[ka]; ++kb)
	{
	  knots2[ka][kb] = Brange[2*ka];
	  knots2[ka][order[ka]+kb] = Brange[2*ka+1];
	}

      Brange[2*ka] = std::max(Brange[2*ka],range[2*ka]);
      Brange[2*ka+1] = std::max(Brange[2*ka],range[2*ka+1]);
      pdel[ka] = (Brange[2*ka+1] - Brange[2*ka])/(double)(order[ka]-1);
   }


  // Sample polynomial
  vector<double> points(order[0]*order[1]*order[2]*dim, 0.0);
  vector<double> parvals1(order[0]), parvals2(order[1]), parvals3(order[2]);
  int ka, kb, kc, kd, kf;
  for (kc=0, w=Brange[4]; kc<order[2]; ++kc, w+=pdel[2])
    {
      parvals3[kc] = w;
      for (kb=0, v=Brange[2]; kb<order[1]; ++kb, v+=pdel[1])
	{
	  parvals2[kb] = v;
	  for (ka=0, u=Brange[0]; ka<order[0]; ++ka, u+=pdel[0])
	    {
	      parvals1[ka] = u;
	      polynomialTerms(u, v, w, degree1, degree2, degree3,
			      tot_degree, tmp);
 
	      for (kd=0; kd<dim; ++kd)
		for (kf=0; kf<num_terms; ++kf)
		  points[dim*((kc*order[2]+kb)*order[0]+ka)+kd] +=
		    pol[kd*num_terms+kf]*tmp[kf];
	    }
	}
    }

  // Interpolate
  shared_ptr<SplineVolume> Bvol;
  interpolateSample(points, dim, parvals1, parvals2,
		    parvals3, order[0],  order[1], order[2],
		    knots2[0], knots2[1], knots2[2], Bvol);

  // Insert inner knots
  for (int ka=0; ka<3; ++ka)
    if (knots[ka].size() > 0)
	Bvol->insertKnot(ka, knots[ka]);

  // Fetch coefficient
  int k1 = order[0] - startmult[0];
  int k2 = order[1] - startmult[1];
  int k3 = order[2] - startmult[2];
  int kk = (k3*Bvol->numCoefs(1)+k2)*Bvol->numCoefs(0) + k1;
  coef = Point(Bvol->ctrl_begin()+kk*dim, Bvol->ctrl_begin()+(kk+1)*dim);

#ifdef DEBUG2
  std::cout << "Coef: " << coef << std::endl;
#endif
  
  return coef;
}

