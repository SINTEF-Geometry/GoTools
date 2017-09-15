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

#include "GoTools/creators/IntCrvEvaluator.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/ClosestPoint.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;

//===========================================================================
IntCrvEvaluator::IntCrvEvaluator(shared_ptr<CurveOnSurface> sfcv1,
				 double start1, double end1,
				 shared_ptr<CurveOnSurface> sfcv2,
				 double start2, double end2,
				 bool same_orientation,
				 int keep_crv)
  : sfcv1_(sfcv1), start1_(start1), end1_(end1), sfcv2_(sfcv2), 
    start2_(start2), end2_(end2), same_orientation_(same_orientation),
    keep_crv_(keep_crv), max_dist_(0.0), max_err_(0.0)
//===========================================================================
{
}

//===========================================================================
IntCrvEvaluator::~IntCrvEvaluator()
//===========================================================================
{
}
 
//===========================================================================
vector<Point> IntCrvEvaluator::eval(double t)
//===========================================================================
{
  vector<Point> result;
  evaluate(t, 0, result);
  return result;
}

//===========================================================================
void IntCrvEvaluator::eval(double t, int n, vector<vector<Point> >& der)
//===========================================================================
{
  vector<Point> result;

  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;
  evaluate(t, n, result);

  der.resize(3);
  for (int ki=0; ki<3; ++ki)
    {
      vector<Point> tmp;
      tmp.insert(tmp.begin(), result.begin()+(n+1)*ki, 
		 result.begin()+(n+1)*(ki+1));
      der[ki] = tmp;
    }
}

//===========================================================================
double IntCrvEvaluator::start()
//===========================================================================
{
  return start1_;
}

//===========================================================================
double IntCrvEvaluator::end()
//===========================================================================
{
  return end1_;
}

//===========================================================================
int IntCrvEvaluator::dim()
//===========================================================================
{
// One geometry curve and two parameter curves
  return sfcv1_->dimension() + 2 + 2;  
}

//===========================================================================
bool IntCrvEvaluator::approximationOK(double par, 
				      const vector<Point>& approxpos,
				      double tol1, double tol2)
//===========================================================================
{
  // Check whether the approximation can ever succeed
  if (max_dist_ > tol1)
    THROW("Too large distance between surfaces");

  // Check geometry curve
  vector<Point> pos = eval(par);
  double dist = approxpos[0].dist(pos[0]);
  max_err_ = std::max(dist, max_err_);

  if (dist > tol1)
    return false;

  // Check parameter curve 1 in geometry space
  Point pp1 = sfcv1_->underlyingSurface()->point(pos[1][0], pos[1][1]);
  Point pp2 = sfcv1_->underlyingSurface()->point(approxpos[1][0], 
						 approxpos[1][1]);
  if (pp1.dist(pp2) > tol1)
    return false;

  // Check parameter curve 2 in geometry space
  pp1 = sfcv2_->underlyingSurface()->point(pos[2][0], pos[2][1]);
  pp2 = sfcv2_->underlyingSurface()->point(approxpos[2][0], 
					   approxpos[2][1]);
  if (pp1.dist(pp2) > tol1)
    return false;

  return true;  // The accuracy is good enough
}

//===========================================================================
int IntCrvEvaluator::nmbCvs()
//===========================================================================
{
  return 3;
}

//===========================================================================
void IntCrvEvaluator::evaluate(double t, int n, vector<Point>& result)
//===========================================================================
{
  double eps = 1.0e-6;
  double ang_tol = 0.05;

  // Evaluation is performed by iterating to the intersection point between
  // the two surfaces defined by the surface curve and a plane being
  // orhtogonal to the intersection curve in the given parameter

  // NB! At most 1. derivitive of curves are computed
  if (n > 1)
    n = 1;
  result.resize(3*(n+1));

  // Evaluate first surface curve
  vector<Point> der1(2);
  sfcv1_->point(der1, t, 1);

  // Find the corresponding point in the other surface curve
  double ts1 = sfcv1_->startparam();
  double te1 = sfcv1_->endparam();
  double ts2 = sfcv2_->startparam();
  double te2 = sfcv2_->endparam();
  double rel = (te2 - ts2)/(te1 - ts1);
  double guess = same_orientation_ ? ts2 + (t-ts1)*rel : te2 - (t-ts1)*rel;
  guess = std::max(ts2, std::min(te2, guess)); // Numerics
  
  double par, clo_dist;
  Point clo_pt;
  sfcv2_->closestPoint(der1[0], ts2, te2, par, clo_pt, clo_dist, &guess);

  // Evaluate second surface curve
  vector<Point> der2(2);
  sfcv2_->point(der2, par, 1);

  // Define plane
  int sgn = same_orientation_ ? 1.0 : -1.0;
   vector<Point> plane_def(2);
   plane_def[0] = 0.5*(der1[0] + der2[0]);
   double len = der1[1].length();
   der1[1].normalize();
   der2[1].normalize();
   plane_def[1] = 0.5*(der1[1] + sgn*der2[1]);

   // Evaluate the first surface
    vector<Point> input_point_1(7, Point(3)); // only entry [0], [1], [2] and [6] will be used
    Point sf_par1 = sfcv1_->faceParameter(t);
    shared_ptr<ParamSurface> sf1 = sfcv1_->underlyingSurface();
    sf1->point(input_point_1, sf_par1[0], sf_par1[1], 1);
    input_point_1[6] = input_point_1[1].cross(input_point_1[2]);
				      
   // Evaluate the second surface
    vector<Point> input_point_2(7, Point(3)); // only entry [1], [1], [2] and [6] will be used
    Point sf_par2 = sfcv2_->faceParameter(par);
    shared_ptr<ParamSurface> sf2 = sfcv2_->underlyingSurface();
    sf2->point(input_point_2, sf_par2[0], sf_par2[1], 1);
    input_point_2[6] = input_point_2[1].cross(input_point_2[2]);

    // Iterate to intersection point
    vector<Point> result_pt_1(7, Point(3));
    vector<Point> result_pt_2(7, Point(3));
    Point surface_1_param, surface_2_param;

    // Distinguish between transversal and tangential points. The latter can move
    // the intersection point too much if we try to find the closest point
    double angle = input_point_1[6].angle(input_point_2[6]);
    angle = std::min(angle, fabs(M_PI-angle));
    if (angle < ang_tol)
      {
	// Project the medium point onto the two surfaces
	Point med = 0.5*(input_point_1[0] + input_point_2[0]);
	double u1, u2, v1, v2, d1, d2;
	Point pt1, pt2;
	sf1->closestPoint(med, u1, v1, pt1, d1, eps, NULL, sf_par1.begin());
	sf1->point(result_pt_1, u1, v1, 2);
	result_pt_1[6] = result_pt_1[1].cross(result_pt_1[2]);
	surface_1_param = Point(u1, v1);
	
	sf2->closestPoint(med, u2, v2, pt2, d2, eps, NULL, sf_par2.begin());
	sf2->point(result_pt_2, u2, v2, 2);
	result_pt_2[6] = result_pt_2[1].cross(result_pt_2[2]);
	surface_2_param = Point(u2, v2);
      }
    else
      {
	int stat = 0;
	double num_tol = 1.0e-15;
	AlgorithmChoice algo = GEOMETRICAL; // FUNCTIONAL; //GEOMETRICAL;
	double dist1 = input_point_1[0].dist(input_point_2[0]);
	ClosestPoint::closestPtSurfSurfPlane(plane_def, 
			       input_point_1, 
			       input_point_2,
			       sf_par1,
			       sf_par2,
			       sfcv1_->underlyingSurface().get(),
			       sfcv2_->underlyingSurface().get(),
			       sqrt(num_tol), // precision in minimization routine cannot
			       // be better than root of computational precision
			       result_pt_1,
			       result_pt_2,
			       surface_1_param,
			       surface_2_param,
			       stat,
			       algo);

	double dist2 = result_pt_1[0].dist(result_pt_2[0]);
	if (dist2 > dist1)
	  {
	    algo = FUNCTIONAL;
	    ClosestPoint::closestPtSurfSurfPlane(plane_def, 
				   input_point_1, 
				   input_point_2,
				   sf_par1,
				   sf_par2,
				   sfcv1_->underlyingSurface().get(),
				   sfcv2_->underlyingSurface().get(),
				   sqrt(num_tol), 
				   result_pt_1,
				   result_pt_2,
				   surface_1_param,
				   surface_2_param,
				   stat,
				   algo);
	  }
      }
    max_dist_ = std::max(max_dist_, result_pt_1[0].dist(result_pt_2[0]));

    if (keep_crv_ == 1)
      result[0] = result_pt_1[0];
    else if (keep_crv_ == 2)
      result[0] = result_pt_2[0];
    else
      result[0] = 0.5*(result_pt_1[0]+result_pt_2[0]);
    result[n+1] = surface_1_param;
    result[2*n+2] = surface_2_param;
    if (n == 1)
      {
	double ang =  result_pt_1[6].angle(result_pt_2[6]);
	if (ang > ang_tol && ang < M_PI-ang_tol)
	  {
	    // The direction of the intersection curve is defined by the
	    // cross product of the surface normal
	    result[1] = result_pt_1[6].cross(result_pt_2[6]);
	    result[1].normalize();
	    if (result[1]*plane_def[1] < 0.0)
	      result[1] *= -1.0;
	    result[1] *= len;
	  }
	else
	  {
	    // Tangential situation. Project average of input tangent
	    // onto the two surfaces
	    double u, v;
	    CoonsPatchGen::blendcoef(&result_pt_1[1][0], &result_pt_1[2][0],
				     &plane_def[1][0], 3, 1, &u, &v);
	    Point proj1 = u*result_pt_1[1] + v*result_pt_1[2];
	    
	    CoonsPatchGen::blendcoef(&result_pt_2[1][0], &result_pt_2[2][0],
				     &plane_def[1][0], 3, 1, &u, &v);
	    Point proj2 = u*result_pt_2[1] + v*result_pt_2[2];
	    proj1.normalize();
	    proj2.normalize();
	    result[1] = 0.5*(proj1 + proj2);
	    result[1] *= len;
	  }

	// The length of the tangent of the intersection curve?

	// Project space tangent into the parameter domain
	double s, t;
	CoonsPatchGen::blendcoef(&result_pt_1[1][0], &result_pt_1[2][0],
				 &result[1][0], 3, 1, &s, &t);
	result[3] = Point(s, t);

	CoonsPatchGen::blendcoef(&result_pt_2[1][0], &result_pt_2[2][0],
				 &result[1][0], 3, 1, &s, &t);
	result[5] = Point(s, t);
      }
}

