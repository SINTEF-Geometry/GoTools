//===========================================================================
//                                                                           
// File: TrimCurve
//                                                                           
// Created: Nov 4, 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/creators/TrimCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
TrimCurve::TrimCurve(CurveOnSurface* sfcv)
  : sfcv_(sfcv)
//===========================================================================
{
  start_ = sfcv->startparam();
  end_ = sfcv->endparam();
  startpt_ = sfcv->ParamCurve::point(start_);
  endpt_ = sfcv->ParamCurve::point(end_);
}

//===========================================================================
TrimCurve::TrimCurve(CurveOnSurface* sfcv, double start, double end)
  : sfcv_(sfcv), start_(start), end_(end)
//===========================================================================
{
  startpt_ = sfcv->ParamCurve::point(start);
  endpt_ = sfcv->ParamCurve::point(end);
}

//===========================================================================
TrimCurve::TrimCurve(Point startpt, Point endpt, CurveOnSurface* sfcv)
  : sfcv_(sfcv), startpt_(startpt), endpt_(endpt)
//===========================================================================
{
  double clo_d;
  Point clo_pt;
  sfcv->closestPoint(startpt, sfcv->startparam(), sfcv->endparam(),
		     start_, clo_pt, clo_d);
  sfcv->closestPoint(endpt, sfcv->startparam(), sfcv->endparam(),
		     end_, clo_pt, clo_d);
}
//===========================================================================
TrimCurve::~TrimCurve()
//===========================================================================
{
}
 
//===========================================================================
vector<Point> TrimCurve::eval(double t)
//===========================================================================
{
  vector<Point> result;
  evaluate(t, 0, result);
  return result;
}

//===========================================================================
void TrimCurve::eval(double t, int n, vector<vector<Point> >& der)
//===========================================================================
{
  vector<Point> result;

  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;
  evaluate(t, n, result);

  der.resize(2);
  for (int ki=0; ki<2; ++ki)
    {
      vector<Point> tmp;
      tmp.insert(tmp.begin(), result.begin()+(n+1)*ki, 
		 result.begin()+(n+1)*(ki+1));
      der[ki] = tmp;
    }
}

//===========================================================================
double TrimCurve::start()
//===========================================================================
{
  return start_;
}

//===========================================================================
double TrimCurve::end()
//===========================================================================
{
  return end_;
}

//===========================================================================
int TrimCurve::dim()
//===========================================================================
{
// One geometry curve and one parameter curve
  return sfcv_->dimension() + 2;  
}

//===========================================================================
bool TrimCurve::approximationOK(double par, 
				const vector<Point>& approxpos,
				double tol1, double tol2)
//===========================================================================
{
  // Check geometry curve
  vector<Point> pos = eval(par);
  if (approxpos[0].dist(pos[0]) > tol1)
    return false;

  // Check the parameter curve in geometry space
  Point pp1 = sfcv_->underlyingSurface()->point(pos[1][0], pos[1][1]);
  Point pp2 = sfcv_->underlyingSurface()->point(approxpos[1][0], 
						approxpos[1][1]);
  if (pp1.dist(pp2) > tol2)
    return false;

  return true;  // The accuracy is good enough
}

//===========================================================================
int TrimCurve::nmbCvs()
//===========================================================================
{
  return 2;
}

//===========================================================================
void TrimCurve::evaluate(double t, int n, vector<Point>& result)
//===========================================================================
{
  // double ang_tol = 0.001;
  double eps = 1.0e-6;
  double num_tol = 1.0e-15;

  // NB! At most 1. derivitive of curves are computed
  if (n > 1)
    n = 1;
  result.resize(2*(n+1));

  // Evaluate first surface curve
  vector<Point> der1(2);
  sfcv_->point(der1, t, 1);

  // Special treatment of the endpoints
  if (fabs(t - sfcv_->startparam()) < num_tol)
    der1[0] = startpt_;
  else if (fabs(t - sfcv_->endparam()) < num_tol)
    der1[0] = endpt_;

  // Project the point onto the surface
  // Compute the closest point to the point on the first curve onto the surface
  // corresponding to the second curve
  double clo_u, clo_v, clo_dist;
  Point clo_pt;
  Point sf_par = sfcv_->faceParameter(t);
  sfcv_->underlyingSurface()->closestPoint(der1[0], clo_u, clo_v, clo_pt,
					   clo_dist, eps, NULL, sf_par.begin());
  
  int write_srf = false;
  if (write_srf)
    {
      std::ofstream out("TrimCurve_sf.g2");
      sfcv_->underlyingSurface()->writeStandardHeader(out);
      sfcv_->underlyingSurface()->write(out);
    }

  // Special treatment of the endpoints
  if (fabs(t - sfcv_->startparam()) < num_tol)
    clo_pt = startpt_;
  else if (fabs(t - sfcv_->endparam()) < num_tol)
    clo_pt = endpt_;

  result[0] = clo_pt;
  result[n+1] = Point(clo_u, clo_v);
  if (n == 1)
    {
      // Evaluate second surface with 1. derivatives
      vector<Point> der2 = sfcv_->underlyingSurface()->point(clo_u, clo_v, 1);
	
      // Project average of input tangent onto the second surface
      double u, v;
      double len = der1[1].length();
      der1[1].normalize();
      CoonsPatchGen::blendcoef(&der2[1][0], &der2[2][0],
			       &der1[1][0], 3, 1, &u, &v);
      Point proj1 = u*der2[1] + v*der2[2];
	    
      proj1.normalize();
      result[1] = proj1;
      result[1] *= len;

      // The space tangent in the parameter domain. Project again
      CoonsPatchGen::blendcoef(&der2[1][0], &der2[2][0],
			       &result[1][0], 3, 1, &u, &v);
      result[3] = Point(u, v);
    }
}

