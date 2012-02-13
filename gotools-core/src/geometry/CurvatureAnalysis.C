//===========================================================================
//                                                                           
// File: CurvatureAnalysis.C                                                 
//                                                                           
// Created: Wed Mar 16 14:34:25 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurvatureAnalysis.C,v 1.10 2009-03-20 07:53:32 vsk Exp $
//                                                                           
//===========================================================================


#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/GeometryTools.h"


using std::min;
using std::max;
using std::vector;

namespace Go
{
//===========================================================================
void CurvatureAnalysis::computeFirstFundamentalForm(const ParamSurface& sf,
				 double u, double v, int derivs,
				 std::vector<double>& form)
//===========================================================================
{
    int sz = (derivs + 1)*(derivs + 2)/2;
    form.resize(3*sz);
    std::vector<Point> pts((derivs + 2)*(derivs + 3)/2);
    // To evaluate k derivatives of E, F, G we need k+1 derivatives
    // of the surface.
    sf.point(pts, u, v, derivs + 1);
    // Compute the first few terms
    form[0] = pts[1]*pts[1];
    form[1] = pts[1]*pts[2];
    form[2] = pts[2]*pts[2];
    if (derivs > 0) {
	form[3] = 2.0*pts[1]*pts[3];
	form[4] = pts[1]*pts[4] + pts[2]*pts[3];
	form[5] = 2.0*pts[2]*pts[4];
	form[6] = 2.0*pts[1]*pts[4];
	form[7] = pts[1]*pts[5] + pts[2]*pts[4];
	form[8] = 2.0*pts[2]*pts[5];
    }
    // Using the 2-variable Leibniz rule,
    // we could eventually compute the higher derivatives.
    for (int i = 2; i <= derivs; ++i) {
	THROW("Higher derivatives of E, F, G not implemented.");
    }
}


//===========================================================================
void CurvatureAnalysis::computeSecondFundamentalForm(const ParamSurface& sf,
				      double u, double v,
				      double form1[3],
				      double form2[3])
//===========================================================================
{
    // To evaluate both fundamental forms, we need the second
    // derivatives and the normal.
    std::vector<Point> pts(6);
    Point normal;
    sf.point(pts, u, v, 2);
    sf.normal(normal, u, v);
    // Compute the first few terms
    form1[0] = pts[1]*pts[1];
    form1[1] = pts[1]*pts[2];
    form1[2] = pts[2]*pts[2];
    form2[0] = normal*pts[3];
    form2[1] = normal*pts[4];
    form2[2] = normal*pts[5];
}


//===========================================================================
void CurvatureAnalysis::curvatures(const ParamSurface& sf,
		double u, double v,
		double& K, double& H)
//===========================================================================
{
    double I[3];
    double II[3];
    computeSecondFundamentalForm(sf, u, v, I, II);
    double denom = I[0]*I[2]-I[1]*I[1];
    K = (II[0]*II[2]-II[1]*II[1])/denom;
    H = (II[0]*I[2]-2*II[1]*I[1]+II[2]*I[0])/(2*denom);
}

//===========================================================================
void CurvatureAnalysis::principalCurvatures(const ParamSurface& sf,
			 double u, double v,
			 double& k1, Point& d1,  // Direction given in par. domain
			 double& k2, Point& d2)
//===========================================================================
{
    double resolution = 1.0e-15;

    // Compute surface derivatives
    int derivs = 1;
    std::vector<Point> pts(3);
    sf.point(pts, u, v, derivs);
    Point Su = pts[1];
    Point Sv = pts[2];

    // Compute 1. and 2. fundamental form
    double I[3];
    double II[3];
    computeSecondFundamentalForm(sf, u, v, I, II);
    double denom = I[0]*I[2]-I[1]*I[1];

    // Calculate the transformation matrix
    double transform[4]; // The sequence is a11, a12, a21, a22.        
    transform[0] = (I[1]*II[1] - I[2]*II[0])/denom;
    transform[1] = (I[1]*II[2] - I[2]*II[1])/denom;
    transform[2] = (I[1]*II[0] - I[0]*II[1])/denom;
    transform[3] = (I[1]*II[1] - I[0]*II[2])/denom;
   
    // Initialize the principal directions to the parameter directions
    d1.setValue(1.0, 0.0);
    d2.setValue(0.0, 1.0);

    // Calculate the principal curvatures
    double a = 1.;
    double b = transform[0] + transform[3];
    double c = transform[0]*transform[3] - transform[1]*transform[2];
      
    double sqrt_arg = b*b - 4.*a*c;
    if (sqrt_arg < resolution)
      {
	 k1 = - b/(2.*a);
	 k2 = k1;
	 return;
      }

    k1 = (- b + sqrt(sqrt_arg))/(2.*a);
    k2 = (- b - sqrt(sqrt_arg))/(2.*a);
      
    // Calculate the curvature directions.
    double ratio, length; 

    // Direction corresponding to the maximal curvature
    if (fabs(transform[0] + k1) < resolution &&
	     fabs(transform[1]) < resolution)
    {
	// Parallel to u direction
	length = 1./sqrt(Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2]);
	d1.setValue(length, 0.0);
    }
    else if (fabs(transform[3] + k1) < resolution &&
		  fabs(transform[2]) < resolution)
    {
	// Parallel to v direction
	length = 1./sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]);
	d1.setValue(0.0, length);
    }
    else if (fabs(transform[0] + k1) < fabs(transform[1]))
    {
	ratio = (transform[0] + k1)/transform[1];
	length = 1./sqrt((Su[0] - ratio*Sv[0])*(Su[0] - ratio*Sv[0]) +
			 (Su[1] - ratio*Sv[1])*(Su[1] - ratio*Sv[1]) +
			 (Su[2] - ratio*Sv[2])*(Su[2] - ratio*Sv[2]));
	d1.setValue(length, -ratio*length);
    }
    else
    {
	ratio = transform[1]/(transform[0] + k1);
	length = 1./sqrt((Sv[0] - ratio*Su[0])*(Sv[0] - ratio*Su[0]) +
			 (Sv[1] - ratio*Su[1])*(Sv[1] - ratio*Su[1]) +
			 (Sv[2] - ratio*Su[2])*(Sv[2] - ratio*Su[2]));
	d1.setValue(-ratio*length, length);
    }

    // Minimal curvature direction
    if (fabs(transform[0] + k2) < resolution &&
	     fabs(transform[1]) < resolution)
    {
	// Parallel to u direction
	length = 1./sqrt(Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2]);
	d2.setValue(length, 0.0);
    }
    else if (fabs(transform[3] + k2) < resolution &&
		  fabs(transform[2]) < resolution)
    {
	// Parallel to v direction
	length = 1./sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]);
	d2.setValue(0.0, length);
    }
    else if (fabs(transform[0] + k2) < fabs(transform[1]))
    {
	ratio = (transform[0] + k2)/transform[1];
	length = 1./sqrt((Su[0] - ratio*Sv[0])*(Su[0] - ratio*Sv[0]) +
			 (Su[1] - ratio*Sv[1])*(Su[1] - ratio*Sv[1]) +
			 (Su[2] - ratio*Sv[2])*(Su[2] - ratio*Sv[2]));
	d2.setValue(length, -ratio*length);
    }
    else
    {
	ratio = transform[1]/(transform[0] + k2);
	length = 1./sqrt((Sv[0] - ratio*Su[0])*(Sv[0] - ratio*Su[0]) +
			 (Sv[1] - ratio*Su[1])*(Sv[1] - ratio*Su[1]) +
			 (Sv[2] - ratio*Su[2])*(Sv[2] - ratio*Su[2]));
	d2.setValue(-ratio*length, length);
    }

	 

}


//===========================================================================
void CurvatureAnalysis::minimalCurvatureRadius(const ParamSurface& sf,
			    double tolerance,
			    double& mincurv,
			    double& pos_u,
			    double& pos_v,
			    double degenerate_eps,
			    double curv_tol)
//===========================================================================
{
    const RectDomain& dom = sf.containingDomain();
    double start_u = dom.umin();
    double end_u = dom.umax();
    double start_v = dom.vmin();
    double end_v = dom.vmax();

    // Modify domain at degenerate boundaries
    bool b, r, t, l;
    bool is_degen = sf.isDegenerate(b, r, t, l, degenerate_eps);
    if (is_degen)
    {
	double par_eps = 1.0e-4;  // @@@ VSK, dec 08. Rather arbitrary
	double fac = 1.0e-3;
	if (b)
	    start_v += std::min(par_eps, fac*(end_v-start_v));
	if (r)
	    end_u -= std::min(par_eps, fac*(end_u-start_u));
	if (t)
	    end_v -= std::min(par_eps, fac*(end_v-start_v));
	if (l)
	    start_u += std::min(par_eps, fac*(end_u-start_u));
    }

  vector<double> param_u;
  vector<double> param_v;
  vector<vector<double> > curvs;

  CurvatureAnalysis::evaluateMinCurvatureRadius(sf,
			     start_u, end_u, start_v, end_v,
			     tolerance,
			     param_u, param_v, curvs,
			     mincurv, pos_u, pos_v,
			     true);

  /*
  // Make an interpolated surface interpolated_surf over the curvatures

  // Intersect by level plane defined by mincurv
  vector<shared_ptr<SplineCurve> > level_curve =
    intersect( aPlane, interpolated_surf);

  // Group curves into a set level_curves in two dimensions, answering to "boundingbox"

  for (int i = 0; i < level_curve.size() ++i)
    {
      BoundingBox b = level_curve[i].boundingBox();
      double box_start_u = max(start_u, b.low()[0] - tol);
      double box_end_u = min(end_u, b.high()[0] - tol);
      double box_start_v = max(start_v, b.low()[1] - tol);
      double box_end_v = min(end_v, b.high()[1] - tol);
      
      param_u.resize(0);
      param_v.resize(0);
      curvs.resize(0);

      local_tol = tolerance * min ((box_end_u - box_start_u) / (start_u-end_u),
				   (box_end_v - box_start_v) / (start_v-end_v));
      local_tol = max(local_tol, 0.1*tolerance);

      evaluateMinCurvatureRadius(sf,
				 box_start_u, box_end_u, box_start_v, box_end_v,
				 local_tol,
				 param_u, param_v, curvs,
				 mincurv, pos_u, pos_v,
				 false);
    }
  */

  // Rerun and make a finer evaluation round local minima
  vector<double> dummy_u;
  vector<double> dummy_v;
  vector<vector<double> > dummy_curvs;

  //double local_tol = 2.0 * tolerance * tolerance / max(end_u-start_u, end_v-start_v);
  double local_tol = tolerance/(double)(std::max(param_u.size(), param_v.size()));
  local_tol = min(local_tol, tolerance);
  double step_u = (end_u-start_u) / ((int)param_u.size() - 1);
  double step_v = (end_v-start_v) / ((int)param_v.size() - 1);

  for (int i = 0; i < (int)param_u.size(); ++i)
    {
      int i_prev = i - (i>0);
      int i_next = i + (i < (int)param_u.size()-1);
      for (int j = 0; j < (int)param_v.size(); ++j)
	{
	  int j_prev = j - (j>0);
	  int j_next = j + (j < (int)param_v.size()-1);
	  double local_curv = curvs[i][j];
	  if (local_curv > curvs[i_next][j] ||
	      local_curv > curvs[i_prev][j] ||
	      local_curv > curvs[i][j_next] ||
	      local_curv > curvs[i][j_prev] ||
	      local_curv > curvs[i_next][j_next] ||
	      local_curv > curvs[i_next][j_prev] ||
	      local_curv > curvs[i_prev][j_next] ||
	      local_curv > curvs[i_prev][j_prev])
	    continue;
	  double limit = 2.0 * local_curv - mincurv + curv_tol*mincurv;
	  bool rand_i = i==i_next || i == i_prev;
	  bool rand_j = j==j_next || j == j_prev;
	  if (!rand_i && curvs[i_prev][j]<limit && curvs[i_next][j]<limit) continue;
	  if (!rand_j && curvs[i][j_prev]<limit && curvs[i][j_next]<limit) continue;
	  if (!rand_i && !rand_j && curvs[i_prev][j_prev]<limit && curvs[i_next][j_next]<limit) continue;
	  if (!rand_i && !rand_j && curvs[i_prev][j_next]<limit && curvs[i_next][j_prev]<limit) continue;

	  dummy_u.resize(0);
	  dummy_v.resize(0);
	  dummy_curvs.resize(0);
	  evaluateMinCurvatureRadius(sf,
				     start_u+double(i_prev)*step_u,
				     start_u+double(i_next)*step_u,
				     start_v+double(j_prev)*step_v,
				     start_v+double(j_next)*step_v,
				     local_tol,
				     dummy_u, dummy_v, dummy_curvs,
				     mincurv, pos_u, pos_v,
				     false);
	}
    }

}


//===========================================================================
void CurvatureAnalysis::evaluateMinCurvatureRadius(const ParamSurface& sf,
				double start_u, double end_u, double start_v, double end_v,
				double tolerance,
				vector<double>& param_u, vector<double>& param_v,
				vector<vector<double> >& curvs,
				double& mincurv,
				double& minpos_u,
				double& minpos_v,
				bool initialize)
//===========================================================================
{
  int minCells = 3;
  int maxCells = 20;

  double area[4];
  area[0] = start_u;
  area[1] = end_u;
  area[2] = start_v;
  area[3] = end_v;

  double len_u, len_v;
  GeometryTools::estimateSurfaceSize(sf, len_u, len_v, area);

  double huge_rad = MAXDOUBLE;
  double tol2 = 10.0*tolerance;
  int pts_u = int (len_u / tol2) + 1;
  int pts_v = int (len_v / tol2) + 1;

  double iso_trim_tol = 0.0001*std::min(end_u - start_u, end_v - start_v);
  iso_trim_tol = std::max(iso_trim_tol, 1.0e-7);

  bool iso_trimmed = sf.isIsoTrimmed(iso_trim_tol);

  pts_u = min(maxCells, max(minCells, pts_u));
  pts_v = min(maxCells, max(minCells, pts_v));

  param_u.resize(pts_u);
  param_v.resize(pts_v);
  curvs.resize(pts_u);

  double step_u = (end_u - start_u)/double(pts_u-1);
  double step_v = (end_v - start_v)/double(pts_v-1);

  double pos_u = start_u;
  for (int i = 0; i < pts_u; pos_u += step_u, ++i)
    {
      curvs[i].resize(pts_v);
      param_u[i] = pos_u;

      double pos_v = start_v;
      for (int j = 0; j < pts_v; pos_v += step_v, ++j)
	{
	  if (i==0) param_v[j] = pos_v;
	  double curveRad;
	  if (iso_trimmed || sf.inDomain(pos_u, pos_v))
	  {
// 	      double gauss_curv, mean_curv;
// 	      curvatures(sf, pos_u, pos_v, gauss_curv, mean_curv);
// 	      if (gauss_curv < 0.0) gauss_curv = -gauss_curv;
// 	      if (gauss_curv < 1e-12) gauss_curv = 1e-12;
//	      curveRad = 1.0 / gauss_curv;

	      // Evaluate principal curvatures
	      Point d1, d2;
	      double k1, k2;
	      CurvatureAnalysis::principalCurvatures(sf, pos_u, pos_v, k1, d1, k2, d2);
	      double kmax = std::max(fabs(k1), fabs(k2));
	      curveRad = (kmax > 1.0e-12) ? 1.0 / kmax : MAXDOUBLE;
	  }
	  else 
	      curveRad = huge_rad;

	  curvs[i][j] = curveRad;
	  if ( (i==0 && j==0 && initialize) || curveRad < mincurv)
	    {
	      mincurv = curveRad;
	      minpos_u = pos_u;
	      minpos_v = pos_v;
	    }
	}

    }

}




} // namespace Go
