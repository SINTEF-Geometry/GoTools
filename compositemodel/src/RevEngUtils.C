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

#include "GoTools/compositemodel/RevEngUtils.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/utils/LUDecomp.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "sislP.h"
#include "newmat.h"
#include "newmatap.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;

//#define DEBUG
//#define DEBUG_BLEND
//#define DEBUG_CONE
//#define DEBUG_APPROX

typedef MatrixXD<double, 3> Matrix3D;

//===========================================================================
void RevEngUtils::principalAnalysis(std::vector<RevEngPoint*>& points, 
				    double lambda[3], double eigenvec[3][3])
//===========================================================================
{
  if (points.size() < 5)
    return;
  vector<Point> remaining(points.size()-1);
  Vector3D xyz = points[0]->getPoint();
  Point curr(xyz[0], xyz[1], xyz[2]);
  for (size_t ki=1; ki<points.size(); ++ki)
    {
      xyz = points[ki]->getPoint();
      Point curr2(xyz[0], xyz[1], xyz[2]);
      remaining[ki-1] = curr2;
    }
  principalAnalysis(curr, remaining, lambda, eigenvec);
}

//===========================================================================
void RevEngUtils::principalAnalysis(Point& curr, vector<Point>& points, 
				    double lambda[3], double eigenvec[3][3])
//===========================================================================
{
  // Compute covariance matrix
  int numpt = (int)points.size();
  double comat[3][3];
  Point mean = curr;
  for (int kr=0; kr<numpt; ++kr)
    mean += points[kr];
  mean /= (double)(numpt+1);

  vector<double> wgt(numpt+1);
  double div = (double)((numpt+1)*(numpt+1));
  wgt[0] = (mean.dist(curr))/div;
  for (int kr=0; kr<numpt; ++kr)
    {
      double dr2 = mean.dist2(points[kr]);
      wgt[kr+1] = exp(-dr2/div);
    }
  
  for (int ki=0; ki<3; ++ki)
    for (int kj=0; kj<3; ++kj)
      {
	double tmp = 0.0;
	tmp += wgt[0]*(curr[ki] - mean[ki])*(curr[kj] - mean[kj]);
	for (int kr=0; kr<numpt; ++kr)
	  {
	    tmp += wgt[kr+1]*(points[kr][ki] - mean[ki])*(points[kr][kj] - mean[kj]);
	  }
	comat[ki][kj] = tmp/(double)(numpt);
      }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ki = 0; ki < 3; ++ki) {
    for (int kj = 0; kj < 3; ++kj) {
      nmat.element(ki, kj) = comat[ki][kj];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    return;
  }
  
  // Singular values
  for (int ki=0; ki<3; ++ki)
    {
      lambda[ki] = diag.element(ki, ki);
      for (int kj=0; kj<3; ++kj)
	eigenvec[ki][kj] = V.element(kj, ki);
    }
}


//===========================================================================
void RevEngUtils::computeLocFunc(Point& curr, std::vector<Point>& points,
			       Point& vec1, Point& vec2, Point& normal, Point& mincvec,
			       double& minc, Point& maxcvec, double& maxc,
			       double& currdist, double& avdist)
//===========================================================================
{
  // Transform points to coordinate system given by vec1 (x-axis) and vec2 (y-axis)
  Matrix3D mat1, mat2, rotmat;
  Vector3D vec1_2(vec1[0], vec1[1], vec1[2]);
  Vector3D vec2_2(vec2[0], vec2[1], vec2[2]);
  Vector3D xaxis(1, 0, 0);
  Vector3D zaxis(0, 0, 1);
  mat1.setToRotation(vec1_2, xaxis);
  Vector3D v1 = mat1*vec1_2;
  Vector3D vec2_3 = mat1*vec2_2;
  mat2.setToRotation(vec2_3, zaxis);
  Vector3D v2 = mat2*vec2_3;
  rotmat = mat2*mat1;
  //rotmat.identity();

  // Perform rotation and sort parameter values and z-value
  int nmbpts = (int)points.size() + 1;
  vector<double> par(2*nmbpts);
  vector<double> zval(nmbpts);
  par[0] = curr[0];
  par[1] = curr[1];
  zval[0] = curr[2];
  for (int ki=1; ki<nmbpts; ++ki)
    {
      Point dv = points[ki-1] - curr;
      Vector3D dv2(dv[0], dv[1], dv[2]);
      Vector3D dvrot = rotmat*dv2;
      //Vector3D dvrot = mat2*dvrot0;
      par[2*ki] = curr[0] + dvrot[0];
      par[2*ki+1] = curr[1] + dvrot[1];
      zval[ki] = curr[2] + dvrot[2];
    }

  // Approximate z-component by biquadratic Bezier function in x and y
  int order = 3;
  shared_ptr<SplineSurface> locsf = RevEngUtils::surfApprox(zval, 1, par, order,
							      order, order, order);

  vector<double> coefs2(3*order*order);
  std::vector<double>::iterator cf = locsf->coefs_begin();
  for (int ka=0; ka<order; ++ka)
    {
      double vpar = locsf->basis_v().grevilleParameter(ka);
      for (int kb=0; kb<order; ++kb, ++cf)
	{
	  double upar = locsf->basis_u().grevilleParameter(kb);
	  coefs2[(ka*order+kb)*3] = upar;
	  coefs2[(ka*order+kb)*3+1] = vpar;
	  coefs2[(ka*order+kb)*3+2] = *cf;
	}
    }
  shared_ptr<SplineSurface> tmp(new SplineSurface(order, order, order, order, 
						  locsf->basis_u().begin(),
						  locsf->basis_v().begin(), &coefs2[0], 3));
#ifdef DEBUG
  int writesurface = 0;
  if (writesurface)
    {
      std::ofstream of("approx_sf.g2");
      tmp->writeStandardHeader(of);
      tmp->write(of);
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << 1 << std::endl;
      of << curr << std::endl;
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << nmbpts << std::endl;
      for (int ka=0; ka<nmbpts; ++ka)
	{
	  Point tmppt(par[2*ka], par[2*ka+1], zval[ka]);
	  of << tmppt << std::endl;
	}
    }
#endif
						  
  // Compute surface normal in curr
  vector<Point> der(3);
  locsf->point(der, par[0], par[1], 1);
  Vector3D norm(-der[1][0], -der[2][0], 1.0);
  norm.normalize();

  // Accuracy of approximation
  currdist = fabs(zval[0] - der[0][0]);
  avdist = currdist;
  for (int ki=1; ki<nmbpts; ++ki)
    {
      Point pos;
      locsf->point(pos, par[2*ki], par[2*ki+1]);
      avdist += fabs(zval[ki] - pos[0]);
    }
  avdist /= (double)nmbpts;
  
  // Compute principal curvatures in curr
  shared_ptr<SISLSurf> sislsf(GoSurf2SISL(*locsf, false));
  int left1 = 0, left2 = 0;
  int stat = 0;
  double k1, k2;
  double d1[2], d2[2];
  s2542(sislsf.get(), 0, 0, 0, &par[0], &left1, &left2, &k1, &k2, d1, d2, &stat);
  Vector3D du(1.0, 0.0, der[1][0]);
  Vector3D dv(0.0, 1.0, der[2][0]);
  Vector3D cvec1 = d1[0]*du + d1[1]*dv;
  Vector3D cvec2 = d2[0]*du + d2[1]*dv;
  minc = k1;
  maxc = k2;

  // Vector3D origin(par[0], par[1], zval[0]);
  // of << "410 1 0 4 0 0 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+norm << std::endl;

  // of << "410 1 0 4 0 55 155 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+cvec1 << std::endl;

  
  // of << "410 1 0 4 155 55 0 255" << std::endl;
  // of << "1" << std::endl;
  // of << origin << " " << origin+cvec2 << std::endl;

  
  
  // Transform results to original coordinate system
  Matrix3D mat3, mat4, rotmat2;
  mat4.setToRotation(zaxis, vec2_3);
  mat3.setToRotation(xaxis, vec1_2);
  rotmat2 = mat3*mat4;
  //rotmat2.identity();
  //Vector3D norm0 = mat4*norm;
  Vector3D norm2 = rotmat2*norm;
  normal = Point(norm2[0], norm2[1], norm2[2]);
  
  Vector3D cvec3 = rotmat2*cvec1;
  mincvec = Point(cvec3[0], cvec3[1],cvec3[2]); 
  Vector3D cvec4 = rotmat2*cvec2;
  maxcvec = Point(cvec4[0], cvec4[1],cvec4[2]); 

  int stop_break = 1;
}

//===========================================================================
void RevEngUtils::smoothSurf(shared_ptr<SplineSurface>& surf, int fixed)
//===========================================================================
{
  SmoothSurf smooth;
  int ncoef1 = surf->numCoefs_u();
  int ncoef2 = surf->numCoefs_v();
  vector<int> coef_known(ncoef1*ncoef2, 0);
  for (int ka=0; ka<ncoef1; ++ka)
    {
      for (int kb=0; kb<fixed; kb++)
	{
	  coef_known[kb*ncoef1+ka] = 1;
	  coef_known[(ncoef2-kb-1)*ncoef1+ka] = 1;
	}
    }
  for (int kb=fixed; kb<ncoef2-fixed; ++kb)
    {
      for (int ka=0; ka<fixed; ++ka)
	{
	  coef_known[kb*ncoef1+ka] = 1;
	  coef_known[(kb+1)*ncoef1-ka-1] = 1;
	}
    }


  int close1 = GeometryTools::analyzePeriodicityDerivs(*surf, 0, 1);
  int close2 = GeometryTools::analyzePeriodicityDerivs(*surf, 1, 1);
  int seem[2];
  seem[0] = (close1 >= 0) ? 1 : 0;
  seem[1] = (close2 >= 0) ? 1 : 0;
  shared_ptr<SplineSurface> surf2(surf->clone());
  smooth.attach(surf2, seem, &coef_known[0]);

  double wgt1 = 0.0, wgt2 = 1.0, wgt3 = 0.0;
  smooth.setOptimize(wgt1, wgt2, wgt3);
  shared_ptr<SplineSurface> surf3;
  smooth.equationSolve(surf3);
  std::swap(surf, surf3);
}

//===========================================================================
shared_ptr<SplineSurface>
RevEngUtils::surfApprox(vector<double>& data, int dim, vector<double>& param,
			int order1, int order2, int ncoef1, int ncoef2,
			bool close1, bool close2,
			int max_iter, double tol, double& maxd, double& avd, 
			int& num_out, vector<double>& parvals, double belt_frac)
//===========================================================================
{
  // Create initial spline space
  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::lowest();
  double vmin = std::numeric_limits<double>::max();
  double vmax = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<param.size(); ki+=2)
    {
      umin = std::min(umin, param[ki]);
      umax = std::max(umax, param[ki]);
      vmin = std::min(vmin, param[ki+1]);
      vmax = std::max(vmax, param[ki+1]);
    }
  double udel = (umax - umin);
  double vdel = (vmax - vmin);
  if (!close1)
    {
      umin -= belt_frac*udel;
      umax += belt_frac*udel;
    }
  if (!close2)
    {
      vmin -= belt_frac*vdel;
      vmax += belt_frac*vdel;
    }
  udel = (umax - umin)/(double)(ncoef1 - order1 + 1);
  vdel = (vmax - vmin)/(double)(ncoef2 - order2 + 1);
  vector<double> et1(order1+ncoef1);
  vector<double> et2(order2+ncoef2);
  vector<double> coef(ncoef1*ncoef2*dim, 0.0);
  for (int ka=0; ka<order1; ++ka)
    {
      et1[ka] = umin;
      et1[ka+ncoef1] = umax;
    }
  for (int ka=0; ka<order2; ++ka)
    {
      et2[ka] = vmin;
      et2[ka+ncoef2] = vmax;
    }
  for (int ka=0; ka<ncoef1-order1; ++ka)
    et1[order1+ka] = umin + (ka+1)*udel; 
  for (int ka=0; ka<ncoef2-order2; ++ka)
    et2[order2+ka] = vmin + (ka+1)*vdel; 

  
 shared_ptr<SplineSurface> surf(new SplineSurface(ncoef1, ncoef2, order1,
						  order2, &et1[0], 
						  &et2[0], &coef[0], dim));
 shared_ptr<SplineSurface> surf1(new SplineSurface(ncoef1, ncoef2, order1,
						   order2, &et1[0], 
						   &et2[0], &coef[0], dim));

 // Approximate
  SmoothSurf approx2;
  vector<int> coef_known(ncoef1*ncoef2, 0);
  int seem[2];
  seem[0] = close1 ? 1 : 0;
  seem[1] = close2 ? 1 : 0;
  int nmbpts = (int)data.size()/dim;
  vector<double> ptwgt(nmbpts, 1.0);
  approx2.attach(surf1, seem, &coef_known[0]);

  //double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  approx2.setOptimize(wgt1, wgt2, wgt3);
  approx2.setLeastSquares(data, param, ptwgt, approxwgt);
  shared_ptr<SplineSurface> surf3;
  approx2.equationSolve(surf3);
#ifdef DEBUG
  std::ofstream of("first_approx.g2");
  surf3->writeStandardHeader(of);
  surf3->write(of);
#endif
  
  ApproxSurf approx(surf3, data, param, dim, tol);
  //ApproxSurf approx(surf3, data, param, dim, tol, 0, false, true, 0, true);
  approx.setMBA(true);
  approx.setFixBoundary(false);
  double acc_frac = 0.6;
  approx.setAccuracyCrit(1, acc_frac);
 shared_ptr<SplineSurface> surf2 = approx.getApproxSurf(maxd, avd, num_out, max_iter);
 if (surf2.get())
   parvals = approx.getParvals();
 return surf2;
}

//===========================================================================
shared_ptr<SplineSurface> RevEngUtils::surfApprox(vector<double>& data, int dim,
						  vector<double>& param, int order1,
						  int order2, int nmb_coef1, int nmb_coef2,
						  double del)
//===========================================================================
{
  // Define spline space
  double umin, umax, vmin, vmax;
  umin = umax = param[0];
  vmin = vmax = param[1];
  for (size_t kj=2; kj<param.size(); kj+=2)
    {
      umin = std::min(umin, param[kj]);
      umax = std::max(umax, param[kj]);
      vmin = std::min(vmin, param[kj+1]);
      vmax = std::max(vmax, param[kj+1]);
    }
  umin -= del;
  umax += del;
  vmin -= del;
  vmax += del;

  return surfApprox(data, dim, param, order1, order2, nmb_coef1, nmb_coef2,
		    umin, umax, vmin, vmax);
}
  
//===========================================================================
shared_ptr<SplineSurface> RevEngUtils::surfApprox(vector<double>& data, int dim,
						  vector<double>& param, int order1,
						  int order2, int nmb_coef1, int nmb_coef2,
						  double umin, double umax,
						  double vmin, double vmax)
//===========================================================================
{
  double udel = (umax - umin)/(double)(nmb_coef1-order1+1);
  double vdel = (vmax - vmin)/(double)(nmb_coef2-order2+1);
  
  vector<double> knots1(order1+nmb_coef1), knots2(order2+nmb_coef2);
  for (int ka=0; ka<order1; ++ka)
    {
      knots1[ka] = umin;
      knots1[nmb_coef1+ka] = umax;
    }
  for (int ka=order1; ka<nmb_coef1; ++ka)
    knots1[ka] = umin + (ka-order1+1)*udel;
  
  for (int ka=0; ka<order2; ++ka)
    {
      knots2[ka] = vmin;
      knots2[nmb_coef2+ka] = vmax;
    }
  for (int ka=order2; ka<nmb_coef2; ++ka)
    knots2[ka] = vmin + (ka-order2+1)*vdel;
  

  vector<double> coefs((nmb_coef1+order1)*(nmb_coef2+order2)*dim, 0.0);
  shared_ptr<SplineSurface> bez(new SplineSurface(nmb_coef1, nmb_coef2, order1, order2, 
						  &knots1[0], &knots2[0], &coefs[0], dim));

  // Approximate
  SmoothSurf approx;
  vector<int> coef_known((nmb_coef1+order1)*(nmb_coef2+order2), 0);
  int seem[2];
  seem[0] = seem[1] = 0;
  int nmbpts = (int)data.size()/dim;
  vector<double> ptwgt(nmbpts, 1.0);
  approx.attach(bez, seem, &coef_known[0]);

  //double wgt1 = 0.0, wgt2 = 0.001, wgt3 = 0.001;
  double wgt1 = 0.0, wgt2 = 0.05, wgt3 = 0.05;
  double approxwgt = 1.0 - wgt1 - wgt2 - wgt3;
  approx.setOptimize(wgt1, wgt2, wgt3);
  approx.setLeastSquares(data, param, ptwgt, approxwgt);

  shared_ptr<SplineSurface> surf;
  approx.equationSolve(surf);
  return surf;
 }

//===========================================================================
bool RevEngUtils::parameterizeOnPrimary(vector<RevEngPoint*>& points,
					shared_ptr<ParamSurface> surf,
					vector<double>& data, vector<double>& param,
					int& inner1, int& inner2, bool& close1,
					bool& close2)
//===========================================================================
{
  double eps = 1.0e-2;
  double angtol = 0.1;

  // Make sure to use an untrimmed surface
  shared_ptr<ParamSurface> sf = surf;
  shared_ptr<BoundedSurface> bdsf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
  if (bdsf.get())
    sf = bdsf->underlyingSurface();
  shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder,ParamSurface>(sf);
  shared_ptr<Cone> cone = dynamic_pointer_cast<Cone,ParamSurface>(sf);
  shared_ptr<Sphere> sph = dynamic_pointer_cast<Sphere,ParamSurface>(sf);
  shared_ptr<Torus> torus = dynamic_pointer_cast<Torus,ParamSurface>(sf);
  close1 = (cyl.get() || cone.get() || sph.get() || torus.get());
  close2 = (sph.get() || torus.get());
  RectDomain dom = sf->containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double udel = u2 - u1;
  double v1 = dom.vmin();
  double v2 = dom.vmax();
  double vdel = v2 - v1;

  // Parameterize
  int dim = surf->dimension();
  param.resize(2*points.size());
  data.reserve(dim*points.size());
  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::lowest();
  double vmin = std::numeric_limits<double>::max();
  double vmax = std::numeric_limits<double>::lowest();
  double sfac = 0.5;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pos(xyz[0], xyz[1], xyz[2]);
      double upar, vpar, dist;
      Point close;
      sf->closestPoint(pos, upar, vpar, close, dist, eps);

      if (sf->isBounded() && dist > 0.1) //eps) TEST
	{
	  // Check boundary
	  double uparbd, vparbd, distbd;
	  Point closebd;
	  double seed[2];
	  seed[0] = upar;
	  seed[1] = vpar;
	  RectDomain dom = sf->containingDomain();
	  sf->closestBoundaryPoint(pos, uparbd, vparbd, closebd, distbd, eps, &dom, seed);
	  if (fabs(distbd - dist) < 0.1*dist)
	    {
	      // Angular check to see if a true closest point is found
	      Point norm;
	      sf->normal(norm, upar, vpar);
	      Point vec = pos - close;
	      double ang = norm.angle(vec);
	      ang = std::min(ang, fabs(M_PI - ang));
	      if (ang > angtol)
		return false;
	    }
	}

      // Check with seam
      // if (ki > 0 && close1 && ((upar-u1 < umin-upar && u2-umax < umin-u1) ||
      // 			       (u2-upar < upar-umax && umin-u1 < u2-umax)))
      // if (ki > 0 && close1 && (upar < umin || upar > umax) &&
      // 	  (fabs(umin-upar+udel) < fabs(upar-umin) || fabs(upar+udel-umax) < fabs(umax-upar)))
	  if (ki > 0 && close1 && (upar < umin || upar > umax) &&
	      std::min(fabs(umin-upar+udel),fabs(upar+udel-umax)) < std::min(fabs(upar-umin),fabs(upar-umax)))
	{
	  // if (ki >= 2)
	  //   {
	  //     // Compare distances
	  //     double d1 = points[ki-2]->pntDist(points[ki-1]);
	  //     double d2 = points[ki-1]->pntDist(points[ki]);
	  //     Point p1(param[2*(ki-2)], param[2*(ki-2)+1]);
	  //     Point p2(param[2*(ki-1)], param[2*(ki-1)+1]);
	  //     Point p3(upar,vpar);
	  //     double dp1 = p1.dist(p2);
	  //     double dp2 = p2.dist(p3);
	  //     if (dp1/dp2 < sfac*d1/d2)
	  // 	{
		  if (upar-u1 < u2-upar)
		    upar += (u2-u1);
		  else
		    upar -= (u2-u1);
	    // 	}
									 
	    //   int stop_break1 = 1;
	    // }
	}
      // if (ki > 0 && close2 && ((vpar-v1 < vmin-vpar && v2-vmax < vmin-v1) ||
      // 			       (v2-vpar < vpar-vmax && vmin-v1 < v2-vmax)))
	  // if (ki > 0 && close2 && (vpar < vmin || vpar > vmax) &&
	  //     (fabs(vmin-vpar+vdel) < fabs(vpar-vmin) || fabs(vpar+vdel-vmax) < fabs(vpar-vmax)))
	  if (ki > 0 && close2 && (vpar < vmin || vpar > vmax) &&
	      std::min(fabs(vmin-vpar+vdel),fabs(vpar+vdel-vmax)) < std::min(fabs(vpar-vmin),fabs(vpar-vmax)))
	{
	  // if (ki >= 2)
	  //   {
	  //     // Compare distances
	  //     double d1 = points[ki-2]->pntDist(points[ki-1]);
	  //     double d2 = points[ki-1]->pntDist(points[ki]);
	  //     Point p1(param[2*(ki-2)], param[2*(ki-2)+1]);
	  //     Point p2(param[2*(ki-1)], param[2*(ki-1)+1]);
	  //     Point p3(upar,vpar);
	  //     double dp1 = p1.dist(p2);
	  //     double dp2 = p2.dist(p3);
	  //     if (dp1/dp2 < sfac*d1/d2)
	  // 	{
		  if (vpar-v1 < v2-vpar)
		    vpar += (v2-v1);
		  else
		    vpar -= (v2-v1);
	    // 	}
									 
	    //   int stop_break2 = 1;
	    // }
	}
      param[2*ki] = upar;
      param[2*ki+1] = vpar;
      umin = std::min(umin, upar);
      umax = std::max(umax, upar);
      vmin = std::min(vmin, vpar);
      vmax = std::max(vmax, vpar);
      data.insert(data.end(), pos.begin(), pos.end());
    }

  if (close1 && umax-umin < 2*M_PI-eps)
    close1 = false;
  if (close2 && vmax-vmin < 2*M_PI-eps)
    close2 = false;
  
  // Reparameterize for surfaces with circular properties
  inner1 = inner2 = 0;
  if (cyl.get())
    {
      double rad = cyl->getRadius();
      for (size_t ki=0; ki<param.size(); ki+=2)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
    }

  if (cone.get())
    {
      double rad1 = cone->radius(0.0, vmin);
      double rad2 = cone->radius(0.0, vmax);
      double rad = 0.5*(rad1 + rad2);
      for (size_t ki=0; ki<param.size(); ki+=2)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
    }
  
  if (sph.get())
    {
      double rad = sph->getRadius();
      for (size_t ki=0; ki<param.size(); ++ki)
	param[ki] *= rad;
      inner1 = 2.0*(umax - umin)/M_PI;
      inner2 = 2.0*(vmax - vmin)/M_PI;
    }

  if (torus.get())
    {
      double rad1 = torus->getMajorRadius();
      double rad2 = torus->getMinorRadius();
      for (size_t ki=0; ki<param.size(); ki+=2)
	{
	  param[ki] *= rad1;
	  param[ki+1] *= rad2;
	}
      inner1 = 2.0*(umax - umin)/M_PI;
      inner2 = 2.0*(vmax - vmin)/M_PI;
    }
  return true;
}

//===========================================================================
void RevEngUtils::parameterizeWithPlane(vector<RevEngPoint*>& pnts,
					const BoundingBox& bbox,
					const Point& vec1, const Point& vec2,
					vector<double>& data, vector<double>& param)
//===========================================================================
{
  vector<Point> pos(pnts.size());
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      Vector3D xyz = pnts[ki]->getPoint();
      pos[ki] = Point(xyz[0], xyz[1], xyz[2]);
    }

  parameterizeWithPlane(pos, bbox, vec1, vec2, data, param);
}

//===========================================================================
void RevEngUtils::parameterizeWithPlane(vector<Point>& pnts, const BoundingBox& bbox,
					const Point& vec1, const Point& vec2,
					vector<double>& data, vector<double>& param)
//===========================================================================
{
  double eps = 1.0e-6;
  Point mid = 0.5*(bbox.low() + bbox.high());
  int dim = mid.dimension();
  double diag = bbox.low().dist(bbox.high());
  int order = 2;
  double et[4];
  et[0] = et[1] = -diag;
  et[2] = et[3] = diag;
  vector<double> coefs;
  int sgn1, sgn2;
  int ka, kb;
  for (kb=0, sgn2=-1; kb<2; ++kb, sgn2=1)
    for (ka=0, sgn1=-1; ka<2; ++ka, sgn1=1)
      {
	Point pos = mid+sgn1*diag*vec1+sgn2*diag*vec2;
	coefs.insert(coefs.end(), pos.begin(), pos.end());
      }

  shared_ptr<SplineSurface> surf(new SplineSurface(order, order, order, order, &et[0], 
						   &et[0], coefs.begin(), dim));
#ifdef DEBUG
  std::ofstream of("parplane.g2");
  surf->writeStandardHeader(of);
  surf->write(of);
#endif
  
  param.resize(2*pnts.size());
  data.reserve(dim*pnts.size());
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(pnts[ki], upar, vpar, close, dist, eps);
      param[2*ki] = upar;
      param[2*ki+1] = vpar;
      data.insert(data.end(), pnts[ki].begin(), pnts[ki].end());
    }
}

//===========================================================================
void RevEngUtils::computeAxis(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  double Cmat[3][3];
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      Cmat[ka][kb] = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (int ka=0; ka<3; ++ka)
	for (int kb=0; kb<3; ++kb)
	  {
	    for (auto it=start; it!=end; ++it)
	      {
		RevEngPoint *pt = *it;
		Point norm1 = pt->getLocFuncNormal();
		Point norm2 = pt->getTriangNormal();
		Point norm = norm1; //0.5*(norm1 + norm2);
		Cmat[ka][kb] += norm[ka]*norm[kb];
		//Cmat[ka][kb] += norm2[ka]*norm2[kb];
	      }
	  }
    }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}


//===========================================================================
void RevEngUtils::computeAxis(vector<Point>& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  double Cmat[3][3];
  for (int ka=0; ka<3; ++ka)
    for (int kb=0; kb<3; ++kb)
      {
	Cmat[ka][kb] = 0.0;
	for (size_t ki=0; ki<points.size(); ++ki)
	  {
	    Cmat[ka][kb] += points[ki][ka]*points[ki][kb];
	  }
      }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}


//===========================================================================
void RevEngUtils::coneAxis(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			      Point& axis, Point& Cx, Point& Cy)
//===========================================================================
{
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  double wgt = 1.0/(double)numpt;

  Point mid(3);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Point norm = pt->getLocFuncNormal();
	  mid += wgt*norm;
	}
    }

  double Cmat[3][3];
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (int ka=0; ka<3; ++ka)
	for (int kb=0; kb<3; ++kb)
	  {
	    Cmat[ka][kb] = 0.0;
	    for (auto it=start; it!=end; ++it)
	      {
		RevEngPoint *pt = *it;
		Point norm = pt->getLocFuncNormal();
		Point vec = norm - mid;
		Cmat[ka][kb] += vec[ka]*vec[kb];
	      }
	  }
    }
  
  // Compute singular values
  NEWMAT::Matrix nmat;
  nmat.ReSize(3, 3);
  for (int ka = 0; ka < 3; ++ka) {
    for (int kb = 0; kb < 3; ++kb) {
      nmat.element(ka, kb) = Cmat[ka][kb];
    }
  }
      
  static NEWMAT::DiagonalMatrix diag;
  static NEWMAT::Matrix V;
  try {
    NEWMAT::SVD(nmat, diag, nmat, V);
  } catch(...) {
    //std::cout << "Exception in SVD" << std::endl;
    exit(-1);
  }
  Cx = Point(V.element(0,0), V.element(1,0), V.element(2,0));
  Cy = Point(V.element(0,1), V.element(1,1), V.element(2,1));
  axis = Point(V.element(0,2), V.element(1,2), V.element(2,2));

}

//===========================================================================
void RevEngUtils::coneApex(vector<pair<vector<RevEngPoint*>::iterator,
			      vector<RevEngPoint*>::iterator> >& points,
			   Point axis, Point& apex, double& phi)
//===========================================================================
{
  double Mmat[3][3], Mi[3][3];
  double bvec[3], bi[3];
  for (int ka=0; ka<3; ++ka)
    {
      bvec[ka] = 0.0;
      for (int kb=0; kb<3; ++kb)
	Mmat[ka][kb] = 0.0;
    }

  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    nmb += (int)(points[ki].second - points[ki].first);

  vector<Point> dird;
  vector<Point> pp;
  dird.reserve(nmb);
  pp.reserve(nmb);
  double wg = 1.0/(double)nmb;
    for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Point norm = pt->getLocFuncNormal();
	  Point tmp = norm.cross(axis);
	  Point di = tmp.cross(norm);
	  di.normalize_checked();
	  dird.push_back(di);
	  Vector3D xyz = pt->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  pp.push_back(pos);

	  for (int ka=0; ka<3; ++ka)
	    for (int kb=0; kb<3; ++kb)
	      {
		if (ka == kb)
		  continue;
		Mi[ka][kb] = -di[ka]*di[kb];
	      }
	  Mi[0][0] = di[1]*di[1] + di[2]*di[2];
	  Mi[1][1] = di[0]*di[0] + di[2]*di[2];
	  Mi[2][2] = di[0]*di[0] + di[1]*di[1];

	  bi[0] = pos[0]*di[1]*di[1] - pos[1]*di[0]*di[1] - pos[2]*di[0]*di[2] + pos[0]*di[2]*di[2];
	  bi[1] = pos[1]*di[2]*di[2] - pos[2]*di[1]*di[2] - pos[0]*di[1]*di[0] + pos[1]*di[0]*di[0];
	  bi[2] = pos[2]*di[0]*di[0] - pos[0]*di[2]*di[0] - pos[1]*di[2]*di[1] + pos[2]*di[1]*di[1];
	  
	  for (int ka=0; ka<3; ++ka)
	    {
	      bvec[ka] += wg*bi[ka];
	      for (int kb=0; kb<3; ++kb)
		Mmat[ka][kb] += wg*Mi[ka][kb];
	    }
	}
    }

#ifdef DEBUG
    std::ofstream of("directions.g2");
    of << "410 1 0 4 155 200 0 255" << std::endl;
    of << dird.size() << std::endl;
    for (size_t ki=0; ki<dird.size(); ++ki)
      of << pp[ki] << " " << pp[ki]+dird[ki] << std::endl;
#endif
    
    double det = 0.0;
    int sgn = 1;
    int ka, kb, kc;
    double ax=0.0, ay=0.0, az=0.0;
    // for (ka=0; ka<3; ++ka, sgn*=-1)
    //   {
    // 	kb = (ka+1)%3;
    // 	kc = (kb+1)%3;
    // 	det += sgn*Mmat[0][ka]*(Mmat[1][kb]*Mmat[2][kc]-Mmat[2][kb]*Mmat[1][kc]);
    // 	ax += sgn*bvec[ka]*(Mmat[1][kb]*Mmat[2][kc]-Mmat[2][kb]*Mmat[1][kc]);
    // 	ay += sgn*Mmat[0][ka]*(bvec[kb]*Mmat[2][kc]-Mmat[2][kb]*bvec[kc]);
    // 	az += sgn*Mmat[0][ka]*(Mmat[1][kb]*bvec[kc]-bvec[kb]*Mmat[1][kc]);
    //   }
    // apex = Point(ax/det, ay/det, az/det);

    double det2 = Mmat[0][0]*(Mmat[1][1]*Mmat[2][2] - Mmat[1][2]*Mmat[2][1]) -
      Mmat[0][1]*(Mmat[1][0]*Mmat[2][2] - Mmat[1][2]*Mmat[2][0]) +
      Mmat[0][2]*(Mmat[1][0]*Mmat[2][1] - Mmat[1][1]*Mmat[2][0]);
    double ax2 = bvec[0]*(Mmat[1][1]*Mmat[2][2] - Mmat[1][2]*Mmat[2][1]) -
      bvec[1]*(Mmat[1][0]*Mmat[2][2] - Mmat[1][2]*Mmat[2][0]) +
      bvec[2]*(Mmat[1][0]*Mmat[2][1] - Mmat[1][1]*Mmat[2][0]);
    double ay2 = Mmat[0][0]*(bvec[1]*Mmat[2][2] - bvec[2]*Mmat[2][1]) -
      Mmat[0][1]*(bvec[0]*Mmat[2][2] - bvec[2]*Mmat[2][0]) +
      Mmat[0][2]*(bvec[0]*Mmat[2][1] - bvec[1]*Mmat[2][0]);
    double az2 = Mmat[0][0]*(Mmat[1][1]*bvec[2] - Mmat[1][2]*bvec[1]) -
      Mmat[0][1]*(Mmat[1][0]*bvec[2] - Mmat[1][2]*bvec[0]) +
      Mmat[0][2]*(Mmat[1][0]*bvec[1] - Mmat[1][1]*bvec[0]);
    apex = (fabs(det2) < 1.0e-6) ? Point(0.0, 0.0, 0.0) : Point(ax2/det2, ay2/det2, az2/det2);
    
    // // std::cout << det << " " << det2 << std::endl;
    // for (int ka=0; ka<3; ++ka)
    //   {
    // 	double tmp = 0.0;
    // 	for (kb=0; kb<3; ++kb)
    // 	  tmp += Mmat[ka][kb]*apex[kb];
    // 	std::cout << tmp << " " << bvec[ka] << std::endl;
    //   }

    double nom=0.0, denom=0.0;
    for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D xyz = pt->getPoint();
	  Point pos(xyz[0], xyz[1], xyz[2]);
	  Point tmp1 = pos - apex;
	  Point tmp2 = tmp1.cross(axis);
	  nom += tmp2.length();
	  denom += tmp1*axis;
	}
    }
    
    double tanphi = nom/denom;
    phi = atan(tanphi);
}

//===========================================================================
void
RevEngUtils::computeSphereProp(vector<pair<vector<RevEngPoint*>::iterator,
			       vector<RevEngPoint*>::iterator> >& points,
			       Point& centre, double& radius)
//===========================================================================
{
  double Amat[4][4];
  double bvec[4];
  // vector<vector<double> > Amat(4, vector<double>(4,0.0));
  // vector<double> bvec(4, 0.0);
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(4);
      
  for (int ka=0; ka<4; ++ka)
    {
      for (int kb=0; kb<4; ++kb)
  	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  double wgt = 1.0/(double)numpt;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      size_t kj = 0;
      for (auto it=start; it!=end; ++it, ++kj)
	{
	  Vector3D curr = (*it)->getPoint();
	  A1[kj][0] = 2*curr[0];
	  A1[kj][1] = 2*curr[1];
	  A1[kj][2] = 2*curr[2];
	  A1[kj][3] = 1.0;
	  b1[kj] = curr.length2();
	}
    }

  for (int ka=0; ka<4; ++ka)
    {
      for (int kb=0; kb<4; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat[ka][kb] += wgt*A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec[ka] += wgt*A1[kr][ka]*b1[kr];
    }

  // double detA = 0.0;
  // double bx[3];
  // bx[0] = bx[1] = bx[2] = 0.0;
  // int sgn1 = 1, sgn2 = 1;
  // for (int kb=0; kb<4; ++kb, sgn1*=(-1))
  //   {
  //     for (int kc=0; kc<4; ++kc)
  // 	{
  // 	  if (kc == kb)
  // 	    continue;
  // 	  int ka1 = (kb == 0);
  // 	  int ka2 = 3 - (kb == 3);
  // 	  detA += sgn1*Amat[0][kb]*
  // 	    (sgn2*Amat[1][kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			       Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[0] += sgn1*bvec[kb]*
  // 	    (sgn2*Amat[1][kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			       Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[1] += sgn1*Amat[0][kb]*
  // 	    (sgn2*bvec[kc]*(Amat[2][ka1]*Amat[3][ka2]-
  // 			    Amat[3][ka1]*Amat[2][ka2]));
  // 	  bx[2] += sgn1*Amat[0][kv]*
  // 	    (sgn2*Amat[1][kc]*(bvec[ka1]Amat[3][ka2]-
  // 			       bvec[ka2]*Amat[3][ka1]));
  // 	  bx[3] += sgn1*Amat[0][kb]*
  // 	    sgn2*Amat[1][kc]*((Amat[2][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
  // 	  sgn2 += -1;
  // 	}
  //   }
  // double sx = bx[0]/detA;
  // double sy = bx[1]/detA;
  // double sz = bx[2]/detA;
  // double r2 = bx[3]/detA;
  LUsolveSystem(Amat, 4, &bvec[0]);
  double sx = bvec[0];
  double sy = bvec[1];
  double sz = bvec[2];
  double r2 = bvec[3];

  centre = Point(sx,sy,sz);

  radius = (r2 + sx*sx + sy*sy + sz*sz < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy + sz*sz);
 }

//===========================================================================
void RevEngUtils::computeCylPosRadius(vector<pair<vector<RevEngPoint*>::iterator,
				      vector<RevEngPoint*>::iterator> >& points,
				      Point& low, Point& high,
				      Point& axis, Point& Cx, Point& Cy,
				      Point& pos, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      numpt += (points[ki].second - points[ki].first);
    }
  Point mid = 0.5*(low+high);
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      size_t kj = 0;
      for (auto it=start; it!=end; ++it, ++kj)
	{
	  Vector3D curr = (*it)->getPoint();
	  Point curr2(curr[0], curr[1], curr[2]);
	  curr2 -= mid;
	  double pxy[3];
	  double px = curr2*Cx;
	  double py = curr2*Cy;
	  pxy[0] = 2*px;
	  pxy[1] = 2*py;
	  pxy[2] = 1.0;
	  double plen2 = px*px + py*py;
	  A1[kj][0] = 2*px;
	  A1[kj][1] = 2*py;
	  A1[kj][2] = 1.0;
	  b1[kj] = plen2;
	  for (int ka=0; ka<3; ++ka)
	    {
	      for (int kb=0; kb<3; ++kb)
		Amat[ka][kb] += pxy[ka]*pxy[kb];
	      bvec[ka] += pxy[ka]*plen2;
	    }
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  Point pos2 = sx*Cx + sy*Cy;
  pos2 += mid;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
  double len = low.dist(high);
  Point vec = mid - pos2;
  Point ax = axis;
  ax.normalize();
  pos = pos2 + (vec*ax)*ax;
 }

 

//===========================================================================
void RevEngUtils::computeCircPosRadius(vector<Point>& points,
				      const Point& axis, const Point& Cx, 
				      const Point& Cy, Point& pos, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = points.size();
  double wgt = 1.0/(double)numpt;
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  Point mid(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<points.size(); ++ki)
      mid += wgt*points[ki];
  
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double pxy[3];
      double px = (points[ki]-mid)*Cx;
      double py = (points[ki]-mid)*Cy;
      pxy[0] = 2*px;
      pxy[1] = 2*py;
      pxy[2] = 1.0;
      double plen2 = px*px + py*py;
      A1[ki][0] = 2*px;
      A1[ki][1] = 2*py;
      A1[ki][2] = 1.0;
      b1[ki] = plen2;
      for (int ka=0; ka<3; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    Amat[ka][kb] += pxy[ka]*pxy[kb];
	  bvec[ka] += pxy[ka]*plen2;
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  if (fabs(detA) < 1.0e-12 || std::isnan(detA))
    THROW("Circle computation fail");
  //std::cout << "Circposradius, detA:" << detA << std::endl;
    
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  pos = sx*Cx + sy*Cy;
  pos += mid;
  pos -= ((pos-mid)*axis)*axis;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
 }



//===========================================================================
void RevEngUtils::computeRadius(vector<Point>& points, Point& axis, 
				Point& Cx, Point& Cy, double& radius)
//===========================================================================
{
  double Amat[3][3];
  double bvec[3];
  size_t numpt = points.size();
  double wgt = 1.0/numpt;
  
  vector<vector<double> > A1(numpt);
  vector<double> b1(numpt);
  for (size_t kj=0; kj<numpt; ++kj)
    A1[kj].resize(3);
      
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat[ka][kb] = 0.0;
      bvec[ka] = 0.0;
    }

  Point mid(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<points.size(); ++ki)
      mid += wgt*points[ki];
  
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double pxy[3];
      double px = (points[ki]-mid)*Cx;
      double py = (points[ki]-mid)*Cy;
      pxy[0] = 2*px;
      pxy[1] = 2*py;
      pxy[2] = 1.0;
      double plen2 = px*px + py*py;
      A1[ki][0] = 2*px;
      A1[ki][1] = 2*py;
      A1[ki][2] = 1.0;
      b1[ki] = plen2;
      for (int ka=0; ka<3; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    Amat[ka][kb] += pxy[ka]*pxy[kb];
	  bvec[ka] += pxy[ka]*plen2;
	}
    }

  double Amat2[3][3];
  double bvec2[3];
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	Amat2[ka][kb] = 0.0;
      bvec2[ka] = 0.0;
    }
  for (int ka=0; ka<3; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	{
	  for (size_t kr=0; kr<numpt; ++kr)
	    Amat2[ka][kb] += A1[kr][ka]*A1[kr][kb];
	}
      for (size_t kr=0; kr<numpt; ++kr)
	bvec2[ka] += A1[kr][ka]*b1[kr];
    }

  double detA = 0.0;
  double bx[3];
  bx[0] = bx[1] = bx[2] = 0.0;
  int sgn = 1;
  for (int kb=0; kb<3; ++kb, sgn*=(-1))
    {
      int ka1 = (kb == 0);
      int ka2 = 2 - (kb == 2);
      detA += sgn*Amat[0][kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[0] += sgn*bvec[kb]*(Amat[1][ka1]*Amat[2][ka2]-Amat[2][ka1]*Amat[1][ka2]);
      bx[1] += sgn*Amat[0][kb]*(bvec[ka1]*Amat[2][ka2]-bvec[ka2]*Amat[2][ka1]);
      bx[2] += sgn*Amat[0][kb]*(Amat[1][ka1]*bvec[ka2]-Amat[1][ka2]*bvec[ka1]);
    }
  if (fabs(detA) < 1.0e-12)
    THROW("Circle with infinite radius");
  
  double sx = bx[0]/detA;
  double sy = bx[1]/detA;
  double r2 = bx[2]/detA;

  radius = (r2 + sx*sx + sy*sy < 0.0) ? 0.0 : sqrt(r2 + sx*sx + sy*sy);
 }


//===========================================================================
void RevEngUtils::computePlane(vector<Point>& points, Point normal,
			       Point mainaxis[3],
			       Point& pos, Point& norm, Point& Cx, Point& Cy)
//===========================================================================
{
  Point pos0(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    pos0 += wgt*points[ki];
  
  ImplicitApprox impl;
  impl.approxPoints(points, 1);
  bool found = impl.projectPoint(pos0, normal, pos, norm);
  if (!found)
    {
      double eps = 1.0e-4;
      double ang_min = 0.25*M_PI;
      pos = pos0;
      Point vec = points[0] - pos;
      norm = Point(0.0, 0.0, 0.0);
      size_t ki, kj;
      for (ki=1; ki<points.size(); ++ki)
	{
	  if (vec.length() > eps)
	    break;
	  vec = points[ki] - pos;
	  for (kj=ki+1; kj<points.size(); ++kj)
	    {
	      Point vec2 = points[kj] - pos;
	      double ang = vec.angle(vec2);
	      if (ang >= ang_min)
		{
		  Point norm2 = vec.cross(vec2);
		  if (normal*norm2 < 0.0)
		    norm2 *= -1;
		  norm2.normalize();
		  norm += norm2;
		  break;
		}
	    }
	}
      norm.normalize();
    }
  if (normal*norm < 0.0)
    norm *= -1.0;

  // Define x-axis
  int ix = -1;
  double minang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(norm);
      ang = std::min(ang, M_PI-ang);
      if (ang < minang)
	{
	  minang = ang;
	  ix = ka;
	}
    }

  Cy = mainaxis[(ix+1)%3].cross(norm);
  Cx = norm.cross(Cy);
}

//===========================================================================
void RevEngUtils::computeLine(vector<Point>& points, Point& pos, Point& dir)
//===========================================================================
{
  size_t num = points.size();
  pos = Point(0.0, 0.0, 0.0);
  dir = Point(0.0, 0.0, 0.0);
  if (num == 0)
    return;
  
  double fac = 1.0/(double)num;
  BoundingBox bb(3);
  for (size_t ki=0; ki<num; ++ki)
    {
      pos += fac*points[ki];
      bb.addUnionWith(points[ki]);
    }

  Point dir0 = bb.high() - bb.low();
  for (size_t ki=0; ki<num; ++ki)
    {
      Point tmp = points[ki] - pos;
      (void)tmp.normalize_checked();
      if (tmp*dir0 < 0.0)
	tmp *= -1;
      dir += fac*tmp;
    }
  (void)dir.normalize_checked();
}

//===========================================================================
void RevEngUtils::projectToPlane(vector<RevEngPoint*>& points,
				 Point& axis, Point& mid, std::vector<Point>& projected,
				 double& maxdist, double& avdist, double dlen)
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > points2;
  points2.push_back(std::make_pair(points.begin(), points.end()));
  projectToPlane(points2, axis, mid, projected, maxdist, avdist, dlen);
}

//===========================================================================
void RevEngUtils::projectToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
				 std::vector<RevEngPoint*>::iterator> >& points,
				 Point& axis, Point& mid, std::vector<Point>& projected,
				 double& maxdist, double& avdist, double dlen)
//===========================================================================
{
  maxdist = 0.0;
  avdist = 0.0;
  int nmb = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Point curr(pnt[0], pnt[1], pnt[2]);
	  Point curr2 = curr - mid;
	  double dd = curr2*axis;
	  if (dlen > 0.0 && dd > dlen)
	    continue;
	  curr2 -= (dd*axis);
	  curr2 += mid;
	  projected.push_back(curr2);
	  double dist = curr.dist(curr2);
	  maxdist = std::max(maxdist, dist);
	  avdist += dist;
	  nmb++;
	}
    }
  avdist /= (double)nmb;
}

//===========================================================================
void RevEngUtils::rotateToPlane(vector<pair<vector<RevEngPoint*>::iterator,
				vector<RevEngPoint*>::iterator> >& points,
				Point& xvec, Point& axis, Point& mid,
				vector<Point>& rotated)
//===========================================================================
{
  Point yvec = xvec.cross(axis);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      vector<RevEngPoint*>::iterator start = points[ki].first;
      vector<RevEngPoint*>::iterator end = points[ki].second;
      for (auto it=start; it!=end; ++it)
	{
	  RevEngPoint *pt = *it;
	  Vector3D pnt = pt->getPoint();
	  Point curr(pnt[0], pnt[1], pnt[2]);
	  curr -= mid;
	  pnt = Vector3D(curr[0], curr[1], curr[2]);
	  curr -= (curr*axis)*axis;
	  double angle = curr.angle(xvec);
	  if (curr*yvec < 0.0)
	    angle *= -1.0;
	  Matrix3D mat;
	  mat.setToRotation(angle, axis[0], axis[1], axis[2]);
	  Vector3D pnt2 = mat*pnt;
	  rotated.push_back(mid + Point(pnt2[0], pnt2[1], pnt2[2]));
	}
    }
 }

//===========================================================================
void RevEngUtils::rotateToPlane(vector<Point>& points,
				Point& xvec, Point& axis, Point& mid,
				vector<Point>& rotated)
//===========================================================================
{
  rotated.resize(points.size());
  Point yvec = xvec.cross(axis);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Point curr = points[ki];
      curr -= mid;
      Vector3D pnt(curr[0], curr[1], curr[2]);
      curr -= (curr*axis)*axis;
      double angle = curr.angle(xvec);
      if (curr*yvec < 0.0)
	angle *= -1.0;
      Matrix3D mat;
      mat.setToRotation(angle, axis[0], axis[1], axis[2]);
      Vector3D pnt2 = mat*pnt;
      rotated[ki] = mid + Point(pnt2[0], pnt2[1], pnt2[2]);
    }
}

//===========================================================================
void RevEngUtils::distToSurf(vector<RevEngPoint*>::iterator start,
			     vector<RevEngPoint*>::iterator end,
			     shared_ptr<ParamSurface> surf, double tol,
			     double& maxdist, double& avdist,
			     int& num_inside, int& num_inside2,
			     vector<RevEngPoint*>& in, vector<RevEngPoint*>& out,
			     vector<double>& parvals,
			     vector<pair<double,double> >& distang,
			     double angtol)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = num_inside2 = 0;
  int num = (int)(end-start);
  parvals.resize(2*num);
  distang.resize(num);
  in.clear();
  out.clear();
  double *seed = 0;
  double seed2[2];
  Point prev;
  double fac = 100.0;
  double upar, vpar, dist;
  Point close;
  Point norm1, norm2, norm3;
  double ang, ang2;
  double dfac = 1.0/(double)num;
  size_t ki=0;
  for (auto it=start; it!=end; ++it, ++ki)
    {
      Vector3D xyz = (*it)->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      if (prev.dimension() == pnt.dimension() && prev.dist(pnt) < fac*tol)
	seed = seed2;
      
      surf->closestPoint(pnt, upar, vpar, close, dist, eps, 0, seed);
      parvals[2*ki] = upar;
      parvals[2*ki+1] = vpar;
      surf->normal(norm1, upar, vpar);
      norm2 = (*it)->getLocFuncNormal();
      norm3 = (*it)->getTriangNormal();
      maxdist = std::max(maxdist, dist);
      avdist += dfac*dist;
      ang = norm1.angle(norm2);
      ang2 = norm1.angle(norm3);
      ang = std::min(std::min(M_PI-ang, ang), std::min(M_PI-ang2,ang2));
      distang[ki] = std::make_pair(dist, ang);
      if (dist <= tol)
	{
	  ++num_inside2;
	  if (angtol < 0.0 || ang <= angtol)
	    {
	      in.push_back(*it);
		++num_inside;
	    }
	  else
	    out.push_back(*it);
	}
      else
	out.push_back(*it);
      seed2[0] = upar;
      seed2[1] = vpar;
      prev = pnt;
    }
}

//===========================================================================
void RevEngUtils::distToSurf(vector<RevEngPoint*>::iterator start,
			     vector<RevEngPoint*>::iterator end,
			     shared_ptr<ParamSurface> surf, double tol,
			     double& maxdist, double& avdist,
			     int& num_inside, int& num_inside2,
			     vector<double>& parvals,
			     vector<pair<double,double> >& distang,
			     double angtol)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = num_inside2 = 0;
  int num = (int)(end-start);
  parvals.resize(2*num);
  distang.resize(num);
  double *seed = 0;
  double seed2[2];
  Point prev;
  double fac = 100.0;
  double upar, vpar, dist;
  Point close;
  Point norm1, norm2, norm3;
  double ang, ang2;
  double dfac = 1.0/(double)num;
  size_t ki=0;
  for (auto it=start; it!=end; ++it, ++ki)
    {
      Vector3D xyz = (*it)->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      if (prev.dimension() == pnt.dimension() && prev.dist(pnt) < fac*tol)
	seed = seed2;
      
      surf->closestPoint(pnt, upar, vpar, close, dist, eps, 0, seed);
      parvals[2*ki] = upar;
      parvals[2*ki+1] = vpar;
      surf->normal(norm1, upar, vpar);
      norm2 = (*it)->getLocFuncNormal();
      norm3 = (*it)->getTriangNormal();
      maxdist = std::max(maxdist, dist);
      avdist += dfac*dist;
      ang = norm1.angle(norm2);
      ang2 = norm1.angle(norm3);
      ang = std::min(std::min(M_PI-ang, ang), std::min(M_PI-ang2,ang2));
      distang[ki] = std::make_pair(dist, ang);
      if (dist <= tol)
	{
	  ++num_inside2;
	  if (angtol < 0.0 || ang <= angtol)
	    ++num_inside;
	}
      seed2[0] = upar;
      seed2[1] = vpar;
      prev = pnt;
    }
}

//===========================================================================
void RevEngUtils::distToSurf(vector<RevEngPoint*>& points,
			     shared_ptr<ParamSurface> surf, double tol,
			     double angtol, double& maxdist, double& avdist,
			     int& inside, int& inside2,
			     vector<pair<double,double> >& dist_ang)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  inside = inside2 = 0;
  int num = (int)points.size();
  dist_ang.resize(num);
  double *seed = 0;
  double seed2[2];
  Point prev;
  double fac = 100.0;
  double upar, vpar, dist;
  Point close;
  Point norm1, norm2, norm3;
  double ang, ang2;
  double dfac = 1.0/(double)num;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      if (prev.dimension() == pnt.dimension() && prev.dist(pnt) < fac*tol)
	seed = seed2;
      
      surf->closestPoint(pnt, upar, vpar, close, dist, eps, 0, seed);
      surf->normal(norm1, upar, vpar);
      norm2 = points[ki]->getLocFuncNormal();
      norm3 = points[ki]->getTriangNormal();
      maxdist = std::max(maxdist, dist);
      avdist += dfac*dist;
      ang = norm1.angle(norm2);
      ang2 = norm1.angle(norm3);
      ang = std::min(std::min(M_PI-ang, ang), std::min(M_PI-ang2,ang2));
      dist_ang[ki] = std::make_pair(dist, ang);
      if (dist <= tol)
	{
	  ++inside2;
	  if (angtol < 0.0 || ang <= angtol)
	    {
	      ++inside;
	    }
	}
      seed2[0] = upar;
      seed2[1] = vpar;
      prev = pnt;
    }
}

//===========================================================================
void RevEngUtils::distToSurf(vector<Point>& points,
			     shared_ptr<ParamSurface> surf, double tol,
			     double& maxdist, double& avdist, int& num_inside,
			     vector<double>& distance)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = (int)points.size();
  distance.resize(num);
  double fac = 1.0/(double)num;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double upar, vpar, dist;
      Point close;
      surf->closestPoint(points[ki], upar, vpar, close, dist, eps);
      maxdist = std::max(maxdist, dist);
      avdist += fac*dist;
      distance[ki] = dist;
      if (dist <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
    }
}

//===========================================================================
void RevEngUtils::distToCurve(vector<Point>& points,
			     shared_ptr<ParamCurve> curve, double tol,
			     double& maxdist, double& avdist, int& num_inside)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = (int)points.size();
  double fac = 1.0/(double)num;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      curve->closestPoint(points[ki], curve->startparam(), curve->endparam(),
			  tpar, close, dist);
      maxdist = std::max(maxdist, dist);
      avdist += fac*dist;
      if (dist <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
    }
}


//===========================================================================
void RevEngUtils::distToCurve(vector<Point>& points,
			      shared_ptr<ParamCurve> curve, double tol,
			      double& maxdist, double& avdist, int& num_inside,
			      vector<double>& dist)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = 0;
  dist.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar;
      Point close;
      curve->closestPoint(points[ki], curve->startparam(), curve->endparam(),
			  tpar, close, dist[ki]);
      maxdist = std::max(maxdist, dist[ki]);
      avdist += dist[ki];
      if (dist[ki] <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
      ++num;
    }
  avdist /= (double)num;
}

//===========================================================================
void RevEngUtils::distToCurve(vector<Point>& points,
			      shared_ptr<ParamCurve> curve, double tol,
			      double& maxdist, double& avdist, int& num_inside,
			      vector<double>& parvals, vector<double>& dist)
//===========================================================================
{
  double eps = 1.0e-6;
  maxdist = avdist = 0.0;
  num_inside = 0;
  int num = 0;
  parvals.resize(points.size());
  dist.resize(points.size());
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar;
      Point close;
      curve->closestPoint(points[ki], curve->startparam(), curve->endparam(),
			  parvals[ki], close, dist[ki]);
      maxdist = std::max(maxdist, dist[ki]);
      avdist += dist[ki];
      if (dist[ki] <= tol)
	++num_inside;
      else
	{
	  int stop_break = 1;
	}
      ++num;
    }
  avdist /= (double)num;
}


//===========================================================================
shared_ptr<ElementarySurface>
RevEngUtils::elemsurfWithAxis(shared_ptr<ElementarySurface> sf_in,
			      vector<RevEngPoint*>& points,
			      Point mainaxis[3], double diag)
//===========================================================================
{
  Point axis = sf_in->direction();
  Point loc = sf_in->location();
  Point Cx = sf_in->direction2();
  if (sf_in->instanceType() == Class_Plane)
    return planeWithAxis(points, axis, loc, mainaxis);
  else if (sf_in->instanceType() == Class_Cylinder)
    return cylinderWithAxis(points, axis, Cx, loc);
  else if (sf_in->instanceType() == Class_Torus)
    return torusWithAxis(points, axis, Cx, loc);
  else if (sf_in->instanceType() == Class_Sphere)
    return sphereWithAxis(points, axis, Cx, loc);
  else if (sf_in->instanceType() == Class_Cone)
    return coneWithAxis(points, axis, Cx, loc, diag);
  else
    {
      shared_ptr<ElementarySurface> dummy;
      return dummy;
    }
}

//===========================================================================
shared_ptr<Plane> RevEngUtils::planeWithAxis(vector<RevEngPoint*>& points,
					     Point axis, Point init_loc,
					     Point mainaxis[3])
//===========================================================================
{
  Point mid(0.0, 0.0, 0.0);
  double wgt = 1.0/(double)points.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D pos0 = points[ki]->getPoint();
      Point pos(pos0[0], pos0[1], pos0[2]);
      Point vec = pos - init_loc;
      Point pos2 = init_loc + (vec*axis)*axis;
      mid += wgt*pos2;
    }

  int ix=-1;
  double min_ang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = axis.angle(mainaxis[ka]);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  ix = ka;
	}
    }
  shared_ptr<Plane> plane(new Plane(mid, axis, mainaxis[(ix+1)%3]));
  return plane;
}

//===========================================================================
shared_ptr<Cylinder> RevEngUtils::cylinderWithAxis(vector<RevEngPoint*>& points,
						   Point axis, Point low, 
						   Point high, Point mainaxis[3])
//===========================================================================
{
  int ix=-1;
  double min_ang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = axis.angle(mainaxis[ka]);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  ix = ka;
	}
    }

  Point Cx = mainaxis[(ix+2)%3].cross(axis);
  Cx.normalize();
  Point Cy = axis.cross(Cx);
  Cy.normalize();

  Point pos;
  double rad;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::computeCylPosRadius(group, low, high, axis, Cx, Cy, pos, rad);

  shared_ptr<Cylinder> cyl(new Cylinder(rad, pos, axis, Cy));
  return cyl;
}

//===========================================================================
shared_ptr<Cylinder> RevEngUtils::cylinderWithAxis(vector<RevEngPoint*>& points,
						   Point axis, Point Cx,
						   Point pos)
//===========================================================================
{
  // Compute radius
  double rad = 0.0;
  double fac = 1.0/(double)points.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      Point pnt2 = pos + ((pnt - pos)*axis)*axis;
      double dd = pnt.dist(pnt2);
      rad += fac*dd;
    }

  shared_ptr<Cylinder> cyl(new Cylinder(rad, pos, axis, Cx));
  return cyl;
}

//===========================================================================
shared_ptr<Torus> RevEngUtils::torusWithAxis(vector<RevEngPoint*>& points,
					     Point axis, Point loc, 
					     Point mainaxis[3])
//===========================================================================
{
  int ix=-1;
  double min_ang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = axis.angle(mainaxis[ka]);
      ang = std::min(ang, M_PI-ang);
      if (ang < min_ang)
	{
	  min_ang = ang;
	  ix = ka;
	}
    }

  Point Cx = mainaxis[(ix+2)%3].cross(axis);
  Cx.normalize();
  Point Cy = axis.cross(Cx);
  Cy.normalize();

  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  rotateToPlane(group, Cy, axis, loc, rotated);
  
  // Approximate rotated points with a circle
  Point centre;
  double radius;
  computeCircPosRadius(rotated, Cx, Cy, axis, centre, radius);

  Point axis_pt = loc + ((centre - loc)*axis)*axis;
  double dist = centre.dist(axis_pt);
  shared_ptr<Torus> torus(new Torus(dist, radius, axis_pt, axis, Cx));

  return torus;
 }

//===========================================================================
shared_ptr<Torus> RevEngUtils::torusWithAxis(vector<RevEngPoint*>& points,
					     Point axis, Point Cx, Point pos)
//===========================================================================
{
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  Point Cy = axis.cross(Cx);
  rotateToPlane(group, Cy, axis, pos, rotated);
  
  // Approximate rotated points with a circle
  Point centre;
  double radius;
  computeCircPosRadius(rotated, Cx, Cy, axis, centre, radius);

  Point axis_pt = pos + ((centre - pos)*axis)*axis;
  double dist = centre.dist(axis_pt);
  shared_ptr<Torus> torus(new Torus(dist, radius, axis_pt, axis, Cx));

  return torus;
 }

//===========================================================================
shared_ptr<Sphere> RevEngUtils::sphereWithAxis(vector<RevEngPoint*>& points,
					       Point axis, 
					       Point mainaxis[3])
//===========================================================================
{
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  Point centre;
  double radius;
  try {
    computeSphereProp(group, centre, radius);
  }
  catch (...)
    {
      shared_ptr<Sphere> dummy;
      return dummy;
    }

  // Define x-axis
  int ix = -1;
  double minang = M_PI;
  for (int ka=0; ka<3; ++ka)
    {
      double ang = mainaxis[ka].angle(axis);
      ang = std::min(ang, M_PI-ang);
      if (ang < minang)
	{
	  minang = ang;
	  ix = ka;
	}
    }

  Point Cy = mainaxis[(ix+1)%3].cross(axis);
  Point Cx = axis.cross(Cy);
  
  shared_ptr<Sphere> sph(new Sphere(radius, centre, axis, Cx));
  return sph;
}

//===========================================================================
shared_ptr<Sphere> RevEngUtils::sphereWithAxis(vector<RevEngPoint*>& points,
						   Point axis, Point Cx,
						   Point pos)
//===========================================================================
{
  // Compute radius
  double rad = 0.0;
  double fac = 1.0/(double)points.size();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      Point pnt(xyz[0], xyz[1], xyz[2]);
      double dd = pnt.dist(pos);
      rad += fac*dd;
    }
  
  shared_ptr<Sphere> sph(new Sphere(rad, pos, axis, Cx));
  return sph;
}

//===========================================================================
shared_ptr<Cone> RevEngUtils::coneWithAxis(vector<RevEngPoint*>& points,
					   Point axis, Point low, 
					   Point high, Point mainaxis[3])
//===========================================================================
{
  shared_ptr<Cylinder> cyl =
    RevEngUtils::cylinderWithAxis(points, axis, low, high, mainaxis);
  Point pnt = cyl->location();
  Point Cx, Cy, Cz;
  cyl->getCoordinateAxes(Cx, Cy, Cz);
  //double rad = cyl->getRadius();

#ifdef DEBUG_CONE
  std::ofstream of("pts_cone.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << points.size() << std::endl;
  for (size_t ki=0; ki<points.size(); ++ki)
    of << points[ki]->getPoint() << std::endl;
#endif
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::rotateToPlane(group, Cy, axis, pnt, rotated);

  double len = low.dist(high);
  shared_ptr<Line> line(new Line(pnt, axis));
  line->setParameterInterval(-len, len);
  
  Point pt1 = line->ParamCurve::point(-len);
  Point pt2 = line->ParamCurve::point(len);
  shared_ptr<SplineCurve> line_cv(new SplineCurve(pt1, -len, pt2, len));
  shared_ptr<SplineCurve> cv1;
  curveApprox(rotated, line_cv, 2, 2, cv1);

  double tclose, dclose;
  Point ptclose;
  cv1->closestPoint(pnt, cv1->startparam(), cv1->endparam(), tclose,
		    ptclose, dclose);
  vector<Point> der(2);
  cv1->point(der, tclose, 1);
  double phi = der[1].angle(axis);
  shared_ptr<Cone> cone(new Cone(dclose, pnt, axis, Cy, phi));

#ifdef DEBUG_CONE
  cone->writeStandardHeader(of);
  cone->write(of);
#endif
  return cone;
}

//===========================================================================
shared_ptr<Cone> RevEngUtils::coneWithAxis(vector<RevEngPoint*>& points,
					   Point axis, Point Cx, Point pos,
					   double len)
//===========================================================================
{
#ifdef DEBUG_CONE
  std::ofstream of("pts_cone.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << points.size() << std::endl;
  for (size_t ki=0; ki<points.size(); ++ki)
    of << points[ki]->getPoint() << std::endl;
#endif
  
  vector<Point> rotated;
  vector<pair<vector<RevEngPoint*>::iterator,
	      vector<RevEngPoint*>::iterator> > group;
  group.push_back(std::make_pair(points.begin(), points.end()));
  RevEngUtils::rotateToPlane(group, Cx, axis, pos, rotated);

  shared_ptr<Line> line(new Line(pos, axis));
  line->setParameterInterval(-len, len);
  
  Point pt1 = line->ParamCurve::point(-len);
  Point pt2 = line->ParamCurve::point(len);
  shared_ptr<SplineCurve> line_cv(new SplineCurve(pt1, -len, pt2, len));
  shared_ptr<SplineCurve> cv1;
  curveApprox(rotated, line_cv, 2, 2, cv1);
#ifdef DEBUG_CONE
  cv1->writeStandardHeader(of);
  cv1->write(of);
#endif

  // Expects one intersection point with plane through point-on-axis
  shared_ptr<Cone> cone;
  double eps = 1.0e-6;
  vector<double> intpar;
  vector<pair<double,double> > intcvs;
  double tpar;
  intersectCurvePlane(cv1.get(), pos, axis, eps, intpar, intcvs);
  if (intpar.size() > 0)
    tpar = intpar[0];
  else if (intcvs.size() > 0)
    tpar = 0.5*(intcvs[0].first + intcvs[0].second);
  else
    return cone;

  vector<Point> der(2);
  cv1->point(der, tpar, 1);
  double phi = der[1].angle(axis);
  double rad = der[0].dist(pos);
  cone = shared_ptr<Cone>(new Cone(rad, pos, axis, Cx, phi));

#ifdef DEBUG_CONE
  if (cone.get())
    {
      cone->writeStandardHeader(of);
      cone->write(of);
    }
#endif
  return cone;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergePlanes(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						    vector<int>& nmbpts,
						    bool set_bound)
//===========================================================================
{
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  Point pos(0.0, 0.0, 0.0);
  Point norm(0.0, 0.0, 0.0);
  vector<RevEngPoint*> all_pts;
  all_pts.reserve(totnmb);
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double wgt = 1.0/(double)totnmb;

      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  Point curr = (*it)->getLocFuncNormal();
	  Vector3D xyz = (*it)->getPoint();
	  pos +=  wgt*Point(xyz[0], xyz[1], xyz[2]);
	  norm += wgt*curr;
	  all_pts.push_back(*it);
	}
    }
  
  // Perform approximation with combined point set
  ImplicitApprox impl;
  impl.approx(points, 1);
  Point pos2, normal2;
  bool found = impl.projectPoint(pos, norm, pos2, normal2);
  if (!found)
    {
      pos2 = pos;
      normal2 = norm;
      normal2.normalize();
    }
   // std::ofstream outviz("implsf_merge.g2");
  // impl->visualize(all_pts, outviz);
 
  shared_ptr<Plane> surf(new Plane(pos2, normal2));
  Point low = bbox.low();
  Point high = bbox.high();
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParameterBounds(-0.5*len, -0.5*len, 0.5*len, 0.5*len);
    }

  return surf;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeCylinders(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						    vector<int>& nmbpts,
						    bool set_bound)
//===========================================================================
{
  // Estimate cylinder axis
  Point axis, Cx, Cy;
  RevEngUtils::computeAxis(points, axis, Cx, Cy);

  // Estimate radius and point on axis
  double rad;
  Point pnt;
  Point low = bbox.low();
  Point high = bbox.high();
  RevEngUtils::computeCylPosRadius(points, low, high,
				   axis, Cx, Cy, pnt, rad);
  shared_ptr<Cylinder> surf(new Cylinder(rad, pnt, axis, Cy));
  if (set_bound)
    {
      double len = low.dist(high);
      surf->setParamBoundsV(-len, len);
    }

  return surf;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeSpheres(vector<pair<vector<RevEngPoint*>::iterator,
						    vector<RevEngPoint*>::iterator> > points,
						    const BoundingBox& bbox,
						     vector<int>& nmbpts, Point& normal)
//===========================================================================
{
  Point centre;
  double radius;
  try {
    RevEngUtils::computeSphereProp(points, centre, radius);
  }
  catch (...)
    {
      shared_ptr<Sphere> dummy;
      return dummy;
    }

  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  vector<Point> pnts;
  pnts.reserve(totnmb);
  for (size_t ki=0; ki<points.size(); ++ki)
    for (auto it=points[ki].first; it!=points[ki].second; ++it)
      {
	Vector3D xyz = (*it)->getPoint();
	pnts.push_back(Point(xyz[0], xyz[1], xyz[2]));
      }

  double eigenvec[3][3];
  double lambda[3];
  Point eigen1, eigen2, eigen3;
  RevEngUtils::principalAnalysis(pnts[0], pnts, lambda, eigenvec);
  Point z_axis = Point(eigenvec[1][0], eigenvec[1][1], eigenvec[1][2]);
  Point x_axis = normal.cross(z_axis);

  shared_ptr<Sphere> sph(new Sphere(radius, centre, z_axis, normal));

  return sph;
}

//===========================================================================
shared_ptr<ParamSurface> RevEngUtils::doMergeTorus(vector<pair<vector<RevEngPoint*>::iterator,
						   vector<RevEngPoint*>::iterator> > points,
						   const BoundingBox& bbox,
						   vector<int>& nmbpts)
//===========================================================================
{
  shared_ptr<ParamSurface> dummy;
  
  int totnmb = 0;
  for (size_t kh=0; kh<nmbpts.size(); ++kh)
    totnmb += nmbpts[kh];
  
  // Compute mean curvature and initial point in plane
  double k2mean = 0.0;
  double wgt = 1.0/(double)totnmb;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it)
	{
	  double kmax = (*it)->maxPrincipalCurvature();
	  k2mean += wgt*kmax;
	}
    }
  double rd = 1.0/k2mean;
  
  vector<Point> centr(totnmb);
  Point mid(0.0, 0.0, 0.0);
  size_t kr = 0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      for (auto it=points[ki].first; it!=points[ki].second; ++it, ++kr)
	{
	  Point norm = (*it)->getLocFuncNormal();
	  Vector3D xyz = (*it)->getPoint();
	  Point xyz2(xyz[0], xyz[1], xyz[2]);
	  centr[kr] = xyz2 + rd*norm;
	  mid += wgt*centr[kr];
	}
    }
  
  ImplicitApprox impl;
  impl.approxPoints(centr, 1);

  double val;
  Point grad;
  impl.evaluate(mid, val, grad);
  grad.normalize_checked();
  Point pos, normal;
  bool found = impl.projectPoint(mid, grad, pos, normal);
  if (!found)
    return dummy;
  double eps1 = 1.0e-8;
  if (normal.length() < eps1)
    return dummy;
  
  Point Cx = centr[0] - mid;
  Cx -= (Cx*normal)*normal;
  Cx.normalize();
  Point Cy = Cx.cross(normal);
  
  double rad;
  Point pnt;
  RevEngUtils::computeCircPosRadius(centr, normal, Cx, Cy, pnt, rad);
  pnt -= ((pnt - pos)*normal)*normal;

  vector<Point> rotated;
  RevEngUtils::rotateToPlane(points, Cx, normal, pnt, rotated);
  Point cpos;
  double crad;
  RevEngUtils::computeCircPosRadius(rotated, Cy, Cx, normal, cpos, crad);
  Point cvec = cpos - pnt;
  double R1 = (cvec - (cvec*normal)*normal).length();
  double R2 = (cvec*normal)*normal.length();
 
  shared_ptr<Torus> surf(new Torus(R1, crad, pnt+R2*normal, normal, Cy));

  return surf;
}


//===========================================================================
shared_ptr<SplineCurve> RevEngUtils::midCurve(shared_ptr<SplineCurve>& cv1,
					      shared_ptr<SplineCurve>& cv2)
//===========================================================================
{
  shared_ptr<SplineCurve> spl1(cv1->clone());
  shared_ptr<SplineCurve> spl2(cv2->clone());

  // Check orientation
  Point pt1 = spl1->ParamCurve::point(spl1->startparam());
  Point pt2 = spl1->ParamCurve::point(spl1->endparam());
  Point pt3 = spl2->ParamCurve::point(spl2->startparam());
  Point pt4 = spl2->ParamCurve::point(spl2->endparam());
  double len1 = pt1.dist(pt3);
  double len2 = pt1.dist(pt4);
  if (len2 < len1)
    spl2->reverseParameterDirection();

  // Ensure same spline room
  spl2->setParameterInterval(spl1->startparam(), spl1->endparam());

  double tol = 1.0e-4;
  vector<shared_ptr<SplineCurve> > curves(2);
  curves[0] = spl1;
  curves[1] = spl2;
  GeometryTools::unifyCurveSplineSpace(curves, tol);

  shared_ptr<SplineCurve> midcv = GeometryTools::curveSum(*curves[0], 0.5,
							  *curves[1], 0.5);
  return midcv;
}


//===========================================================================
void  RevEngUtils::curveApprox(vector<Point>& points,
			       shared_ptr<ParamCurve> cvin,
			       int ik, int in, 
			       shared_ptr<SplineCurve>& curve)
//===========================================================================
{
  vector<double> pts;
  vector<double> param;
  double tmin = cvin->startparam();
  double tmax = cvin->endparam();
  double tmin2 = tmax;
  double tmax2 = tmin;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      cvin->closestPoint(points[ki], tmin, tmax, tpar, close, dist);
      pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      param.push_back(tpar);
      tmin2 = std::min(tmin2, tpar);
      tmax2 = std::max(tmax2, tpar);
    }

  double tdel = (tmax2 - tmin2)/(double)(in - ik + 1);
  vector<double> et(ik+in);
  for (int ka=0; ka<ik; ++ka)
    {
      et[ka] = tmin2;
      et[in+ka] = tmax2;
    }
  for (int ka=ik; ka<in; ++ka)
    et[ka] = tmin2 + (ka-ik+1)*tdel;

  vector<double> ecoef(3*in, 0.0);
  shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, &et[0], &ecoef[0], 3));

  SmoothCurve smooth(3);
  vector<int> cfn(in, 0);
  vector<double> wgts(param.size(), 1.0);
  smooth.attach(cv, &cfn[0]);

  smooth.setOptim(0.0, 0.001, 0.001);
  smooth.setLeastSquares(pts, param, wgts, 0.998);

  smooth.equationSolve(curve);
  int stop_break = 1;
}

//===========================================================================
void  RevEngUtils::curveApprox(vector<Point>& points,
			       vector<double>& param,
			       int ik, int in, 
			       shared_ptr<SplineCurve>& curve)
//===========================================================================
{
  if (points.size() == 0)
    return;
  if (points.size() != param.size())
    return;
  vector<double> pts;
  double tmin = std::numeric_limits<double>::max();
  double tmax = std::numeric_limits<double>::lowest();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      pts.insert(pts.end(), points[ki].begin(), points[ki].end());
      tmin = std::min(tmin, param[ki]);
      tmax = std::max(tmax, param[ki]);
    }

  double tdel = (tmax - tmin)/(double)(in - ik + 1);
  vector<double> et(ik+in);
  for (int ka=0; ka<ik; ++ka)
    {
      et[ka] = tmin;
      et[in+ka] = tmax;
    }
  for (int ka=ik; ka<in; ++ka)
    et[ka] = tmin + (ka-ik+1)*tdel;

  vector<double> ecoef(3*in, 0.0);
  shared_ptr<SplineCurve> cv(new SplineCurve(in, ik, &et[0], &ecoef[0], 3));

  SmoothCurve smooth(3);
  vector<int> cfn(in, 0);
  vector<double> wgts(param.size(), 1.0);
  smooth.attach(cv, &cfn[0]);

  smooth.setOptim(0.0, 0.001, 0.001);
  smooth.setLeastSquares(pts, param, wgts, 0.998);

  smooth.equationSolve(curve);
  int stop_break = 1;
}

//===========================================================================
shared_ptr<SplineCurve> RevEngUtils::createCurve(vector<RevEngPoint*>& points,
						 int degree, double tol, int maxiter)
//===========================================================================
{
  shared_ptr<SplineCurve> cv;
  if (points.size() < 2)
    return cv;
  
  // Parameterize curves and fetch data points
  vector<double> param(points.size(), 0.0);
  vector<double> pts;
  pts.reserve(3*points.size());
  Vector3D prev = points[0]->getPoint();
  double tmp = 0.0;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      Vector3D xyz = points[ki]->getPoint();
      pts.insert(pts.end(), xyz.begin(), xyz.end());
      param[ki] = tmp + prev.dist(xyz);
      prev = xyz;
      tmp = param[ki];
    }

  double smoothwgt = 0.1;
  ApproxCurve approx(pts, param, 3, tol, degree+1, degree+1);
  approx.setSmooth(smoothwgt);

  double maxdist, avdist;
  cv = approx.getApproxCurve(maxdist, avdist, maxiter);

  return cv;
}

//===========================================================================
void RevEngUtils::extractLinearPoints(vector<RevEngPoint*>& points,
				      vector<Point>& rotated, double len,
				      Point& pos, Point& axis, double rad,
				      Point& axis2, bool plane,
				      double tol, double angtol,
				      vector<pair<double,double> >& dist_ang,
				      vector<RevEngPoint*>& linear, bool start,
				      vector<RevEngPoint*>& remaining)
//===========================================================================
{
  // Parametarize the points according to axis
  vector<double> param(points.size());
  vector<double> distance(points.size());
  vector<int> perm(points.size());
  shared_ptr<Line> line(new Line(pos, axis));
  line->setParameterInterval(-len, len);

  double tmin = -len;
  double tmax = len;
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      line->closestPoint(rotated[ki], -len, len, tpar, close, dist);
      param[ki] = tpar;
      distance[ki] = dist;
      perm[ki] = (int)ki;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
    }

  // Sort
  for (size_t ki=0; ki<perm.size(); ++ki)
    for (size_t kj=ki+1; kj<perm.size(); ++kj)
      if (param[perm[kj]] < param[perm[ki]])
	std::swap(perm[ki], perm[kj]);

  // Identify linear points
  double pihalf = 0.5*M_PI;
  int ix1 = (start) ? 0 : (int)perm.size()-1;
  int ix2 = (start) ? (int)perm.size() : -1;
  int sgn = (start) ? 1 : -1;
  int lim1 = std::min(10, (int)points.size()/200);
  int lim2 = std::min(100, (int)points.size()/50);
  int ka, kb;
  int num_out = 0;
  double fac = 0.9;
  double dfac = 0.5;
  for (ka=ix1; ka!=ix2; ka=kb)
    {
      Point norm = points[perm[ka]]->getLocFuncNormal();
      Point norm2 = points[perm[ka]]->getTriangNormal();
      double ang = norm.angle(axis2);
      double ang2 = norm2.angle(axis2);
      if (plane)
	ang = std::min(std::min(ang,M_PI-ang), std::min(ang2,M_PI-ang2));
      else
	ang = std::min(fabs(pihalf-ang), fabs(pihalf-ang2));
      double dd = fabs(distance[perm[ka]]-rad);
      double anglim = (dist_ang.size() > 0) ? dist_ang[perm[ka]].second : 0.0;
      double angtol2 = std::max(fac*anglim, angtol);
      double tol2 = (dist_ang.size() > 0) ?
	std::max(tol, dfac*dist_ang[perm[ka]].first) : tol;;
      if (dd <= tol2 && ang <= angtol2)
	kb = ka+sgn;
      else
	{
	  num_out++;
	  for (kb=ka+sgn; kb!=ix2; kb+=sgn)
	    {
	      norm = points[perm[kb]]->getLocFuncNormal();
	      norm2 = points[perm[kb]]->getTriangNormal();
	      ang = norm.angle(axis2);
	       ang2 = norm2.angle(axis2);
	      ang = std::min(fabs(pihalf-ang), fabs(pihalf-ang2));
	      dd = fabs(distance[perm[kb]]-rad);
	      anglim = (dist_ang.size() > 0) ? dist_ang[perm[kb]].second : 0.0;
	      angtol2 = std::max(fac*anglim, angtol);
	      tol2 = (dist_ang.size() > 0) ?
		std::max(tol, dfac*dist_ang[perm[kb]].first) : tol;;
	      if (dd <= tol2 && ang <= angtol2)
		break;
	      num_out++;
	      if (abs(kb-ka) > lim1 || num_out > lim2)
		break;
	    }
	}
      if (abs(kb-ka) > lim1 || num_out > lim2)
	break;
    }

  // Extract linear points
  if (ka == ix2)
    ka+=sgn;
  for (kb=ix1; kb!=ka && kb<(int)perm.size() && kb>=0; kb+=sgn)
    linear.push_back(points[perm[kb]]);
  for (; kb!=ix2 && kb<(int)perm.size() && kb>=0; kb+=sgn)
    remaining.push_back(points[perm[kb]]);
  int stop_break = 1;
}

struct LinInfo
{
  double t1_, t2_;
  double mind_, maxd_, avd_, avd2_;
  int nmb_;

  LinInfo()
  {
    t1_ = t2_ = mind_ = maxd_ = avd_ = avd2_ = 0.0;
    nmb_ = 0;
  }

  LinInfo(double t1, double t2, int nmb, double mind, double maxd,
	  double avd, double avd2)
  {
    t1_ = t1;
    t2_ = t2;
    nmb_ = nmb;
    mind_ = mind;
    maxd_ = maxd;
    avd_ = avd;
    avd2_ = avd2;
  }
};

//===========================================================================
bool RevEngUtils::extractLinearPoints(vector<Point>& points,
				      shared_ptr<Line>& line, Point norm,
				      double tol, bool start, double& splitpar,
				      vector<Point>& linear, 
				      vector<Point>& remaining)
//===========================================================================
{
  // Parametarize the points according to axis
  vector<double> param(points.size());
  vector<double> distance(points.size());
  vector<int> perm(points.size());

  double tmin = line->startparam();
  double tmax = line->endparam();
  for (size_t ki=0; ki<points.size(); ++ki)
    {
      double tpar, dist;
      Point close;
      line->closestPoint(points[ki], line->startparam(), line->endparam(),
			 tpar, close, dist);
      param[ki] = tpar;
      int sgn = ((close-points[ki])*norm > 0.0) ? 1 : -1;
      distance[ki] = sgn*dist;
      perm[ki] = (int)ki;
      tmin = std::min(tmin, tpar);
      tmax = std::max(tmax, tpar);
    }

  // Sort
  for (size_t ki=0; ki<perm.size(); ++ki)
    for (size_t kj=ki+1; kj<perm.size(); ++kj)
      if (param[perm[kj]] < param[perm[ki]])
	std::swap(perm[ki], perm[kj]);

  int ndiv = 100;
  double tdel = (tmax - tmin)/(double)ndiv;
  vector<LinInfo> info(ndiv);
  double t1, t2;
  int ka;
  size_t kj=0;
  for (ka=0, t1=tmin, t2=t1+tdel; ka<ndiv; ++ka, t1+=tdel, t2+=tdel)
    {
      int nmb = 0;
      double maxd = std::numeric_limits<double>::lowest();
      double mind = std::numeric_limits<double>::max();
      double avd = 0.0, avd2 = 0.0;
      for (; kj<perm.size() && param[perm[kj]] <= t2; ++kj)
	{
	  avd += fabs(distance[perm[kj]]);
	  avd2 += distance[perm[kj]];
	  maxd = std::max(maxd, distance[perm[kj]]);
	  mind = std::min(mind, distance[perm[kj]]);
	  ++nmb;
	}
      if (nmb > 0)
	{
	  avd /= (double)nmb;
	  avd2 /= (double)nmb;
	}
      else
	{
	  mind = maxd = avd = avd2 = 0.0;
	}
      LinInfo curr_info(t1, t2, nmb, mind, maxd, avd, avd2);
      info[ka] = curr_info;
    }

  // Count number of information entities without any points
  int num_zero = 0;
  for (size_t ki=0; ki<info.size(); ++ki)
    if (info[ki].nmb_ == 0)
      num_zero++;

  int zero_lim = ndiv/3;
  if (num_zero > zero_lim)
    return false;
  
  vector<double> par;
  vector<double> range;
  vector<double> avdist1;
  vector<double> avdist2;
  for (int ka=0; ka<ndiv; ++ka)
    {
      if (info[ka].nmb_ > 0)
	{
	  par.push_back(0.5*(info[ka].t1_ + info[ka].t2_));
	  range.push_back(info[ka].maxd_ - info[ka].mind_);
	  avdist1.push_back(info[ka].avd_);
	  avdist2.push_back(info[ka].avd2_);
	}
    }

  int in = 6;
  int ik = 4;
  double smoothwgt = 1.0e-5;
  int maxiter = 4;
  double tol2 = 1.0e-6;
  double maxdist, avdist;
  ApproxCurve approx2(avdist1, par, 1, tol2, in, ik);
  shared_ptr<SplineCurve> cv2 = approx2.getApproxCurve(maxdist, avdist, maxiter);
#ifdef DEBUG_BLEND
  double maxdistd, avdistd;
  ApproxCurve approx1(range, par, 1, tol2, in, ik);
  shared_ptr<SplineCurve> cv1 = approx1.getApproxCurve(maxdistd, avdistd, maxiter);
  ApproxCurve approx3(avdist2, par, 1, tol2, in, ik);
  shared_ptr<SplineCurve> cv3 = approx3.getApproxCurve(maxdistd, avdistd, maxiter);
  std::ofstream of("dist_cvs.g2");
  SplineDebugUtils::writeSpace1DCurve(*cv1, of);
  SplineDebugUtils::writeSpace1DCurve(*cv2, of);
  SplineDebugUtils::writeSpace1DCurve(*cv3, of);
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " 0 0 " << tmax << " 0 0 " << std::endl;
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " " << tol << " 0 " << tmax << " " << tol << " 0 " << std::endl;
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " " << 2*tol << " 0 " << tmax << " " << 2*tol << " 0 " << std::endl;
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " " << -tol << " 0 " << tmax << " " << -tol << " 0 " << std::endl;

  shared_ptr<SplineCurve> cv1_2(cv1->derivCurve(1));
  shared_ptr<SplineCurve> cv2_2(cv2->derivCurve(1));
  shared_ptr<SplineCurve> cv3_2(cv3->derivCurve(1));

  Point zero(1);
  zero[0] = 0.0;
  vector<double> intpar1, intpar2, intpar3;
  vector<pair<double,double> > intcv1, intcv2, intcv3;
  intersectCurvePoint(cv1_2.get(), zero, tol2, intpar1, intcv1);
  intersectCurvePoint(cv2_2.get(), zero, tol2, intpar2, intcv2);
  intersectCurvePoint(cv3_2.get(), zero, tol2, intpar3, intcv3);
 for (auto it=cv1_2->coefs_begin(); it != cv1_2->coefs_end(); ++it)
    (*it) *= 0.01;
  for (auto it=cv2_2->coefs_begin(); it != cv2_2->coefs_end(); ++it)
    (*it) *= 0.01;
  for (auto it=cv3_2->coefs_begin(); it != cv3_2->coefs_end(); ++it)
    (*it) *= 0.01;
  std::ofstream of2("dist_cvs_2.g2");
  SplineDebugUtils::writeSpace1DCurve(*cv1_2, of2);
  SplineDebugUtils::writeSpace1DCurve(*cv2_2, of2);
  SplineDebugUtils::writeSpace1DCurve(*cv3_2, of2);
#endif
  // Set limitation of average distance for the line to be an accepted
  // approximation
  int asgn = (start) ? 1 : -1;
  double av1 = 0.0, minav2=2.0*tol, maxav2=0.0;
  int an = 0;
  double afac = 2.0;
  for (int ka=(start)?0:ndiv-1; ka!=ndiv/2; ka+=asgn)
    {
      double av2=0.0;
      double fac = 0.1;
      for (int kb=0; kb<10; ++kb)
	av2 += fac*info[ka+asgn*kb].avd_;

      if (an>0 && av2>afac*av1/(double)an)
	break;
      av1 += av2;
      ++an;
      minav2 = std::min(minav2, av2);
      maxav2 = std::max(maxav2, av2);
    }
  if (an > 0)
    av1 /= (double)an;
  double av2del = maxav2 - minav2;
  av2del = std::max(av2del, avdist);

  splitpar = (start) ? tmin - tmax : 2*tmax;
  vector<double> intpar, intpar_n;
  vector<pair<double,double> > intcv, intcv_n;
  Point tolpt(1);
  double tolfac = 3;
  double toldel = tolfac*av2del;
  toldel = std::max(toldel, 0.25*tol);
  tolpt[0] = std::min(av1 + toldel, tol);
#ifdef DEBUG_BLEND
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " " << tolpt << " 0 " << tmax << " " << tolpt << " 0 " << std::endl;
  of << "410 1 0 0" << std::endl;
  of << "1" << std::endl;
  of << tmin << " " << av1-tolfac*av2del << " 0 " << tmax << " " << av1-tolfac*av2del << " 0 " << std::endl;
#endif
  if (an > 0)
    {
      intersectCurvePoint(cv2.get(), tolpt, tol2, intpar, intcv);
      tolpt[0] = av1 - std::min(toldel, tol-av1);
      intersectCurvePoint(cv2.get(), tolpt, tol2, intpar_n, intcv_n);
      if (intpar_n.size() > 0)
	intpar.insert(intpar.end(), intpar_n.begin(), intpar_n.end());
    }
      
  if (intpar.size() > 0)
    {
      std::sort(intpar.begin(), intpar.end());
      splitpar = (start) ? intpar[0] : intpar[intpar.size()-1];
    }
#ifdef DEBUG_BLEND  
  else
    {
  //int nmbsplit = 0;
  int splitsgn = (start) ? 1 : -1;
  vector<Point> der(2);
  for (size_t kr=0; kr<intpar1.size(); ++kr)
    {
      cv1_2->point(der, intpar1[kr], 1);
      if (splitsgn*der[1][0] > 0.0)
	{
	  // splitpar += intpar1[kr];
	  // nmbsplit++;
	  splitpar = (start) ? std::max(splitpar, intpar1[kr]) :
	    std::min(splitpar, intpar1[kr]);
	}
    }
  for (size_t kr=0; kr<intpar2.size(); ++kr)
    {
      cv2_2->point(der, intpar2[kr], 1);
      if (splitsgn*der[1][0] > 0.0)
	{
	  // splitpar += intpar2[kr];
	  // nmbsplit++;
	  splitpar = (start) ? std::max(splitpar, intpar2[kr]) :
	    std::min(splitpar, intpar2[kr]);
	}
    }
  for (size_t kr=0; kr<intpar3.size(); ++kr)
    {
      cv3_2->point(der, intpar3[kr], 1);
      if (splitsgn*der[1][0] > 0.0)
	{
	  // splitpar += intpar3[kr];
	  // nmbsplit++;
	  splitpar = (start) ? std::max(splitpar, intpar3[kr]) :
	    std::min(splitpar, intpar3[kr]);
	}
    }
    }
#endif
  if (splitpar > tmin && splitpar < tmax) //nmbsplit > 0)
    {
      //splitpar /= (double)nmbsplit;
      if (start)
	{
	  size_t kr;
	  for (kr=0; kr<perm.size() && param[perm[kr]] <= splitpar; ++kr)
	    linear.push_back(points[perm[kr]]);
	  for (; kr<perm.size(); ++kr)
	    remaining.push_back(points[perm[kr]]);
	}
      else
	{
	  size_t kr;
	  for (kr=0; kr<perm.size() && param[perm[kr]] < splitpar; ++kr)
	    remaining.push_back(points[perm[kr]]);
	  for (; kr<perm.size(); ++kr)
	    linear.push_back(points[perm[kr]]);
	}
    }
  else
    remaining.insert(remaining.end(), points.begin(), points.end());
  
  return true;
}


//===========================================================================
void RevEngUtils::identifyEndPoints(vector<RevEngPoint*> edge_pts, shared_ptr<CurveOnSurface>& sfcv,
				    RevEngPoint*& first_pt, double& t1,
				    RevEngPoint*& last_pt, double& t2)
//===========================================================================
{
  shared_ptr<ParamCurve> pcrv = sfcv->parameterCurve();
  shared_ptr<ParamCurve> spacecrv = sfcv->spaceCurve();
  bool parpref = sfcv->parPref();
  t1 = sfcv->endparam();
  t2 = sfcv->startparam();
  double tpar, dist;
  Point close;
  double t3 = sfcv->startparam();
  double t4 = sfcv->endparam();
   if (pcrv && parpref)
    {
      for (size_t ki=0; ki<edge_pts.size(); ++ki)
	{
	  Vector2D uv = edge_pts[ki]->getPar();
	  Point ppt(uv[0], uv[1]);
	  pcrv->closestPoint(ppt, t3, t4, tpar, close, dist);
	  if (tpar < t1)
	    {
	      t1 = tpar;
	      first_pt = edge_pts[ki];
	    }
	  if (tpar > t2)
	    {
	      t2 = tpar;
	      last_pt = edge_pts[ki];
	    }
	}
    }
  else
    {
      for (size_t ki=0; ki<edge_pts.size(); ++ki)
	{
	  Vector3D xyz = edge_pts[ki]->getPoint();
	  Point ppt(xyz[0], xyz[1], xyz[2]);
	  spacecrv->closestPoint(ppt, t3, t4, tpar, close, dist);
	  if (tpar < t1)
	    {
	      t1 = tpar;
	      first_pt = edge_pts[ki];
	    }
	  if (tpar > t2)
	    {
	      t2 = tpar;
	      last_pt = edge_pts[ki];
	    }
	}
     }
}


//===========================================================================
void RevEngUtils::setLoopSeq(vector<shared_ptr<CurveOnSurface> >& cvs)
//===========================================================================
{
  if (cvs.size() <= 1)
    return;
  
  Point pos1 = cvs[0]->ParamCurve::point(cvs[0]->endparam());
  for (size_t ki=1; ki<cvs.size(); ++ki)
    {
      Point pos2 = cvs[ki]->ParamCurve::point(cvs[ki]->startparam());
      Point pos3 = cvs[ki]->ParamCurve::point(cvs[ki]->endparam());
      double dd2 = pos1.dist(pos2);
      double dd3 = pos1.dist(pos3);
      for (size_t kj=ki+1; kj<cvs.size(); ++kj)
	{
	  Point pos4 = cvs[kj]->ParamCurve::point(cvs[kj]->startparam());
	  Point pos5 = cvs[kj]->ParamCurve::point(cvs[kj]->endparam());
	  double dd4 = pos1.dist(pos4);
	  double dd5 = pos1.dist(pos5);
	  if (std::min(dd4,dd5) < std::min(dd2,dd3))
	    {
	      std::swap(cvs[ki], cvs[kj]);
	      std::swap(dd2, dd4);
	      std::swap(dd3, dd5);
	    }
	}
      if (dd3 < dd2)
	cvs[ki]->reverseParameterDirection();
      pos1 = cvs[ki]->ParamCurve::point(cvs[ki]->endparam());
    }
}


//===========================================================================
void RevEngUtils::identifyConGroups(vector<RevEngPoint*>& init,
				    vector<vector<RevEngPoint*> >& groups)
//===========================================================================
{
  // The points may belong to different regions and have different
  // classifications. Mark to identify
  int mark = 1;
  for (size_t ki=0; ki<init.size(); ++ki)
    init[ki]->setMarkIx(mark);

  for (size_t ki=0; ki<init.size(); ++ki)
    {
      if (init[ki]->visited())
	continue;
      vector<RevEngPoint*> curr_group;
      init[ki]->fetchConnectedMarked(mark, curr_group);
       groups.push_back(curr_group);
     }

  for (size_t ki=0; ki<init.size(); ++ki)
    {
      init[ki]->unsetMarkIx();
      init[ki]->unsetVisited();
    }

  // Sort groups according to size
  for (size_t ki=0; ki<groups.size(); ++ki)
    for (size_t kj=ki+1; kj<groups.size(); ++kj)
      if (groups[kj].size() > groups[ki].size())
	std::swap(groups[ki], groups[kj]);
}


