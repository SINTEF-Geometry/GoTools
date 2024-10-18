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

#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SurfaceOfRevolution.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/SurfaceOfLinearExtrusion.h"

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;
using std::swap;


namespace Go
{

//===========================================================================
CurveLoop SurfaceTools::outerBoundarySfLoop(shared_ptr<ParamSurface> surf,
					    double degenerate_epsilon) 
//===========================================================================
{
  // It is convenient to let boundary loops be described as CurveOnSurface
  // to store as much information as possible. Due to problems with shared_ptr, 
  // this is not possible from within SplineSurface.
  // This function is implemented to get around this problem

  shared_ptr<SplineSurface> spline_sf = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);

  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);

  if (bd_sf.get())
    return bd_sf->outerBoundaryLoop(degenerate_epsilon);
  else 
    {
      // Test for degeneracy.
      bool deg[4];
      if (degenerate_epsilon < 0.0 || surf->dimension() == 1)
	deg[0] = deg[1] = deg[2] = deg[3] = false;
      else
	surf->isDegenerate(deg[0], deg[1], deg[2], deg[3], degenerate_epsilon);

      RectDomain dom = surf->containingDomain();
      vector< shared_ptr< ParamCurve > >  vec;
      int pardir[] = {2, 1, 2, 1};
      int boundary[] = {2, 1, 3, 0};
      double parval[] = {dom.vmin(), dom.umax(), dom.vmax(), dom.umin()};
      if (spline_sf.get())
	{
	  // Spline surface
	  for (int edgenum = 0; edgenum < 4; ++edgenum) {
	    if (!deg[edgenum])
	      {
		// Fetch geometry curve
		shared_ptr<ParamCurve> edgecurve(spline_sf->edgeCurve(edgenum));

		// Construct curve on surface with knowledge about what it is
		shared_ptr<ParamCurve> sfcv = 
		  shared_ptr<ParamCurve>(new CurveOnSurface(surf, edgecurve, 
							    pardir[edgenum], 
							    parval[edgenum], 
							    boundary[edgenum]));
		if (edgenum == 2 || edgenum == 3)
		  sfcv->reverseParameterDirection();
		vec.push_back(sfcv);
	      }
	  }
	}
      else
	{
	  // The boundary loop of non bounded surfaces misses information
	  // about the surface and parameter curves. This information must
	  // be added
	  CurveLoop cv_loop = surf->outerBoundaryLoop(degenerate_epsilon);
	  int nmb_cvs = cv_loop.size();
	  if (nmb_cvs == 0)
	    return cv_loop;

	  shared_ptr<CurveOnSurface> cv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv_loop[0]);
	  if (cv.get())
	    return cv_loop; // Already curve on surface curves
	  
	  // Make new loop with curve-on-surface curves 
	  int ki, kj;
	  for (ki=0, kj=0; ki<nmb_cvs; ++kj)
	    {
	      if (deg[kj])
		continue;
	      shared_ptr<ParamCurve> sfcv = 
		shared_ptr<ParamCurve>(new CurveOnSurface(surf, cv_loop[ki], 
							  pardir[kj],
							  parval[kj],
							  boundary[kj]));
	      vec.push_back(sfcv);
	      ++ki;
	    }
	}

      return CurveLoop(vec, (degenerate_epsilon < 0.0) ? DEFAULT_SPACE_EPSILON :
		       degenerate_epsilon);
    }
}

//===========================================================================
std::vector<CurveLoop> SurfaceTools::allBoundarySfLoops(shared_ptr<ParamSurface> surf,
					  double degenerate_epsilon) 
//===========================================================================
{
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);

  if (bd_sf.get())
    {
      //return bd_sf->absolutelyAllBoundaryLoops();
      return bd_sf->allBoundaryLoops();
    }
  else
    {
      // There is only one boundary loop...
      std::vector<CurveLoop> cvloopvec;
      cvloopvec.push_back(outerBoundarySfLoop(surf, degenerate_epsilon));
      return cvloopvec;
    }
}

//===========================================================================
std::vector<CurveLoop> 
SurfaceTools::absolutelyAllBoundarySfLoops(shared_ptr<ParamSurface> surf,
			     double degenerate_epsilon) 
//===========================================================================
{
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);

  if (bd_sf.get())
    {
      return bd_sf->absolutelyAllBoundaryLoops();
    }
  else
    {
      // There is only one boundary loop...
      // Use a negative degeneracy tolarance to tell that also degenerate
      // boundaries must be included in the loop
      std::vector<CurveLoop> cvloopvec;
      //cvloopvec.push_back(SurfaceTools::outerBoundarySfLoop(surf, 0.0));
      cvloopvec.push_back(SurfaceTools::outerBoundarySfLoop(surf, -0.1));
      return cvloopvec;
    }
}

//===========================================================================
void 
SurfaceTools::iterateCornerPos(Point& vertex, 
		 vector<pair<shared_ptr<ParamSurface>, Point> > sfs,
		 double tol)
//===========================================================================
{
  Point prev, curr;
  curr = vertex;
  double wgt_fac = 10.0;
  double wgt = 1.0;
  double wgt_sum;
  int max_iter = 10;
  int kr = 0;

  // Iterate until the vertex point doesn't move
  do
    {
      prev = curr;
      curr.setValue(0.0);
      wgt_sum = 0.0;

      for (size_t ki=0; ki<sfs.size(); ++ki)
	{
	  double seed[2];
	  double clo_u, clo_v, clo_dist;
	  Point clo_pt;
	  seed[0] = sfs[ki].second[0];
	  seed[1] = sfs[ki].second[1];
	  sfs[ki].first->closestPoint(prev, clo_u, clo_v, clo_pt,
				      clo_dist, 0.001*tol, NULL, seed);

	  // Check if the surface is elementary
	  shared_ptr<ElementarySurface> elemsf = 
	    dynamic_pointer_cast<ElementarySurface, ParamSurface>(sfs[ki].first);
	  double curr_wgt = (elemsf.get()) ? wgt*wgt_fac : wgt;
	  curr += curr_wgt*clo_pt;
	  wgt_sum += curr_wgt;

	  sfs[ki].second = Point(clo_u, clo_v);
	}
      curr /= wgt_sum;
      kr++;
      if (kr > max_iter)
	break;
    }
  while (prev.dist(curr) > tol);

  vertex = curr;
}

//===========================================================================
  bool SurfaceTools::cornerToCornerSfs(shared_ptr<ParamSurface> sf1,
			 shared_ptr<CurveOnSurface> sf_cv1,
			 shared_ptr<ParamSurface> sf2,
			 shared_ptr<CurveOnSurface> sf_cv2,
			 double tol)
//===========================================================================
  {
  int bd1, bd2;  // Specifies the surface boundaries corresponding to 
  // the current edges
  // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax
  bool same_orient1, same_orient2;
  bd1 = sf_cv1->whichBoundary(tol, same_orient1);
  bd2 = sf_cv2->whichBoundary(tol, same_orient2);
  if (bd1 < 0 || bd2 < 0)
    return false;  // Adjacency not along boundary

  // Get surface parameters at corners
  RectDomain dom1 = sf1->containingDomain();
  RectDomain dom2 = sf2->containingDomain();
  double corn1_1[2], corn1_2[2], corn2_1[2], corn2_2[2];
  if (bd1 == 0 || bd1 == 1)
    {
      if (bd1 == 0)
	corn1_1[0] = corn1_2[0] = dom1.umin();
      else
	corn1_1[0] = corn1_2[0] = dom1.umax();
      corn1_1[1] = dom1.vmin();
      corn1_2[1] = dom1.vmax();
    }
  else if (bd1 == 2 || bd1 == 3)
    {
      if (bd1 == 2)
	corn1_1[1] = corn1_2[1] = dom1.vmin();
      else
	corn1_1[1] = corn1_2[1] = dom1.vmax();
      corn1_1[0] = dom1.umin();
      corn1_2[0] = dom1.umax();
    }
  if (bd2 == 0 || bd2 == 1)
    {
      if (bd2 == 0)
	corn2_1[0] = corn2_2[0] = dom2.umin();
      else
	corn2_1[0] = corn2_2[0] = dom2.umax();
      corn2_1[1] = dom2.vmin();
      corn2_2[1] = dom2.vmax();
    }
  else if (bd2 == 2 || bd2 == 3)
    {
      if (bd2 == 2)
	corn2_1[1] = corn2_2[1] = dom2.vmin();
      else
	corn2_1[1] = corn2_2[1] = dom2.vmax();
      corn2_1[0] = dom2.umin();
      corn2_2[0] = dom2.umax();
    }
   
  // Evaluate surface corners
  Point pt1 = sf1->point(corn1_1[0], corn1_1[1]);
  Point pt2 = sf1->point(corn1_2[0], corn1_2[1]);
  Point pt3 = sf2->point(corn2_1[0], corn2_1[1]);
  Point pt4 = sf2->point(corn2_2[0], corn2_2[1]);
  
  if (pt1.dist(pt3) > tol && pt1.dist(pt4) > tol)
    return false;
  if (pt2.dist(pt3) > tol && pt2.dist(pt4) > tol)
    return false;
  if (pt3.dist(pt1) > tol && pt3.dist(pt2) > tol)
    return false;
  if (pt4.dist(pt1) > tol && pt4.dist(pt2) > tol)
    return false;
  
  return true;
  }

//===========================================================================
  bool SurfaceTools::getSfAdjacencyInfo(shared_ptr<ParamSurface> sf1,
			  shared_ptr<CurveOnSurface> sf_cv1,
			  shared_ptr<ParamSurface> sf2,
			  shared_ptr<CurveOnSurface> sf_cv2,
			  double tol,
			  int& bd1, int& bd2, bool& same_orient)
//===========================================================================
  {
    // bd1, bd2:
    // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax
    bool same_orient1, same_orient2;
    bd1 = sf_cv1->whichBoundary(tol, same_orient1);
    bd2 = sf_cv2->whichBoundary(tol, same_orient2);
    if (bd1 < 0 || bd2 < 0)
      return false;  // Adjacency not along boundary

    Point f1_p1 = sf_cv1->faceParameter(sf_cv1->startparam());
    Point f1_p2 = sf_cv1->faceParameter(sf_cv1->endparam());
    Point f2_p1 = sf_cv2->faceParameter(sf_cv2->startparam());
    Point f2_p2 = sf_cv2->faceParameter(sf_cv2->endparam());
    bool opposite = false;
    Point p1 = sf1->ParamSurface::point(f1_p1[0], f1_p1[1]);
    Point p2 = sf1->ParamSurface::point(f1_p2[0], f1_p2[1]);
    Point p3 = sf2->ParamSurface::point(f2_p1[0], f2_p1[1]);
    Point p4 = sf2->ParamSurface::point(f2_p2[0], f2_p2[1]);
    if ((p2 - p1)*(p4 -p3) < 0.0)
      opposite = true;
    if ((same_orient1 && !same_orient2) || (!same_orient1 && same_orient2))
      opposite = !opposite;
    same_orient = !opposite;
    return true;
   }

//===========================================================================
  bool SurfaceTools::getCorrCoefEnum(shared_ptr<SplineSurface> sf1,
		       shared_ptr<SplineSurface> sf2,
		       int bd1, int bd2, bool same_orient,
		       vector<pair<int,int> >& enumeration)
//===========================================================================
  {
    int kn1 = sf1->numCoefs_u();
    int kn2 = sf1->numCoefs_v();
    int kn3 = sf2->numCoefs_u();
    int kn4 = sf2->numCoefs_v();

    int nmb1 = (bd1 == 0 || bd1 == 1) ? kn2 : kn1;
    int nmb2 = (bd2 == 0 || bd2 == 1) ? kn4 : kn3;
    if (nmb1 != nmb2)
      return false;  // No correspondence

    enumeration.resize(nmb1);
    int start1 = (bd1 == 0 || bd1 == 2) ? 0 :
      ((bd1 == 1) ? kn1-1 : kn1*(kn2-1));
    int del1 = (bd1 == 0 || bd1 == 1) ? kn1 : 1;

    int start2 = (bd2 == 0 || bd2 == 2) ? 0 :
      ((bd2 == 1) ? kn3-1 : kn3*(kn4-1));
    int del2 = (bd2 == 0 || bd2 == 1) ? kn3 : 1;
    if (!same_orient)
      {
	start2 += (nmb2-1)*del2;
	del2 *= -1;
      }

    int ki, idx1, idx2;
    for (ki=0, idx1=start1, idx2=start2; ki<nmb1; ++ki, idx1+=del1, idx2+=del2)
      enumeration[ki] = make_pair(idx1, idx2);
 
    return true;
 }

//===========================================================================
  bool SurfaceTools::getCoefEnumeration(shared_ptr<SplineSurface> sf, int bd,
			  std::vector<int>& enumeration)
//===========================================================================
  {
    if (bd < 0 || bd > 3)
      return false;

    int kn1 = sf->numCoefs_u();
    int kn2 = sf->numCoefs_v();

    int nmb = (bd == 0 || bd == 1) ? kn2 : kn1;
    enumeration.resize(nmb);
    int start = (bd == 0 || bd == 2) ? 0 :
      ((bd == 1) ? kn1-1 : kn1*(kn2-1));
    int del = (bd == 0 || bd == 1) ? kn1 : 1;

    int ki, idx;
    for (ki=0, idx=start; ki<nmb; ++ki, idx+=del)
      enumeration[ki] = idx;

    return true;
  }

//===========================================================================
  bool SurfaceTools::getCoefEnumeration(shared_ptr<SplineSurface> sf, int bd,
					std::vector<int>& enumeration_bd,
					std::vector<int>& enumeration_bd2)
//===========================================================================
  {
    if (bd < 0 || bd > 3)
      return false;

    int kn1 = sf->numCoefs_u();
    int kn2 = sf->numCoefs_v();

    int nmb = (bd == 0 || bd == 1) ? kn2 : kn1;
    enumeration_bd.resize(nmb);
    enumeration_bd2.resize(nmb);
    int start = (bd == 0 || bd == 2) ? 0 :
      ((bd == 1) ? kn1-1 : kn1*(kn2-1));
    int del = (bd == 0 || bd == 1) ? kn1 : 1;

    int step = (bd == 0 || bd == 1) ? 1 : kn1;
    int sign = (bd == 0 || bd == 2) ? 1 : -1;
    int ki, idx;
    for (ki=0, idx=start; ki<nmb; ++ki, idx+=del)
    {
	enumeration_bd[ki] = idx;
	enumeration_bd2[ki] = idx + sign*step;
    }

    return true;
  }

//===========================================================================
  bool SurfaceTools::getCornerCoefEnum(shared_ptr<SplineSurface> sf, int bd1, int bd2,
			  int& enumeration)
//===========================================================================
  {
    if (bd1 < 0 || bd1 > 3 || bd2 < 0 || bd2 > 3)
      return false;
    if (bd1 == bd2)
      return false;  // No corner specified

    int kn1 = sf->numCoefs_u();
    int kn2 = sf->numCoefs_v();

    if ((bd1 == 0 || bd1 == 2) && (bd2 == 0 || bd2 == 2))
      enumeration = 0;
    else if  ((bd1 == 2 || bd1 == 1) && (bd2 == 2 || bd2 == 1))
      enumeration = kn1-1;
    else if  ((bd1 == 0 || bd1 == 3) && (bd2 == 0 || bd2 == 3))
      enumeration = (kn2-1)*kn1;
    else
      enumeration = kn1*kn2-1;

    return true;
  }

//===========================================================================
  int SurfaceTools::checkCoefCoLinearity(shared_ptr<SplineSurface> sf1,
			    shared_ptr<SplineSurface> sf2,
			    int bd1, int bd2, bool same_orient,
			    double tol, double ang_tol,
			    vector<vector<int> >& enumeration)
//===========================================================================
  {
    int colinear = 1;

    int dim = sf1->dimension();
    int kn1 = sf1->numCoefs_u();
    int kn2 = sf1->numCoefs_v();
    int kn3 = sf2->numCoefs_u();
    int kn4 = sf2->numCoefs_v();
    vector<double>::iterator c1 = sf1->coefs_begin();
    vector<double>::iterator c2 = sf2->coefs_begin();

    int nmb1 = (bd1 == 0 || bd1 == 1) ? kn2 : kn1;
    int nmb2 = (bd2 == 0 || bd2 == 1) ? kn4 : kn3;
    if (nmb1 != nmb2)
      return 0;  // No correspondence in number of coefficients

    enumeration.resize(nmb1);
    int start1 = (bd1 == 0 || bd1 == 2) ? 0 :
      ((bd1 == 1) ? kn1-1 : kn1*(kn2-1));
    int del1 = (bd1 == 0 || bd1 == 1) ? kn1 : 1;
    int del3 = (bd1 == 0 || bd1 == 1) ? 1 : kn1;
    if (bd1 == 1 || bd1 == 3)
      del3 *= -1;

    int start2 = (bd2 == 0 || bd2 == 2) ? 0 :
      ((bd2 == 1) ? kn3-1 : kn3*(kn4-1));
    int del2 = (bd2 == 0 || bd2 == 1) ? kn3 : 1;
    int del4 = (bd2 == 0 || bd2 == 1) ? 1 : kn3;
    if (bd2 == 1 || bd2 == 3)
      del4 *= -1;
    if (!same_orient)
      {
	start2 += (nmb2-1)*del2;
	del2 *= -1;
      }

    int ki, idx1, idx2;
    for (ki=0, idx1=start1, idx2=start2; ki<nmb1; ++ki, idx1+=del1, idx2+=del2)
      {
	// Check colinearity
	Point pos1(c1+idx1*dim, c1+(idx1+1)*dim);
	Point pos2(c2+idx2*dim, c2+(idx2+1)*dim);
	Point pos3(c1+(idx1+del3)*dim, c1+(idx1+del3+1)*dim);
	Point pos4(c2+(idx2+del4)*dim, c2+(idx2+del4+1)*dim);
	double ang = (pos1-pos3).angle(pos4-pos2);

	if (pos1.dist(pos2) > tol || ang > ang_tol)
	  colinear = 2;

	// Store coefficient indices
	vector<int> coef_enum(4);
	coef_enum[0] = idx1+del3;
	coef_enum[1] = idx1;
	coef_enum[2] = idx2;
	coef_enum[3] = idx2+del4;
	enumeration[ki] = coef_enum;
      }
 
    return colinear;
 }

//===========================================================================
// find a good seed for closest point computation
void SurfaceTools::surface_seedfind(const Point& pt, 
		      const ParamSurface& sf, 
		      const RectDomain* rd,
		      double& u,
		      double& v)
//===========================================================================
{
  // Evaluate a number of points in a rectangular grid and find the
  // closest one
  int nmb_sample = 7;
  vector<Point> samples(7*7);

  // Domain
  RectDomain dom;
  if (rd)
    dom = *rd;
  else
    dom = sf.containingDomain();

  int ki, kj, idx;
  double udel = (dom.umax() - dom.umin())/(int)(nmb_sample-1);
  double vdel = (dom.vmax() - dom.vmin())/(int)(nmb_sample-1);
  double upar, vpar;
  for (idx=0, kj=0, vpar=dom.vmin(); kj<nmb_sample; ++kj, vpar+=vdel)
    for (ki=0, upar=dom.umin(); ki<nmb_sample; ++idx, ++ki, upar+=udel)
	samples[idx] = sf.point(upar, vpar);

  // Find closest sampling point
  int min_idx=0;
  double min_dist = pt.dist(samples[0]);
  double dist;
  for (idx=1; idx<(int)samples.size(); ++idx)
    {
      dist = pt.dist(samples[idx]);
      if (dist < min_dist)
	{
	  min_dist = dist;
	  min_idx = idx;
	}
    }

  // Could check distance to neighbouring sample points and make the
  // seed more precise, but return only the corresponding parameter 
  // value for the time being
  kj = min_idx/nmb_sample;
  ki = (min_idx % nmb_sample);
  u = dom.umin() + ki*udel;
  v = dom.vmin() + kj*vdel;
}

//===========================================================================
// 
void SurfaceTools::parameterizeByBaseSurf(const  ParamSurface& sf, 
			    const vector<double>& points,
			    vector<double>& parvals)
//===========================================================================
{
  double eps = 1.0e-6;
  int dim = sf.dimension();
  int nmb_pts = (int)points.size()/dim;
  RectDomain dom = sf.containingDomain();

  for (int ki=0; ki<nmb_pts; ++ki)
    {
      double seed[2];
      double u, v, dist;
      Point close;
      Point curr(points.begin()+ki*dim, points.begin()+(ki+1)*dim);
      SurfaceTools::surface_seedfind(curr, sf, &dom, seed[0], seed[1]);
      sf.closestPoint(curr, u, v, close, dist, eps, &dom, seed);
      parvals.push_back(u);
      parvals.push_back(v);
    }
}

//===========================================================================
  double SurfaceTools::estimateTangentLength(ParamSurface *surf, int pardir, 
			       bool at_start)
//===========================================================================
  {
    int nmb_sample = 5;
    vector<Point> pts(3);
    double len = 0.0;
    RectDomain dom = surf->containingDomain();
    if (pardir == 1)
      {
	double upar = (at_start) ? dom.umin() : dom.umax();
	double vpar = dom.vmin(); 
	double del = (dom.vmax() - vpar)/(double)(nmb_sample-1);
	for (int ki=0; ki<nmb_sample; ++ki, vpar+=del)
	  {
	    surf->point(pts, upar, vpar, 1);
	    len += pts[1].length();
	  }
      }
    else
      {
	double upar = dom.umin();
	double vpar = (at_start) ? dom.vmin() : dom.vmax();
	double del = (dom.umax() - upar)/(double)(nmb_sample-1);
	for (int ki=0; ki<nmb_sample; ++ki, upar+=del)
	  {
	    surf->point(pts, upar, vpar, 1);
	    len += pts[2].length();
	  }
      }
    len /= (double)nmb_sample;
    return len;
  }

//===========================================================================
void SurfaceTools::checkSurfaceClosed(const ParamSurface& sf,
				      bool& closed_dir_u, bool& closed_dir_v,
				      double closed_tol_geo)
//===========================================================================
{
    closed_dir_u = false;
    closed_dir_v = false;

    // We convert from the geometric tolerance to the corresponding parametric tolerance.
    const Point sf_epspar = SurfaceTools::getParEpsilon(sf, closed_tol_geo);

    if (sf.instanceType() == Class_SplineSurface) {
        const SplineSurface& spline_sf =
            dynamic_cast<const SplineSurface&>(sf);
        // Making sure the sf is k-regular.
        shared_ptr<SplineSurface> spline_under_sf(
            spline_sf.subSurface(spline_sf.startparam_u(),
                                spline_sf.startparam_v(),
                                spline_sf.endparam_u(),
                                spline_sf.endparam_v()));
        surfaceClosed(*spline_under_sf, 
            closed_dir_u, closed_dir_v, closed_tol_geo);
        return;
    }
    else if (sf.instanceType() == Class_Cone ||
        sf.instanceType() == Class_Cylinder ||
        sf.instanceType() == Class_Sphere ||
        sf.instanceType() == Class_Torus) {
            const ElementarySurface& es
                = dynamic_cast<const ElementarySurface&>(sf);
            es.isClosed(closed_dir_u, closed_dir_v);
            return;
    }
    else if (sf.instanceType() == Class_SurfaceOfRevolution) {
        const SurfaceOfRevolution& sor
            = dynamic_cast<const SurfaceOfRevolution&>(sf);

        // u-direction
        const RectDomain& domain =
            dynamic_cast<const RectDomain&>(sor.parameterDomain());
        double interval_length_u = domain.umax() - domain.umin();
        double interval_length_v = domain.vmax() - domain.vmin();
        if (sor.isSwapped())
            swap(interval_length_u, interval_length_v);
        if (fabs(interval_length_u - 2.0*M_PI) < sf_epspar[0])
            closed_dir_u = true;

        // v-direction
        shared_ptr<ParamCurve> pc 
            = dynamic_pointer_cast<ParamCurve>(sor.getCurve());
        closed_dir_v = pc->isClosed();

        if (sor.isSwapped())
            swap(closed_dir_u, closed_dir_v);

        return;
    }


    return;
}


//===========================================================================
void SurfaceTools::surfaceClosed(const SplineSurface& sf,
				 bool& closed_dir_u, bool& closed_dir_v,
				 double closed_tol_geo)
//===========================================================================
{
  // Assuming k-regular surface, summing dist (L1-norm) between
  // corresponding row coefs.
  double num_tol = 1e-12;
  int ik1 = sf.order_u();
  int ik2 = sf.order_v();
  int in1 = sf.numCoefs_u();
  int in2 = sf.numCoefs_v();
  vector<double>::const_iterator et1 = sf.basis_u().begin();
  vector<double>::const_iterator et2 = sf.basis_v().begin();

  if (!((et1[ik1-1] - et1[0] < num_tol) &&
         (et1[in1+ik1-1] - et1[in1] < num_tol) &&
         (et2[ik2-1] - et2[0] < num_tol) &&
         (et2[in2+ik2-1] - et2[in2] < num_tol))) {
             THROW("Surface not k-regular");
  }
  //ASSERT((et1[ik1-1] - et1[0] < num_tol) &&
  //       (et1[in1+ik1-1] - et1[in1] < num_tol) &&
  //       (et2[ik2-1] - et2[0] < num_tol) &&
  //       (et2[in2+ik2-1] - et2[in2] < num_tol));

  // Current case is not rational ...
  //ASSERT(!sf.rational());

  // We first check vmin vs vmax.
  double sum_dist = 0.0;
  int dim = sf.dimension();
  vector<double>::const_iterator coefs;
//   if (sf.rational())
//     {
//       coefs = sf.rcoefs_begin();
//       dim++;
//     }
//   else
    coefs = sf.coefs_begin();
  for (int ki = 0; ki < in1*dim; ++ki)
    sum_dist += fabs(coefs[in1*dim*(in2-1)+ki] - coefs[ki]);
  closed_dir_v = (sum_dist < closed_tol_geo);

  // Then umin vs umax.
  sum_dist = 0.0;
  for (int ki = 0; ki < in2; ++ki)
    for (int kj = 0; kj < dim; ++kj)
      sum_dist += fabs(coefs[(ki*in1+in1-1)*dim+kj] - coefs[ki*in1*dim+kj]);
  closed_dir_u = (sum_dist < closed_tol_geo);

  return;
}


//===========================================================================
Point SurfaceTools::getParEpsilon(const ParamSurface& sf, double epsgeo)
//===========================================================================
{
    Point sf_epspar(2);

    // A tolerance is typically shared between two neighbours. Also we
    // do not want to approximate close to the exact tolerance.
    const double scaling = 0.5;

    double sf_length_u, sf_length_v; // Average values, sampled.
    // For a cylinder the size is infinite. We assume that such cases
    // are curve length parametrized in that direction.
    if (sf.instanceType() == Class_Cylinder)
    {
	const Cylinder& cyl = dynamic_cast<const Cylinder&>(sf);
	sf_epspar[0] = epsgeo/cyl.getRadius();
	sf_epspar[1] = epsgeo;
	if (cyl.isSwapped())
	{
	    std::swap(sf_epspar[0], sf_epspar[1]);
	}
    }
    else if (sf.instanceType() == Class_SurfaceOfLinearExtrusion)
    {
	const SurfaceOfLinearExtrusion& sole = dynamic_cast<const SurfaceOfLinearExtrusion&>(sf);
        double cv_length = sole.getCurve()->estimatedCurveLength();
        double umin = sole.getCurve()->startparam();
        double umax = sole.getCurve()->endparam();
	sf_epspar[0] = epsgeo*(umax - umin)/cv_length;
	sf_epspar[1] = epsgeo;
	if (sole.isSwapped())
	{
	    std::swap(sf_epspar[0], sf_epspar[1]);
	}
    }
    else if (sf.instanceType() == Class_Cone)
    {
	const Cone& cone = dynamic_cast<const Cone&>(sf);
	sf_epspar[0] = epsgeo/cone.getRadius();
	sf_epspar[1] = epsgeo;
	if (cone.isSwapped())
	{
	    std::swap(sf_epspar[0], sf_epspar[1]);
	}
    }
    else
    {
      // Set number of sampling points based on surface size
      BoundingBox bbox = sf.boundingBox();
      double len = bbox.low().dist(bbox.high());
      double fac = 10000.0;
      int num_samples = (len/epsgeo < fac) ? 5 : 10;
      sf.ParamSurface::estimateSfSize(sf_length_u, sf_length_v, 
				      num_samples, num_samples);

	RectDomain rect_dom = sf.containingDomain();
	double dom_length_u = rect_dom.umax() - rect_dom.umin();
	double dom_length_v = rect_dom.vmax() - rect_dom.vmin();

	// Subtracting 1.0 to avoid numerical issues.
	bool u_dom_inf = (dom_length_u > MAXDOUBLE - 1.0);
	bool v_dom_inf = (dom_length_v > MAXDOUBLE - 1.0);

	sf_epspar[0] = (u_dom_inf) ? epsgeo : epsgeo*dom_length_u/sf_length_u;
	sf_epspar[1] = (v_dom_inf) ? epsgeo : epsgeo*dom_length_v/sf_length_v;

	// Avoid extremely large tolerances
	double lim = 0.001*std::max(dom_length_u, dom_length_v);
	sf_epspar[0] = std::min(lim, sf_epspar[0]);
	sf_epspar[1] = std::min(lim, sf_epspar[1]);
    }

    sf_epspar[0] *= scaling;
    sf_epspar[1] *= scaling;

    return sf_epspar;
}

  //===========================================================================
  void 
  SurfaceTools::setResolutionFromDensity(shared_ptr<ParamSurface> surf,
					 double density, 
					 int min_nmb, int max_nmb,
					 int& u_res, int& v_res)
  //===========================================================================
  {
    // Estimate size of surface/underlying surface
    RectDomain dom = surf->containingDomain();
	
    double len_u, len_v;
    surf->estimateSubSfSize(dom.umin(), dom.umax(), len_u, 
			    dom.vmin(), dom.vmax(), len_v);

    u_res = (int)(len_u/density);
    v_res = (int)(len_v/density);
    double fac = len_u/len_v;
    u_res = std::max(min_nmb, std::min(u_res, (int)(fac*max_nmb)));
    v_res = std::max(min_nmb, std::min(v_res, (int)(max_nmb/fac)));

  }



} // end namespace Go
