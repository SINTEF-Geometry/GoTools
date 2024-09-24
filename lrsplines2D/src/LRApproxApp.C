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


#include "GoTools/lrsplines2D/LRApproxApp.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSplineMBA.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/creators/Eval1D3DSurf.h"
#include "GoTools/utils/omp.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;
using std::string;

//#define DEBUG


//=============================================================================
void LRApproxApp::pointCloud2Spline(vector<double>& points, int dim,
				    double domain[], double reduced_domain[],
				    double eps, int max_iter,
				    shared_ptr<LRSplineSurface>& surf,
				    double& maxdist, double& avdist, 
				    double& avdist_out, int& nmb_out,
				    int mba, int initmba, int tomba)
//=============================================================================
{
  // Define parameters
  double smoothwg = 1.0e-10; 
  if (tomba == 5)
    tomba = std::min(tomba, max_iter-1);    // Turn to the mba method at 
  // iteration level 5 or in the last iteration

  // Translate data points to origo
  int del = 2+dim;
  int nmb_points = (int)points.size()/del;

  int ki, kj;
  double mid[2];
  mid[0] = 0.5*(domain[0] + domain[1]);
  mid[1] = 0.5*(domain[2] + domain[3]);
  double low = std::numeric_limits<double>::max();
  double high = std::numeric_limits<double>::lowest();
  for (ki=0; ki<nmb_points; ++ki)
    {
      for (kj=del-3; kj<del-1; ++kj)
	{
	  points[del*ki+kj] -= mid[kj-del+3];
	}
      low = std::min(low, points[del*ki+kj]);
      high = std::max(high, points[del*ki+kj]);
    }

#ifdef DEBUG
  // Write translated points to g2 format
  vector<double> points2;
  points2.reserve(nmb_points*3);
  for (ki=0, kj=0; ki<nmb_points; ++ki, kj+=del)
    points2.insert(points2.end(), points.begin()+kj, points.begin()+kj+3);
  PointCloud3D cloud(points2.begin(), nmb_points);

  std::ofstream of2("translated_points.g2");
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
  
  std::cout << std::endl;
  std::cout << "Input points read and pre processed. Ready for surface creation.";
  std::cout << std::endl << std::endl;
#endif
  
  // Make approximation engine
  // First make initial tensor-product spline space
  double ntol = 1.0e-15;
  int nmb_coef = 14;
  int order = 3; 
  vector<double> knots_u;
  vector<double> knots_v;
  knots_u.reserve(nmb_coef+order);
  knots_v.reserve(nmb_coef+order);
  double del_u = (reduced_domain[1]-reduced_domain[0])/(double)(nmb_coef-order-1);
  double del_v = (reduced_domain[3]-reduced_domain[2])/(double)(nmb_coef-order-1);
  for (ki=0; ki<order; ++ki)
    knots_u.push_back(domain[0]-mid[0]);
  if (reduced_domain[0] > domain[0]+ntol)
    knots_u.push_back(reduced_domain[0]-mid[0]);
  for (ki=1; ki<nmb_coef-order-1; ++ki)
    knots_u.push_back(reduced_domain[0]-mid[0]+ki*del_u);
  if (reduced_domain[1] < domain[1]-ntol)
    knots_u.push_back(reduced_domain[1]-mid[0]);
  for (ki=0; ki<order; ++ki)
    knots_u.push_back(domain[1]-mid[0]);

  for (ki=0; ki<order; ++ki)
    knots_v.push_back(domain[2]-mid[1]);
  if (reduced_domain[2] > domain[2]+ntol)
    knots_v.push_back(reduced_domain[2]-mid[1]);
  for (ki=1; ki<nmb_coef-order-1; ++ki)
    knots_v.push_back(reduced_domain[2]-mid[1]+ki*del_v);
  if (reduced_domain[3] < domain[3]-ntol)
    knots_v.push_back(reduced_domain[3]-mid[1]);
  for (ki=0; ki<order; ++ki)
    knots_v.push_back(domain[3]-mid[1]);

  double mba_coef = 0.0;
  if (initmba)
    mba_coef = 0.5*(low+high);
  LRSurfApprox approx(order, knots_u, order, knots_v, points, del-2, 
		      eps, initmba ? true : false, mba_coef, true, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  if (mba)
    approx.setUseMBA(true);
  else
    {
      approx.setSwitchToMBA(tomba);
      approx.setMakeGhostPoints(true);
    }
  //approx.setVerbose(true);
  if (del == 3)
    {
      double zrange = high - low;
      double zfac = std::max(eps, 0.005*zrange);
      approx.addLowerConstraint(low - zfac);
      approx.addUpperConstraint(high + zfac);
      approx.setLocalConstraint(zfac);
    }

  // Approximate
  surf = approx.getApproxSurf(maxdist, avdist,avdist_out, nmb_out, max_iter);

#ifdef DEBUG
  std::cout << std::endl;
  std::cout << "Approximation completed. " << std::endl;
 
  std::cout << "Total number of points: " << nmb_points << std::endl;
  std::cout << "Number of elements: " << surf->numElements() << std::endl;
  std::cout << "Maximum distance: " << maxdist << std::endl;
  std::cout << "Average distance: " << avdist << std::endl;
  std::cout << "Average distance for points outside of the tolerance: " << avdist_out << std::endl;
  std::cout << "Number of points outside the tolerance: " << nmb_out << std::endl;
#endif
  if (surf.get())
    {
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (surf->dimension() == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	}
#endif
	  
      // Translate
      if (surf->dimension() == 3)
	{
	  Point tmp(mid[0], mid[1], 0.0);
	  surf->translate(tmp);
	}
      else
	{
	  // Update parameter domain
	  double umin = surf->paramMin(XFIXED);
	  double umax = surf->paramMax(XFIXED);
	  double vmin = surf->paramMin(YFIXED);
	  double vmax = surf->paramMax(YFIXED);

	  surf->setParameterDomain(umin + mid[0], umax + mid[0],
				   vmin + mid[1], vmax + mid[1]);
	}  
    }
}

//=============================================================================
void LRApproxApp::pointCloud2Spline(vector<double>& points, 
				    shared_ptr<LRSplineSurface>& init_surf,
				    vector<double>& extent,
				    double eps, int max_iter,
				    shared_ptr<LRSplineSurface>& surf,
				    double& maxdist, double& avdist, 
				    double& avdist_out, int& nmb_out,
				    int mba, int tomba)
//=============================================================================
{
  // Define parameters
  double smoothwg = 1.0e-10; 
  if (tomba >= max_iter)
    tomba = std::min(tomba, max_iter-1);    

  // Translate data points to origo
  int dim = init_surf->dimension();
  int del = 2+dim;
  int nmb_points = (int)points.size()/del;

 // Move point cloud to origo
  double umin = init_surf->paramMin(XFIXED);
  double umax = init_surf->paramMax(XFIXED);
  double vmin = init_surf->paramMin(YFIXED);
  double vmax = init_surf->paramMax(YFIXED);
  Point mid;
  int ki, kj;
  if (dim == 1)
    mid.setValue(0.5*(umin+umax), 0.5*(vmin+vmax), 0.0);
  else
    {
      mid = Point(0.5*(extent[2*(del-3) + 2*(del-3)+1]),
		  0.5*(extent[2*(del-2) + 2*(del-2)+1]), 0.0);
    }
  for (ki=0; ki<nmb_points; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      {
	points[del*ki+kj] -= mid[kj-del+3];
      }

  // Move surface to origo
  if (dim == 1)
    {
      init_surf->setParameterDomain(umin - mid[0], umax - mid[0],
				vmin - mid[1], vmax - mid[1]);
    }
  else
    init_surf->translate(-mid);
      
#ifdef DEBUG
  // Write translated points to g2 format
  vector<double> points2;
  points2.reserve(nmb_points*3);
  for (ki=0, kj=0; ki<nmb_points; ++ki, kj+=del)
    points2.insert(points2.end(), points.begin()+kj, points.begin()+kj+3);
  PointCloud3D cloud(points2.begin(), nmb_points);

  std::ofstream of2("translated_points.g2");
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
  
  std::cout << std::endl;
  std::cout << "Input points read and pre processed. Ready for surface creation.";
  std::cout << std::endl << std::endl;
#endif
  
  // Make approximation engine
  bool repar = true;
  LRSurfApprox approx(init_surf, points, eps, true, repar, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  if (mba)
    approx.setUseMBA(true);
  else
    {
      approx.setSwitchToMBA(tomba);
      approx.setMakeGhostPoints(true);
    }
  //approx.setVerbose(true);
  if (del == 3)
    {
      double zrange = extent[5] - extent[4];
      double zfac = std::max(eps, 0.005*zrange);
      approx.addLowerConstraint(extent[4] - zfac);
      approx.addUpperConstraint(extent[5] + zfac);
      approx.setLocalConstraint(zfac);
    }

  // Approximate
  surf = approx.getApproxSurf(maxdist, avdist,avdist_out, nmb_out, max_iter);

#ifdef DEBUG
  std::cout << std::endl;
  std::cout << "Approximation completed. " << std::endl;
 
  std::cout << "Total number of points: " << nmb_points << std::endl;
  std::cout << "Number of elements: " << surf->numElements() << std::endl;
  std::cout << "Maximum distance: " << maxdist << std::endl;
  std::cout << "Average distance: " << avdist << std::endl;
  std::cout << "Average distance for points outside of the tolerance: " << avdist_out << std::endl;
  std::cout << "Number of points outside the tolerance: " << nmb_out << std::endl;
#endif
  if (surf.get())
    {
#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      if (dim == 3)
	{
	  surf->writeStandardHeader(of1);
	  surf->write(of1);
	}
      else
	{
	  shared_ptr<LRSplineSurface> surf2(surf->clone());
	  surf2->to3D();
	  surf2->writeStandardHeader(of1);
	  surf2->write(of1);
	}
#endif
	  
      // Translate
      if (dim == 3)
	{
	  Point tmp(mid[0], mid[1], 0.0);
	  surf->translate(tmp);
	  init_surf->translate(tmp);
	}
      else
	{
	  // Update parameter domain
	  double umin = surf->paramMin(XFIXED);
	  double umax = surf->paramMax(XFIXED);
	  double vmin = surf->paramMin(YFIXED);
	  double vmax = surf->paramMax(YFIXED);

	  surf->setParameterDomain(umin + mid[0], umax + mid[0],
				   vmin + mid[1], vmax + mid[1]);
	  init_surf->setParameterDomain(umin + mid[0], umax + mid[0],
					vmin + mid[1], vmax + mid[1]);
	}  
    }
}

static int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

static int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//=============================================================================
void LRApproxApp::computeDistPointSpline(vector<double>& points,
					 shared_ptr<LRSplineSurface>& surf,
					 double& max_above, double& max_below, 
					 double& avdist, int& nmb_points,
					 vector<double>& pointsdist,
					 int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height

  pointsdist.reserve(nmb_pts*4);

  // Get all knot values in the u-direction
  const double* const uknots = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));
  double aeps = 0.001;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = std::numeric_limits<double>::lowest();
  max_below = std::numeric_limits<double>::max();
  avdist = 0.0;
  nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  int ki, kj, kr, ka;
  double *curr;
  double dist;

  int pp0, pp1;
  Element2D* elem = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)points.size() && points[pp0+1] < (*knotv); 
       pp0+=3);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)points.size() && points[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)points.size();

      // Sort the current sub set of points according to the u-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant points and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && points[pp2] < (*knotu); pp2+=3);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	  elem = elements[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      surf->point(pos, curr[0], curr[1], elem);
	      dist = curr[2]-pos[0];

	      if (use_proj)
		{
		  Point clo_pt;
		  double clo_u, clo_v, clo_dist;
		  double seed[2];
		  seed[0] = curr[0];
		  seed[1] = curr[1];
		  Point pt(curr[0], curr[1], curr[2]);
		  surf->setCurrentElement(elem);
		  evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					clo_dist, aeps, 1, seed);
		  if (clo_dist < fabs(dist))
		    dist = (dist < 0.0) ? -clo_dist : clo_dist;
		}

	      max_above = std::max(max_above, dist);
	      max_below = std::min(max_below, dist);
	      avdist += fabs(dist);
	      pointsdist.push_back(curr[0]);
	      pointsdist.push_back(curr[1]);
	      pointsdist.push_back(curr[2]);
	      pointsdist.push_back(dist);

	      nmb_points++;
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }
  if (nmb_points > 0)
    avdist /= nmb_points;
}

//=============================================================================
void LRApproxApp::computeDistPointSpline_omp(vector<double>& points,
					     shared_ptr<LRSplineSurface>& surf,
					     double& max_above, double& max_below, 
					     double& avdist, int& nmb_points,
					     vector<double>& pointsdist,
					     int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height

  pointsdist.reserve(nmb_pts*4);

  // Get all knot values in the u-direction
  const double* const uknots_begin = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  const int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots_begin = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = max_below = avdist = 0.0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  int pp0, pp1;

  int num_kj = vknots_end - vknots_begin - 1; // Threshold for the number of elements in the v-dir.
  vector<int> num_pts(num_kj, 0);
  vector<vector<double> > pts_dist(num_kj);
  int ki, kj, kr;
#pragma omp parallel default(none) private(ki, kj, kr, knotv, knotu, pp0, pp1) \
  shared(surf, points, num_pts, pts_dist, num_kj, elements, evalsrf, uknots_begin, uknots_end, nmb_knots_u, vknots_begin, vknots_end)
  {
      Point pos;
      int nump;
      int pp2, pp3;
      Element2D* elem = NULL;
      double *curr;
      double dist;
      double aeps = 0.001;
#pragma omp for OMP_SCHEDULE_AUTO
      for (kj=0; kj < num_kj; ++kj)
      {
	  knotv = vknots_begin + kj; // Left side of global element.
	  // We locate first index of point belonging to element.
	  for (pp0=0; pp0<(int)points.size() && points[pp0+1] < knotv[0]; pp0+=3);

	  // We locate last index of point belonging to element.
	  for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < knotv[1]; pp1+=3);

	  if (knotv+2 == vknots_end)
	      for (; pp1<(int)points.size() && points[pp1+1] <= knotv[1]; pp1+=3);

	  // Sort the current sub set of points according to the u-parameter
	  qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

	  // Traverse the relevant points and identify the associated element
	  for (pp2=pp0, knotu=uknots_begin; pp2<pp1 && points[pp2] < (*knotu); pp2+=3);
	  for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	  {
	      for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	      if (knotu+1 == uknots_end)
		  for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	      //   pp3 = pp1;
	  
	      // Fetch associated element
	      elem = elements[kj*(nmb_knots_u-1)+ki];

	      nump = (pp3 - pp2)/3;
	      for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	      {
		  // Evaluate
		  surf->point(pos, curr[0], curr[1], elem);
		  dist = curr[2]-pos[0];

		  if (evalsrf.get())
		    {
		      Point clo_pt;
		      double clo_u, clo_v, clo_dist;
		      double seed[2];
		      seed[0] = curr[0];
		      seed[1] = curr[1];
		      Point pt(curr[0], curr[1], curr[2]);
		      surf->setCurrentElement(elem);
		      evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					    clo_dist, aeps, 1, seed);
		      if (clo_dist < fabs(dist))
			dist = (dist < 0.0) ? -clo_dist : clo_dist;
		    }

		  pts_dist[kj].push_back(curr[0]);
		  pts_dist[kj].push_back(curr[1]);
		  pts_dist[kj].push_back(curr[2]);
		  pts_dist[kj].push_back(dist);

		  num_pts[kj]++;
	      }
	      pp2 = pp3;
	  }
	  pp0 = pp1;
      } // End of parallel loop.
  } // End of parallel region.

  nmb_points = 0;
  double dist;
  for (kj = 0; kj < (int)num_pts.size(); ++kj)
  {
      nmb_points += num_pts[kj];
      for (ki = 0; ki < num_pts[kj]; ++ki)
      {
	  dist = pts_dist[kj][ki*4+3];
	  max_above = std::max(max_above, dist);
	  max_below = std::min(max_below, dist);
	  avdist += fabs(dist);
      }
      pointsdist.insert(pointsdist.end(), pts_dist[kj].begin(), pts_dist[kj].end());
  }
  if (nmb_points > 0)
    avdist /= nmb_points;
}


//=============================================================================
void LRApproxApp::classifyCloudFromDist(vector<double>& points,
					shared_ptr<LRSplineSurface>& surf,
					vector<double>& limits,
					double& max_above, double& max_below, 
					double& avdist, int& nmb_points,
					vector<vector<double> >& level_points,
					vector<int>& nmb_group,
					int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height

  // Get all knot values in the u-direction
  const double* const uknots = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));
  double aeps = 0.001;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = max_below = avdist = 0.0;
  nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  int ki, kj, kr, ka;
  double *curr;
  double dist;

  int pp0, pp1;
  Element2D* elem = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)points.size() && points[pp0+1] < (*knotv); 
       pp0+=3);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)points.size() && points[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)points.size();

      // Sort the current sub set of points according to the u-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant points and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && points[pp2] < (*knotu); pp2+=3);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	  elem = elements[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      surf->point(pos, curr[0], curr[1], elem);
	      dist = curr[2]-pos[0];

	      if (evalsrf.get())
		{
		  Point clo_pt;
		  double clo_u, clo_v, clo_dist;
		  double seed[2];
		  seed[0] = curr[0];
		  seed[1] = curr[1];
		  Point pt(curr[0], curr[1], curr[2]);
		  surf->setCurrentElement(elem);
		  evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					clo_dist, aeps, 1, seed);
		  if (clo_dist < fabs(dist))
		    dist = (dist < 0.0) ? -clo_dist : clo_dist;
		}

	      max_above = std::max(max_above, dist);
	      max_below = std::min(max_below, dist);
	      avdist += fabs(dist);
	      nmb_points++;
	  
	      // Find classification
	      for (ka=0; ka<(int)limits.size(); ++ka)
		if (dist < limits[ka])
		  {
		    level_points[ka].push_back(curr[0]);
		    level_points[ka].push_back(curr[1]);
		    level_points[ka].push_back(curr[2]);
		    break;
		  }
	      if (ka == (int)limits.size())
		{
		  level_points[ka].push_back(curr[0]);
		  level_points[ka].push_back(curr[1]);
		  level_points[ka].push_back(curr[2]);
		}
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }
  avdist /= nmb_points;

  nmb_group.resize(level_points.size());
  for (size_t kk=0; kk<nmb_group.size(); ++kk)
    nmb_group[kk] = (int)level_points[kk].size()/3;
}



//=============================================================================
void LRApproxApp::classifyCloudFromDist_omp(vector<double>& points,
					    shared_ptr<LRSplineSurface>& surf,
					    vector<double>& limits,
					    double& max_above, double& max_below, 
					    double& avdist, int& nmb_points,
					    vector<vector<double> >& level_points,
					    vector<int>& nmb_group,
					    int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height

  // Get all knot values in the u-direction
  const double* const uknots_begin = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  const int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);

  // Get all knot values in the v-direction
  const double* const vknots_begin = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = max_below = avdist = 0.0;
  nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  const int num_kj = vknots_end - vknots_begin - 1; // Threshold for the number of elements in the v-dir.

  vector<double> all_max_above(num_kj, 0.0);
  vector<double> all_max_below(num_kj, 0.0);
  vector<double> all_dist(num_kj, 0.0);
  vector<int> all_nmb_points(num_kj, 0);
  vector<vector<vector<double> > > all_level_points(num_kj, vector<vector<double> >(level_points.size()));;

  int kj;
#pragma omp parallel default(none) private(kj) \
  shared(surf, points, limits, elements, all_max_above, all_max_below, all_dist, all_nmb_points, all_level_points, evalsrf, uknots_begin, uknots_end, nmb_knots_u, vknots_begin, vknots_end, num_kj)
  {
      Element2D* elem = NULL;
      int ki, kr, ka;
      int pp0, pp1;
      int pp2, pp3;
      Point pos;
      int nump;
      double *curr;
      double dist;
      const double* knotu;
      const double* knotv;
      double aeps = 0.001;
#pragma omp for OMP_SCHEDULE_AUTO
      for (kj = 0; kj < num_kj; ++kj)
      {
      	  knotv = vknots_begin + kj; // Left side of global element.

	  for (pp0=0; pp0<(int)points.size() && points[pp0+1] < knotv[0]; pp0+=3);

	  for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < knotv[1]; pp1+=3);

	  if (knotv+2 == vknots_end)
	      for (; pp1<(int)points.size() && points[pp1+1] <= knotv[1]; pp1+=3);
	  // 	pp1 = (int)points.size();

	  // Sort the current sub set of points according to the u-parameter
	  qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

	  // Traverse the relevant points and identify the associated element
	  for (pp2=pp0, knotu=uknots_begin; pp2<pp1 && points[pp2] < (*knotu); pp2+=3);
	  for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	  {
	      for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	      if (knotu+1 == uknots_end)
		  for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	      //   pp3 = pp1;
	  
	      // Fetch associated element
	      elem = elements[kj*(nmb_knots_u-1)+ki];

	      nump = (pp3 - pp2)/3;
	      for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	      {
		  // Evaluate
		  surf->point(pos, curr[0], curr[1], elem);
		  dist = curr[2]-pos[0];

		  if (evalsrf.get())
		    {
		      Point clo_pt;
		      double clo_u, clo_v, clo_dist;
		      double seed[2];
		      seed[0] = curr[0];
		      seed[1] = curr[1];
		      Point pt(curr[0], curr[1], curr[2]);
		      surf->setCurrentElement(elem);
		      evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					    clo_dist, aeps, 1, seed);
		      if (clo_dist < fabs(dist))
			dist = (dist < 0.0) ? -clo_dist : clo_dist;
		    }

		  all_max_above[kj] = std::max(all_max_above[kj], dist);
		  all_max_below[kj] = std::min(all_max_below[kj], dist);
		  all_dist[kj] += fabs(dist);
		  all_nmb_points[kj]++;
	  
		  // Find classification
		  for (ka=0; ka<(int)limits.size(); ++ka)
		      if (dist < limits[ka])
		      {
			  all_level_points[kj][ka].push_back(curr[0]);
			  all_level_points[kj][ka].push_back(curr[1]);
			  all_level_points[kj][ka].push_back(curr[2]);
			  break;
		      }
		  if (ka == (int)limits.size())
		  {
		      all_level_points[kj][ka].push_back(curr[0]);
		      all_level_points[kj][ka].push_back(curr[1]);
		      all_level_points[kj][ka].push_back(curr[2]);
		  }
	      }
	      pp2 = pp3;
	  }
	  pp0 = pp1;
      } // End of parallel for loop.
  } // End of parallel region.


  for (kj = 0; kj < num_kj; ++kj)
  {
      max_above = std::max(all_max_above[kj], max_above);
      max_below = std::min(all_max_below[kj], max_below);
      avdist += all_dist[kj];
      nmb_points += all_nmb_points[kj];

      for (size_t ki = 0; ki < all_level_points[kj].size(); ++ki)
      {
	  level_points[ki].insert(level_points[ki].end(),
				  all_level_points[kj][ki].begin(), all_level_points[kj][ki].end());
      }
  }

  avdist /= nmb_points;

  nmb_group.resize(level_points.size());
  for (size_t kk=0; kk<nmb_group.size(); ++kk)
      nmb_group[kk] = (int)level_points[kk].size()/3;
}


//=============================================================================
void LRApproxApp::categorizeCloudFromDist(vector<double>& points,
					  shared_ptr<LRSplineSurface>& surf,
					  vector<double>& limits,
					  double& max_above, double& max_below, 
					  double& avdist, int& nmb_points,
					  vector<int>& classification,
					  vector<int>& nmb_group,
					  int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height
  classification.clear();
  classification.reserve(nmb_pts);
  nmb_group.resize(limits.size()+1);
  for (size_t kk=0; kk<nmb_group.size(); ++kk)
    nmb_group[kk] = 0;

  // Get all knot values in the u-direction
  const double* const uknots = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));
  double aeps = 0.001;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = max_below = avdist = 0.0;
  nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  int ki, kj, kr, ka;
  double *curr;
  double dist;

  int pp0, pp1;
  Element2D* elem = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)points.size() && points[pp0+1] < (*knotv); 
       pp0+=3)
    classification.push_back(-1);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)points.size() && points[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)points.size();

      // Sort the current sub set of points according to the u-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant points and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && points[pp2] < (*knotu); pp2+=3)
	classification.push_back(-1);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	  elem = elements[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      surf->point(pos, curr[0], curr[1], elem);
	      dist = curr[2]-pos[0];

	      if (evalsrf.get())
		{
		  Point clo_pt;
		  double clo_u, clo_v, clo_dist;
		  double seed[2];
		  seed[0] = curr[0];
		  seed[1] = curr[1];
		  Point pt(curr[0], curr[1], curr[2]);
		  surf->setCurrentElement(elem);
		  evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					clo_dist, aeps, 1, seed);
		  if (clo_dist < fabs(dist))
		    dist = (dist < 0.0) ? -clo_dist : clo_dist;
		}

	      max_above = std::max(max_above, dist);
	      max_below = std::min(max_below, dist);
	      avdist += fabs(dist);
	      nmb_points++;
	  
	      // Find classification
	      for (ka=0; ka<(int)limits.size(); ++ka)
		if (dist < limits[ka])
		  {
		    nmb_group[ka]++;
		    classification.push_back(ka);
		    break;
		  }
	      if (ka == (int)limits.size())
		{
		  nmb_group[ka]++;
		  classification.push_back(ka);
		}
	    }
	  pp2 = pp3;
	}
      for (; pp2<pp1; pp2+=3)
	classification.push_back(-1);
      pp0 = pp1;
    }
  for (; pp0<(int)points.size(); pp0+=3)
    classification.push_back(-1);

  avdist /= nmb_points;
}


//=============================================================================
void LRApproxApp::categorizeCloudFromDist_omp(vector<double>& points,
					      shared_ptr<LRSplineSurface>& surf,
					      vector<double>& limits,
					      double& max_above, double& max_below, 
					      double& avdist, int& nmb_points,
					      vector<int>& classification,
					      vector<int>& nmb_group,
					      int use_proj)
//=============================================================================
{
  if (surf->dimension() != 1)
    return;   // Not handled
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height
  classification.clear();
  classification.reserve(nmb_pts);
  nmb_group.resize(limits.size()+1);
  for (size_t kk=0; kk<nmb_group.size(); ++kk)
    nmb_group[kk] = 0;

  // Get all knot values in the u-direction
  const double* const uknots_begin = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  const int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);

  // Get all knot values in the v-direction
  const double* const vknots_begin = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);

  shared_ptr<Eval1D3DSurf> evalsrf;
  if (use_proj)
    evalsrf = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(surf));

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  max_above = max_below = avdist = 0.0;
  nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);

  const int num_kj = vknots_end - vknots_begin - 1; // Threshold for the number of elements in the v-dir.

  vector<double> all_max_above(num_kj, 0.0);
  vector<double> all_max_below(num_kj, 0.0);
  vector<double> all_dist(num_kj, 0.0);
  vector<int> all_nmb_points(num_kj, 0);
  vector<vector<int> > all_classification(num_kj);
  vector<vector<int> > all_nmb_group(num_kj, vector<int>(limits.size() + 1, 0));

  int kj;

#pragma omp parallel default(none) private(kj) \
  shared(surf, points, limits, elements, all_max_above, all_max_below, all_dist, all_nmb_points, all_classification, all_nmb_group, evalsrf, uknots_begin, uknots_end, nmb_knots_u, vknots_begin, vknots_end, num_kj)
  {
      Element2D* elem = NULL;
      int ki, kr, ka;
      int pp0, pp1;
      int pp2, pp3;
      Point pos;
      int nump;
      double *curr;
      double dist;
      const double* knotu;
      const double* knotv;
      double aeps = 0.001;
#pragma omp for OMP_SCHEDULE_AUTO
      for (kj = 0; kj < num_kj; ++kj)
      {
      	  knotv = vknots_begin + kj; // Left side of global element.

	  for (pp0=0; pp0<(int)points.size() && points[pp0+1] < knotv[0]; pp0+=3)
	  {
	      if (kj == 0)
	      {
		  all_classification[kj].push_back(-1);
	      }
	  }

	  for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < knotv[1]; pp1+=3);

	  if (knotv+2 == vknots_end)
	      for (; pp1<(int)points.size() && points[pp1+1] <= knotv[1]; pp1+=3);

	  // Sort the current sub set of points according to the u-parameter
	  qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

	  // Traverse the relevant points and identify the associated element
	  for (pp2=pp0, knotu=uknots_begin; pp2<pp1 && points[pp2] < (*knotu); pp2+=3)
	  {
	      all_classification[kj].push_back(-1);
	  }

	  for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	  {
	      for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	      if (knotu+1 == uknots_end)
	      {
		  for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	      }
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	      elem = elements[kj*(nmb_knots_u-1)+ki];

	      nump = (pp3 - pp2)/3;
	      for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	      {
		  // Evaluate
		  surf->point(pos, curr[0], curr[1], elem);
		  dist = curr[2]-pos[0];

		  if (evalsrf.get())
		    {
		      Point clo_pt;
		      double clo_u, clo_v, clo_dist;
		      double seed[2];
		      seed[0] = curr[0];
		      seed[1] = curr[1];
		      Point pt(curr[0], curr[1], curr[2]);
		      surf->setCurrentElement(elem);
		      evalsrf->closestPoint(pt, clo_u, clo_v, clo_pt,
					    clo_dist, aeps, 1, seed);
		      if (clo_dist < fabs(dist))
			dist = (dist < 0.0) ? -clo_dist : clo_dist;
		    }

		  all_max_above[kj] = std::max(all_max_above[kj], dist);
		  all_max_below[kj] = std::min(all_max_below[kj], dist);
		  all_dist[kj] += fabs(dist);
		  all_nmb_points[kj]++;
	  
		  // Find classification
		  for (ka=0; ka<(int)limits.size(); ++ka)
		      if (dist < limits[ka])
		      {
			  all_nmb_group[kj][ka]++;
			  all_classification[kj].push_back(ka);
			  break;
		      }
		  if (ka == (int)limits.size())
		  {
		      all_nmb_group[kj][ka]++;
		      all_classification[kj].push_back(ka);
		  }
	      }
	      pp2 = pp3;
	  }
	  for (; pp2<pp1; pp2+=3)
	  all_classification[kj].push_back(-1);
	  pp0 = pp1;

	  if (kj == num_kj - 1)
	  {
	      for (kr = pp0; kr < (int)points.size(); kr += 3)
	      {
		  all_classification[kj].push_back(-1);
	      }
	  }
      } // End of parallel for loop.
  } // End of parallel region.

  for (kj = 0; kj < num_kj; ++kj)
  {
      max_above = std::max(all_max_above[kj], max_above);
      max_below = std::min(all_max_below[kj], max_below);
      avdist += all_dist[kj];
      nmb_points += all_nmb_points[kj];
      classification.insert(classification.end(),
			    all_classification[kj].begin(), all_classification[kj].end());
      for (size_t ki = 0; ki < all_nmb_group[kj].size(); ++ki)
      {
	  nmb_group[ki] += all_nmb_group[kj][ki];
      }
  }

  avdist /= nmb_points;
}

//=============================================================================
double getSignedMaxDist(LRBSpline2D* bspline, int sgn)
//=============================================================================
{
  double res = 0.0;
  const int dim = bspline->dimension();
  int del = 3 + dim;
  const vector<Element2D*>& elements = bspline->supportedElements();
  for (size_t ki=0; ki<elements.size(); ++ki)
    {
      vector<double>& points = elements[ki]->getDataPoints();
      for (size_t kj=0; kj<points.size(); kj+=del)
	{
	  double dist = points[kj+del-1];
	  if (sgn*dist > 0.0)
	    {
	      if (sgn < 0)
		res = std::min(res, dist);
	      else
		res = std::max(res, dist);
	    }
	}
	     }
      return res;
}

//=============================================================================
void LRApproxApp::limitingSurfs(vector<double>& points,  // The points are modified!!!
				shared_ptr<ParamSurface>& surf,
				int nmb_iter,
				shared_ptr<ParamSurface>& limsf1,
				shared_ptr<ParamSurface>& limsf2)
//=============================================================================
{
  shared_ptr<BoundedSurface> bdsf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  shared_ptr<ParamSurface> parsf = surf;
  if (bdsf.get())
    {
      parsf = bdsf->underlyingSurface();
    }

  shared_ptr<LRSplineSurface> sf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(parsf);
 
  if (!sf.get())
   {
     THROW("Input surface is not of type LR B-spline:");
   }
 
  shared_ptr<LRSplineSurface> lim1, lim2;
  limitingSurfs(points, sf, nmb_iter, lim1, lim2);

  if (bdsf.get())
    {
      // The input surface is trimmed. Trim limit surfaces accordingly
      vector<CurveLoop> loops = bdsf->allBoundaryLoops();
      vector<CurveLoop> loops1(loops.size());
      vector<CurveLoop> loops2(loops.size());
      for (size_t ki=0; ki<loops.size(); ++ki)
	{
	  int nmb = loops[ki].size();
	  vector<shared_ptr<ParamCurve> > loop_cvs1(nmb);
	  vector<shared_ptr<ParamCurve> > loop_cvs2(nmb);
	  for (int kj=0; kj<nmb; ++kj)
	    {
	      shared_ptr<ParamCurve> curr = loops[ki][kj];
	      shared_ptr<CurveOnSurface> sfcv = 
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr);
	      if (!sfcv.get())
		{
		  THROW("Trimming curve of wrong type");
		}
	      shared_ptr<ParamCurve> tmp_par1(sfcv->parameterCurve()->clone());
	      shared_ptr<ParamCurve> cv1(new CurveOnSurface(lim1, 
							    tmp_par1,
							    true));
	      shared_ptr<ParamCurve> tmp_par2(sfcv->parameterCurve()->clone());
	      shared_ptr<ParamCurve> cv2(new CurveOnSurface(lim2, 
							    tmp_par2,
							    true));
	      loop_cvs1[kj] = cv1;
	      loop_cvs2[kj] = cv2;
	    }
	  double eps = std::max(1.0e-12, loops[ki].getSpaceEpsilon());
	  loops1[ki].setCurves(loop_cvs1);
	  loops1[ki].setSpaceEpsilon(eps);
	  loops2[ki].setCurves(loop_cvs2);
	  loops2[ki].setSpaceEpsilon(eps);
	}
      limsf1 = shared_ptr<ParamSurface>(new BoundedSurface(lim1, loops1));
      limsf2 = shared_ptr<ParamSurface>(new BoundedSurface(lim2, loops2));
    }
  else
    {
      limsf1 = lim1;
      limsf2 = lim2;
    }
 }

//=============================================================================
void LRApproxApp::limitingSurfs(vector<double>& points,  // The points are modified!!!
				shared_ptr<LRSplineSurface>& surf,
				int nmb_iter,
				shared_ptr<LRSplineSurface>& limsf1,
				shared_ptr<LRSplineSurface>& limsf2)
//=============================================================================
{
  // Define parameters
  int mba = 1;
  double eps = 1.0e-6;  // Not really used

  int dim = surf->dimension();
  if (dim != 1)
    THROW("Dimension different from 1 in computation of limit surfaces");

  // Create initial error limit surfaces having the same knots as the initial
  // surface and level value zero
  limsf1 = shared_ptr<LRSplineSurface>(new LRSplineSurface(*surf));
  limsf2 = shared_ptr<LRSplineSurface>(new LRSplineSurface(*surf));

  Point coef(dim);
  coef.setValue(0.0);
  for (LRSplineSurface::BSplineMap::const_iterator it1 = limsf1->basisFunctionsBegin();
       it1 != limsf1->basisFunctionsEnd(); ++it1)
    limsf1->setCoef(coef, it1->second.get());

  for (LRSplineSurface::BSplineMap::const_iterator it1 = limsf2->basisFunctionsBegin();
       it1 != limsf2->basisFunctionsEnd(); ++it1)
    limsf2->setCoef(coef, it1->second.get());

  // Distribute points to limit surfaces and elements
  // Get all knot values in the u-direction
  const double* const uknots = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  const int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  // Construct meshes of element pointers
  vector<Element2D*> elements0;
  vector<Element2D*> elements1;
  vector<Element2D*> elements2;
  surf->constructElementMesh(elements0);
  limsf1->constructElementMesh(elements1);
  limsf2->constructElementMesh(elements2);

  // Sort points in v-direction
  const int nmb_pts = (int)points.size()/3;    // Parameter value + height
  qsort(&points[0], nmb_pts, 3*sizeof(double), compare_v_par);
  int ki, kj, kr;
  double *curr;
  double dist;

  int pp0, pp1;
  Element2D* elem0 = NULL;
  Element2D* elem1 = NULL;
  Element2D* elem2 = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)points.size() && points[pp0+1] < (*knotv); 
       pp0+=3);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)points.size() && points[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)points.size();

      // Sort the current sub set of points according to the u-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant points and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && points[pp2] < (*knotu); pp2+=3);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && points[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated elements
	  elem0 = elements0[kj*(nmb_knots_u-1)+ki];
	  elem1 = elements1[kj*(nmb_knots_u-1)+ki];
	  elem2 = elements2[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&points[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      surf->point(pos, curr[0], curr[1], elem0);
	      dist = curr[2]-pos[0];
	      curr[2] = dist;  // Change to difference field setting
	      if (dist < 0.0)
		elem1->addDataPoints(points.begin()+pp2+3*kr, 
				     points.begin()+pp2+3*(kr+1),
				     3, false);
	      else if (dist > 0.0)
		elem2->addDataPoints(points.begin()+pp2+3*kr, 
				     points.begin()+pp2+3*(kr+1),
				     3, false);
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }

  // Perform approximation
  bool init_mba = true; //false;
  double mba_level = 0.0;
  vector<double> points_dummy;
  LRSurfApprox approx1(limsf1, points_dummy, eps, init_mba, mba_level, false, false);
  approx1.setUseMBA(true);
  approx1.addUpperConstraint(0.0);
  approx1.setMBAiter(nmb_iter);
  approx1.setMBASign(-1);

  int max_iter = 0; //1;
  double maxdist1, avdist1, avdist_total1; // will be set below
  int nmb_out_eps1;        // will be set below
  shared_ptr<LRSplineSurface> surf1;
  try {
  surf1 = approx1.getApproxSurf(maxdist1, avdist_total1, avdist1, nmb_out_eps1, max_iter);
      }
  catch (...)
    {
  std::cout << "ERROR: Failed creating limit surface ";
}

  LRSurfApprox approx2(limsf2, points_dummy, eps, init_mba, mba_level, false, false);
  approx2.setUseMBA(true);
  approx2.addLowerConstraint(0.0);
  approx2.setMBAiter(nmb_iter);
  approx2.setMBASign(1);

  double maxdist2, avdist2, avdist_total2; // will be set below
  int nmb_out_eps2;        // will be set below
  shared_ptr<LRSplineSurface> surf2;
  try {
  surf2 = approx2.getApproxSurf(maxdist2, avdist_total2, avdist2, nmb_out_eps2, max_iter);
      }
  catch (...)
    {
  std::cout << "ERROR: Failed creating limit surface ";
}

  LRSplineSurface::BSplineMap::const_iterator it1;
  LRSplineSurface::BSplineMap::const_iterator it2;
#ifdef DEBUG
  std::ofstream of1("limsf1_0.g2");
  std::ofstream of2("limsf2_0.g2");
  surf1->writeStandardHeader(of1);
  surf1->write(of1);
  surf2->writeStandardHeader(of2);
  surf2->write(of2);
  int nmb_basis1 = 0;
  int nmb_basis2 = 0;
  double max_dist1 = 0.0, max_dist2 = 0.0;

  double maxdiff0 = 0.0;
  double mindiff0 = std::numeric_limits<double>::max();
  double avdiff0 = 0.0;
  double maxd0 = 0.0;
  double mind0 = std::numeric_limits<double>::max();
  double avd0 = 0.0;
  int nmb_basis = surf1->numBasisFunctions();
  it1 = surf1->basisFunctionsBegin();
  it2 = surf2->basisFunctionsBegin();
  for (; it1 != surf1->basisFunctionsEnd(); ++it1, ++it2)
    {
      Point coef1 = it1->second->Coef();
      Point coef2 = it2->second->Coef();
      double diff = coef2[0] - coef1[0];
      Point greville = it1->second->getGrevilleParameter();
      Point pos1, pos2;
      surf1->point(pos1, greville[0], greville[1]);
      surf2->point(pos2, greville[0], greville[1]);
      double ptdist = pos1.dist(pos2);
      maxdiff0 = std::max(maxdiff0, diff);
      mindiff0 = std::min(mindiff0, diff);
      avdiff0 += diff;
      maxd0 = std::max(maxd0, ptdist);
      mind0 = std::min(mind0, ptdist);
      avd0 += ptdist;
    }
  avdiff0 /= (double)nmb_basis;
  avd0 /= (double)nmb_basis;
  std::cout << "Number of basis functions: " << nmb_basis << std::endl;
  std::cout << "Maximum width (coefficient): " << maxdiff0 << ", minimum width: " << mindiff0 << ", average width: " << avdiff0 << std::endl;
  std::cout << "Maximum width (greville point): " << maxd0 << ", minimum width: " << mind0 << ", average width: " << avd0 << std::endl;
#endif

  double fac = 1.0;
  shared_ptr<LRSplineSurface> levelsf1(new LRSplineSurface(*surf1));
  shared_ptr<LRSplineSurface> levelsf2(new LRSplineSurface(*surf2));

  it1 = levelsf1->basisFunctionsBegin();  
  it2 = surf1->basisFunctionsBegin();  
  for (;it1 != levelsf1->basisFunctionsEnd(); ++it1, ++it2)
    {
      double dist1 = getSignedMaxDist(it2->second.get(), -1);
      Point cf(1);
      cf[0] = dist1;
#ifdef DEBUG
      if (fabs(dist1) > 1.0e-10)
	nmb_basis1++;
      max_dist1 = std::min(dist1, max_dist1);
#endif
      levelsf1->setCoef(cf, it1->second.get());
    }

  it1 = levelsf2->basisFunctionsBegin();
  it2 = surf2->basisFunctionsBegin();
  for (; it1 != levelsf2->basisFunctionsEnd(); ++it1, ++it2)
    {
      double dist2 = getSignedMaxDist(it2->second.get(), 1);
#ifdef DEBUG
      if (fabs(dist2) > 1.0e-10)
	nmb_basis2++;
      max_dist2 = std::max(dist2, max_dist2);
#endif
      Point cf(1);
      cf[0] = dist2;
      levelsf2->setCoef(cf, it1->second.get());
    }

  surf1->addSurface(*levelsf1, fac);
  surf2->addSurface(*levelsf2, fac);
#ifdef DEBUG
  std::ofstream of3("limsf1_1.g2");
  std::ofstream of4("limsf2_1.g2");
  surf1->writeStandardHeader(of3);
  surf1->write(of3);
  surf2->writeStandardHeader(of4);
  surf2->write(of4);

  std::cout << "Number of updated coefficients, limit surface 1: " << nmb_basis1 << std::endl;
  std::cout << "Maximum distance, surface 1: " << max_dist1 << std::endl;
  std::cout << "Number of updated coefficients, limit surface 2: " << nmb_basis2 << std::endl;
  std::cout << "Maximum distance, surface 2: " << max_dist2 << std::endl;

  double maxdiff = 0.0;
  double mindiff = std::numeric_limits<double>::max();
  double avdiff = 0.0;
  double maxd1 = 0.0;
  double mind1 = std::numeric_limits<double>::max();
  double avd1 = 0.0;
  it1 = surf1->basisFunctionsBegin();
  it2 = surf2->basisFunctionsBegin();
  for (; it1 != surf1->basisFunctionsEnd(); ++it1, ++it2)
    {
      Point coef1 = it1->second->Coef();
      Point coef2 = it2->second->Coef();
      double diff = coef2[0] - coef1[0];
      Point greville = it1->second->getGrevilleParameter();
      Point pos1, pos2;
      surf1->point(pos1, greville[0], greville[1]);
      surf2->point(pos2, greville[0], greville[1]);
      double ptdist = pos1.dist(pos2);
      maxdiff = std::max(maxdiff, diff);
      mindiff = std::min(mindiff, diff);
      avdiff += diff;
      maxd1 = std::max(maxd1, ptdist);
      mind1 = std::min(mind1, ptdist);
      avd1 += ptdist;
    }
  avdiff /= (double)nmb_basis;
  avd1 /= (double)nmb_basis;
  std::cout << "Maximum width: " << maxdiff << ", minimum width: " << mindiff << ", average width: " << avdiff << std::endl;
  std::cout << "Maximum width (greville point): " << maxd1 << ", minimum width: " << mind1 << ", average width: " << avd1 << std::endl;
#endif

  surf1->addSurface(*surf, fac);
  surf2->addSurface(*surf, fac);

  limsf1 = surf1;
  limsf2 = surf2;
}

//=============================================================================
double getMaxDistSignificantPoints(LRBSpline2D* bspline, int& nmb_overlap)
//=============================================================================
{
  double res = 0.0;
  const int dim = bspline->dimension();
  int del = 3 + dim;
  const vector<Element2D*>& elements = bspline->supportedElements();
  double maxdist = 0.0;
  for (size_t ki=0; ki<elements.size(); ++ki)
    {
      vector<double>& points = elements[ki]->getSignificantPoints();
      for (size_t kj=0; kj<points.size(); kj+=del)
	{
	  double dist = points[kj+del-1];
	  if (fabs(dist) > maxdist)
	    {
	      maxdist = fabs(dist);
	      nmb_overlap = elements[ki]->nmbSupport();
	      res = dist;
	    }
	}
    }
  return res;
}

//=============================================================================
void  
LRApproxApp::updateSurfWithSignificantPts(shared_ptr<LRSplineSurface>& surf,
					  double tol, double tol_sign,
					  double fac_pos, double fac_neg,
					  double& maxdist, double& avdist,
					  double& avdist_out, int& nmb_out,
					  double& maxsign, double& avsign,
					  int& nmb_out_sign)
//=============================================================================
{
  if (surf->dimension() != 1)
    {
      MESSAGE("Surface dimension different from 1, no significant point update");
      return;
    }

  int nmb_changed = 0;
  double min_tol = 0.01;
  double min_tol_sign = min_tol*tol_sign/tol;
  const int dim = surf->dimension();
  // vector<LRBSpline2D*> to_change;
  // vector<Point> coef_change;
  LRSplineSurface::BSplineMap::const_iterator it1;
  for (it1=surf->basisFunctionsBegin();
       it1 != surf->basisFunctionsEnd(); ++it1)
    {
      // int nmb_bsplines = 1;
      // double dist = getMaxDistSignificantPoints(it1->second.get(), 
      // 						nmb_bsplines);
      const vector<Element2D*>& elements = it1->second->supportedElements();
      size_t ki;
      for (ki=0; ki<elements.size(); ++ki)
	{
	  int del = elements[ki]->getNmbValPrPoint();
	  if (del == 0)
	    del = dim+3;  // Parameter pair, point and distance
	  vector<double>& points = elements[ki]->getSignificantPoints();
	  size_t kj;
	  for (kj=0; kj<points.size(); kj+=del)
	    {
	      double dist = points[kj+dim+2];
	      double height = points[kj+dim+1];
	      double tol2 = (height < 0.0) ? tol_sign - fac_neg*height :
			    tol_sign + fac_pos*height;
	      if (fabs(dist) > tol2)
		{
		  // // Modify coefficient
		  // Point coef = it1->second->Coef();
		  // double diff = fabs(dist) - tol_sign/(double)nmb_bsplines;
		  // if (dist < 0)
		  //   diff *= -1;
		  
		  // Point coef2 = coef;
		  // coef2[0] += diff;
		  // to_change.push_back(it1->second.get());
		  // coef_change.push_back(coef2);
		  it1->second->setFixCoef(0);
		  ++nmb_changed;
		  break;
		}
	    }
	  if (kj < points.size())
	    break;
	}
      if (ki == elements.size())
	it1->second->setFixCoef(1);
    }

  // for (size_t kr=0; kr<to_change.size(); ++kr)
  //   surf->setCoef(coef_change[kr], to_change[kr]);
    
  if (nmb_changed > 0)
    LRSplineMBA::MBAUpdate(surf.get(), 1.0, 0.0, true);

  if (nmb_changed > 0)
    {
      // Update accuracy information
      // Initiate accuracy information
      maxdist = 0.0;
      avdist = 0.0;
      avdist_out = 0.0;
      nmb_out = 0;
      nmb_out_sign = 0;
      maxsign = 0.0;
      avsign = 0.0;

      LRSplineSurface::ElementMap::const_iterator it2;
      int num = surf->numElements();
      int ki, kj;
      double bval, sfval;

      int nmb_pts_all = 0, nmb_sign_all = 0;
      double dist;
      double *curr;
      for (it2=surf->elementsBegin(), kj=0; kj<num; ++it2, ++kj)
	{
	  if (!(it2->second->hasDataPoints() || 
		it2->second->hasSignificantPoints()))
	    {
	      // Reset accuracy information in element
	      it2->second->resetAccuracyInfo();
	      continue;   // No points in which to check accuracy
	    }

	  double av_dist = 0.0, max_dist = 0.0, acc_err = 0.0, acc_out = 0.0;
	  double av_dist_sign = 0.0, max_dist_sign = 0.0;
	  int nmb_out_el = 0, nmb_out_sign_el = 0;
	  double tol2;

	  vector<double>& points = it2->second->getDataPoints();
	  vector<double>& sign_points = it2->second->getSignificantPoints();
	  int nmb_pts = it2->second->nmbDataPoints();
	  int del = it2->second->getNmbValPrPoint();
	  if (del == 0)
	    del = dim+3;  // Parameter pair, point and distance
	  int nmb_sign = (int)sign_points.size()/del;
	  int nmb_all = nmb_pts + nmb_sign;

	  nmb_sign_all += nmb_sign;
	  nmb_pts_all += nmb_pts;

	  // Check if the accuracy can have changed
	  const vector<LRBSpline2D*>& bsplines = it2->second->getSupport();
	  size_t nb;
	  for (nb=0; nb<bsplines.size(); ++nb)
	    if (!bsplines[nb]->coefFixed())
	      break;

	  if (nb < bsplines.size())
	    {
	      // Recompute accuracy
	      for (ki=0, curr=(nmb_pts>0) ? &points[0] : &sign_points[0]; 
		   ki<nmb_all; ++ki)
		{
		  bool outlier = (del > dim+3 && curr[del-1] < 0.0);
		  sfval = 0.0;
		  for (size_t kr=0; kr<bsplines.size(); ++kr)
		    {
		      bsplines[kr]->evalpos(curr[0], curr[1], &bval);
		      sfval += bval;
		    }
	      
		  dist = curr[dim+1] - sfval;
		  curr[dim+2] = dist;
		  if (outlier)
		    nmb_pts_all--;
		  else
		    {
		      double dist2 = fabs(dist);
		      double height = curr[dim+1];
		      max_dist = std::max(max_dist, dist2);
		      acc_err += dist2;
		      tol2 = (ki < nmb_pts) ? tol : tol_sign;
		      tol2 = (height < 0.0) ? tol2 - fac_neg*height :
			tol2 + fac_pos*height;
		      tol = std::max(tol2, (ki<nmb_pts) ? min_tol : 
				     min_tol_sign);
		      if (ki >= nmb_pts)
			{
			  av_dist_sign += dist2;
			  max_dist_sign = std::max(max_dist_sign, dist2);
			}
		      if (dist2 > tol2)
			{
			  av_dist += dist2;
			  nmb_out_el++;
			  acc_out += (dist2 - tol2);
			  if (ki >= nmb_pts)
			    {
			      nmb_out_sign_el++;
			    }
			}
		    }
		  if (ki == nmb_pts-1 && nmb_all > nmb_pts)
		    curr = &sign_points[0];
		  else
		    curr += del;
		}
	      it2->second->setAccuracyInfo(acc_err, av_dist, max_dist,
					   nmb_out_el, nmb_out_sign_el,
					   acc_out);
	    }
	  else
	    {
	      it2->second->getAccuracyInfo(av_dist, max_dist, nmb_out_el, 
					   nmb_out_sign_el);
	      acc_err = it2->second->getAccumulatedError();
	      acc_out = it2->second->getAccumulatedOutside();

	      int nmb_outliers = it2->second->getNmbOutliers();
	      nmb_pts_all -= nmb_outliers;

	      int nmb_sign2;
	      it2->second->getInfoSignificantPoints(dim, max_dist_sign,
						    av_dist_sign, nmb_sign2);
	      av_dist *= nmb_out_el;
	      av_dist_sign *= nmb_sign;
	    }
	  maxdist = std::max(maxdist, max_dist);
	  avdist += acc_err;
	  avdist_out += av_dist;
	  nmb_out += nmb_out_el;
	  maxsign = std::max(maxsign, max_dist_sign);
	  avsign += av_dist_sign;
	  nmb_out_sign += nmb_out_sign_el;
	}
      avdist_out /= (double)nmb_out;
      avdist /= (double)(nmb_pts_all);
      avsign /= (double)(nmb_sign_all);
    }
}
