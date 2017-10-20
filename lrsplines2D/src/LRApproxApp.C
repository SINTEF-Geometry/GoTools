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
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/Utils.h"
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
  double low = HUGE;
  double high = -HUGE;
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
      approx.addLowerConstraint(low - 0.1*(high-low));
      approx.addUpperConstraint(high + 0.1*(high-low));
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
      approx.addLowerConstraint(extent[4] - 0.1*(zrange));
      approx.addUpperConstraint(extent[5] + 0.1*(zrange));
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

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
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
					 vector<double>& pointsdist)
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
					     vector<double>& pointsdist)
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
  int ki, kj, kr, ka;
#pragma omp parallel default(none) private(ki, kj, kr, knotv, knotu, pp0, pp1) \
  shared(surf, points, num_pts, pts_dist, num_kj, elements)
  {
      Point pos;
      int nump;
      int pp2, pp3;
      Element2D* elem = NULL;
      double *curr;
      double dist;
#pragma omp for schedule(auto)
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
  for (kj = 0; kj < num_pts.size(); ++kj)
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
					vector<int>& nmb_group)
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
					    vector<int>& nmb_group)
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
    shared(surf, points, limits, elements, all_max_above, all_max_below, all_dist, all_nmb_points, all_level_points)
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
#pragma omp for schedule(auto)
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
					  vector<int>& nmb_group)
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
					      vector<int>& nmb_group)
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
    shared(surf, points, limits, elements, all_max_above, all_max_below, all_dist, all_nmb_points, all_classification, all_nmb_group)
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
#pragma omp for schedule(auto)
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
