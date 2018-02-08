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

//#define DEBUG_ADAPT

#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/compositemodel/ftSmoothSurf.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"
#include "GoTools/compositemodel/ttlTriang.h"
#include "GoTools/compositemodel/ttlPoint.h"
#include "GoTools/compositemodel/cmUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/parametrization/PrPlanarGraph_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrOrganizedPoints.h"

using std::vector;
using std::list;


namespace Go
{
//===========================================================================
  shared_ptr<SplineSurface>
  AdaptSurface::approxInSplineSpace(shared_ptr<ParamSurface> surf,
				    shared_ptr<SplineSurface> surf2,
				    double tol)
//===========================================================================
  {
    shared_ptr<SplineSurface> updated_sf;

    // Compute correspondance of boundary curves.
    int idx = 0;  // Index of corner in surf corresponding to first 
    // corner in surf2
    bool turned = false; // Whether or not the boundary loop is turned
    bool correspond = getCornerCorrespondance(surf2, surf, idx, turned);
    if (!correspond)
      return updated_sf;  // Different number of corners

    // Get boundary loop
    CurveLoop loop = surf->outerBoundaryLoop(0);
    
    // Approximate the boundary curves given the corresponding spline spaces 
    // of the spline surface as the initial space
#ifdef DEBUG_ADAPT
    std::ofstream of0("cv0_out.g2");
#endif
    vector<shared_ptr<ParamCurve> > bd_cvs(4);
    for (int dir=0; dir<2; ++dir)
      {
	// Fetch spline space
	BsplineBasis basis = surf2->basis(dir);

	// Select boundary curves
	shared_ptr<ParamCurve> cvs[2];
	cvs[0] = shared_ptr<ParamCurve>(loop[(dir+idx)%4]->clone());
	cvs[1] = shared_ptr<ParamCurve>(loop[(dir+idx+2)%4]->clone());
	if (turned)
	  cvs[0]->reverseParameterDirection();
	else
	  cvs[1]->reverseParameterDirection();

#ifdef DEBUG_ADAPT
    for (int kr=0; kr<2; ++kr)
      {
	cvs[kr]->writeStandardHeader(of0);
	cvs[kr]->write(of0);
      }
#endif
	// Perform approximation
	vector<shared_ptr<SplineCurve> > app_cvs =
	  curveApprox(cvs, 2, basis, tol);
	bd_cvs[dir] = app_cvs[0];
	bd_cvs[dir+2] = app_cvs[1];
      }

    // Make initial spline surface as a Coons patch
    bd_cvs[2]->reverseParameterDirection();
    bd_cvs[3]->reverseParameterDirection();

 #ifdef DEBUG_ADAPT
   std::ofstream of("cv_out.g2");
    for (int kr=0; kr<4; ++kr)
      {
	bd_cvs[kr]->writeStandardHeader(of);
	bd_cvs[kr]->write(of);
      }
#endif

    // Just to be sure
    Point pos0 = bd_cvs[0]->point(bd_cvs[0]->startparam());
    Point pos1 = bd_cvs[0]->point(bd_cvs[0]->endparam());
    for (size_t kj=1; kj<bd_cvs.size(); ++kj)
      {
	Point pos2 = bd_cvs[kj]->point(bd_cvs[kj]->startparam());
	Point pos3 = bd_cvs[kj]->point(bd_cvs[kj]->endparam());
	if (kj == 1)
	  {
	    double d0 = std::min(pos0.dist(pos2), pos0.dist(pos3));
	    double d1 = std::min(pos1.dist(pos2), pos1.dist(pos3));
	    if (d0 < d1)
	      {
		// First curve has inconsistent orientation
		bd_cvs[0]->reverseParameterDirection();
		std::swap(pos0, pos1);
	      }
	  }
	      
	if (pos1.dist(pos2) < pos1.dist(pos3))
	  pos1 = pos3;
	else
	  {
	    bd_cvs[kj]->reverseParameterDirection();
	    pos1 = pos2;
	  }
      }
   
    CurveLoop boundary(bd_cvs, loop.getSpaceEpsilon());
    shared_ptr<SplineSurface> init_sf(CoonsPatchGen::createCoonsPatch(boundary));
    // Adapt the initial surface to the given input surface within the
    // given tolerance
    // The spline space of the initial surface may be refined
    updated_sf = adaptSurface(surf, init_sf, tol);
    
    return updated_sf;
  }

//===========================================================================
  vector<shared_ptr<SplineSurface> >
  AdaptSurface::expressInSameSplineSpace(shared_ptr<ParamSurface> surf1,
					 shared_ptr<ParamSurface> surf2,
					 double tol)
//===========================================================================
  {
    vector<shared_ptr<SplineSurface> > updated_sf;

    // Compute correspondance of boundary curves.
    int idx = 0;  // Index of corner in surf corresponding to first 
    // corner in surf2
    bool turned = false; // Whether or not the boundary loop is turned
    bool correspond = getCornerCorrespondance(surf1, surf2, idx, turned);
    if (!correspond)
      return updated_sf;  // Different number of corners

    // Get boundary loops
    CurveLoop loop1 = surf1->outerBoundaryLoop();
    CurveLoop loop2 = surf2->outerBoundaryLoop();
    
    // Approximate the boundary curves given the corresponding spline spaces 
    // of the spline surface as the initial space
#ifdef DEBUG_ADAPT
    std::ofstream of0("cv0_out.g2");
#endif
    vector<shared_ptr<ParamCurve> > bd_cvs1(4);
    vector<shared_ptr<ParamCurve> > bd_cvs2(4);
    for (int dir=0; dir<2; ++dir)
      {
	// Select boundary curves
	shared_ptr<ParamCurve> cvs[4];
	cvs[0] = shared_ptr<ParamCurve>(loop1[(dir)%4]->clone());
	cvs[1] = shared_ptr<ParamCurve>(loop1[(dir+2)%4]->clone());
	cvs[2] = shared_ptr<ParamCurve>(loop2[(dir+idx)%4]->clone());
	cvs[3] = shared_ptr<ParamCurve>(loop2[(dir+idx+2)%4]->clone());

	cvs[1]->reverseParameterDirection();
	if (turned)
	  cvs[2]->reverseParameterDirection();
	else
	  cvs[3]->reverseParameterDirection();

 #ifdef DEBUG_ADAPT
   for (int kr=0; kr<4; ++kr)
      {
	cvs[kr]->writeStandardHeader(of0);
	cvs[kr]->write(of0);
      }
#endif
	// Perform approximation
	vector<shared_ptr<SplineCurve> > app_cvs =
	  curveApprox(cvs, 4, tol);
	bd_cvs1[dir] = app_cvs[0];
	bd_cvs1[dir+2] = app_cvs[1];
	bd_cvs2[dir] = app_cvs[2];
	bd_cvs2[dir+2] = app_cvs[3];
      }

    // Make initial spline surface as a Coons patch
    bd_cvs1[2]->reverseParameterDirection();
    bd_cvs1[3]->reverseParameterDirection();
    bd_cvs2[2]->reverseParameterDirection();
    bd_cvs2[3]->reverseParameterDirection();

#ifdef DEBUG_ADAPT
    std::ofstream of("cv_out.g2");
    for (int kr=0; kr<4; ++kr)
      {
	bd_cvs1[kr]->writeStandardHeader(of);
	bd_cvs1[kr]->write(of);
	bd_cvs2[kr]->writeStandardHeader(of);
	bd_cvs2[kr]->write(of);
      }
#endif

    CurveLoop boundary1(bd_cvs1, loop1.getSpaceEpsilon());
    shared_ptr<SplineSurface> init_sf1(CoonsPatchGen::createCoonsPatch(boundary1));
    // Adapt the initial surface to the given input surface within the
    // given tolerance
    // The spline space of the initial surface may be refined
    shared_ptr<SplineSurface> update1_sf = adaptSurface(surf1, init_sf1, tol);
    updated_sf.push_back(update1_sf);
    
    CurveLoop boundary2(bd_cvs2, loop2.getSpaceEpsilon());
    shared_ptr<SplineSurface> init_sf2(CoonsPatchGen::createCoonsPatch(boundary2));
    // Adapt the initial surface to the given input surface within the
    // given tolerance
    // The spline space of the initial surface may be refined
    shared_ptr<SplineSurface> update2_sf = adaptSurface(surf2, init_sf2, tol);
    updated_sf.push_back(update2_sf);
    
    return updated_sf;
  }

//===========================================================================
  bool
  AdaptSurface::getCornerCorrespondance(shared_ptr<ParamSurface> surf1,
					shared_ptr<ParamSurface> surf2,
					int& idx, bool& turned)
//===========================================================================
  {
    // Fetch the corner points of the two surfaces
    vector<pair<Point, Point> > corners1;
    vector<pair<Point, Point> > corners2;

    surf1->getCornerPoints(corners1);
    surf2->getCornerPoints(corners2);
    if (corners1.size() != corners2.size())
      return false;  // Different number of corners

    // Compute distances between corresponding corners in all possible
    // corner configuarations.
    // The corner sequences are expected to loop around the surfaces, but
    // the sequences may be turned compared to each other and the position
    // of the first point in the sequence may differ
    int nmb = (int)corners1.size();
    vector<double> dists(2*nmb);
    int ki, kj, kr;
    double dd;
    for (ki=0; ki<nmb; ++ki)
      {
	for (kj=0, dd=0; kj<nmb; ++kj)
	  {
	    kr = (ki + kj) % nmb;
	    dd += corners1[kj].first.dist(corners2[kr].first);
	  }
	dists[ki] = dd;
      }

    for (ki=0; ki<nmb; ++ki)
      {
	for (kj=0, dd=0; kj<nmb; ++kj)
	  {
	    kr = (ki + kj) % nmb;
	    kr = nmb - kr - 1;
	    dd += corners1[kj].first.dist(corners2[kr].first);
	  }
	dists[nmb+ki] = dd;
      }

    // Find smallest distance
    double min_dist = dists[0];
    int min_idx = 0;
    for (ki=1; ki<2*nmb; ++ki)
      {
	if (dists[ki] < min_dist)
	  {
	    min_dist = dists[ki];
	    min_idx = ki;
	  }
      }
    idx = min_idx % nmb;
    turned = (min_idx >= nmb);
    return true;
   }

//===========================================================================
  vector<shared_ptr<SplineCurve> > 
  AdaptSurface::curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
			    const BsplineBasis& init_basis, double tol,
			    int nmb_init_sample_pr_seg)
//===========================================================================
  {
    return CurveCreators::curveApprox(cvs, nmb_cvs, init_basis, tol,
				      nmb_init_sample_pr_seg);
  }

//===========================================================================
  vector<shared_ptr<SplineCurve> > 
  AdaptSurface::curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
			    double tol, double degree,
			    int nmb_init_sample_pr_seg)
//===========================================================================
  {
    return CurveCreators::curveApprox(cvs, nmb_cvs, tol, degree,
				      nmb_init_sample_pr_seg);
  }

//===========================================================================
  void AdaptSurface::createTriangulation(shared_ptr<ParamSurface> surf, 
					 const RectDomain& dom,
					 shared_ptr<ftPointSet>& points, 
					 vector<int>& corner,
					 bool consider_joint, int nmb)
//===========================================================================
  {
    int nmb_sample = (nmb <= 0) ? 10 : nmb; //50 20; //10;// Number of pts to sample in one direction. 
    //double bd_fac = 1.5;
    getBoundaryData(surf, dom, nmb_sample, points, corner);
 
#ifdef DEBUG_ADAPT
    std::ofstream of2("points1.g2");
    points->write(of2);

    std::ofstream of2p("par1.g2");
    points->write2D(of2p);
#endif
    // Sample points in the inner of the surface
    getInnerData(surf, dom, nmb_sample, points, consider_joint);

#ifdef DEBUG_ADAPT
    std::ofstream of3("points2.g2");
    points->write(of3);

    std::ofstream of3p("par2.g2");
    points->write2D(of3p);
#endif
    // Complete point set topology
    updatePointTopology(surf, *points);
    
#ifdef DEBUG_ADAPT
    std::ofstream of4("points3.g2");
    points->write(of4);

    std::ofstream of4p("par3.g2");
    points->write2D(of4p);
#endif
  }

//===========================================================================
  shared_ptr<SplineSurface> 
  AdaptSurface::adaptSurface(shared_ptr<ParamSurface> surf, 
			     shared_ptr<SplineSurface> init_surf, double tol)
//===========================================================================
  {
    shared_ptr<SplineSurface> result;

#ifdef DEBUG_ADAPT
    std::ofstream of1("coons.g2");
    init_surf->writeStandardHeader(of1);
    init_surf->write(of1);
    surf->writeStandardHeader(of1);
    surf->write(of1);

    std::cout << "Init surf: " << init_surf->numCoefs_u() << ", ";
    std::cout << init_surf->numCoefs_v() << std::endl;
#endif

    // Get domain
    RectDomain dom;
    try {
      dom = surf->containingDomain();
    }
    catch (...)
      {
	dom = init_surf->containingDomain();
      }

    BoundingBox bbox = surf->boundingBox();
    double len = bbox.low().dist(bbox.high());
    double fac1 = 10000.0;
    double fac2 = 100.0;
    int nmb = (len/tol < fac2) ? 5 : ((len/tol < fac1) ? 15 : 30);

    // Sample points and create triangulation
    vector<int> corner;
    shared_ptr<ftPointSet> points = shared_ptr<ftPointSet>(new ftPointSet());
    createTriangulation(surf, dom, points, corner, false, nmb);

    // Make sure that the initial surface has a parameterization that
    // corresponds to the geometry
    double len_u, len_v;
    GeometryTools::estimateSurfaceSize(*init_surf, len_u, len_v);
    init_surf->setParameterDomain(0.0, len_u, 0.0, len_v);

     // Parameterize the sampling points
    double max_error1 = 1.0e8;
    try {
      max_error1 = parameterizePoints(init_surf, points, corner);
    }
    catch (...)
      {
	points->computeDist(init_surf);
	max_error1 = points->getMaxDist();
	result = init_surf;
#ifdef DEBUG_ADAPT
    std::cout << "Result surf: " << result->numCoefs_u() << ", ";
    std::cout << result->numCoefs_v() << ", max error: " << max_error1 << std::endl;
#endif
    return result;
      }

#ifdef DEBUG_ADAPT
    std::ofstream of5p("par4.g2");
    points->write2D(of5p);
#endif
    if (max_error1 < tol)
      {
	result = init_surf;
      }
    else
      {
	// Approximate point set, i.e. modify the initial surface with
	// respect to the point set
	int maxiter = 3; //2; //5;
	double max_error2, mean_error;
	result = doApprox(init_surf, maxiter, points, tol, max_error2, mean_error);
	

#ifdef DEBUG_ADAPT
	std::cout << "Approximation errors: " << max_error2 << " " << mean_error << std::endl;
#endif
	if (max_error2 > max_error1)
	  {
#ifdef DEBUG_ADAPT
	    std::cout << "Max init error: " << max_error1 << ". Interchanges surfaces" << std::endl;
#endif
	    result = init_surf;
	  }
	else
	  max_error1 = max_error2;
      }
    //result->setParameterDomain(0.0, 1.0, 0.0, 1.0);

#ifdef DEBUG_ADAPT
    std::cout << "Result surf: " << result->numCoefs_u() << ", ";
    std::cout << result->numCoefs_v() << ", max error: " << max_error1 << std::endl;
#endif
    return result;
  }

//===========================================================================
  double 
  AdaptSurface::parameterizePoints(shared_ptr<SplineSurface> init_surf,
				   shared_ptr<ftPointSet> points,
				   vector<int> corner)
//===========================================================================
  {
    double epsge = 1.0e-6; 

    // We must make sure that the ftPointSet has the neighbour
    // structure the way the parametrization code expects it.
    points->orderNeighbours();

    // Parametrize boundary points.
    // The parametrization routine seems to work best with the unit domain.
    double u1 = init_surf->startparam_u();
    double u2 = init_surf->endparam_u();
    double v1 = init_surf->startparam_v();
    double v2 = init_surf->endparam_v();
    init_surf->setParameterDomain(0.0, 1.0, 0.0, 1.0);
    // Parameterize
    PrPrmUniform par;
    PrParametrizeBdy bdy;
    shared_ptr<PrOrganizedPoints> op = shared_ptr<PrOrganizedPoints>(points);

    if (corner.size() < 4)
      {
	// Modify corner information to handle degenerate surfaces.
	CurveLoop loop = init_surf->outerBoundaryLoop(-1);

	// Identify degenerate boundaries
	int nmb = loop.size();
	vector<Point> mid;
	vector<double> len;
	for (int ki=0; ki<nmb; ++ki)
	  {
	    shared_ptr<ParamCurve> cv = loop[ki];
	    Point pt1 = cv->point(cv->startparam());
	    Point pt2 = cv->point(cv->endparam());
	    len.push_back(pt1.dist(pt2));
	    mid.push_back(0.5*(pt1+pt2));
	  }
	for (int ki=0; ki<nmb; ++ki)
	  for (int kj=ki+1; kj<nmb; ++kj)
	    if (len[kj] < len[ki])
	      {
		std::swap(len[ki], len[kj]);
		std::swap(mid[ki], mid[kj]);
	      }
      
	// Identify missing corners
	for (int ki=(int)corner.size(); ki<nmb; ++ki)
	  {
	    int ix = 0;
	    Vector3D pos = points->get3dNode(corner[0]);
	    Point pos2(pos[0], pos[1], pos[2]);
	    double mindist = mid[ki-(int)corner.size()].dist(pos2);
	    for (size_t kr=1; kr<corner.size(); ++kr)
	      {
		Vector3D pos = points->get3dNode(corner[kr]);
		Point pos2(pos[0], pos[1], pos[2]);
		double dist = mid[ki-(int)corner.size()].dist(pos2);
		if (dist < mindist)
		  {
		    ix = (int)kr;
		    mindist = dist;
		  }
	      }
	    corner.insert(corner.begin()+ix, corner[ix]);
	    ++ki;
	  }
	int stop_break = 1;
      }

    try {
      bdy.attach(op);
      //    bdy.setParamKind(PrUNIFBDY);
      double umin = init_surf->startparam_u();
      double umax = init_surf->endparam_u();
      double vmin = init_surf->startparam_v();
      double vmax = init_surf->endparam_v();
      bdy.parametrize(corner[0], corner[1], corner[2], corner[3], 
		      umin, umax, vmin, vmax);

      par.attach(op);
      par.setBiCGTolerance(0.00001);
      par.parametrize();
    } catch(...) {
      THROW("Parameterization failed");
    }
#ifdef DEBUG_ADAPT
    std::ofstream of5p("par_m.g2");
    points->write2D(of5p);
#endif

    double dist1 = points->reparBdy(init_surf, false/*true*/); // Second parametrization, we use existing uv.
    double dist2 = points->reparInnerPoints(init_surf, false/*true*/);

    // Reparameterize surface and points based on initial surface domain
    init_surf->setParameterDomain(u1, u2, v1, v2);
    
    double u_frac = u2 - u1;
    double v_frac = v2 - v1;
    for (int ki = 0; ki < (int)points->size(); ++ki) {
      double u = points->getU(ki) * u_frac;
      double v = points->getV(ki) * v_frac;
      if ((*points)[ki]->isOnBoundary()) { // If bnd point, make sure it stays on bnd.
	if (fabs(u - u1) < epsge)
	  u = u1;
	else if (fabs(u - u2) < epsge)
	  u = u2;
	if (fabs(v - v1) < epsge)
	  v = v1;
	else if (fabs(v - v2) < epsge)
	  v = v2;
      }
      points->setU(ki, u);
      points->setV(ki, v);
    }
    return std::max(dist1, dist2);
  }

//===========================================================================
  double
  AdaptSurface::projectPoints(shared_ptr<SplineSurface> surf, 
			      shared_ptr<ftPointSet> points)
//===========================================================================
  {
    double eps = 1e-13;
    double u, v, dist;
    Point pt(3);
    Point clo_pnt(3);
    double error = 0.0;

    int nmb = points->size();
    for (int ki=0; ki<nmb; ++ki)
      {
	ftSamplePoint *curr = (*points)[ki];
	pt.setValue(curr->getPoint().begin());
	if (curr->isOnBoundary())
	  surf->closestBoundaryPoint(pt, u, v, clo_pnt, dist, eps);
	else
	  surf->closestPoint(pt, u, v, clo_pnt, dist, eps);
	
	curr->setDist(dist);
	curr->setPar(Vector2D(u,v));
	error = std::max(error, dist);
      }
    return error;
  }

//===========================================================================
  shared_ptr<SplineSurface>
  AdaptSurface::doApprox(shared_ptr<SplineSurface> init_surf, int max_iter,
			 shared_ptr<ftPointSet> points, double tol,
			 double& max_error, double& mean_error)
//===========================================================================
  {
    bool isOK;
    double prevmax, prevmean;
    double prevapprox;

    max_error = points->getMaxDist();
    mean_error = points->getMeanDist();
#ifdef DEBUG_ADAPT
    std::cout << "-1: " << max_error << " " << mean_error << std::endl;
#endif

    // Make instance for surface modification and smoothing.
    double approx_orig_tol = -1.0; // We will not perform inner smoothing.
    vector<int> edge_derivs(4, 1);
    shared_ptr<SplineSurface> surf(init_surf->clone());
    int nmb_coef = surf->numCoefs_u()*surf->numCoefs_v();
    int nmb_pnts = points->size();
    double approxweight = 0.5*((double)(nmb_pnts/nmb_coef));
    approxweight = std::max(std::min(0.99, approxweight), 0.8);
    approxweight = 1.0 - 1.0e-6;  // TEST
    int max_iter_init = 6;
    ftSmoothSurf smoothsrf(surf, tol, approx_orig_tol, edge_derivs, 
			   max_iter_init);
    smoothsrf.setApproxWeight(approxweight);

    bool reparam = true; 
    isOK = smoothsrf.update(*points, tol, reparam);
 #ifdef DEBUG_ADAPT
    std::ofstream of0("adapt1.g2");
    surf->writeStandardHeader(of0);
    surf->write(of0);
#endif
   

  // Iterate to make sure that the surface is accurate enough.
    smoothsrf.getError(max_error, mean_error);
    prevmax = max_error;
    prevmean = mean_error;
    int iter=0;
#ifdef DEBUG_ADAPT
    std::cout << "0: " << max_error << " " << mean_error << std::endl;
#endif

  // Make a copy of the previous version of the surface.
    shared_ptr<SplineSurface> prev_surf;
    prev_surf = init_surf;

    approxweight = smoothsrf.getApproxWeight();
    while (!isOK && iter < max_iter) {
      iter++;

    // Add more degrees of freedom to the surface where the
    // error is too large.
      // TEST
      // if (iter > 1)
      // 	{
      try {
	smoothsrf.refineSurf(*points);
      }
      catch(...)
      {
	break;
      }
	// }

      // Update the surface with the new points and refined data
      // size.
      try {
	isOK = smoothsrf.update(*points, tol, reparam);
      }
      catch(...)
	{
	  break;
	}

      prevapprox = approxweight;
      approxweight = smoothsrf.getApproxWeight();
      smoothsrf.getError(max_error, mean_error);
#ifdef DEBUG_ADAPT
      std::cout << iter << ": " << max_error << " " << mean_error << std::endl;
#endif
      if (/*iter > 1 && */(!(max_error < 0.95*prevmax || mean_error < 0.95*prevmean) ||
		       (max_error >= 0.99*prevmax && (mean_error >= 0.5*prevmean ||
						      prevmean < tol))))
	break;
      prevmax = max_error;
      prevmean = mean_error;
      prev_surf  = shared_ptr<SplineSurface>(surf->clone());
    }

    if (iter > 0 && ((max_error >= 0.99*prevmax &&
		      !(mean_error < 0.95*prevmean && 
			max_error < 1.01*prevmax)) ||
		     (max_error < 1.01*prevmax && prevmean < tol)))
      {
	surf = prev_surf;  // Revert to the previous version
	approxweight = prevapprox;
	max_error = prevmax;
	mean_error = prevmean;
#ifdef DEBUG_ADAPT
	std::cout << "final: " << max_error << " " << mean_error << std::endl;
#endif
      }
  
#ifdef DEBUG_ADAPT
    std::ofstream of("adapt_surf.g2");
    prev_surf->writeStandardHeader(of);
    prev_surf->write(of);
    surf->writeStandardHeader(of);
    surf->write(of);
#endif
    return surf;
  }

  //===========================================================================
  void
  AdaptSurface::getBoundaryData(shared_ptr<ParamSurface> surf, 
				const RectDomain& dom,
				int nmb_sample,shared_ptr<ftPointSet> points, 
				vector<int>& corner)
//===========================================================================
  {
    // Bounding domain
    double pardom[4];
    pardom[0] = dom.vmin();
    pardom[1] = dom.umax();
    pardom[2] = dom.vmax();
    pardom[3] = dom.umin();

    double eps = 1.0e-6;  // Tolerance in closest point
    bool set_second = false;

    // Fetch all curves
    vector<CurveLoop> loops = 
      SurfaceTools::absolutelyAllBoundarySfLoops(surf, 
					       DEFAULT_SPACE_EPSILON); //surf->allBoundaryLoops();

    // To get a close to uniform distribution of points, estimate the
    // average edge length
    size_t ki, kr;
    int nmb_cvs = 0;
    double av_len = 0.0;
    vector<double> cv_len;
    for (ki=0; ki<loops.size(); ++ki)
      {
	int size = loops[ki].size();
	nmb_cvs += size;
	for (int kj=0; kj<size; ++kj)
	  {
	    double len = loops[ki][kj]->estimatedCurveLength();
	    av_len += len;
	    cv_len.push_back(len);
	  }
      }
    av_len /= (double)nmb_cvs;

    // For each boundary curve, evaluate an appropriate number of points
    // and add them to the point set
    shared_ptr<ftFaceBase> dummy;
    int bd = 1;  // Indicates outer boundary point
    vector<double> corner_ang;
    for (ki=0, kr=0; ki<loops.size(); ++ki)
      {
	PointIter prevpt = 0;
	PointIter firstpt = 0;
	int size = loops[ki].size();
	int kj;
	for (kj=0; kj<size; ++kj, ++kr)
	  {
	    shared_ptr<ParamCurve> curr = loops[ki][kj];
	    shared_ptr<CurveOnSurface> bd_cv =
	      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr);

	    // Compute number of points
	    int nsample = (int)(nmb_sample*cv_len[kr]/av_len);

	    double t1 = curr->startparam();
	    double t2 = curr->endparam();
	    double tdel = (t2 - t1)/(double)(nsample);
	    double tol = std::max(1.0e-7, 1.0e-5*(t2 - t1));
	    double tpar;
	    
	    // Evaluate the start point of a curve which is a corner point
	    Point pos = curr->point(t1);
	    corner.push_back((points->size() == 0) ? 0 :
			     points->lastAdded()->getIndex() + 1);

	    // Check tangent in corner
	    vector<Point> der1 = curr->point(t1, 1);
	    shared_ptr<ParamCurve> prev = loops[ki][(kj>0) ? kj-1 : size-1];
	    vector<Point> der2 = prev->point(prev->endparam(), 1);
	    double ang = der1[1].angle(der2[1]);
	    corner_ang.push_back(ang);
	    
	    Vector3D pnt3D(pos[0], pos[1], pos[2]);
	    Vector2D par;
	    if (bd_cv)
	      {
		Point face_par = bd_cv->faceParameter(t1);
		par = Vector2D(face_par[0], face_par[1]);
	      }
	    else if (size == 4)
	      {
		// Expecting rectangular surface
		par[kj%2] = (kj < 2) ? t1 : t2;
		par[1 - (kj%2)] = pardom[kj];
	      }
	    else
	      {
		// Closest point
		double clo_u, clo_v, clo_dist;
		Point clo_pt;
		surf->closestPoint(pos, clo_u, clo_v, clo_pt, clo_dist, eps);
		par = Vector2D(clo_u, clo_v);
	      }
	    shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, 
								      bd,
								      dummy,
								      par));
	    ftpnt->setPar(par);
	    PointIter latestpt = points->addEntry(ftpnt);
	    if (prevpt)
	      {
		prevpt->addNeighbour(latestpt);
		latestpt->addNeighbour(prevpt);
		if (set_second)
		  {
		    points->setSecond(latestpt);
		    set_second = false;
		  }
	      }
	    else
	      {
		firstpt = latestpt;
		if (ki == 0)
		  {
		    points->setFirst(firstpt);
		    set_second = true;
		  }
	      }
	    prevpt = latestpt;
	    
	    // Evaluate the inner points of the boundary curve
	    //int kh;
	    tpar = std::min(t1+tdel, curr->nextSegmentVal(t1, true, tol));
	    //for (kh=1, tpar=t1+tdel; kh<nsample; ++kh, tpar+=tdel)
	    while (tpar < t2)
	      {
		pos = curr->point(tpar);
		pnt3D = Vector3D (pos[0], pos[1], pos[2]);

		// Compute face parameter
		if (bd_cv)
		  {
		    Point face_par = bd_cv->faceParameter(tpar);
		    par = Vector2D(face_par[0], face_par[1]);
		  }
		else if (size == 4)
		  {
		    // Expecting rectangular surface
		    par[kj%2] = (kj < 2) ? tpar : t2 - tpar;
		    par[1 - (kj%2)] = pardom[kj];
		  }
		else
		  {
		    // Closest point
		    double clo_u, clo_v, clo_dist;
		    Point clo_pt;
		    surf->closestPoint(pos, clo_u, clo_v, clo_pt, clo_dist, eps);
		    par = Vector2D(clo_u, clo_v);
		  }
		shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, 
									  bd,
									  dummy,
									  par));
		ftpnt->setPar(par);
		latestpt = points->addEntry(ftpnt);
		prevpt->addNeighbour(latestpt);
		latestpt->addNeighbour(prevpt);
		prevpt = latestpt;

		tpar = std::min(tpar+tdel, curr->nextSegmentVal(tpar, true, tol));
	      }
	  }
	// Close loop
	prevpt->addNeighbour(firstpt);
	firstpt->addNeighbour(prevpt);
	bd = 2;  // Inner boundary points
      }
      
  }

  //===========================================================================
  void
  AdaptSurface::getInnerData(shared_ptr<ParamSurface> surf, 
			     const RectDomain& dom,
			     int nmb_sample, shared_ptr<ftPointSet> points,
			     bool consider_joint)
//===========================================================================
  {
    // Preparatory computations
    // Get estimated length of surface sides
    double len_u, len_v;
    try {
    surf->estimateSfSize(len_u, len_v);
    }
    catch (...)
      {
	len_u = len_v = 1.0;
      }

    int cv_dir = 0;
    int pt_dir = 1;

    // Check ratio of boundaries of underlying surface to decide which
    // parameter direction should give rise to iso-parametric curves
    shared_ptr<ParamSurface> tmp_sf = surf;
    shared_ptr<BoundedSurface> bd_sf = 
      dynamic_pointer_cast<BoundedSurface,ParamSurface>(tmp_sf);
    if (bd_sf.get())
      tmp_sf = bd_sf->underlyingSurface();
    RectDomain tmp_dom = tmp_sf->containingDomain();
    double len_u1, len_u2, len_v1, len_v2;
    GeometryTools::estimateIsoCurveLength(*tmp_sf, false, tmp_dom.umin(), len_u1);
    GeometryTools::estimateIsoCurveLength(*tmp_sf, false, tmp_dom.umax(), len_u2);
    GeometryTools::estimateIsoCurveLength(*tmp_sf, true, tmp_dom.vmin(), len_v1);
    GeometryTools::estimateIsoCurveLength(*tmp_sf, true, tmp_dom.vmax(), len_v2);
    double frac1 = std::min(len_u1, len_u2)/std::max(len_u1, len_u2);
    double frac2 = std::min(len_v1, len_v2)/std::max(len_v1, len_v2);
    if (frac1 > frac2)
      {
	cv_dir = 1;
	pt_dir = 0;
      }

    shared_ptr<ftFaceBase> dummy;
    int bd = 0;  // Indicates inner point

    // Number of points to sample in each parameter direction
    // int nmb_u = (int)(nmb_sample*len_u/len_v);
    // int nmb_v = (int)(nmb_sample*len_v/len_u);
    double fac = len_u/len_v;
    double len = (double)nmb_sample/sqrt(fac);
    int nmb_u = (int)(fac*len);
    int nmb_v = (int)(len);
    int min_samples = 1; //3;
    nmb_u = std::max(nmb_u, min_samples);
    nmb_v = std::max(nmb_v, min_samples);
    nmb_u = std::min(nmb_u, 2*nmb_sample);
    nmb_v = std::min(nmb_v, 2*nmb_sample);

    int nmb_cv = (cv_dir == 0) ? nmb_u : nmb_v;
    int nmb_pt = (pt_dir == 0) ? nmb_u : nmb_v;
    double pt_len = (cv_dir == 0) ? len_v : len_u;

    // Fetch constant parameter curves in the selected curve direction
    double u1 = (cv_dir == 0) ? dom.umin() : dom.vmin();
    double u2 = (cv_dir == 0) ? dom.umax() : dom.vmax();
    double udel = (u2 - u1)/(double)(nmb_cv+1);
    double tol1 = std::max(1.0e-7, 1.0e-5*(u2-u1));
    double par[2];
    //int ki, kj;
    size_t kr;
    if (consider_joint)
      par[cv_dir] = std::min(u1+udel, surf->nextSegmentVal(cv_dir, u1, true, tol1));
    else
      par[cv_dir] = u1+udel;
    //for (ki=0, par[cv_dir]=u1+udel; ki<nmb_cv; ++ki, par[cv_dir]+=udel)
    while (par[cv_dir] < u2)
      {
	vector<shared_ptr<ParamCurve> > crvs = 
	  surf->constParamCurves(par[cv_dir], cv_dir /*false*/);

// #ifdef DEBUG_ADAPT
// 	std::ofstream of("conscrv.g2");
// 	for (size_t ka=0; ka<crvs.size(); ++ka)
// 	  {
// 	    shared_ptr<CurveOnSurface> sf_cv = 
// 	      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(crvs[ka]);
// 	    if (sf_cv.get())
// 	      sf_cv->ensureSpaceCrvExistence(0.0001);
// 	    crvs[ka]->geometryCurve()->writeStandardHeader(of);
// 	    crvs[ka]->geometryCurve()->write(of);
// 	  }
// #endif
	if (crvs.size() == 0)
	  {
	    if (consider_joint)
	      par[cv_dir] = std::min(par[cv_dir]+udel, surf->nextSegmentVal(cv_dir, par[cv_dir], 
									    true, tol1));
	    else
	      par[cv_dir] += udel;
	    continue;  // Outside domain of surface
	  }

	// Distribute sampling points
	double av_len = 0.0;
	vector<double> cv_len(crvs.size());
	for (kr=0; kr<crvs.size(); ++kr)
	  {
	    double len = crvs[kr]->estimatedCurveLength();
	    av_len += len;
	    cv_len[kr] = len;
	  }
	double curr_len = av_len;
	av_len /= (double)crvs.size();

	// Evaluate sampling points
	int curr_nmb = (int)(nmb_pt*(curr_len/pt_len)) + 1;
	for (kr=0; kr<crvs.size(); ++kr)
	  {
	    int nmb = (int)(curr_nmb*cv_len[kr]/av_len);
	    nmb = std::max(nmb, min_samples);
	    double v1 = crvs[kr]->startparam();
	    double v2 = crvs[kr]->endparam();
	    double vdel = (v2 - v1)/(double)(nmb+1);
	    double tol2 = std::max(1.0e-7, 1.0e-5*(v2-v1));
	    if (consider_joint)
	      par[pt_dir] = std::min(v1+vdel, surf->nextSegmentVal(pt_dir, v1, true, tol2));
	    else
	      par[pt_dir] = v1 + vdel;
	    //for (kj=0, par[pt_dir]=v1+vdel; kj<nmb; ++kj, par[pt_dir]+=vdel)
	    while (par[pt_dir] < v2)
	      {
		Point pos = crvs[kr]->point(par[pt_dir]);
		Vector3D pnt3D(pos[0], pos[1], pos[2]);
		Vector2D pntpar(par[0], par[1]);
		shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, 
									  bd,
									  dummy,
									  pntpar));
		ftpnt->setPar(pntpar);
		points->addEntry(ftpnt);
		
		if (consider_joint)
		  par[pt_dir] = std::min(par[pt_dir]+vdel, surf->nextSegmentVal(pt_dir, par[pt_dir], 
								  true, tol2));
		else
		  par[pt_dir] += vdel;
	      }
	  }
	if (consider_joint)
	  par[cv_dir] = std::min(par[cv_dir]+udel, surf->nextSegmentVal(cv_dir, par[cv_dir], 
									true, tol1));
	else
	  par[cv_dir] += udel;
      }
  }

//===========================================================================
void AdaptSurface::updatePointTopology(shared_ptr<ParamSurface> surf, 
				       ftPointSet& points)
//===========================================================================
{
  int j, m;
    vector<ttlPoint*>  tri_points;

    // Parameter values are then mapped accordingly.
    RectDomain rect_dom = surf->containingDomain();
    RectDomain new_dom = cmUtils::geometricParamDomain(surf.get());
    for (j = 0; j < points.size(); ++j) {
      ftSurfaceSetPoint* sspnt = dynamic_cast<ftSurfaceSetPoint*>(points[j]);
      if (sspnt == 0) 
	THROW("Expecting ftSurfaceSetPoint");

      double u = sspnt->parValue(0)[0];
      double v = sspnt->parValue(0)[1];
      double new_u = (new_dom.umax() - new_dom.umin())*(u - rect_dom.umin())/
	(rect_dom.umax() - rect_dom.umin()) + new_dom.umin();
      double new_v = (new_dom.vmax() - new_dom.vmin())*(v - rect_dom.vmin())/
	(rect_dom.vmax() - rect_dom.vmin()) + new_dom.vmin();
      tri_points.push_back(new ttlPoint(sspnt, new_u, new_v));
    }

    hetriang::Triangulation triang;
    triang.createDelaunay(tri_points.begin(),
			  tri_points.end());

    for (j = 0; j < (int)tri_points.size(); ++j)
      delete tri_points[j]; // No more need for object.

    // We run through the triangulation, updating the structure 
    const list<hetriang::Edge*>& l_edges = triang.getLeadingEdges();
    // For the triangulation, we run through all triangles.
    list<hetriang::Edge*>::const_iterator leading_edge_it = l_edges.begin();
    while (leading_edge_it != l_edges.end()) 
      {
	hetriang::Edge* curr_edge = *leading_edge_it;
	shared_ptr<hetriang::Node> source_node = 
	  curr_edge->getSourceNode();
	shared_ptr<hetriang::Node> target_node = 
	  curr_edge->getTargetNode();

	for (m = 0; m < 3; ++m) 
	  {
	    std::vector<PointIter> neighbours =
	      source_node->pointIter()->getNeighbours();
	    for (j = 0; j < (int)neighbours.size(); ++j)
	      if (target_node->pointIter() == neighbours[j])
		break;
	    // If break was executed, connection already exists.
	    // If both nodes are on boundary we dont't make the points they
	    // are referring to neighbours (as all boundary points already
	    // have got their maximum of two boundary neighbours).
	    // We could have allowed two points on a subsurfaceboundary
	    // to be neighbours, but it would reault in a conflict when
	    // topology is to be used in the context of a graph.
	    if (j == (int)neighbours.size())
	      if (!((source_node->pointIter()->isOnSubSurfaceBoundary()
		     && target_node->pointIter()->isOnSubSurfaceBoundary())))
		// Add neighbour (in one direction at the time).
		(source_node->pointIter())->
		  addNeighbour(target_node->pointIter());
	    if (m == 2)
	      break;
	    curr_edge = curr_edge->getNextEdgeInFace();
	    source_node = curr_edge->getSourceNode();
	    target_node = curr_edge->getTargetNode();
	  }
	++leading_edge_it; // We iterate to next triangle.
      }
}

} // namespace Go

