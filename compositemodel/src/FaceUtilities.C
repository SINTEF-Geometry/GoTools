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

#include "GoTools/compositemodel/FaceUtilities.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/creators/ModifySurf.h"

using namespace Go;
using std::vector;
using std::pair;

//===========================================================================
  void
  FaceUtilities::getBoundaryData(ftSurface* face, int nmb_sample, 
				 vector<SamplePointData>& sample_points)
//===========================================================================
  {
    double angtol = 0.01;
    double curvtol = 0.01;

    // Fetch associated surface
    shared_ptr<ParamSurface> surf = face->surface();

    // Fetch all edges
    vector<shared_ptr<ftEdge> > edges = face->getAllEdges();

    // Compute average number of sampling points pr edge
    int av_sample = nmb_sample/(int)edges.size();

    // To get a close to uniform distribution of points, estimate the
    // average edge length
    size_t ki, kr;
    int nmb_cvs = (int)edges.size();
    double av_len = 0.0;
    vector<double> cv_len(edges.size());
    for (ki=0; ki<edges.size(); ++ki)
      {
	double len = edges[ki]->estimatedCurveLength(edges[ki]->tMin(),
						     edges[ki]->tMax());
	av_len += len;
	cv_len[ki] = len;
      }
    av_len /= (double)nmb_cvs;

    // For each boundary curve, evaluate an appropriate number of points
    // and add them to the point set
    for (ki=0; ki<edges.size(); ++ki)
      {
	// Compute number of points
	int nsample = (int)(av_sample*cv_len[ki]/av_len+1);
	nsample = std::max(nsample, 3);

	double t1 = edges[ki]->tMin();
	double t2 = edges[ki]->tMax();
	double tdel = (t2 - t1)/(double)(nsample);
	double tpar;
	    
	// Evaluate the start point of a curve which is a corner point
	Point pos = edges[ki]->point(t1);
	Point par = edges[ki]->faceParameter(t1);
	Point normal = face->normal(par[0], par[1]);
	Point prev_par = par;

	double Kcurv, Hcurv;
	CurvatureAnalysis::curvatures(*surf, par[0], par[1], Kcurv, Hcurv);

	// Fetch vertex
	shared_ptr<Vertex> vx1 = edges[ki]->getVertex(true);

	// Get all associated faces and compute normals and curvatures
	vector<pair<ftSurface*, Point> > faces_par = vx1->getFaces();
	bool normal_OK = true;
	bool curv_OK = true;
	for (kr=0; kr<faces_par.size(); ++kr)
	  {
	    if (faces_par[kr].first == face)
	      continue;  // Already computed

	    Point pos2 = faces_par[kr].first->point(faces_par[kr].second[0],
						    faces_par[kr].second[1]);
	    Point normal2 = faces_par[kr].first->normal(faces_par[kr].second[0],
							faces_par[kr].second[1]);
	    
	    double Kcurv2, Hcurv2;
	    CurvatureAnalysis::curvatures(*faces_par[kr].first->surface(), 
		       faces_par[kr].second[0], faces_par[kr].second[1],
		       Kcurv2, Hcurv2);

	    // Check if the normal and curvature information is unique
	    double ang = normal.angle(normal2);
	    if (ang > angtol)
	      normal_OK = false;

	    if (fabs(Hcurv-Hcurv2) > curvtol)
	      curv_OK = false;
	  }

	if (!normal_OK)
	  normal.setValue(MAXDOUBLE);
	if (!curv_OK)
	  Hcurv = MAXDOUBLE;

	sample_points.push_back(SamplePointData(pos, normal, Hcurv,
						face, par[0], par[1],
						edges[ki].get(), t1));
	    
	

	// Evaluate the inner points of the edge
	int kh;
	Point dummy_vec;
	for (kh=1, tpar=t1+tdel; kh<nsample; ++kh, tpar+=tdel)
	  {
	    pos = edges[ki]->point(tpar);
	    par = edges[ki]->faceParameter(tpar, prev_par.begin());
	    normal = face->normal(par[0], par[1]);
	    CurvatureAnalysis::curvatures(*surf, par[0], par[1], Kcurv, Hcurv);

	    // Check adjacent surface
	    if (edges[ki]->twin())
	      {
		double clo_u, clo_v, clo_dist, clo_par;
		Point clo_pos;
		ftSurface *face2 = 
		  edges[ki]->twin()->geomEdge()->face()->asFtSurface();
		//ftEdgeBase *edge2 = 
		(void)face2->closestBoundaryPoint(pos, dummy_vec, clo_u, clo_v, 
						  clo_pos, clo_dist, clo_par);

		Point normal2 = face2->normal(clo_u, clo_v);
		double Kcurv2, Hcurv2;
		CurvatureAnalysis::curvatures(*face2->surface(), 
			   clo_u, clo_v, Kcurv2, Hcurv2);
		double ang = normal.angle(normal2);
		if (ang > angtol)
		  normal.setValue(MAXDOUBLE);

		if (fabs(Hcurv-Hcurv2) > curvtol)
		  Hcurv = MAXDOUBLE;
	      }
	    sample_points.push_back(SamplePointData(pos, normal, Hcurv,
						    face, par[0], par[1],
						    edges[ki].get(), t1));
	  }
      }
  }
      

  //===========================================================================
  void
  FaceUtilities::getInnerData(ftSurface* face, int nmb_sample_u, 
			      int nmb_sample_v, 
			      vector<SamplePointData>& sample_points)
//===========================================================================
  {
    // Fetch associated surface
    shared_ptr<ParamSurface> surf = face->surface();
    shared_ptr<ParamSurface> sf;

    // Check for bounded surface
    double tol2d = 1.0e-4;
    bool use_domain = true;
    shared_ptr<BoundedSurface> bd_surf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    if (bd_surf.get())
      {
	// A trimmed surface is found
	// Get underlying surface 
	sf = bd_surf->underlyingSurface();
	if (bd_surf->isIsoTrimmed(tol2d))
	  {
	    RectDomain domain = bd_surf->containingDomain();
	    RectDomain dom2 = sf->containingDomain();
	    double umin = std::max(domain.umin(), dom2.umin());
	    double umax = std::min(domain.umax(), dom2.umax());
	    double vmin = std::max(domain.vmin(), dom2.vmin());
	    double vmax = std::min(domain.vmax(), dom2.vmax());
    
	    vector<shared_ptr<ParamSurface> > sfs = sf->subSurfaces(umin, vmin, umax, vmax);
	    sf = sfs[0];
	  }
	else
	  use_domain = false;
      }
    else 
      sf = surf;

    RectDomain domain = sf->containingDomain();
    if (use_domain)
      {
	double udel = (domain.umax() - domain.umin())/(double)(nmb_sample_u-1);
	double vdel = (domain.vmax() - domain.vmin())/(double)(nmb_sample_v-1);
	int ki, kj;
	double upar, vpar;
	for (kj=0, vpar=domain.vmin()+vdel; kj<nmb_sample_v-2; ++kj, vpar+=vdel)
	  for (ki=0, upar=domain.umin()+udel; ki<nmb_sample_u-2; ++ki, upar+=udel)
	    {
	      vector<Point> der(3);
	      sf->point(der, upar, vpar, 1);
	      Point normal = der[1].cross(der[2]);

	      double Kcurv, Hcurv;
	      CurvatureAnalysis::curvatures(*sf, upar, vpar, Kcurv, Hcurv);

	      sample_points.push_back(SamplePointData(der[0], normal, Hcurv,
						      face, upar, vpar));
	    }
      }
    else
      {
	// Fetch constant parameter curves in the u-direction
	int min_samples = 1;
	double u1 = domain.umin();
	double u2 = domain.umax();
	double udel = (u2 - u1)/(double)(nmb_sample_u-1);
	double upar;
	int ki, kj;
	size_t kr;
	for (ki=0, upar=u1+udel; ki<nmb_sample_u-2; ++ki, upar+=udel)
	  {
	    vector<shared_ptr<ParamCurve> > crvs = 
	      surf->constParamCurves(upar, false);
	    if (crvs.size() == 0)
	      continue;  // Outside domain of surface

	    // Distribute sampling points
	    double v1 = domain.vmin();
	    double v2 = domain.vmax();
	    double av_len = 0.0;
	    vector<double> cv_len(crvs.size());
	    double par_len = 0.0;
	    for (kr=0; kr<crvs.size(); ++kr)
	      {
		double len = crvs[kr]->estimatedCurveLength();
		av_len += len;
		cv_len[kr] = len;
		par_len += (crvs[kr]->endparam() - crvs[kr]->startparam());
	      }
	    av_len /= (double)crvs.size();

	    // Evaluate sampling points
	    int curr_nmb = (int)((nmb_sample_v-2)*par_len/(v2 - v1) + 1);
	    for (kr=0; kr<crvs.size(); ++kr)
	      {
		int nmb = (int)(curr_nmb*cv_len[kr]/av_len);
		nmb = std::max(nmb, min_samples);
		v1 = crvs[kr]->startparam();
		v2 = crvs[kr]->endparam();
		double vdel = (v2 - v1)/(double)(nmb+1);
		double vpar;
		for (kj=0, vpar=v1+vdel; kj<nmb; ++kj, vpar+=vdel)
		  {
		    vector<Point> der(3);
		    sf->point(der, upar, vpar, 1);
		    Point normal = der[1].cross(der[2]);

		    double Kcurv, Hcurv;
		    CurvatureAnalysis::curvatures(*sf, upar, vpar, Kcurv, Hcurv);
		    
		    sample_points.push_back(SamplePointData(der[0], normal, 
							    Hcurv, face, 
							    upar, vpar));
		  }

	      }
	  }
      }
  }

//===========================================================================
bool
FaceUtilities::enforceCoLinearity(ftSurface *face1, ftEdge *edge1,
				  ftSurface *face2, 
				  double tol, double ang_tol)
//===========================================================================
{
  // Check input
  if ((!face1->isSpline()) || (!face2->isSpline()))
    return false;  // Associated surfaces are not spline surfaces, cannot
                   // modify coefficients

  if (!edge1->hasConnectivityInfo())
    return false;  // No tangency information

  // Check tangency
  int status = edge1->getConnectivityInfo()->WorstStatus();
  if (status > 1)
    return false;  // A corner curve

  // Fetch information about common boundary
  AdjacencyInfo adj_info = face1->getAdjacencyInfo(edge1, face2, tol);
  if (adj_info.adjacency_found_ == false)
    return false;

    // Fetch surface geometry
  shared_ptr<ParamSurface> srf1 = face1->surface();
  shared_ptr<ParamSurface> srf2 = face2->surface();
  shared_ptr<SplineSurface> splsf1 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
  shared_ptr<SplineSurface> splsf2 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
  if (!splsf1.get() || !splsf2.get())
    return false;

  // Check if the surface coefficients at the common boundary are almost
  // co-linear and fetch the local enumeration of the associated coefficients
  vector<vector<int> > coef_enum;
  int colinear = SurfaceTools::checkCoefCoLinearity(splsf1, splsf2, adj_info.bd_idx_1_, 
				      adj_info.bd_idx_2_, adj_info.same_orient_,
				       tol, ang_tol, coef_enum);
  if (coef_enum.size() == 0)
    return false;  // No information computed
  if (colinear == 0)
    return false;

  // Peform surface modification
  bool smoothed = ModifySurf::enforceCoefCoLinearity(splsf1, adj_info.bd_idx_1_, 
						     splsf2, adj_info.bd_idx_2_, 
						     tol, coef_enum);
  if (!smoothed)
    return false;

  // Update edge information if possible
  (void)edge1->updateEdgeInfo(tol);
  ftEdge *edge2 = edge1->twin()->geomEdge();
  if (edge2)
    (void)edge2->updateEdgeInfo(tol);
  
  return true;
}

//===========================================================================
bool
FaceUtilities::enforceVxCoLinearity(shared_ptr<Vertex> vx, 
				    double tol, double ang_tol)
//===========================================================================
{
  // Fetch all associated faces
  vector<pair<ftSurface*, Point> > faces = vx->getFaces();

  if (faces.size() != 4)
    return false;  // Not a regular corner

  vector<ftEdge*> edges;
  vector<Point> norms(faces.size());
  vector<shared_ptr<SplineSurface> > sfs(faces.size());
  vector<int> vx_enumeration(faces.size());
  size_t ki, kj;
  for (ki=0; ki<faces.size(); ++ki)
    {
      if (!faces[ki].first->isSpline())
	return false;   // Cannot modify

      shared_ptr<ParamSurface> curr_sf = faces[ki].first->surface();
      shared_ptr<SplineSurface> splsf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(curr_sf);
      if (!splsf.get())
	return false;  // Something wrong
      sfs[ki] = splsf;

      vector<ftEdge*> curr_edges = vx->getFaceEdges(faces[ki].first);
      if (curr_edges.size() != 2)
	return false;  // Unclear how to handle this

      if (!(curr_edges[0]->twin() && curr_edges[1]->twin()))
	return false;  // At the face set boundary. Do not change
      
      // Collect face normal
      norms[ki] = faces[ki].first->normal(faces[ki].second[0], faces[ki].second[1]);

      // Fetch coefficient enumeration of surface corner related to vertex
      // First fetch edge enumeration
      int idx1 = curr_edges[0]->getCurveIndex();
      int idx2 = curr_edges[1]->getCurveIndex();
      if (idx1 < 0 || idx2 < 0)
	return false;   // Edges meeting in vertex are not iso trimmed boundaries
      bool vx_enum_found = SurfaceTools::getCornerCoefEnum(splsf, idx1, idx2, vx_enumeration[ki]);
      if (!vx_enum_found)
	return false;
      
      edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
    }
  
  // Check tangent plane continuity at vertex
  for (ki=1; ki<norms.size(); ++ki)
    if (norms[0].angle(norms[ki]) > ang_tol)
      return false;
  
  // Compute adjacency info
  vector<pair<vector<int>, pair<int,int> > > coef_cond;
  for (ki=0; ki<faces.size(); ++ki)
    {
      ftEdge* twin1 = edges[2*ki]->twin()->geomEdge();
      ftEdge* twin2 = edges[2*ki+1]->twin()->geomEdge();
      if (!(twin1 && twin2))
	return false;

      for (kj=ki+1; kj<faces.size(); ++kj)
	{
	  ftEdge *curr_edge = NULL;
	  if (twin1 == edges[2*kj] || twin1 == edges[2*kj+1])
	    curr_edge = edges[2*ki];
	  else if (twin2 == edges[2*kj] || twin2 == edges[2*kj+1])
	    curr_edge = edges[2*ki+1];
	  if (curr_edge)
	    {
	      AdjacencyInfo adj_info = 
		faces[ki].first->getAdjacencyInfo(curr_edge, faces[kj].first, tol);
	      vector<vector<int> > coef_enum;
	      if (adj_info.adjacency_found_ == true)
		{
		  int colinear = SurfaceTools::checkCoefCoLinearity(sfs[ki], sfs[kj], 
						       adj_info.bd_idx_1_, 
		  				       adj_info.bd_idx_2_, 
		  				       adj_info.same_orient_,
		  				       tol, ang_tol, coef_enum);
		  if (colinear > 0)
		    {
		      for (size_t kr=0; kr<coef_enum.size(); ++kr)
			coef_cond.push_back(make_pair(coef_enum[kr], 
						      std::make_pair(ki,kj)));
		    }
		}
	    }
	}
    }

  // Enforce colinearity
  bool smoothed = ModifySurf::enforceVxCoefCoLinearity(sfs, vx_enumeration, coef_cond,
						       tol);

  return smoothed;
}



