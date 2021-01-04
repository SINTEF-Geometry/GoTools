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

#include "GoTools/lrsplines2D/LRMinMax.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ExtremalPoint.h"
#include "GoTools/utils/BoundingBox.h"
#include <fstream>
#include <iostream>

//#define DEBUG

using std::vector;
using std::pair;
using std::make_pair;
using namespace Go;

namespace {
  void addContourLoop(shared_ptr<ParamCurve> crv, double isoval,
		     const CurveLoop& loop, double epsge,
		      vector<pair<vector<shared_ptr<ParamCurve> >, double> >& crv_loops,
		     vector<BoundingBox>& bbox);

  void extractInnerCurves(vector<pair<vector<shared_ptr<ParamCurve> >, double> >& cvs,
			  vector<BoundingBox>& bbox,
			  double eps);
}; // end anonymous namespace 

//===========================================================================
void 
LRMinMax::computeMinMaxPoints(shared_ptr<ParamSurface> surface,
			      vector<pair<shared_ptr<ParamCurve>, double> >& contour_crvs,
			      double tol, double epsge,
			      vector<pair<Point, Point> >& minpoints,
			      vector<pair<Point, Point> >& maxpoints)
//===========================================================================
{
  // Check for height surface
  if (surface->dimension() != 1)
    THROW("Surface dimension different from 1");

  // Create curve loop 
  CurveLoop loop = SurfaceTools::outerBoundarySfLoop(surface, epsge);

  // Domain of surface
  RectDomain dom = surface->containingDomain();

  // Fetch parameter curves
  int nmb = loop.size();
  vector<shared_ptr<ParamCurve> > par_cvs(nmb);
  for (int ka=0; ka<nmb; ++ka)
    {
      if (loop[ka]->dimension() == 2)
	par_cvs[ka] = loop[ka];
      else
	{
	  shared_ptr<CurveOnSurface> sf_cv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loop[ka]);
	  if (!sf_cv.get())
	    {
	      THROW("Parameter curve not available");
	    }

	  if (!sf_cv->hasParameterCurve())
	    {
	      sf_cv->ensureParCrvExistence(tol);
	    }
	  par_cvs[ka] = sf_cv->parameterCurve();
	  if (!par_cvs[ka].get())
	    {
	      THROW("Parameter loop not available");
	    }
	}
    }
  CurveLoop parloop(par_cvs, epsge);

  shared_ptr<ParamSurface> surf;
  shared_ptr<BoundedSurface> bd_surf = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surface);
  if (bd_surf.get())
    surf = bd_surf->underlyingSurface();
  else
    surf = surface;

  // Ensure closed contour loops
  vector<pair<vector<shared_ptr<ParamCurve> >,double> > cvs;
  vector<BoundingBox> bbox;
  for (size_t ki=0; ki<contour_crvs.size(); ++ki)
    {
      addContourLoop(contour_crvs[ki].first, contour_crvs[ki].second, 
		     parloop, epsge, cvs, bbox);
    }

  // Ensure consistent direction of curves in loop
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      if (cvs[ki].first.size() == 1)
	continue;

      Point pos1 = cvs[ki].first[0]->point(cvs[ki].first[0]->startparam());
      Point pos2 = cvs[ki].first[0]->point(cvs[ki].first[0]->endparam());
      Point pos3 = cvs[ki].first[1]->point(cvs[ki].first[1]->startparam());
      if (pos3.dist(pos1) < pos3.dist(pos2))
	cvs[ki].first[0]->reverseParameterDirection();
    }

  // For each combination of curves, check if one curve can possibly lie 
  // inside the other. If so, check if it is the case. Curves that are
  // external to another curve is removed.
  extractInnerCurves(cvs, bbox, epsge);

#ifdef DEBUG
  std::ofstream of("inner_cvs.g2");
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      for (size_t kj=0; kj<cvs[ki].first.size(); ++kj)
	{
	  cvs[ki].first[kj]->writeStandardHeader(of);
	  cvs[ki].first[kj]->write(of);
	}
    }
#endif

  // Make sure that the curves are oriented counter clockwise and
  // create trimmed surfaces for each inner curve loop.
  // Sort according to whether they contain minimum points or maximum
  // points.
#ifdef DEBUG
  std::ofstream ofmx("max_sfs.g2");
  std::ofstream ofmn("min_sfs.g2");
#endif
  int nmb_search = 0;
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      bool done = LoopUtils::makeLoopCCW(cvs[ki].first, epsge);

      // Pick part of surface containing loop
      // First make sure that the sub surface domain does not exceed the
      // domain of the surface
      double umin = std::max(bbox[ki].low()[0], dom.umin());
      double umax = std::min(bbox[ki].high()[0], dom.umax());
      double vmin = std::max(bbox[ki].low()[1], dom.vmin());
      double vmax = std::min(bbox[ki].high()[1], dom.vmax());
      if (umax-umin < tol)
	{
	  umax = std::min(dom.umax(), umax+tol);
	  umin = std::max(dom.umin(), umin-tol);
	}
      if (vmax - vmin < tol)
	{
	  vmax = std::min(dom.vmax(), vmax+tol);
	  vmin = std::max(dom.vmin(), vmin-tol);
	}
      if (umax < umin || vmax < vmin)
	THROW("Mismatch between curves and surface");
      
      vector<shared_ptr<ParamSurface> > sub_sfs =
	surf->subSurfaces(umin, vmin, umax, vmax, epsge);

      // Create trimmed surface
      // First make CurveOnSurfaceCurves
      vector<shared_ptr<CurveOnSurface> > trim_cvs(cvs[ki].first.size());
      for (size_t kj=0; kj<cvs[ki].first.size(); ++kj)
	trim_cvs[kj] = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(sub_sfs[0],
							cvs[ki].first[kj],
							true));

      shared_ptr<BoundedSurface> sub_bd(new BoundedSurface(sub_sfs[0],
							   trim_cvs, epsge));

      // Check surface configuration
      double upar, vpar;
      Point innerpt = sub_bd->getInternalPoint(upar, vpar);
      double tpar = 0.5*(cvs[ki].first[0]->startparam()+
			 cvs[ki].first[0]->endparam());
      Point bdpar = cvs[ki].first[0]->point(tpar);
      Point bdpt = sub_bd->ParamSurface::point(bdpar[0], bdpar[1]);
      int sgn = 0;
      if (innerpt[0] >= bdpt[0])
	{
#ifdef DEBUG
	  sub_bd->writeStandardHeader(ofmx);
	  sub_bd->write(ofmx);
#endif
	  sgn = 1;
	}
      else
	{
#ifdef DEBUG
	  sub_bd->writeStandardHeader(ofmx);
	  sub_bd->write(ofmx);
#endif
	  sgn = -1;
	}

	// Find global extremal points on sub surface
	vector<pair<Point, Point> > extpoints;
	int nmb_tp = 
	  computeExtremalPoints(sub_bd, sgn, tol, epsge, extpoints);
	nmb_search += nmb_tp;
	if (sgn < 0)
	  minpoints.insert(minpoints.end(), extpoints.begin(), extpoints.end());
	else
	  maxpoints.insert(maxpoints.end(), extpoints.begin(), extpoints.end());
	int stop_break0 = 1;
    }

#ifdef DEBUG
  std::cout << "Number of searches: " << nmb_search << std::endl;
#endif
  int stop_break = 1;
}

//===========================================================================
int LRMinMax::computeExtremalPoints(shared_ptr<ParamSurface> surface,
				    int sgn, double tol, double epsge,
				    vector<pair<Point, Point> >& extpoints)
//===========================================================================
{
  // Check for height surface
  if (surface->dimension() != 1)
    THROW("Surface dimension different from 1");

  shared_ptr<LRSplineSurface> lrsurf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(surface);
  shared_ptr<BoundedSurface> bdsurf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
  CurveBoundedDomain bddom;
  CurveLoop loop = SurfaceTools::outerBoundarySfLoop(surface, epsge);
  int nmb = loop.size();

  if (bdsurf.get())
    {
      lrsurf = dynamic_pointer_cast<LRSplineSurface, ParamSurface>(bdsurf->underlyingSurface());
      bddom = bdsurf->parameterDomain();
    }

  if (!lrsurf.get())
    {
      THROW("No LR B-spline surface is found");
    }
#ifdef DEBUG
  std::ofstream of1("curr_sf.g2");
  surface->writeStandardHeader(of1);
  surface->write(of1);

  std::ofstream of2("curr_el.g2");
  LineCloud lines = lrsurf->getElementBds();
  lines.writeStandardHeader(of2);
  lines.write(of2);
#endif

  //break LR surface up into individual patches
  int threshold_missing = 100; //50;
  const vector<pair<shared_ptr<LRSplineSurface>,
		    LRSplineSurface::PatchStatus> > surf_fragments = 
    lrsurf->subdivideIntoSimpler(threshold_missing, tol, &bddom);

  // Collect tensor-product spline patches
  vector<shared_ptr<ParamSurface> > tpsfs;
  if (surf_fragments.size() == 1)
    {
      // Change the complete surface into a tensor-product surface
      shared_ptr<SplineSurface> curr_tp(lrsurf->asSplineSurface());
      vector<shared_ptr<CurveOnSurface> > curr_loop;
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Fetch parameter curve
	  shared_ptr<ParamCurve> cv = loop[kr];
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
	  shared_ptr<ParamCurve> pcv;
	  if (cv->dimension() == 2)
	    pcv = cv;
	  else if (sfcv.get())
	    {
	      if (!sfcv->hasParameterCurve())
		sfcv->ensureParCrvExistence(epsge);
	      pcv = sfcv->parameterCurve();
	    }
	  else
	    continue;  // No parameter curve information

	  // Make CurveOnSurface curve with respect to current
	  // surface
	  shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(curr_tp,
							      pcv, true));
	  curr_loop.push_back(sfcv2);
	}
      shared_ptr<BoundedSurface> bdsurf2(new BoundedSurface(curr_tp, curr_loop,
							    epsge));
#ifdef DEBUG
      std::ofstream of0("tp_bd.g2");
      bdsurf2->writeStandardHeader(of0);
      bdsurf2->write(of0);
#endif

      tpsfs.push_back(bdsurf2);
    }
  else
    {
      for (size_t ki=0; ki<surf_fragments.size(); ++ki)
	{
	  LRSplineSurface::PatchStatus stat = surf_fragments[ki].second;
	  if (stat != LRSplineSurface::OUTSIDE)
	    {
#ifdef DEBUG
	      std::ofstream of3("curr_sub.g2");
	      surf_fragments[ki].first->writeStandardHeader(of3);
	      surf_fragments[ki].first->write(of3);
#endif
	      shared_ptr<SplineSurface> curr_tp(surf_fragments[ki].first->asSplineSurface());
	      if (stat == LRSplineSurface::INTERSECT)
		{
		  // Create trimmed surface
		  // First represent spline surface as trimmed
		  shared_ptr<BoundedSurface> curr_bd(new BoundedSurface(curr_tp,
									epsge));

#ifdef DEBUG
		  std::ofstream of4("par_crvs_trim.g2");
		  CurveLoop loop2 = SurfaceTools::outerBoundarySfLoop(curr_bd, epsge);
		  for (int ka=0; ka<loop2.size(); ++ka)
		    {
		      shared_ptr<ParamCurve> tmp = loop2[ka];
		      shared_ptr<CurveOnSurface> bdtmp = 
			dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
		      if (bdtmp.get())
			{
			  shared_ptr<ParamCurve> tmp2 = bdtmp->parameterCurve();
			  if (tmp2.get())
			    {
			      tmp2->writeStandardHeader(of4);
			      tmp2->write(of4);
			    }
			}
		    }
#endif
		  // Collect trimming loop pieces which is inside the
		  // current surface
		  vector<shared_ptr<CurveOnSurface> > cv_pieces;
		  for (int kr=0; kr<nmb; ++kr)
		    {
		      // Fetch parameter curve
		      shared_ptr<ParamCurve> cv = loop[kr];
		      shared_ptr<CurveOnSurface> sfcv =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		      shared_ptr<ParamCurve> pcv;
		      if (cv->dimension() == 2)
			pcv = cv;
		      else if (sfcv.get())
			{
			  if (!sfcv->hasParameterCurve())
			    sfcv->ensureParCrvExistence(epsge);
			  pcv = sfcv->parameterCurve();
			}
		      else
			continue;  // No parameter curve information
#ifdef DEBUG
		      of4 << "100 1 0 4 255 0 0 255" << std::endl;
		      pcv->write(of4);
#endif
		      // Make CurveOnSurface curve with respect to current
		      // surface
		      shared_ptr<CurveOnSurface> sfcv2(new CurveOnSurface(curr_tp,
									  pcv, true));

		      // Intersect curve with the trimming loop of the surface
		      vector<shared_ptr<CurveOnSurface> > trim_cvs = 
			BoundedUtils::intersectWithSurface(*sfcv2, *curr_bd,
							   epsge, true);
#ifdef DEBUG
		      std::ofstream of5("curve_piece.g2");
		      for (size_t kh=0; kh<trim_cvs.size(); ++kh)
			{
			  of5 << "100 1 0 4 0 255 0 255" << std::endl;
			  trim_cvs[kh]->parameterCurve()->write(of5);
			}
#endif
		      cv_pieces.insert(cv_pieces.end(), trim_cvs.begin(), 
				       trim_cvs.end());
		    }
	      
		  // Split surface with respect to trim segments
		  vector<shared_ptr<BoundedSurface> > trim_sfs;
		  if (cv_pieces.size() > 0)
		    trim_sfs =
		      BoundedUtils::splitWithTrimSegments(curr_bd, cv_pieces,
							  epsge);
	      
		  // Select inside surface
		  for (size_t kj=0; kj<trim_sfs.size(); ++kj)
		    {
		      double upar, vpar;
		      Point inner = trim_sfs[kj]->getInternalPoint(upar, vpar);
		      Vector2D parpt(upar, vpar);
		      if (bddom.isInDomain(parpt, epsge))
			tpsfs.push_back(trim_sfs[kj]);
		    }
		}
	      else
		tpsfs.push_back(curr_tp);
	    }
	}
    }

  // Compute extremal points
  Point vec(1);
  vec[0] = (sgn < 0) ? -1.0 : 1.0;
  ExtremalPoint::computeExtremalPoints(tpsfs, vec, epsge, extpoints);
  int stop_break = 1;
  return tpsfs.size();
}


namespace {

  void addContourLoop(shared_ptr<ParamCurve> crv, double isoval,
		      const CurveLoop& parloop, double epsge,
		      vector<pair<vector<shared_ptr<ParamCurve> >, double> >& crv_loops,
		      vector<BoundingBox>& bbox)
  {
    // Compute distance between curve endpoints
    Point pos1 = crv->point(crv->startparam());
    Point pos2 = crv->point(crv->endparam());
    
    // Assign bounding boxes
    BoundingBox bb = crv->boundingBox();
    vector<shared_ptr<ParamCurve> > contour_loop;
    contour_loop.push_back(crv);
    double dist0 = pos1.dist(pos2);
    if (dist0 > epsge)
      {
	// Close the curve loop
	int ind1, ind2;
	double par1, par2, dist1, dist2;
	Point pt1, pt2;
	parloop.closestPoint(pos1, ind1, par1, pt1, dist1);
	parloop.closestPoint(pos2, ind2, par2, pt2, dist2);

	if (dist0 < std::min(dist1, dist2))
	    {
	      // Close original loop
	      shared_ptr<ParamCurve> crv2(new SplineCurve(pos1, pos2));
	      BoundingBox bb2 = crv2->boundingBox();
	      contour_loop.push_back(crv2);
	      bb.addUnionWith(bb2);
	      crv_loops.push_back(make_pair(contour_loop,isoval));
	      bbox.push_back(bb);
	    }
	  else
	    {
	      // Include part of surface loop in the curve loop
	      // Compute curve lengths
	      double len0 = crv->estimatedCurveLength();
	      int nmb = parloop.size();
	      if (ind1 > ind2 || (ind1 == ind2 && par1 > par2))
		{
		  std::swap(ind1, ind2);
		  std::swap(par1, par2);
		}
	      double len1 = 
		parloop[ind1]->estimatedCurveLength(par1, 
						    (ind2 != ind1) ? 
						    parloop[ind1]->endparam() : par2);
	      for (int ka=ind1+1; ka<ind2; ++ka)
		len1 += parloop[ka]->estimatedCurveLength();
	      if (ind1 != ind2)
		len1 += parloop[ind2]->estimatedCurveLength(parloop[ind2]->startparam(),
							 par2);

	      double len2 = 
		parloop[ind2]->estimatedCurveLength(par2, parloop[ind2]->endparam());
	      for (int ka=(ind2+1)%nmb; ka!=ind1; ka=(ka+1)%nmb)
		len2 += parloop[ka]->estimatedCurveLength();
	      len2 += parloop[ind1]->estimatedCurveLength(parloop[ind1]->startparam(),
						       par1);
	      
	      vector<shared_ptr<ParamCurve> > contour_loop2;
	      contour_loop2.push_back(shared_ptr<ParamCurve>(crv->clone()));
	      BoundingBox bb3 = bb;
	      if (len2 >= len0 + len1)
		{
		  for (int ka=ind1; ka<=ind2; ++ka)
		    {
		      double t1 = (ka == ind1) ? par1 : parloop[ka]->startparam();
		      double t2 = (ka == ind2) ? par2 : parloop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> crv2;
		      if (ka != ind1 && ka != ind2) 
			crv2 = shared_ptr<ParamCurve>(parloop[ka]->clone());
		      else if (t2 - t1 >= epsge)
			crv2 =
			  shared_ptr<ParamCurve>(parloop[ka]->subCurve(t1, t2));
		      else
			continue;
		      BoundingBox bb2 = crv2->boundingBox();
		      contour_loop.push_back(crv2);
		      bb.addUnionWith(bb2);
		    }
		  crv_loops.push_back(make_pair(contour_loop,isoval));
		  bbox.push_back(bb);
		}

	      if (len1 >= len0 + len2)
		{
		  int ka, kb;
		  int nmb2 = nmb - (ind2-ind1);
		  for (ka=ind2, kb=0; kb<=nmb2; ka=(ka+1)%nmb, ++kb)
		    {
		      double t1 = (ka == ind2) ? par2 : parloop[ka]->startparam();
		      double t2 = (ka == ind1) ? par1 : parloop[ka]->endparam();
		      
		      shared_ptr<ParamCurve> crv2;
		      if (ka != ind1 && ka != ind2) 
			crv2 = shared_ptr<ParamCurve>(parloop[ka]->clone());
		      else if (t2 - t1 >= epsge)
			crv2 =
			  shared_ptr<ParamCurve>(parloop[ka]->subCurve(t1, t2));
		      else
			continue;
		      BoundingBox bb2 = crv2->boundingBox();
		      contour_loop2.push_back(crv2);
		      bb3.addUnionWith(bb2);
		    }
		  crv_loops.push_back(make_pair(contour_loop2,isoval));
		  bbox.push_back(bb3);
		}
	    }
	}
      else
	{
	  crv_loops.push_back(make_pair(contour_loop,isoval));
	  bbox.push_back(bb);
	}  
}

void 
extractInnerCurves(vector<pair<vector<shared_ptr<ParamCurve> >, double> >& cvs,
		   vector<BoundingBox>& bbox,
		   double eps)
{
  size_t ki, kj;
  for (ki=0; ki<cvs.size(); )
    {
      for (kj=ki+1; kj<cvs.size(); )
	{
	  if (!(bbox[ki].overlaps(bbox[kj])))
	    {
	      ++kj;
	      continue;  // No possibility of one curve lying inside the other
	    }

	  if (fabs(cvs[ki].second - cvs[kj].second) < eps)
	    {
	      ++kj;
	      continue;  // Same height of contour curves. Should mean both
	      // minimum and maximum point
	    }

	  double size1 = bbox[ki].low().dist(bbox[ki].high());
	  double size2 = bbox[kj].low().dist(bbox[kj].high());
	  size_t ix1 = (size1 < size2) ? ki : kj;
	  size_t ix2 = (ix1 == kj) ? ki : kj;

	  int ka = 0;
	  for (ka=0; ka<2; ++ka)
	    {
	      // Check if curve ix1 lies inside curve ix2
	      // Compute one point on the curve
	      Point pos = 
		cvs[ix1].first[0]->ParamCurve::point(0.5*(cvs[ix1].first[0]->startparam()+
							  cvs[ix1].first[0]->endparam()));
	      shared_ptr<CurveLoop> cvloop(new CurveLoop(cvs[ix2].first, eps, false));
	      CurveBoundedDomain cvdom(cvloop);

	      Vector2D pos2(pos[0], pos[1]);
	      bool inside = cvdom.isInDomain(pos2, eps);
	      if (inside)
		{
		  break;
		}
	      std::swap(ix1, ix2);
	    }

	  if (ka < 2)
	    {
	      // An "inside" curve is found. Remove the other one
	      cvs.erase(cvs.begin()+ix2);
	      bbox.erase(bbox.begin()+ix2);
	      if (ix2 == ki)
		break;
	    }
	  else
	    ++kj;
	}
      if (kj == cvs.size())
	++ki;
    }
}

}; // end anonymous namespace
