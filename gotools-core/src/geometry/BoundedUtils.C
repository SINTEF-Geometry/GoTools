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

//#define DEBUG1

//#define SBR_DBG

#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>
#include <utility>
#include "sisl.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Values.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/SurfaceOfLinearExtrusion.h"
#include "GoTools/creators/CoonsPatchGen.h"


using namespace Go;
using std::vector;
using std::pair;
using std::swap;
using CoonsPatchGen::blendcoef;

namespace {
// Given input of curves which may start in end pt of curr_cv (within space_eps), return
// index of the one with the greatest angle (in the interval (-pi, pi]).
// If there is no such curve, return -1.
int leftMostCurve(CurveOnSurface& cv,
		  vector<shared_ptr<CurveOnSurface> >& other_cvs,
		  double space_eps, bool par_cv, double& angle);

  // Used by leftMostCurve in tangential situations
void recomputeAngles(vector<Point>& curr_val, double ang_tol,
		     vector<shared_ptr<CurveOnSurface> >& cvs,
		     int idx1, double& ang1, int idx2, double& ang2);

  double getSeed(Point space_pt, CurveOnSurface& cv_on_sf, bool par_cv);

  // Given a point at the seam, we use a marching approach to decide which side to choose.
  // If the curve follows the seam at all values this test is inconclusive.
  void marchOutSeamPoint(const ParamSurface& surface, const ParamCurve& space_cv,
			 double tpar, bool to_the_right, bool at_u_seam, bool at_v_seam,
			 double epsgeo, Point& par_pt, bool& success);

  double maxCurveGap(shared_ptr<CurveOnSurface> curves[], int nmb_curves,
		     bool par_crv, double eps, vector<Point>& free_endpt);

  void splitLoopCvs(const BoundedSurface& sf,
		    vector<shared_ptr<CurveOnSurface> >& old_loop_cvs,
		    vector<shared_ptr<CurveOnSurface> >& part_bd_cvs,
		    vector<Point>& part_bd_endpt, double min_loop_tol,
		    double eps, double epspar, double knot_diff_tol,
		    int last_split, bool par_cv);
}; // end anonymous namespace 

namespace Go {

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::intersectWithSurface(CurveOnSurface& curve,
				   BoundedSurface& bounded_surf,
				   double epsge, bool only_inner_curves)
//===========================================================================
{
  double int_tol = 0.1*epsge; //1e-06;
  shared_ptr<const ParamSurface> under_sf = bounded_surf.underlyingSurface();
  double ptol = getParEps(epsge, under_sf.get());
  int_tol = ptol;
  double tdel = curve.endparam() - curve.startparam();
  double delfac = 0.01;
  //double fuzzy = 1.0e-8;

    // We extract boundary loop, and check for intersections.
    // We do not handle trimming of trimmed surfaces with holes.
    //vector<CurveLoop> boundary_loops = bounded_surf.allBoundaryLoops(epsge);
  vector<CurveLoop> boundary_loops = bounded_surf.allBoundaryLoops(DEFAULT_SPACE_EPSILON);

    // Knowing that bounded_surf is BoundedSurface, we extract CurveOnSurface's.
    vector<shared_ptr<CurveOnSurface> > loop_curves;
    int i, j, kk;
    for (i = 0; i < int(boundary_loops.size()); ++i) {
	for (j = 0; j < boundary_loops[i].size(); ++j) {
	    loop_curves.push_back(dynamic_pointer_cast<CurveOnSurface, ParamCurve>(boundary_loops[i][j]));
	}
    }

    // We make sure that all calculations rely on the parameter curves.
    curve.ensureParCrvExistence(epsge);
    if (!curve.parameterCurve())
      THROW("Not possible to generate parameter curve");
    // if (!curve.parPref()) {
    // 	curve = CurveOnSurface(curve.underlyingSurface(), curve.parameterCurve(),
    // 			       curve.spaceCurve(), true);
    // }
    for (i = 0; i < int(loop_curves.size()); ++i) {
      loop_curves[i]->ensureParCrvExistence(epsge);
      if (!loop_curves[i]->parameterCurve())
	THROW("Not possible to generate parameter curve");
// 	if (!loop_curves[i]->parPref()) {
// 	    loop_curves[i] = shared_ptr<CurveOnSurface>
// 		(new CurveOnSurface(loop_curves[i]->underlyingSurface(),
// 				    loop_curves[i]->parameterCurve(),
// 				    loop_curves[i]->spaceCurve(), 
// 				    loop_curves[i]->parPref()));
// 	}
    }

    const CurveBoundedDomain& domain = bounded_surf.parameterDomain();

    vector<pair<double, double> > int_params;
    vector<pair<pair<double, double>,pair<double, double> > > int_crvs;
    vector<double> all_int_params;
    vector<int> pretop;
    shared_ptr<SplineCurve> first_curve =
	dynamic_pointer_cast<SplineCurve, ParamCurve>(curve.parameterCurve());
    if (first_curve.get() == 0)
      first_curve = 
	shared_ptr<SplineCurve>(curve.parameterCurve()->geometryCurve());
    ALWAYS_ERROR_IF(first_curve.get() == 0,
		    "Intersection routine not implemented for general curves.");

    for (j = 0; j < int(loop_curves.size()); ++j) {
	shared_ptr<SplineCurve> second_curve =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(loop_curves[j]->parameterCurve());
	if (second_curve.get() == 0)
	  second_curve = 
	    shared_ptr<SplineCurve>(loop_curves[j]->parameterCurve()->geometryCurve());
	ALWAYS_ERROR_IF(second_curve.get() == 0,
			"Intersect. routine not implemented for general curves.");

	int_params.clear();
	// We use SISL for finding curve-curve intersection.
	intersect2Dcurves(first_curve.get(), second_curve.get(), int_tol,
			  int_params, pretop, int_crvs);
	for (kk = 0; kk < int(int_params.size()); ++kk) 
	  {
	    // Evaluate
	    Point pt1, pt2;
	    curve.point(pt1, int_params[kk].first);
	    loop_curves[j]->point(pt2, int_params[kk].second);
	    double init_dist = pt1.dist(pt2);

	    // Iterate to take the curve preferences into account
	    double par1, par2, dist;
	    Point ptc1, ptc2;
	    int max_passes = 10;
	    ClosestPoint::closestPtCurves(&curve, loop_curves[j].get(), curve.startparam(),
					  curve.endparam(), loop_curves[j]->startparam(),
					  loop_curves[j]->endparam(),
					  int_params[kk].first, int_params[kk].second,
					  par1, par2, dist, ptc1, ptc2, max_passes);
	  
	    if (dist < init_dist)
	      all_int_params.push_back(par1);
	    else
	      all_int_params.push_back(int_params[kk].first);
	}
    }

    // We sort the intersection parameters in ascending order.
    double fuzzy = 1.0e-10;
    fuzzy *= (first_curve->endparam() - first_curve->startparam());
    fuzzy = std::max(1.0e-10, std::min(fuzzy, 1.0e-6));
    sort(all_int_params.begin(), all_int_params.end());
    if (all_int_params.size() > 0 && 
	all_int_params[0]-first_curve->startparam() < fuzzy)
      all_int_params[0] = first_curve->startparam();
    else if ((all_int_params.size() == 0) || 
	     (all_int_params.front() != first_curve->startparam())) {
    	all_int_params.insert(all_int_params.begin(), first_curve->startparam());
    }
    if (all_int_params.size() > 0 && 
	first_curve->endparam()-all_int_params[all_int_params.size()-1] < fuzzy)
      all_int_params[all_int_params.size()-1] = first_curve->endparam();
    else if (all_int_params.back() != first_curve->endparam()) {
    	all_int_params.insert(all_int_params.end(), first_curve->endparam());
    }

    // We extract parts of curve inside trimmed domain.
    // First mark segments
    vector<int> int_seg_type(all_int_params.size()-1, -1);
    // 0 = short curve, 1 = outside, 2 = inside, 3 = on the boundary
    vector<shared_ptr<CurveOnSurface> > inside_segments;
    double knot_diff_tol = 0.01*getParEps(epsge, under_sf.get()); // We may not trust pcv to repr space_cv.
    for (j = 0; j < int(all_int_params.size()) - 1; ++j) {
	double from_par = all_int_params[j];
	double to_par = all_int_params[j+1];
	Point tmp1 = first_curve->ParamCurve::point(from_par);
	Point tmp2 = first_curve->ParamCurve::point(0.5*(from_par+to_par));
	Point tmp3 = first_curve->ParamCurve::point(to_par);
	double len = tmp1.dist(tmp2) + tmp2.dist(tmp3);
	if (to_par - from_par < knot_diff_tol || len < epsge) {
	    // all_int_params.erase(all_int_params.begin() + j + 1);
	    // --j;
	    // continue;
	  int_seg_type[j] = 0;  // Short curve
	  continue;
	}
	if (from_par < curve.startparam())
	  from_par = curve.startparam();
	if (to_par > curve.endparam())
	  to_par = curve.endparam();
	double med_par = 0.5*(from_par + to_par);
	Point med_pt = first_curve->ParamCurve::point(med_par);
	int is_in_domain = -1;
	if (to_par-from_par < delfac*tdel)
	  {
	    // Extra testing
	    Point med3 = curve.ParamCurve::point(med_par);
	    Point med_sf = curve.underlyingSurface()->point(med_pt[0],med_pt[1]);
	    double med_dist = med3.dist(med_sf);
	    if (true) //med_dist > epsge)
	      {
		// Perform test in geometry space
		double u1, u2, v1, v2, d1, d2;
		Point clo1, clo2;
		bounded_surf.closestBoundaryPoint(med3, u1, v1, clo1, 
						  d1, epsge);
		if (d1 < epsge)
		  is_in_domain = 2;
		else
		  {
		    bounded_surf.closestPoint(med3, u2, v2, clo2, 
					      d2, epsge);
		    if (d2 < epsge)
		      is_in_domain = 1;
		    else
		      is_in_domain = 0;
		  }
	      }
	  }
	if (is_in_domain < 0)
	  {
	    try {
	      is_in_domain = domain.isInDomain2(Vector2D(med_pt[0], med_pt[1]), 
						ptol);
	      if (is_in_domain == 2)
		{
		  int is_in_domain2 = 
		    domain.isInDomain2(Vector2D(med_pt[0], med_pt[1]), 
				       knot_diff_tol);
		  if (is_in_domain2 == 1)
		    is_in_domain = 1;
		}
	    }
	    catch (...)
	      {
		if (len > 10.0*epsge)
		  {
		    // Use geometric test
		    Point med3 = curve.ParamCurve::point(med_par);
		    Point med_sf = 
		      curve.underlyingSurface()->point(med_pt[0],med_pt[1]);
		    double med_dist = med3.dist(med_sf);

		    double u1, u2, v1, v2, d1, d2;
		    Point clo1, clo2;
		    bounded_surf.closestBoundaryPoint(med3, u1, v1, clo1, 
						      d1, epsge);
		    if (d1 < epsge)
		      is_in_domain = 2;
		    else
		      {
			bounded_surf.closestPoint(med3, u2, v2, clo2, 
						  d2, epsge);
			if (d2 < epsge)
			  is_in_domain = 1;
			else
			  is_in_domain = 0;
		      }
		  }
		else
		  is_in_domain = 0;
	      }
	  }
	int_seg_type[j] = /*(is_in_domain >= 1) ? 2 : 1; /*/is_in_domain + 1;
	// if (is_in_domain == 1)
	//   {
	//     // inside_segments.push_back(shared_ptr<CurveOnSurface>
	//     // 			    (dynamic_cast<CurveOnSurface*>
	//     // 			     (curve.subCurve(from_par, to_par))));
	//     int_seg_type[j] = 2;  // Inside segment
	//   }
	// else
	//   int_seg_type[j] = 1;    // Outside segment
    }

    // Simplify segmentation by joining segments of the same type (inside/outside)
    for (j=0; j<int(int_seg_type.size()); ++j)
      {
	int k;
	for (k=j+1; k<int(int_seg_type.size()); ++k)
	  if (int_seg_type[k] != int_seg_type[j])
	    break;
	if (k > j+1)
	  {
	    // Simplify
	    int type = int_seg_type[j];
	    if (int_seg_type[j] == 0)
	      {
		// Check if the segment is still small
		double from_par = all_int_params[j];
		double to_par = all_int_params[k];
		Point tmp1 = first_curve->ParamCurve::point(from_par);
		Point tmp2 = first_curve->ParamCurve::point(0.5*(from_par+to_par));
		Point tmp3 = first_curve->ParamCurve::point(to_par);
		double len = tmp1.dist(tmp2) + tmp2.dist(tmp3);
		if (to_par - from_par > knot_diff_tol && len > epsge) 
		  {
		    Point med_pt = first_curve->ParamCurve::point(0.5*(from_par+to_par));
		    int is_in_domain = 0;
		    try {
		      is_in_domain = 
			domain.isInDomain2(Vector2D(med_pt[0], med_pt[1]), 
					   knot_diff_tol);
		    }
		    catch (...)
		      {
			if (len > 10.0*epsge)
			  THROW("Could not trim intersection curve with surface boundary");
			else
			  is_in_domain = 0;
		      }
		    type = (is_in_domain >= 1) ? 2 : 1;
		  }
		else
		  {
		    all_int_params.erase(all_int_params.begin()+j+1, all_int_params.begin()+k);
		    int_seg_type.erase(int_seg_type.begin()+j+1, int_seg_type.begin()+k);
		  }
	      }
	    else
	      {
		all_int_params.erase(all_int_params.begin()+j+1, all_int_params.begin()+k);
		int_seg_type.erase(int_seg_type.begin()+j+1, int_seg_type.begin()+k);
	      }
	    int_seg_type[j] = type;
	  }
      }

    for (j=1; j<int(int_seg_type.size()-1); ++j)
      {
	if (int_seg_type[j] == 0)
	  {
	    if (int_seg_type[j-1] == int_seg_type[j+1])
	      {
		// Remove small segment
		all_int_params.erase(all_int_params.begin()+j, all_int_params.begin()+j+1);
		int_seg_type.erase(int_seg_type.begin()+j, int_seg_type.begin()+j+1);
		--j;
	      }
	  }
      }

    // Extract segments
    for (j=1; j<int(all_int_params.size()); ++j)
      {
	if (int_seg_type[j-1] == 2 || 
	    (only_inner_curves == false && int_seg_type[j-1] == 3))
	  {
	    inside_segments.push_back(shared_ptr<CurveOnSurface>
				      (dynamic_cast<CurveOnSurface*>
				       (curve.subCurve(all_int_params[j-1], 
						       all_int_params[j],
						       fuzzy))));
	  }
      }

    return inside_segments;
}


//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::trimWithPlane(const shared_ptr<ParamSurface>& surf,
			      Point point, Point normal, double epsge)
//===========================================================================
{
    int ki;
    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    shared_ptr<BoundedSurface> bounded_sf;
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	    vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));
	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();

    vector<shared_ptr<CurveOnSurface> > int_segments;
    try {
      int_segments =
	  intersectWithPlane(underlying_sf, point, normal, epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with plane.");
    }
    vector<shared_ptr<BoundedSurface> > return_sfs;
    if (int_segments.size() == 0)
	return return_sfs; // Empty intersection, return empty vector.

    // We extract those parts of int_segments which are actually on bounded surface.
    vector<shared_ptr<CurveOnSurface> > interior_segments;
    for (ki = 0; ki < int(int_segments.size()); ++ki) {
	vector<shared_ptr<CurveOnSurface> > new_segments =
	    intersectWithSurface(*int_segments[ki],
				 *bounded_sf, 0.1*epsge);
	interior_segments.insert(interior_segments.begin(),
				 new_segments.begin(), new_segments.end());
    }

    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
    try {
      loop_curves = getBoundaryLoops(*bounded_sf, interior_segments, epsge);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }

    return_sfs = createTrimmedSurfs(loop_curves, bounded_sf->underlyingSurface(), epsge);

    return return_sfs;
}

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::getPlaneIntersections(const shared_ptr<ParamSurface>& surf,
				    Point point, Point normal, double epsge,
				    shared_ptr<BoundedSurface>& bounded_sf)
//===========================================================================
{
    int ki;
    vector<shared_ptr<CurveOnSurface> > interior_segments;

    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();

    vector<shared_ptr<CurveOnSurface> > int_segments;
    try {
      int_segments =
	  intersectWithPlane(underlying_sf, point, normal, 0.7*epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with plane.");
    }

#ifdef DEBUG1
    std::ofstream debug("int_crv.g2");
    bounded_sf->writeStandardHeader(debug);
    bounded_sf->write(debug);
    for (ki=0; ki<int(int_segments.size()); ++ki)
      {
	int_segments[ki]->spaceCurve()->writeStandardHeader(debug);
	int_segments[ki]->spaceCurve()->write(debug);
      }
#endif
    
    // We extract those parts of int_segments which are actually on bounded surface.
    for (ki = 0; ki < int(int_segments.size()); ++ki) {
	vector<shared_ptr<CurveOnSurface> > new_segments =
	    intersectWithSurface(*int_segments[ki],
				 *bounded_sf, 0.1*epsge);
	interior_segments.insert(interior_segments.begin(),
				 new_segments.begin(), new_segments.end());
    }
    return interior_segments;
}

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::getCylinderIntersections(const shared_ptr<ParamSurface>& surf,
				       Point point, Point axis, double radius, 
				       double epsge,
				       shared_ptr<BoundedSurface>& bounded_sf)
//===========================================================================
{
    int ki;
    vector<shared_ptr<CurveOnSurface> > interior_segments;

    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();

    vector<shared_ptr<CurveOnSurface> > int_segments;
    try {
      int_segments =
	intersectWithCylinder(underlying_sf, point, axis, radius, 0.7*epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with cylinder.");
    }

#ifdef DEBUG1
    std::ofstream debug("int_crv.g2");
    bounded_sf->writeStandardHeader(debug);
    bounded_sf->write(debug);
    for (ki=0; ki<int(int_segments.size()); ++ki)
      {
	int_segments[ki]->spaceCurve()->writeStandardHeader(debug);
	int_segments[ki]->spaceCurve()->write(debug);
      }
#endif
    
    // We extract those parts of int_segments which are actually on bounded surface.
    for (ki = 0; ki < int(int_segments.size()); ++ki) {
	vector<shared_ptr<CurveOnSurface> > new_segments =
	    intersectWithSurface(*int_segments[ki],
				 *bounded_sf, 0.1*epsge);
	interior_segments.insert(interior_segments.begin(),
				 new_segments.begin(), new_segments.end());
    }
    return interior_segments;
}

//===========================================================================
void 
BoundedUtils::getSurfaceIntersections(const shared_ptr<ParamSurface>& surf1,
				      const shared_ptr<ParamSurface>& surf2,
				      double epsge,
				      vector<shared_ptr<CurveOnSurface> >& int_cv1,
				      shared_ptr<BoundedSurface>& bounded_sf1,
				      vector<shared_ptr<CurveOnSurface> >& int_cv2,
				      shared_ptr<BoundedSurface>& bounded_sf2,
				      bool only_inner_curves)
//===========================================================================
{
    vector<shared_ptr<CurveOnSurface> > int_segments1;
    vector<shared_ptr<CurveOnSurface> > int_segments2;

    // We convert the input surfaces to bounded surfaces (given that input is
    // a bounded or spline surface).
    if (surf1->instanceType() == Class_BoundedSurface)
	bounded_sf1 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf1);
    else //if (surf1->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf1,
								 epsge);
	    bounded_sf1 = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf1, loops));
#ifdef DEBUG1
	    int state;
	    bounded_sf1->analyzeLoops();
// 	    bool valid = bounded_sf1->isValid(state);
// 	    if (!valid)
// 	      std::cout << "Surface not valid: " << state << std::endl;
#endif
	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    if (surf2->instanceType() == Class_BoundedSurface)
	bounded_sf2 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf2);
    else //if (surf1->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf2,
								 epsge);
	    bounded_sf2 = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf2, loops));
#ifdef DEBUG1
	    int state;
	    bounded_sf2->analyzeLoops();
// 	    bool valid = bounded_sf2->isValid(state);
// 	    if (!valid)
// 	      std::cout << "Surface not valid: " << state << std::endl;
#endif
	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

//     shared_ptr<SplineSurface> under_sf1 =
// 	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf1->underlyingSurface());
//     shared_ptr<SplineSurface> under_sf2 =
// 	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf2->underlyingSurface());
//     ALWAYS_ERROR_IF(under_sf1.get() == 0 || under_sf2.get() == 0,
// 		    "Requiring that underlying surface is a SplineSurface.");
    shared_ptr<ParamSurface> under_sf1 = bounded_sf1->underlyingSurface();
    shared_ptr<ParamSurface> under_sf2 = bounded_sf2->underlyingSurface();

    // Avoid major extension of intersection domain
    double domainfac = 1.5;
    double redfac = 0.05;
    RectDomain dom1 = bounded_sf1->containingDomain();
    RectDomain dom2 = under_sf1->containingDomain();
    double ll1 = dom1.diagLength();
    double ll2 = dom2.diagLength();
    if (ll2 > domainfac*ll1)
      {
	Point par_eps = SurfaceTools::getParEpsilon(*under_sf2, epsge);
	double u1 = std::max(dom1.umin()-5.0*std::max(par_eps[0],epsge), 
			     dom2.umin());
	double u2 = std::min(dom1.umax()+5.0*std::max(par_eps[0],epsge), 
			     dom2.umax());
	double v1 = std::max(dom1.vmin()-5.0*std::max(par_eps[1],epsge), 
			     dom2.vmin());
	double v2 = std::min(dom1.vmax()+5.0*std::max(par_eps[1],epsge), 
			     dom2.vmax());
	vector<shared_ptr<ParamSurface> > sub =
	  under_sf1->subSurfaces(u1, v1, u2, v2);
	if (sub.size() == 1)
	  under_sf1 = sub[0];
      }
    RectDomain dom3 = bounded_sf2->containingDomain();
    RectDomain dom4 = under_sf2->containingDomain();
    double ll3 = dom3.diagLength();
    double ll4 = dom4.diagLength();
    if (ll4 > domainfac*ll3)
      {
	Point par_eps = SurfaceTools::getParEpsilon(*under_sf2, epsge);
	double u1 = std::max(dom3.umin()-5.0*std::max(par_eps[0],epsge), 
			     dom4.umin());
	double u2 = std::min(dom3.umax()+5.0*std::max(par_eps[0],epsge), 
			     dom4.umax());
	double v1 = std::max(dom3.vmin()-5.0*std::max(par_eps[1],epsge), 
			     dom4.vmin());
	double v2 = std::min(dom3.vmax()+5.0*std::max(par_eps[1],epsge), 
			     dom4.vmax());
	vector<shared_ptr<ParamSurface> > sub =
	  under_sf2->subSurfaces(u1, v1, u2, v2);
	if (sub.size() == 1)
	  under_sf2 = sub[0];
      }

#ifdef DEBUG1
    std::ofstream debug0("int_crv_surf0.g2");
    under_sf1->writeStandardHeader(debug0);
    under_sf1->write(debug0);
    under_sf2->writeStandardHeader(debug0);
    under_sf2->write(debug0);
#endif
    try {
      getIntersectionCurve(under_sf1, under_sf2, int_segments1, int_segments2, 
			   epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with plane.");
    }

    // Ensure correct surface pointer
    for (size_t kj=0; kj<int_segments1.size(); ++kj)
      {
	int_segments1[kj]->setUnderlyingSurface(bounded_sf1->underlyingSurface());
	int_segments2[kj]->setUnderlyingSurface(bounded_sf2->underlyingSurface());
      }

#ifdef DEBUG1
    std::ofstream debug("int_crv_surf.g2");
    bounded_sf1->writeStandardHeader(debug);
    bounded_sf1->write(debug);
    bounded_sf2->writeStandardHeader(debug);
    bounded_sf2->write(debug);
    for (int ki=0; ki<int(int_segments1.size()); ++ki)
      {
	int_segments1[ki]->spaceCurve()->writeStandardHeader(debug);
	int_segments1[ki]->spaceCurve()->write(debug);
      }
#endif
    
    // We extract those parts of int_segments which are actually on bounded surfaces.
    if (int_segments1.size() > 0)
      {
	try {
	  intersectWithSurfaces(int_segments1, bounded_sf1, int_segments2, 
				bounded_sf2, /*0.1**/epsge, only_inner_curves);
	} catch (...) {
	  //THROW("Failed intersecting the two spline surfaces.");
	}
      }

#ifdef DEBUG1
    std::ofstream debug2("int_crv_surf2.g2");
    for (int ki=0; ki<int(int_segments1.size()); ++ki)
      {
	int_segments1[ki]->spaceCurve()->writeStandardHeader(debug2);
	int_segments1[ki]->spaceCurve()->write(debug2);

      }
#endif

     if (int_segments1.size() > 0)
      {
	int_cv1.insert(int_cv1.end(), int_segments1.begin(), 
		       int_segments1.end());
	int_cv2.insert(int_cv2.end(), int_segments2.begin(), 
		       int_segments2.end());
      }
}

//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::splitWithPlane(const shared_ptr<ParamSurface>& surf,
			     Point point, Point normal, double epsge)
//===========================================================================
{
    int ki;
    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    shared_ptr<BoundedSurface> bounded_sf;
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();

    vector<shared_ptr<CurveOnSurface> > int_segments;
    try {
      int_segments =
	  intersectWithPlane(underlying_sf, point, normal, epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with plane.");
    }
    vector<shared_ptr<BoundedSurface> > return_sfs;
    if (int_segments.size() == 0)
	return return_sfs; // Empty intersection, return empty vector.

    // We extract those parts of int_segments which are actually on bounded surface.
    vector<shared_ptr<CurveOnSurface> > interior_segments;
    for (ki = 0; ki < int(int_segments.size()); ++ki) {
	vector<shared_ptr<CurveOnSurface> > new_segments =
	    intersectWithSurface(*int_segments[ki],
				 *bounded_sf, 0.1*epsge);
	interior_segments.insert(interior_segments.begin(),
				 new_segments.begin(), new_segments.end());
    }

    // Insert the intersection segment also in the opposite direction
    size_t nmb_seg = interior_segments.size();
    for (size_t kr=0; kr<nmb_seg; ++kr)
      {
	shared_ptr<CurveOnSurface> tmp_seg = 
	  shared_ptr<CurveOnSurface>(interior_segments[kr]->clone());
	tmp_seg->reverseParameterDirection();
	interior_segments.push_back(tmp_seg);
      }
	    
    

    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
    try {
      loop_curves = getBoundaryLoops(*bounded_sf, interior_segments, 
				     epsge, (int)nmb_seg);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }

    return_sfs = createTrimmedSurfs(loop_curves, bounded_sf->underlyingSurface(), epsge);

    return return_sfs;
}


//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::splitBetweenParams(const shared_ptr<ParamSurface>& surf,
				 Point parval1, Point parval2, double epsge)
//===========================================================================
{
    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    shared_ptr<BoundedSurface> bounded_sf;
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<SplineSurface> underlying_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf->underlyingSurface());
    ALWAYS_ERROR_IF(underlying_sf.get() == 0,
		    "Requiring that underlying surface is a SplineSurface.");

    // Make trimming curve
    double ptol = 1.0e-8;
    shared_ptr<CurveOnSurface> trimcrv;
    if (fabs(parval1[0]-parval2[0]) < ptol)
      {
	trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 1,
							parval1[0], parval1[1],
							parval2[1], -1));
      }
    else if (fabs(parval1[1]-parval2[1]) < ptol)
      {
	trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 2,
							parval1[1], parval1[0],
							parval2[0], -1));
      }
    else
      {
	// Make parameter curve between parameter values
	shared_ptr<ParamCurve> pcrv = 
	  shared_ptr<ParamCurve>(new SplineCurve(parval1, parval2));

	// Make trimming curve
	trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, pcrv, true));
	bool updated;
	updated = trimcrv->ensureSpaceCrvExistence(epsge);
      }
    
    // Define the new trimming loop segment in both directions
    vector<shared_ptr<CurveOnSurface> > segments;
    segments.push_back(trimcrv);
    shared_ptr<CurveOnSurface> tmp_seg = 
      shared_ptr<CurveOnSurface>(trimcrv->clone());
    tmp_seg->reverseParameterDirection();
    segments.push_back(tmp_seg);

    // Define boundary loops of the new surfaces
    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
    try {
      loop_curves = getBoundaryLoops(*bounded_sf, segments, epsge, 1);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }

    // Make surfaces
    vector<shared_ptr<BoundedSurface> > return_sfs = 
      createTrimmedSurfs(loop_curves, bounded_sf->underlyingSurface(), epsge);

    return return_sfs;
}

//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::splitBetweenParPairs(const shared_ptr<ParamSurface>& surf,
				   vector<pair<Point,Point> > parvals, 
				   double epsge)
//===========================================================================
{
    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    shared_ptr<BoundedSurface> bounded_sf;
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<SplineSurface> underlying_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf->underlyingSurface());
    ALWAYS_ERROR_IF(underlying_sf.get() == 0,
		    "Requiring that underlying surface is a SplineSurface.");

    vector<shared_ptr<CurveOnSurface> > segments;
    for (size_t ki=0; ki<parvals.size(); ++ki)
      {
	// Make parameter curve between parameter values
	shared_ptr<ParamCurve> pcrv = 
	  shared_ptr<ParamCurve>(new SplineCurve(parvals[ki].first, 
						 parvals[ki].second));

	// Make trimming curve
	shared_ptr<CurveOnSurface> trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, pcrv, true));
	bool updated;
	updated = trimcrv->ensureSpaceCrvExistence(epsge);
    
	// Define the new trimming loop segment in both directions
	segments.push_back(trimcrv);
      }

    int nmb_seg = (int)segments.size();
    for (int kj=0; kj<nmb_seg; ++kj)
      {
	shared_ptr<CurveOnSurface> tmp_seg = 
	  shared_ptr<CurveOnSurface>(segments[kj]->clone());
	tmp_seg->reverseParameterDirection();
	segments.push_back(tmp_seg);
      }

    // Define boundary loops of the new surfaces
    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
    try {
      loop_curves = getBoundaryLoops(*bounded_sf, segments, epsge, nmb_seg);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }

    // Make surfaces
    vector<shared_ptr<BoundedSurface> > return_sfs = 
      createTrimmedSurfs(loop_curves, bounded_sf->underlyingSurface(), epsge);

    return return_sfs;
}


//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::getTrimCrvsParam(const shared_ptr<ParamSurface>& surf,
			       Point parval1, Point parval2, double epsge,
			       shared_ptr<BoundedSurface>& bounded_sf)
//===========================================================================
{
    // We convert the input surface to a bounded surface (given that input is
    // a bounded or spline surface).
    if (surf->instanceType() == Class_BoundedSurface)
	bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    else //if (surf->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf,
								 epsge);
	    bounded_sf = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();

    // Make trimming curve
    double ptol = 1.0e-8;
    shared_ptr<CurveOnSurface> trimcrv;
    if (fabs(parval1[0]-parval2[0]) < ptol)
      {
	trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 1,
							parval1[0], 
							std::min(parval1[1],parval2[1]),
							std::max(parval1[1],parval2[1]), 
							-1));
	if (parval2[1] < parval1[1])
	  trimcrv->reverseParameterDirection();
      }
    else if (fabs(parval1[1]-parval2[1]) < ptol)
      {
	trimcrv = 
	  shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 2,
							parval1[1], 
							std::min(parval1[0],parval2[0]),
							std::max(parval1[0],parval2[0]), 
							-1));
	if (parval2[0] < parval1[0])
	  trimcrv->reverseParameterDirection();
      }
    else
      {    
	// Make parameter curve between parameter values
	shared_ptr<ParamCurve> pcrv = 
	  shared_ptr<ParamCurve>(new SplineCurve(parval1, parval2));

	trimcrv = shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 
								pcrv, true));
      }

    bool updated;
    updated = trimcrv->ensureSpaceCrvExistence(epsge);
    
    // Define the new trimming loop segment in both directions
    vector<shared_ptr<CurveOnSurface> > segments =
      intersectWithSurface(*trimcrv, *bounded_sf, 0.1*epsge);

    return segments;
}

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::getTrimCrvsPcrv(const shared_ptr<ParamSurface>& surf,
			      shared_ptr<ParamCurve>& pcurve, double epsge,
			      shared_ptr<BoundedSurface>& bounded_sf)
//===========================================================================
{
  if (surf->instanceType() == Class_BoundedSurface)
    bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  else //if (surf->instanceType() == Class_SplineSurface)
    try {
      vector<CurveLoop> loops = 
	SurfaceTools::absolutelyAllBoundarySfLoops(surf,epsge);
      bounded_sf = 
	shared_ptr<BoundedSurface>(new BoundedSurface(surf, loops));
		  
    } catch (...) {
      THROW("Something went wrong, returning.");
    }
  shared_ptr<ParamSurface> underlying_sf = bounded_sf->underlyingSurface();
  shared_ptr<CurveOnSurface> trimcrv = 
    shared_ptr<CurveOnSurface>(new CurveOnSurface(underlying_sf, 
						  pcurve, true));
  bool updated;
  updated = trimcrv->ensureSpaceCrvExistence(epsge);
    
  // Define the new trimming loop segment in both directions
  vector<shared_ptr<CurveOnSurface> > segments =
    intersectWithSurface(*trimcrv, *bounded_sf, 0.1*epsge);

  return segments;
}
//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::trimSurfWithSurf(const shared_ptr<ParamSurface>& sf1,
			       const shared_ptr<ParamSurface>& sf2,
			       double epsge)
//===========================================================================
{
    // Transform the two surfaces to bounded surfaces.
    shared_ptr<BoundedSurface> bounded_sf1, bounded_sf2;
    if (sf1->instanceType() == Class_BoundedSurface)
	bounded_sf1 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf1);
    else //if (sf1->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(sf1,
								 epsge);
	    bounded_sf1 = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(sf1, loops));

	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");
    if (sf2->instanceType() == Class_BoundedSurface)
	bounded_sf2 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf2);
    else //if (sf2->instanceType() == Class_SplineSurface)
	try {
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(sf2,
								 epsge);
	    bounded_sf2 = 
	      shared_ptr<BoundedSurface>(new BoundedSurface(sf2, loops));


	} catch (...) {
	    THROW("Something went wrong, returning.");
	}
//     else
// 	THROW("Surface type not supported.");

    shared_ptr<SplineSurface> under_sf1 =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf1->underlyingSurface());
    shared_ptr<SplineSurface> under_sf2 =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf2->underlyingSurface());
    ALWAYS_ERROR_IF(under_sf1.get() == 0 || under_sf2.get() == 0,
		"Underlying surfaces must both be of type SplineSurface!");

#ifdef SISL_DEPENDENT_DEBUG
    // We write to file parametric trim cvs, as well as surfs and spatial trim cvs.
    vector<shared_ptr<BoundedSurface> > both_sfs;
    both_sfs.push_back(bounded_sf1);
    both_sfs.push_back(bounded_sf2);
    for (int ki = 0; ki < int(both_sfs.size()); ++ki) {
	std::ofstream debug("data/debug.g2");
	both_sfs[ki]->underlyingSurface()->writeStandardHeader(debug);
	both_sfs[ki]->underlyingSurface()->write(debug);
	vector<CurveLoop> all_bd_loops = both_sfs[ki]->allBoundaryLoops(epsge);
	for (int kj = 0; kj < int(all_bd_loops.size()); ++kj)
	    for (int kk = 0; kk < all_bd_loops[kj].size(); ++kk) {
		shared_ptr<CurveOnSurface> trim_cv =
		    dynamic_pointer_cast<CurveOnSurface>(all_bd_loops[kj][kk]);
		ASSERT(trim_cv.get() != 0); // All surfs on GoBOundedSurface is of this type.
		shared_ptr<SplineCurve> pcv =
		    dynamic_pointer_cast<SplineCurve, ParamCurve>(trim_cv->parameterCurve());
		ASSERT((pcv.get() != 0) && (pcv->dimension() == 2));
		writeSpaceParamCurve(*pcv, debug);
		trim_cv->spaceCurve()->writeStandardHeader(debug);
		trim_cv->spaceCurve()->write(debug);
	    }
    }
#endif // SISL_DEPENDENT_DEBUG

    vector<shared_ptr<CurveOnSurface> > int_segments1, int_segments2;
    double tol = 1e-04*epsge;
    shared_ptr<ParamSurface> under1 = under_sf1;
    shared_ptr<ParamSurface> under2 = under_sf2;
    try {
	getIntersectionCurve(under1, under2, int_segments1, int_segments2, tol);
    } catch (...) {
	THROW("Failed intersecting the two spline surfaces.");
    }
    vector<shared_ptr<BoundedSurface> > trimmed_sfs;
    if (int_segments1.size() == 0) { // Size of the two vectors is equal.
	return trimmed_sfs; // Empty intersection, return empty vector.
    }

    intersectWithSurfaces(int_segments1, bounded_sf1, int_segments2, bounded_sf2, 
			  0.1*epsge);

#ifdef SISL_DEPENDENT_DEBUG
    std::ofstream debug("data/debug.g2");
    under_sf1->writeStandardHeader(debug);
    under_sf1->write(debug);
    under_sf2->writeStandardHeader(debug);
    under_sf2->write(debug);
    for (int kj = 0; kj < int(int_segments1.size()); ++kj) {
	shared_ptr<CurveOnSurface> cv_on_sf1 = int_segments1[kj];
	cv_on_sf1->spaceCurve()->writeStandardHeader(debug);
	cv_on_sf1->spaceCurve()->write(debug);
	shared_ptr<SplineCurve> spline_cv =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf1->parameterCurve());
	writeSpaceParamCurve(*spline_cv, debug);
	shared_ptr<CurveOnSurface> cv_on_sf2 = int_segments2[kj];
	cv_on_sf2->spaceCurve()->writeStandardHeader(debug);
	cv_on_sf2->spaceCurve()->write(debug);
	spline_cv =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf2->parameterCurve());
	writeSpaceParamCurve(*spline_cv, debug);
	// @@sbr This is to removed as soon as possible!!!
	if (false) {
	    int erase_no = 1;
	    int_segments1.erase(int_segments1.begin() + erase_no);
	    int_segments2.erase(int_segments2.begin() + erase_no);
	}
    }
#endif // SISL_DEPENDENT_DEBUG

    // We pass the intersection curves through a function generating trim curve.
    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves1, loop_curves2;
    try {
      loop_curves1 = getBoundaryLoops(*bounded_sf1, int_segments1, epsge);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }
    vector<shared_ptr<BoundedSurface> > trim_sfs1 =
       createTrimmedSurfs(loop_curves1, under_sf1, epsge);
    trimmed_sfs.insert(trimmed_sfs.end(), trim_sfs1.begin(), trim_sfs1.end());

    try {
      loop_curves2 = getBoundaryLoops(*bounded_sf2, int_segments2, epsge);
    } catch (...) {
	MESSAGE("Failed extracting loops. Suspecting input curve ended in middle of surface.");
    }
    vector<shared_ptr<BoundedSurface> > trim_sfs2 =
       createTrimmedSurfs(loop_curves2, under_sf2, epsge);
    trimmed_sfs.insert(trimmed_sfs.end(), trim_sfs2.begin(), trim_sfs2.end());

    return trimmed_sfs;
}




//===========================================================================
BoundedSurface* BoundedUtils::convertToBoundedSurface(const SplineSurface& surf,
							  double space_epsilon)
//===========================================================================
{
    shared_ptr<ParamSurface> under_surf;
    under_surf = shared_ptr<ParamSurface>(surf.clone());
    vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(under_surf,
							 space_epsilon);

    return new BoundedSurface(under_surf, loops);
}

//===========================================================================
  shared_ptr<BoundedSurface> 
  BoundedUtils::convertToBoundedSurface(shared_ptr<ParamSurface> surf,
					double space_epsilon)
//===========================================================================
{
  shared_ptr<BoundedSurface> bounded_sf;
  if (surf->instanceType() == Class_BoundedSurface) 
    {
      bounded_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    }
  else
    {
      shared_ptr<ParamSurface> surf2(surf->clone());
      vector<CurveLoop> loops = 
	SurfaceTools::absolutelyAllBoundarySfLoops(surf2, space_epsilon);
      bounded_sf = shared_ptr<BoundedSurface>(new BoundedSurface(surf2, loops));
    }
  return bounded_sf;
}

//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::splitWithTrimSegments(shared_ptr<BoundedSurface> surf,
				    vector< shared_ptr< CurveOnSurface > >& bnd_cvs,
				    double eps)
//===========================================================================
{
  // Ensure that the input curve loop has got a parametric representation
  size_t ki, kj;
  for (ki=0; ki<bnd_cvs.size(); ++ki)
    bnd_cvs[ki]->ensureParCrvExistence(eps);

  // Remove duplicates
  for (ki=0; ki<bnd_cvs.size(); ++ki)
    for (kj=ki+1; kj<bnd_cvs.size(); ++kj)
      if (checkCurveCoincidence(bnd_cvs[ki], bnd_cvs[kj], eps, false))
	{
	  bnd_cvs.erase(bnd_cvs.begin()+kj);
	  kj--;
	}
  
  // Represent input curves with both orientations
  size_t nmb_crvs = bnd_cvs.size();
  for (ki=0; ki<nmb_crvs; ++ki)
    {
      shared_ptr<CurveOnSurface> tmp_crv = 
	shared_ptr<CurveOnSurface>(bnd_cvs[ki]->clone());
      tmp_crv->reverseParameterDirection();
      bnd_cvs.push_back(tmp_crv);
    }

  // First get boundary loops corresponding to the split surfaces
  vector< vector< shared_ptr< CurveOnSurface > > > new_loops = 
      getBoundaryLoops(*(surf.get()), bnd_cvs, eps, (int)nmb_crvs);

  // Make trimmed surfaces
  vector<shared_ptr<BoundedSurface> > trim_sfs =
    createTrimmedSurfs(new_loops, surf->underlyingSurface(), eps);

  return trim_sfs;
 }

//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::subtractSfPart(shared_ptr<BoundedSurface> surf,
			     vector< shared_ptr< CurveOnSurface > >& bnd_cvs,
			     double eps)
//===========================================================================
{
  // Ensure that the input curve loop has got a parametric representation
  size_t ki;
  for (ki=0; ki<bnd_cvs.size(); ++ki)
    bnd_cvs[ki]->ensureParCrvExistence(eps);

  // Compute parameter value inside loop
  Point midpar(0.0, 0.0);
  for (ki=0; ki<bnd_cvs.size(); ++ki)
    {
      double t1 = 0.5*(bnd_cvs[ki]->startparam() + bnd_cvs[ki]->endparam());
      Point mid = bnd_cvs[ki]->parameterCurve()->point(t1);
      midpar += mid;
    }
  midpar /= (int)(bnd_cvs.size());

  // To ensure that we get the correct part of the surface, we make all.
  // Represent input loop with both orientations
  size_t nmb_crvs = bnd_cvs.size();
  for (ki=0; ki<nmb_crvs; ++ki)
    {
      shared_ptr<CurveOnSurface> tmp_crv = 
	shared_ptr<CurveOnSurface>(bnd_cvs[ki]->clone());
      tmp_crv->reverseParameterDirection();
      bnd_cvs.push_back(tmp_crv);
    }

  // First get boundary loops corresponding to the split surfaces
  vector< vector< shared_ptr< CurveOnSurface > > > new_loops = 
      getBoundaryLoops(*(surf.get()), bnd_cvs, eps, (int)nmb_crvs);

  // Make trimmed surfaces
  vector<shared_ptr<BoundedSurface> > trim_sfs =
    createTrimmedSurfs(new_loops, surf->underlyingSurface(), eps);

  // Remove specified surface part

  for (ki=0; ki<trim_sfs.size(); ++ki)
  {
      if (trim_sfs[ki]->inDomain(midpar[0], midpar[1]))
      {
          trim_sfs.erase(trim_sfs.begin()+ki);
          break;
      }
  }

  return trim_sfs;
}


//===========================================================================
vector< vector< shared_ptr< CurveOnSurface > > >
BoundedUtils::getBoundaryLoops(const BoundedSurface& sf,
			       vector< shared_ptr< CurveOnSurface > >& part_bd_cvs,
			       double eps, int last_split)
//===========================================================================
{
    double a_tol = 1.0e-8;  
    int dim = sf.dimension();

    vector<vector<shared_ptr<CurveOnSurface> > > new_loops;

    if (part_bd_cvs.size() == 0)
	THROW("Curve input vector was empty!");
    int ki, kj;

    // Compute minimum length of the part boundary curves
    double min_part_bd_len = std::numeric_limits<double>::max();
    for (ki=0; ki<(int)part_bd_cvs.size(); ++ki)
      min_part_bd_len = std::min(min_part_bd_len, 
				 part_bd_cvs[ki]->estimatedCurveLength());
    min_part_bd_len *= 0.98; // Reduce slightly

    // Compute largest gap between part boundary curves
    // Fetch also endpoints of chains of part boundary curves
    vector<Point> part_bd_endpt;
    int nmb_part_bd_cvs = (int)part_bd_cvs.size()/2;  // Represented twice,
    // one for each orientation
    double part_bd_gap = maxCurveGap(&part_bd_cvs[0], nmb_part_bd_cvs,
				     (dim == 1), eps, part_bd_endpt);
    

    // We next must extract boundary loops and check for intersections.
    // Function returns all curves, including degenerate ones.
    vector<CurveLoop> boundary_loops = sf.absolutelyAllBoundaryLoops();

    // Knowing that surf is a BoundedSurface, we extract CurveOnSurface's.
    vector<shared_ptr<CurveOnSurface> > old_loop_cvs; // All loop curves stored in common vector.
    //double min_loop_tol = 1.5*boundary_loops[0].getSpaceEpsilon();
    double min_loop_tol = 1.1*boundary_loops[0].getMaxCurveDist();
    min_loop_tol = std::max(min_loop_tol, eps);

    shared_ptr<const ParamSurface> under_sf = sf.underlyingSurface();
//     double knot_diff_tol =
//       0.01*getParEps(min_loop_tol, under_sf.get()); // We may not trust pcv to repr space_cv.
    double knot_diff_tol =
      getParEps(min_loop_tol, under_sf.get()); // We may not trust pcv to repr space_cv.
    for (ki = 0; ki < int(boundary_loops.size()); ++ki) {
      min_loop_tol = std::min(min_loop_tol, 1.5*boundary_loops[ki].getSpaceEpsilon());
      min_loop_tol = std::max(min_loop_tol, 
			      1.1*boundary_loops[ki].getMaxCurveDist());
	for (kj = 0; kj < boundary_loops[ki].size(); ++kj)
	  {
	    // Ensure a certain degree of correspondence between curve lengths
	    // and parameterization
	    shared_ptr<ParamCurve> tmp_crv(boundary_loops[ki][kj]->clone());
	    shared_ptr<CurveOnSurface> tmp_crv2 = 
	      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(tmp_crv);
	    tmp_crv2->ensureParCrvExistence(min_loop_tol);
	    double cv_len = (dim == 1) ?
	      tmp_crv2->parameterCurve()->estimatedCurveLength() :
	      tmp_crv2->estimatedCurveLength();
	    double t1 = tmp_crv2->startparam();
	    double t2 = tmp_crv2->endparam();
	    double tdel = t2 - t1;
	    if (cv_len > tdel && cv_len > min_loop_tol) //&& t2 - t1 <= knot_diff_tol)
	      {
		tmp_crv2->setParameterInterval(t1, t1+cv_len);
		tdel = cv_len;
	      }
	    old_loop_cvs.push_back(tmp_crv2);
	    knot_diff_tol = std::min(knot_diff_tol, 0.2*tdel);
	  }
    }
    knot_diff_tol = std::max(knot_diff_tol, DEFAULT_PARAMETER_EPSILON);

#ifdef DEBUG1
    if (dim > 1)
      {
	std::ofstream out0_1("old_loop_cvs1.g2");
	for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
	  {
	    out0_1 << "100 1 0 4 255 0 0 255" << std::endl;
	    shared_ptr<SplineCurve> tmp1 = 
	      shared_ptr<SplineCurve>(old_loop_cvs[ki]->geometryCurve());
	    tmp1->write(out0_1);
	  }
      }
#endif

    double min_tool = 1.0e-5;
    min_loop_tol = std::max(min_tool, min_loop_tol);
    min_loop_tol = std::max(min_loop_tol, eps);
    double epspar = min_loop_tol; // Assuming parameter domain reflects the geometry...
    // Since we do not know anything about this and it is mainly a test additional
    // to the geometry space test to avoid confusing seam curves, we make it a little bigger
    Point par_eps = SurfaceTools::getParEpsilon(sf, min_loop_tol);
    epspar = std::max(min_loop_tol, 0.5*(par_eps[0]+par_eps[1]));
    epspar *= 10.0;

    // Check part_bd_cvs agains old_loop_cvs
    for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
      for (kj=0; kj<(int)part_bd_cvs.size(); ++kj)
	{
	  shared_ptr<ParamCurve> loop_cv = (dim == 1) ?
	    old_loop_cvs[ki]->parameterCurve() : old_loop_cvs[ki];
	  shared_ptr<ParamCurve> part_cv = (dim == 1) ?
	    part_bd_cvs[kj]->parameterCurve() : part_bd_cvs[kj];
	  int coincidence = checkCurveCoinc(loop_cv, part_cv,
					    /*1.1**/min_loop_tol);
	  if (coincidence == 1 || coincidence == 3)
	    {
	      // Coincidence or the part boundary curve is embedded in the
	      // loop curve. Remove part boundary curve
	      part_bd_cvs.erase(part_bd_cvs.begin() + kj);
	      if (last_split > 0 && kj < last_split)
		--last_split;
	      kj--;
	    }
	}

    if (part_bd_cvs.size() == 0)
      {
	// No need for rearrangement. Return existing boundary loop
	for (size_t kf=0; kf<boundary_loops.size(); ++kf)
	  {
	    vector<shared_ptr<CurveOnSurface> > bd_loop;
	    int bd_size = boundary_loops[kf].size();
	    for (int kk=0; kk<bd_size; ++kk)
	      {
		shared_ptr<CurveOnSurface> tmp_sfcv = 
		  dynamic_pointer_cast<CurveOnSurface,ParamCurve>(boundary_loops[kf][kk]);
		  bd_loop.push_back(tmp_sfcv);
	      }
	    new_loops.push_back(bd_loop);
	  }
	return new_loops;
      }

    // We run part_bd_cvs, splitting if they start/end in inner part of a cv in old_loop_cvs.
    splitLoopCvs(sf, old_loop_cvs, part_bd_cvs, part_bd_endpt, min_loop_tol,
		 eps, epspar, knot_diff_tol, last_split, (dim == 1));

#ifdef DEBUG1
    std::ofstream out1("loop_cvs1.g2");
    for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
      {
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(old_loop_cvs[ki]->geometryCurve());
	if (tmp1.get())
	  {
	    out1 << "100 1 0 4 255 0 0 255" << std::endl;
	    tmp1->write(out1);
	  }
      }
    for (ki=0; ki<(int)part_bd_cvs.size(); ++ki)
      {
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(part_bd_cvs[ki]->geometryCurve());
	if (tmp1.get())
	  {
	    out1 << "100 1 0 4 0 255 0 255" << std::endl;
	    tmp1->write(out1);
	  }
      }
    std::ofstream out1_2("loop_cvs_par.g2");
    for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
      {
	out1_2 << "100 1 0 4 255 0 0 255" << std::endl;
	shared_ptr<ParamCurve> tmp1 = 
	  shared_ptr<ParamCurve>(old_loop_cvs[ki]->parameterCurve());
	tmp1->write(out1_2);
      }
    for (ki=0; ki<(int)part_bd_cvs.size(); ++ki)
      {
	out1_2 << "100 1 0 4 0 255 0 255" << std::endl;
	shared_ptr<ParamCurve> tmp1 = 
	  shared_ptr<ParamCurve>(part_bd_cvs[ki]->parameterCurve());
	tmp1->write(out1_2);
      }
#endif // DEBUG1

    if (part_bd_cvs.size() == 0)
      {
	// No need for rearrangement. Return existing boundray loop
	for (size_t kf=0; kf<boundary_loops.size(); ++kf)
	  {
	    vector<shared_ptr<CurveOnSurface> > bd_loop;
	    int bd_size = boundary_loops[kf].size();
	    for (int kk=0; kk<bd_size; ++kk)
	      bd_loop.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(boundary_loops[kf][kk]));
	    new_loops.push_back(bd_loop);
	  }
	return new_loops;
      }

    // We are now left with the task of arranging meeting cv parts. Starting w/first new cv.
    // Curves may be given in any order.
    shared_ptr<CurveOnSurface> curr_crv = part_bd_cvs[0];
    vector<shared_ptr<CurveOnSurface> > curr_loop;
    curr_loop.push_back(curr_crv);
    part_bd_cvs.erase(part_bd_cvs.begin());
    double t1 = curr_crv->startparam();
    double t2 = curr_crv->endparam();
    Point total_par_start_pt = curr_crv->parameterCurve()->point(t1);
    Point total_space_start_pt = (dim == 1) ? total_par_start_pt :
      curr_crv->ParamCurve::point(t1);
    Point curr_par_end_pt = curr_crv->parameterCurve()->point(t2);
    Point curr_space_end_pt = (dim == 1) ? curr_par_end_pt :
      curr_crv->ParamCurve::point(t2);
    double space_end_dist = total_space_start_pt.dist(curr_space_end_pt);
    double par_end_dist = total_par_start_pt.dist(curr_par_end_pt);

#ifdef DEBUG1
    if (space_end_dist >= min_loop_tol && par_end_dist < epspar)
      {
	shared_ptr<const ParamSurface> psf = sf.underlyingSurface();
	Point geom1 = psf->point(total_par_start_pt[0], total_par_start_pt[1]);
	Point geom2 = psf->point(curr_par_end_pt[0], curr_par_end_pt[1]);
	double gdist = geom1.dist(geom2);
	int stop_break = 1;
      }
#endif
//     bool cw_loop;
//     double max_angle = -4.0; // Angles are measured on (-pi, pi].
    // We must locate max angle among remaining parts, both new and old.
//     bool loop_locally_closed = false; // We may want to continue after loop is closed.
//     double local_angle = -2*M_PI; // Illegal value.
    // double min_ang = 0.01;

    double min_loop_tol2 = min_loop_tol;
    int prev_part_ind = -1, prev_old_ind = -1;
    bool first_is_part_bd = true;
    while (true) { 
      vector<shared_ptr<CurveOnSurface> > tmp_vec;
      if (nmb_part_bd_cvs > 1 || prev_old_ind >= 0)
	tmp_vec.insert(tmp_vec.end(), part_bd_cvs.begin(), part_bd_cvs.end());
      int curr_part_size = (int)tmp_vec.size();
      tmp_vec.insert(tmp_vec.end(), old_loop_cvs.begin(), old_loop_cvs.end());
      double angtol = 0.01;
      double tmp_ang;
      int part_ind = -1, old_ind = -1;
      int tmp_ind;
#ifdef DEBUG1
      std::ofstream out_tmp("loop_cvs_tmp.g2");
      for (size_t ix=0; ix<curr_loop.size(); ++ix)
	{
	  if (curr_loop[ix]->spaceCurve().get() && dim != 1)
	    {
	      curr_loop[ix]->spaceCurve()->writeStandardHeader(out_tmp);
	      curr_loop[ix]->spaceCurve()->write(out_tmp);
	    }
	  else
	    {
	      curr_loop[ix]->parameterCurve()->writeStandardHeader(out_tmp);
	      curr_loop[ix]->parameterCurve()->write(out_tmp);
	    }
	}
#endif // DEBUG1
      tmp_ind = leftMostCurve(*curr_crv, tmp_vec, min_loop_tol2, 
			      (dim == 1), tmp_ang);
      if ((fabs(2*M_PI-tmp_ang) < angtol || tmp_ind < 0) && tmp_vec.size() > 0)
	{
	  // Try again with a larger space tolerance
	  tmp_ind = leftMostCurve(*curr_crv, tmp_vec, 
				  std::min(10.0*min_loop_tol2, min_part_bd_len), 
				  (dim == 1), tmp_ang);
	}
	  
      double part_angle = (tmp_ind < curr_part_size) ? tmp_ang : 8.0;
      if (tmp_ind < curr_part_size)
	part_ind = tmp_ind;
      double old_angle = (tmp_ind >= curr_part_size) ? tmp_ang : 8.0;
      if (tmp_ind >= curr_part_size)
	old_ind = tmp_ind - curr_part_size;

	if (space_end_dist >= min_loop_tol && 
	    space_end_dist < 50.0*min_loop_tol && 
	    par_end_dist < epspar /*&& part_ind >= 0*/)
	  {
	    // Make an extra check to avoid bypassing short connecting curves
	    double cv_len = std::numeric_limits<double>::max(); 
	    if (tmp_ind >= 0) 
	      cv_len = (dim == 1) ? 
		tmp_vec[tmp_ind]->parameterCurve()->estimatedCurveLength() :
		tmp_vec[tmp_ind]->estimatedCurveLength();

	    // Check if the start- and endpoints correspond to loose ends
	    // of the part boundary curves
	    size_t kr, kh;
	    for (kr=0; kr<part_bd_endpt.size(); ++kr)
	      {
		if (total_space_start_pt.dist(part_bd_endpt[kr]) < min_loop_tol)
		  break;
	      }
	    for (kh=0; kh<part_bd_endpt.size(); ++kh)
	      {
		if (curr_space_end_pt.dist(part_bd_endpt[kh]) < min_loop_tol)
		  break;
	      }
	    
	    double cand_end_dist = std::numeric_limits<double>::max();
	    if (cv_len < 50.0*min_loop_tol)
	      {
		// Check endpoint of candidate curve
		Point cand_end = (dim == 1) ?
		  tmp_vec[tmp_ind]->parameterCurve()->point(tmp_vec[tmp_ind]->endparam()) :
		  tmp_vec[tmp_ind]->ParamCurve::point(tmp_vec[tmp_ind]->endparam());
		cand_end_dist = total_space_start_pt.dist(cand_end);
	      }
	    if ((space_end_dist < min_loop_tol + part_bd_gap ||
		space_end_dist < 2.0*min_loop_tol) && 
		cand_end_dist > space_end_dist)
	      {
		// Probably a closed loop despite a big gap
		min_loop_tol2 = std::min(min_part_bd_len,
					 space_end_dist + a_tol);
	      }
	    else if (space_end_dist < 10.0*min_loop_tol &&
		     (kr < part_bd_endpt.size() || kh < part_bd_endpt.size()))
	      
	      {
		// A missing connection in the part boundary curves,
		// assumes a closed loop
		min_loop_tol2 = std::min(min_part_bd_len,
					 space_end_dist + a_tol);
	      }
	    else if (cand_end_dist > space_end_dist)
	      {
		// Make a test run for the next curve
		vector<shared_ptr<CurveOnSurface> > tmp_vec2;
		if (nmb_part_bd_cvs > 1 || prev_old_ind >= 0)
		  tmp_vec2.insert(tmp_vec2.end(), part_bd_cvs.begin(), 
				  part_bd_cvs.end());
		int curr_part_size2 = (int)tmp_vec2.size();
		tmp_vec2.insert(tmp_vec2.end(), old_loop_cvs.begin(), 
				old_loop_cvs.end());
		double tmp_ang2;
		int tmp_ind2 = leftMostCurve(*curr_crv, tmp_vec2, 
					     min_loop_tol, (dim == 1), tmp_ang2);
		if (tmp_ind2 < curr_part_size2 && 
		    fabs(2.0*M_PI - tmp_ang2) < angtol)
		  {
		    // Probably a closed loop despite a big gap
		    min_loop_tol2 = std::min(min_part_bd_len,
					     space_end_dist + a_tol);
		  }
	      }
	    int stop_break = 1;
	  }

	if (space_end_dist < min_loop_tol2 && par_end_dist < epspar) 
	  {
	    // Check if next segment is to the left of the first in curr_loop.
	    //tmp_vec.push_back(curr_loop.front());
	    tmp_vec.insert(tmp_vec.begin(), curr_loop.front());
	    double angle;
	    tmp_ind = leftMostCurve(*curr_crv, tmp_vec, min_loop_tol2, 
				    (dim == 1), angle);
	    //if (tmp_ind == (int)tmp_vec.size() - 1)
	    if (tmp_ind == 0)
	      part_ind = old_ind = -1;
	    else
	      {
		// A new segment to the left of the first segment is found.
		// Be careful not to join two loops in tangential cases
		// Check also curve lengths
		double len = (dim == 1) ?
		  tmp_vec[tmp_ind]->parameterCurve()->estimatedCurveLength() :
		  tmp_vec[tmp_ind]->estimatedCurveLength();
		if ((first_is_part_bd && prev_old_ind >= 0) ||
		    len < min_loop_tol2)
		  part_ind = old_ind = -1;
	      }
	  }

	if (part_ind != -1 && old_ind != -1) 
	  { // We choose the smallest angle.
	    if (part_angle <= old_angle)
		old_ind = -1;
	    else
		part_ind = -1;
	  } 
	else if (part_ind == -1 && old_ind == -1) { // No new segment.
	  if (space_end_dist < min_loop_tol2 /*&& par_end_dist < epspar*/ ||
	      (space_end_dist < std::min(10.0*min_loop_tol2, min_part_bd_len) 
	       && par_end_dist < epspar)) 
	    {
#ifdef DEBUG1
		std::ofstream out2("loop_cvs2.g2");
		for (size_t ix=0; ix<curr_loop.size(); ++ix)
		  {
		    if (curr_loop[ix]->spaceCurve().get())
		      {
			curr_loop[ix]->spaceCurve()->writeStandardHeader(out_tmp);
			curr_loop[ix]->spaceCurve()->write(out_tmp);
		      }
		    else
		      {
			curr_loop[ix]->parameterCurve()->writeStandardHeader(out_tmp);
			curr_loop[ix]->parameterCurve()->write(out_tmp);
		      }
		  }
#endif // DEBUG1

		new_loops.push_back(curr_loop);

	    } // If not this was a dead end ...
	    curr_loop.clear(); // No more segments in loop.
	    // We then look for start segment of next loop.
	    if (part_bd_cvs.size() != 0) {
		part_ind = 0;
		first_is_part_bd = true;
	    } else if (old_loop_cvs.size() != 0) {
	       part_ind = -1;
	       old_ind = 0;
	       first_is_part_bd = false;
	    } else {
		break; // No more cvs. 
	    }
	}
	if (part_ind != -1) {
	    curr_crv = part_bd_cvs[part_ind];
	    part_bd_cvs.erase(part_bd_cvs.begin() + part_ind);
	} else {
	    curr_crv = old_loop_cvs[old_ind];
	    old_loop_cvs.erase(old_loop_cvs.begin() + old_ind);
	}
	t1 = curr_crv->startparam();
	t2 = curr_crv->endparam();
	if (curr_loop.size() == 0) {
	  total_par_start_pt = curr_crv->parameterCurve()->point(t1);
	  total_space_start_pt = (dim == 1) ? total_par_start_pt :
	    curr_crv->ParamCurve::point(t1);
	}

	curr_loop.push_back(curr_crv);
	curr_par_end_pt = curr_crv->parameterCurve()->point(t2);
	curr_space_end_pt = (dim == 1) ? curr_par_end_pt :
	  curr_crv->ParamCurve::point(t2);
	space_end_dist = total_space_start_pt.dist(curr_space_end_pt);
	par_end_dist = total_par_start_pt.dist(curr_par_end_pt);

	//min_loop_tol2 = min_loop_tol;
	prev_part_ind = part_ind;
	prev_old_ind = old_ind;
    }

    return new_loops; // We should be done
}


//===========================================================================
int BoundedUtils::checkCurveCoinc(shared_ptr<ParamCurve> cv1, 
				  shared_ptr<ParamCurve> cv2, 
				  double tol)
//===========================================================================
{
  int coinc = 0;  // Initially no coincidence

  // Check if two curves are identical, or one is embedded in the other.
  // Compare endpoints
  int dim = cv1->dimension();
  
  shared_ptr<ParamCurve> tmp1 = cv1;
  shared_ptr<ParamCurve> tmp2 = cv2;
  if (dim == 1)
    {
      // Check existence of parameter curve
      shared_ptr<CurveOnSurface> sf_cv1 = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
      shared_ptr<CurveOnSurface> sf_cv2 = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
      if ((!sf_cv1.get()) || (!sf_cv2.get()))
	return coinc;
      if ((!sf_cv1->hasParameterCurve()) || (!sf_cv2->hasParameterCurve()))
	return coinc;  // Not possible to check
      tmp1 = sf_cv1->parameterCurve();
      tmp2 = sf_cv2->parameterCurve();
    }
  Point pt1 = tmp1->point(tmp1->startparam());
  Point pt2 = tmp1->point(tmp1->endparam());
  Point pt3 = tmp2->point(tmp2->startparam());
  Point pt4 = tmp2->point(tmp2->endparam());
  double d1 = pt1.dist(pt3);
  double d2 = pt1.dist(pt4);
  double d3 = pt2.dist(pt3);
  double d4 = pt2.dist(pt4);

  double param1[2], param2[2];  // Endpoint of coincidence interval
  if ((d1 < tol && d4 < tol) || (d2 < tol && d3 < tol))
    {
      // Potential full coincidence
      param1[0] = tmp1->startparam();
      param1[1] = tmp1->endparam();
      param2[0] = tmp2->startparam();
      param2[1] = tmp2->endparam();
      coinc = 1;
    }
  else
    {
      // Check if cv1 is embedded in cv2
      Point clo1, clo2, clo3;
      double par1, par2, par3, dist1, dist2, dist3;
      Point mid = tmp1->point(0.5*(tmp1->startparam()+tmp1->endparam()));
      tmp2->closestPoint(pt1, tmp2->startparam(), tmp2->endparam(),
			par1, clo1, dist1);
      tmp2->closestPoint(pt2, tmp2->startparam(), tmp2->endparam(),
			par2, clo2, dist2);
      tmp2->closestPoint(mid, tmp2->startparam(), tmp2->endparam(),
			par3, clo3, dist3);
      if (dist1 < tol && dist2 < tol && dist3 < tol)
	{
	  param1[0] = tmp1->startparam();
	  param1[1] = tmp1->endparam();
	  param2[0] = std::min(par1, par2);
	  param2[1] = std::max(par1, par2);
	  coinc = 2;
	}
      else
	{
	  // Check of cv2 is embedded in cv1
	  mid = tmp2->point(0.5*(tmp2->startparam()+tmp2->endparam()));
	  tmp1->closestPoint(pt3, tmp1->startparam(), tmp1->endparam(),
			    par1, clo1, dist1);
	  tmp1->closestPoint(pt4, tmp1->startparam(), tmp1->endparam(),
			    par2, clo2, dist2);
	  tmp1->closestPoint(mid, tmp1->startparam(), tmp1->endparam(),
			    par3, clo3, dist3);
 	  if (dist1 < tol && dist2 < tol && dist3 < tol)
	    {
	      param1[0] = std::min(par1, par2);
	      param1[1] = std::max(par1, par2);
	      param2[0] = tmp2->startparam();
	      param2[1] = tmp2->endparam();
	      coinc = 3;
	    }
	}
    }
 	  
  if (coinc == 0)
    return coinc;  // No coincidence


  // Check the inner of the curve
  int nmb_sample = 3;
  double tdel = (param1[1] - param1[0])/(double)(nmb_sample+1);
  double tpar;
  int ki;
  for (ki=0, tpar=param1[0]+tdel; ki<nmb_sample; ++ki, tpar+=tdel)
    {
      Point clo;
      double par, dist;
      Point pos = tmp1->point(tpar);
      tmp2->closestPoint(pos, param2[0], param2[1], par, clo, dist);
      if (dist >= tol)
	return 0;
    }

  return coinc;
}

//===========================================================================
int BoundedUtils::checkCurveCoincidence(shared_ptr<CurveOnSurface> cv1, 
					shared_ptr<CurveOnSurface> cv2, double tol,
					bool same_orient)
//===========================================================================
{
  int nmb_sample = 2;  // Number of points to check in the inner of a candidate curve

  int dim = cv1->dimension();
  
  shared_ptr<ParamCurve> tmp1 = cv1;
  shared_ptr<ParamCurve> tmp2 = cv2;
  if (dim == 1)
    {
      if ((!cv1->hasParameterCurve()) || (!cv2->hasParameterCurve()))
	return 0;  // Not possible to check
      tmp1 = cv1->parameterCurve();
      tmp2 = cv2->parameterCurve();
    }

  // Evaluate endpoints
  Point pt1 = tmp1->ParamCurve::point(tmp1->startparam());
  Point pt2 = tmp1->ParamCurve::point(tmp1->endparam());
  Point pt3 = tmp2->ParamCurve::point(tmp2->startparam());
  Point pt4 = tmp2->ParamCurve::point(tmp2->endparam());
  double d1 = pt1.dist(pt3);
  double d2 = pt1.dist(pt4);
  double d3 = pt2.dist(pt3);
  double d4 = pt2.dist(pt4);
  if (same_orient)
    {
      if (d1 >= tol || d4 >= tol)
	return 0;
    }
  else
    {
      if (!((d1 < tol && d4 < tol) ||
	    (d2 < tol && d3 < tol)))
	return 0; // Not a candidate for coincidence
    }

  // Check the inner of the curves
  int kr;
  double t1 = tmp1->startparam();
  double t2 = tmp1->endparam();
  double tdel = (t2 - t1)/(double)(nmb_sample+1);
  double tpar;
  for (kr=0, tpar=t1+tdel; kr<nmb_sample; ++kr, tpar+=tdel)
    {
      Point pt5 = tmp1->ParamCurve::point(tpar);
      double clo_par, clo_dist;
      Point clo_pt;
      tmp2->closestPoint(pt5, tmp2->startparam(), tmp2->endparam(), 
			clo_par, clo_pt, clo_dist);
      if (clo_dist >= tol)
	return 0;
    }

  return 1;
}

//===========================================================================
vector<shared_ptr<BoundedSurface> >
BoundedUtils::createTrimmedSurfs(vector<vector<shared_ptr<CurveOnSurface> > >&
				 loops,
				 shared_ptr<ParamSurface> under_sf,
				 double epsgeo)
//===========================================================================
{
   int ki, kj;
   vector<shared_ptr<BoundedSurface> > return_sfs;
   vector<vector<shared_ptr<CurveOnSurface> > > cw_loops, ccw_loops;
   // double int_tol = 1e-08;
   for (ki = 0; ki < int(loops.size()); ++ki) {
     // Check for degenerate loops. Loops that represent lines or 
     // points are dismissed
     double deg_eps = epsgeo;
     if (loops[ki].size() == 3)
       {
	 // Use a small tolerance to avoid dismissing loops with a small angle
	 deg_eps = 1.0e-6;
       }
     if (loopIsDegenerate(loops[ki], deg_eps))
       continue;

     // Sort according to loop orientation
     if (LoopUtils::paramIsCCW(loops[ki], epsgeo, epsgeo)) {
     //if (LoopUtils::paramIsCCW(loops[ki], epsgeo, int_tol)) {
	 ccw_loops.push_back(loops[ki]);
      } else {
	 cw_loops.push_back(loops[ki]);
      }
   }

   // We run through both types, making sure that the outer trim
   // curves are ordered by inclusion (i.e. the smallest one is first).
   for (ki = 0; ki < int(ccw_loops.size()); ++ki) {
       for (kj = ki + 1; kj < int(ccw_loops.size()); ++kj) {
           if (!LoopUtils::firstLoopInsideSecond(ccw_loops[ki], ccw_loops[kj],
						 epsgeo, epsgeo /*int_tol*/)) {
               vector<shared_ptr<CurveOnSurface> > tmp_loop = ccw_loops[kj];
	       ccw_loops.insert(ccw_loops.begin() + ki, tmp_loop);
	       ccw_loops.erase(ccw_loops.begin() + kj + 1);
// 	       ccw_loops.erase(ccw_loops.begin() + kj);
// 	       --kj;
	   }//  else if (LoopUtils::firstLoopInsideSecond(ccw_loops[kj],
// 						       ccw_loops[ki],
// 						       epsgeo, int_tol)) {
// 	       ccw_loops.erase(ccw_loops.begin() + ki);
// 	       --ki;
// 	       break;
// 	   }
       }
   }
   // The same applies to these loops.
   for (ki = 0; ki < int(cw_loops.size()); ++ki) {
       for (kj = ki + 1; kj < int(cw_loops.size()); ++kj) {
	   if (!LoopUtils::firstLoopInsideSecond(cw_loops[ki], cw_loops[kj],
						 epsgeo, epsgeo /*int_tol*/)) {
               vector<shared_ptr<CurveOnSurface> > tmp_loop = cw_loops[kj];
	       cw_loops.insert(cw_loops.begin() + ki, tmp_loop);
	       cw_loops.erase(cw_loops.begin() + kj + 1);
// 	       cw_loops.erase(cw_loops.begin() + kj);
// 	       --kj;
	   }//  else if (LoopUtils::firstLoopInsideSecond(cw_loops[ki],
// 						       cw_loops[kj],
// 						       epsgeo, int_tol)) {
// 	       cw_loops.erase(cw_loops.begin() + ki);
// 	       --ki;
// 	       break;
// 	   }
       }
   }


#ifdef DEBUG1
   std::ofstream out("trim_sfs.g2");
#endif // DEBUG1

   // For each member of ccw_loops, we extract those cw_loops which lie inside.
   for (ki = 0; ki < int(ccw_loops.size()); ++ki) {
      vector<vector<shared_ptr<CurveOnSurface> > > loops2;
      loops2.push_back(ccw_loops[ki]);
      for (kj = 0; kj < int(cw_loops.size()); ++kj) {
	bool found;
	try {
	  found = LoopUtils::firstLoopInsideSecond(cw_loops[kj], ccw_loops[ki],
						   epsgeo, epsgeo /*int_tol*/);
	}
	catch (...)
	  {
	    continue;
	  }
	 if (found) {
	    loops2.push_back(cw_loops[kj]);
	    cw_loops.erase(cw_loops.begin() + kj);
	    --kj;
	 }
      }

      shared_ptr<BoundedSurface> surf =
	shared_ptr<BoundedSurface>(new BoundedSurface(under_sf, loops2, epsgeo));
#ifdef DEBUG1
      surf->writeStandardHeader(out);
      surf->write(out);
#endif // DEBUG1

      return_sfs.push_back(surf);
   }

   if (cw_loops.size() != 0) {
      MESSAGE("Not all input (cw) loops were used.");
   }

   return return_sfs;
}


//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::intersectWithPlane(shared_ptr<ParamSurface>& surf,
				   Point pnt, Point normal, double geom_tol)
//===========================================================================
{
    vector<shared_ptr<CurveOnSurface> > curves;

    // Convert the surface to a SISLSurf in order to use SISL functions
    // on it. The "false" argument dictates that the SISLSurf will only
    // copy pointers to arrays, not the arrays themselves.
    shared_ptr<SplineSurface> tmp_sf;
    SplineSurface *splinesf = surf->getSplineSurface();
    if (!splinesf)
      {
	tmp_sf = shared_ptr<SplineSurface>(surf->asSplineSurface());
	splinesf = tmp_sf.get();
      }
    if (splinesf == 0)
      THROW("Requiring surface to be a SplineSurface.");
    SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
    int dim = 3;
    double epsco = 1e-15; // Not used
//     double epsge = 1e-6;
    int numintpt;
    double* pointpar = 0;
    int numintcr;
    SISLIntcurve** intcurves = 0;
    int stat;
    // Find the topology of the intersection
    s1851(sislsf, pnt.begin(), normal.begin(), dim, epsco, geom_tol,
	  &numintpt, &pointpar, &numintcr, &intcurves, &stat);
    // @@sbr Not sure this is the right solution. Maybe stat!=0 because of warning.
    ALWAYS_ERROR_IF(stat<0,
		"s1851 returned code: " << stat);
#ifdef DEBUG1
    if (stat > 0)
      {
	std::cout << "s1851: " << stat << std::endl;
      }
#endif
    // pointpar is not used any further
    free(pointpar);
    double maxstep = 0.0;
    int makecurv = 2;     // Make both geometric and parametric curves
    int graphic = 0;      // Do not draw the curve
//     epsge = tol_.neighbour;
    for (int i = 0; i < numintcr; ++i) {
	// March out the intersection curves
	s1314(sislsf, pnt.begin(), normal.begin(), dim, epsco, geom_tol,
	      maxstep, intcurves[i], makecurv, graphic, &stat);
	SISLCurve* sc = intcurves[i]->pgeom;
	if (sc == 0) {
	    MESSAGE("s1314 returned code: " << stat << ", returning.");
	    continue;
	    // freeIntcrvlist(intcurves, numintcr);
	    // freeSurf(sislsf);
	    // return curves;
	}
	double* t = sc->et;
	double* c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	// Convert the geometric curve to Go format
	shared_ptr<ParamCurve> gcv(new SplineCurve(sc->in, sc->ik, t, c, 3));
	sc = intcurves[i]->ppar1;
	t = sc->et;
	c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	shared_ptr<ParamCurve> pcv(new SplineCurve(sc->in, sc->ik, t, c, 2));
	// We prefer parameter curves.
	// curves.push_back(shared_ptr<CurveOnSurface>
	// 		 (new CurveOnSurface(surf, pcv, gcv, true)));
	curves.push_back(shared_ptr<CurveOnSurface>
			 (new CurveOnSurface(surf, pcv, gcv, false)));

	// We would like the curve to have direction in consistency with trimming.
	vector<Point> derivs(2);
	gcv->point(derivs, gcv->startparam(), 1);
	Point par_pt = pcv->point(pcv->startparam());
	Point surf_normal;
	surf->normal(surf_normal, par_pt[0], par_pt[1]);
	Point dir_ind = surf_normal%derivs[1]; // Cross product.
	// If angle(dir_ind, normal) < 90, turn curve (should be in (90, 180).
	double inner_product = dir_ind*normal;
	if (inner_product > 0)
	  curves[curves.size()-1]->reverseParameterDirection();
    }
    freeIntcrvlist(intcurves, numintcr);
    freeSurf(sislsf);

    return curves;
}

//===========================================================================
vector<pair<Point, Point> >
BoundedUtils::intersectWithLine(shared_ptr<ParamSurface>& surf,
				Point pnt, Point dir, double geom_tol)
//===========================================================================
{
  vector<pair<Point,Point> > result;

  shared_ptr<SplineSurface> tmp_sf;
  SplineSurface *splinesf = surf->getSplineSurface();
  if (!splinesf)
    {
      tmp_sf = shared_ptr<SplineSurface>(surf->asSplineSurface());
      splinesf = tmp_sf.get();
    }
  if (splinesf == 0)
    {
      MESSAGE("Requiringsurface to be a SplineSurface.");
      return result;
    }
  SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
  int dim = 3;
  double epsco = 1e-15; // Not used
  //     double epsge = 1e-6;
  int numintpt;
  double* pointpar = 0; // array containing the parameter values of single intersect. pt.
  int numintcr; // number of intersection curves
  SISLIntcurve** intcurves = 0;
  int stat;

  // Find the intersection points
  s1856(sislsf, pnt.begin(), dir.begin(), dim, epsco, geom_tol,
	&numintpt, &pointpar, &numintcr, &intcurves, &stat);
  MESSAGE_IF(stat!=0, "s1856 returned code: " << stat);

  int i;
  for (i = 0; i < numintpt; ++i)
    {
      double u = pointpar[2*i];
      double v = pointpar[2*i+1];
      Point pt = surf->point(u, v);

      result.push_back(std::make_pair(Point(u,v), pt));
    }

  for (i=0; i<numintcr; ++i)
    {
      int npt = intcurves[i]->ipoint;
      double u = 0.5*(intcurves[i]->epar1[0]+intcurves[i]->epar1[2*(npt-1)]);
      double v = 0.5*(intcurves[i]->epar1[1]+intcurves[i]->epar1[2*npt-1]);
      Point pt = surf->point(u, v);

      result.push_back(std::make_pair(Point(u,v), pt));
    }
  if (sislsf)
    freeSurf(sislsf);
  if (pointpar)
    free(pointpar);
  if (intcurves)
    freeIntcrvlist(intcurves, numintcr);

  return result;
}

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::intersectWithCylinder(shared_ptr<ParamSurface>& surf,
				    Point pnt, Point vec, double radius, 
				    double geom_tol)
//===========================================================================
{
    vector<shared_ptr<CurveOnSurface> > curves;

    // Convert the surface to a SISLSurf in order to use SISL functions
    // on it. The "false" argument dictates that the SISLSurf will only
    // copy pointers to arrays, not the arrays themselves.
    shared_ptr<SplineSurface> tmp_sf;
    SplineSurface *splinesf = surf->getSplineSurface();
    if (!splinesf)
      {
	tmp_sf = shared_ptr<SplineSurface>(surf->asSplineSurface());
	splinesf = tmp_sf.get();
      }
    ALWAYS_ERROR_IF(splinesf == 0,
		    "Requiringsurface to be a SplineSurface.");
    SISLSurf* sislsf = GoSurf2SISL(*splinesf, false);
    int dim = 3;
    double epsco = 1e-15; // Not used
//     double epsge = 1e-6;
    int numintpt;
    double* pointpar = 0;
    int numintcr;
    SISLIntcurve** intcurves = 0;
    int stat;
    // Find the topology of the intersection
    s1853(sislsf, pnt.begin(), vec.begin(), radius, dim, epsco, geom_tol,
	  &numintpt, &pointpar, &numintcr, &intcurves, &stat);
    // @@sbr Not sure this is the right solution. Maybe stat!=0 because of warning.
    ALWAYS_ERROR_IF(stat!=0,
		"s1851 returned code: " << stat);
    // pointpar is not used any further
    free(pointpar);
    double maxstep = 0.0;
    int makecurv = 2;     // Make both geometric and parametric curves
    int graphic = 0;      // Do not draw the curve
//     epsge = tol_.neighbour;
     for (int i = 0; i < numintcr; ++i) {
	// March out the intersection curves
       s1316(sislsf,pnt.begin(), vec.begin(), radius, dim, epsco, geom_tol,
	      maxstep, intcurves[i], makecurv, graphic, &stat);
	SISLCurve* sc = intcurves[i]->pgeom;
	if (sc == 0) {
	    MESSAGE("s1316 returned code: " << stat << ", returning.");
	    continue;
	    // freeIntcrvlist(intcurves, numintcr);
	    // freeSurf(sislsf);
	    // return curves;
	}
	double* t = sc->et;
	double* c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	// Convert the geometric curve to Go format
	shared_ptr<ParamCurve> gcv(new SplineCurve(sc->in, sc->ik, t, c, 3));
	sc = intcurves[i]->ppar1;
	t = sc->et;
	c = (sc->ikind==2 || sc->ikind==4)? sc->rcoef : sc->ecoef;
	shared_ptr<ParamCurve> pcv(new SplineCurve(sc->in, sc->ik, t, c, 2));
	// We prefer parameter curves.
	// curves.push_back(shared_ptr<CurveOnSurface>
	// 		 (new CurveOnSurface(surf, pcv, gcv, true)));
	curves.push_back(shared_ptr<CurveOnSurface>
			 (new CurveOnSurface(surf, pcv, gcv, false)));

	// We would like the curve to have direction in consistency with trimming.
	vector<Point> derivs(2);
	gcv->point(derivs, gcv->startparam(), 1);
	Point par_pt = pcv->point(pcv->startparam());
	Point surf_normal;
	surf->normal(surf_normal, par_pt[0], par_pt[1]);
	Point dir_ind = surf_normal%derivs[1]; // Cross product.
	// If angle(dir_ind, normal) < 90, turn curve (should be in (90, 180).
	double inner_product = dir_ind*vec;
	if (inner_product > 0)
	  curves[curves.size()-1]->reverseParameterDirection();
    }
    freeIntcrvlist(intcurves, numintcr);
    freeSurf(sislsf);

    return curves;
}

//===========================================================================
void
BoundedUtils::getIntersectionCurve(shared_ptr<ParamSurface>& sf1,
				   shared_ptr<ParamSurface>& sf2,
				   vector<shared_ptr<CurveOnSurface> >& int_segments1,
				   vector<shared_ptr<CurveOnSurface> >& int_segments2,
				   double epsgeo)
//===========================================================================
{
    // @@sbr epsgeo should be close to 1e-06;
    ALWAYS_ERROR_IF((int_segments2.size() != 0) || (int_segments2.size() != 0),
		"Segment vectors must be empty!");

    shared_ptr<SplineSurface> tmp_sf1, tmp_sf2;
    SplineSurface *surf1 = sf1->getSplineSurface();
    SplineSurface *surf2 = sf2->getSplineSurface();
    if (!surf1)
      {
	tmp_sf1 = shared_ptr<SplineSurface>(sf1->asSplineSurface());
	surf1 = tmp_sf1.get();
      }
    if (!surf2)
      {
	tmp_sf2 = shared_ptr<SplineSurface>(sf2->asSplineSurface());
	surf2 = tmp_sf2.get();
      }
    if (!(surf1 && surf2))
      THROW("No underlying spline surface");

    // Check if the surface is closed
    bool close1_u, close1_v, close2_u, close2_v;
    SurfaceTools::checkSurfaceClosed(*surf1, close1_u, close1_v);
    SurfaceTools::checkSurfaceClosed(*surf2, close2_u, close2_v);

    SISLSurf* sisl_sf1 = GoSurf2SISL(*surf1);
    SISLSurf* sisl_sf2 = GoSurf2SISL(*surf2);
    if (close1_u)
      sisl_sf1->cuopen_1 = 0;
    if (close1_v)
      sisl_sf1->cuopen_2 = 0;
    if (close2_u)
      sisl_sf2->cuopen_1 = 0;
    if (close2_v)
      sisl_sf2->cuopen_2 = 0;

    double epsco = epsgeo; // Actually not used.
    int nmb_int_pts;
    double* pointpar1 = 0;
    double* pointpar2 = 0;
    int nmb_int_cvs;
    SISLIntcurve** intcurves = 0;
    int status = 0;
    int ki;
    s1859(sisl_sf1, sisl_sf2, epsco, epsgeo, &nmb_int_pts,
	  &pointpar1, &pointpar2, &nmb_int_cvs, &intcurves, &status);
    if (status < 0)
      {
	if (sisl_sf1) freeSurf(sisl_sf1);
	if (sisl_sf2) freeSurf(sisl_sf2);
	if (pointpar1) free(pointpar1);
	if (pointpar2) free(pointpar2);
	if (intcurves)
	  freeIntcrvlist(intcurves, nmb_int_cvs);
	THROW("Failed intersecting surfs.");
      }
    MESSAGE_IF(status != 0, "Returned status value: " << status);
    if (nmb_int_cvs == 0)
      {
	if (nmb_int_pts > 0)
	  {
	    // Check accuracy in intersection points
	    double mindist = std::numeric_limits<double>::max();
	    for (ki=0; ki<nmb_int_pts; ++ki)
	      {
		Point pos1 = surf1->ParamSurface::point(pointpar1[2*ki],
							pointpar1[2*ki+1]);
		Point pos2 = surf2->ParamSurface::point(pointpar2[2*ki],
							pointpar2[2*ki+1]);
		double dd = pos1.dist(pos2);
		mindist = std::min(mindist, dd);
	      }

	    // Rerun with a reduced tolerance
	    double epsgeo2 = std::max(0.5*mindist, 1.e-6);
	    s1859(sisl_sf1, sisl_sf2, epsco, epsgeo2, &nmb_int_pts,
		  &pointpar1, &pointpar2, &nmb_int_cvs, &intcurves, &status);
	    if (status < 0)
	      {
		if (sisl_sf1) freeSurf(sisl_sf1);
		if (sisl_sf2) freeSurf(sisl_sf2);
		if (pointpar1) free(pointpar1);
		if (pointpar2) free(pointpar2);
		if (intcurves)
		  freeIntcrvlist(intcurves, nmb_int_cvs);
		THROW("Failed intersecting surfs.");
	      }
	  }

	if (nmb_int_cvs == 0)
	  {
	    if (sisl_sf1) freeSurf(sisl_sf1);
	    if (sisl_sf2) freeSurf(sisl_sf2);
	    if (pointpar1) free(pointpar1);
	    if (pointpar2) free(pointpar2);
	    return; // Surfaces did not intersect.
	  }
      }

    double maxstep = (double)0;
    int makecurv = 2; // Make both geometric and parametric curves.
    int draw = 0;
    for (ki = 0; ki < nmb_int_cvs; ++ki) {
      // Set marching tolerance
      double march_eps = epsgeo;

      // Allow for some slack if the distance in the guide points are
      // close to the tolerance
      int nmb_guide = intcurves[ki]->ipoint;
      double max_dist = 0.0;
      for (int kj=0; kj<nmb_guide; ++kj)
	{
	  Point pos1 = surf1->ParamSurface::point(intcurves[ki]->epar1[2*kj],
						intcurves[ki]->epar1[2*kj+1]);
	  Point pos2 = surf2->ParamSurface::point(intcurves[ki]->epar2[2*kj],
						intcurves[ki]->epar2[2*kj+1]);
	  double dist_guide = pos1.dist(pos2);

	  // If the guide point lies at a boundary of both surfaces (only
	  // relevant for the first and the last guide point), it might be
	  // possible to improve the accuracy by iterating the point in one
	  // surface
	  if (kj==0 || kj<nmb_guide-1)
	    {
	      double u1, u2, v1, v2, d1, d2;
	      Point clo1, clo2;
	      double *seed1 = intcurves[ki]->epar1+2*kj;
	      double *seed2 = intcurves[ki]->epar2+2*kj;
	      surf1->closestBoundaryPoint(pos1, u1, v1, clo1, d1, epsgeo, NULL,
					  seed1);
	      surf2->closestBoundaryPoint(pos2, u2, v2, clo2, d2, epsgeo, NULL, 
					  seed2);
	      if (d1 < epsgeo && d2 < epsgeo)
		{
		  // Iterate with respect to the other surface
		  double u3, u4, v3, v4, d3, d4;
		  Point clo3, clo4;
		  surf1->closestPoint(pos2, u3, v3, clo3, d3, epsgeo, NULL, seed1);
		  surf2->closestPoint(pos1, u4, v4, clo4, d4, epsgeo, NULL, seed2);
		  if (d3 < dist_guide && d3 < d4)
		    {
		      intcurves[ki]->epar1[2*kj] = u3;
		      intcurves[ki]->epar1[2*kj+1] = v3;
		      dist_guide = d3;
		    }
		  else if (d4 < dist_guide)
		    {
		      intcurves[ki]->epar2[2*kj] = u4;
		      intcurves[ki]->epar2[2*kj+1] = v4;
		      dist_guide = d4;
		    }
		}
	    }
	  max_dist = std::max(max_dist, dist_guide);
	}
      if (max_dist > 0.5*march_eps)
	march_eps += 0.5*max_dist;
 
      if (intcurves[ki]->pgeom)
	{
	  // The geometry curve exists already. Use GoTools to create parameter curves
	  shared_ptr<SplineCurve> space_curve1(SISLCurve2Go(intcurves[ki]->pgeom));
	  space_curve1->makeKnotStartRegular();
	  space_curve1->makeKnotEndRegular();
	  shared_ptr<SplineCurve> space_curve2(space_curve1->clone());
	  shared_ptr<CurveOnSurface> sf_cv1, sf_cv2;
	  if (intcurves[ki]->ppar1)
	    {
	      shared_ptr<SplineCurve> pcurve1(SISLCurve2Go(intcurves[ki]->ppar1));
	      pcurve1->makeKnotStartRegular();
	      pcurve1->makeKnotEndRegular();
	      sf_cv1 = shared_ptr<CurveOnSurface>(new CurveOnSurface(sf1, pcurve1, 
								     space_curve1, false));
	    }
	  else
	    {
	      sf_cv1 = shared_ptr<CurveOnSurface>(new CurveOnSurface(sf1, space_curve1, 
								     false));
	      Point par1(intcurves[ki]->epar1[0], intcurves[ki]->epar1[1]);
	      int npar = intcurves[ki]->ipar1;
	      Point par2(intcurves[ki]->epar1[2*npar-2], intcurves[ki]->epar1[2*npar-1]);
	      sf_cv1->ensureParCrvExistence(march_eps, NULL, &par1, &par2);
	    }

	  if (intcurves[ki]->ppar2)
	    {
	      shared_ptr<SplineCurve> pcurve2(SISLCurve2Go(intcurves[ki]->ppar2));
	      pcurve2->makeKnotStartRegular();
	      pcurve2->makeKnotEndRegular();
	      sf_cv2 = shared_ptr<CurveOnSurface>(new CurveOnSurface(sf2, pcurve2, 
								     space_curve2, false));
	    }
	  else
	    {
	      sf_cv2 = shared_ptr<CurveOnSurface>(new CurveOnSurface(sf2, space_curve2, 
								     false));
	      Point par1(intcurves[ki]->epar2[0], intcurves[ki]->epar2[1]);
	      int npar = intcurves[ki]->ipoint;
	      Point par2(intcurves[ki]->epar2[2*npar-2], intcurves[ki]->epar2[2*npar-1]);
	      sf_cv2->ensureParCrvExistence(march_eps, NULL, &par1, &par2);
	    }

	  // We make sure the intersection cvs have the right direction (area below sfs
	  // to be trimmed away).
	  consistentIntersectionDir(sf_cv1, sf_cv2, epsgeo);
	  consistentIntersectionDir(sf_cv2, sf_cv1, epsgeo);

	  int_segments1.push_back(sf_cv1);
	  int_segments2.push_back(sf_cv2);
	}
      else
	{
	  s1310(sisl_sf1, sisl_sf2, intcurves[ki], march_eps, maxstep, makecurv, draw, &status);
	  SISLCurve* sc = intcurves[ki]->pgeom;
	  if (sc == 0) {
	    // Turn the sequence of guide points and try again
	    int kpoint = intcurves[ki]->ipoint;
	    int kpar1 = intcurves[ki]->ipar1;
	    int kpar2 = intcurves[ki]->ipar2;
	    for (int ka=0; ka<kpoint/2; ++ka)
	      {
	  	for (int kb=0; kb<kpar1; ++kb)
	  	  std::swap(intcurves[ki]->epar1[kpar1*ka+kb],
	  		    intcurves[ki]->epar1[kpar1*(kpoint-ka-1)+kb]);
	  	for (int kb=0; kb<kpar2; ++kb)
	  	  std::swap(intcurves[ki]->epar2[kpar2*ka+kb],
	  		    intcurves[ki]->epar2[kpar2*(kpoint-ka-1)+kb]);
	      }
	    s1310(sisl_sf1, sisl_sf2, intcurves[ki], march_eps, maxstep, makecurv, draw, &status);
	    sc = intcurves[ki]->pgeom;
	    if (sc == 0) {
	      
	      MESSAGE("s1310 returned code: " << status << ", returning.");
	      continue;
	    }

	  }

	  // Check distances between first and last guide points and the 
	  // endpoints of the generated curve
	  Point par1_1(intcurves[ki]->epar1[0],intcurves[ki]->epar1[1]);
	  Point par1_2(intcurves[ki]->epar2[0],intcurves[ki]->epar2[1]);
	  Point par2_1(intcurves[ki]->epar1[2*nmb_guide-2],
		       intcurves[ki]->epar1[2*nmb_guide-1]);
	  Point par2_2(intcurves[ki]->epar2[2*nmb_guide-2],
		       intcurves[ki]->epar2[2*nmb_guide-1]);
	  Point pos1_1 = surf1->ParamSurface::point(par1_1[0],par1_1[1]);
	  Point pos1_2 = surf2->ParamSurface::point(par1_2[0],par1_2[1]);
	  Point pos2_1 = surf1->ParamSurface::point(par2_1[0],par2_1[1]);
	  Point pos2_2 = surf2->ParamSurface::point(par2_2[0],par2_2[1]);
	  shared_ptr<SplineCurve> tmp_space(SISLCurve2Go(sc));
	  Point pos1_3 = tmp_space->ParamCurve::point(tmp_space->startparam());
	  Point pos2_3 = tmp_space->ParamCurve::point(tmp_space->endparam());
	  Point mid1 = 0.5*(pos1_1+pos1_2);
	  Point mid2 = 0.5*(pos2_1+pos2_2);
	  double dd1 = std::min(pos1_3.dist(mid1), pos1_3.dist(mid2));
	  double dd2 = std::min(pos2_3.dist(mid1), pos2_3.dist(mid2));
	  if (dd1 > epsgeo || dd2 > epsgeo)
	    {
	      // Check if the missing curve piece follows a boundary
	      // curve in one of the surfaces
	      double ptol = 1.0e-6;  // Arbitrary. Should be related to the
	      // geometry tolerance
	      RectDomain dom[2];
	      dom[0] = surf1->containingDomain();
	      dom[1] = surf2->containingDomain();

	      // Find parameter pair(s) limiting the curve piece
	      vector<pair<Point,Point> > lim_par;
	      shared_ptr<SplineCurve> tmp_par1(SISLCurve2Go(intcurves[ki]->ppar1));
	      shared_ptr<SplineCurve> tmp_par2(SISLCurve2Go(intcurves[ki]->ppar2));
	      if (dd1 > epsgeo)
		{
		  Point pp1_1 = tmp_par1->ParamCurve::point(tmp_par1->startparam());
		  Point pp1_2 = tmp_par2->ParamCurve::point(tmp_par2->startparam());
		  Point pp2_1, pp2_2;
		  if (pos1_3.dist(mid1) < pos1_3.dist(mid2))
		    {
		      pp2_1 = par1_1;
		      pp2_2 = par1_2;
		    }
		  else
		    {
		      pp2_1 = par2_1;
		      pp2_2 = par2_2;
		    }
		  lim_par.push_back(std::make_pair(pp1_1, pp2_1));
		  lim_par.push_back(std::make_pair(pp1_2, pp2_2));
		}
	      if (dd2 > epsgeo)
		{
		  Point pp1_1 = tmp_par1->ParamCurve::point(tmp_par1->endparam());
		  Point pp1_2 = tmp_par2->ParamCurve::point(tmp_par2->endparam());
		  Point pp2_1, pp2_2;
		  if (pos2_3.dist(mid1) < pos2_3.dist(mid2))
		    {
		      pp2_1 = par1_1;
		      pp2_2 = par1_2;
		    }
		  else
		    {
		      pp2_1 = par2_1;
		      pp2_2 = par2_2;
		    }
		  lim_par.push_back(std::make_pair(pp1_1, pp2_1));
		  lim_par.push_back(std::make_pair(pp1_2, pp2_2));
		}
	      for (size_t kr=0; kr<lim_par.size(); kr+=2)
		{
		  int kh;
		  for (kh=0; kh<2; ++kh)
		    {
		      SplineSurface *surf = (kh == 0) ? surf1 : surf2;
		      int bd = dom[kh].whichBoundary(Vector2D(lim_par[kr+kh].first[0],
							      lim_par[kr+kh].first[1]),
						     Vector2D(lim_par[kr+kh].second[0],
							      lim_par[kr+kh].second[1]),
						     ptol);
		      if (bd >= 0)
			{
			  // Pick associated boundary curve
			  int edg_ix[4] = {3, 1, 0, 2};
			  double epar1[3], epar2[3];
			  shared_ptr<SplineCurve> edg_cv(surf->edgeCurve(edg_ix[bd]));
			  epar1[0] = (bd <= 1) ? lim_par[kr+kh].first[1] : 
			    lim_par[kr+kh].first[0];
			  epar2[0] = (bd <= 1) ? lim_par[kr+kh].second[1] : 
			    lim_par[kr+kh].second[0];
			  epar1[1] = lim_par[kr+1-kh].first[0];
			  epar1[2] = lim_par[kr+1-kh].first[1];
			  epar2[1] = lim_par[kr+1-kh].second[0];
			  epar2[2] = lim_par[kr+1-kh].second[1];
			  
			  // Check coincidence
			  SISLCurve *qc = Curve2SISL(*edg_cv);
			  s1785(qc, (kh==0) ? sisl_sf2 : sisl_sf1, epsgeo, 
				epar1, epar2, 1, &status);
			  freeCurve(qc); // Not needed anymore
			  
			  // status < 0: error
			  // status == 0: no coincidence
			  // status == 1: coincidence
			  if (status <= 0)
			    continue;

			  double t1 = std::min(epar1[0],epar2[0]);
			  double t2 = std::max(epar1[0],epar2[0]);
			  Point par1 = (epar1[0] < epar2[0]) ? lim_par[kr+kh].first :
			    lim_par[kr+kh].second;
			  Point par2 = (epar1[0] < epar2[0]) ? lim_par[kr+kh].second :
			    lim_par[kr+kh].first;
			  shared_ptr<SplineCurve> space_cv(edg_cv->subCurve(t1,t2));
			  shared_ptr<SplineCurve> par_cv1(new SplineCurve(par1, t1,
									  par2, t2));
			  shared_ptr<CurveOnSurface> sf_cv1(new CurveOnSurface((kh==0) ? sf1 : sf2,
									       par_cv1,
									       space_cv,
									       false));
			  shared_ptr<SplineCurve> space_cv2(space_cv->clone());
			  shared_ptr<CurveOnSurface> sf_cv2(new CurveOnSurface((kh==0) ? sf2 : sf1,
									       space_cv2,
									       false));
			  sf_cv2->ensureParCrvExistence(march_eps);
			  consistentIntersectionDir(sf_cv1, sf_cv2, epsgeo);
			  consistentIntersectionDir(sf_cv2, sf_cv1, epsgeo);

			  int_segments1.push_back((kh==0) ? sf_cv1 : sf_cv2);
			  int_segments2.push_back((kh==0) ? sf_cv2 : sf_cv1);

			  break;
			}
		    }
		  if (kh < 2)
		    break;
		}
	    }

	  // ALWAYS_ERROR_IF(status < 0,
	  // 	    "Failed intersecting surfs.");
	  // MESSAGE_IF(status != 0, "Returned status value: " << status);

	  shared_ptr<SplineCurve> pcurve1(SISLCurve2Go(intcurves[ki]->ppar1));
	  shared_ptr<SplineCurve> pcurve2(SISLCurve2Go(intcurves[ki]->ppar2));
	  shared_ptr<SplineCurve> space_curve1(SISLCurve2Go(intcurves[ki]->pgeom));
	  //shared_ptr<SplineCurve> space_curve2(new SplineCurve(*space_curve1));
	  shared_ptr<SplineCurve> space_curve2(space_curve1->clone());

	  // Make sure that the curves are k-regular
	  pcurve1->makeKnotStartRegular();
	  pcurve1->makeKnotEndRegular();
	  pcurve2->makeKnotStartRegular();
	  pcurve2->makeKnotEndRegular();
	  space_curve1->makeKnotStartRegular();
	  space_curve1->makeKnotEndRegular();
	  space_curve2->makeKnotStartRegular();
	  space_curve2->makeKnotEndRegular();

	  // We make sure the intersection cvs have the right direction (area below sfs
	  // to be trimmed away).
	  consistentIntersectionDir(*pcurve1, *space_curve1, *sf1,
				    *pcurve2, *space_curve2, *sf2, epsgeo);
	  consistentIntersectionDir(*pcurve2, *space_curve2, *sf2,
				    *pcurve1, *space_curve1, *sf1, epsgeo);

	  // int_segments1.push_back(shared_ptr<CurveOnSurface>
	  // 			(new CurveOnSurface(sf1, pcurve1, space_curve1, true)));
	  // int_segments2.push_back(shared_ptr<CurveOnSurface>
	  // 			(new CurveOnSurface(sf2, pcurve2, space_curve2, true)));
	  int_segments1.push_back(shared_ptr<CurveOnSurface>
				  (new CurveOnSurface(sf1, pcurve1, space_curve1, false)));
	  int_segments2.push_back(shared_ptr<CurveOnSurface>
				  (new CurveOnSurface(sf2, pcurve2, space_curve2, false)));
	}
    }

    // As we're using SISL we must not forget to release memory.
    if (sisl_sf1)
      freeSurf(sisl_sf1);
    if (sisl_sf2)
      freeSurf(sisl_sf2);
    if (pointpar1)
      free(pointpar1);
    if (pointpar2)
      free(pointpar2);
    freeIntcrvlist(intcurves, nmb_int_cvs);
}


#ifdef __BORLANDC__
} // namespace Go
// in the ::<anonymous> namespace below, C++Builder did not "see" the
// Go::<anonymous>::comparepair_second template
#endif

namespace {
    template <typename PairType>
    struct comparepair_second
    {
	bool operator() (const PairType& p1, const PairType& p2)
	{
	    return p1.second < p2.second;
	}
    };
}

#ifdef __BORLANDC__
namespace Go
{
#endif

//===========================================================================
void BoundedUtils::translateBoundedSurf(Point trans_vec, BoundedSurface& bd_sf,
					  double deg_eps)
//===========================================================================
{
    int ki, kj;
    shared_ptr<SplineSurface> spline_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bd_sf.underlyingSurface());
    ASSERT(spline_sf.get() != 0);
    GeometryTools::translateSplineSurf(trans_vec, *spline_sf);
    vector<CurveLoop> all_bd_loops = bd_sf.allBoundaryLoops(deg_eps);
    for (ki = 0; ki < int(all_bd_loops.size()); ++ki) {
	for (kj = 0; kj < all_bd_loops[ki].size(); ++kj) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(all_bd_loops[ki][kj]);
	    ASSERT(cv_on_sf.get() != 0);
	    shared_ptr<SplineCurve> space_cv =
		dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf->spaceCurve());
	    ASSERT(space_cv.get() != 0);
	    GeometryTools::translateSplineCurve(trans_vec, *space_cv);
	}
    }
}

//===========================================================================
void BoundedUtils::rotateBoundedSurf(Point rot_axis, double alpha,
				       BoundedSurface& bd_sf, double deg_eps)
//===========================================================================
{
    int ki, kj;
    shared_ptr<SplineSurface> spline_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(bd_sf.underlyingSurface());
    ASSERT(spline_sf.get() != 0);
    GeometryTools::rotateSplineSurf(rot_axis, alpha, *spline_sf);
    vector<CurveLoop> all_bd_loops = bd_sf.allBoundaryLoops(deg_eps);
    for (ki = 0; ki < int(all_bd_loops.size()); ++ki) {
	for (kj = 0; kj < all_bd_loops[ki].size(); ++kj) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(all_bd_loops[ki][kj]);
	    ASSERT(cv_on_sf.get() != 0);
	    shared_ptr<SplineCurve> space_cv =
		dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf->spaceCurve());
	    ASSERT(space_cv.get() != 0);
	    GeometryTools::rotateSplineCurve(rot_axis, alpha, *space_cv);
	}
    }
}

//===========================================================================
void BoundedUtils::intersectWithSurfaces(vector<shared_ptr<CurveOnSurface> >& cvs1,
					   shared_ptr<BoundedSurface>& bd_sf1,
					   vector<shared_ptr<CurveOnSurface> >& cvs2,
					   shared_ptr<BoundedSurface>& bd_sf2,
					 double epsge, bool only_inner_curves)
//===========================================================================
{
    ASSERT(cvs1.size() == cvs2.size());
    int ki, kj;
    // The cvs may be defined in opposite directions.
    vector<Point> start_pt1 = cvs1[0]->ParamCurve::point(cvs1[0]->startparam(), 1);
    vector<Point> start_pt2 = cvs2[0]->ParamCurve::point(cvs2[0]->startparam(), 1);
    bool opp_dir = false;
    if ((start_pt1[0].dist(start_pt2[0]) > epsge) || // If cv is closed we must check tangent.
	((start_pt1[1].dist(-start_pt2[1])) < (start_pt1[1].dist(start_pt2[1])))) {
	opp_dir = true;
    }

    // We run through the parallell vectors extracting parts lying on sf.
    double fac = 1.5;
    double angtol = 0.1;
    double pihalf = 0.5*M_PI;
    vector<shared_ptr<CurveOnSurface> > new_cvs1, new_cvs2; //(cvs2.size());
    for (ki = 0; ki < int(cvs1.size()); ++ki) {
      // First make sure that the parameteriazation of geometry curves and
      // parameter curves correspond 
      if (!cvs1[ki]->sameOrientation())
	cvs1[ki]->enableSameOrientation();
      if (!cvs1[ki]->sameCurve(epsge))
	{
	  if (cvs1[ki]->parPref())
	    {
	      shared_ptr<ParamCurve> tmp_cv = cvs1[ki]->spaceCurve();
	      cvs1[ki]->unsetSpaceCurve();
	      cvs1[ki]->ensureSpaceCrvExistence(epsge);
	      DirectionCone cone1 = tmp_cv->directionCone();
	      DirectionCone cone2 = cvs1[ki]->spaceCurve()->directionCone();
	      double ang = cone1.centre().angle(cone2.centre());
	      if (cone2.greaterThanPi() ||
		  (cone2.angle() > fac*cone1.angle() && cone2.angle() > angtol) ||
		  cone2.angle() + 0.5*ang >= pihalf)
		{
#ifdef DEBUG
		  std::cout << "Check curve" << std::endl;
#endif
		  cvs1[ki]->setSpaceCurve(tmp_cv);
		}
	    }
	  else
	    {
	      shared_ptr<ParamCurve> tmp_cv = cvs1[ki]->parameterCurve();
	      cvs1[ki]->unsetParameterCurve();
	      bool found = cvs1[ki]->ensureParCrvExistence(epsge);
	      if (!found)
		cvs1[ki]->setParameterCurve(tmp_cv);
	      else
		{
		  DirectionCone cone1 = tmp_cv->directionCone();
		  DirectionCone cone2 = cvs1[ki]->parameterCurve()->directionCone();
		  double ang = cone1.centre().angle(cone2.centre());
		  if (cone2.greaterThanPi() ||
		      (cone2.angle() > fac*cone1.angle() && cone2.angle() > angtol) ||
		      cone2.angle() + 0.5*ang >= pihalf)
		    {
#ifdef DEBUG
		      std::cout << "Check curve" << std::endl;
#endif
		      cvs1[ki]->setParameterCurve(tmp_cv);
		    }
		}
	    }
	}

	vector<shared_ptr<CurveOnSurface> > new_cvs = 
	  intersectWithSurface(*cvs1[ki], *bd_sf1, /*0.1**/epsge,
			       only_inner_curves);
	new_cvs1.insert(new_cvs1.end(), new_cvs.begin(), new_cvs.end());
	int other_ind = opp_dir ? (int)cvs1.size() - 1 - ki : ki;

	if (!cvs2[other_ind]->sameOrientation())
	  cvs2[other_ind]->enableSameOrientation();
	if (!cvs2[other_ind]->sameCurve(epsge))
	{
	  if (cvs2[other_ind]->parPref())
	    {
	      shared_ptr<ParamCurve> tmp_cv = cvs2[other_ind]->spaceCurve();
	      cvs2[other_ind]->unsetSpaceCurve();
	      bool found = cvs2[other_ind]->ensureSpaceCrvExistence(epsge);
	      DirectionCone cone1, cone2;
	      if (found)
		{
		  cone1 = tmp_cv->directionCone();
		  cone2 = cvs2[other_ind]->spaceCurve()->directionCone();
		}
	      if (!found || cone2.greaterThanPi() ||
		  (cone2.angle() > fac*cone1.angle() && cone2.angle() > angtol))
		{
#ifdef DEBUG
		  std::cout << "Check curve" << std::endl;
#endif
		  cvs2[other_ind]->setSpaceCurve(tmp_cv);
		}
	    }
	  else
	    {
	      shared_ptr<ParamCurve> tmp_cv = cvs2[other_ind]->parameterCurve();
	      cvs2[other_ind]->unsetParameterCurve();
	      bool found = cvs2[other_ind]->ensureParCrvExistence(epsge);
	      DirectionCone cone1, cone2;
	      if (found)
		{
		  cone1 = tmp_cv->directionCone();
		  cone2 = cvs1[ki]->parameterCurve()->directionCone();
		}
	      if (!found || cone2.greaterThanPi() ||
		  (cone2.angle() > fac*cone1.angle() && cone2.angle() > angtol))
		{
#ifdef DEBUG
		  std::cout << "Check curve" << std::endl;
#endif
		  cvs2[other_ind]->setParameterCurve(tmp_cv);
		}
	    }
	}
	new_cvs = intersectWithSurface(*cvs2[other_ind], *bd_sf2, 
				       /*0.1**/epsge, only_inner_curves);
	new_cvs2.insert(new_cvs2.end(), new_cvs.begin(), new_cvs.end());
    }

    // It is more convenient to work with the cvs if they are going in the same direction.
    if (opp_dir) {
	reverse(new_cvs2.begin(), new_cvs2.end());
	for (ki = 0; ki < int(new_cvs2.size()); ++ki) {
	    new_cvs2[ki]->reverseParameterDirection();
	}
    }

    // We then make sure that the resulting cvs span the same spatial pts.
    // First split curves which end in the interior of another curve
    double epsge2 = 1.5*epsge; //1.1*epsge;  // Allow for some slack due to approximate
    double knot_diff_tol = std::min(1e-05, 0.1*epsge2);
    // intersection curves
    splitIntersectingCurves(new_cvs1, new_cvs2, epsge2, knot_diff_tol);

    // We then remove parts with a space cv lying on one sf only.
    // We extract matches to input vectors.
    cvs1.clear();
    cvs2.clear();
    for (ki = 0; ki < int(new_cvs1.size()); ++ki) {
	Point start1 = new_cvs1[ki]->ParamCurve::point(new_cvs1[ki]->startparam());
	Point end1 = new_cvs1[ki]->ParamCurve::point(new_cvs1[ki]->endparam());
	for (kj = 0; kj < int(new_cvs2.size()); ++kj) {
	    Point start2 = new_cvs2[kj]->ParamCurve::point(new_cvs2[kj]->startparam());
	    Point end2 = new_cvs2[kj]->ParamCurve::point(new_cvs2[kj]->endparam());
	    if ((start1.dist(start2) < epsge2) && (end1.dist(end2) < epsge2)) {
		cvs1.push_back(new_cvs1[ki]);
		cvs2.push_back(new_cvs2[kj]);
		break;
	    }
	}
    }

    // If second cv was reversed we must reverse it again (as it may be significant wrt trim area).
    if (opp_dir) {
	reverse(cvs2.begin(), cvs2.end());
	for (ki = 0; ki < int(cvs2.size()); ++ki) {
	    cvs2[ki]->reverseParameterDirection();
	}
    }
}

  void 
  BoundedUtils::splitIntersectingCurves(vector<shared_ptr<CurveOnSurface> >& cvs1, 
					vector<shared_ptr<CurveOnSurface> >& cvs2, 
					double epsge, double knot_diff_tol)
  {
    int ki, kj;
    for (ki = 0; ki < int(cvs1.size()); ++ki) {
      double t1 = cvs1[ki]->startparam();
      double t2 = cvs1[ki]->endparam();
	Point start1 = cvs1[ki]->ParamCurve::point(t1);
	Point end1 = cvs1[ki]->ParamCurve::point(t2);
	if (start1.dist(end1) < epsge)
	  {
	    // Check mid point
	    Point mid1 = cvs1[ki]->ParamCurve::point(0.5*(t1+t2));
	    if (mid1.dist(start1) < epsge && mid1.dist(end1) < epsge)
	      continue;
	  }
	// We split cvs which end in the interior of another cv.
	for (kj = 0; kj < int(cvs2.size()); ++kj) {
	  if (cvs1[ki].get() == cvs2[kj].get())
	    continue;  // The same curve

	  double t3 = cvs2[kj]->startparam();
	  double t4 = cvs2[kj]->endparam();
	  Point start2 = cvs2[kj]->ParamCurve::point(t3);
	  Point end2 = cvs2[kj]->ParamCurve::point(t4);
	    if (start2.dist(end2) < epsge)
	      {
		// Check mid point
		Point mid2 = cvs2[kj]->ParamCurve::point(0.5*(t3+t4));
		if (mid2.dist(start2) < epsge && mid2.dist(end2) < epsge)
		  continue;
	      }

	    // We then perform closest pt calculations for all end pts.
	    double clo_t, clo_dist;
	    Point clo_pt;
	    double len1 = cvs1[ki]->estimatedCurveLength();
	    cvs1[ki]->ParamCurve::closestPoint(end2, clo_t, clo_pt, clo_dist);
	    double len1_1 = 
	      cvs1[ki]->estimatedCurveLength(cvs1[ki]->startparam(), clo_t);
	    double len1_2 = 
	      cvs1[ki]->estimatedCurveLength(clo_t, cvs1[ki]->endparam());
	    if ((clo_dist < epsge && len1 > epsge && len1_1 > epsge && len1_2 > epsge) && 
		((clo_t - knot_diff_tol > cvs1[ki]->startparam()) &&
		 (clo_t + knot_diff_tol < cvs1[ki]->endparam()))) 
	      {
		shared_ptr<CurveOnSurface> new_cv(cvs1[ki]->subCurve(clo_t, cvs1[ki]->endparam()));
		cvs1.insert(cvs1.begin() + ki + 1, new_cv);
		cvs1[ki] = shared_ptr<CurveOnSurface>
		    (cvs1[ki]->subCurve(cvs1[ki]->startparam(), clo_t));
	    }
	    cvs1[ki]->ParamCurve::closestPoint(start2, clo_t, clo_pt, clo_dist);
	    double len1_3 = 
	      cvs1[ki]->estimatedCurveLength(cvs1[ki]->startparam(), clo_t);
	    double len1_4 = 
	      cvs1[ki]->estimatedCurveLength(clo_t, cvs1[ki]->endparam());	    
	    if ((clo_dist < epsge && len1 > epsge &&  
		 len1_3 > epsge && len1_4 > epsge) && 
		((clo_t - knot_diff_tol > cvs1[ki]->startparam()) &&
		 (clo_t + knot_diff_tol < cvs1[ki]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(cvs1[ki]->subCurve(clo_t, cvs1[ki]->endparam()));
		cvs1.insert(cvs1.begin() + ki + 1, new_cv);
		cvs1[ki] = shared_ptr<CurveOnSurface>
		    (cvs1[ki]->subCurve(cvs1[ki]->startparam(), clo_t));
	    }
	    double len2 = cvs2[kj]->estimatedCurveLength();
	    if (len2 < epsge)
	      continue;
	    cvs2[kj]->ParamCurve::closestPoint(end1, clo_t, clo_pt, clo_dist);
	    double len2_1 = 
	      cvs2[kj]->estimatedCurveLength(cvs2[kj]->startparam(), clo_t);
	    double len2_2 = 
	      cvs2[kj]->estimatedCurveLength(clo_t, cvs2[kj]->endparam());
	    if ((clo_dist < epsge && len2_1 > epsge && len2_2 > epsge) && 
		((clo_t - knot_diff_tol > cvs2[kj]->startparam()) &&
				       (clo_t + knot_diff_tol < cvs2[kj]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(cvs2[kj]->subCurve(clo_t, cvs2[kj]->endparam()));
		cvs2.insert(cvs2.begin() + kj + 1, new_cv);
		cvs2[kj] = shared_ptr<CurveOnSurface>
		    (cvs2[kj]->subCurve(cvs2[kj]->startparam(), clo_t));
	    }
	    cvs2[kj]->ParamCurve::closestPoint(start1, clo_t, clo_pt, clo_dist);
	    double len2_3 = 
	      cvs2[kj]->estimatedCurveLength(cvs2[kj]->startparam(), clo_t);
	    double len2_4 = 
	      cvs2[kj]->estimatedCurveLength(clo_t, cvs2[kj]->endparam());

	    if ((clo_dist < epsge && len2_3 > epsge && len2_4 > epsge) && 
		((clo_t - knot_diff_tol > cvs2[kj]->startparam()) &&
		 (clo_t + knot_diff_tol < cvs2[kj]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(cvs2[kj]->subCurve(clo_t, cvs2[kj]->endparam()));
		cvs2.insert(cvs2.begin() + kj + 1, new_cv);
		cvs2[kj] = shared_ptr<CurveOnSurface>
		    (cvs2[kj]->subCurve(cvs2[kj]->startparam(), clo_t));
	    }
	}
    }
  }

vector<vector<shared_ptr<BoundedSurface> > >
BoundedUtils::trimSurfsWithSurfs(const vector<shared_ptr<ParamSurface> >& sfs1,
				   const vector<shared_ptr<ParamSurface> >& sfs2, double epsge)
{
    int ki, kj;
    int nmb1 = (int)sfs1.size();
    int nmb2 = (int)sfs2.size();
 
    vector<vector<shared_ptr<BoundedSurface> > > trimmed_sfs(2);
    // Transform the two sets of surfaces to bounded surfaces.
    vector<shared_ptr<BoundedSurface> > bounded_sfs1, bounded_sfs2;
    vector<shared_ptr<ParamSurface> > under_sfs1, under_sfs2;
    for (ki = 0; ki < int(sfs1.size()); ++ki)
	if (sfs1[ki]->instanceType() == Class_BoundedSurface) {
	    shared_ptr<BoundedSurface> bounded_sf =
		dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs1[ki]);
	    bounded_sfs1.push_back(bounded_sf);
	    shared_ptr<SplineSurface> under_sf =
		dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf->underlyingSurface());
	    ALWAYS_ERROR_IF(under_sf.get() == 0,
			"Expecting underlying surface to be a spline surface.");
	    under_sfs1.push_back(under_sf);
	} else //if (sfs1[ki]->instanceType() == Class_SplineSurface)
	    try {
	      vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(sfs1[ki],
								     epsge);
	      shared_ptr<BoundedSurface> bounded_sf = 
		shared_ptr<BoundedSurface>(new BoundedSurface(sfs1[ki], loops));
		bounded_sfs1.push_back(bounded_sf);
	    shared_ptr<ParamSurface> under_sf = bounded_sf->underlyingSurface();
	    under_sfs1.push_back(under_sf);
	    } catch (...) {
		THROW("Something went wrong, returning.");
	    }
// 	else
// 	    THROW("Surface type not supported.");
    for (ki = 0; ki < int(sfs2.size()); ++ki)
	if (sfs2[ki]->instanceType() == Class_BoundedSurface) {
	    shared_ptr<BoundedSurface> bounded_sf =
		dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs2[ki]);
	    bounded_sfs2.push_back(bounded_sf);
	    shared_ptr<SplineSurface> under_sf =
		dynamic_pointer_cast<SplineSurface, ParamSurface>(bounded_sf->underlyingSurface());
	    ALWAYS_ERROR_IF(under_sf.get() == 0,
			"Expecting underlying surface to be a spline surface.");
	    under_sfs2.push_back(under_sf);
	} else //if (sfs2[ki]->instanceType() == Class_SplineSurface)
	    try {
	      vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(sfs2[ki],
								     epsge);
	      shared_ptr<BoundedSurface> bounded_sf = 
		shared_ptr<BoundedSurface>(new BoundedSurface(sfs2[ki], loops));
		bounded_sfs2.push_back(bounded_sf);
	    shared_ptr<ParamSurface> under_sf = bounded_sf->underlyingSurface();
	    under_sfs2.push_back(under_sf);
	    } catch (...) {
		THROW("Something went wrong, returning.");
	    }
// 	else
// 	    THROW("Surface type not supported.");


    // We next perform pairwise intersection against members of the two sets.
    vector<vector<shared_ptr<CurveOnSurface> > > all_int_segments1(nmb1), all_int_segments2(nmb2);
    for (ki = 0; ki < int(under_sfs1.size()); ++ki)
	for (kj = 0; kj < int(under_sfs2.size()); ++kj) {
	    vector<shared_ptr<CurveOnSurface> > int_segments1, int_segments2;
	    double tol = 1e-04*epsge;
	    try {
	      shared_ptr<ParamSurface> under1 = under_sfs1[ki];
	      shared_ptr<ParamSurface> under2 = under_sfs2[kj];
#ifdef DEBUG1
	      std::ofstream debug0("int_crv_surf0.g2");
	      under1->writeStandardHeader(debug0);
	      under1->write(debug0);
	      under2->writeStandardHeader(debug0);
	      under2->write(debug0);
#endif
	      getIntersectionCurve(under1, under2, int_segments1, int_segments2, tol);
	    } catch (...) {
		THROW("Failed intersecting the two spline surfaces.");
	    }
	    if (int_segments1.size() != 0) {
		// We must then extract parts lying on both the actual sfs.
		intersectWithSurfaces(int_segments1, bounded_sfs1[ki],
				      int_segments2, bounded_sfs2[kj], 0.1*epsge);

		all_int_segments1[ki].insert(all_int_segments1[ki].end(),
					     int_segments1.begin(), int_segments1.end());
		all_int_segments2[kj].insert(all_int_segments2[kj].end(),
					     int_segments2.begin(), int_segments2.end());
	    }
	}

    // Finally we extract BoundedSurface's from segments.
    for (ki = 0; ki < int(all_int_segments1.size()); ++ki) {
	vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	try {
	    loop_curves =
	      getBoundaryLoops(*bounded_sfs1[ki], all_int_segments1[ki], epsge);
	} catch (...) {
	    MESSAGE("Failed extracting boundary loop.");
	}
	vector<shared_ptr<BoundedSurface> > trim_sfs =
	   createTrimmedSurfs(loop_curves, under_sfs1[ki], epsge);
	trimmed_sfs[0].insert(trimmed_sfs[0].end(), trim_sfs.begin(), trim_sfs.end());
    }
    for (ki = 0; ki < int(all_int_segments2.size()); ++ki) {
	vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
	try {
	    loop_curves =
	      getBoundaryLoops(*bounded_sfs2[ki], all_int_segments2[ki], epsge);
	} catch (...) {
	    MESSAGE("Failed extracting boundary loop.");
	}
	vector<shared_ptr<BoundedSurface> > trim_sfs =
	   createTrimmedSurfs(loop_curves, under_sfs2[ki], epsge);
	trimmed_sfs[1].insert(trimmed_sfs[1].end(), trim_sfs.begin(), trim_sfs.end());

    }

    return trimmed_sfs;
}

  bool compare_ints(std::pair<double,Point> p1, std::pair<double,Point> p2)
  {
    return (p1.first < p2.first);
  }

//===========================================================================
vector<shared_ptr<ParamCurve> > 
BoundedUtils::trimCurveWithSurf(const shared_ptr<ParamCurve>& curve,
				const shared_ptr<ParamSurface>& surf,
				Point midpoint, double epsge, double min_len)
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > result;

  // Compute intersections
  vector<pair<double,Point> > int_pts;
  vector<int> pretop;
  vector<pair<pair<double,Point>, pair<double,Point> > > int_cvs;
  intersectParamCurveSurf(curve.get(), surf.get(), epsge, int_pts,
			  pretop, int_cvs);

  // Extend endpoints of curve pieces with intersection curves
  for (size_t ki=0; ki<int_cvs.size(); ++ki)
    {
      int_pts.push_back(int_cvs[ki].first);
      int_pts.push_back(int_cvs[ki].second);
    }

  if (int_pts.size() == 0)
    {
      result.push_back(curve);
      return result;
    }

  // Check if the inside test must be negated
  bool negate = false;
  if (midpoint.dimension() > 0)
    {
      double upar, vpar, dist;
      Point clo_pt;
      surf->closestPoint(midpoint, upar, vpar, clo_pt, dist, epsge);
      Point norm0;
      surf->normal(norm0, upar, vpar);
      Point vec0 = midpoint - clo_pt;
      if (norm0*vec0 > 0.0)
	negate = true;
    }

  // Sort intersection points
  std::sort(int_pts.begin(), int_pts.end(), compare_ints);
  double t1 = curve->startparam();
  double t2 = curve->endparam();
  double del = 0.01*(t2 - t1);
  
  // Check start segment
  if (curve->estimatedCurveLength(t1, int_pts[0].first) > min_len)
    {
      double par = std::max(0.5*(t1+int_pts[0].first), 
			       int_pts[0].first-del);
      Point pos = curve->point(par);
      Point int_pos = curve->point(int_pts[0].first);
      Point normal, sf_pos;
      surf->normal(normal, int_pts[0].second[0], int_pts[0].second[1]);
      Point vec = pos - int_pos;
      if (negate)
	vec *= -1;
      if (vec*normal < 0.0)
	{
	  shared_ptr<ParamCurve> sub(curve->subCurve(t1, int_pts[0].first));
	  result.push_back(sub);
	}
    }

  // Internal segments
  double a_tol = 1.0e-12;
  for (size_t ki=1; ki<int_pts.size(); ++ki)
    {
      if (curve->estimatedCurveLength(int_pts[ki-1].first, 
				      int_pts[ki].first) > min_len)
	{
	  // Check first coincidence intervals
	  size_t kj;
	  for (kj=0; kj<int_cvs.size(); ++kj)
	    {
	      if (fabs(int_pts[ki-1].first - int_cvs[kj].first.first) < a_tol &&
		  fabs(int_pts[ki].first - int_cvs[kj].second.first) < a_tol)
		break;
	    }
	  if (kj < int_cvs.size())
	    continue;  // Coincidence interval

	  double par = std::max(0.5*(int_pts[ki-1].first+int_pts[ki].first), 
				int_pts[ki].first-del);
	  Point pos = curve->point(par);
	  Point int_pos = curve->point(int_pts[ki].first);
	  Point normal, sf_pos;
	  surf->normal(normal, int_pts[ki].second[0], int_pts[ki].second[1]);
	  Point vec = pos - int_pos;
	  if (negate)
	    vec *= -1;
	  if (vec*normal < 0.0)
	    {
	      shared_ptr<ParamCurve> sub(curve->subCurve(int_pts[ki-1].first,
							 int_pts[ki].first));
	      result.push_back(sub);
	    }
	}
    }
  
  // Last segment
  size_t ix = int_pts.size()-1;
  if (curve->estimatedCurveLength(int_pts[ix].first, t2) > min_len)
    {
      double par = std::min(0.5*(int_pts[ix].first+t2), 
			    int_pts[ix].first+del);
      Point pos = curve->point(par);
      Point int_pos = curve->point(int_pts[ix].first);
      Point normal, sf_pos;
      surf->normal(normal, int_pts[0].second[ix], int_pts[ix].second[1]);
      Point vec = pos - int_pos;
      if (negate)
	vec *= -1;
      if (vec*normal < 0.0)
	{
	  shared_ptr<ParamCurve> sub(curve->subCurve(int_pts[ix].first, t2));
	  result.push_back(sub);
	}
    }
  
  // Check coincidence with larger tolerance
  if (min_len > epsge)
    {
    }

  return result;
}

//===========================================================================
vector<shared_ptr<BoundedSurface> > 
BoundedUtils::trimSurfWithSurfs(shared_ptr<ParamSurface>& sf,
				const vector<shared_ptr<ParamSurface> >& sfs2, 
				double epsge, double minlen)
//===========================================================================
{
  size_t ki, kj;
    int nmb2 = (int)sfs2.size();

    if (minlen <= 0.0)
      minlen = epsge;
 
    // Represent the surfaces as bounded
    shared_ptr<BoundedSurface> bd_sf = convertToBoundedSurface(sf, epsge);
    shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
    vector<shared_ptr<BoundedSurface> > bd_sfs2(nmb2);
    vector<shared_ptr<ParamSurface> > under_sfs2(nmb2);
    for (kj = 0; kj < sfs2.size(); ++kj) 
      {
	bd_sfs2[kj] = convertToBoundedSurface(sfs2[kj], epsge);
	under_sfs2[kj] = bd_sfs2[kj]->underlyingSurface();
      }
    

    // We intersect the given surface with all members of the other set
    vector<shared_ptr<CurveOnSurface> >  all_int_segments;
    for (kj = 0; kj < sfs2.size(); ++kj) 
      {
#ifdef DEBUG1
	std::ofstream debug0("int_crv_surf0.g2");
	under_sf->writeStandardHeader(debug0);
	under_sf->write(debug0);
	under_sfs2[kj]->writeStandardHeader(debug0);
	under_sfs2[kj]->write(debug0);
#endif
	vector<shared_ptr<CurveOnSurface> > int_segments1, int_segments2;
	try {
	  getIntersectionCurve(under_sf, under_sfs2[kj], 
			       int_segments1, int_segments2, epsge);
	} catch (...) {
	  THROW("Failed intersecting the two spline surfaces.");
	}
	if (int_segments1.size() != 0) 
	  {
	  // We must then extract parts lying on both the actual sfs.
	  intersectWithSurfaces(int_segments1, bd_sf,
				int_segments2, bd_sfs2[kj], epsge);

	  all_int_segments.insert(all_int_segments.end(),
				  int_segments1.begin(), int_segments1.end());
	  }
      }

#ifdef DEBUG1
    std::ofstream out1("int_seg1.g2");
    for (ki=0; ki<all_int_segments.size(); ++ki)
      {
	out1 << "100 1 0 4 255 0 0 255" << std::endl;
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(all_int_segments[ki]->geometryCurve());
	tmp1->write(out1);
      }
#endif

    // Split intersection curves where they intersect (in case the surfaces
    // in the set is a possibly intersecting collection)
    vector<vector<double> > par(all_int_segments.size());
    for (ki=0; ki<all_int_segments.size(); ++ki)
      {
	for (kj=ki+1; kj<all_int_segments.size(); ++kj)
	  {
	    vector<pair<double,double> > int_pars;
	    bool use_par_crvs = (all_int_segments[ki]->hasParameterCurve() &&
				 all_int_segments[kj]->hasParameterCurve());
	    if (use_par_crvs)
	      intersectParamCurves(all_int_segments[ki]->parameterCurve().get(), 
				   all_int_segments[kj]->parameterCurve().get(), 
				   epsge, int_pars);
	    else
	      intersectParamCurves(all_int_segments[ki].get(), 
				   all_int_segments[kj].get(), epsge, int_pars);
	    if (int_pars.size() > 0)
	      {
		// Split intersection curves
		// First represent in independent vectors and add endpoints
		for (size_t kr=0; kr<int_pars.size(); ++kr)
		  {
		    if (use_par_crvs)
		      {
			// Iterate to a better position
			double par1, par2, dist;
			Point ptc1, ptc2;
			double start1 = all_int_segments[ki]->startparam();
			double end1 = all_int_segments[ki]->endparam();
			double start2 = all_int_segments[kj]->startparam();
			double end2 = all_int_segments[kj]->endparam();
			ClosestPoint::closestPtCurves(all_int_segments[ki].get(),
						      all_int_segments[kj].get(),
						      start1, end1,
						      start2, end2,
						      int_pars[kr].first, 
						      int_pars[kr].second,
						      par1, par2, dist, ptc1, ptc2);
			par[ki].push_back(par1);
			par[kj].push_back(par2);
		      }
		    else
		      {
			par[ki].push_back(int_pars[kr].first);
			par[kj].push_back(int_pars[kr].second);
		      }
		  }
	      }
	  }
      }

     for (ki=0; ki<par.size();)
      {
	if (par[ki].size() == 0)
	  {
	    ++ki;
	    continue;
	  }

	std::sort(par[ki].begin(), par[ki].end());

	// Extend with endpoints if necessary
	if (fabs(par[ki][0]-all_int_segments[ki]->startparam()) < epsge)
	  par[ki][0] = all_int_segments[ki]->startparam();
	else
	  par[ki].insert(par[ki].begin(), all_int_segments[ki]->startparam());
	if (fabs(all_int_segments[ki]->endparam()-par[ki][par[ki].size()-1]) < epsge)
	  par[ki][par[ki].size()-1] = all_int_segments[ki]->endparam();
	else
	  par[ki].push_back(all_int_segments[ki]->endparam());

	// Split intersection curves in found parameters
	vector<shared_ptr<CurveOnSurface> > sub_int;
	for (size_t kr=1; kr<par[ki].size(); ++kr)
	  {
	    if (par[ki][kr]-par[ki][kr-1] > epsge)
	      {
		shared_ptr<CurveOnSurface> sub(all_int_segments[ki]->subCurve(par[ki][kr-1], par[ki][kr]));
		sub_int.push_back(sub);
	      }
	  }
	if (sub_int.size() > 0)
	  {
	    all_int_segments.erase(all_int_segments.begin()+ki);
	    all_int_segments.insert(all_int_segments.end(), 
				    sub_int.begin(), sub_int.end());
	    par.erase(par.begin()+ki);
	  }
	else
	  ki++;
      }

#ifdef DEBUG1
    std::ofstream out2("int_seg2.g2");
    for (ki=0; ki<all_int_segments.size(); ++ki)
      {
	out2 << "100 1 0 4 0 255 0 255" << std::endl;
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(all_int_segments[ki]->geometryCurve());
	tmp1->write(out2);
      }
#endif

    // Remove intersection curves with loose ends
    // First collect the boundary curves of the initial surface
    vector<CurveLoop> boundary_loops = bd_sf->absolutelyAllBoundaryLoops();
    vector<shared_ptr<ParamCurve> > bd_cvs;
    for (size_t kr=0; kr<boundary_loops.size(); ++kr)
      {
	vector<shared_ptr<ParamCurve> > curr_cvs = boundary_loops[kr].getCurves();
	bd_cvs.insert(bd_cvs.end(), curr_cvs.begin(), curr_cvs.end());
      }

    // Check endpoints of intersection curves
    for (ki=0; ki<all_int_segments.size(); )
      {
	// First check with itself
	Point pos[2];
	all_int_segments[ki]->point(pos[0], all_int_segments[ki]->startparam());
	all_int_segments[ki]->point(pos[1], all_int_segments[ki]->endparam());
	if (pos[0].dist(pos[1]) < minlen)
	  {
	    ++ki;
	    continue;
	  }

	int ka;
	for (ka=0; ka<2; ++ka)
	  {
	    // First check if the endpoint coincides with the endpoint of 
	    // another intersection curve
	    size_t kj;
	    for (kj=0; kj<all_int_segments.size(); kj++)
	      {
		if (ki == kj)
		  continue;
		double endpar[2];
		endpar[0] = all_int_segments[kj]->startparam();	
		endpar[1] = all_int_segments[kj]->endparam();
		int kb;
		for (kb=0; kb<2; ++kb)
		  {
		    Point pos2 = all_int_segments[kj]->ParamCurve::point(endpar[kb]);
		    if (pos[ka].dist(pos2) < minlen)
		      break;
		  }
		if (kb < 2)
		  break;
	      }
	    if (kj < all_int_segments.size())
	      continue;

	    // Check if the current intersection curve ends up in a 
	    // boundary curve
	    size_t kr;
	    for (kr=0; kr<bd_cvs.size(); ++kr)
	      {
		double cpar, dist;
		Point ptc;
		bd_cvs[kr]->closestPoint(pos[ka], bd_cvs[kr]->startparam(),
					 bd_cvs[kr]->endparam(), cpar, ptc,
					     dist);
		if (dist < minlen)
		  break;
	      }
	    if (kr == bd_cvs.size())
	      break;
	  }

	if (ka < 2)
	  {
	    // Loose end found, remove curve
	    all_int_segments.erase(all_int_segments.begin()+ki);
	  }
	else
	  ++ki;
      }
    
 #ifdef DEBUG1
    std::ofstream out3("int_seg3.g2");
    for (ki=0; ki<all_int_segments.size(); ++ki)
      {
	out3 << "100 1 0 4 0 255 0 255" << std::endl;
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(all_int_segments[ki]->geometryCurve());
	tmp1->write(out3);
      }
#endif

   // Extend intersection curve pool with oppositely oriented curves to
    // ensure that a closed loop can be extracted
    size_t nmb_seg = all_int_segments.size();
    for (kj=0; kj<nmb_seg; ++kj)
      {
	shared_ptr<CurveOnSurface> cv(all_int_segments[kj]->clone());
	cv->reverseParameterDirection();
	all_int_segments.push_back(cv);
      }

    // Finally we extract BoundedSurface's from segments.
    vector<vector<shared_ptr<CurveOnSurface> > > loop_curves;
    try {
      loop_curves =
	getBoundaryLoops(*bd_sf, all_int_segments, epsge);
    } catch (...) {
      MESSAGE("Failed extracting boundary loop.");
    }
    vector<shared_ptr<BoundedSurface> > trim_sfs =
      createTrimmedSurfs(loop_curves, under_sf, epsge);
    return trim_sfs;
}


//===========================================================================
void BoundedUtils::trimSurfaceKinks(const BoundedSurface& sf, double max_normal_angle,
				     std::vector<double>& g1_disc_u, 
				     std::vector<double>& g1_disc_v,
				     bool compute_g1_disc)
//===========================================================================
{
    // First compute kinks in the underlying surface
    const BoundedSurface *parsf = &sf;
    shared_ptr<const SplineSurface> surf;
    shared_ptr<const BoundedSurface> bd;
    while (parsf)
    {
	shared_ptr<const ParamSurface> under_sf = parsf->underlyingSurface();
	surf = dynamic_pointer_cast<const SplineSurface, const ParamSurface>(under_sf);
	parsf = 0;
	if (!surf.get())
	{
	    bd = dynamic_pointer_cast<const BoundedSurface, const ParamSurface>(under_sf);
	    parsf = bd.get();
	}
    }
    if (!surf.get())
    {
	// The underlying surface is not a spline. Nothing to do
	return;
    }
    
    // Compute kinks
    GeometryTools::surfaceKinks(*surf, max_normal_angle, g1_disc_u, g1_disc_v, compute_g1_disc);

    // Get parameter domain
    const CurveBoundedDomain domain = sf.parameterDomain();

    // For all kinks, check if the corresponding constant parameter curve lies in the
    // trimmed surface
    double int_tol = 1.0e-6;
    double u_start = surf->startparam_u();
    double u_end = surf->endparam_u();
    double v_start = surf->startparam_v();
    double v_end = surf->endparam_v();
    int ki;
    for (ki=0; ki<(int)g1_disc_u.size(); ki++)
    {
	Point pnt1(g1_disc_u[ki], v_start), pnt2(g1_disc_u[ki], v_end);
	SplineCurve parcrv(pnt1, pnt2);

	vector<double> trim_params;
	domain.findPcurveInsideSegments(parcrv, int_tol, trim_params);

	if (trim_params.size() == 0 || trim_params.size() == 1)
	{
	  // Not intersecting or non-recognized tangential intersection. 
	  // Remove discontinuity parameter
	  g1_disc_u.erase(g1_disc_u.begin()+ki);
	  ki--;
	}
    }

    for (ki=0; ki<(int)g1_disc_v.size(); ki++)
    {
	Point pnt1(u_start, g1_disc_v[ki]), pnt2(u_end, g1_disc_v[ki]);
	SplineCurve parcrv(pnt1, pnt2);

	vector<double> trim_params;
	domain.findPcurveInsideSegments(parcrv, int_tol, trim_params);
	if (trim_params.size() == 0 || trim_params.size() == 1)
	{
	  // Not intersecting or non-recognized tangential intersection. 
	  //Remove discontinuity parameter
	  g1_disc_v.erase(g1_disc_v.begin()+ki);
	  ki--;
	}
    }

}	
	
//===========================================================================
int
BoundedUtils::checkAndFixLoopOrientation(shared_ptr<BoundedSurface> surf)
//===========================================================================
{
    // First check if this computation is performed already
  int fix = 0;
    if (surf->orientationIsSet())
      {
	fix = surf->orientationOK();
	return (fix) ? 1 : 0;
      }

    int nmb_loops = surf->numberOfLoops();
    for (int ki=0; ki<nmb_loops; ++ki)
    {
	// Get curves belonging to this loop
	vector<shared_ptr<ParamCurve> > curves = surf->loop(ki)->getCurves();

	vector<shared_ptr<CurveOnSurface> > sf_curves;
	LoopUtils::representAsSurfaceCurves(curves, surf, sf_curves);

	// Check orientation of current loop
	double tol = surf->loop(ki)->getSpaceEpsilon();  // Tolerance used in intersection
	bool is_CCW;
	try {
	  is_CCW = LoopUtils::paramIsCCW(sf_curves, tol, tol);
	}
	catch (...)
	  {
	    MESSAGE("Non-existing parameter curves. Cannot check orientation of loop");
	    return 0;
	  }

	if (is_CCW && ki > 0)
	{
	    // Inner trimming curve. Counter clockwise orientation.  Turn orientation
	    surf->turnLoopOrientation(ki);
	    fix = 2;
	}
	if (!is_CCW && ki == 0)
	{
	    // Outer trimming curve. Not counter clockwise orientation.  Turn orientation
	    surf->turnLoopOrientation(ki);
	    fix = 2;
	}
    }
    if (!fix)
	surf->setOrientationOK();

    return fix;
}




//-----------------------------------------------------------------------------
shared_ptr<SplineSurface>
BoundedUtils::makeTrimmedPlane(shared_ptr<Plane>& plane, 
			       vector<shared_ptr<ParamCurve> >& space_crvs)
//-----------------------------------------------------------------------------
{
    // Make rotated bounding box of trimming curves
    vector<Point> axis(3);
    plane->getSpanningVectors(axis[0],axis[1]);  // Two vectors
						 // spanning the plane
    BoundingBox box = space_crvs[0]->boundingBox();
    for (int i=1; i<(int)space_crvs.size(); ++i) {
	box.addUnionWith(space_crvs[i]->boundingBox());
    }

    // Get corners of bounding box 
    Point low = box.low();
    Point high = box.high();

    // Extending the box with 10 % in each dir (to handle a flat box).
    // Extend the trimmed plane with 10% in each dir.
    Point vec = high - low;
    double zero_tol = 1e-06;
    int ki;
    for (ki = 0; ki < 3; ++ki) {
      if (vec[ki] < zero_tol) {
	// We then add a randomly chosen value.
	double random_value = 1.0;
	low[ki] -= random_value;
	high[ki] += random_value;
      } else {
	low[ki] -= 0.1*vec[ki];
	high[ki] += 0.1*vec[ki];
      }
    }

//     // The box is set from the space curves. For a spline
//     // curve that is the coefficients.
//     Point rot_low, rot_high;
//     if (space_crvs[0]->instanceType() == Class_SplineCurve) {
// 	axis[2] = axis[0]%axis[1];
// 	vector<double> all_cv_coefs;
// 	for (size_t ki = 0; ki < space_crvs.size(); ++ki) {
// 	    shared_ptr<SplineCurve> spline_cv =
// 		dynamic_pointer_cast<SplineCurve, ParamCurve>(space_crvs[ki]);
// 	    all_cv_coefs.insert(all_cv_coefs.end(),
// 				spline_cv->coefs_begin(),
// 				spline_cv->coefs_end());
// 	}
// 	int dim = space_crvs[0]->dimension();
// 	int nmb_coefs = all_cv_coefs.size()/dim;
// 	RotatedBox rot_box(all_cv_coefs.begin(), dim,
// 			   nmb_coefs, 1, &axis[0]);
// 	rot_low = rot_box.low();
// 	rot_high = rot_box.high();
// 	Point proj_rot_low = plane->projectPoint(rot_low);
// 	Point proj_rot_high = plane->projectPoint(rot_high);
// 	double debug_val = 0.0;
//     }

    Point normal = plane->getNormal();

    MatrixXD<double, 3> cs;
    Vector3D from(normal[0],normal[1],normal[2]);
    Vector3D to(0.0, 0.0, 1.0);
    cs.setToRotation(to, from);
//     low = cs*low;
//     high = cs*high;
//     Point y2(low[0],high[1],low[2]); // Box corner.
//     Point x2(high[0],low[1],low[2]); // Box corner.

    MatrixXD<double, 3> cs1, cs2; // Rotation for the two axis directions.
    Vector3D from1(axis[0][0],axis[0][1],axis[0][2]);
    Vector3D to1(1.0, 0.0, 0.0);
    Vector3D from2(axis[1][0],axis[1][1],axis[1][2]);
    Vector3D to2(0.0, 1.0, 0.0);
    cs1.setToRotation(to1, from1);
    cs2.setToRotation(to2, from2);
    low = cs1*low;
    low = cs2*low;
    high = cs1*high;
    high = cs2*high;
    Point y2(low[0],high[1],low[2]); // Box corner.
    Point x2(high[0],low[1],low[2]); // Box corner.

    // Rotate back to xyz coordinate system
    MatrixXD<double, 3> cs3;
    cs3.setToRotation(from, to);
//     low = cs3*low;
//     high = cs3*high;
//     y2 = cs3*y2;
//     x2 = cs3*x2;

    MatrixXD<double, 3> cs1b, cs2b;
    cs1b.setToRotation(from1, to1);
    cs2b.setToRotation(from2, to2);
    low = cs2b*low;
    low = cs1b*low;
    high = cs2b*high;
    high = cs1b*high;
    y2 = cs2b*y2;
    y2 = cs1b*y2;
    x2 = cs2b*x2;
    x2 = cs1b*x2;
    
    // Project the box corners onto the plane
    low = plane->projectPoint(low);
    high = plane->projectPoint(high);
    y2 = plane->projectPoint(y2);
    x2 = plane->projectPoint(x2);

//     // @@@ VSK 0209. This is a hack. It seems that either some iges
//     // files have a misplaced plane or the plane definition is
//     // buggy. Move the plane position to fit with the trimming curve
//     int ki;
//     double dist;
//     double maxdist = -1;
//     double mindist = 1.0e8;
//     double frac = 0.001;
//     for (ki=0; ki<(int)space_crvs.size(); ++ki) 
//     {
// 	shared_ptr<SplineCurve> c1 =
// 	    (space_crvs[ki]->instanceType() == Class_SplineCurve) ?
// 	    dynamic_pointer_cast<SplineCurve>(space_crvs[ki]) :
// 	    shared_ptr<SplineCurve>(space_crvs[ki]->geometryCurve());
// 	if (c1)
// 	{
// 	    vector<double>::iterator start = c1->coefs_begin();
// 	    vector<double>::iterator end = c1->coefs_end();
// 	    for (; start != end; start+=3)
// 	    {
// 		dist = plane->distance(Point(start, start+3));
// 		maxdist = std::max(maxdist, dist);
// 		mindist = std::min(mindist, dist);
// 	    }
// 	}
//     }
//     if (maxdist >= 0.0 && maxdist - mindist < frac*maxdist)
//     {
// 	// Move trimmed plane
// 	Point pnt = space_crvs[0]->point(space_crvs[0]->startparam()); 
// 	Point vec = pnt - plane->projectPoint(pnt);
// 	low += vec;
// 	high += vec;
// 	x2 += vec;
// 	y2 += vec;
//     }

    // We project low & high onto plane to get parameter values, as we
    // want to keep parametrization of plane.
    double umin, umax, vmin, vmax;
    double clo_dist_low, clo_dist_high;
    double epsgeo = 1e-06;
    Point clo_low, clo_high;
    plane->closestPoint(low, umin, vmin, clo_low, clo_dist_low, epsgeo);
    plane->closestPoint(high, umax, vmax, clo_high, clo_dist_high, epsgeo);
    if (umax < umin) {
      swap(low, x2);
      swap(y2, high);
      swap(umin, umax);
    }
    if (vmax < vmin) {
      swap(low, y2);
      swap(x2, high);
      swap(vmin, vmax);
    }
	
    // Make bilinear spline surface
//     double knots[4];
//     knots[0] = knots[1] = 0.0;
//     knots[2] = knots[3] = 1.0;
    double knots_u[4];
    knots_u[0] = knots_u[1] = umin;
    knots_u[2] = knots_u[3] = umax;
    double knots_v[4];
    knots_v[0] = knots_v[1] = vmin;
    knots_v[2] = knots_v[3] = vmax;
    double coefs[12];
    for (ki=0; ki<3; ++ki)
	coefs[ki] = low[ki];
    for (ki=0; ki<3; ++ki)
	coefs[3+ki] = x2[ki];
    for (ki=0; ki<3; ++ki)
	coefs[6+ki] = y2[ki];
    for (ki=0; ki<3; ++ki)
	coefs[9+ki] = high[ki];
    shared_ptr<SplineSurface> surf = 
	shared_ptr<SplineSurface>
      (new SplineSurface(2, 2, 2, 2, knots_u, knots_v, coefs, 3));
	
		
    return surf;
}


//===========================================================================
void BoundedUtils::translatePlaneToCurves(shared_ptr<Go::Plane>& plane,
					  std::vector<shared_ptr
					  <Go::ParamCurve> >&
					  space_crvs)
//===========================================================================
{

#if 0
    // Get corners of bounding box 
    Point low = box.low();
    Point high = box.high();

    // Extending the box with 10 % in each dir (to handle a flat box).
    // Extend the trimmed plane with 10% in each dir.
    Point vec = high - low;
    double zero_tol = 1e-06;
    int ki;
    for (ki = 0; ki < 3; ++ki) {
      if (vec[ki] < zero_tol) {
	// We then add a randomly chosen value.
	double random_value = 1.0;
	low[ki] -= random_value;
	high[ki] += random_value;
      } else {
	low[ki] -= 0.1*vec[ki];
	high[ki] += 0.1*vec[ki];
      }
    }

    Point normal = plane->getNormal();

    MatrixXD<double, 3> cs;
    Vector3D from(normal[0],normal[1],normal[2]);
    Vector3D to(0.0, 0.0, 1.0);
    cs.setToRotation(to, from);
//     low = cs*low;
//     high = cs*high;
//     Point y2(low[0],high[1],low[2]); // Box corner.
//     Point x2(high[0],low[1],low[2]); // Box corner.

    MatrixXD<double, 3> cs1, cs2; // Rotation for the two axis directions.
    Vector3D from1(axis[0][0],axis[0][1],axis[0][2]);
    Vector3D to1(1.0, 0.0, 0.0);
    Vector3D from2(axis[1][0],axis[1][1],axis[1][2]);
    Vector3D to2(0.0, 1.0, 0.0);
    cs1.setToRotation(to1, from1);
    cs2.setToRotation(to2, from2);
    low = cs1*low;
    low = cs2*low;
    high = cs1*high;
    high = cs2*high;
    Point y2(low[0],high[1],low[2]); // Box corner.
    Point x2(high[0],low[1],low[2]); // Box corner.

    // Rotate back to xyz coordinate system
    MatrixXD<double, 3> cs3;
    cs3.setToRotation(from, to);

#endif

    int ki;
    double dist;
    double maxdist = -1;
    double mindist = 1.0e8;
    double frac = 0.001;
    for (ki=0; ki<(int)space_crvs.size(); ++ki) 
    {
	shared_ptr<SplineCurve> c1 =
	    (space_crvs[ki]->instanceType() == Class_SplineCurve) ?
	    dynamic_pointer_cast<SplineCurve>(space_crvs[ki]) :
	    shared_ptr<SplineCurve>(space_crvs[ki]->geometryCurve());
	if (c1)
	{
	    vector<double>::iterator start = c1->coefs_begin();
	    vector<double>::iterator end = c1->coefs_end();
	    for (; start != end; start+=3)
	    {
		dist = plane->distance(Point(start, start+3));
		maxdist = std::max(maxdist, dist);
		mindist = std::min(mindist, dist);
	    }
	}
    }
    if (maxdist >= 0.0 && maxdist - mindist < frac*maxdist)
    {
	MESSAGE("maxdist: " << maxdist);

	// Make rotated bounding box of trimming curves
	vector<Point> axis(3);
	plane->getSpanningVectors(axis[0],axis[1]);  // Two vectors
						 // spanning the plane
	BoundingBox box = space_crvs[0]->boundingBox();
	for (int i=1; i<(int)space_crvs.size(); ++i) {
	    box.addUnionWith(space_crvs[i]->boundingBox());
	}

	// Method does not handle rotation. We check whether the
	// normal for the curve-plane coincides with normal for
	// plane. Which is the case if box diagonal is orthogonal to
	// the plane (within a tolerance).
	Point low = box.low();
	Point high = box.high();
	Point diag = high - low;
	Point plane_normal = plane->getNormal();
	double ang = plane_normal.angle_smallest(diag);
	// Move trimmed plane
	Point pnt = space_crvs[0]->point(space_crvs[0]->startparam()); 
	Point proj_pt = plane->projectPoint(pnt);
	Point vec = pnt - proj_pt;
	Point loc = plane->getPoint();
	Point transl_loc = loc + vec;

	if (fabs(0.5*M_PI - ang) < 0.005)
	{
	    *plane = Plane(transl_loc, plane_normal);
// 	    std::cout << "vec: (" << vec[0] << ", " << vec[1] << ", " <<
// 		vec[2] << ")" << std::endl;
	}
	else
	{
	    MESSAGE("We needed to alter normal for plane.");
	    // We use spanning vectors and diag to create the new normal.
	    Point axis1, axis2;
	    plane->getSpanningVectors(axis1, axis2);
	    Point avg_axis = 0.5*(axis1 + axis2);
	    Point new_normal;
	    if (avg_axis.angle_smallest(diag) < 0.01)
		new_normal = axis1.cross(diag);
	    else
		new_normal = avg_axis.cross(diag);
	    *plane = Plane(transl_loc, new_normal);
	}
#if 0
	low += vec;
	high += vec;
	x2 += vec;
	y2 += vec;
#endif
    }

#if 0
    // We project low & high onto plane to get parameter values, as we
    // want to keep parametrization of plane.
    double umin, umax, vmin, vmax;
    double clo_dist_low, clo_dist_high;
    double epsgeo = 1e-06;
    Point clo_low, clo_high;
    plane->closestPoint(low, umin, vmin, clo_low, clo_dist_low, epsgeo);
    plane->closestPoint(high, umax, vmax, clo_high, clo_dist_high, epsgeo);
    if (umax < umin) {
      swap(low, x2);
      swap(y2, high);
      swap(umin, umax);
    }
    if (vmax < vmin) {
      swap(low, y2);
      swap(x2, high);
      swap(vmin, vmax);
    }
#endif

}


//===========================================================================
void BoundedUtils::fixInvalidBoundedSurface(shared_ptr<BoundedSurface>& bd_sf,
					    double max_tol_mult)
//===========================================================================
{
    if (max_tol_mult < 1.0)
    {
	max_tol_mult = 1.0;
    }

    int init_state = 0;
    bool sf_ok = bd_sf->isValid(init_state);

    if (sf_ok)
	return;

    // We allow to increase the tolerance for the boundary loops by
    // a factor.
    int nmb_seg_samples = 20;//100;
    double min_epsgeo = 1e-07;

    int nmb_loops = bd_sf->numberOfLoops();
    vector<double> init_loop_tol(nmb_loops, -1.0);
    vector<double> init_loop_sf_dist(nmb_loops, -1.0);
    vector<double> new_loop_sf_dist(nmb_loops, -1.0);
    for (int kj = 0; kj < nmb_loops; ++kj)
    {
	shared_ptr<CurveLoop> loop = bd_sf->loop(kj);
	init_loop_tol[kj] = loop->getSpaceEpsilon();
	init_loop_sf_dist[kj] =
	    bd_sf->maxLoopSfDist(kj, nmb_seg_samples);
    }

#ifdef SBR_DBG
    {
        std::ofstream outfile_curr_bd_sf("tmp/curr_bd_sf.g2");
        std::ofstream outfile_failures("tmp/bd_sf_failures.g2");
        SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_curr_bd_sf, 0.0);
    }
#endif
    if (init_state == 0)
    {
	MESSAGE("State: Failed analyzing input surface! "
		"Nothing more to be done. ");
#if 0
	bd_sf->writeStandardHeader(outfile_valid);
	bd_sf->write(outfile_valid);
#endif
    }
    else
    {
	// We first fix the trim cvs (i.e. project space cv).
	ASSERT(init_state < 0);
	// We first make sure that loop tol is set
	// appropriately.
	vector<double> allowed_tol(bd_sf->numberOfLoops(), -1.0);
	for (int kj = 0; kj < bd_sf->numberOfLoops(); ++kj)
	{
	    double loop_sf_dist =  bd_sf->maxLoopSfDist(kj);
	    // As a result of removing invalid parameter
	    // curves, the loop may seem further away, in
	    // which case we should make sure the
	    // tolerance is suitable.
	    double space_eps = bd_sf->loop(kj)->getSpaceEpsilon();
	    allowed_tol[kj] = space_eps*max_tol_mult;
	    if (loop_sf_dist > allowed_tol[kj])
	    {
		MESSAGE("Large dist from curve to sf! dist = " <<
			loop_sf_dist);
		if (allowed_tol[kj] > space_eps)
		{
		    MESSAGE("Altering to " << allowed_tol[kj]);
		    bd_sf->loop(kj)->setSpaceEpsilon(allowed_tol[kj]);
		    new_loop_sf_dist[kj] = allowed_tol[kj];
		}
	    }
	    else if (1.1*loop_sf_dist > space_eps)
	    {
		double new_space_eps = 1.1*loop_sf_dist;
                const double num_tol = 1.0e-14;
                if (new_space_eps - space_eps > num_tol)
                {
                    MESSAGE("Altering loop tol from space_eps = " <<
                            space_eps << ", to new_space_eps = " <<
                            new_space_eps);
                    bd_sf->loop(kj)->setSpaceEpsilon(new_space_eps);
                    new_loop_sf_dist[kj] = new_space_eps;
                }
	    }
	    else if (space_eps < min_epsgeo)
	    {
		MESSAGE("Altering loop tol from space_eps = " <<
			space_eps << ", to min_epsgeo = " <<
			min_epsgeo);
		bd_sf->loop(kj)->setSpaceEpsilon(min_epsgeo);
		new_loop_sf_dist[kj] = min_epsgeo;
	    }
	}

	int pos_state = -init_state;
	if ((pos_state%2 == 1))
	{// || (pos_state%4 > 1)) {
	    // Remove curves with mismatch between par and
	    // space cv,
	    MESSAGE("State: Mismatch for cvs, "
		    "trying to fix!");
	    // @@sbr072009 We should only a remove curve if
	    // it is invalid using allowed_tol
	    // (loop_tol*max_mult_tol).
	    double max_tol_mult = 1.0;
	    bd_sf->removeMismatchCurves(max_tol_mult);

	    // Make sure that tol for loop makes sense.
	    for (int kj = 0; kj < bd_sf->numberOfLoops(); ++kj)
	    {
		new_loop_sf_dist[kj] = bd_sf->maxLoopSfDist(kj);
		// As a result of removing invalid parameter
		// curves, the loop may seem further away, in
		// which case we should make sure the
		// tolerance is suitable.
		double space_eps = bd_sf->loop(kj)->getSpaceEpsilon();
		if (1.1*new_loop_sf_dist[kj] > allowed_tol[kj])
		{
		    MESSAGE("Large dist from curve to sf!, dist = " <<
			    new_loop_sf_dist[kj]);
		}
		else if (1.1*new_loop_sf_dist[kj] > space_eps)
		{
		    double new_space_eps = 1.1*new_loop_sf_dist[kj];
		    MESSAGE("Altering loop tol from space_eps = " <<
			    space_eps << ", to new_space_eps = " <<
			    new_space_eps);
		    bd_sf->loop(kj)->setSpaceEpsilon(new_space_eps);
		}
	    }

	    bd_sf->analyzeLoops();
	    int bd_sf_state = 0;
	    sf_ok = bd_sf->isValid(bd_sf_state);

	    pos_state = -bd_sf_state;
	    if (!sf_ok && pos_state%2 == 1)
	    {
		MESSAGE("State: Failed removing inconsistent curves!");
#ifdef SBR_DBG
                {
                    std::ofstream outfile_failures("tmp/bd_sf_failures.g2");
                    SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
                }
#endif
// 			  bd_sf->writeStandardHeader(outfile_failures);
// 			  bd_sf->write(outfile_failures);
// 		++nmb_failures;
		return;
	    }
	}
	if (!sf_ok && (pos_state%4 > 1))
	{

#ifdef SBR_DBG
            {
                std::ofstream outfile_not_ok("tmp/bd_sf_not_ok.g2");
                SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_not_ok, 0.0);
            }
#endif
	    // Project missing parameter curves.
	    MESSAGE("State: Missing par cv, trying to fix!");
	    // There is no point in projecting missing parameter curves
	    // if existing curves are not within input tolerance.
	    // We check if we need to enlarge epsgeo.
            const double gap = 1.0e-03;
            const double neighbour = 1.0e-02;
            const double kink = 1.0e-02;
	    CreatorsUtils::fixTrimCurves(bd_sf, 1.0, gap, neighbour, kink);
	    bd_sf->analyzeLoops();
	    int bd_sf_state = 0;
	    sf_ok = bd_sf->isValid(bd_sf_state);

	    pos_state = -bd_sf_state;
	    if (!sf_ok && pos_state%4 > 1) {
		MESSAGE("State: Failed projecting (classType: " <<
			bd_sf->underlyingSurface()->instanceType()
			<< ")!");
// 			  bd_sf->writeStandardHeader(outfile_failures);
// 			  bd_sf->write(outfile_failures);
// 		writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
		return;
	    }
	    else if (pos_state%2 > 1)
	    {
		MESSAGE("State: Failed. Projection not a valid loop!");
#ifdef SBR_DBG
                {
                    std::ofstream outfile_failures("tmp/bd_sf_failures.g2");
                    SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
                }
#endif
		return;
	    }
	}
	if ((pos_state%8 > 1) || (pos_state%16 > 1))
	{
	    // Loops not closed or in wrong order/orientation.
	    double max_gap = -1.0;
	    bool success = bd_sf->fixInvalidSurface(max_gap);
	    if (!success)
	    {
		MESSAGE("max_gap = " << max_gap);
	    }
	}

	// We are done, time to sum up.
	bd_sf->analyzeLoops();
	int bd_sf_state = 0;
	bool sf_valid = bd_sf->isValid(bd_sf_state);
	if (!sf_valid)
	{
	    MESSAGE("State: Obj not valid after fixing! "
		    "sf_state: " << bd_sf_state);
#ifdef SBR_DBG
            {
                std::ofstream outfile_failures("tmp/bd_sf_failures.g2");
                SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
            }
#endif
	}
	else
	{
	    // Writing to file the fixed surfaces.
	    MESSAGE("State: Surface valid after fixing "
		    "trim curves! bd_sf_state = " << bd_sf_state);
#ifdef SBR_DBG
            std::ofstream outfile_fixed("tmp/bd_sf_fixed.g2");
	    bd_sf->writeStandardHeader(outfile_fixed);
	    bd_sf->write(outfile_fixed);
#endif
	}
    }
}

//===========================================================================
bool BoundedUtils::loopIsDegenerate(vector<shared_ptr<CurveOnSurface> >& loop,
				    double epsgeo)
//===========================================================================
{
  // Assumes that no curve in the loop intersect itself
  // For each curve in the loop, compute a number of sampling points in
  // the inner and try to project these points onto the other loop curves
  int max_nmb_sample = 5;
  int min_nmb_sample = 1;
  bool degen = true; // Initially
  double len = 0.0;
  Point prev;
  if (loop.size() == 0)
    return true;

  int dim = loop[0]->underlyingSurface()->dimension();
  for (size_t ki=0; ki<loop.size(); ++ki)
    {
      double cv_len = (dim == 1) ?
	loop[ki]->parameterCurve()->estimatedCurveLength() :
	loop[ki]->estimatedCurveLength();
      int nmb_sample = (int)(0.05*cv_len/epsgeo);
      nmb_sample = std::min(std::max(nmb_sample, min_nmb_sample), max_nmb_sample);
      double t1 = loop[ki]->startparam();
      double t2 = loop[ki]->endparam();
      double tdel = (t2 - t1)/(double)(nmb_sample+1);
      double tpar;
      int kj;
      for (kj=0, tpar=t1+tdel; kj<nmb_sample; ++kj, tpar+=tdel)
	{
	  // Evaluate sampling point
	  Point pos;
	  if (dim == 1)
	    pos = loop[ki]->faceParameter(tpar);
	  else
	    loop[ki]->point(pos, tpar);

	  if (prev.dimension() == pos.dimension())
	    len += prev.dist(pos);

	  // Project
	  size_t kr;
	  for (kr=0; kr<loop.size(); ++kr)
	    {
	      if (kr == ki)
		continue;  // Same curve

	      double tpar2, dist;
	      Point pos2;
	      if (dim == 1)
		loop[kr]->parameterCurve()->closestPoint(pos, 
							 loop[kr]->startparam(),
							 loop[kr]->endparam(), 
							 tpar2, pos2, dist);
	      else
		loop[kr]->closestPoint(pos, loop[kr]->startparam(),
				       loop[kr]->endparam(), tpar2,
				       pos2, dist);

	      if (dist < epsgeo)
		break;
	    }

	  if (kr == loop.size())
	    {
	      // The loop is not degenerate
	      degen = false;
	    }
	  prev = pos;
	}
    }
  if (len < epsgeo)
    degen = true;

  // All sampling points coincide with some other part of the loop
  return degen;
}


//==========================================================================
bool BoundedUtils::createMissingParCvs(Go::BoundedSurface& bd_sf)
//==========================================================================
{
    bool all_par_cvs_ok = true;

    double max_gap = bd_sf.maxLoopGap();
    double epsgeo = bd_sf.getEpsGeo();
    if (max_gap > epsgeo)
    {
	MESSAGE("The epgeo should be increased! epsgeo = " << epsgeo << ", max_gap = " << max_gap);
    }
    else
    {
	;//MESSAGE("All OK, epsgeo = " << epsgeo << ", max_gap = " << max_gap);
    }

#ifndef NDEBUG
    {
	std::ofstream debug("tmp/debug_pre.g2");
	SplineDebugUtils::writeTrimmedInfo(bd_sf, debug);
	std::ofstream debug2("tmp/seam_info.g2");
	SplineDebugUtils::writeSeamInfo(bd_sf, debug2);
	double debug_val = 0.0;
    }
#endif // NDEBUG

    if (bd_sf.underlyingSurface()->instanceType() == Class_SurfaceOfLinearExtrusion) {
        // Unbounded surface of linear extrusion is not handled well in the messy projection routine.
        // Unbounded handling is restricted to elementary surfaces only.
        shared_ptr<SurfaceOfLinearExtrusion> surf_of_lin_extr =
            dynamic_pointer_cast<SurfaceOfLinearExtrusion>(bd_sf.underlyingSurface());
        RectDomain cont_dom = surf_of_lin_extr->containingDomain();
        double max_domain_val = 1.0e08;
        double umin = cont_dom.umin();
        double umax = cont_dom.umax();
        double vmin = cont_dom.vmin();
        double vmax = cont_dom.vmax();
        RectDomain bounded_cont_dom;
        if ((umin < -max_domain_val) || (umax > max_domain_val) ||
            (vmin < -max_domain_val) || (vmax > max_domain_val))
        {
            umin = std::max(umin, -max_domain_val);
            vmin = std::max(vmin, -max_domain_val);
            umax = std::min(umax, max_domain_val);
            vmax = std::min(vmax, max_domain_val); 
            surf_of_lin_extr->setParameterBounds(umin, vmin, umax, vmax);
        }
    }

    for (int ki = 0; ki < bd_sf.numberOfLoops(); ++ki)
    {
        shared_ptr<CurveLoop> loop = bd_sf.loop(ki);
	const bool loop_is_ccw = (ki == 0);
        bool loop_par_cvs_ok = createMissingParCvs(*loop, loop_is_ccw);
        if (!loop_par_cvs_ok)
        {
            all_par_cvs_ok = false;
        }
    }

#ifndef NDEBUG
    {
	if (!all_par_cvs_ok)
	{
	    MESSAGE("all_par_cvs_ok: " << all_par_cvs_ok);
	    std::ofstream debug("tmp/debug_post.g2");
	    Go::SplineDebugUtils::writeTrimmedInfo(bd_sf, debug);
	    double debug_val = 0.0;
	}
    }
#endif // NDEBUG
 
    return all_par_cvs_ok;
}


//==========================================================================
bool BoundedUtils::createMissingParCvs(vector<CurveLoop>& bd_loops)
//==========================================================================
{
    bool all_par_cvs_ok = true;

    int loop_id = -1;
    for ( auto& bd_loop : bd_loops )
    {
        ++loop_id;
	const bool loop_is_ccw = (loop_id == 0);        
        bool loop_cvs_ok = createMissingParCvs(bd_loop, loop_is_ccw);
        if (!loop_cvs_ok)
        {
            all_par_cvs_ok = false;
        }
    }

    return all_par_cvs_ok;
}

//==========================================================================
bool BoundedUtils::createMissingParCvs(CurveLoop& bd_loop, bool loop_is_ccw)
//==========================================================================
{
    bool all_par_cvs_ok = true;

    vector<pair<shared_ptr<Point>, shared_ptr<Point> > > loop_end_par_pts = getEndParamPoints(bd_loop, loop_is_ccw);

    // If an end point is in a singularity and we are missing a deg space curve we must add one.
    // We are assuming that all projected end points are within epsgeo. If the corresponding par pts
    // differ by more than epspar the point is assumed to be a singularity in need of a degenerete
    // segment in space and a corresponding linear segment in the paremeter domain.
    shared_ptr<CurveOnSurface> first_cv = dynamic_pointer_cast<CurveOnSurface>(bd_loop[0]);
    ASSERT(first_cv.get() != NULL);
    shared_ptr<ParamSurface> sf = first_cv->underlyingSurface();
    const double epsgeo = bd_loop.getSpaceEpsilon();
    Point par_eps = SurfaceTools::getParEpsilon(*sf, epsgeo);
    const double epspar = std::min(par_eps[0], par_eps[1]);
    const double epstang = 1.0e-06;
    std::vector<shared_ptr<ParamCurve> > loop_cvs = bd_loop.getCurves();
    for (size_t ki = 0; ki < loop_end_par_pts.size(); ++ki)
    {
        size_t next_ind = (ki + 1)%(loop_end_par_pts.size());
        if ((loop_end_par_pts[ki].second.get() != nullptr) && (loop_end_par_pts[next_ind].first.get() != nullptr))
        {
            shared_ptr<Point> end_par_pt = loop_end_par_pts[ki].second;
            vector<Point> end_sf_pt = sf->point((*end_par_pt)[0], (*end_par_pt)[1], 1);
            shared_ptr<Point> start_par_pt = loop_end_par_pts[next_ind].first;
            vector<Point> start_sf_pt = sf->point((*start_par_pt)[0], (*start_par_pt)[1], 1);
            double pardist = end_par_pt->dist(*start_par_pt);
            double space_dist = end_sf_pt[0].dist(start_sf_pt[0]);
            if (pardist > epspar)
            {
                if (space_dist < epsgeo)
                {   // This warning means one of the following:

                    // If the surface point is a singularity we add a degenerate segment.
                    const double par_dist_0 = fabs((*end_par_pt)[0] - (*start_par_pt)[0]);
                    const double par_dist_1 = fabs((*end_par_pt)[1] - (*start_par_pt)[1]);
                    if ((par_dist_0 < epspar) || (par_dist_1 < epspar))
                    {
                        // If the tangents for the index wich the deviating parameter value is less
                        // than a small tolerance we have found a surface singularity and we need to
                        // insert a degenerate segment.
                        const int diff_ind = (par_dist_0 < epspar) ? 1 : 0;
                        const double end_length = end_sf_pt[1 + diff_ind].length();
                        const double start_length = start_sf_pt[1 + diff_ind].length();
                        if ((end_length < epstang) && (start_length < epstang))
                        {
                            // 1) We miss a degenerate edge (i.e. the surface point is a singularity, the par
                            //    points share parameter value in one of the directions).
                            // std::cout << "WARNING: Suspecting: Add a degenerate edge! pardist = " <<
                            //     pardist << " (epspar = " << epspar << "). end_length: " << end_length <<
                            //     ", start_length: " << start_length << ". UPDATE BD_SF WITH LOOP!" << std::endl;
                            // The Line object relies on a non-zero directional vector, hence we use a SplineCurve.
#if 1
                            shared_ptr<Line> deg_line(new Line(end_sf_pt[0], start_sf_pt[0], 0.0, 1.0));
#else
                            shared_ptr<SplineCurve> deg_line(new SplineCurve(end_sf_pt[0], 0.0, start_sf_pt[0], 1.0));
#endif
                            const double par_pref = false;
                            shared_ptr<CurveOnSurface> deg_cv_on_sf(new CurveOnSurface(sf, deg_line, par_pref));
                            loop_end_par_pts.insert(loop_end_par_pts.begin() + ki + 1,
                                                    std::make_pair(loop_end_par_pts[ki].second,
                                                                   loop_end_par_pts[next_ind].first));
                            loop_cvs.insert(loop_cvs.begin() + ki + 1, deg_cv_on_sf);
                            continue;
                        }
                        else
                        {
                            // 2) We failed projection onto the seam of a closed surface (like sphere, cylinder, torus).
                            std::cout << "WARNING: Suspecting: Failed projecting onto seam of closed surface." <<
                                std::endl;
                        }
                    }
                    else
                    {
                        // 3) The projection routine is not accurate enough.
                        std::cout << "WARNING: Suspecting: Projection is inaccurate." << " par_dist_0: " <<
                            par_dist_0 << ", par_dist_1: " << par_dist_1 << ", epspar: " << epspar <<
                            ", space_dist: " << space_dist << ", epsgeo: " << epsgeo << std::endl;
                    }
                }
                else
                {
                    // 4) The space curve is too far from the surface.
                    std::cout << "WARNING: Suspecting: The loop is not connected! space_dist = " <<
                        space_dist << ", epsgeo = " << epsgeo << ")" << std::endl;
                }
            }
        }
    }

#if 1
    if (bd_loop.size() != loop_cvs.size())
    {
        bd_loop.setCurves(loop_cvs, false);
    }
#endif

    int num_segments = bd_loop.size();
    vector<int> loop_cv_ind(num_segments);
    vector<bool> failed_once(num_segments, false);
    for (int kr = 0; kr < bd_loop.size(); ++kr)
    {
        loop_cv_ind[kr] = kr;
    }

    // We start by creating parameter end points.
    // If that fails for a curve there is no need to even try to project ...
    // We may also increase the epsgeo if the lifted end par pts are too far apart.

    // If a parameter curve is missing we try to project from the space curve.
    for (size_t ki=0; ki< loop_cv_ind.size(); ++ki) {
        // Try to generate the parameter curve if it does not
        // exist already
        const int curr_cv_ind = loop_cv_ind[ki];
        shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(bd_loop[curr_cv_ind]);
        ASSERT(cv_on_sf.get() != NULL);
        if (cv_on_sf->parameterCurve().get() != NULL)
        {   // Parameter curve is present.
            continue;
        }

        shared_ptr<Point> start_pt = loop_end_par_pts[curr_cv_ind].first;
        shared_ptr<Point> end_pt = loop_end_par_pts[curr_cv_ind].second;

        int num_loop_cvs = bd_loop.size();
        const int prev_cv_ind = (curr_cv_ind-1+num_loop_cvs)%num_loop_cvs;
        shared_ptr<ParamCurve> prev_cv = bd_loop[prev_cv_ind];
        shared_ptr<CurveOnSurface> prev_cos = dynamic_pointer_cast<CurveOnSurface>(prev_cv);
        if ((start_pt.get() == NULL) && (prev_cos->parameterCurve()))
        {
            start_pt = shared_ptr<Point>
                (new Point(prev_cos->parameterCurve()->point(prev_cos->parameterCurve()->endparam())));
        }
        const int next_cv_ind = (curr_cv_ind+1)%num_loop_cvs;
        shared_ptr<ParamCurve> next_cv = bd_loop[next_cv_ind];
        shared_ptr<CurveOnSurface> next_cos = dynamic_pointer_cast<CurveOnSurface>(next_cv);
        if ((end_pt.get() == NULL) && (next_cos->parameterCurve()))
        {
            end_pt = shared_ptr<Point>
                (new Point(next_cos->parameterCurve()->point(next_cos->parameterCurve()->startparam())));
        }

        if ((start_pt.get() == NULL) || (end_pt.get() == NULL))
        {
            if (!failed_once[curr_cv_ind])
            {
                failed_once[curr_cv_ind] = true;
                loop_cv_ind.erase(loop_cv_ind.begin() + ki);
                loop_cv_ind.push_back(curr_cv_ind);
                --ki;
                continue;
            }
        }

        shared_ptr<ParamSurface> under_sf = cv_on_sf->underlyingSurface();
        shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();

#ifndef NDEBUG
        {
            std::ofstream debug3("tmp/undersf_outer_loop.g2");
            Go::SplineDebugUtils::writeOuterBoundaryLoop(*under_sf, debug3);
            double debug_val = 0.0;
        }
#endif // NDEBUG

        RectDomain cont_dom = under_sf->containingDomain();
        double max_domain_val = 1.0e06;
        double umin = cont_dom.umin();
        double umax = cont_dom.umax();
        double vmin = cont_dom.vmin();
        double vmax = cont_dom.vmax();
        RectDomain* domain_of_interest = NULL;
        RectDomain bounded_cont_dom;
        if ((umin < -max_domain_val) || (umax > max_domain_val) ||
            (vmin < -max_domain_val) || (vmax > max_domain_val))
        {
            Array<double, 2> ll, ur;
            ll[0] = std::max(umin, -max_domain_val);
            ll[1] = std::max(vmin, -max_domain_val);
            ur[0] = std::min(umax, max_domain_val);
            ur[1] = std::min(vmax, max_domain_val); 
            bounded_cont_dom = RectDomain(ll, ur);
            domain_of_interest = &bounded_cont_dom;
        }
        bool cv_ok = cv_on_sf->ensureParCrvExistence(epsgeo, domain_of_interest, start_pt.get(), end_pt.get());

#ifndef PROJECT_CURVES_DEBUG
        const int num_samples = 1000;
        bool same_trace = cv_on_sf->sameTrace(epsgeo, num_samples);
        double max_trace_diff = cv_on_sf->maxTraceDiff(num_samples);
        double debug_val = 0.0;
#endif

#if 1
        if (!cv_ok)
        {
            if (failed_once[curr_cv_ind])
            {
                all_par_cvs_ok = false;
            }
            else
            {
                failed_once[curr_cv_ind] = true;
                loop_cv_ind.erase(loop_cv_ind.begin() + ki);
                loop_cv_ind.push_back(curr_cv_ind);
                --ki;
                continue;
            }
        }
#else
        if (!cv_ok)
        {
            all_par_cvs_ok = false;
        }
#endif
    }

    if (all_par_cvs_ok)
    {
        bd_loop.analyze();
    }

    return all_par_cvs_ok;
}


//===========================================================================
vector<pair<shared_ptr<Point>, shared_ptr<Point> > >
BoundedUtils::getEndParamPoints(const Go::CurveLoop& bd_loop, bool ccw_loop)
//===========================================================================
{
    const int num_segs = (int)bd_loop.size();
    vector<pair<shared_ptr<Point>, shared_ptr<Point> > > bd_par_pts(num_segs); // The bd of the curve, i.e. start and end.

    const double epsgeo = bd_loop.getSpaceEpsilon();

    const bool loop_is_ccw = ccw_loop;
    const bool loop_is_cw = !ccw_loop;

    // We try to find end param points for all segments.  It should be straightforward except for cases
    // where the segment has an end point at a surface seam (closed surface), with a tangent parallel to
    // the seam.
    for (int ki = 0; ki < num_segs; ++ki)
    {
	// We expect the bd_loop to consist of CurveOnSurface's.
	shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(bd_loop[ki]);
        int prev_ind = (ki - 1 + num_segs)%num_segs;
	shared_ptr<CurveOnSurface> prev_cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(bd_loop[prev_ind]);
        int next_ind = (ki + 1)%num_segs;
	shared_ptr<CurveOnSurface> next_cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(bd_loop[next_ind]);
	assert((cv_on_sf.get() != NULL) && (prev_cv_on_sf.get() != NULL) && (next_cv_on_sf.get() != NULL));
	shared_ptr<ParamCurve> pcv = cv_on_sf->parameterCurve();
	shared_ptr<ParamCurve> prev_pcv = prev_cv_on_sf->parameterCurve();
	shared_ptr<ParamCurve> next_pcv = next_cv_on_sf->parameterCurve();
	if (pcv.get() != NULL)
	{
	    bd_par_pts[ki] = std::make_pair(shared_ptr<Point>(new Point(pcv->point(pcv->startparam()))),
                                            shared_ptr<Point>(new Point(pcv->point(pcv->endparam()))));
	}
	else
	{
	    double* seed = NULL;
	    bool check_bd = true;
	    bd_par_pts[ki].first =
                cv_on_sf->projectSpacePoint(cv_on_sf->startparam(), epsgeo,
                                            seed,
                                            loop_is_ccw, loop_is_cw,
                                            check_bd);
	    bd_par_pts[ki].second =
                cv_on_sf->projectSpacePoint(cv_on_sf->endparam(), epsgeo,
                                            seed,
                                            loop_is_ccw, loop_is_cw,
                                            check_bd);

	    double dist = -1.0;
	    if (bd_par_pts[ki].first.get() != NULL)
	    {
		double dist_start = cv_on_sf->spaceDist(cv_on_sf->startparam(), *(bd_par_pts[ki].first));
		if (dist_start > dist)
		{
		    dist = dist_start;
		}
	    }
	    if (bd_par_pts[ki].second.get() != NULL)
	    {
		double dist_end = cv_on_sf->spaceDist(cv_on_sf->endparam(), *(bd_par_pts[ki].second));
		if (dist_end > dist)
		{
		    dist = dist_end;
		}
	    }
	    if (dist > epsgeo)
	    {
		MESSAGE("Inconsistent input to curve approximation: dist = "
			<< dist << ", epsgeo = " << epsgeo);// << ". Altering epsgeo to: " << 1.1*dist);
		//epsgeo = 1.1*dist;
	    }
	}
    }

    // We take another round to see if we failed for any curves. Using
    // neighbour curve pt as seed if it exists. If that also is
    // missing we must use a marching approach.
    // We start with a segment for which there exists info at start (possibly from end of previous segment).
    int start_ind = 0;
    for (int ki = 0; ki < num_segs; ++ki)
    {
	int prev_ind = (ki - 1 + num_segs)%num_segs;
	if ((bd_par_pts[ki].first.get() != NULL) || (bd_par_pts[prev_ind].second.get() != NULL))
	{
	    start_ind = ki;
	    break;
	}
    }

    for (int ki = 0; ki < num_segs; ++ki)
    {
      int curr_ind = (start_ind + ki)%num_segs;
	shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(bd_loop[curr_ind]);
	int prev_ind = (curr_ind - 1 + num_segs)%num_segs;
	int next_ind = (curr_ind + 1)%num_segs;
	if (bd_par_pts[curr_ind].first.get() == NULL)
	{
	    if (bd_par_pts[prev_ind].second.get() != NULL)
	    {
		Point seed_pt = *(bd_par_pts[prev_ind].second);
		bool check_bd = true;
		//int follows_seem_ind = -1;
		bd_par_pts[curr_ind].first = cv_on_sf->projectSpacePoint(cv_on_sf->startparam(), epsgeo,
									 &seed_pt[0],
//									 follows_seem_ind,
									 loop_is_ccw, loop_is_cw,
									 check_bd);
		// if (follows_seem_ind != -1)
		// {
		//     MESSAGE("We should look for the next/prev segment with a well defined par pt!");
		// }
	    }
	    else
	    {
		;//MESSAGE("Case not handled, time to add the marching approach!");
	    }
	}

	if (bd_par_pts[curr_ind].second.get() == NULL)
	{
	    if (bd_par_pts[next_ind].first.get() != NULL)
	    {
		Point seed_pt = *(bd_par_pts[next_ind].first);
		bool check_bd = true;
		//int follows_seem_ind = -1;
		bd_par_pts[curr_ind].second = cv_on_sf->projectSpacePoint(cv_on_sf->endparam(), epsgeo,
									  &seed_pt[0],
//									  follows_seem_ind,
									  loop_is_ccw, loop_is_cw,
									  check_bd);
		// if (follows_seem_ind != -1)
		// {
		//     MESSAGE("We should look for the next/prev segment with a well defined par pt!");
		// }
	    }
	    else
	    {
		;//MESSAGE("Case not handled, time to add the marching approach!");
	    }
	}
    }

    return bd_par_pts;
}


shared_ptr<Point> BoundedUtils::projectSpacePoint(const ParamSurface& sf,
                                                  const ParamCurve& space_cv,
                                                  double tpar, double epsgeo,
                                                  double* seed,
                                                  bool ccw_loop,
                                                  bool cw_loop ,
                                                  bool check_bd)
{

    // If we give info on loop orientation, we must choose a direction.
    if (ccw_loop && cw_loop)
    {
        THROW("The caller can not choose both ccw and cw direction for the loop.");
    }

    const double ang_tol = 1e-02;

    bool closed_dir_u, closed_dir_v;
    Go::SurfaceTools::checkSurfaceClosed(sf, closed_dir_u, closed_dir_v, epsgeo);
    bool is_closed = (closed_dir_u || closed_dir_v);

    vector<Point> cv_pt = space_cv.point(tpar, 1);
    double clo_u, clo_v, clo_dist;
    Point clo_pt;
    const double eps = 1e-10;
    sf.closestPoint(cv_pt[0], clo_u, clo_v, clo_pt, clo_dist, eps, NULL, seed);
    vector<Point> sf_pt = sf.point(clo_u, clo_v, 1);
    const double deg_tol = 1.0e-06;
    const double length_cv_der = cv_pt[1].length();
    const bool deg_cv_pt = (length_cv_der < deg_tol);
    const double length_uder = sf_pt[1].length();
    const double length_vder = sf_pt[2].length();
    const bool deg_uder = (length_uder < deg_tol);
    const bool deg_vder = (length_vder < deg_tol);
    if (deg_uder || deg_vder)
    {
        if (deg_cv_pt) // If the curve is degenerate in the point the projection must be handled by other
                       // methods (like using parameters of neighbouring curves if this is a bd point).
        {
            return shared_ptr<Point>(NULL);
        }

        if (seed)
        { // @@sbr201711 We must decide on how to handle this situation. Trust seed to give satisfactory param value?
            std::cout << "DEBUG: We were given a seed!" << std::endl;
        }

        // We need to use a marching approach to find the correct parameter. Or use the space
        // tangent. For the cone case this should suffice. The same with the sphere.
        double ang_rad = (deg_uder) ? cv_pt[1].angle(sf_pt[2]) : cv_pt[1].angle(sf_pt[1]);
        std::cout << "DEBUG: The surface is degenerate in this point! ang_rad = " << ang_rad << std::endl;
        double tstep = 1.0e-03;
        double tpar2 = (tpar + tstep < space_cv.endparam()) ? tpar + tstep : tpar - tstep;
        vector<Point> cv_pt2 = space_cv.point(tpar2, 1);
        double clo_u2, clo_v2, clo_dist2;
        Point clo_pt2;
        sf.closestPoint(cv_pt2[0], clo_u2, clo_v2, clo_pt2, clo_dist2, eps, NULL, seed);
        std::cout << "DEBUG: clo_u: " << clo_u << ", clo_u2: " << clo_u2 << ", clo_v: " << clo_v <<
            ", clo_v2: " << clo_v2 << std::endl;

        // @@sbr201801 Another option is to use the angle with the end tangents at both ends of the deg
        // edge. Or even better we may actually search for the parameter with the corresponding
        // direction! If the end tangents are the same we must use a marching approach to find the
        // correct parameter along the edge.
        const RectDomain& rect_dom = sf.containingDomain();
        if (deg_uder)
        {
            Point par_pt_min(rect_dom.umin(), clo_v);
            Point par_pt_max(rect_dom.umax(), clo_v);
            // We then compute the directional derivatives.
            vector<Point> pt_min = sf.point(par_pt_min[0], par_pt_min[1], 1);
            double ang_min = cv_pt[1].angle(pt_min[2]);
            bool min_parallel = ((fabs(ang_min) < ang_tol) || (fabs(fabs(ang_min) - M_PI) < ang_tol));
            vector<Point> pt_max = sf.point(par_pt_max[0], par_pt_max[1], 1);
            double ang_max = cv_pt[1].angle(pt_max[2]);
            bool max_parallel = ((fabs(ang_max) < ang_tol) || (fabs(fabs(ang_max) - M_PI) < ang_tol));
            if (min_parallel != max_parallel)
            {
                clo_u = (min_parallel) ? par_pt_min[0] : par_pt_max[0];
                return shared_ptr<Point>(new Point(clo_u, clo_v));
            }
            else if (min_parallel && max_parallel)
            {
                // We use cw/ccw info if it exists.
                if (ccw_loop || cw_loop)
                {
                    bool same_dir = (fabs(ang_max) < 0.5*M_PI);
                    clo_u = ((same_dir && ccw_loop) || ((!same_dir) && cw_loop)) ? par_pt_max[0] : par_pt_min[0];
                    return shared_ptr<Point>(new Point(clo_u, clo_v));
                }
            }
            else if ((!min_parallel) && (!max_parallel))
            {
                std::cout << "WARNING: Tangent to/from degenerate point is not along the min/max value!" << std::endl;
                // We should add a search in the tanget space of the surface. We are not guaranteed to follow an iso
                // line, making the search more complex.
            }
        }
        else
        {
            Point par_pt_min(clo_u, rect_dom.vmin());
            Point par_pt_max(clo_u, rect_dom.vmax());
            // We then compute the directional derivatives.
            vector<Point> pt_min = sf.point(par_pt_min[0], par_pt_min[1], 1);
            double ang_min = cv_pt[1].angle(pt_min[1]);
            bool min_parallel = ((fabs(ang_min) < ang_tol) || (fabs(fabs(ang_min) - M_PI) < ang_tol));
            vector<Point> pt_max = sf.point(par_pt_max[0], par_pt_max[1], 1);
            double ang_max = cv_pt[1].angle(pt_max[1]);
            bool max_parallel = ((fabs(ang_max) < ang_tol) || (fabs(fabs(ang_max) - M_PI) < ang_tol));
            if (min_parallel != max_parallel)
            {
                clo_v = (min_parallel) ? par_pt_min[1] : par_pt_max[1];
                return shared_ptr<Point>(new Point(clo_u, clo_v));
            }
            else if (min_parallel && max_parallel)
            {
                // We use cw/ccw info if it exists.
                if (ccw_loop || cw_loop)
                {
                    bool same_dir = (fabs(ang_max) < 0.5*M_PI);
                    clo_v = ((same_dir && ccw_loop) || ((!same_dir) && cw_loop)) ? par_pt_min[1] : par_pt_max[1];
                    return shared_ptr<Point>(new Point(clo_u, clo_v));
                }
            }
            else if ((!min_parallel) && (!max_parallel))
            {
                std::cout << "WARNING: Tangent to/from degenerate point is not along the min/max value!" << std::endl;
                // We should add a search in the tanget space of the surface. We are not guaranteed to follow an iso
                // line, making the search more complex.
            }
        }
        double upar = (deg_uder) ? clo_u2 : clo_u;
        double vpar = (deg_vder) ? clo_v2 : clo_v;
        Point sf_pt_deg = sf.point(upar, vpar);
        double sf_pt_deg_dist = sf_pt_deg.dist(cv_pt[0]); // We compare with the projection, not with the
        // curve (which may be relatively far away).
        if (sf_pt_deg_dist < epsgeo) // We only accept the point if it is within epsgeo.
        {
            std::cout << "DEBUG: Degenerate point, enabling special handling! clo_dist: " << clo_dist <<
                ", sf_pt_deg_dist: " << sf_pt_deg_dist << std::endl;
            if (deg_uder)
            {
                clo_u = clo_u2;
            }
            else
            {
                clo_v = clo_v2;
            }
            return shared_ptr<Point>(new Point(clo_u, clo_v));
        }
    }

    bool sf_is_bounded = sf.isBounded();
    if (sf_is_bounded && check_bd)
    {
	try {
	    double clo_u_bd, clo_v_bd, clo_dist_bd;
	    Point clo_pt_bd;
	    sf.closestBoundaryPoint(cv_pt[0], clo_u_bd, clo_v_bd, clo_pt_bd, clo_dist_bd, eps, NULL, seed);
	    if (clo_dist_bd < clo_dist)
	    {
		clo_dist = clo_dist_bd;
		clo_u = clo_u_bd;
		clo_v = clo_v_bd;
		clo_pt = clo_pt_bd;
	    }
	}
	catch (...)
	{
	    MESSAGE("Suspecting the surface is not bounded.");
	}
    }

    const Point sf_epspar = SurfaceTools::getParEpsilon(sf, epsgeo);
    const double epspar = std::min(sf_epspar[0], sf_epspar[1]);
    const double knot_diff_tol = epspar;//1e-08;
    const RectDomain rect_dom = sf.containingDomain();
    const double umin = rect_dom.umin();
    const double umax = rect_dom.umax();
    const double vmin = rect_dom.vmin();
    const double vmax = rect_dom.vmax();
    const bool at_u_start = (fabs(clo_u - umin) < knot_diff_tol);
    const bool at_u_end = (fabs(clo_u - umax) < knot_diff_tol);
    const bool at_v_start = (fabs(clo_v - vmin) < knot_diff_tol);
    const bool at_v_end = (fabs(clo_v - vmax) < knot_diff_tol);

    // By at_u_bd we mean that the seam corresponds to a u-parameter.
    const bool at_u_bd = (at_u_start || at_u_end);
    const bool at_v_bd = (at_v_start || at_v_end);

    if ((seed != NULL) || ((!(closed_dir_u && at_u_bd) && !(closed_dir_v && at_v_bd))))// && (!deg_uder && !deg_vder)))
    {
	// Simple case.
	shared_ptr<Point> par_pt(new Point(clo_u, clo_v));
	return par_pt;
    }
    else // We are at the seam. We use the direction of the projected tangent to choose the side.
    {
	// We consider an angle of more than ang_tol to be non-tangential.
	const bool handle_u_seam = (at_u_bd && closed_dir_u);
	const bool handle_v_seam = (at_v_bd && closed_dir_v);

	// If we are at the end of the curve, our test differs slightly.
	const bool at_cv_end = (fabs(tpar - space_cv.endparam()) < knot_diff_tol);
	const bool at_cv_start = (fabs(tpar - space_cv.startparam()) < knot_diff_tol);

#if 0
	// If we cross the seem the task is impossible. Use a seed.
	// The caller can get a seed on both sides by picking a tpar
	// slightly larger and smaller.
	if (!at_cv_start && !at_cv_end && !cw_loop && !ccw_loop)
	{
	    // std::cout << "WARNING: Case requires a seed." << std::endl;
	    return shared_ptr<Point>(NULL);
	}
#endif

	vector<Point> sf_pt2 = sf.point(clo_u, clo_v, 1);
	// vector<Point> cv_pt = space_cv.point(tpar, 1);
	double ang_u_space = cv_pt[1].angle(sf_pt2[1]);
	double ang_v_space = cv_pt[1].angle(sf_pt2[2]);

	// The range of angle2() is [0, 2*M_PI).
	Point u_dir(1.0, 0.0);
	Point v_dir(0.0, 1.0);

	// We project the space curve tangent onto the surface.
	Point proj_par_pt(clo_u, clo_v);
	Point par_tangent = projectSpaceCurveTangent(sf, space_cv, proj_par_pt, tpar);

	double ang_u = u_dir.angle2(par_tangent); // This is in the parameter domain.
	// But we want [-M_PI, M_PI).
	if (ang_u >= M_PI)
	{
	    ang_u -= 2*M_PI;
	}
	// To fix problems with uneven scaling of domain directions we check angle for space tangents.
	if (ang_u_space < ang_tol)
	{
	    ang_u = 0.0;
	}
	double ang_v = v_dir.angle2(par_tangent);
	if (ang_v >= M_PI)
	{
	    ang_v -= 2*M_PI;
	}
	if (ang_v_space < ang_tol)
	{
	    ang_v = 0.0;
	}
        bool u_parallel = ((fabs(ang_u) < ang_tol) || (fabs(fabs(ang_u) - M_PI) < ang_tol));
        bool v_parallel = ((fabs(ang_v) < ang_tol) || (fabs(fabs(ang_v) - M_PI) < ang_tol));

        // By using the direction of the tangent in combination with the cw_loop/ccw_loop information we
        // should be able to handle points at a double seam.
        if (handle_u_seam && handle_v_seam)
	{   
            // We must compare the tangent with the surface tangent. Based on the direction and position on the curve (begin/end)
            // we can conclude on which side to project.
            if ((!cw_loop) && (!ccw_loop))
            {
                std::cout << "WARNING: The method branch expects the input to be part of a loop!" << std::endl;
            }

            if ((!u_parallel) && (!v_parallel))
            {
                // It is trivial to extend the method to support this case.
                std::cout << "WARNING: Double seam, non-tangential, case not handled!" << 
                    " ang_u: " << ang_u << ", ang_v: " << ang_v << std::endl;
                return shared_ptr<Point>(NULL);
            }
            else if (u_parallel)
            {   // We use the u-dir to pick the u param. We use the ccw/cw info to pick the v param.
                if (fabs(ang_u) < ang_tol)
                {
                    clo_u = (at_cv_start) ? umin : umax;
                    clo_v = (ccw_loop) ? vmin : vmax;
                }
                else if (fabs(fabs(ang_u) - M_PI) < ang_tol)
                {
                    clo_u = (at_cv_start) ? umax : umin;
                    clo_v = (ccw_loop) ? vmax : vmin;
                }
            }
            else // v_parallel
            {   // We use the u-dir to pick the v param. We use the ccw/cw info to pick the u param.
                if (fabs(ang_v) < ang_tol)
                {
                    clo_v = (at_cv_start) ? vmin : vmax;
                    clo_u = (ccw_loop) ? umax : umin;
                }
                else if (fabs(fabs(ang_v) - M_PI) < ang_tol)
                {
                    clo_v = (at_cv_start)? vmax : vmin;
                    clo_u = (ccw_loop) ? umin : umax;
                }
            }

            return shared_ptr<Point>(new Point(clo_u, clo_v));
        }

        // ElementarySurface* elem_sf = surface_->elementarySurface();
        // const double sign = (elem_sf == nullptr) ? 1.0 : (elem_sf->isSwapped() ? -1.0 : 1.0);
	else if (handle_u_seam)
	{
	    if (v_parallel)//(fabs(ang_v) < ang_tol) || ((fabs(ang_v) - M_PI) < ang_tol))// || ((fabs(ang_v - M_PI) < ang_tol)))
	    {   // We are following the seam, with no seed given, not handled currently.
		// @@sbr201506 We could use a marching approach to handle some cases.
//		MESSAGE("Following the seam! constdir_: " << constdir_ << ", constval_: " << constval_);
                if (cw_loop || ccw_loop)
                {
                    if (fabs(ang_v) < ang_tol)
                    {
//                    clo_v = (at_cv_start) ? vmin : vmax;
                        clo_u = (ccw_loop) ? umax : umin;
                    }
                    else if (fabs(fabs(ang_v) - M_PI) < ang_tol)
                    {
//                    clo_v = (at_cv_start)? vmax : vmin;
                        clo_u = (ccw_loop) ? umin : umax;
                    }
                }
                else
                {
                    bool march_left_success = false;
                    Point march_left_pt = Point(clo_u, clo_v);
                    marchOutSeamPoint(sf, space_cv, tpar, false, true, false, epsgeo,
                                      march_left_pt, march_left_success);
                    bool march_right_success = false;
                    Point march_right_pt = Point(clo_u, clo_v);
                    marchOutSeamPoint(sf, space_cv, tpar, true, true, false, epsgeo,
                                      march_right_pt, march_right_success);
                    if (march_left_success && march_right_success)
                    {
                        double dist = march_left_pt.dist(march_right_pt);
                        if (dist < knot_diff_tol)
                        {
                            // The result should be the same.
                            return shared_ptr<Point>(new Point(march_left_pt));
                        }
                        else
                        {
                            MESSAGE("Marching ended in mismatch.");
                        }
                    }
                    else if (march_left_success)
                    {
                        return shared_ptr<Point>(new Point(march_left_pt));
                    }
                    else if (march_right_success)
                    {
                        return shared_ptr<Point>(new Point(march_right_pt));
                    }
                    else
                    {
                        // This should mean that the whole curve is following the seam. We use loop orientation
                        // to choose side.
                        if (ccw_loop)
                        {
                            // If tangent is increasing we choose umax, otherwise umin.
                            clo_u = (par_tangent[1] > 0.0) ? umax : umin;
                            return shared_ptr<Point>(new Point(clo_u, clo_v));
                        }
                        else if (cw_loop)
                        {
                            clo_u = (par_tangent[1] > 0.0) ? umin : umax;
                            return shared_ptr<Point>(new Point(clo_u, clo_v));
                        }
//		    follows_seem_dir = 2;
                        MESSAGE("Marching failed. ccw_loop: " << ccw_loop << ", cw_loop: " << cw_loop);
                        return shared_ptr<Point>(NULL);
                    }
                }
	    }
	    else
	    {
		if (at_cv_end) // The end of the space curve.
		{
		    if (((ang_v < 0.0) && at_u_start) ||
			((ang_v > 0.0) && at_u_end))
		    {
			clo_u = (at_u_start) ? umax : umin;
		    }
		    return shared_ptr<Point>(new Point(clo_u, clo_v));
		}
		else if (at_cv_start)
		{
		    if (((ang_v < 0.0) && at_u_end) ||
			((ang_v > 0.0) && at_u_start))
		    {
			clo_u = (at_u_start) ? umax : umin;
		    }
		    return shared_ptr<Point>(new Point(clo_u, clo_v));
		}
		else
		{
		    MESSAGE("This routine does not handle curves crossing the seam!");
		    return shared_ptr<Point>(NULL);
		}
	    }
	}
	// else // at_v_bd
	else if (handle_v_seam) // Constant v parameter for the seam.
	{
	    if (u_parallel)//(fabs(ang_u) < ang_tol) || ((fabs(ang_u) - M_PI) < ang_tol))// || ((fabs(ang_u - M_PI) < ang_tol)))
	    { // We are along the seam.
//		MESSAGE("Following the seam! constdir_: " << constdir_ << ", constval_: " << constval_);
                if (cw_loop || ccw_loop)
                {
                    if (fabs(ang_u) < ang_tol)
                    {
//                    clo_u = (at_cv_start) ? umin : umax;
                        clo_v = (ccw_loop) ? vmin : vmax;
                    }
                    else if (fabs(fabs(ang_u) - M_PI) < ang_tol)
                    {
//                    clo_u = (at_cv_start) ? umax : umin;
                        clo_v = (ccw_loop) ? vmax : vmin;
                    }
                }
                else
                {
                    bool march_left_success = false;
                    Point march_left_pt = Point(clo_u, clo_v);
                    marchOutSeamPoint(sf, space_cv, tpar, false, false, true, epsgeo,
                                      march_left_pt, march_left_success);
                    bool march_right_success = false;
                    Point march_right_pt = Point(clo_u, clo_v);
                    marchOutSeamPoint(sf, space_cv, tpar, true, false, true, epsgeo,
                                      march_right_pt, march_right_success);
                    if (march_left_success && march_right_success)
                    {
                        double dist = march_left_pt.dist(march_right_pt);
                        if (dist < knot_diff_tol)
                        {
                            // The result should be the same.
                            return shared_ptr<Point>(new Point(march_left_pt));
                        }
                        else
                        {
                            MESSAGE("Marching ended in mismatch.");
                        }
                    }
                    else if (march_left_success)
                    {
                        return shared_ptr<Point>(new Point(march_left_pt));
                    }
                    else if (march_right_success)
                    {
                        return shared_ptr<Point>(new Point(march_right_pt));
                    }
                    else
                    {
//		    follows_seem_dir = 2;
                        // This should mean that the whole curve is following the seam. We use loop orientation
                        // to choose side.
                        if (ccw_loop)
                        {
                            // If tangent is increasing we choose umax, otherwise umin.
                            clo_v = (par_tangent[0] > 0.0) ? vmin : vmax;
                            return shared_ptr<Point>(new Point(clo_u, clo_v));
                        }
                        else if (cw_loop)
                        {
                            clo_v = (par_tangent[0] > 0.0) ? vmax : vmin;
                            return shared_ptr<Point>(new Point(clo_u, clo_v));
                        }
                        MESSAGE("Marching failed. ccw_loop: " << ccw_loop << ", cw_loop: " << cw_loop);
                        return shared_ptr<Point>(NULL);
                    }
                }
		// MESSAGE("Following the seam!");
		// return shared_ptr<Point>(NULL);
	    }
	    else
	    {
		if (at_cv_end) // The end of the space curve.
		{
		    if (((ang_u > 0.0) && at_v_start) ||
			((ang_u < 0.0) && at_v_end))
		    {
			clo_v = (at_v_start) ? vmax : vmin;
		    }
		    return shared_ptr<Point>(new Point(clo_u, clo_v));
		}
		else if (at_cv_start)
		{
		    if (((ang_u > 0.0) && at_v_end) ||
			((ang_u < 0.0) && at_v_start))
		    {
			clo_v = (at_v_start) ? vmax : vmin;
		    }
		    return shared_ptr<Point>(new Point(clo_u, clo_v));
		}
		else
		{
		    MESSAGE("This routine does not handle curves crossing the seam!");
		    return shared_ptr<Point>(NULL);
		}
	    }
	}
    }

    return shared_ptr<Point>(new Point(clo_u, clo_v));
}


//===========================================================================
Point BoundedUtils::projectSpaceCurveTangent(const ParamSurface& sf, const ParamCurve& space_cv,
                                             const Point& par_pt, double tpar)
//===========================================================================
{
    vector<Point> space_cv_pt = space_cv.point(tpar, 1);
    const int dim = space_cv.dimension();
    // We then project the tangent into the parameter domain.
    vector<Point> sf_pt = sf.ParamSurface::point(par_pt[0], par_pt[1], 1);
    // We describe the curve tangent as a linear combination of the partial derivs.
    double coef1, coef2;
    blendcoef(&sf_pt[1][0], &sf_pt[2][0],
	      &space_cv_pt[1][0], dim, 1, &coef1, &coef2);

    Point pt_dir(2);
    pt_dir[0] = coef1;
    pt_dir[1] = coef2;

    return pt_dir;
}

//===========================================================================
  void BoundedUtils::consistentIntersectionDir(ParamCurve& inters_pcv,
					       ParamCurve& inters_space_cv,
					       const ParamSurface& sf,
					       const ParamCurve& other_inters_pcv,
					       const ParamCurve& other_inters_space_cv,
					       const ParamSurface& other_sf,
					       double epsgeo)
//===========================================================================
{
    Point par_pt = inters_pcv.point(inters_pcv.startparam());
    Point sf_normal;
    sf.normal(sf_normal, par_pt[0], par_pt[1]);
    vector<Point> space_pt = inters_space_cv.point(inters_space_cv.startparam(), 1);
    Point tangent = space_pt[1];
    tangent.normalize_checked();
    Point other_par_pt1 = other_inters_pcv.point(other_inters_pcv.startparam());
    Point other_par_pt2 = other_inters_pcv.point(other_inters_pcv.endparam());
    vector<Point> other_space_pt1 =
	other_inters_space_cv.point(other_inters_space_cv.startparam(), 1);
    vector<Point> other_space_pt2 =
	other_inters_space_cv.point(other_inters_space_cv.endparam(), 1);
    double dist1 = other_space_pt1[0].dist(space_pt[0]);
    double dist2 = other_space_pt2[0].dist(space_pt[0]);
    Point other_sf_normal;
    if ((dist1 < epsgeo) && (dist2 < epsgeo)) {
	// Both end pts of other_inters_pcv match that of inters_pcv, we must choose based on
	// tangents of space curve.
	tangent.normalize_checked();
	Point other_tangent1 = other_space_pt1[1];
	other_tangent1.normalize_checked();
	Point other_tangent2 = other_space_pt2[1];
	other_tangent2.normalize_checked();
	double min_dist1 = std::min(tangent.dist(other_tangent1), tangent.dist(-other_tangent1));
	double min_dist2 = std::min(tangent.dist(other_tangent2), tangent.dist(-other_tangent2));
	if (min_dist1 < min_dist2) {
	    other_sf.normal(other_sf_normal, other_par_pt1[0], other_par_pt1[1]);
	} else {
	    other_sf.normal(other_sf_normal, other_par_pt2[0], other_par_pt2[1]);
	}
    } else if (dist1 < epsgeo) {
	other_sf.normal(other_sf_normal, other_par_pt1[0], other_par_pt1[1]);
    } else if (dist2 < epsgeo) {
	other_sf.normal(other_sf_normal, other_par_pt2[0], other_par_pt2[1]);
    } else {
	THROW("Input cvs must match in end pts!");
    }

    // Direction of inters_space_cv (in start par) should be that of other_sf_normal X sf_normal.
    Point cross_prod;
    cross_prod.setToCrossProd(other_sf_normal, sf_normal);
    cross_prod.normalize_checked();
    if (cross_prod.dist(-tangent) < cross_prod.dist(tangent)) {
	inters_pcv.reverseParameterDirection();
	inters_space_cv.reverseParameterDirection();
    }
}

//===========================================================================
  void 
  BoundedUtils::consistentIntersectionDir(shared_ptr<CurveOnSurface> sf_cv1,
					  shared_ptr<CurveOnSurface> sf_cv2,
					  double epsgeo)
//===========================================================================
{
  if (!sf_cv1->hasParameterCurve())
    return; // Nothing to do
  Point par_pt = sf_cv1->parameterCurve()->point(sf_cv1->startparam());
  Point sf_normal;
  sf_cv1->underlyingSurface()->normal(sf_normal, par_pt[0], par_pt[1]);
  vector<Point> space_pt = sf_cv1->spaceCurve()->point(sf_cv1->startparam(), 1);
  Point tangent = space_pt[1];
  tangent.normalize_checked();
  if (!sf_cv2->hasParameterCurve())
    return; // Nothing to do
  Point other_par_pt1 = sf_cv2->parameterCurve()->point(sf_cv2->startparam());
  Point other_par_pt2 = sf_cv2->parameterCurve()->point(sf_cv2->endparam());
  vector<Point> other_space_pt1 = 
    sf_cv2->spaceCurve()->point(sf_cv2->startparam(), 1);
  vector<Point> other_space_pt2 = 
    sf_cv2->spaceCurve()->point(sf_cv2->endparam(), 1);
  double dist1 = other_space_pt1[0].dist(space_pt[0]);
  double dist2 = other_space_pt2[0].dist(space_pt[0]);
  Point other_sf_normal;
  if ((dist1 < epsgeo) && (dist2 < epsgeo)) {
    // Both end pts of other_inters_pcv match that of inters_pcv, we must choose based on
    // tangents of space curve.
    tangent.normalize_checked();
    Point other_tangent1 = other_space_pt1[1];
    other_tangent1.normalize_checked();
    Point other_tangent2 = other_space_pt2[1];
    other_tangent2.normalize_checked();
    double min_dist1 = std::min(tangent.dist(other_tangent1), tangent.dist(-other_tangent1));
    double min_dist2 = std::min(tangent.dist(other_tangent2), tangent.dist(-other_tangent2));
    if (min_dist1 < min_dist2) {
      sf_cv2->underlyingSurface()->normal(other_sf_normal, other_par_pt1[0], other_par_pt1[1]);
    } else {
      sf_cv2->underlyingSurface()->normal(other_sf_normal, other_par_pt2[0], other_par_pt2[1]);
    }
  } else if (dist1 < epsgeo) {
    sf_cv2->underlyingSurface()->normal(other_sf_normal, other_par_pt1[0], other_par_pt1[1]);
  } else if (dist2 < epsgeo) {
    sf_cv2->underlyingSurface()->normal(other_sf_normal, other_par_pt2[0], other_par_pt2[1]);
  } else {
    THROW("Input cvs must match in end pts!");
  }

  // Direction of inters_space_cv (in start par) should be that of other_sf_normal X sf_normal.
  Point cross_prod;
  cross_prod.setToCrossProd(other_sf_normal, sf_normal);
  cross_prod.normalize_checked();
  if (cross_prod.dist(-tangent) < cross_prod.dist(tangent)) {
    sf_cv1->reverseParameterDirection();
  }
}

//===========================================================================
  double BoundedUtils::getParEps(double space_eps, const ParamSurface *sf)
//===========================================================================
{
    // Compute average length of boundary loop curves and length of
    // (diagonal of) parameter domain. The parameter eps is then
    // the space eps times the ratio of these lengths.

    int nmb_samples = 10;
    CurveLoop crvs = sf->outerBoundaryLoop();
    int nmb_crvs = crvs.size();

    // Average length of curves in space
    double space_len = 0.0;
    for (int i = 0; i < nmb_crvs; ++i)
	space_len += crvs[i]->estimatedCurveLength(nmb_samples);
    space_len /= double(nmb_crvs);

    // Average length of parameter domain
    RectDomain dom = sf->containingDomain();
    double par_len = (dom.lowerLeft()).dist(dom.upperRight());

    double par_eps = space_eps * par_len / space_len;

    // Control for unbalanced parameter domains
    double fac = 0.01;
    par_eps = std::min(par_eps, fac*(dom.umax()-dom.umin()));
    par_eps = std::min(par_eps, fac*(dom.vmax()-dom.vmin()));

    return par_eps;
}



} // end namespace Go


namespace {

//===========================================================================
  double getSeed(Point space_pt, CurveOnSurface& cv_on_sf, bool par_cv)
//===========================================================================
{
   shared_ptr<SplineCurve> spline_cv;
   if (cv_on_sf.parPref() || par_cv) {
      spline_cv =
	 dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf.parameterCurve());
   } else {
      spline_cv =
	 dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf.spaceCurve());
   }

   int nmb_samples;
   if (spline_cv.get())
     nmb_samples = std::max(3, spline_cv->numCoefs());
   else
     nmb_samples = 10;

   double tmin = cv_on_sf.startparam();
   double tmax = cv_on_sf.endparam();
   double tstep = (tmax - tmin)/(nmb_samples - 1);
   std::vector<std::pair<double, double> > par_and_dist;
   par_and_dist.reserve(nmb_samples);
   int ki;
   for (ki = 0; ki < nmb_samples; ++ki) {
      double tpar = tmin + ki*tstep;
      Point cv_pt;
      if (par_cv)
	cv_pt = cv_on_sf.parameterCurve()->point(tpar);
      else
	cv_on_sf.point(cv_pt, tpar);
      double dist = space_pt.dist(cv_pt);
      par_and_dist.push_back(std::make_pair(tpar, dist));
   }

   std::sort(par_and_dist.begin(), par_and_dist.end(),
	     comparepair_second< std::pair<double, double> >());
   double guess_param = par_and_dist[0].first;

   return guess_param;
}



//===========================================================================
int leftMostCurve(CurveOnSurface& cv,
		  vector<shared_ptr<CurveOnSurface> >& other_cvs,
		  double space_eps, bool par_cv, double& angle)
//===========================================================================
{
    // A closed surface needs special care, we use cvs starting in the same pt in parametric domain.
    // Tolerance may be somewhat larger than the corresponding spatial value as it is solely used
    // to handle closed surfaces.
  double a_tol = 1.0e-10;  // To include equality in angular test
    RectDomain domain = cv.underlyingSurface()->containingDomain();
    double par_closed_eps = 0.5*std::min(domain.umax() - domain.umin(), domain.vmax() - domain.vmin());

    ALWAYS_ERROR_IF(cv.parameterCurve().get() == 0, "could not find parameter curve");
    double endpar = cv.endparam();
    vector<Point> par_end_pt =
	cv.parameterCurve()->ParamCurve::point(endpar, 1);
    Point end_pt = par_cv ? par_end_pt[0] : cv.ParamCurve::point(endpar);
    int left_most_ind = -1;
    int zero_ind = -1;
    bool opposite = false;
    double min_angle = 8.0; // More than 2 pi
    const double ANG_TOL = 1e-06; // We do not want to return edge returning in the same direction.
    for (int ki = 0; ki < int(other_cvs.size()); ++ki) {
	ALWAYS_ERROR_IF(other_cvs[ki]->parameterCurve().get() == 0, "missing parameter curve");
	double startpar = other_cvs[ki]->startparam();
	vector<Point> par_start_pt = 
	  other_cvs[ki]->parameterCurve()->ParamCurve::point(startpar, 1);
	Point start_pt = par_cv ? par_start_pt[0] :
	  other_cvs[ki]->ParamCurve::point(startpar);

	if ((end_pt.dist(start_pt) < space_eps) &&
	    (par_end_pt[0].dist(par_start_pt[0]) < par_closed_eps)) {
	    // We now must compute the angle.
	    double local_angle = par_end_pt[1].angle(par_start_pt[1]);
	    Point vec2 = Point(par_end_pt[1][1],-par_end_pt[1][0]);

	    if (vec2*par_start_pt[1] < -ANG_TOL) 
	      local_angle *= -1.0; 
	    local_angle += M_PI;

	    if (fabs(local_angle - min_angle) < a_tol)
	      {
		// A tangential situation. Recompute angles
		recomputeAngles( par_end_pt, ANG_TOL, other_cvs, 
				left_most_ind, min_angle, ki, local_angle);
	      }

	    if (fabs(vec2*par_start_pt[1]) <= ANG_TOL && zero_ind >= 0 &&
		!(opposite && par_end_pt[1]*par_start_pt[1] >= 0.0))
	      {
		// A tangential situation. Recompute angles
		double ang2 = 2.0*M_PI;
		recomputeAngles( par_end_pt, ANG_TOL, other_cvs, 
				zero_ind, ang2, ki, local_angle);
		if (ang2 < min_angle)
		  {
		    min_angle = ang2;
		    left_most_ind = zero_ind;
		    zero_ind = -1;
		  }
		if (local_angle < min_angle+a_tol)
		  zero_ind = -1;
	      }
	    if (fabs(vec2*par_start_pt[1]) <= ANG_TOL && (zero_ind < 0 || opposite))
	      {
		local_angle = 2.0*M_PI;
		zero_ind = ki;
		opposite =  (par_end_pt[1]*par_start_pt[1] < 0.0);
	      }

	    if (local_angle < min_angle+a_tol) {
		min_angle = local_angle;
		left_most_ind = ki;
	    }
	}
    }

    if (min_angle > M_PI && (!opposite) &&  zero_ind >= 0)
      {
	min_angle = 0.0;
	left_most_ind = zero_ind;
      }

    angle = min_angle;
    return left_most_ind;
}

//===========================================================================
  void recomputeAngles(vector<Point>& curr_val, double ang_tol,
		       vector<shared_ptr<CurveOnSurface> >& cvs,
		       int idx1, double& ang1, int idx2, double& ang2)
//===========================================================================
  {
    // (Approximate) tangential situation. Recompute angles based on
    // less local information to get a more reliable decision
    // NB! The current choice of evaluation point is quite arbitrary
    // and can be improved (VSK, 201004)
    double t1 = cvs[idx1]->startparam() + 0.1*(cvs[idx1]->endparam()  - 
					       cvs[idx1]->startparam());
    double t2 = cvs[idx2]->startparam() + 0.1*(cvs[idx2]->endparam()  - 
	   cvs[idx2]->startparam());
    Point pos1_1 = cvs[idx1]->parameterCurve()->ParamCurve::point(cvs[idx1]->startparam());
    Point pos2_1 = cvs[idx2]->parameterCurve()->ParamCurve::point(cvs[idx2]->startparam());
    Point pos1_2 = cvs[idx1]->parameterCurve()->ParamCurve::point(t1);
    Point pos2_2 = cvs[idx2]->parameterCurve()->ParamCurve::point(t2);
    Point vec = Point(curr_val[1][1], -curr_val[1][0]);
    Point vec1 = pos1_2 - pos1_1;
    Point vec2 = pos2_2 - pos2_1;
    ang1 = curr_val[1].angle(vec1);
    if (vec*vec1 < -ang_tol)
      ang1 *= -1;
    ang1 += M_PI;

    ang2 = curr_val[1].angle(vec2);
    if (vec*vec2 < -ang_tol)
      ang2 *= -1;
    ang2 += M_PI;
  }



//===========================================================================
void marchOutSeamPoint(const ParamSurface& surface, const ParamCurve& space_cv,
                       double tpar, bool to_the_right, bool at_u_seam, bool at_v_seam,
                       double epsgeo, Point& par_pt, bool& success)
//===========================================================================
{
    const Point sf_epspar = SurfaceTools::getParEpsilon(surface, epsgeo);
    const double epspar = std::min(sf_epspar[0], sf_epspar[1]);
    const double knot_diff_tol = epspar;//1e-08;
    Point space_pt = space_cv.ParamCurve::point(tpar);
    double tmin = space_cv.startparam();
    double tmax = space_cv.endparam();

    RectDomain rect_dom = surface.containingDomain();

    if (at_u_seam && at_v_seam)
    {   // We must extend our method to keep on going until ok_left_u and ok_left_v are ok etc.
	MESSAGE("Function is not really prepared for a point in a torus corner ...");
    }

    const double march_limit = (to_the_right) ? tmax : tmin;

    if (fabs(march_limit - tpar) < knot_diff_tol)
    {   // We are at the end going in that direction.
	success = false;
	return;
    }

//    Point par_pt_march = par_pt;
    double range = fabs(tpar - march_limit);
    int sign = (to_the_right) ? 1 : -1;
    double tstep_frac = 1e-04; // 14 steps to get past the start param.
    int mult = 1;
    double step_tpar = tpar + sign*tstep_frac*mult*range;
    Point seed_pt = par_pt;
    while ((step_tpar > tmin) && (step_tpar < tmax))
    {
	double clo_u, clo_v, clo_dist;
	Point clo_pt;
        Point step_space_pt = space_cv.ParamCurve::point(step_tpar);
	surface.closestPoint(step_space_pt, clo_u, clo_v, clo_pt, clo_dist, epsgeo, NULL, &seed_pt[0]);
	if ((at_u_seam) && (fabs(par_pt[0] - clo_u) > knot_diff_tol))
	{
	    par_pt[0] = (clo_u < 0.5*(rect_dom.umin() + rect_dom.umax())) ? rect_dom.umin() : rect_dom.umax();
	    success = true;
	}
	if ((at_v_seam) && (fabs(par_pt[1] - clo_v) > knot_diff_tol))
	{
	    par_pt[1] = (clo_v < 0.5*(rect_dom.vmin() + rect_dom.vmax())) ? rect_dom.vmin() : rect_dom.vmax();
	    success = true;
	}
	if (success)
	{
	    break;
	}
	else
	{
	    seed_pt = Point(clo_u, clo_v);
	    mult *= 2;
	    step_tpar = tpar + sign*tstep_frac*mult*range;
	}
    }
}

//===========================================================================
  double maxCurveGap(shared_ptr<CurveOnSurface> curves[], int nmb_curves,
		     bool par_crv, double eps, vector<Point>& free_endpt)
//===========================================================================
{
  if (nmb_curves == 0)
    return 0.0;

  double curve_gap = 0.0;
  for (int ki=0; ki<nmb_curves; ++ki)
    {
      double t1 = curves[ki]->startparam();
      double t2 = curves[ki]->endparam();
      Point pos1 = par_crv ? curves[ki]->parameterCurve()->point(t1) :
	curves[ki]->ParamCurve::point(t1);
      Point pos2 = par_crv ? curves[ki]->parameterCurve()->point(t2) :
	curves[ki]->ParamCurve::point(t2);
      double mindist = std::numeric_limits<double>::max();
      double min1 = std::numeric_limits<double>::max(); 
      double min2 = std::numeric_limits<double>::max();
      if (nmb_curves < 2)
	{
	  min1 = min2 = mindist = 0.0;
	}
      else
	{
	  for (int kj=0; kj<nmb_curves; ++kj)
	    {
	      if (ki == kj)
		continue;
	      double t3 = curves[kj]->startparam();
	      double t4 = curves[kj]->endparam();
	      Point pos3 = par_crv ? curves[kj]->parameterCurve()->point(t3) :
		curves[kj]->ParamCurve::point(t3);
	      Point pos4 = par_crv ? curves[kj]->parameterCurve()->point(t4) :
		curves[kj]->ParamCurve::point(t4);
	      double dd1 = std::min(pos1.dist(pos3), pos1.dist(pos4));
	      double dd2 = std::min(pos2.dist(pos3), pos2.dist(pos4));
	      double dist = std::min(dd1, dd2);
	      mindist = std::min(mindist, dist);
	      min1 = std::min(min1, dd1);
	      min2 = std::min(min2, dd2);
	    }
	  if (mindist <= eps)
	    curve_gap = std::max(curve_gap, mindist);
	  if (min1 > eps)
	    free_endpt.push_back(pos1);
	  if (min2 > eps)
	    free_endpt.push_back(pos2);
	}
    }
      return curve_gap;
  }

//===========================================================================
  void splitLoopCvs(const BoundedSurface& sf,
		    vector<shared_ptr<CurveOnSurface> >& old_loop_cvs,
		    vector<shared_ptr<CurveOnSurface> >& part_bd_cvs,
		    vector<Point>& part_bd_endpt, double min_loop_tol,
		    double eps, double epspar, double knot_diff_tol,
		    int last_split, bool par_cv)
//===========================================================================
{
  double a_tol = 1.0e-8;
  int ki, kj;
  double epspar2 = (par_cv) ? min_loop_tol : epspar;
  if (last_split < 0)
    last_split = (int)part_bd_cvs.size();
  for (ki = 0; ki < last_split; ++ki) {
    if (!part_bd_cvs[ki]->parameterCurve().get()) {
      // Suppose we should set curve to prefer parametric part.
      THROW("Only space curves, method uses parameter curves...");
    }
    double t1 = part_bd_cvs[ki]->startparam();
    double t2 = part_bd_cvs[ki]->endparam();
    Point par_start_pt =  part_bd_cvs[ki]->parameterCurve()->point(t1);
    Point space_start_pt = par_cv ? par_start_pt :
      part_bd_cvs[ki]->ParamCurve::point(t1);
    Point par_end_pt =  part_bd_cvs[ki]->parameterCurve()->point(t2);
    Point space_end_pt = par_cv ? par_end_pt :
      part_bd_cvs[ki]->ParamCurve::point(t2);

    // Check if a part boundary curve endpoint is a loose end that should be connected
    // to an old loop curve
    bool loose_end1 = false, loose_end2 = false;
    for (size_t kr=0; kr<part_bd_endpt.size(); ++kr)
      {
	if (space_start_pt.dist(part_bd_endpt[kr]) < min_loop_tol)
	  loose_end1 = true;
	if (space_end_pt.dist(part_bd_endpt[kr]) < min_loop_tol)
	  loose_end2 = true;
      }

    for (kj = 0; kj < int(old_loop_cvs.size()); ++kj) {

      if (!old_loop_cvs[kj]->parameterCurve().get())
	// Suppose we should set curve to prefer parametric part.
	THROW("Only space curves, method uses parameter curves...");
      
      double t3 = old_loop_cvs[kj]->startparam();
      double t4 = old_loop_cvs[kj]->endparam();
      Point par_start = 
	old_loop_cvs[kj]->parameterCurve()->point(t3);
      Point par_end = 
	old_loop_cvs[kj]->parameterCurve()->point(t4);
      Point geom_start = par_cv ? par_start :
	old_loop_cvs[kj]->ParamCurve::point(t3);
      Point geom_end = par_cv ? par_end :
	old_loop_cvs[kj]->ParamCurve::point(t4);
      shared_ptr<ParamCurve> pcv = old_loop_cvs[kj]->parameterCurve();
      if (pcv->endparam() - pcv->startparam() < knot_diff_tol) {
	// Check length of curve in geometry space and parameter space
	if (geom_start.dist(geom_end) < min_loop_tol && 
	    par_start.dist(par_end) < epspar2)
	  {
	    MESSAGE("Loop segment smaller than loop tolerance, removing segment.");
	    old_loop_cvs.erase(old_loop_cvs.begin() + kj);
	    --kj;
	    continue;
	  }
      }

      // We compute distance from start and end pt of part_bd_cv.
      double start_t, end_t, clo_dist_start, clo_dist_end;
      Point clo_start_pt, clo_end_pt;
      double seed = getSeed(space_start_pt, *old_loop_cvs[kj], par_cv);
      pcv->closestPoint(par_start_pt, pcv->startparam(), pcv->endparam(),
			start_t, clo_start_pt, clo_dist_start, &seed);

      // Check against total_start_pt.
      seed = getSeed(space_end_pt, *old_loop_cvs[kj], par_cv);
      pcv->closestPoint(par_end_pt, pcv->startparam(), pcv->endparam(),
			end_t, clo_end_pt, clo_dist_end, &seed);

      double space_start_dist, space_end_dist;
      if (!par_cv)
	{
	  // Modify endpoint distances with respect to space curve information
	  Point clo_start_pt_space = 
	    sf.ParamSurface::point(clo_start_pt[0], clo_start_pt[1]);
	  Point clo_end_pt_space = 
	    sf.ParamSurface::point(clo_end_pt[0], clo_end_pt[1]);

	  double par_close_start, par_close_end, dist_close_start, dist_close_end;
	  Point point_close_start, point_close_end;
	  old_loop_cvs[kj]->closestPoint(space_start_pt, 
					 old_loop_cvs[kj]->startparam(),
					 old_loop_cvs[kj]->endparam(),
					 par_close_start, point_close_start,
					 dist_close_start, &start_t);
      
	  old_loop_cvs[kj]->closestPoint(space_end_pt, 
					 old_loop_cvs[kj]->startparam(),
					 old_loop_cvs[kj]->endparam(),
					 par_close_end, point_close_end,
					 dist_close_end, &end_t);

	  space_start_dist = clo_start_pt_space.dist(space_start_pt);
	  space_end_dist = clo_end_pt_space.dist(space_end_pt);

	  if (std::min(space_start_dist, dist_close_start) < min_loop_tol)
	    {
	      space_start_dist = dist_close_start;
	      start_t = par_close_start;
	      min_loop_tol = std::max(min_loop_tol, dist_close_start+a_tol);
	    }
	  if (std::min(space_end_dist, dist_close_end) < min_loop_tol)
	    {
	      space_end_dist = dist_close_end;
	      end_t = par_close_end;
	      min_loop_tol = std::max(min_loop_tol, dist_close_end+a_tol);
	    }
	}
      else
	{
	  space_start_dist = clo_dist_start;
	  space_end_dist = clo_dist_end;
	}

      if (space_start_dist < min_loop_tol)
	{
	  loose_end1 = false;

	  // Check if a loose endpoint of a part boundary curve	
	  // is coincident with the start point
	  for (size_t kr=0; kr<part_bd_endpt.size();)
	    {
	      if (space_start_pt.dist(part_bd_endpt[kr]) < min_loop_tol)
		part_bd_endpt.erase(part_bd_endpt.begin()+kr);
	      else
		++kr;
	    }
	}

      if (space_end_dist < min_loop_tol)
	{
	  loose_end2 = false;
	  for (size_t kr=0; kr<part_bd_endpt.size();)
	    {
	      if (space_end_pt.dist(part_bd_endpt[kr]) < min_loop_tol)
		part_bd_endpt.erase(part_bd_endpt.begin()+kr);
	      else
		++kr;
	    }
	}

      bool is_split = false;
      if (space_start_dist < min_loop_tol  &&
	  ((start_t - knot_diff_tol > old_loop_cvs[kj]->startparam() &&
	    start_t + knot_diff_tol < old_loop_cvs[kj]->endparam()) ||
	   (geom_start.dist(space_start_pt) > min_loop_tol &&
	    geom_end.dist(space_start_pt) > min_loop_tol)))
	{
	  is_split = true;
	  vector<shared_ptr<ParamCurve> > sub_cvs = 
	    old_loop_cvs[kj]->split(start_t);
	  for (size_t k2=0; k2<sub_cvs.size(); ++k2)
	    {
	      shared_ptr<CurveOnSurface> sf_cv = 
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[k2]);
	      if (sf_cv.get())
		sf_cv->ensureParCrvExistence(eps);
	    }
	  shared_ptr<CurveOnSurface> sub1 =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[0]);
	  shared_ptr<CurveOnSurface> sub2 =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[1]);
	  double len1 = par_cv ? 
	    sub1->parameterCurve()->estimatedCurveLength(3) : 
	    sub1->estimatedCurveLength(3);
	  double len2 = par_cv ? 
	    sub2->parameterCurve()->estimatedCurveLength(3) : 
	    sub2->estimatedCurveLength(3);
	  if (len1 > min_loop_tol || len1 > len2)
	    old_loop_cvs[kj] = sub1;
	  else
	    old_loop_cvs[kj] = sub2;

	  if (len1 > min_loop_tol && len2 > min_loop_tol)
	    {
	      old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub2);
	      if (space_end_dist < min_loop_tol &&
		  old_loop_cvs[kj]->startparam() < end_t &&
		  old_loop_cvs[kj]->endparam() > end_t)
		std::swap(old_loop_cvs[kj], old_loop_cvs[kj+1]);
	      kj++;
	    }
	}
      if (space_end_dist < min_loop_tol &&
	  ((end_t - knot_diff_tol > old_loop_cvs[kj]->startparam() &&
	    end_t + knot_diff_tol < old_loop_cvs[kj]->endparam())||
	   (geom_start.dist(space_end_pt) > min_loop_tol &&
	    geom_end.dist(space_end_pt) > min_loop_tol && (!is_split))))
	{
	  vector<shared_ptr<ParamCurve> > sub_cvs = 
	    old_loop_cvs[kj]->split(end_t);
	  for (size_t k2=0; k2<sub_cvs.size(); ++k2)
	    {
	      shared_ptr<CurveOnSurface> sf_cv = 
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[k2]);
	      if (sf_cv.get())
		sf_cv->ensureParCrvExistence(eps);
	    }
	  shared_ptr<CurveOnSurface> sub1 =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[0]);
	  shared_ptr<CurveOnSurface> sub2 =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[1]);
	  double len1 = par_cv ? 
	    sub1->parameterCurve()->estimatedCurveLength(3) : 
	    sub1->estimatedCurveLength(3);
	  double len2 = par_cv ? 
	    sub2->parameterCurve()->estimatedCurveLength(3) : 
	    sub2->estimatedCurveLength(3);

	  if (len1 > min_loop_tol || len1 > len2)
	    old_loop_cvs[kj] = sub1;
	  else
	    old_loop_cvs[kj] = sub2; 

	  if (len1 > min_loop_tol && len2 > min_loop_tol)
	    {
	      old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub2);
	      kj++;
	    }
	}
    }
    double min_loop_fac = 10.0;
    vector<Point> loose_endpt;
    if (loose_end1)
      loose_endpt.push_back(space_start_pt);
    if (loose_end2)
      loose_endpt.push_back(space_end_pt);

    for (size_t kr=0; kr<loose_endpt.size(); ++kr)
      {
	// Part boundary curve ends in nothing. Find closest point on old boundary
	// curves
	int ix = -1;
	double mindist = std::numeric_limits<double>::max();
	double loop_par;
	for (kj = 0; kj < int(old_loop_cvs.size()); ++kj) 
	  {
	    double t1 = old_loop_cvs[kj]->startparam();
	    double t2 = old_loop_cvs[kj]->endparam();
	    double par, dist;
	    Point close;
	    if (par_cv)
	      old_loop_cvs[kj]->parameterCurve()->closestPoint(loose_endpt[kr],
							       t1, t2, par,
							       close, dist);
	    else
	      old_loop_cvs[kj]->closestPoint(loose_endpt[kr], t1, t2,
					     par, close, dist);
	    if (dist < mindist)
	      {
		mindist = dist;
		ix = kj;
		loop_par = par;
	      }
	  }
	if (mindist < min_loop_fac*min_loop_tol)
	  {
	    // Split curve anyway and adjust tolerance
	    if (loop_par > old_loop_cvs[ix]->startparam()+knot_diff_tol &&
		loop_par < old_loop_cvs[ix]->endparam()-knot_diff_tol)
	      {
		vector<shared_ptr<ParamCurve> > sub_cvs = 
		  old_loop_cvs[ix]->split(loop_par);
		for (size_t k2=0; k2<sub_cvs.size(); ++k2)
		  {
		    shared_ptr<CurveOnSurface> sf_cv = 
		      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[k2]);
		    if (sf_cv.get())
		      sf_cv->ensureParCrvExistence(eps);
		  }

		shared_ptr<CurveOnSurface> sub1 =
		  dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[0]);
		shared_ptr<CurveOnSurface> sub2 =
		  dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[1]);
		double len1 = par_cv ? 
		  sub1->parameterCurve()->estimatedCurveLength(3) : 
		  sub1->estimatedCurveLength(3);
		double len2 = par_cv ? 
		  sub2->parameterCurve()->estimatedCurveLength(3) : 
		  sub2->estimatedCurveLength(3);

		if (len1 > min_loop_tol)
		  old_loop_cvs[ix] = sub1;
		else
		  old_loop_cvs[ix] = sub2;

		if (len1 > min_loop_tol && len2 > min_loop_tol)
		  {
		    old_loop_cvs.insert(old_loop_cvs.begin() + ix, sub2);
		  }
	      }
	    min_loop_tol = std::max(min_loop_tol,mindist+a_tol);
	  }
      }
  }
}

}; // end anonymous namespace
