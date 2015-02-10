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

#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>
#include <utility>
#include "GoTools/geometry/SISL_code.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/utils/Values.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ClosestPoint.h"


using namespace Go;
using std::vector;
using std::pair;
using std::swap;

namespace {
// Given input of curves which may start in end pt of curr_cv (within space_eps), return
// index of the one with the greatest angle (in the interval (-pi, pi]).
// If there is no such curve, return -1.
int leftMostCurve(CurveOnSurface& cv,
		  vector<shared_ptr<CurveOnSurface> >& other_cvs,
		  double space_eps, double& angle);

  // Used by leftMostCurve in tangential situations
void recomputeAngles(vector<Point>& curr_val, double ang_tol,
		       vector<shared_ptr<CurveOnSurface> >& cvs,
		     int idx1, double& ang1, int idx2, double& ang2);

// Given input of intersection curves, we make sure orientation of intersection curve is in
// correpondence with orientation of surfaces (part below other_sf is to be trimmed away).
void consistentIntersectionDir(ParamCurve& inters_pcv,
			       ParamCurve& inters_space_cv,
			       const Go::ParamSurface& sf,
			       const ParamCurve& other_inters_pcv,
			       const ParamCurve& other_inters_space_cv,
			       const Go::ParamSurface& other_sf,
			       double epsgeo);

double getSeed(Point space_pt, CurveOnSurface& cv_on_sf);

// Based on parameter domain and lengths of sf's edges, return corresponding par_eps.
  double getParEps(double space_eps, const ParamSurface *sf);

}; // end anonymous namespace 

namespace Go {

//===========================================================================
vector<shared_ptr<CurveOnSurface> >
BoundedUtils::intersectWithSurface(CurveOnSurface& curve,
				     BoundedSurface& bounded_surf,
				     double epsge)
//===========================================================================
{
  double int_tol = 0.1*epsge; //1e-06;

    // We extract boundary loop, and check for intersections.
    // We do not handle trimming of trimmed surfaces with holes.
    vector<CurveLoop> boundary_loops = bounded_surf.allBoundaryLoops(epsge);

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
    // if (!curve.parPref()) {
    // 	curve = CurveOnSurface(curve.underlyingSurface(), curve.parameterCurve(),
    // 			       curve.spaceCurve(), true);
    // }
    for (i = 0; i < int(loop_curves.size()); ++i) {
      loop_curves[i]->ensureParCrvExistence(epsge);
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
    ALWAYS_ERROR_IF(first_curve.get() == 0,
		    "Intersection routine not implemented for general curves.");

    for (j = 0; j < int(loop_curves.size()); ++j) {
	shared_ptr<SplineCurve> second_curve =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(loop_curves[j]->parameterCurve());
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
	    ClosestPoint::closestPtCurves(&curve, loop_curves[j].get(), curve.startparam(),
			    curve.endparam(), loop_curves[j]->startparam(),
			    loop_curves[j]->endparam(),
			    int_params[kk].first, int_params[kk].second,
			    par1, par2, dist, ptc1, ptc2);
	  
	    if (dist < init_dist)
	      all_int_params.push_back(par1);
	    else
	      all_int_params.push_back(int_params[kk].first);
	}
    }

    // We sort the intersetcion parameters in ascending order.
    sort(all_int_params.begin(), all_int_params.end());
    if ((all_int_params.size() == 0) || (all_int_params.front() != first_curve->startparam())) {
    	all_int_params.insert(all_int_params.begin(), first_curve->startparam());
    }
    if (all_int_params.back() != first_curve->endparam()) {
    	all_int_params.insert(all_int_params.end(), first_curve->endparam());
    }

    // We extract parts of curve inside trimmed domain.
    vector<shared_ptr<CurveOnSurface> > inside_segments;
    shared_ptr<const ParamSurface> under_sf = bounded_surf.underlyingSurface();
    double knot_diff_tol = 0.01*getParEps(epsge, under_sf.get()); // We may not trust pcv to repr space_cv.
    for (j = 0; j < int(all_int_params.size()) - 1; ++j) {
	double from_par = all_int_params[j];
	double to_par = all_int_params[j+1];
	Point tmp1 = first_curve->ParamCurve::point(from_par);
	Point tmp2 = first_curve->ParamCurve::point(0.5*(from_par+to_par));
	Point tmp3 = first_curve->ParamCurve::point(to_par);
	double len = tmp1.dist(tmp2) + tmp2.dist(tmp3);
	if (to_par - from_par < knot_diff_tol || len < epsge) {
	    all_int_params.erase(all_int_params.begin() + j + 1);
	    --j;
	    continue;
	}
	if (from_par < curve.startparam())
	  from_par = curve.startparam();
	if (to_par > curve.endparam())
	  to_par = curve.endparam();
	double med_par = 0.5*(from_par + to_par);
	Point med_pt = first_curve->ParamCurve::point(med_par);
	if (domain.isInDomain(Vector2D(med_pt[0], med_pt[1]), int_tol))
	  inside_segments.push_back(shared_ptr<CurveOnSurface>
				    (dynamic_cast<CurveOnSurface*>
				     (curve.subCurve(from_par, to_par))));
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
				      shared_ptr<BoundedSurface>& bounded_sf2)
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
	    bool valid = bounded_sf1->isValid(state);
	    if (!valid)
	      std::cout << "Surface not valid: " << state << std::endl;
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
	    bool valid = bounded_sf2->isValid(state);
	    if (!valid)
	      std::cout << "Surface not valid: " << state << std::endl;
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

    try {
      getIntersectionCurve(under_sf1, under_sf2, int_segments1, int_segments2, 
			   epsge);
    } catch (...) {
      THROW("Failed intersecting spline surface with plane.");
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
				bounded_sf2, /*0.1**/epsge);
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
//   int i;

//     // We extract both parametric and geometric boundary curves.
// // #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
// //     const RectDomain& domain = static_cast<const RectDomain&>(surf.parameterDomain());
// // #else
//     const RectDomain& domain = surf.parameterDomain();
// // #endif

//     CurveLoop space_loop = surf.outerBoundaryLoop(space_epsilon);
//     vector<shared_ptr<ParamCurve> > space_curves;
//     for (i = 0; i < space_loop.size(); ++i)
// 	space_curves.push_back(space_loop[i]);

//     vector<Point> par_pts;
//     par_pts.push_back(Point(domain.umin(), domain.vmin()));
//     par_pts.push_back(Point(domain.umax(), domain.vmin()));
//     par_pts.push_back(Point(domain.umax(), domain.vmax()));
//     par_pts.push_back(Point(domain.umin(), domain.vmax()));

//     // If surface is degenerate, we must create degenerate space_curve.
//     if (space_curves.size() < 4) {
// 	vector<Point> space_pts;
// 	for (i = 0; i < int(par_pts.size()); ++i)
// 	    space_pts.push_back(surf.ParamSurface::point
// 				(par_pts[i][0], par_pts[i][1]));
// 	bool deg[4];
// 	surf.isDegenerate(deg[0], deg[1], deg[2], deg[3], space_epsilon);
// 	for (i = 0; i < 4; ++i)
// 	    if (deg[i]) {
// // 	int nmb_missing_crvs = 4 - space_curves.size();
// // 	for (int i = 0; i < space_pts.size(); ++i)
// // 	    if (space_pts[i].dist(space_pts[(i+1)%4]) < space_epsilon)
// 		space_curves.insert(space_curves.begin() + i,
// 				    shared_ptr<ParamCurve>
// 				    (new SplineCurve(space_pts[i],
// 						       space_pts[(i+1)%4])));
// 	    }
//     }
//     ASSERT(space_curves.size() == 4);

//     vector<shared_ptr<SplineCurve> > param_curves;
//     for (i = 0; i < int(par_pts.size()); ++i)
// 	param_curves.push_back(shared_ptr<SplineCurve>
// 			       (new SplineCurve(par_pts[i], 
// 						  space_curves[i]->startparam(),
// 						  par_pts[(i+1)%4],
// 						  space_curves[i]->endparam())));

    shared_ptr<ParamSurface> under_surf;
// #if _MSC_VER > 0 && _MSC_VER < 1300
//     under_surf = shared_ptr<SplineSurface>
//       (dynamic_cast<SplineSurface*>(surf.clone()));
// #else
    under_surf = shared_ptr<ParamSurface>(surf.clone());
// #endif
    vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(under_surf,
							 space_epsilon);
//     vector<shared_ptr<CurveOnSurface> > surf_curves;
//     for (i = 0; i < 4; ++i)
// 	surf_curves.push_back // We prefer parameter curves.
// 	    (shared_ptr<CurveOnSurface>(new CurveOnSurface(under_surf,
// 							       param_curves[i],
// 							       space_curves[i],
// 							       true)));

    return new BoundedSurface(under_surf, loops);
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

    vector<vector<shared_ptr<CurveOnSurface> > > new_loops;
//===========================================================================
vector< vector< shared_ptr< CurveOnSurface > > >
BoundedUtils::getBoundaryLoops(const BoundedSurface& sf,
			       vector< shared_ptr< CurveOnSurface > >& part_bd_cvs,
			       double eps, int last_split)
//===========================================================================
{
  double a_tol = 1.0e-8;  

    vector<vector<shared_ptr<CurveOnSurface> > > new_loops;

    if (part_bd_cvs.size() == 0)
	THROW("Curve input vector was empty!");
    int ki, kj;

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
    knot_diff_tol = std::max(knot_diff_tol, DEFAULT_PARAMETER_EPSILON);
    for (ki = 0; ki < int(boundary_loops.size()); ++ki) {
      min_loop_tol = std::min(min_loop_tol, 1.5*boundary_loops[ki].getSpaceEpsilon());
      min_loop_tol = std::max(min_loop_tol, 
			      1.1*boundary_loops[ki].getMaxCurveDist());
	for (kj = 0; kj < boundary_loops[ki].size(); ++kj)
	  {
	    shared_ptr<CurveOnSurface> tmp_crv = 
	      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(boundary_loops[ki][kj]);
	    tmp_crv->ensureParCrvExistence(min_loop_tol);
	    old_loop_cvs.push_back(tmp_crv);
	  }
    }

    double min_tool = 1.0e-5;
    min_loop_tol = std::max(min_tool, min_loop_tol);
    min_loop_tol = std::max(min_loop_tol, eps);
    double epspar = min_loop_tol; // Assuming parameter domain reflects the geometry...
    // Since we do not know anything about this and it is mainly a test additional
    // to the geometry space test to avoid confusing seam curves, we make it a little bigger
    epspar *= 10.0;

    // We run part_bd_cvs, splitting if they start/end in inner part of a cv in old_loop_cvs.
    if (last_split < 0)
	last_split = (int)part_bd_cvs.size();
    for (ki = 0; ki < last_split; ++ki) {
	if (!part_bd_cvs[ki]->parameterCurve().get()) {
	    // Suppose we should set curve to prefer parametric part.
	    THROW("Only space curves, method uses parameter curves...");
	}
	Point space_start_pt = part_bd_cvs[ki]->ParamCurve::point(part_bd_cvs[ki]->startparam());
	Point par_start_pt =  part_bd_cvs[ki]->parameterCurve()->point
	    (part_bd_cvs[ki]->parameterCurve()->startparam());
	Point space_end_pt = part_bd_cvs[ki]->ParamCurve::point(part_bd_cvs[ki]->endparam());
	Point par_end_pt =  part_bd_cvs[ki]->parameterCurve()->point
	    (part_bd_cvs[ki]->parameterCurve()->endparam());
	for (kj = 0; kj < int(old_loop_cvs.size()); ++kj) {

	    if (!old_loop_cvs[kj]->parameterCurve().get())
		// Suppose we should set curve to prefer parametric part.
		THROW("Only space curves, method uses parameter curves...");
	    shared_ptr<SplineCurve> pcv =
		dynamic_pointer_cast<SplineCurve, ParamCurve>(old_loop_cvs[kj]->parameterCurve());
	    if (pcv->endparam() - pcv->startparam() < knot_diff_tol) {
	      // Check length of curve in geometry space and parameter space
	      Point geom_start = 
		old_loop_cvs[kj]->ParamCurve::point(old_loop_cvs[kj]->startparam());
	      Point geom_end = 
		old_loop_cvs[kj]->ParamCurve::point(old_loop_cvs[kj]->endparam());
	      Point par_start = 
		old_loop_cvs[kj]->parameterCurve()->point(old_loop_cvs[kj]->startparam());
	      Point par_end = 
		old_loop_cvs[kj]->parameterCurve()->point(old_loop_cvs[kj]->endparam());
	      if (geom_start.dist(geom_end) < min_loop_tol && 
		  par_start.dist(par_end) < epspar)
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
 	    double seed = getSeed(space_start_pt, *old_loop_cvs[kj]);
	    pcv->closestPoint(par_start_pt, pcv->startparam(), pcv->endparam(),
			      start_t, clo_start_pt, clo_dist_start, &seed);
	    //pcv->basis().knotIntervalFuzzy(start_t, knot_diff_tol);
	    // Check against total_start_pt.
	    seed = getSeed(space_end_pt, *old_loop_cvs[kj]);
	    pcv->closestPoint(par_end_pt, pcv->startparam(), pcv->endparam(),
			      end_t, clo_end_pt, clo_dist_end, &seed);
	    //pcv->basis().knotIntervalFuzzy(end_t, knot_diff_tol);
	    Point clo_start_pt_space = sf.ParamSurface::point(clo_start_pt[0], clo_start_pt[1]);
	    Point clo_end_pt_space = sf.ParamSurface::point(clo_end_pt[0], clo_end_pt[1]);

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

	    double space_start_dist = clo_start_pt_space.dist(space_start_pt);
	    double space_end_dist = clo_end_pt_space.dist(space_end_pt);
	    if (std::min(space_start_dist, dist_close_start) < min_loop_tol
		/*true*/ /*dist_close_start < space_start_dist*/)
	    // if (!(space_start_dist < min_loop_tol && 
	    // 	  dist_close_start >= min_loop_tol))
	      {
		space_start_dist = dist_close_start;
		start_t = par_close_start;
		min_loop_tol = std::max(min_loop_tol, dist_close_start+a_tol);
	      }
	    if (std::min(space_end_dist, dist_close_end) < min_loop_tol
		/*true*/ /*dist_close_end < space_end_dist*/)
	    // if (!(space_end_dist < min_loop_tol &&
	    // 	  dist_close_end >= min_loop_tol))
	      {
		space_end_dist = dist_close_end;
		end_t = par_close_end;
		min_loop_tol = std::max(min_loop_tol, dist_close_end+a_tol);
	      }

	    if ((space_start_dist < min_loop_tol) &&
		(start_t - knot_diff_tol > old_loop_cvs[kj]->startparam()) &&
		(start_t + knot_diff_tol < old_loop_cvs[kj]->endparam())) 
	      {
		vector<shared_ptr<ParamCurve> > sub_cvs = 
		  old_loop_cvs[kj]->split(start_t);
		for (size_t k2=0; k2<sub_cvs.size(); ++k2)
		  {
		    shared_ptr<CurveOnSurface> sf_cv = 
		      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[k2]);
		    if (sf_cv.get())
		      sf_cv->ensureParCrvExistence(eps);
		  }
		old_loop_cvs[kj] =  
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[0]);
		shared_ptr<CurveOnSurface> sub_cv =  
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[1]);
		old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub_cv);
		kj++;

		// Check
		Point tmp_pos1 = sub_cvs[0]->point(start_t);
		Point tmp_pos2 = sub_cvs[1]->point(start_t);
		Point tmp_par1 = old_loop_cvs[kj-1]->parameterCurve()->point(start_t);
		Point tmp_par2 = old_loop_cvs[kj]->parameterCurve()->point(start_t);
		Point tmp_sf1 = sf.ParamSurface::point(tmp_par1[0],tmp_par1[1]);
		Point tmp_sf2 = sf.ParamSurface::point(tmp_par2[0],tmp_par2[1]);
		int stop_break;
		stop_break = 1;
		

		// shared_ptr<ParamCurve> sub1(old_loop_cvs[kj]->subCurve(old_loop_cvs[kj]->startparam(), 
		// 						       start_t));
		// Point space_div_1 = sub1->point(sub1->endparam());

		// shared_ptr<CurveOnSurface> sub_cv =  
		//   dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub1);
		
		// sub1 = shared_ptr<ParamCurve>
		//     (old_loop_cvs[kj]->subCurve(start_t, old_loop_cvs[kj]->endparam()));

		// Point space_div_2 = sub1->point(sub1->startparam());

		// old_loop_cvs[kj] =  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub1);
		    
		// old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub_cv);
		// ++kj;
	    }
	    if ((space_end_dist < min_loop_tol) &&
		(end_t - knot_diff_tol > old_loop_cvs[kj]->startparam()) &&
		(end_t + knot_diff_tol < old_loop_cvs[kj]->endparam())) 
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
		old_loop_cvs[kj] =  
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[0]);
		shared_ptr<CurveOnSurface> sub_cv =  
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cvs[1]);
		old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub_cv);
		kj++;

		// Check
		Point tmp_pos1 = sub_cvs[0]->point(end_t);
		Point tmp_pos2 = sub_cvs[1]->point(end_t);
		Point tmp_par1 = old_loop_cvs[kj-1]->parameterCurve()->point(end_t);
		Point tmp_par2 = old_loop_cvs[kj]->parameterCurve()->point(end_t);
		Point tmp_sf1 = sf.ParamSurface::point(tmp_par1[0],tmp_par1[1]);
		Point tmp_sf2 = sf.ParamSurface::point(tmp_par2[0],tmp_par2[1]);
		int stop_break;
		stop_break = 1;
		
		// shared_ptr<ParamCurve> sub1
		//     (old_loop_cvs[kj]->subCurve(old_loop_cvs[kj]->startparam(), end_t));
		// Point space_div_1 = sub1->point(sub1->endparam());

		// shared_ptr<CurveOnSurface> sub_cv =  
		//   dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub1);
		// sub1 = shared_ptr<ParamCurve>
		//     (old_loop_cvs[kj]->subCurve(end_t, old_loop_cvs[kj]->endparam()));
		// Point space_div_2 = sub1->point(sub1->startparam());

		// old_loop_cvs[kj] =  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub1);
		// old_loop_cvs.insert(old_loop_cvs.begin() + kj, sub_cv);
		// ++kj;
	    }
	}
    }

    // // Remove coincident curves. Check first part_bd_cvs
    // for (ki=0; ki<(int)part_bd_cvs.size(); ++ki)
    //   for (kj=ki+1; kj<(int)part_bd_cvs.size(); ++kj)
    // 	{
    // 	  int coincidence = checkCurveCoincidence(part_bd_cvs[ki], part_bd_cvs[kj], 
    // 						  min_loop_tol, true);
    // 	  if (coincidence)
    // 	    {
    // 	      // Coincidence. Remove last curve
    // 	      part_bd_cvs.erase(part_bd_cvs.begin() + kj);
    // 	      kj--;
    // 	    }
    // 	}

    // Check part_bd_cvs agains old_loop_cvs
    
    for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
      for (kj=0; kj<(int)part_bd_cvs.size(); ++kj)
	{
	  int coincidence = checkCurveCoincidence(old_loop_cvs[ki], part_bd_cvs[kj], 
						  min_loop_tol, false);
	  if (coincidence)
	    {
	      // Coincidence. Remove last curve
	      part_bd_cvs.erase(part_bd_cvs.begin() + kj);
	      kj--;
	    }
	}
	    
#ifdef DEBUG1
    std::ofstream out1("loop_cvs1.g2");
    for (ki=0; ki<(int)old_loop_cvs.size(); ++ki)
      {
	out1 << "100 1 0 4 255 0 0 255" << std::endl;
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(old_loop_cvs[ki]->geometryCurve());
	tmp1->write(out1);
      }
    for (ki=0; ki<(int)part_bd_cvs.size(); ++ki)
      {
	out1 << "100 1 0 4 0 255 0 255" << std::endl;
	shared_ptr<SplineCurve> tmp1 = 
	  shared_ptr<SplineCurve>(part_bd_cvs[ki]->geometryCurve());
	tmp1->write(out1);
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
    Point total_space_start_pt = curr_crv->ParamCurve::point(curr_crv->startparam());
    Point curr_space_end_pt = curr_crv->ParamCurve::point(curr_crv->endparam());
    double space_end_dist = total_space_start_pt.dist(curr_space_end_pt);
    Point total_par_start_pt = curr_crv->parameterCurve()->point(curr_crv->startparam());
    Point curr_par_end_pt = curr_crv->parameterCurve()->point(curr_crv->endparam());
    double par_end_dist = total_par_start_pt.dist(curr_par_end_pt);
//     bool cw_loop;
//     double max_angle = -4.0; // Angles are measured on (-pi, pi].
    // We must locate max angle among remaining parts, both new and old.
//     bool loop_locally_closed = false; // We may want to continue after loop is closed.
//     double local_angle = -2*M_PI; // Illegal value.
    // double min_ang = 0.01;

    while (true) { //((end_dist > min_loop_tol) || (part_bd_cvs.size() != 0)) {
// 	double part_angle;
// 	int part_ind = leftMostCurve(*curr_crv, part_bd_cvs, min_loop_tol, 
// 				     part_angle, 2*M_PI);
// 	double old_angle;
// 	int old_ind = leftMostCurve(*curr_crv, old_loop_cvs, min_loop_tol, 
// 				    old_angle, 2*M_PI);
      vector<shared_ptr<CurveOnSurface> > tmp_vec;
      tmp_vec.insert(tmp_vec.end(), part_bd_cvs.begin(), part_bd_cvs.end());
      tmp_vec.insert(tmp_vec.end(), old_loop_cvs.begin(), old_loop_cvs.end());
      double tmp_ang;
      int part_ind = -1, old_ind = -1;
      int tmp_ind = leftMostCurve(*curr_crv, tmp_vec, min_loop_tol, tmp_ang);
      double part_angle = (tmp_ind < (int)part_bd_cvs.size()) ? tmp_ang : 8.0;
      if (tmp_ind < (int)part_bd_cvs.size())
	part_ind = tmp_ind;
      double old_angle = (tmp_ind >= (int)part_bd_cvs.size()) ? tmp_ang : 8.0;
      if (tmp_ind >= (int)part_bd_cvs.size())
	old_ind = tmp_ind - (int)part_bd_cvs.size();

      if (space_end_dist < min_loop_tol && par_end_dist < epspar) {
	// Check if next segment is to the left of the first in curr_loop.
	// 	    vector<shared_ptr<CurveOnSurface> > dummy_vec;
	// 	    dummy_vec.push_back(curr_loop.front());
	tmp_vec.push_back(curr_loop.front());
	double angle;
	tmp_ind = leftMostCurve(*curr_crv, tmp_vec, min_loop_tol, angle);
	if (tmp_ind == (int)tmp_vec.size() - 1)
	  part_ind = old_ind = -1;
// 	if (angle < part_angle + a_tol) {
// 	  part_ind = -1;
// 	}
// 	if (angle < old_angle + a_tol) {
// 	  old_ind = -1;
// 	}
// 	  part_ind = -1;
// 	  old_ind = -1;
      }
	if (part_ind != -1 && old_ind != -1) 
	  { // We choose the smallest angle.
// 	  if (fabs(part_angle - old_angle) < min_ang)
// 	    {
// 	      // @@@ VSK 0111. This part of code is not updated with
// 	      // regard to a modified use of leftMostCurve. Must be updated

// 	      // (Approximate) tangential situation. Recompute angles based on
// 	      // less local information to get a more reliable decision
// 	      // NB! The current choice of evaluation point is quite arbitrary
// 	      // and can be improved (VSK, 201004)
// 	      double t1 = old_loop_cvs[old_ind]->startparam() +
// 		0.1*(old_loop_cvs[old_ind]->endparam()  - 
// 		     old_loop_cvs[old_ind]->startparam());
// 	      double t2 = part_bd_cvs[part_ind]->startparam() +
// 		0.1*(part_bd_cvs[part_ind]->endparam()  - 
// 		     part_bd_cvs[part_ind]->startparam());
// 	      Point tmp1 = old_loop_cvs[old_ind]->parameterCurve()->ParamCurve::point(t1);
// 	      Point tmp2 = part_bd_cvs[part_ind]->parameterCurve()->ParamCurve::point(t2);
// 	      Point vec1 = tmp1 - curr_par_end_pt;
// 	      Point vec2 = tmp2 - curr_par_end_pt;
// 	      vector<Point> curr_val =
// 		curr_crv->parameterCurve()->ParamCurve::point(curr_crv->endparam(), 1);
// 	      old_angle = M_PI - curr_val[1].angle(vec1);
// 	      if (curr_val[1][0]*vec1[1] - vec1[0]*curr_val[1][1] < min_ang)
// 		old_angle += M_PI;
// // 	      if (old_angle + min_ang > M_PI)
// // 		old_angle -= 2.0*M_PI;

// 	      part_angle = M_PI - curr_val[1].angle(vec2);
// 	      if (curr_val[1][0]*vec2[1] - vec2[0]*curr_val[1][1] < min_ang)
// 		part_angle += M_PI;
// // 	      if (part_angle + min_ang > M_PI)
// // 		part_angle -= 2.0*M_PI;
// 	    }

	    if (part_angle <= old_angle)
		old_ind = -1;
	    else
		part_ind = -1;
	  } 
	else if (part_ind == -1 && old_ind == -1) { // No new segment.
	  if (space_end_dist < min_loop_tol /*&& par_end_dist < epspar*/) {
#ifdef DEBUG1
		std::ofstream out2("loop_cvs2.g2");
		for (size_t ix=0; ix<curr_loop.size(); ++ix)
		  {
		    curr_loop[ix]->spaceCurve()->writeStandardHeader(out2);
		    curr_loop[ix]->spaceCurve()->write(out2);
		  }
#endif // DEBUG1

		new_loops.push_back(curr_loop);

	    } // If not this was a dead end ...
	    curr_loop.clear(); // No more segments in loop.
	    // We then look for start segment of next loop.
	    if (part_bd_cvs.size() != 0) {
		part_ind = 0;
	    } else if (old_loop_cvs.size() != 0) {
	       part_ind = -1;
	       old_ind = 0;
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
	if (curr_loop.size() == 0) {
	    total_space_start_pt = curr_crv->ParamCurve::point(curr_crv->startparam());
	    total_par_start_pt = curr_crv->parameterCurve()->point(curr_crv->startparam());
	}

	curr_loop.push_back(curr_crv);
	curr_space_end_pt = curr_crv->ParamCurve::point(curr_crv->endparam());
	curr_par_end_pt = curr_crv->parameterCurve()->point(curr_crv->endparam());
	space_end_dist = total_space_start_pt.dist(curr_space_end_pt);
	par_end_dist = total_par_start_pt.dist(curr_par_end_pt);
    }

    return new_loops; // We should be done
}


//===========================================================================
int BoundedUtils::checkCurveCoincidence(shared_ptr<CurveOnSurface> cv1, 
					shared_ptr<CurveOnSurface> cv2, double tol,
					bool same_orient)
//===========================================================================
{
  int nmb_sample = 2;  // Number of points to check in the inner of a candidate curve

  // Evaluate endpoints
  Point pt1 = cv1->ParamCurve::point(cv1->startparam());
  Point pt2 = cv1->ParamCurve::point(cv1->endparam());
  Point pt3 = cv2->ParamCurve::point(cv2->startparam());
  Point pt4 = cv2->ParamCurve::point(cv2->endparam());
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
  double t1 = cv1->startparam();
  double t2 = cv1->endparam();
  double tdel = (t2 - t1)/(double)(nmb_sample+1);
  double tpar;
  for (kr=0, tpar=t1+tdel; kr<nmb_sample; ++kr, tpar+=tdel)
    {
      Point pt5 = cv1->ParamCurve::point(tpar);
      double clo_par, clo_dist;
      Point clo_pt;
      cv2->closestPoint(pt5, cv2->startparam(), cv2->endparam(), 
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
     if (loopIsDegenerate(loops[ki], epsgeo))
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
    SplineSurface *splinesf = surf->asSplineSurface();
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
    SplineSurface *splinesf = surf->asSplineSurface();
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

    SplineSurface *surf1 = sf1->asSplineSurface();
    SplineSurface *surf2 = sf2->asSplineSurface();
    if (!(surf1 && surf2))
      THROW("No underlying spline surface");

    SISLSurf* sisl_sf1 = GoSurf2SISL(*surf1);
    SISLSurf* sisl_sf2 = GoSurf2SISL(*surf2);
    double epsco = epsgeo; // Actually not used.
    int nmb_int_pts;
    double* pointpar1 = 0;
    double* pointpar2 = 0;
    int nmb_int_cvs;
    SISLIntcurve** intcurves = 0;
    int status = 0;
    int ki;
    //double march_eps = std::min(0.01,100.0*epsgeo); //0.01;
    double march_eps = 0.75*epsgeo; //std::min(0.001,10.0*epsgeo); //0.01;
    //double march_eps = std::min(0.0005,5.0*epsgeo); //0.01;
    s1859(sisl_sf1, sisl_sf2, epsco, epsgeo, &nmb_int_pts,
	  &pointpar1, &pointpar2, &nmb_int_cvs, &intcurves, &status);
    ALWAYS_ERROR_IF(status < 0,
		"Failed intersecting surfs.");
    MESSAGE_IF(status != 0, "Returned status value: " << status);
    if (nmb_int_cvs == 0)
	return; // Surfaces did not intersect.

    double maxstep = (double)0;
    int makecurv = 2; // Make both geometric and parametric curves.
    int draw = 0;
    for (ki = 0; ki < nmb_int_cvs; ++ki) {
	s1310(sisl_sf1, sisl_sf2, intcurves[ki], march_eps, maxstep, makecurv, draw, &status);
	ALWAYS_ERROR_IF(status < 0,
		    "Failed intersecting surfs.");
	MESSAGE_IF(status != 0, "Returned status value: " << status);

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
					   double epsge)
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
    vector<shared_ptr<CurveOnSurface> > new_cvs1, new_cvs2; //(cvs2.size());
    for (ki = 0; ki < int(cvs1.size()); ++ki) {
	vector<shared_ptr<CurveOnSurface> > new_cvs = 
	  intersectWithSurface(*cvs1[ki], *bd_sf1, 0.1*epsge);
	new_cvs1.insert(new_cvs1.end(), new_cvs.begin(), new_cvs.end());
	int other_ind = opp_dir ? (int)cvs1.size() - 1 - ki : ki;
	new_cvs = intersectWithSurface(*cvs2[other_ind], *bd_sf2, 0.1*epsge);
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
    double knot_diff_tol = 1e-05;
    for (ki = 0; ki < int(new_cvs1.size()); ++ki) {
	Point start1 = new_cvs1[ki]->ParamCurve::point(new_cvs1[ki]->startparam());
	Point end1 = new_cvs1[ki]->ParamCurve::point(new_cvs1[ki]->endparam());
	// We split cvs which end in the interior of another cv.
	for (kj = 0; kj < int(new_cvs2.size()); ++kj) {
	    Point start2 = new_cvs2[kj]->ParamCurve::point(new_cvs2[kj]->startparam());
	    Point end2 = new_cvs2[kj]->ParamCurve::point(new_cvs2[kj]->endparam());
	    // We then perform closest pt calculations for all end pts.
	    double clo_t, clo_dist;
	    Point clo_pt;
	    new_cvs1[ki]->ParamCurve::closestPoint(end2, clo_t, clo_pt, clo_dist);
	    if ((clo_dist < epsge) && ((clo_t - knot_diff_tol > new_cvs1[ki]->startparam()) &&
				       (clo_t + knot_diff_tol < new_cvs1[ki]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(new_cvs1[ki]->subCurve(clo_t, new_cvs1[ki]->endparam()));
		new_cvs1.insert(new_cvs1.begin() + ki + 1, new_cv);
		new_cvs1[ki] = shared_ptr<CurveOnSurface>
		    (new_cvs1[ki]->subCurve(new_cvs1[ki]->startparam(), clo_t));
	    }
	    new_cvs1[ki]->ParamCurve::closestPoint(start2, clo_t, clo_pt, clo_dist);
	    if ((clo_dist < epsge) && ((clo_t - knot_diff_tol > new_cvs1[ki]->startparam()) &&
				       (clo_t + knot_diff_tol < new_cvs1[ki]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(new_cvs1[ki]->subCurve(clo_t, new_cvs1[ki]->endparam()));
		new_cvs1.insert(new_cvs1.begin() + ki + 1, new_cv);
		new_cvs1[ki] = shared_ptr<CurveOnSurface>
		    (new_cvs1[ki]->subCurve(new_cvs1[ki]->startparam(), clo_t));
	    }
	    new_cvs2[kj]->ParamCurve::closestPoint(end1, clo_t, clo_pt, clo_dist);
	    if ((clo_dist < epsge) && ((clo_t - knot_diff_tol > new_cvs2[kj]->startparam()) &&
				       (clo_t + knot_diff_tol < new_cvs2[kj]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(new_cvs2[kj]->subCurve(clo_t, new_cvs2[kj]->endparam()));
		new_cvs2.insert(new_cvs2.begin() + kj + 1, new_cv);
		new_cvs2[kj] = shared_ptr<CurveOnSurface>
		    (new_cvs2[kj]->subCurve(new_cvs2[kj]->startparam(), clo_t));
	    }
	    new_cvs2[kj]->ParamCurve::closestPoint(start1, clo_t, clo_pt, clo_dist);
	    if ((clo_dist < epsge) && ((clo_t - knot_diff_tol > new_cvs2[kj]->startparam()) &&
				       (clo_t + knot_diff_tol < new_cvs2[kj]->endparam()))) {
		shared_ptr<CurveOnSurface> new_cv(new_cvs2[kj]->subCurve(clo_t, new_cvs2[kj]->endparam()));
		new_cvs2.insert(new_cvs2.begin() + kj + 1, new_cv);
		new_cvs2[kj] = shared_ptr<CurveOnSurface>
		    (new_cvs2[kj]->subCurve(new_cvs2[kj]->startparam(), clo_t));
	    }
	}
    }

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
	    if ((start1.dist(start2) < epsge) && (end1.dist(end2) < epsge)) {
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
    std::ofstream outfile_curr_bd_sf("tmp/curr_bd_sf.g2");
    SplineDebugUtils::writeTrimmedInfo(*bd_sf, outfile_curr_bd_sf, 0.0);
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
		MESSAGE("Altering loop tol from space_eps = " <<
			space_eps << ", to new_space_eps = " <<
			new_space_eps);
		bd_sf->loop(kj)->setSpaceEpsilon(new_space_eps);
		new_loop_sf_dist[kj] = new_space_eps;
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
#if 0
		writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
#endif
// 			  bd_sf->writeStandardHeader(outfile_failures);
// 			  bd_sf->write(outfile_failures);
// 		++nmb_failures;
		return;
	    }
	}
	if (!sf_ok && (pos_state%4 > 1))
	{
	    // Project missing parameter curves.
	    MESSAGE("State: Missing par cv, trying to fix!");
	    // There is no point in projecting missing parameter curves
	    // if existing curves are not within input tolerance.
	    // We check if we need to enlarge epsgeo.
	    CreatorsUtils::fixTrimCurves(bd_sf);
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
		writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
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
	    writeTrimmedInfo(*bd_sf, outfile_failures, 0.0);
#endif
	}
	else
	{
	    // Writing to file the fixed surfaces.
	    MESSAGE("State: Surface valid after fixing "
		    "trim curves! bd_sf_state = " << bd_sf_state);
#ifdef SBR_DBG
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
  int nmb_sample = 5;
  for (size_t ki=0; ki<loop.size(); ++ki)
    {
      double t1 = loop[ki]->startparam();
      double t2 = loop[ki]->endparam();
      double tdel = (t2 - t1)/(double)(nmb_sample+1);
      double tpar;
      int kj;
      for (kj=0, tpar=t1+tdel; kj<nmb_sample; ++kj, tpar+=tdel)
	{
	  // Evaluate sampling point
	  Point pos;
	  loop[ki]->point(pos, tpar);

	  // Project
	  size_t kr;
	  for (kr=0; kr<loop.size(); ++kr)
	    {
	      if (kr == ki)
		continue;  // Same curve

	      double tpar2, dist;
	      Point pos2;
	      loop[kr]->closestPoint(pos, loop[kr]->startparam(),
				      loop[kr]->endparam(), tpar2,
				      pos2, dist);
	      if (dist < epsgeo)
		break;
	    }

	  if (kr == loop.size())
	    {
	      // The loop is not degenerate
	      return false;
	    }
	}
    }

  // All sampling points coincide with some other part of the loop
  return true;
}


} // end namespace Go

namespace {

//===========================================================================
double getSeed(Point space_pt, CurveOnSurface& cv_on_sf)
//===========================================================================
{
   shared_ptr<SplineCurve> spline_cv;
   if (cv_on_sf.parPref()) {
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
		  double space_eps, double& angle)
//===========================================================================
{
    // A closed surface needs special care, we use cvs starting in the same pt in parametric domain.
    // Tolerance may be somewhat larger than the corresponding spatial value as it is solely used
    // to handle closed surfaces.
  double a_tol = 1.0e-10;  // To include equality in angular test
    RectDomain domain = cv.underlyingSurface()->containingDomain();
    double par_closed_eps = 0.5*std::min(domain.umax() - domain.umin(), domain.vmax() - domain.vmin());

    ALWAYS_ERROR_IF(cv.parameterCurve().get() == 0, "could not find parameter curve");
    Point end_pt = cv.ParamCurve::point(cv.endparam());
    vector<Point> geom_end_pt = cv.ParamCurve::point(cv.endparam(), 1);
    vector<Point> par_end_pt =
	cv.parameterCurve()->ParamCurve::point(cv.parameterCurve()->endparam(), 1);
    int left_most_ind = -1;
    int zero_ind = -1;
    bool opposite = false;
    double min_angle = 8.0; // More than 2 pi
    const double ANG_TOL = 1e-06; // We do not want to return edge returning in the same direction.
    for (int ki = 0; ki < int(other_cvs.size()); ++ki) {
	ALWAYS_ERROR_IF(other_cvs[ki]->parameterCurve().get() == 0, "missing parameter curve");
	Point start_pt = other_cvs[ki]->ParamCurve::point(other_cvs[ki]->startparam());
	vector<Point> geom_start_pt = 
	  other_cvs[ki]->ParamCurve::point(other_cvs[ki]->startparam(), 1);
	vector<Point> par_start_pt = other_cvs[ki]->parameterCurve()->ParamCurve::point
	    (other_cvs[ki]->parameterCurve()->startparam(), 1);
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
void consistentIntersectionDir(ParamCurve& inters_pcv,
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
double getParEps(double space_eps, const ParamSurface *sf)
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



// }; // end anonymous namespace

//     /// Assuming input sets of surfaces are oriented consistently within the set, we return
//     /// the surfaces resulting from pairwise trimming.
//     vector<vector<shared_ptr<BoundedSurface> > >
//       trimSurfsWithSurfs(const vector<shared_ptr<Go::ParamSurface> >& sfs1,
// 			 const vector<shared_ptr<Go::ParamSurface> >& sfs2, double epsge);
    

}; // end anonymous namespace
