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

#include "GoTools/compositemodel/cmUtils.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/utils/CurvatureUtils.h"
#include "GoTools/compositemodel/ftEdge.h"


using std::vector;


namespace Go
{

//===========================================================================
double cmUtils::estimatedCurveLength(ftEdgeBase* edge, int nmb_samples)
//===========================================================================
{
    ALWAYS_ERROR_IF(nmb_samples < 2,
		"At least 2 points needed to estimate curve length!");

    double total_length = 0.0;
    double tmin = edge->tMin();
    double tmax = edge->tMax();
    double tstep = (tmax - tmin) / (nmb_samples - 1);
    for (int ki = 0; ki < nmb_samples - 1; ++ki) {
	Point left_pt = edge->point(tmin + ki*tstep);
	Point right_pt = edge->point(tmin + (ki + 1)*tstep);
	total_length += (right_pt - left_pt).length();
    }

    return total_length;
}

//===========================================================================
RectDomain cmUtils::geometricParamDomain(ParamSurface* sf)
//===========================================================================
{
  double len_u, len_v;
  GeometryTools::estimateSurfaceSize(*sf, len_u, len_v);

  RectDomain domain = sf->containingDomain();

  double u_from = domain.umin();
  double u_to = u_from + len_u;
  double v_from = domain.vmin();
  double v_to = v_from + len_v;
  
  Vector2D ll(u_from, v_from);
  Vector2D ur(u_to, v_to);
  
  return RectDomain(ll, ur);
}

//===========================================================================
bool cmUtils::abovePlane(Point pt, Point plane_pt, Point normal)
//===========================================================================
{
    Point vector = pt - plane_pt;
    double inner_product = vector*normal;

    return (inner_product > 0); // Simple test, but nonetheless true
}

    //===========================================================================
void 
cmUtils::extendWithDegBd(vector<int>& corner, 
			 vector< shared_ptr<ParamCurve> >& bd_curves, 
			 vector< shared_ptr<ParamCurve> >& cross_curves, 
			 int idxmin)
//---------------------------------------------------------------------------
//
// Purpose: Extend the set of boundary curves with degenerate curves if
//          a degenerate surface is to be created.
//
//===========================================================================
{
  shared_ptr<ParamCurve> dummy;  // Empty curve
  if (bd_curves.size() == 2)
    {
      // Every second curve must be degenerate. First evaluate the
      // existing curves in their start- and endpoint.
      Point startp1, startp2, endp1, endp2;
      bd_curves[0]->point(startp2, bd_curves[0]->startparam());
      bd_curves[0]->point(endp2, bd_curves[0]->endparam());
      bd_curves[1]->point(startp1, bd_curves[1]->startparam());
      bd_curves[1]->point(endp1, bd_curves[1]->endparam());

      // Make degenerate curves and insert them into the 
      // vector of boundary curves.
      shared_ptr<ParamCurve> degcv1, degcv2;
      SplineCurve *degspline1 = new SplineCurve(endp1, startp2);
      degspline1->raiseOrder(4-degspline1->order());
      SplineCurve *degspline2 = new SplineCurve(endp2, startp1);
      degspline2->raiseOrder(4-degspline2->order());
      degcv1 = shared_ptr<ParamCurve>(degspline1);
      degcv2 = shared_ptr<ParamCurve>(degspline2);
      bd_curves.insert(bd_curves.begin(), 1, degcv1);
      bd_curves.insert(bd_curves.begin() + 2, 1, degcv2);
      cross_curves.insert(cross_curves.begin(), 1, dummy);
      cross_curves.insert(cross_curves.begin() + 2, 1, dummy);

      // Update the corner indices to indicate the existence of
      // a degenerate boundary.
      int corner0 = corner[0];
      int corner2 = corner[2];
      corner.insert(corner.begin(), 1, corner0);
      corner.insert(corner.begin() + 2, 1, corner2);
    }
  else if (bd_curves.size() == 3)
    {
      // Insert a degenerate curve opposite to the boundary
      // curve of minimum length.
      int idx1, idx2;
      idx1 = (idxmin + 1) % 3;
      idx2 = (idx1 + 1) % 3;
      Point pnt1, pnt2;
      bd_curves[idx1]->point(pnt1, bd_curves[idx1]->endparam());
      bd_curves[idx2]->point(pnt2, bd_curves[idx2]->startparam());
      shared_ptr<ParamCurve> degcv 
	= shared_ptr<ParamCurve>(new SplineCurve(pnt1, pnt2));
      bd_curves.insert(bd_curves.begin() + idx2, 1, degcv);
      cross_curves.insert(cross_curves.begin() + idx2, 1, dummy);
      int corner_idx2 = corner[idx2];
      corner.insert(corner.begin() + idx2, 1, corner_idx2);
    }
}

    
//===========================================================================
void cmUtils::updateWithNewCorners(const vector<ftEdgeBase*>& outer_loop, vector<int>& corners,
				   const vector<Point>& add_corner_pts, double eps_geo)
//===========================================================================
{
    for (size_t ki = 0; ki < add_corner_pts.size(); ++ki) {
	int clo_edge_ind = -1;
	double clo_par;
	double global_clo_dist = 1.0e10;
	for (size_t kj = 0; kj < outer_loop.size(); ++kj) {
	    double clo_t, clo_dist;
	    Point clo_pt;
	    outer_loop[kj]->closestPoint(add_corner_pts[ki], clo_t, clo_pt, clo_dist);
	    if (clo_dist < global_clo_dist) {
		global_clo_dist = clo_dist;
		clo_edge_ind = (int)kj;
		clo_par = clo_t;
	    }
	}
	if (global_clo_dist < eps_geo) {
	    double knot_diff_tol = 1e-05;
	    if (fabs(clo_par - outer_loop[clo_edge_ind]->tMin()) < knot_diff_tol) {
		size_t kj;
		for (kj = 0; kj < corners.size() - 1; ++kj) {
		    if ((corners[kj] < clo_edge_ind) && (corners[kj+1] > clo_edge_ind)) {
			corners.insert(corners.begin() + kj + 1, clo_edge_ind);
			break;
		    }
		}
		if (kj == corners.size() - 1) {
		    corners.push_back(clo_edge_ind);
		}
	    } else if (fabs(clo_par - outer_loop[clo_edge_ind]->tMax()) < knot_diff_tol) {
		int new_cn = (int)((clo_edge_ind + 1)%(outer_loop.size()));
		size_t kj;
		for (kj = 0; kj < corners.size() - 1; ++kj) {
		    if ((corners[kj] < new_cn) && (corners[kj+1] > new_cn)) {
			corners.insert(corners.begin() + kj + 1, new_cn);
			break;
		    }
		}
		if (kj == corners.size() - 1) {
		    corners.push_back(new_cn);
		}
	    } else {
		// @@sbr Seems to me that the edges were not split prior to function call...
		MESSAGE("Should not happen!");
	    }
	}
    }

    std::sort(corners.begin(), corners.end());
}


    
//===========================================================================
bool cmUtils::insideBdTrain(const Vector2D& sample_pt, const vector<Vector2D>& bd_train)
//===========================================================================
{
   ASSERT(bd_train.size() > 1);

   int nmb_bd = (int)bd_train.size();
   vector<double> intersections;
   int ki, kj;
   int par_ind;
   int other_par_ind;
   for (ki = 0; ki < 2; ++ki) {
       par_ind = ki;
       other_par_ind = 1 - ki;
       for (kj = 0; kj < nmb_bd; ++kj) {
	   Vector2D seg_from = bd_train[kj];
	   Vector2D seg_to = bd_train[(kj + 1)%nmb_bd];
	   // We check whether the edge crosses line through node.
	   // @@sbr Not taking care of scenario w/two segments lying after one another, on the same side.
	   if (((seg_from[ki] <= sample_pt[ki]) && (seg_to[ki] > sample_pt[ki])) ||
	       ((seg_from[ki] >= sample_pt[ki]) && (seg_to[ki] < sample_pt[ki]))) {
	       // We must extract intersection parameter for param_node-edge.
	       intersections.push_back(seg_from[other_par_ind] + (sample_pt[ki] - seg_from[ki])*
				       (seg_to[other_par_ind] - seg_from[other_par_ind])/
				       (seg_to[ki] - seg_from[ki]));
	   }
       }
       if ((intersections.size())%2 == 0) {
	   break;
       }
   }

   std::sort(intersections.begin(), intersections.end());
   std::vector<double>::const_iterator iter = 
     std::lower_bound(intersections.begin(),intersections.end(), 
		      sample_pt[other_par_ind]);
   int nmb_to_left = (int)(iter - intersections.begin());
   int nmb_to_right = (int)((std::vector<double>::const_iterator) 
			    intersections.end()-iter);

   if ((nmb_to_left + nmb_to_right)%2 != 0) {
       MESSAGE("Intersection method expecting even number of intersections... Fix method!");
       if ((nmb_to_left == 0) || (nmb_to_right == 0)) {
	   return false;
       }
   }

   return ((nmb_to_left%2) != 0);
}


    //===========================================================================
void cmUtils::reparametrizeBdCvs(const vector<shared_ptr<SplineCurve> >& bd_cvs,
				 double appr_tol,
				 vector<shared_ptr<SplineCurve> >& new_bd_cvs)
//===========================================================================
{
  int nmb_sample_pts = 1000;
 vector<int> rep(4, 1);
 int fix_ind = 1;
 for (size_t ki = 0; ki < bd_cvs.size(); ++ki)
   if (rep[ki])
     {
       double tmin = bd_cvs[ki]->startparam();
       double tmax = bd_cvs[ki]->endparam();
       double tstep = (tmax - tmin)/((double)(nmb_sample_pts - 1));
       // First we sample points.
       int der = 2;
       vector<vector<Point> > sample_der_pts;
       vector<double> sample_pts;
       vector<double> pars; // @@sbr Are these values really used?
       vector<double> curvatures;
       double sum_curv = 0.0;
       for (int kj = 0; kj < nmb_sample_pts; ++kj)
	 {
	   double tpar = tmin + kj*tstep;
	   vector<Point> der_pt(3);
	   bd_cvs[ki]->point(der_pt, tpar, der);
	   sample_der_pts.push_back(der_pt);
	   sample_pts.insert(sample_pts.end(), der_pt[0].begin(), der_pt[0].end());
	   pars.push_back(tpar);
	   vector<Point> unit_der;
	   double curv = 1.0/curvatureRadius(der_pt, unit_der);

	   if (int(ki) == fix_ind) // @@sbr Rather case specific...
	     { // We want Length to 0.25 in start of cv, 1.0 in end (as fraction of orig length).
	       if (kj > 980)
		 curv = curvatures.back();
	     }
	   sum_curv += curv;
	   curvatures.push_back(curv);
	 }

       // We then assign parameter values to the points, according to curvature.
       // Get information about tangent lengts and parameter interval
       vector<double>::const_iterator min_curv = min_element(curvatures.begin(), curvatures.end());
       vector<double>::const_iterator max_curv = max_element(curvatures.begin(), curvatures.end());
       vector<double> new_pars;
       new_pars.push_back(pars[0]);
       for (size_t kj = 0; kj < sample_der_pts.size() - 1; ++kj)
	 {
	   double parint, len1, len2; // @@sbr len1 & len2 not used.
	   getHermiteData(sample_der_pts[kj], sample_der_pts[kj+1], parint, len1, len2);

	   double new_par;
	   // 	  new_par = new_pars.back() + curvatures[kj]*(tmax - tmin)/sum_curv;
	   double a = 2.0; //1.75; //1.9; //1.5;
	   double b = 0.4; //0.25; //0.1; //0.5;
	   new_par = new_pars.back() +
	     (a*(curvatures[kj] - *min_curv) + b*(*max_curv - curvatures[kj]))/
	     (*max_curv - *min_curv);
	   new_pars.push_back(new_par);
	 }

       // The next step is to create an approximation of the curve.
       int dim = 3;
       vector<double> knots(bd_cvs[ki]->basis().begin(), bd_cvs[ki]->basis().end());
       ApproxCurve appr_cv(sample_pts, new_pars, dim, appr_tol,
			     bd_cvs[ki]->numCoefs(), bd_cvs[ki]->order()); //, knots);

       double maxdist, avdist;
       int max_iter = 10;
       new_bd_cvs.push_back(appr_cv.getApproxCurve(maxdist, avdist, max_iter));
     }
   else
     {
       new_bd_cvs.push_back(bd_cvs[ki]);
     }
}

//===========================================================================
void cmUtils::reparametrizeBdCvs2(const vector<shared_ptr<SplineCurve> >& bd_cvs,
				  double appr_tol,
				  vector<shared_ptr<SplineCurve> >& new_bd_cvs)
//===========================================================================
{
    int nmb_sample_pts = 1000;
    int nmb_avg = 5; // Number to use for avg curvature, in each direction.
    vector<int> rep(4, 1);
    int fix_ind1 = 1;
    //  int fix_ind2 = 2;
    for (size_t ki = 0; ki < bd_cvs.size(); ++ki)
	if (rep[ki]) {
	    double tmin = bd_cvs[ki]->startparam();
	    double tmax = bd_cvs[ki]->endparam();
	    double tstep = (tmax - tmin)/((double)(nmb_sample_pts - 1));
	    // First we sample points.
	    int der = 2;
	    vector<vector<Point> > sample_der_pts;
	    vector<double> sample_pts;
	    vector<double> pars; // @@sbr Are these values really used?
	    vector<double> curvatures;
	    vector<double> avg_curvatures;
	    double sum_curv = 0.0;
	    for (int kj = 0; kj < nmb_sample_pts; ++kj) {
		double tpar = tmin + kj*tstep;
		vector<Point> der_pt(3);
		bd_cvs[ki]->point(der_pt, tpar, der);
		sample_der_pts.push_back(der_pt);
		sample_pts.insert(sample_pts.end(), der_pt[0].begin(), der_pt[0].end());
		pars.push_back(tpar);
		vector<Point> unit_der;
		double curv = 1.0/curvatureRadius(der_pt, unit_der);

		if (int(ki) == fix_ind1) { // @@sbr Rather case specific...
		    // We want Length to 0.25 in start of cv, 1.0 in end (as fraction of orig length).
		    if (kj > 980)
			curv = curvatures.back();
		}
		sum_curv += curv;
		curvatures.push_back(curv);
	    }

	    //        if (ki == fix_ind2)
	    // 	 {
	    // 	   for (kj = 20; kj < curvatures.size(); ++kj)
	    // 	     curvatures[kj] = curvatures.back();
	    // 	 }

	    for (size_t kj = 0; kj < curvatures.size(); ++kj) {
		double sum = curvatures[ki];
		int nmb = 1;
		for (int kk = 1; kk < nmb_avg + 1; ++kk) {
		    if (int(kj - kk) > -1) {
			sum += curvatures[kj-kk];
			++nmb;
		    }
		    if (int(kj) + kk < nmb_sample_pts) {
			sum += curvatures[kj+kk];
			++nmb;
		    }
		}
		avg_curvatures.push_back(sum/(double)nmb);
	    }

	    // We then assign parameter values to the points, according to curvature.
	    // Get information about tangent lengts and parameter interval
	    vector<double>::const_iterator min_curv = min_element(avg_curvatures.begin(), avg_curvatures.end());
	    vector<double>::const_iterator max_curv = max_element(avg_curvatures.begin(), avg_curvatures.end());
	    vector<double> new_pars;
	    new_pars.push_back(pars[0]);
	    for (size_t kj = 0; kj < sample_der_pts.size() - 1; ++kj) {
		double parint, len1, len2; // @@sbr len1 & len2 not used.
		getHermiteData(sample_der_pts[kj], sample_der_pts[kj+1], parint, len1, len2);

		double new_par;
		// 	  new_par = new_pars.back() + curvatures[kj]*(tmax - tmin)/sum_curv;
		double a = 2.0; //1.75; //1.9; //1.5;
		double b = 0.4; //0.25; //0.1; //0.5;
		new_par = new_pars.back() +
		    (a*(avg_curvatures[kj] - *min_curv) + b*(*max_curv - avg_curvatures[kj]))/
		    (*max_curv - *min_curv);
		new_pars.push_back(new_par);
	    }

	    // The next step is to create an approximation of the curve.
	    int dim = 3;
	    vector<double> knots(bd_cvs[ki]->basis().begin(), bd_cvs[ki]->basis().end());
	    ApproxCurve appr_cv(sample_pts, new_pars, dim, appr_tol,
				  bd_cvs[ki]->numCoefs(), bd_cvs[ki]->order()); //, knots);

	    double maxdist, avdist;
	    int max_iter = 10;
	    new_bd_cvs.push_back(appr_cv.getApproxCurve(maxdist, avdist, max_iter));
	} else {
	    new_bd_cvs.push_back(bd_cvs[ki]);
	}
}


    //===========================================================================
void cmUtils::splitEdgesInCorners(vector<shared_ptr<ftEdgeBase> >& edges,
				  const vector<Point>& corner_pts, double eps_geo)
//===========================================================================
{
    for (size_t ki = 0; ki < corner_pts.size(); ++ki) {
	int clo_edge_ind = -1;
	double clo_par;
	double global_clo_dist = 1.0e10;
	for (size_t kj = 0; kj < edges.size(); ++kj) {
	    double clo_t, clo_dist;
	    Point clo_pt;
	    edges[kj]->closestPoint(corner_pts[ki], clo_t, clo_pt, clo_dist);
	    if (clo_dist < global_clo_dist) {
		global_clo_dist = clo_dist;
		clo_edge_ind = (int)kj;
		clo_par = clo_t;
	    }
	}
	if (global_clo_dist < eps_geo) {
	    double knot_diff_tol = 1e-05;
	    if ((fabs(clo_par - edges[clo_edge_ind]->tMin()) > knot_diff_tol) &&
		((fabs(clo_par - edges[clo_edge_ind]->tMax()) > knot_diff_tol))) {
		shared_ptr<ftEdgeBase> new_edge(edges[clo_edge_ind]->split(clo_par));
		edges.insert(edges.begin() + clo_edge_ind + 1, new_edge);
	    }
	}
    }
}


    //===========================================================================
void cmUtils::cwOrientation(vector<ftEdgeBase*>& meeting_edges, vector<bool>& start,
			    double angle_tol)
//===========================================================================
{
#ifdef FANTASTIC_DEBUG
    std::ofstream debug("data/debug.g2");
    double knot_diff_tol = 1e-05;
    for (size_t ki = 0; ki < meeting_edges.size(); ++ki) {
	shared_ptr<ParamCurve> edge_cv;
	try {
	    edge_cv = shared_ptr<ParamCurve>
		(meeting_edges[ki]->geomEdge()->geomCurve()->subCurve
		 (meeting_edges[ki]->tMin(), meeting_edges[ki]->tMax(), knot_diff_tol));
	    SplineCurve* edge_space_cv = edge_cv->geometryCurve();
	    if (edge_space_cv) {
		edge_space_cv->writeStandardHeader(debug);
		edge_space_cv->write(debug);
	    }
	} catch (...) {
	    MESSAGE("Failed extracting sub curve.");
	}
    }
#endif // FANTASTIC_DEBUG

    vector<pair<double, int> > angles; // We store index in meeting_edges and tangent in common vertex.
    double tpar = start[0] ? meeting_edges[0]->tMin() : meeting_edges[0]->tMax();
    Point ref_tangent = meeting_edges[0]->tangent(tpar);
    if (!start[0]) {
	ref_tangent *= -1.0;
    }
    Point ref_normal(0.0, 0.0, 0.0);
    // We must middle normal to handle edges meeting sharply (angle less than 90 degrees).
    for (size_t ki = 0; ki < meeting_edges.size(); ++ki) {
	tpar = start[ki] ? meeting_edges[ki]->tMin() : meeting_edges[ki]->tMax();
	ref_normal += meeting_edges[ki]->normal(tpar);
    }
    ref_normal.normalize();

    // And then we project the ref_tangent onto plane defined by ref_normal.
    // We first find one vec orthogonal to normal.
    Point vec1 = ref_normal%ref_tangent;
    vec1.normalize();
    Point vec2 = ref_normal%vec1;
    vec2.normalize();
    double coef1, coef2;
    CoonsPatchGen::blendcoef(vec1.begin(), vec2.begin(), ref_tangent.begin(),
			       3, 1, &coef1, &coef2);
    ref_tangent = coef1*vec1 + coef2*vec2;

    for (size_t ki = 1; ki < meeting_edges.size(); ++ki) {
	tpar = start[ki] ? meeting_edges[ki]->tMin() : meeting_edges[ki]->tMax();
	Point tangent = meeting_edges[ki]->tangent(tpar);
	if (!start[ki]) {
	    tangent *= -1.0;
	}
	// We must then project the angle onto plane defined by ref_normal.
	CoonsPatchGen::blendcoef(vec1.begin(), vec2.begin(), tangent.begin(),
				   3, 1, &coef1, &coef2);
	tangent = coef1*vec1 + coef2*vec2;

	double ccw_angle = ccwAngle(ref_tangent, tangent, &ref_normal);
	double cw_angle = 2*M_PI - ccw_angle;
	if (fabs(2*M_PI - cw_angle) < angle_tol) {
	    cw_angle = 0.0;
	}
	angles.push_back(std::make_pair(cw_angle, ki));
    }
    sort(angles.begin(), angles.end());

    vector<ftEdgeBase*> sorted_meeting_edges;
    vector<bool> sorted_start;
    sorted_meeting_edges.push_back(meeting_edges[0]);
    sorted_start.push_back(start[0]);
    for (size_t ki = 0; ki < angles.size(); ++ki) {
	// We must check if current angle and the next are almost equal. If so we may reorder.
	// (Almost equal angles is obviously quite normal in a graph.)
	if ((ki < angles.size() - 1) && (fabs(angles[ki].first - angles[ki+1].first) < angle_tol) &&
	    (sorted_meeting_edges.back()->face()->getId() ==
	     meeting_edges[angles[ki+1].second]->face()->getId())) {
	    swap(angles[ki], angles[ki+1]);
	}
	sorted_meeting_edges.push_back(meeting_edges[angles[ki].second]);
	sorted_start.push_back(start[angles[ki].second]);
    }

    // We then set our return values.
    meeting_edges = sorted_meeting_edges;
    start = sorted_start;
}

//===========================================================================
void cmUtils::cwOrientation2(vector<ftEdgeBase*>& meeting_edges, vector<bool>& start)
//===========================================================================
{
    vector<ftEdgeBase*> sorted_edges;
    vector<bool> sorted_start;
    vector<int> twin_less_edges;
    // If any of the input edges lack a twin, we choose the one ending in vertex.
    for (size_t ki = 0; ki < meeting_edges.size(); ++ki) {
	if (meeting_edges[ki]->twin() == 0) {
	    twin_less_edges.push_back((int)ki);
	}
    }
    if (twin_less_edges.size() != 0 && twin_less_edges.size() != 2) {
	// @@sbr If this is to happen we may solve the situation by dividing into
	// two separate groups.
	// Currently not really using input edges when sorting...
	THROW("Two or no edges should miss their twin!");
    }

    int first_edge;
    if (twin_less_edges.size() == 0) {
	size_t ki;
	for (ki = 0; ki < start.size(); ++ki) {
	    if (!start[ki]) {
		break;
	    }
	}
	ASSERT(ki < start.size());
	first_edge = (int)ki;
    } else {
	if (start[twin_less_edges[0]]) {
	    first_edge = twin_less_edges[1];
	} else {
	    first_edge = twin_less_edges[1];
	}
    }
    sorted_edges.push_back(meeting_edges[first_edge]);
    sorted_start.push_back(false);

    int nmb_faces = (int)meeting_edges.size()/2;
    for (int ki = 0; ki < nmb_faces - 1; ++ki) {
	sorted_edges.push_back(sorted_edges.back()->next());
	sorted_start.push_back(true);
	sorted_edges.push_back(sorted_edges.back()->twin());
	sorted_start.push_back(false);
    }
    sorted_edges.push_back(sorted_edges.back()->next());
    sorted_start.push_back(true);

    // As first edge is to be start in vertex we move first last.
    sorted_edges.push_back(sorted_edges.front());
    sorted_edges.erase(sorted_edges.begin());
    sorted_start.push_back(sorted_start.front());
    sorted_start.erase(sorted_start.begin());

    meeting_edges = sorted_edges;
    start = sorted_start;
}


    //===========================================================================
vector<pair<shared_ptr<ParamCurve>, ftFaceBase*> >
cmUtils::getG1FaceCurves(ftCurve& bd_curve)
//===========================================================================
{
    int nmb_segs = bd_curve.numSegments();
    vector<vector<ftCurveSegment> > g1_segments;
    vector<ftCurveSegment> curr_vec;
    curr_vec.push_back(bd_curve.segment(0));
    if ((curr_vec.back().jointAfter() >= JOINT_G0) ||
	((nmb_segs > 1) && (curr_vec.back().face(0) != bd_curve.segment(1).face(0)))) {
	g1_segments.push_back(curr_vec);
	curr_vec.clear();
    }
    for (int ki = 1; ki < nmb_segs; ++ki) {
	const ftCurveSegment& curr_seg = bd_curve.segment(ki);
	curr_vec.push_back(curr_seg);
	if ((curr_vec.back().jointAfter() >= JOINT_G0) ||
	    ((ki + 1) < nmb_segs && (curr_vec.back().face(0) != bd_curve.segment(ki+1).face(0)))) {
	    g1_segments.push_back(curr_vec);
	    curr_vec.clear();
	}
    }

    // If first and last elements in vector are of same type, meeting G1, elements are put
    // in a common vector.
    if ((g1_segments.back().back().face(0) == g1_segments.front().front().face(0)) &&
	(g1_segments.back().back().jointAfter() < JOINT_G0)) {
	g1_segments.front().insert(g1_segments.front().begin(),
				   g1_segments.back().begin(), g1_segments.back().end());
    }

    vector<pair<shared_ptr<ParamCurve>, ftFaceBase*> > g1_face_curves;
    // Finally we make continuous curves.
    for (size_t ki = 0; ki < g1_segments.size(); ++ki) {
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	shared_ptr<ParamCurve> edge_cv((ParamCurve*)(g1_segments[ki][0].spaceCurve()->clone()));
#else
	shared_ptr<ParamCurve> edge_cv(g1_segments[ki][0].spaceCurve()->clone());
#endif
	shared_ptr<ParamCurve> total_space_cv;
	if (edge_cv->instanceType() == Class_SplineCurve) {
	    total_space_cv = edge_cv;
	} else if (edge_cv->instanceType() == Class_CurveOnSurface) {
	    total_space_cv = (dynamic_pointer_cast<CurveOnSurface, ParamCurve>(edge_cv))->spaceCurve();
	} else {
	    THROW("Unexpected curve type!");
	}
	for (size_t kj = 1; kj < g1_segments[ki].size(); ++kj) {
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	    edge_cv = shared_ptr<ParamCurve>((ParamCurve*)(g1_segments[ki][kj].spaceCurve()->clone()));
#else
	    edge_cv = shared_ptr<ParamCurve>(g1_segments[ki][kj].spaceCurve()->clone());
#endif
	    shared_ptr<ParamCurve> space_cv;
	    if (edge_cv->instanceType() == Class_SplineCurve) {
		space_cv = edge_cv;
	    } else if (edge_cv->instanceType() == Class_CurveOnSurface) {
		space_cv = (dynamic_pointer_cast<CurveOnSurface, ParamCurve>(edge_cv))->spaceCurve();
	    } else {
		THROW("Unexpected curve type!");
	    }
	    total_space_cv->appendCurve(space_cv.get());
	}
	g1_face_curves.push_back(std::make_pair(total_space_cv, g1_segments[ki][0].face(0)));
    }

    return g1_face_curves;
}

    
//===========================================================================
double cmUtils::ccwAngle(Point first_leg, Point second_leg, Point* normal)
//===========================================================================
{
    if (first_leg.dimension() == 3) {
	ASSERT(normal != NULL);
    }

    double angle = first_leg.angle(second_leg);
    if (first_leg.dimension() == 2) {
	if (first_leg[0]*second_leg[1] - first_leg[1]*second_leg[0] < 0.0) {
	    angle = 2*M_PI - angle;
	}
    } else {
	if ((first_leg%second_leg).dist(-*normal) < (first_leg%second_leg).dist(*normal)) {
	    angle = 2*M_PI - angle;
	}
    }

    return angle;
}


//===========================================================================
std::vector<int> cmUtils::removeInnerCorners(const std::vector<ftEdgeBase*>& outer_loop, std::vector<int>& corners)
//===========================================================================
{
    std::vector<int> new_corners;
    // By allowing inner corners, defined as corner points attached to with multiple surfaces, we may
    // handle such cases. This approach will (currently) not handle cases where less than 4 of the
    // corners are affiliated with multiple surfaces. We expect all surfaces in the set to be smooth and
    // should not experience cases with more than 4 resulting corner points.
    // Another approach would be to let the user define the 4 corner points.
    const int loop_size = outer_loop.size();
    for (size_t ki = 0; ki < corners.size(); ++ki)
    {
        ftEdgeBase* curr_edge = outer_loop[corners[ki]];
        int prev_ind = (corners[ki] + loop_size - 1)%(loop_size);
        ftEdgeBase* prev_loop_edge = outer_loop[prev_ind];
        ftFaceBase* curr_face = curr_edge->geomEdge()->face();
        ftFaceBase* prev_face = prev_loop_edge->geomEdge()->face();
        if (curr_face == prev_face)
        {
        // vector<ftEdgeBase*> adj_edges;
        // vector<bool> at_start;
        // curr_edge->adjacentEdges(true, adj_edges, at_start);
        // if (adj_edges.size() == 1)
        // {
            new_corners.push_back(corners[ki]);
        }
    }

    return new_corners;
}

    
} // namespace Go
