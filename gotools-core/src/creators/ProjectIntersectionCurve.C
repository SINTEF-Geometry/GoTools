//===========================================================================
//                                                                           
// File: ProjectIntersectionCurve.C                                                    
//                                                                           
// Created: Tue Nov 26 09:00:22 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: ProjectIntersectionCurve.C,v 1.4 2005-11-17 08:55:30 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/creators/ProjectIntersectionCurve.h"

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/PointCloud.h"

#include <fstream>

using namespace Go;
using std::vector;
using std::shared_ptr;


//===========================================================================
ProjectIntersectionCurve::ProjectIntersectionCurve(std::shared_ptr<SplineCurve>& inters_crv,
						   std::shared_ptr<SplineCurve>& p_crv,
						   std::shared_ptr<SplineCurve>& other_p_crv,
						   std::shared_ptr<ParamSurface>& surf,
						   std::shared_ptr<ParamSurface>& other_surf,
						   double offset_dist, double other_offset_dist,
						   double epsgeo)
    : inters_crv_(inters_crv), p_crv_(p_crv), other_p_crv_(other_p_crv),
      surf_(surf), other_surf_(other_surf),
      offset_dist_(offset_dist), other_offset_dist_(other_offset_dist), epsgeo_(epsgeo)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(inters_crv_.get() == 0 || p_crv.get() == 0 || other_p_crv.get() == 0 ||
		surf_.get() == 0 || other_surf_.get() == 0,
		"Missing input data.");
    ALWAYS_ERROR_IF(inters_crv_->dimension() != 3 ||
		inters_crv_->dimension() != surf_->dimension() ||
		surf_->dimension() != other_surf_->dimension() ||
		p_crv->dimension() != 2 || other_p_crv->dimension() != 2,
		"Dimension mismatch.");
}

//===========================================================================
ProjectIntersectionCurve::~ProjectIntersectionCurve()
//===========================================================================
{
}

//===========================================================================
Point ProjectIntersectionCurve::eval(double t) const
//===========================================================================
{
    Point space_pt = inters_crv_->ParamCurve::point(t);
    // We find parameter values of closest points on respective surfs.
    double clo_u, clo_v, clo_dist, other_clo_u, other_clo_v, other_clo_dist;
    Point clo_pt, other_clo_pt;
    double epsilon = 1e-10;
    vector<double> seed;
    //if (p_crv_.get() != 0)
    seed = getSuggestedSurfaceParameter(*inters_crv_, t, *p_crv_, false);
    surf_->closestBoundaryPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist, epsgeo_, NULL, &seed[0]);
    vector<double> max_seed_error(2);
    if (seed.size() == 2) {
	RectDomain under_domain = surf_->containingDomain();
	max_seed_error[0] = 0.1*(under_domain.umax() - under_domain.umin());
	max_seed_error[1] = 0.1*(under_domain.vmax() - under_domain.vmin());
	if ((fabs(clo_u - seed[0]) > max_seed_error[0]) || (fabs(clo_v - seed[1]) > max_seed_error[1])) {
	    surf_->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
				epsgeo_, NULL, &seed[0]);
	}
    }
    //if (other_p_crv_.get() != 0)
    seed = getSuggestedSurfaceParameter(*inters_crv_, t, *other_p_crv_, true);
    // else
    // seed.clear();
    other_surf_->closestBoundaryPoint(space_pt, other_clo_u, other_clo_v, other_clo_pt,
				      other_clo_dist, epsgeo_, NULL, &seed[0]);
    if (seed.size() == 2) {
	RectDomain under_domain = other_surf_->containingDomain();
	max_seed_error[0] = 0.1*(under_domain.umax() - under_domain.umin());
	max_seed_error[1] = 0.1*(under_domain.vmax() - under_domain.vmin());
	if ((fabs(other_clo_u - seed[0]) > max_seed_error[0]) ||
	    (fabs(other_clo_v - seed[1]) > max_seed_error[1])) {
	    other_surf_->closestPoint(space_pt, other_clo_u, other_clo_v, other_clo_pt,
				      other_clo_dist, epsgeo_, NULL, &seed[0]);
	}
    }
    // @@ We ought to test dist to actual point.
    Point normal, other_normal;
    surf_->normal(normal, clo_u, clo_v);
    other_surf_->normal(other_normal, other_clo_u, other_clo_v);
    Point offset_space_pt = space_pt + offset_dist_*normal + other_offset_dist_*other_normal;
    // We next must project onto surf_. @@ Should use a sensible seed.
    // Our best guess for a seed so far is (clo_u, clo_v). Suppose we could project the offset curve
    // ont surface and thus approximate the seed.
    seed.resize(2);
    seed[0] = clo_u;
    seed[1] = clo_v;
    surf_->closestPoint(offset_space_pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, NULL, &seed[0]);

#ifdef DEBUG_CREATORS
    Point normal_to = space_pt + offset_dist_*normal;
    SplineCurve normal_cv(space_pt, normal_to);
    Point other_normal_to = space_pt + other_offset_dist_*other_normal;
    SplineCurve other_normal_cv(space_pt, other_normal_to);
    SplineCurve offset_cv(space_pt, offset_space_pt); 
    vector<Point> pts;
    pts.push_back(clo_pt);
    PointCloud3D pt_cloud(pts[0].begin(), pts.size());
    std::ofstream of("data/debug.g2");
    normal_cv.writeStandardHeader(of);
    normal_cv.write(of);
    other_normal_cv.writeStandardHeader(of);
    other_normal_cv.write(of);
    offset_cv.writeStandardHeader(of);
    offset_cv.write(of);
    pt_cloud.writeStandardHeader(of);
    pt_cloud.write(of);
#endif // DEBUG_CREATORS

    return clo_pt;
}

//===========================================================================
void ProjectIntersectionCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
      der[0] = eval(t);
  else {
      int dim = inters_crv_->dimension();

      // To find the tangent vector we project direction of space curve onto tangent plane in proj pt.
      vector<Point> space_pt = inters_crv_->ParamCurve::point(t, 1);
      // We find parameter values of closest points on respective surfs.
      double clo_u, clo_v, clo_dist, other_clo_u, other_clo_v, other_clo_dist;
      Point clo_pt, other_clo_pt;
      double epsilon = 1e-10;
      vector<double> seed;
      if (p_crv_.get() != 0)
	  seed = getSuggestedSurfaceParameter(*inters_crv_, t, *p_crv_, false);
      // Given a boundary loop, the input of startpt may be interpreted as the endpt.
      // Could be fatal if parametric representation of boundary loop is not a loop.
      // Hence we try with the closestPoint routine, for which the seed is more reliable.
      surf_->closestBoundaryPoint(space_pt[0], clo_u, clo_v, clo_pt, clo_dist,
				  epsgeo_, NULL, &seed[0]);
      vector<double> max_seed_error(2);
      if (seed.size() == 2) {
	  RectDomain under_domain = surf_->containingDomain();
	  max_seed_error[0] = 0.1*(under_domain.umax() - under_domain.umin());
	  max_seed_error[1] = 0.1*(under_domain.vmax() - under_domain.vmin());
	  if ((fabs(clo_u - seed[0]) > max_seed_error[0]) || (fabs(clo_v - seed[1]) > max_seed_error[1])) {
	      surf_->closestPoint(space_pt[0], clo_u, clo_v, clo_pt, clo_dist,
				  epsgeo_, NULL, &seed[0]);
	  }
      }
      if (p_crv_.get() != 0)
	  seed = getSuggestedSurfaceParameter(*inters_crv_, t, *other_p_crv_, true);
      else
	  seed.clear();
      other_surf_->closestBoundaryPoint(space_pt[0], other_clo_u, other_clo_v, other_clo_pt,
					other_clo_dist, epsgeo_, NULL, &seed[0]);
      if (seed.size() == 2) {
	  RectDomain under_domain = other_surf_->containingDomain();
	  max_seed_error[0] = 0.1*(under_domain.umax() - under_domain.umin());
	  max_seed_error[1] = 0.1*(under_domain.vmax() - under_domain.vmin());
	  if ((fabs(other_clo_u - seed[0]) > max_seed_error[0]) ||
	      (fabs(other_clo_v - seed[1]) > max_seed_error[1])) {
	      other_surf_->closestPoint(space_pt[0], other_clo_u, other_clo_v, other_clo_pt,
					other_clo_dist, epsgeo_, NULL, &seed[0]);
	  }
      }

      // @@ We ought to test dist to actual point.
      Point normal, other_normal;
      surf_->normal(normal, clo_u, clo_v);
      other_surf_->normal(other_normal, other_clo_u, other_clo_v);
      Point offset_space_pt = space_pt[0] + offset_dist_*normal + other_offset_dist_*other_normal;
      // We next must project onto surf_. Should use a sensible seed.
      // A better approach would be to project offset pt onto tangent plane in sampled pt.
      seed.resize(2);
      seed[0] = clo_u;
      seed[1] = clo_v;
      surf_->closestPoint(offset_space_pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, NULL, &seed[0]);
      der[0] = clo_pt;

#ifdef DEBUG_CREATORS
      Point normal_to = space_pt[0] + offset_dist_*normal;
      SplineCurve normal_cv(space_pt[0], normal_to);
      Point other_normal_to = space_pt[0] + other_offset_dist_*other_normal;
      SplineCurve other_normal_cv(space_pt[0], other_normal_to);
      SplineCurve offset_cv(space_pt[0], offset_space_pt); 
      vector<Point> pts;
      pts.push_back(clo_pt);
      PointCloud3D pt_cloud(pts[0].begin(), pts.size());
      std::ofstream of("data/debug.g2");
      normal_cv.writeStandardHeader(of);
      normal_cv.write(of);
      other_normal_cv.writeStandardHeader(of);
      other_normal_cv.write(of);
      offset_cv.writeStandardHeader(of);
      offset_cv.write(of);
      pt_cloud.writeStandardHeader(of);
      pt_cloud.write(of);
#endif // DEBUG_CREATORS

      // From the surface we create the derivatives going in the u- and v-direction.
      vector<Point> surf_pts(3);
      surf_->point(surf_pts, clo_u, clo_v, 1);
      // We next describe the dir of space_crv as linear combination of the partial derivs.
      // space_pt[1] = s*surf_pts[1] + t*surf_pts[2].
      double s, t;
      CoonsPatchGen::blendcoef(&surf_pts[1][0], &surf_pts[2][0],
				 &space_pt[1][0], dim, 1, &s, &t);
      der[1] = s*surf_pts[1] + t*surf_pts[2];
  }
}

//===========================================================================
double ProjectIntersectionCurve::start() const
//===========================================================================
{
  return inters_crv_->startparam();
}


//===========================================================================
double ProjectIntersectionCurve::end() const
//===========================================================================
{
  return inters_crv_->endparam();
}

//===========================================================================
int ProjectIntersectionCurve::dim() const
//===========================================================================
{
    return inters_crv_->dimension();
}

//===========================================================================
bool ProjectIntersectionCurve::approximationOK(double par, Point approxpos,
					       double tol1, double tol2) const
//===========================================================================
{
    // Only first tolerance is used.
//   double tol3 = 0.000001*tol1;
  Point pos = eval(par);
  double dist = pos.dist(approxpos);

  // @@sbr Currently no input tolerance is used, only the epsgeo_.
  return (dist < epsgeo_);
}

//===========================================================================
std::vector<double>
ProjectIntersectionCurve::getSuggestedSurfaceParameter(const SplineCurve& space_cv, 
						       double t,
						       const SplineCurve& param_cv,
						       bool pcv_turned) const
//===========================================================================
{
    ALWAYS_ERROR_IF(param_cv.dimension() != 2,
		"Parameter curve expected to have dimension 2...");

    // We compute the distance from start_param as a fraction of total distance.
    double frac = (t - space_cv.startparam())/(space_cv.endparam() - space_cv.startparam());
    double param_t = (pcv_turned) ?
	param_cv.endparam() - (param_cv.endparam() - param_cv.startparam())*frac :
	param_cv.startparam() + (param_cv.endparam() - param_cv.startparam())*frac;
    Point param_pt = param_cv.ParamCurve::point(param_t);

    vector<double> return_vec(param_pt.begin(), param_pt.end());

    return return_vec;
}
