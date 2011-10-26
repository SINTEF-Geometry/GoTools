//===========================================================================
//                                                                           
// File: SmoothTransition.C                                                
//                                                                           
// Created: Mon Dec  9 11:21:48 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SmoothTransition.C,v 1.5 2005-11-17 08:55:31 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/creators/SmoothTransition.h"

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/LineCloud.h"

#include <fstream>

using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::vector;
using std::max;
using std::min;



//===========================================================================
SmoothTransition::SmoothTransition(std::shared_ptr<const SplineCurve>& inters_crv,
				   std::shared_ptr<const SplineCurve>& p_crv1,
				   std::shared_ptr<const SplineCurve>& p_crv2,
				   std::shared_ptr<const ParamSurface> surf1,
				   std::shared_ptr<const ParamSurface> surf2,
				   double offset_dist1, double offset_dist2,
				   double epsgeo)
    : inters_crv_(inters_crv), p_crv1_(p_crv1), p_crv2_(p_crv2),
      surf1_(surf1), surf2_(surf2), offset_dist1_(offset_dist1),
      offset_dist2_(offset_dist2), epsgeo_(epsgeo), kinktol_(1e-02)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(inters_crv_.get() == 0 || p_crv1.get() == 0 || p_crv2.get() == 0 ||
		surf1_.get() == 0 || surf2_.get() == 0,
		"Missing input data.");
    ALWAYS_ERROR_IF(inters_crv_->dimension() != 3 ||
		inters_crv_->dimension() != surf1_->dimension() ||
		surf1_->dimension() != surf2_->dimension() ||
		p_crv1->dimension() != 2 || p_crv2->dimension() != 2,
		"Dimension mismatch.");

    // We create sisl surfaces from the input surfs. May be skipped in the future.
    shared_ptr<SplineSurface> sf1, sf2;
    if (surf1_->instanceType() == Class_SplineSurface)
	under_surf1_ = dynamic_pointer_cast<const SplineSurface, const ParamSurface>(surf1_);
    else if (surf1_->instanceType() == Class_BoundedSurface)
	under_surf1_ = dynamic_pointer_cast<const SplineSurface, const ParamSurface>
	    (dynamic_pointer_cast<const BoundedSurface, const ParamSurface>(surf1_)->underlyingSurface());
    if (surf2_->instanceType() == Class_SplineSurface)
	under_surf2_ = dynamic_pointer_cast<const SplineSurface, const ParamSurface>(surf2_);
    else if (surf2_->instanceType() == Class_BoundedSurface)
	under_surf2_ = dynamic_pointer_cast<const SplineSurface, const ParamSurface>
	    (dynamic_pointer_cast<const BoundedSurface, const ParamSurface>(surf2_)->underlyingSurface());

    int nmb_samples = 10;
    tangent_lengths_.resize(6);
    tangent_lengths_[0] = inters_crv_->ParamCurve::estimatedCurveLength(nmb_samples)/
	(inters_crv_->endparam() - inters_crv_->startparam());
    tangent_lengths_[1] = p_crv1_->ParamCurve::estimatedCurveLength(nmb_samples)/
	(p_crv1_->endparam() - p_crv1_->startparam());
    tangent_lengths_[3] = tangent_lengths_[0];
    tangent_lengths_[4] = p_crv2_->ParamCurve::estimatedCurveLength(nmb_samples)/
	(p_crv2_->endparam() - p_crv2_->startparam());

}

//===========================================================================
SmoothTransition::~SmoothTransition()
//===========================================================================
{
}

//===========================================================================
vector<Point> SmoothTransition::eval(double t)
//===========================================================================
{
    // @@sbr We should somehow figure out direction of cross tangent curves.

    const int nmb_return_pts = 6;
    vector<Point> return_pts(nmb_return_pts);
    int derivs = 1; //2;

    // We compute parameter point on the two underlying surfaces.
    vector<Point> curve_pt = inters_crv_->ParamCurve::point(t, derivs);

    // Modify the plane normal if we are at a surface boundary to avoid
    // iterating to far off the boundary. THIS IS A HACK, vsk 0102.
    if (t == inters_crv_->startparam() || t == inters_crv_->endparam())
      {
	RectDomain domain1 = surf1_->containingDomain();
	RectDomain domain2 = surf2_->containingDomain();

	// Check if the parameter point lies at a boundary of surface 1
	Point ptpar;
	vector<Point> gpt;
	Point normal, vec;
	bool change = false;
	if (p_crv1_.get())
	  {
	    ptpar = p_crv1_->ParamCurve::point(t);
	    surf1_->point(gpt, ptpar[0], ptpar[1], 1);
	    if (fabs(ptpar[0]-domain1.umin()) || 
		fabs(ptpar[0]-domain1.umax()))
	      {
		normal = gpt[1] % curve_pt[1];
		vec = gpt[1] % normal;
		curve_pt[1] = vec;
		change = true;
	      }
	    else if (fabs(ptpar[1]-domain1.vmin()) || 
		     fabs(ptpar[1]-domain1.vmax()))
	      {
		normal = gpt[0] % curve_pt[1];
		vec = gpt[0] % normal;
		curve_pt[1] = vec;
		change = true;
	      }
	  }

	if (!change && p_crv2_.get())
	  {
	    ptpar = p_crv2_->ParamCurve::point(t);
	    surf2_->point(gpt, ptpar[0], ptpar[1], 1);
	    if (fabs(ptpar[0]-domain2.umin()) || 
		fabs(ptpar[0]-domain2.umax()))
	      {
		normal = gpt[1] % curve_pt[1];
		vec = gpt[1] % normal;
		curve_pt[1] = vec;
		change = true;
	      }
	    else if (fabs(ptpar[1]-domain2.vmin()) || 
		     fabs(ptpar[1]-domain2.vmax()))
	      {
		normal = gpt[0] % curve_pt[1];
		vec = gpt[0] % normal;
		curve_pt[1] = vec;
		change = true;
	      }
	  }
      }
    // End hack

    // We next calculate parameter points given by intersection of both surfs and plane.
    Point param_guess_pt1, param_guess_pt2;
    guessParameterPoints(curve_pt[0], t, *inters_crv_, *surf1_, *surf2_, *p_crv1_, *p_crv2_,
			 offset_dist1_, offset_dist2_, param_guess_pt1, param_guess_pt2);

    // We compute derivatives in both surfaces.
    vector<Point> guess_pt1 = surf1_->point(param_guess_pt1[0], param_guess_pt1[1], 2);
    vector<Point> guess_pt2 = surf2_->point(param_guess_pt2[0], param_guess_pt2[1], 2);
    Point normal1, normal2;
    surf1_->normal(normal1, param_guess_pt1[0], param_guess_pt1[1]);
    surf2_->normal(normal2, param_guess_pt2[0], param_guess_pt2[1]);
    guess_pt1.push_back(normal1);
    guess_pt2.push_back(normal2);
    // Given initial values for iteration we locate intersection between sfs and plane.
    // The curve_pt defines the mentioned plane.
    double tstep = 0.0; //= DNULL; // DNULL?
    vector<Point> sipnt1(7), sipnt2(7), sioffpnt1(7), sioffpnt2(7);
    Point sipar1, sipar2;
    offsetIntersectionIterate(offset_dist1_, offset_dist2_, curve_pt,
			      guess_pt1, guess_pt2,
			      param_guess_pt1, param_guess_pt2,
			      *under_surf1_, *under_surf2_, tstep, epsgeo_,
			      sipnt1, sipnt2, sioffpnt1, sioffpnt2,
			      sipar1, sipar2);


#ifdef CREATORS_DEBUG
    std::ofstream of("data/debug.g2");
    SplineCurve normal_cv1(sipnt1[0], sioffpnt1[0]);
    normal_cv1.writeStandardHeader(of);
    normal_cv1.write(of);
    SplineCurve normal_cv2(sipnt2[0], sioffpnt2[0]);
    normal_cv2.writeStandardHeader(of);
    normal_cv2.write(of);
#endif // CREATORS_DEBUG

//     int maxinf = 100;
    vector<Point> s3dinf(3), s3dbb1(3), s3dbb2(3), s3dcrt1(3), s3dcrt2(3);
    vector<Point> sp1inf(3), sp2inf(3);
    vector<double> curv_radis(7); // @@sbr Is radius of curvature really used?
    offsetIntersectionPoints(sipnt1, sipnt2, sioffpnt1, sioffpnt2,
			     sipar1, sipar2, s3dinf, s3dbb1, s3dbb2,
			     s3dcrt1, s3dcrt2, sp1inf, sp2inf, curv_radis);

    // From these parameter points we extract spatial and cross tangent points.
    return_pts[0] = s3dbb1[0]; //surf_pts1[0];
    return_pts[1] = sp1inf[0];
    return_pts[2] = s3dcrt1[0]; //cross_pts1[0];
    return_pts[3] = s3dbb2[0]; //surf_pts2[0];
    return_pts[4] = sp2inf[0];
    return_pts[5] = s3dcrt2[0]; //cross_pts2[0];
//     return_pts[5] *= -1.0; // @@sbr This is not the final solution!!!

#ifdef CREATORS_DEBUG
    std::ofstream of2("data/debug.g2");
    Point cross_to1 = return_pts[0] + return_pts[2];
    SplineCurve cross_cv1(return_pts[0], cross_to1);
    Point cross_to2 = return_pts[3] + return_pts[5];
    SplineCurve cross_cv2(return_pts[3], cross_to2);
    cross_cv1.writeStandardHeader(of2);
    cross_cv1.write(of2);
    cross_cv2.writeStandardHeader(of2);
    cross_cv2.write(of2);
#endif // CREATORS_DEBUG

    return return_pts;
}

//===========================================================================
void SmoothTransition::eval(double t, int n, std::vector<std::vector<Point> >& ders)
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  int derivs = 1;

  const int nmb_return_pts = 6;
  if (int(ders.size()) != nmb_return_pts)
      ders.resize(nmb_return_pts);
  if (n == 0) {
      vector<Point> evals = eval(t);
      for (size_t ki = 0; ki < evals.size(); ++ki) {
	  ders.resize(1);
	  ders[ki].push_back(evals[ki]);
      }
  } else {
      for (size_t ki = 0; ki < ders.size(); ++ki)
	  ders[ki].resize(n+1); // We will compute pos and n der for each curve.

      // @@sbr We should somehow figure out direction of cross tangent curves.

      // We compute parameter point on the two underlying surfaces.
      vector<Point> curve_pt = inters_crv_->ParamCurve::point(t, derivs);

    // Modify the plane normal if we are at a surface boundary to avoid
    // iterating to far off the boundary. THIS IS A HACK, vsk 0102.
    if (t == inters_crv_->startparam() || t == inters_crv_->endparam())
      {
	RectDomain domain1 = surf1_->containingDomain();
	RectDomain domain2 = surf2_->containingDomain();

	// Check if the parameter point lies at a boundary of surface 1
	double eps = 0.000001;
	Point ptpar;
	vector<Point> gpt(3);
	Point normal, vec;
	bool change = false;
	if (p_crv1_.get())
	  {
	    ptpar = p_crv1_->ParamCurve::point(t);
	    surf1_->point(gpt, ptpar[0], ptpar[1], 1);
	    normal = gpt[1] % gpt[2];
	    if (fabs(ptpar[0]-domain1.umin())<eps || 
		fabs(ptpar[0]-domain1.umax())<eps)
	      {
		vec = normal % gpt[2];
		vec.normalize();
		curve_pt[1] = vec;
		change = true;
	      }
	    else if (fabs(ptpar[1]-domain1.vmin())<eps || 
		     fabs(ptpar[1]-domain1.vmax())<eps)
	      {
		vec = normal % gpt[1];
		vec.normalize();
		curve_pt[1] = vec;
		change = true;
	      }
	  }

	if (!change && p_crv2_.get())
	  {
	    ptpar = p_crv2_->ParamCurve::point(t);
	    surf2_->point(gpt, ptpar[0], ptpar[1], 1);
	    normal = gpt[1] % gpt[2];
	    if (fabs(ptpar[0]-domain2.umin())<eps || 
		fabs(ptpar[0]-domain2.umax())<eps)
	      {
		vec = normal % gpt[2];
		vec.normalize();
		curve_pt[1] = vec;
		change = true;
	      }
	    else if (fabs(ptpar[1]-domain2.vmin())<eps || 
		     fabs(ptpar[1]-domain2.vmax())<eps)
	      {
		vec = normal % gpt[1];
		vec.normalize();
		curve_pt[1] = vec;
		change = true;
	      }
	  }
      }
    // End hack

#ifdef CREATORS_DEBUG
    std::ofstream of6("data/debug.g2");
    vector<double> pts;
    pts.insert(pts.end(), curve_pt[0].begin(), curve_pt[0].end());
    Point to_pt3 = curve_pt[0] + curve_pt[1];
    pts.insert(pts.end(), to_pt3.begin(), to_pt3.end());
    LineCloud line_cloud3(&pts[0], 1);
    line_cloud3.writeStandardHeader(of6);
    line_cloud3.write(of6);
#endif // CREATORS_DEBUG

      // We next calculate parameter points given by intersection of both surfs and plane.
      Point param_guess_pt1, param_guess_pt2;
      guessParameterPoints(curve_pt[0], t, *inters_crv_, *surf1_, *surf2_, *p_crv1_, *p_crv2_,
			   offset_dist1_, offset_dist2_, param_guess_pt1, param_guess_pt2);

      // We compute derivatives in both surfaces.
      vector<Point> guess_pt1 = surf1_->point(param_guess_pt1[0], param_guess_pt1[1], 2);
      vector<Point> guess_pt2 = surf2_->point(param_guess_pt2[0], param_guess_pt2[1], 2);
      Point normal1, normal2;
      surf1_->normal(normal1, param_guess_pt1[0], param_guess_pt1[1]);
      surf2_->normal(normal2, param_guess_pt2[0], param_guess_pt2[1]);
      guess_pt1.push_back(normal1);
      guess_pt2.push_back(normal2);
      // Given initial values for iteration we locate intersection between sfs and plane.
      // The curve_pt defines the mentioned plane.
      double tstep = 0.0; //= DNULL; // DNULL?
      vector<Point> sipnt1(7), sipnt2(7), sioffpnt1(7), sioffpnt2(7);
      Point sipar1, sipar2;
      offsetIntersectionIterate(offset_dist1_, offset_dist2_, curve_pt,
				guess_pt1, guess_pt2,
				param_guess_pt1, param_guess_pt2,
				*under_surf1_, *under_surf2_, 
				tstep, epsgeo_, sipnt1, sipnt2, sioffpnt1, 
				sioffpnt2, sipar1, sipar2);

#ifdef CREATORS_DEBUG
      std::ofstream of("data/debug.g2");
      SplineCurve normal_cv1(sipnt1[0], sioffpnt1[0]);
      normal_cv1.writeStandardHeader(of);
      normal_cv1.write(of);
      SplineCurve normal_cv2(sipnt2[0], sioffpnt2[0]);
      normal_cv2.writeStandardHeader(of);
      normal_cv2.write(of);
#endif // CREATORS_DEBUG

      //     int maxinf = 100;
      vector<Point> s3dinf(3), s3dbb1(3), s3dbb2(3), s3dcrt1(3), s3dcrt2(3);
      vector<Point> sp1inf(3), sp2inf(3);
      vector<double> curv_radis(7); // @@sbr Is radius of curvature really used?
      offsetIntersectionPoints(sipnt1, sipnt2, sioffpnt1, sioffpnt2,
			       sipar1, sipar2, s3dinf, s3dbb1, s3dbb2,
			       s3dcrt1, s3dcrt2, sp1inf, sp2inf, curv_radis);

#ifdef CREATORS_DEBUG
      std::ofstream of3("data/debug.g2");
      Point tangent_to1 = s3dbb1[0] + s3dbb1[1];
      SplineCurve tangent_cv1(s3dbb1[0], tangent_to1);
      tangent_cv1.writeStandardHeader(of3);
      tangent_cv1.write(of3);
      Point tangent_to2 = s3dbb2[0] + s3dbb2[1];
      SplineCurve tangent_cv2(s3dbb2[0], tangent_to2);
      tangent_cv2.writeStandardHeader(of3);
      tangent_cv2.write(of3);
#endif // CREATORS_DEBUG

      // From these parameter points we extract spatial and cross tangent points.
      copy(s3dbb1.begin(), s3dbb1.begin() + derivs + 1, ders[0].begin());
      copy(sp1inf.begin(), sp1inf.begin() + derivs + 1, ders[1].begin());
      copy(s3dcrt1.begin(), s3dcrt1.begin() + derivs + 1, ders[2].begin());
      copy(s3dbb2.begin(), s3dbb2.begin() + derivs + 1, ders[3].begin());
      copy(sp2inf.begin(), sp2inf.begin() + derivs + 1, ders[4].begin());
      copy(s3dcrt2.begin(), s3dcrt2.begin() + derivs + 1, ders[5].begin());
//       ders[5][0] *= -1.0; // @@sbr This is not the final solution!!!

      // We need to scale the tangents in the parameter domain. Tangents already nomalized.
      for (int ki = 0; ki < 6; ki += 3) {
	  ders[ki][1] *= tangent_lengths_[ki];
	  ders[ki+1][1] *= tangent_lengths_[ki+1];
      }
      // We test direction of tangents to see if any need to be turned. Expecting this is the case.
      if (curve_pt[1]*ders[0][1] < 0.0) {
	  ders[0][1] *= -1.0; // @@sbr This is not the final solution!!!
	  ders[1][1] *= -1.0; // @@sbr This is not the final solution!!!
	  ders[3][1] *= -1.0; // @@sbr This is not the final solution!!!
	  ders[4][1] *= -1.0; // @@sbr This is not the final solution!!!
      } else {
	  MESSAGE("Was expecting tangent to point in the other direction...");
      }

#ifdef CREATORS_DEBUG
      vector<double> line_pts0, line_pts1, line_pts2;
      line_pts0.insert(line_pts0.end(), s3dinf[0].begin(), s3dinf[0].end());
      Point to_pt = s3dinf[0] + s3dinf[1];
      line_pts0.insert(line_pts0.end(), to_pt.begin(), to_pt.end());
      line_pts1.insert(line_pts1.end(), ders[0][0].begin(), ders[0][0].end());
      to_pt = ders[0][0] + ders[0][1];
      line_pts1.insert(line_pts1.end(), to_pt.begin(), to_pt.end());
      line_pts2.insert(line_pts2.end(), ders[3][0].begin(), ders[3][0].end());
      to_pt = ders[3][0] + ders[3][1];
      line_pts2.insert(line_pts2.end(), to_pt.begin(), to_pt.end());
      LineCloud line_cloud0(&line_pts0[0], 1);
      LineCloud line_cloud1(&line_pts1[0], 1);
      LineCloud line_cloud2(&line_pts2[0], 1);
      std::ofstream of2("data/debug.g2");
      line_cloud0.writeStandardHeader(of2);
      line_cloud0.write(of2);
      line_cloud1.writeStandardHeader(of2);
      line_cloud1.write(of2);
      line_cloud2.writeStandardHeader(of2);
      line_cloud2.write(of2);
#endif // CREATORS_DEBUG

  }
}

//===========================================================================
double SmoothTransition::start()
//===========================================================================
{
  return inters_crv_->startparam();
}


//===========================================================================
double SmoothTransition::end()
//===========================================================================
{
  return inters_crv_->endparam();
}

//===========================================================================
int SmoothTransition::dim()
//===========================================================================
{
    return inters_crv_->dimension();
}

//===========================================================================
bool SmoothTransition::approximationOK(double par, 
				       const std::vector<Point>& approxpos,
				       double tol1, 
				       double tol2)
//===========================================================================
{
    // Only first tolerance is used.
//   double tol3 = 0.000001*tol1;
    vector<Point> pos = eval(par);
    for (size_t ki = 0; ki < approxpos.size(); ki += 3) {
	double dist1 = pos[ki].dist(approxpos[ki]);
	double dist2 = 0.0;//pos[ki+1].dist(approxpos[ki+1]); //@@sbr Fix!
	if ((dist1 > epsgeo_) || (dist2 > epsgeo_))
	    return false;
    }
    // We next test cross tangent curves to see if kink is acceptable.
    for (size_t ki = 2; ki < approxpos.size(); ki += 3) {
	double dist = 0.0; //pos[ki].dist(approxpos[ki]); //@@sbr Fix!
	if (dist > kinktol_)
	    return false;
    }

    // @@sbr Currently no input tolerance is used, only the epsgeo_.
    return true; //(dist < epsgeo_);
}

// Given space point, we project onto surface, returning parameter values
//===========================================================================
Point SmoothTransition::projectPoint(const Point& space_pt, const ParamSurface& surf,
				     std::vector<double>& seed, bool boundary_pt,
				     double epsgeo, double& dist)
//===========================================================================
{
    // We find parameter values of closest point on surf.
    double clo_u, clo_v, clo_dist;
    Point clo_pt, other_clo_pt;
    double epsilon = 1e-5;

    if (boundary_pt) {
	surf.closestBoundaryPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
				  epsgeo_, NULL, &seed[0]);
	if (seed.size() == 2) {
	    // Closest boundary point is not that faithful to seed.
	    RectDomain under_domain = surf.containingDomain();
	    double max_seed_error_u = 0.1*(under_domain.umax() - under_domain.umin());
	    double max_seed_error_v = 0.1*(under_domain.vmax() - under_domain.vmin());
	    if ((fabs(clo_u - seed[0]) > max_seed_error_u) ||
		(fabs(clo_v - seed[1]) > max_seed_error_v)) {
		surf.closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
				  epsgeo_, NULL, &seed[0]);
	    }
	}
    } else
	surf.closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
			  epsilon, NULL, &seed[0]);
    dist = space_pt.dist(clo_pt);

#ifdef CREATORS_DEBUG
    Point cp_pt = space_pt;
    SplineCurve proj_cv(cp_pt, clo_pt);
    std::ofstream of("data/debug.g2");
    proj_cv.writeStandardHeader(of);
    proj_cv.write(of);
#endif // CREATORS_DEBUG

    return Point(clo_u, clo_v);
}




// }

//===========================================================================
std::vector<Point>
SmoothTransition::computeCrosstangentValues(std::vector<Point>& cv_space_pt,
					    std::vector<Point>& cv_local_pt,
					    int derivs)
//===========================================================================
{
    vector<Point> return_vec;

    if (derivs > 1)
	MESSAGE("Only first derivative will be computed.");

    Point diff = cv_space_pt[0] - cv_local_pt[0];
    // We start by computing the cross tangent vector as cross product of tangent and normal:
    Point cross_tan = cv_space_pt[1]%diff; // We will later normalize this vector.

    if (derivs > 0) {
	// We next compute the derivative of the cross tangent curve.
	double cross_length = cross_tan.length();
	double cl3 = cross_length*cross_length*cross_length;
	Point diff2 = cv_space_pt[1] - cv_local_pt[1];
	Point vec1 = diff%cv_space_pt[1]; // Hmm, is this what we want?
	Point vec2 = diff2%cv_space_pt[1]; // Some reparametrization involved here. See blend_s1304.
	vec1 += vec2;
	double tdot = cross_tan*vec1;

	// Finally we compute the tangent of the cross tangent curve.
	Point cross_tan_tan = vec1/cross_length - cross_tan*tdot/cl3;
	cross_tan_tan.normalize(); // @@sbr Suppose cord length parametrization is adequate.
	return_vec.push_back(cross_tan_tan);
    }

    cross_tan.normalize();
    return_vec.insert(return_vec.begin(), cross_tan);

    return return_vec;
}

//===========================================================================
void SmoothTransition::guessParameterPoints(const Point& inters_cv_pt, double t,
					      const SplineCurve& inters_cv,
					      const ParamSurface& sf1,
					      const ParamSurface& sf2,
					      const SplineCurve& p_inters_cv1,
					      const SplineCurve& p_inters_cv2,
					      double offset_dist1, double offset_dist2,
					      Point& guess_pt1, Point& guess_pt2)
//===========================================================================
{
    vector<double> seed;
    Point par_pt1(2), par_pt2(2);
    double dist;
    if (p_crv1_.get() != 0)
      {
	seed = getSuggestedSurfaceParameter(inters_cv, t, *p_crv1_, false);
	par_pt1.setValue(&seed[0]);
      }
    else
      {
    // @@sbr Currently we require the curve to be a bd curve. However the situation
    // may be more general. In addition seed will not used when iterating on bd curves in
    // a closest bd point evaluation.
	par_pt1 = projectPoint(inters_cv_pt, sf1, seed, false, epsgeo_, dist);
	if (dist > epsgeo_) {
	  //     MESSAGE("Projection of point seems to have failed.");
	}
      }
    seed.clear();
    if (p_crv2_.get() != 0)
      {
	seed = getSuggestedSurfaceParameter(inters_cv, t, *p_crv2_, true);
	par_pt2.setValue(&seed[0]);
      }
    else
      {
	par_pt2 = projectPoint(inters_cv_pt, sf2, seed, false, epsgeo_, dist);
	if (dist > epsgeo_) {
	  //     MESSAGE("Projection of point seems to have failed.");
	}
      }

#ifdef CREATORS_DEBUG
    Point debug_pt = sf1.point(par_pt1[0], par_pt1[1]);
    Point debug_normal;
    sf1.normal(debug_normal, par_pt1[0], par_pt1[1]);
    Point normal_to = debug_pt + offset_dist1_*debug_normal;
    SplineCurve normal_cv(debug_pt, normal_to);
    std::ofstream of1("data/debug.g2");
    normal_cv.writeStandardHeader(of1);
    normal_cv.write(of1);
    debug_pt = sf2.point(par_pt2[0], par_pt2[1]);
    sf2.normal(debug_normal, par_pt2[0], par_pt2[1]);
    normal_to = debug_pt + offset_dist2_*debug_normal;
    normal_cv = SplineCurve(debug_pt, normal_to);
    normal_cv.writeStandardHeader(of1);
    normal_cv.write(of1);
#endif // CREATORS_DEBUG

    Point normal1, normal2;
    sf1.normal(normal1, par_pt1[0], par_pt1[1]);
    sf2.normal(normal2, par_pt2[0], par_pt2[1]);
    Point offset_space_pt = inters_cv_pt + offset_dist1_*normal1 + offset_dist2_*normal2;

//     for (ki = 1; ki < offset_space_pt.size(); ++ki)
// 	offset_space_pt[ki] = curve_pt[ki]; // @@ Hmm, is this a valid assumption...?
    // By using the partial derivatives in the respective surfs, we make a qualified guess on
    // new parameter value.
    double tolerance = 1e-05;
    seed = getSuggestedSurfaceParameter(par_pt1, sf1, offset_space_pt, tolerance);
    guess_pt1 = projectPoint(offset_space_pt, sf1,
			     seed, false, epsgeo_, dist);
    seed = getSuggestedSurfaceParameter(par_pt2, sf2, offset_space_pt, tolerance);
    guess_pt2 = projectPoint(offset_space_pt, sf2,
			     seed, false, epsgeo_, dist);



}

//===========================================================================
vector<double>
SmoothTransition::getSuggestedSurfaceParameter(const SplineCurve& space_cv, double t,
						 const SplineCurve& param_cv, bool pcv_turned)
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

//===========================================================================
vector<double>
SmoothTransition::getSuggestedSurfaceParameter(Point& surf_par_pt,
						 const ParamSurface& surf,
						 Point& space_pt,
						 double tolerance)
//===========================================================================
{
    // We first compute the derivatives of the surface in input point.
    vector<Point> surf_pt = surf.point(surf_par_pt[0], surf_par_pt[1], 1);
    int dim = surf.dimension();
    double s, t;
    Point diff_vec = space_pt - surf_pt[0];
    CoonsPatchGen::blendcoef(&surf_pt[1][0], &surf_pt[2][0],
			       &diff_vec[0], dim, 1, &s, &t);

    Vector2D par_vec(surf_par_pt[0] + s, surf_par_pt[1] + t);
    // We next must make sure that we stay inside valid parameter domain.
    const Domain& domain = surf.parameterDomain();
    Vector2D clo_vec;
    domain.closestInDomain(par_vec, clo_vec, tolerance);

#ifdef CREATORS_DEBUG
    Point to_pt = surf.point(clo_vec[0], clo_vec[1]);
    SplineCurve offset_cv(surf_pt[0], space_pt);
    SplineCurve proj_offset_cv(surf_pt[0], to_pt);
    std::ofstream of("data/debug.g2");
    offset_cv.writeStandardHeader(of);
    offset_cv.write(of);
    proj_offset_cv.writeStandardHeader(of);
    proj_offset_cv.write(of);
#endif // CREATORS_DEBUG

    vector<double> return_vec(clo_vec.begin(), clo_vec.end());
    return return_vec;
}
